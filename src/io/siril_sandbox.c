// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * siril_sandbox — Linux sandbox for the Python script run phase.
 *
 * Constrains a Python script's ability to do damage *around* the sirilpy
 * interface using raw CPython (open(), os.system, socket, subprocess, ctypes).
 * Two kernel-enforced defenses are applied to the spawned interpreter and can
 * NOT be lifted from inside the interpreter (native/ctypes code bypasses any
 * Python-level hook):
 *
 *   1. Landlock: confine WRITES to the project working dir + venv + temp; reads
 *      are left broadly allowed (see the read-confinement note below).
 *   2. seccomp-bpf: block internet egress by denying socket(2) for AF_INET /
 *      AF_INET6, while allowing AF_UNIX (Siril's IPC socket is AF_UNIX).
 *
 * The ruleset and seccomp filter are built entirely in the PARENT
 * (siril_sandbox_prepare); the post-fork/pre-exec child_setup only performs
 * async-signal-safe syscalls (prctl / landlock_restrict_self / seccomp /
 * _exit).
 *
 * NOTE on read-confinement being out of scope for this cut: with egress cut
 * (seccomp), a script can still read e.g. ~/.ssh but has nowhere to send it.
 * Read-confinement (a Landlock allow-list of readable paths) is a deliberate
 * later refinement, not an oversight.
 */

/* =====================================================================
 * LIMITATIONS — what this Landlock + seccomp approach CANNOT close.
 * Do NOT mistake this sandbox for full containment. The following holes
 * are known and out of scope for this cut:
 *
 *  - AF_UNIX reach to privileged local daemons. We MUST allow AF_UNIX
 *    for Siril's own IPC socket, and neither seccomp (it cannot see the
 *    connect(2) path/sockaddr) nor Landlock (it does not mediate unix
 *    socket connect) can stop a script from connecting to e.g.
 *    /var/run/docker.sock (→ root), the ssh-agent socket (→ the user's
 *    private keys) or the D-Bus session/system bus. Closing this needs a
 *    mount/network namespace or an LSM (AppArmor) profile.
 *
 *  - X11 keylogging / screenshot for GUI (PyQt) scripts. Any access to
 *    $DISPLAY is full session snooping and cannot be closed while still
 *    handing the script an X11 GUI. Needs Wayland or a nested X server.
 *
 *  - Trusted computing base. Siril's own C IPC handlers (CMD_*, the SHM
 *    transport) run UN-sandboxed; a memory-safety bug there that is
 *    reachable from a script is a full escape. Hardening those is a
 *    separate, ongoing effort.
 *
 *  - Resource limits (RLIMIT_NPROC / RLIMIT_FSIZE below) are generous by
 *    necessity (busy-desktop UID process count; large astro images) and
 *    are DoS-hardening only, NOT a security boundary. They do not contain
 *    total disk or CPU exhaustion — that needs cgroups/quota.
 *
 *  - No kernel-attack-surface reduction. The seccomp filter is a set of
 *    targeted denials, not a syscall allow-list, so a kernel LPE reached
 *    through an otherwise-allowed syscall is NOT mitigated.
 * ===================================================================== */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE   /* O_PATH */
#endif

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/* Windows AppContainer APIs (used in the _WIN32 branch far below) need the
 * Win8/10 API level. The project defines _WIN32_WINNT=0x0600 (Vista) globally,
 * so we MUST raise it and pull in <windows.h> HERE, before "core/siril.h"
 * transitively includes windows.h at the lower level (after which its include
 * guard would freeze the API surface). This mirrors io/FITS_symlink.c, which
 * sets the level before its first include for the same reason. */
#ifdef _WIN32
#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00
#undef WINVER
#define WINVER 0x0A00
#include <windows.h>
#endif

#include "core/siril.h"
#include "io/siril_sandbox.h"

#ifdef __linux__

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/prctl.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <linux/filter.h>
#include <linux/seccomp.h>
#include <linux/audit.h>
#include <linux/unistd.h>

#ifdef HAVE_LANDLOCK
#include <linux/landlock.h>
#endif

#include "core/siril_log.h"
#include "core/siril_app_dirs.h"

/* Fail-open policy switch. When 1, an unavailable sandbox (old kernel, Landlock
 * disabled) results in a plain unsandboxed spawn plus a one-line warning. Flip
 * to 0 to fail-closed once the MVP has bedded in. */
#define SANDBOX_FAIL_OPEN 1

/* ------------------------------------------------------------------ */
/* seccomp BPF: (a) deny process-inspection syscalls (ptrace family),  */
/*     which are a genuine sandbox ESCAPE (same-uid ptrace/process_vm   */
/*     into an un-sandboxed sibling), and (b) in the network-off filter */
/*     restrict socket() to AF_UNIX only (Siril IPC).                   */
/* Built as static const filters selected per build arch so we never   */
/* ship a filter that silently no-ops on an unexpected architecture.   */
/* Two filters are provided and selected in child_setup on the stored  */
/* allow_network flag; seccomp is now installed WHENEVER the arch is    */
/* supported (regardless of allow_network) because it also performs the */
/* process-inspection denial that is independent of networking.        */
/* ------------------------------------------------------------------ */

#if defined(__x86_64__)
#  define SANDBOX_AUDIT_ARCH AUDIT_ARCH_X86_64
#  define SANDBOX_HAVE_ARCH 1
#elif defined(__aarch64__)
#  define SANDBOX_AUDIT_ARCH AUDIT_ARCH_AARCH64
#  define SANDBOX_HAVE_ARCH 1
#else
#  define SANDBOX_HAVE_ARCH 0
#endif

#ifndef AF_UNIX
#define AF_UNIX   1
#endif

/* Syscall-number fallbacks. To keep the BPF instruction layout FIXED (so
 * #ifdefs never shift jump offsets), every compare instruction is always
 * emitted; a syscall number absent from the build headers is given a value
 * that can never match a real nr (0xffffffffu), so its compare is a harmless
 * no-op instead of vanishing and moving every jt/jf below it. */
#ifdef __NR_ptrace
#  define SBX_NR_PTRACE            __NR_ptrace
#else
#  define SBX_NR_PTRACE            0xffffffffu
#endif
#ifdef __NR_process_vm_readv
#  define SBX_NR_PROCESS_VM_READV  __NR_process_vm_readv
#else
#  define SBX_NR_PROCESS_VM_READV  0xffffffffu
#endif
#ifdef __NR_process_vm_writev
#  define SBX_NR_PROCESS_VM_WRITEV __NR_process_vm_writev
#else
#  define SBX_NR_PROCESS_VM_WRITEV 0xffffffffu
#endif
#ifdef __NR_process_madvise
#  define SBX_NR_PROCESS_MADVISE   __NR_process_madvise
#else
#  define SBX_NR_PROCESS_MADVISE   0xffffffffu
#endif
#ifdef __NR_kcmp
#  define SBX_NR_KCMP              __NR_kcmp
#else
#  define SBX_NR_KCMP              0xffffffffu
#endif
#ifdef __NR_pidfd_getfd
#  define SBX_NR_PIDFD_GETFD       __NR_pidfd_getfd
#else
#  define SBX_NR_PIDFD_GETFD       0xffffffffu
#endif

#if SANDBOX_HAVE_ARCH
/* ---- network-off filter (default policy) --------------------------
 * Layout (indices are load-bearing for the jt/jf offsets — count them):
 *   [0]  load arch
 *   [1]  arch == build arch ? ->[3] : ->[2]                (jt=1, jf=0)
 *   [2]  RET_KILL                                         (arg-layout guard)
 *   [3]  load nr
 *   [4]  nr == ptrace            ? ->EPERM[15] : next      (jt=10)
 *   [5]  nr == process_vm_readv  ? ->EPERM[15] : next      (jt=9)
 *   [6]  nr == process_vm_writev ? ->EPERM[15] : next      (jt=8)
 *   [7]  nr == process_madvise   ? ->EPERM[15] : next      (jt=7)
 *   [8]  nr == kcmp              ? ->EPERM[15] : next      (jt=6)
 *   [9]  nr == pidfd_getfd       ? ->EPERM[15] : next      (jt=5)
 *   [10] nr == socket ? next[11] : ->ALLOW[14]             (jt=0, jf=3)
 *   [11] load args[0] (domain)
 *   [12] domain == AF_UNIX ? ->ALLOW[14] : ->EAFNOSUPPORT[13] (jt=1, jf=0)
 *   [13] RET_ERRNO(EAFNOSUPPORT)
 *   [14] RET_ALLOW
 *   [15] RET_ERRNO(EPERM)
 * socket() is a direct syscall on both x86_64 and aarch64.
 */
static const struct sock_filter sandbox_seccomp_insns[] = {
	/* [0] */ BPF_STMT(BPF_LD  | BPF_W | BPF_ABS, offsetof(struct seccomp_data, arch)),
	/* [1] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SANDBOX_AUDIT_ARCH, 1, 0),
	/* [2] */ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_KILL_PROCESS),
	/* [3] */ BPF_STMT(BPF_LD  | BPF_W | BPF_ABS, offsetof(struct seccomp_data, nr)),
	/* [4] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PTRACE,            10, 0),
	/* [5] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PROCESS_VM_READV,   9, 0),
	/* [6] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PROCESS_VM_WRITEV,  8, 0),
	/* [7] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PROCESS_MADVISE,    7, 0),
	/* [8] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_KCMP,               6, 0),
	/* [9] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PIDFD_GETFD,        5, 0),
	/* [10]*/ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, __NR_socket,               0, 3),
	/* [11]*/ BPF_STMT(BPF_LD  | BPF_W | BPF_ABS, offsetof(struct seccomp_data, args[0])),
	/* [12]*/ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, AF_UNIX,                   1, 0),
	/* [13]*/ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_ERRNO | (EAFNOSUPPORT & SECCOMP_RET_DATA)),
	/* [14]*/ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_ALLOW),
	/* [15]*/ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_ERRNO | (EPERM & SECCOMP_RET_DATA)),
};

/* ---- network-on filter (future allow_network == TRUE) -------------
 * Keeps the process-inspection (ptrace family) denials but does NOT
 * restrict socket() families (inet permitted).
 *   [0]  load arch
 *   [1]  arch == build arch ? ->[3] : ->[2]                (jt=1, jf=0)
 *   [2]  RET_KILL
 *   [3]  load nr
 *   [4]  nr == ptrace            ? ->EPERM[11] : next      (jt=6)
 *   [5]  nr == process_vm_readv  ? ->EPERM[11] : next      (jt=5)
 *   [6]  nr == process_vm_writev ? ->EPERM[11] : next      (jt=4)
 *   [7]  nr == process_madvise   ? ->EPERM[11] : next      (jt=3)
 *   [8]  nr == kcmp              ? ->EPERM[11] : next      (jt=2)
 *   [9]  nr == pidfd_getfd       ? ->EPERM[11] : next      (jt=1)
 *   [10] RET_ALLOW
 *   [11] RET_ERRNO(EPERM)
 */
static const struct sock_filter sandbox_seccomp_insns_net[] = {
	/* [0] */ BPF_STMT(BPF_LD  | BPF_W | BPF_ABS, offsetof(struct seccomp_data, arch)),
	/* [1] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SANDBOX_AUDIT_ARCH, 1, 0),
	/* [2] */ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_KILL_PROCESS),
	/* [3] */ BPF_STMT(BPF_LD  | BPF_W | BPF_ABS, offsetof(struct seccomp_data, nr)),
	/* [4] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PTRACE,            6, 0),
	/* [5] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PROCESS_VM_READV,  5, 0),
	/* [6] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PROCESS_VM_WRITEV, 4, 0),
	/* [7] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PROCESS_MADVISE,   3, 0),
	/* [8] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_KCMP,              2, 0),
	/* [9] */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, SBX_NR_PIDFD_GETFD,       1, 0),
	/* [10]*/ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_ALLOW),
	/* [11]*/ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_ERRNO | (EPERM & SECCOMP_RET_DATA)),
};
#endif /* SANDBOX_HAVE_ARCH */

/* Resource limits (DoS hardening only, NON-fatal — a failed setrlimit is
 * ignored, never _exit). Tunable via these #defines. */
/* RLIMIT_NPROC caps the real-UID's TOTAL live process count (not just this
 * subtree), so it must be generous or the child may be unable to fork at all
 * on a busy desktop. */
#define SANDBOX_MAX_PROCS       4096
/* RLIMIT_FSIZE caps a single file's size; exceeding it raises SIGXFSZ and
 * kills the writer. Astro images/stacks are large, so this is VERY generous
 * (64 GiB). It only stops the single-giant-file case; total-disk exhaustion
 * needs cgroups/quota (out of scope). */
#define SANDBOX_MAX_FSIZE_BYTES (64ULL << 30)

/* PORTABILITY NOTE — OpenBSD (future, not implemented).
 * OpenBSD's unveil(2) is a close analogue of the Landlock write-confinement
 * below (a per-path filesystem allow-list that survives execve), and pledge(2)
 * — with its execpromises argument — is the analogue of the seccomp syscall
 * restriction. Both are plain syscalls callable from a GSpawnChildSetupFunc
 * exactly like the Linux child_setup here, so an #elif defined(__OpenBSD__)
 * branch could provide an equivalent per-child sandbox with modest effort.
 * Other non-Linux/macOS POSIX systems currently fall through to the
 * unconfined generic branch (see the bottom of this file): FreeBSD Capsicum is
 * all-or-nothing capability mode that needs a capsicum-aware program (would
 * break a stock interpreter); jails are root-only/heavyweight; NetBSD/HURD have
 * no practical per-process mechanism. Left unimplemented deliberately — could
 * not be built or tested in the dev environment. */

/* ------------------------------------------------------------------ */
/* Landlock thin syscall wrappers (glibc provides no wrappers).       */
/* ------------------------------------------------------------------ */

#ifdef HAVE_LANDLOCK

#ifndef landlock_create_ruleset
static inline long sandbox_landlock_create_ruleset(
		const struct landlock_ruleset_attr *attr,
		size_t size, __u32 flags) {
	return syscall(__NR_landlock_create_ruleset, attr, size, flags);
}
#endif

#ifndef landlock_add_rule
static inline long sandbox_landlock_add_rule(
		int ruleset_fd, enum landlock_rule_type rule_type,
		const void *rule_attr, __u32 flags) {
	return syscall(__NR_landlock_add_rule, ruleset_fd, rule_type,
	               rule_attr, flags);
}
#endif

#ifndef landlock_restrict_self
static inline long sandbox_landlock_restrict_self(int ruleset_fd, __u32 flags) {
	return syscall(__NR_landlock_restrict_self, ruleset_fd, flags);
}
#endif

/* Handled (restricted) write-ish access rights, masked to ABI in prepare(). */
#define SANDBOX_ACCESS_WRITE_BASE ( \
	LANDLOCK_ACCESS_FS_WRITE_FILE | \
	LANDLOCK_ACCESS_FS_REMOVE_DIR | \
	LANDLOCK_ACCESS_FS_REMOVE_FILE | \
	LANDLOCK_ACCESS_FS_MAKE_CHAR | \
	LANDLOCK_ACCESS_FS_MAKE_DIR | \
	LANDLOCK_ACCESS_FS_MAKE_REG | \
	LANDLOCK_ACCESS_FS_MAKE_SOCK | \
	LANDLOCK_ACCESS_FS_MAKE_FIFO | \
	LANDLOCK_ACCESS_FS_MAKE_BLOCK | \
	LANDLOCK_ACCESS_FS_MAKE_SYM)

/* Write rights that are meaningful on a NON-directory path (regular file,
 * device node, ...). The directory-scoped rights (REMOVE_*, MAKE_*, REFER)
 * only apply when the rule's target is a directory; adding them on a file
 * target makes landlock_add_rule fail with EINVAL. WRITE_FILE and TRUNCATE
 * (ABI≥3) are the only file-applicable write rights. */
#define SANDBOX_ACCESS_FILE_ONLY \
	(LANDLOCK_ACCESS_FS_WRITE_FILE | LANDLOCK_ACCESS_FS_TRUNCATE)

#endif /* HAVE_LANDLOCK */

/* ------------------------------------------------------------------ */
/* Sandbox object — POD only, read by the post-fork child.            */
/*                                                                     */
/* IMPORTANT: the Landlock ruleset is built in the CHILD (child_setup) */
/* and NOT in the parent. GLib closes every descriptor >= 3 that it    */
/* does not manage BEFORE calling child_setup (see the                 */
/* GSpawnChildSetupFunc contract), so a parent-created ruleset fd      */
/* would already be closed by the time child_setup runs, and           */
/* landlock_restrict_self() would fail with EBADF. Instead the parent  */
/* only computes the access mask and the list of writable paths; the   */
/* child creates the ruleset, adds the rules, and restricts itself,    */
/* using only async-signal-safe syscalls (open/fstat/close/syscall).   */
/* ------------------------------------------------------------------ */

/* Internal to the Linux implementation (no longer exposed in the header, which
 * now offers only siril_sandbox_spawn()). */
typedef struct _SirilSandbox SirilSandbox;
struct _SirilSandbox {
	int install_seccomp;               /* 1 if a seccomp filter applies     */
	int allow_network;                 /* select net-on vs net-off filter   */
#ifdef HAVE_LANDLOCK
	guint64 landlock_access;           /* handled write mask, 0 => disabled */
	char **paths;                      /* heap array of g_strdup'd writable  */
	int n_paths;                       /* paths (built-in roots + extras)    */
#endif
};

/* Number of built-in writable roots the Linux ruleset always grants (wd, venv,
 * user-data, <config>/siril, $TMPDIR, /var/tmp, /dev/{null,zero,full}). */
#define SANDBOX_N_BUILTIN_PATHS 9

#ifdef HAVE_LANDLOCK
/* Add `path` to the ruleset with the handled write rights. ASYNC-SIGNAL-SAFE:
 * only open/fstat/close/syscall. Missing paths are skipped silently. A failed
 * add is skipped (non-fatal): the path simply stays non-writable rather than
 * discarding the whole ruleset. Called only from the post-fork child. */
static void sandbox_child_add_path(int ruleset_fd, guint64 access,
                                   const char *path) {
	int fd;

	if (!path || !*path)
		return;

	fd = open(path, O_PATH | O_CLOEXEC);
	if (fd < 0)
		return;   /* does not exist / not accessible: nothing to grant */

	/* Directory-scoped rights (REMOVE_*, MAKE_*, REFER) only apply when the
	 * rule target is a directory; adding them on a regular file or device
	 * node (e.g. /dev/null) makes landlock_add_rule fail with EINVAL. For a
	 * non-directory target, keep only the file-applicable write rights. */
	{
		struct stat st;
		if (fstat(fd, &st) == 0 && !S_ISDIR(st.st_mode))
			access &= (guint64) SANDBOX_ACCESS_FILE_ONLY;
	}

	{
		struct landlock_path_beneath_attr pb = {
			.allowed_access = access,
			.parent_fd = fd,
		};
		/* Ignore failure: a single path that cannot be granted is skipped,
		 * it must not tear down the whole ruleset. */
		syscall(__NR_landlock_add_rule, ruleset_fd,
		        LANDLOCK_RULE_PATH_BENEATH, &pb, 0);
	}
	close(fd);
}
#endif /* HAVE_LANDLOCK */

static SirilSandbox *sandbox_prepare(const SirilSandboxPolicy *policy) {
	int have_seccomp = 0;
	int have_landlock = 0;
	SirilSandbox *sb;

#if SANDBOX_HAVE_ARCH
	/* seccomp is available on this build arch. It installs REGARDLESS of
	 * allow_network, because it also performs the process-inspection
	 * (ptrace family) denial that is independent of networking; the
	 * allow_network flag only selects which of the two filters is used. */
	have_seccomp = 1;
#else
	/* Unknown build arch: do not attempt to ship a filter. */
	have_seccomp = 0;
#endif

	sb = g_malloc0(sizeof(SirilSandbox));
	sb->allow_network = policy->allow_network ? 1 : 0;

#ifdef HAVE_LANDLOCK
	{
		int abi = (int) syscall(__NR_landlock_create_ruleset, NULL, 0,
		                        LANDLOCK_CREATE_RULESET_VERSION);
		if (abi >= 1) {
			guint64 access = SANDBOX_ACCESS_WRITE_BASE;
			const char *tmpdir;
			int n = 0, n_extra = 0;

			/* Mask rights down to what the detected ABI supports:
			 * requesting unsupported bits yields EINVAL. */
			if (abi >= 2)
				access |= LANDLOCK_ACCESS_FS_REFER;
			if (abi >= 3)
				access |= LANDLOCK_ACCESS_FS_TRUNCATE;
			sb->landlock_access = access;

			/* Size the path array for the built-in roots plus any per-script
			 * extra writable paths from the policy/manifest. */
			if (policy->extra_write_paths)
				for (const char * const *p = policy->extra_write_paths; *p; p++)
					n_extra++;
			sb->paths = g_new0(char *, SANDBOX_N_BUILTIN_PATHS + n_extra);

			/* Assemble the writable path list (resolved here in the parent
			 * so the child needs no g_getenv). Writes are confined to the
			 * project working dir, the venv (python writes __pycache__/.pyc
			 * there), temp dirs and the standard write-safe char devices.
			 * Reads are NOT handled and so remain allowed everywhere. */
			if (policy->wd && *policy->wd)
				sb->paths[n++] = g_strdup(policy->wd);
			if (policy->venv_path && *policy->venv_path)
				sb->paths[n++] = g_strdup(policy->venv_path);

			/* Siril's own writable directories, via the app-dirs accessors
			 * (platform/XDG-correct — do not hardcode): user data
			 * (~/.local/share/siril: prefs, venv metadata) and the per-app
			 * config dir (<config>/siril — NOT the whole XDG config base, which
			 * would expose every app's config). The scripts repo and SPCC data
			 * dir are deliberately NOT granted: they are git-managed by Siril's
			 * own C code, and since reads are unrestricted (read rights are
			 * unhandled) scripts can still read them — they just cannot write.
			 * A dir that does not exist yet is simply skipped when the rule is
			 * added in the child. */
			{
				const char *d;
				gchar *cfg;
				d = siril_get_user_data_dir();
				if (d && *d)
					sb->paths[n++] = g_strdup(d);
				cfg = g_build_filename(siril_get_config_dir(), PACKAGE, NULL);
				if (cfg && *cfg)
					sb->paths[n++] = cfg;   /* heap-owned; freed in finish() */
				else
					g_free(cfg);
			}

			tmpdir = g_getenv("TMPDIR");
			sb->paths[n++] = g_strdup((tmpdir && *tmpdir) ? tmpdir : "/tmp");
			sb->paths[n++] = g_strdup("/var/tmp");
			sb->paths[n++] = g_strdup("/dev/null");
			sb->paths[n++] = g_strdup("/dev/zero");
			sb->paths[n++] = g_strdup("/dev/full");

			/* Per-script extra writable paths (from the granted manifest).
			 * extra_read_paths is a no-op on Linux: reads are unrestricted. */
			if (policy->extra_write_paths)
				for (const char * const *p = policy->extra_write_paths; *p; p++)
					if (**p)
						sb->paths[n++] = g_strdup(*p);

			sb->n_paths = n;

			have_landlock = 1;
		}
		/* abi < 0 => ENOSYS/EOPNOTSUPP (kernel < 5.13 or Landlock off). */
	}
#endif /* HAVE_LANDLOCK */

	if (!have_landlock && !have_seccomp) {
		/* Neither confinement could be prepared. */
		g_free(sb);
#if SANDBOX_FAIL_OPEN
		siril_log_message(_("Script sandbox unavailable on this kernel; "
		                    "running Python without filesystem/network "
		                    "confinement.\n"));
#endif
		return NULL;
	}

	sb->install_seccomp = have_seccomp;
	return sb;
}

/* Runs in the child after fork, before exec. ASYNC-SIGNAL-SAFE: only
 * open/fstat/close/prctl/syscall/seccomp/_exit. No g_*, malloc, logging, or
 * errno message formatting. Failures surface as the child exiting 127. The
 * ruleset is built HERE (not in the parent) because GLib has already closed
 * any parent-created fd by this point. */
static void sandbox_child_setup(gpointer user_data) {
	const SirilSandbox *sb = (const SirilSandbox *) user_data;

	if (!sb)
		return;

	/* PR_SET_NO_NEW_PRIVS is required for both unprivileged seccomp and
	 * Landlock. */
	if (prctl(PR_SET_NO_NEW_PRIVS, 1, 0, 0, 0) != 0)
		_exit(127);

#ifdef HAVE_LANDLOCK
	if (sb->landlock_access != 0) {
		struct landlock_ruleset_attr attr = {
			.handled_access_fs = sb->landlock_access,
		};
		int ruleset_fd = (int) syscall(__NR_landlock_create_ruleset,
		                               &attr, sizeof(attr), 0);
		if (ruleset_fd < 0)
			_exit(127);

		for (int i = 0; i < sb->n_paths; i++)
			sandbox_child_add_path(ruleset_fd, sb->landlock_access,
			                       sb->paths[i]);

		if (syscall(__NR_landlock_restrict_self, ruleset_fd, 0) != 0)
			_exit(127);
		close(ruleset_fd);
	}
#endif

	/* Resource limits (DoS hardening only). NON-fatal: a failed setrlimit
	 * must never break a legit script, so its return is ignored and we do
	 * NOT _exit. See the SANDBOX_MAX_* comments for the generous caps. */
	{
		struct rlimit rl;
		rl.rlim_cur = SANDBOX_MAX_PROCS;
		rl.rlim_max = SANDBOX_MAX_PROCS;
		(void) setrlimit(RLIMIT_NPROC, &rl);

		rl.rlim_cur = SANDBOX_MAX_FSIZE_BYTES;
		rl.rlim_max = SANDBOX_MAX_FSIZE_BYTES;
		(void) setrlimit(RLIMIT_FSIZE, &rl);
	}

#if SANDBOX_HAVE_ARCH
	if (sb->install_seccomp) {
		struct sock_fprog prog;
		/* Select the network-off (AF_UNIX-only + ptrace denials) or
		 * network-on (ptrace denials only) filter. Both are static const
		 * arrays in parent memory that survive fork; reading them here is
		 * async-signal-safe. */
		if (sb->allow_network) {
			prog.len = (unsigned short) (sizeof(sandbox_seccomp_insns_net) /
			                             sizeof(sandbox_seccomp_insns_net[0]));
			prog.filter = (struct sock_filter *) sandbox_seccomp_insns_net;
		} else {
			prog.len = (unsigned short) (sizeof(sandbox_seccomp_insns) /
			                             sizeof(sandbox_seccomp_insns[0]));
			prog.filter = (struct sock_filter *) sandbox_seccomp_insns;
		}
		if (syscall(__NR_seccomp, SECCOMP_SET_MODE_FILTER, 0, &prog) != 0)
			_exit(127);
	}
#endif
}

static void sandbox_finish(SirilSandbox *sb) {
	if (!sb)
		return;
#ifdef HAVE_LANDLOCK
	if (sb->paths) {
		for (int i = 0; i < sb->n_paths; i++)
			g_free(sb->paths[i]);
		g_free(sb->paths);
	}
#endif
	g_free(sb);
}

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const SirilSandboxPolicy *policy, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	/* Explicit, granted opt-out (glue scripts spawning opaque third-party
	 * software): no confinement at all. */
	if (policy->unsandboxed)
		return g_spawn_async_with_pipes(working_dir, argv, envp,
			G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
			NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);

	/* Build the ruleset spec in the parent; the actual Landlock ruleset and
	 * seccomp filter are installed by sandbox_child_setup() post-fork. If
	 * preparation fails entirely (unsupported kernel) sb is NULL and we spawn
	 * without a child_setup (fail-open). */
	SirilSandbox *sb = sandbox_prepare(policy);
	gboolean ok = g_spawn_async_with_pipes(
		working_dir, argv, envp,
		G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
		sb ? sandbox_child_setup : NULL, sb,
		child_pid, NULL, stdout_fd, stderr_fd, error);
	sandbox_finish(sb);
	return ok;
}

#elif defined(__APPLE__)

/* ==================================================================== */
/* macOS: confine via /usr/bin/sandbox-exec with a generated Seatbelt   */
/* (SBPL) profile prepended to the argv.                                */
/*                                                                       */
/* Why argv-wrap and not sandbox_init() in a child_setup: sandbox_init() */
/* is not async-signal-safe (it parses/compiles the profile and         */
/* allocates), so it is unsafe to call between fork() and exec() in the  */
/* multithreaded Siril process. sandbox-exec applies the profile in a    */
/* separate exec, side-stepping that entirely.                          */
/*                                                                       */
/* NB sandbox-exec is Apple-deprecated but remains the only mechanism    */
/* for a DYNAMIC, per-child, path-scoped profile — App Sandbox is        */
/* app-wide and static (set at code-signing), so it cannot express       */
/* "confine just this Python child with these dirs". The mechanism       */
/* (Seatbelt) underpins App Sandbox itself and is stable in practice.    */
/*                                                                       */
/* LIMITATION: this SBPL profile is a best-effort draft that could not   */
/* be validated on macOS in the build environment; it needs a runtime    */
/* test on a Mac (the CI macos job only packages, it does not run a      */
/* script under the profile). Same residual limits as Linux apply        */
/* (X11/Wayland n/a on mac; the C IPC is still the TCB; local IPC        */
/* sockets to privileged daemons are not filtered here).                 */
/* ==================================================================== */

#include "core/siril_log.h"
#include "core/siril_app_dirs.h"

#define SANDBOX_EXEC_PATH "/usr/bin/sandbox-exec"

/* Append `s` to the profile as an SBPL double-quoted string literal, escaping
 * backslashes and double quotes. */
static void sb_append_quoted(GString *p, const char *s) {
	g_string_append_c(p, '"');
	for (const char *c = s; *c; c++) {
		if (*c == '\\' || *c == '"')
			g_string_append_c(p, '\\');
		g_string_append_c(p, *c);
	}
	g_string_append_c(p, '"');
}

/* Append (allow <op> (subpath "PATH")) if path is non-empty. */
static void sb_allow_subpath(GString *p, const char *op, const char *path) {
	if (!path || !*path)
		return;
	g_string_append_printf(p, "(allow %s (subpath ", op);
	sb_append_quoted(p, path);
	g_string_append(p, "))");
}

/* Build a self-contained Seatbelt profile string (passed via sandbox-exec -p).
 * Allow-by-default (so the interpreter reads its libs freely), then confine:
 * deny all writes except the granted roots, and deny IP network egress while
 * leaving AF_UNIX (Siril IPC) untouched. */
static gchar *build_seatbelt_profile(const SirilSandboxPolicy *policy) {
	GString *p = g_string_new("(version 1)(allow default)");

	/* Write confinement: global deny, then re-allow the granted roots (SBPL
	 * is last-match-wins, so the allows must follow the deny). */
	g_string_append(p, "(deny file-write*)");
	sb_allow_subpath(p, "file-write*", policy->wd);
	sb_allow_subpath(p, "file-write*", policy->venv_path);
	sb_allow_subpath(p, "file-write*", siril_get_user_data_dir());
	{
		gchar *cfg = g_build_filename(siril_get_config_dir(), PACKAGE, NULL);
		sb_allow_subpath(p, "file-write*", cfg);
		g_free(cfg);
	}
	sb_allow_subpath(p, "file-write*", g_get_tmp_dir());
	sb_allow_subpath(p, "file-write*", "/dev/null");
	sb_allow_subpath(p, "file-write*", "/dev/random");
	sb_allow_subpath(p, "file-write*", "/dev/urandom");

	/* Per-script extra writable paths from the granted manifest. extra_read
	 * paths are a near no-op here (reads are already allowed by default) but
	 * are granted explicitly for symmetry with the Windows read-confined case. */
	if (policy->extra_write_paths)
		for (const char * const *e = policy->extra_write_paths; *e; e++)
			sb_allow_subpath(p, "file-write*", *e);
	if (policy->extra_read_paths)
		for (const char * const *e = policy->extra_read_paths; *e; e++)
			sb_allow_subpath(p, "file-read*", *e);

	/* Network: deny IP egress (exfil). AF_UNIX stays allowed by (allow
	 * default), so the Siril IPC socket is unaffected. */
	if (!policy->allow_network)
		g_string_append(p, "(deny network-outbound (remote ip \"*:*\"))");

	/* Deny reading other processes' state (ptrace/task_for_pid analogue). */
	g_string_append(p, "(deny process-info*)");

	return g_string_free(p, FALSE);
}

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const SirilSandboxPolicy *policy, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	/* Explicit, granted opt-out (glue scripts spawning opaque third-party
	 * software): no confinement at all. */
	if (policy->unsandboxed)
		return g_spawn_async_with_pipes(working_dir, argv, envp,
			G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
			NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);

	/* If sandbox-exec is missing (should never happen on macOS), fall back to
	 * an unsandboxed spawn with a warning rather than failing outright. */
	if (!g_file_test(SANDBOX_EXEC_PATH, G_FILE_TEST_IS_EXECUTABLE)) {
		siril_log_message(_("Script sandbox unavailable (%s missing); running "
		                    "Python without confinement.\n"), SANDBOX_EXEC_PATH);
		return g_spawn_async_with_pipes(working_dir, argv, envp,
			G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
			NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);
	}

	gchar *profile = build_seatbelt_profile(policy);

	/* Prepend: sandbox-exec -p <profile> <original argv...> */
	GPtrArray *wrapped = g_ptr_array_new();
	g_ptr_array_add(wrapped, (gpointer) SANDBOX_EXEC_PATH);
	g_ptr_array_add(wrapped, (gpointer) "-p");
	g_ptr_array_add(wrapped, profile);
	for (gchar **a = argv; a && *a; a++)
		g_ptr_array_add(wrapped, *a);
	g_ptr_array_add(wrapped, NULL);

	gboolean ok = g_spawn_async_with_pipes(working_dir,
		(gchar **) wrapped->pdata, envp,
		G_SPAWN_DO_NOT_REAP_CHILD,   /* absolute path, no SEARCH_PATH needed */
		NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);

	g_ptr_array_free(wrapped, FALSE);
	g_free(profile);
	return ok;
}

#elif defined(_WIN32)

/* ==================================================================== */
/* Windows: AppContainer confinement.                                   */
/*                                                                       */
/* AppContainer cannot be applied through g_spawn: it needs a           */
/* SECURITY_CAPABILITIES attribute passed to CreateProcess at creation  */
/* time. So this path owns process creation and reproduces g_spawn's    */
/* contract: it fills *child_pid (the process HANDLE cast to GPid),     */
/* *stdout_fd / *stderr_fd (CRT fds via _open_osfhandle), and returns   */
/* TRUE / FALSE + g_set_error.                                          */
/*                                                                       */
/* Defenses:                                                            */
/*   - No network capability in the container => no internetClient etc.,  */
/*     so network egress is denied by construction (nothing to un-set).  */
/*     removableStorage is a FILE capability and does not change this.   */
/*   - AppContainer is deny-by-default for the user's files, so the      */
/*     child cannot even load its interpreter until we add explicit      */
/*     access-allowed ACEs for the container SID. To MATCH the Linux/    */
/*     macOS model (reads permitted everywhere; exfiltration blocked by  */
/*     the network cut, not by read confinement), we grant read/execute  */
/*     broadly — the user profile, every fixed drive root, and removable */
/*     media via the removableStorage capability. Writes stay confined.  */
/*                                                                       */
/* HONESTY / LIMITATIONS (this file could NOT be compiled or run in the  */
/* dev environment — no MinGW / Windows):                                */
/*   - The broad read grants (user profile + fixed drive roots) should    */
/*     now cover the interpreter/stdlib wherever a uv-managed Python lives */
/*     (under the user profile or a fixed drive), so the earlier "base    */
/*     Python install not readable" risk should be resolved — but this,   */
/*     like everything here, is UNVERIFIED without a Windows toolchain.   */
/*   - The C:\ root DACL edit typically fails without elevation and is     */
/*     skipped; the user-profile grant covers C:\Users\<user>, but reads   */
/*     of C:\ locations OUTSIDE the profile, and files with protected/     */
/*     explicit ACLs that do not inherit, may still be denied.            */
/*   - The IPC pipe SD (companion change in siril_pythonmodule.c) grants  */
/*     ALL APPLICATION PACKAGES; that too wants runtime confirmation.     */
/*   - The granted ACEs PERSIST on disk (they name only this one          */
/*     container SID). The MVP leaves them in place across runs.          */
/* ==================================================================== */

/* <windows.h> is already included at the top of this file at the required
 * Win8/10 API level (see the _WIN32_WINNT bump there). Only the extra
 * headers this branch needs are pulled in here. */
#include <userenv.h>   /* CreateAppContainerProfile / DeriveAppContainerSid... */
#include <aclapi.h>    /* GetNamedSecurityInfo / SetEntriesInAcl / SetNamed...  */
#include <sddl.h>      /* ConvertStringSidToSid (internetClient capability SID) */
#include <io.h>        /* _open_osfhandle                                       */
#include <fcntl.h>     /* _O_RDONLY                                             */

#include "core/siril_log.h"
#include "core/siril_app_dirs.h"

/* Fixed container name/identity. A stable name means the profile is created
 * once and re-derived on subsequent runs. */
#define SANDBOX_APPCONTAINER_NAME L"org.siril.script_sandbox"
#define SANDBOX_APPCONTAINER_DESC L"Siril Python script sandbox"

/* ---- MinGW-w64 header workarounds ---------------------------------------
 * Some AppContainer bits are absent or inconsistent across MinGW-w64 header
 * vintages. Declare/guard each locally so the CI build does not depend on a
 * particular MSYS2 header snapshot. Every workaround here is listed in the
 * report so it can be confirmed against the actual CI toolchain.
 */

/* PROC_THREAD_ATTRIBUTE_SECURITY_CAPABILITIES / _HANDLE_LIST are usually
 * present in <processthreadsapi.h>, but older MinGW-w64 lacked them. They are
 * fixed constants (ProcThreadAttributeValue-encoded), so define if missing. */
#ifndef PROC_THREAD_ATTRIBUTE_SECURITY_CAPABILITIES
/* ProcThreadAttributeSecurityCapabilities == 9; Input|Additive, no Thread. */
#define PROC_THREAD_ATTRIBUTE_SECURITY_CAPABILITIES \
	ProcThreadAttributeValue(9, FALSE, TRUE, FALSE)
#endif
#ifndef PROC_THREAD_ATTRIBUTE_HANDLE_LIST
/* ProcThreadAttributeHandleList == 2; Input, no Thread/Additive. */
#define PROC_THREAD_ATTRIBUTE_HANDLE_LIST \
	ProcThreadAttributeValue(2, FALSE, TRUE, FALSE)
#endif

/* SECURITY_CAPABILITIES and SID_AND_ATTRIBUTES are declared in <winnt.h> when
 * _WIN32_WINNT >= 0x0602 (bumped above). CreateAppContainerProfile /
 * DeriveAppContainerSidFromAppContainerName are prototyped in <userenv.h> on
 * MinGW-w64 >= ~8.x. If the CI toolchain is older and the compile fails with
 * "unknown type SECURITY_CAPABILITIES" or an implicit-declaration warning for
 * those functions, a local fallback declaration must be added here (we could
 * not feature-test it without a Windows toolchain in this environment). This
 * is called out in the report as an item to confirm on CI. */

/* ---- UTF-16 command line assembly ---------------------------------------
 * Build the full UTF-16 command line from a UTF-8 argv[], quoting each element
 * per the CommandLineToArgvW rules: an element is wrapped in double quotes only
 * when it is empty or contains whitespace or a quote; a run of backslashes is
 * doubled only when it immediately precedes a double quote (or the closing
 * quote), and a literal double quote is escaped as \". This is the canonical
 * MS "Everyone quotes command line arguments the wrong way" algorithm.
 *
 * Returns a newly-allocated, NUL-terminated wchar_t* built via a GArray; free
 * it with g_free() (it is the detached g_array data pointer). Returns NULL on
 * a UTF-8 -> UTF-16 conversion failure. */
static wchar_t *sandbox_win_build_cmdline(gchar **argv) {
	GArray *buf = g_array_new(TRUE, FALSE, sizeof(wchar_t)); /* zero-terminated */

	for (gchar **a = argv; a && *a; a++) {
		if (a != argv) {
			wchar_t sp = L' ';
			g_array_append_val(buf, sp);
		}

		glong wlen = 0;
		gunichar2 *warg = g_utf8_to_utf16(*a, -1, NULL, &wlen, NULL);
		if (!warg) {
			g_array_free(buf, TRUE);
			return NULL;
		}

		/* Decide whether quoting is needed. */
		gboolean need_quotes = (warg[0] == L'\0');
		for (glong i = 0; !need_quotes && i < wlen; i++) {
			if (warg[i] == L' ' || warg[i] == L'\t' || warg[i] == L'"')
				need_quotes = TRUE;
		}

		if (!need_quotes) {
			for (glong i = 0; i < wlen; i++) {
				wchar_t c = (wchar_t) warg[i];
				g_array_append_val(buf, c);
			}
		} else {
			wchar_t q = L'"';
			g_array_append_val(buf, q);
			for (glong i = 0; i < wlen; i++) {
				/* Count a run of backslashes. */
				glong nbs = 0;
				while (i < wlen && warg[i] == L'\\') { nbs++; i++; }
				if (i == wlen) {
					/* Escape all backslashes so the closing quote is literal. */
					for (glong k = 0; k < nbs * 2; k++) {
						wchar_t bs = L'\\';
						g_array_append_val(buf, bs);
					}
					break;
				} else if (warg[i] == L'"') {
					/* Escape the backslashes AND the quote. */
					for (glong k = 0; k < nbs * 2 + 1; k++) {
						wchar_t bs = L'\\';
						g_array_append_val(buf, bs);
					}
					wchar_t dq = L'"';
					g_array_append_val(buf, dq);
				} else {
					/* Backslashes not before a quote stay literal. */
					for (glong k = 0; k < nbs; k++) {
						wchar_t bs = L'\\';
						g_array_append_val(buf, bs);
					}
					wchar_t c = (wchar_t) warg[i];
					g_array_append_val(buf, c);
				}
			}
			g_array_append_val(buf, q);
		}
		g_free(warg);
	}

	/* g_array has a trailing zero element (TRUE above). Detach the buffer. */
	wchar_t *line = (wchar_t *) g_array_free(buf, FALSE);
	return line;
}

/* Build the UTF-16, double-NUL-terminated environment block from a UTF-8
 * envp[] (each "KEY=VALUE"). Returns NULL to mean "inherit parent env" when
 * envp is NULL. On conversion failure returns NULL and sets *failed. */
static wchar_t *sandbox_win_build_env_block(gchar **envp, gboolean *failed) {
	*failed = FALSE;
	if (!envp)
		return NULL;   /* NULL => CreateProcessW inherits the parent env. */

	GArray *buf = g_array_new(FALSE, FALSE, sizeof(wchar_t));
	for (gchar **e = envp; *e; e++) {
		glong wlen = 0;
		gunichar2 *w = g_utf8_to_utf16(*e, -1, NULL, &wlen, NULL);
		if (!w) {
			*failed = TRUE;
			g_array_free(buf, TRUE);
			return NULL;
		}
		for (glong i = 0; i < wlen; i++) {   /* wlen excludes the terminator */
			wchar_t c = (wchar_t) w[i];
			g_array_append_val(buf, c);
		}
		wchar_t z = L'\0';
		g_array_append_val(buf, z);          /* terminate this entry */
		g_free(w);
	}
	/* A completely empty block must still be "\0\0"; add the final NUL. */
	wchar_t z = L'\0';
	g_array_append_val(buf, z);
	return (wchar_t *) g_array_free(buf, FALSE);
}

/* Add an inheritable access-allowed ACE for `sid` granting `access` on `path`.
 * Best-effort: a failure is logged and ignored (a missing/ungrantable path
 * must not abort the whole spawn — the child will simply be denied that path).
 * Uses the ANSI Get/SetNamedSecurityInfo since `path` is a UTF-8 g-string; we
 * convert to the local codepage via g_win32 helpers is avoided by using the
 * wide variants with a UTF-16 path. */
static void grant_appcontainer(const char *path, DWORD access, PSID sid) {
	if (!path || !*path || !sid)
		return;

	wchar_t *wpath = (wchar_t *) g_utf8_to_utf16(path, -1, NULL, NULL, NULL);
	if (!wpath) {
		siril_log_debug("sandbox: grant path utf16 conversion failed\n");
		return;
	}

	PACL old_dacl = NULL, new_dacl = NULL;
	PSECURITY_DESCRIPTOR sd = NULL;
	DWORD rc;

	rc = GetNamedSecurityInfoW(wpath, SE_FILE_OBJECT,
	                           DACL_SECURITY_INFORMATION,
	                           NULL, NULL, &old_dacl, NULL, &sd);
	if (rc != ERROR_SUCCESS) {
		siril_log_debug("sandbox: GetNamedSecurityInfo(%s) failed: %lu\n",
		                path, rc);
		g_free(wpath);
		return;
	}

	EXPLICIT_ACCESSW ea;
	ZeroMemory(&ea, sizeof(ea));
	ea.grfAccessPermissions = access;
	ea.grfAccessMode        = GRANT_ACCESS;   /* merge with existing ACEs */
	ea.grfInheritance       = SUB_CONTAINERS_AND_OBJECTS_INHERIT;
	ea.Trustee.TrusteeForm  = TRUSTEE_IS_SID;
	ea.Trustee.TrusteeType  = TRUSTEE_IS_WELL_KNOWN_GROUP;
	ea.Trustee.ptstrName    = (LPWSTR) sid;

	rc = SetEntriesInAclW(1, &ea, old_dacl, &new_dacl);
	if (rc != ERROR_SUCCESS) {
		siril_log_debug("sandbox: SetEntriesInAcl(%s) failed: %lu\n", path, rc);
		if (sd) LocalFree(sd);
		g_free(wpath);
		return;
	}

	rc = SetNamedSecurityInfoW(wpath, SE_FILE_OBJECT,
	                           DACL_SECURITY_INFORMATION,
	                           NULL, NULL, new_dacl, NULL);
	if (rc != ERROR_SUCCESS)
		siril_log_debug("sandbox: SetNamedSecurityInfo(%s) failed: %lu\n",
		                path, rc);

	if (new_dacl) LocalFree(new_dacl);
	if (sd) LocalFree(sd);
	g_free(wpath);
}

/* Fall back to a plain unsandboxed g_spawn (fail-open policy, matches Linux). */
static gboolean sandbox_win_spawn_unsandboxed(const char *working_dir,
                                              gchar **argv, gchar **envp,
                                              GPid *child_pid, gint *stdout_fd,
                                              gint *stderr_fd, GError **error) {
	return g_spawn_async_with_pipes(working_dir, argv, envp,
		G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
		NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);
}

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const SirilSandboxPolicy *policy, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	/* Explicit, granted opt-out (glue scripts spawning opaque third-party
	 * software): no AppContainer at all — a plain g_spawn. NB the child is then
	 * a normal-user process; the IPC pipe DACL companion change grants ALL
	 * APPLICATION PACKAGES, so the same fail-open pipe-access caveat noted there
	 * applies (a token-user ACE on the pipe is needed for non-AppContainer
	 * clients). */
	if (policy->unsandboxed)
		return sandbox_win_spawn_unsandboxed(working_dir, argv, envp,
		                                     child_pid, stdout_fd, stderr_fd, error);

	/* ---- 1. Container SID ------------------------------------------------ */
	PSID container_sid = NULL;
	HRESULT hr = CreateAppContainerProfile(SANDBOX_APPCONTAINER_NAME,
	                                       SANDBOX_APPCONTAINER_NAME,
	                                       SANDBOX_APPCONTAINER_DESC,
	                                       NULL, 0, &container_sid);
	if (hr == HRESULT_FROM_WIN32(ERROR_ALREADY_EXISTS)) {
		/* Profile already exists (created on a previous run): re-derive. */
		hr = DeriveAppContainerSidFromAppContainerName(
			SANDBOX_APPCONTAINER_NAME, &container_sid);
	}
	if (FAILED(hr) || !container_sid) {
		/* Fail-open: cannot build the container SID (too old an OS, policy,
		 * ...). Warn once and run the script unsandboxed rather than block it. */
		siril_log_message(_("Script sandbox unavailable on this system "
		                    "(AppContainer SID could not be created); running "
		                    "Python without confinement.\n"));
		return sandbox_win_spawn_unsandboxed(working_dir, argv, envp,
		                                     child_pid, stdout_fd, stderr_fd,
		                                     error);
	}

	/* ---- 2. Filesystem grants ------------------------------------------- */
	/* Read/write on the writable roots. */
	const DWORD RW = GENERIC_READ | GENERIC_WRITE | GENERIC_EXECUTE;
	const DWORD RX = GENERIC_READ | GENERIC_EXECUTE;

	grant_appcontainer(policy->wd, RW, container_sid);
	grant_appcontainer(policy->venv_path, RW, container_sid);
	grant_appcontainer(siril_get_user_data_dir(), RW, container_sid);
	{
		gchar *cfg = g_build_filename(siril_get_config_dir(), PACKAGE, NULL);
		grant_appcontainer(cfg, RW, container_sid);
		g_free(cfg);
	}

	/* Per-script extra paths from the granted manifest: extra_write RW,
	 * extra_read RX. With the broad read grants below, Windows is read-open like
	 * Linux/macOS, so extra_read_paths now only matter for locations the broad
	 * grant cannot reach (e.g. network shares, protected-ACL system paths). */
	if (policy->extra_write_paths)
		for (const char * const *e = policy->extra_write_paths; *e; e++)
			grant_appcontainer(*e, RW, container_sid);
	if (policy->extra_read_paths)
		for (const char * const *e = policy->extra_read_paths; *e; e++)
			grant_appcontainer(*e, RX, container_sid);

	/* Interpreter tree: grant read/execute on the directory of argv[0] (the
	 * venv Scripts dir) AND its parent. NB the BASE python install (real
	 * stdlib for a uv-managed interpreter) is NOT known here — this grant set
	 * WILL likely need dev iteration once a real script runs and Windows
	 * reports which reads are denied. See the limitations block above. */
	if (argv && argv[0] && *argv[0]) {
		gchar *interp_dir = g_path_get_dirname(argv[0]);
		if (interp_dir && *interp_dir) {
			grant_appcontainer(interp_dir, RX, container_sid);
			gchar *interp_parent = g_path_get_dirname(interp_dir);
			if (interp_parent && *interp_parent)
				grant_appcontainer(interp_parent, RX, container_sid);
			g_free(interp_parent);
		}
		g_free(interp_dir);
	}

	/* Read-anywhere, to match the Linux/macOS model: reads are permitted
	 * broadly and exfiltration is prevented by the network cut, not by read
	 * confinement (so a script can open images outside the working dir).
	 * AppContainer is deny-by-default for the filesystem, so read/execute must
	 * be granted explicitly:
	 *   - the user profile (the user owns it, so the DACL edit succeeds); this
	 *     also covers Pictures/Documents/Desktop, so the known-folder
	 *     capabilities are unnecessary;
	 *   - every FIXED drive root, best-effort — works on user-owned data drives
	 *     (e.g. D:\); the C:\ root DACL edit typically fails without elevation
	 *     and is skipped, but the profile grant already covers C:\Users\<user>;
	 *   - removable media is handled by the removableStorage capability added at
	 *     the spawn step below (a dynamic mount cannot be pre-ACLed).
	 * Writes remain confined to the roots granted earlier. */
	grant_appcontainer(g_get_home_dir(), RX, container_sid);
	{
		DWORD drives = GetLogicalDrives();
		for (int i = 0; i < 26; i++) {
			if (!(drives & (1u << i)))
				continue;
			char root[4] = { (char) ('A' + i), ':', '\\', '\0' };
			if (GetDriveTypeA(root) == DRIVE_FIXED)
				grant_appcontainer(root, RX, container_sid);
		}
	}

	/* ---- 4. Spawn -------------------------------------------------------- */
	/* All locals used after any `goto cleanup` are declared up front (no
	 * initializer that a goto could skip) to stay clean under
	 * -Wjump-misses-init on MinGW GCC. */
	gboolean result = FALSE;
	HANDLE out_r = INVALID_HANDLE_VALUE, out_w = INVALID_HANDLE_VALUE;
	HANDLE err_r = INVALID_HANDLE_VALUE, err_w = INVALID_HANDLE_VALUE;
	LPPROC_THREAD_ATTRIBUTE_LIST attr_list = NULL;
	wchar_t *cmdline = NULL, *env_block = NULL, *wworkdir = NULL;
	SECURITY_CAPABILITIES sec_caps;
	SID_AND_ATTRIBUTES cap_attrs[2];   /* removableStorage [+ internetClient] */
	PSID inet_cap_sid = NULL;          /* internetClient, only if allow_network */
	PSID removable_cap_sid = NULL;     /* removableStorage, for read-anywhere    */
	STARTUPINFOEXW si;
	PROCESS_INFORMATION pi;
	HANDLE inherit_handles[2];
	SIZE_T attr_size = 0;
	DWORD create_flags = EXTENDED_STARTUPINFO_PRESENT | CREATE_UNICODE_ENVIRONMENT;
	gboolean env_failed = FALSE;
	ZeroMemory(&sec_caps, sizeof(sec_caps));
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));

	/* stdout/stderr pipes: write ends inheritable, read ends NOT. */
	SECURITY_ATTRIBUTES pipe_sa;
	pipe_sa.nLength = sizeof(pipe_sa);
	pipe_sa.lpSecurityDescriptor = NULL;
	pipe_sa.bInheritHandle = TRUE;

	if (!CreatePipe(&out_r, &out_w, &pipe_sa, 0) ||
	    !CreatePipe(&err_r, &err_w, &pipe_sa, 0)) {
		g_set_error(error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
		            "CreatePipe failed: %lu", GetLastError());
		goto cleanup;
	}
	/* The read ends must NOT be inherited by the child. */
	SetHandleInformation(out_r, HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(err_r, HANDLE_FLAG_INHERIT, 0);

	/* Command line + environment + working dir (UTF-16). */
	cmdline = sandbox_win_build_cmdline(argv);
	if (!cmdline) {
		g_set_error(error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
		            "failed to build UTF-16 command line");
		goto cleanup;
	}
	env_block = sandbox_win_build_env_block(envp, &env_failed);
	if (env_failed) {
		g_set_error(error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
		            "failed to build UTF-16 environment block");
		goto cleanup;
	}
	if (working_dir && *working_dir)
		wworkdir = (wchar_t *) g_utf8_to_utf16(working_dir, -1, NULL, NULL, NULL);

	/* STARTUPINFOEX with two proc-thread attributes:
	 *   - SECURITY_CAPABILITIES (the container SID, zero capabilities)
	 *   - HANDLE_LIST (exactly the two pipe write ends we allow to inherit) */
	si.StartupInfo.cb = sizeof(STARTUPINFOEXW);
	si.StartupInfo.dwFlags = STARTF_USESTDHANDLES;
	/* WARNING (untested, verify on Windows): with PROC_THREAD_ATTRIBUTE_HANDLE_LIST
	 * set below, ONLY handles in that list are inherited. hStdInput here is the
	 * parent's console stdin, which is (a) not in the list and (b) likely not
	 * marked inheritable — under STARTF_USESTDHANDLES that may make CreateProcessW
	 * fail (ERROR_INVALID_HANDLE) on every launch. If so, the fix is to open an
	 * inheritable handle to "NUL" for stdin and add it to inherit_handles[], or to
	 * add the std handles to the list. Scripts do not need real stdin. */
	si.StartupInfo.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
	si.StartupInfo.hStdOutput = out_w;
	si.StartupInfo.hStdError = err_w;

	InitializeProcThreadAttributeList(NULL, 2, 0, &attr_size);
	attr_list = (LPPROC_THREAD_ATTRIBUTE_LIST) g_malloc0(attr_size);
	if (!InitializeProcThreadAttributeList(attr_list, 2, 0, &attr_size)) {
		g_set_error(error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
		            "InitializeProcThreadAttributeList failed: %lu",
		            GetLastError());
		goto cleanup;
	}

	sec_caps.AppContainerSid = container_sid;
	sec_caps.Capabilities = cap_attrs;
	sec_caps.CapabilityCount = 0;
	/* removableStorage (well-known cap SID S-1-15-3-10): read access to removable
	 * media — the dynamic-mount counterpart of the fixed-drive read ACEs above.
	 * It is a FILE-access capability, not a network one, so it does not weaken
	 * the no-network default. */
	if (ConvertStringSidToSidW(L"S-1-15-3-10", &removable_cap_sid)) {
		cap_attrs[sec_caps.CapabilityCount].Sid = removable_cap_sid;
		cap_attrs[sec_caps.CapabilityCount].Attributes = SE_GROUP_ENABLED;
		sec_caps.CapabilityCount++;
	} else {
		siril_log_debug("sandbox: removableStorage capability SID failed: %lu\n",
		                GetLastError());
	}
	/* internetClient (S-1-15-3-1) only when the policy allows network. With NO
	 * network capability an AppContainer cannot use WinSock at all — that is
	 * exactly how the default (no-network) policy denies egress by construction. */
	if (policy->allow_network) {
		if (ConvertStringSidToSidW(L"S-1-15-3-1", &inet_cap_sid)) {
			cap_attrs[sec_caps.CapabilityCount].Sid = inet_cap_sid;
			cap_attrs[sec_caps.CapabilityCount].Attributes = SE_GROUP_ENABLED;
			sec_caps.CapabilityCount++;
		} else {
			siril_log_debug("sandbox: internetClient capability SID failed: %lu\n",
			                GetLastError());
		}
	}
	if (sec_caps.CapabilityCount == 0)
		sec_caps.Capabilities = NULL;   /* both derivations failed */
	if (!UpdateProcThreadAttribute(attr_list, 0,
	        PROC_THREAD_ATTRIBUTE_SECURITY_CAPABILITIES,
	        &sec_caps, sizeof(sec_caps), NULL, NULL)) {
		g_set_error(error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
		            "UpdateProcThreadAttribute(SECURITY_CAPABILITIES) failed: %lu",
		            GetLastError());
		goto cleanup;
	}

	inherit_handles[0] = out_w;
	inherit_handles[1] = err_w;
	if (!UpdateProcThreadAttribute(attr_list, 0,
	        PROC_THREAD_ATTRIBUTE_HANDLE_LIST,
	        inherit_handles, sizeof(inherit_handles), NULL, NULL)) {
		g_set_error(error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
		            "UpdateProcThreadAttribute(HANDLE_LIST) failed: %lu",
		            GetLastError());
		goto cleanup;
	}
	si.lpAttributeList = attr_list;

	if (!CreateProcessW(
	        NULL,             /* application name: taken from cmdline */
	        cmdline,          /* mutable UTF-16 command line */
	        NULL, NULL,       /* process/thread security */
	        TRUE,             /* inherit handles (restricted by HANDLE_LIST) */
	        create_flags,
	        env_block,        /* NULL => inherit parent env */
	        wworkdir,         /* current dir (NULL => inherit) */
	        &si.StartupInfo,
	        &pi)) {
		g_set_error(error, G_SPAWN_ERROR, G_SPAWN_ERROR_FAILED,
		            "CreateProcessW failed: %lu", GetLastError());
		goto cleanup;
	}

	/* ---- 5. Return: process HANDLE as child_pid, read ends as CRT fds ---- */
	CloseHandle(pi.hThread);
	/* The parent must close the write ends so it sees EOF when the child
	 * exits; ownership of the read ends transfers to the CRT fds below. */
	CloseHandle(out_w); out_w = INVALID_HANDLE_VALUE;
	CloseHandle(err_w); err_w = INVALID_HANDLE_VALUE;

	if (stdout_fd) {
		*stdout_fd = _open_osfhandle((intptr_t) out_r, _O_RDONLY);
		if (*stdout_fd >= 0)
			out_r = INVALID_HANDLE_VALUE;   /* fd owns it now */
	}
	if (stderr_fd) {
		*stderr_fd = _open_osfhandle((intptr_t) err_r, _O_RDONLY);
		if (*stderr_fd >= 0)
			err_r = INVALID_HANDLE_VALUE;   /* fd owns it now */
	}
	if (child_pid)
		*child_pid = (GPid) pi.hProcess;   /* HANDLE, per g_spawn on Windows */

	result = TRUE;

cleanup:
	if (attr_list) {
		DeleteProcThreadAttributeList(attr_list);
		g_free(attr_list);
	}
	if (out_w != INVALID_HANDLE_VALUE) CloseHandle(out_w);
	if (err_w != INVALID_HANDLE_VALUE) CloseHandle(err_w);
	if (!result) {
		if (out_r != INVALID_HANDLE_VALUE) CloseHandle(out_r);
		if (err_r != INVALID_HANDLE_VALUE) CloseHandle(err_r);
	}
	g_free(cmdline);
	g_free(env_block);
	g_free(wworkdir);
	if (inet_cap_sid)
		LocalFree(inet_cap_sid);        /* allocated by ConvertStringSidToSid */
	if (removable_cap_sid)
		LocalFree(removable_cap_sid);   /* allocated by ConvertStringSidToSid */
	if (container_sid)
		FreeSid(container_sid);
	return result;
}

#else /* other POSIX — no confinement available; spawn plainly. */

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const SirilSandboxPolicy *policy, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	(void) policy;
	return g_spawn_async_with_pipes(working_dir, argv, envp,
		G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
		NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);
}

#endif /* platform */
