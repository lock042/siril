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

#define SANDBOX_MAX_PATHS 16

/* Internal to the Linux implementation (no longer exposed in the header, which
 * now offers only siril_sandbox_spawn()). */
typedef struct _SirilSandbox SirilSandbox;
struct _SirilSandbox {
	int install_seccomp;               /* 1 if a seccomp filter applies     */
	int allow_network;                 /* select net-on vs net-off filter   */
#ifdef HAVE_LANDLOCK
	guint64 landlock_access;           /* handled write mask, 0 => disabled */
	char *paths[SANDBOX_MAX_PATHS];    /* g_strdup'd writable paths          */
	int n_paths;
#endif
};

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

static SirilSandbox *sandbox_prepare(const char *wd,
                                     const char *venv_path,
                                     gboolean allow_network) {
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
	sb->allow_network = allow_network ? 1 : 0;

#ifdef HAVE_LANDLOCK
	{
		int abi = (int) syscall(__NR_landlock_create_ruleset, NULL, 0,
		                        LANDLOCK_CREATE_RULESET_VERSION);
		if (abi >= 1) {
			guint64 access = SANDBOX_ACCESS_WRITE_BASE;
			const char *tmpdir;
			int n = 0;

			/* Mask rights down to what the detected ABI supports:
			 * requesting unsupported bits yields EINVAL. */
			if (abi >= 2)
				access |= LANDLOCK_ACCESS_FS_REFER;
			if (abi >= 3)
				access |= LANDLOCK_ACCESS_FS_TRUNCATE;
			sb->landlock_access = access;

			/* Assemble the writable path list (resolved here in the parent
			 * so the child needs no g_getenv). Writes are confined to the
			 * project working dir, the venv (python writes __pycache__/.pyc
			 * there), temp dirs and the standard write-safe char devices.
			 * Reads are NOT handled and so remain allowed everywhere. */
			if (wd && *wd)
				sb->paths[n++] = g_strdup(wd);
			if (venv_path && *venv_path)
				sb->paths[n++] = g_strdup(venv_path);

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
			sb->n_paths = n;   /* <= SANDBOX_MAX_PATHS */

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
	for (int i = 0; i < sb->n_paths; i++)
		g_free(sb->paths[i]);
#endif
	g_free(sb);
}

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const char *wd, const char *venv_path,
                             gboolean allow_network, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	/* Build the ruleset spec in the parent; the actual Landlock ruleset and
	 * seccomp filter are installed by sandbox_child_setup() post-fork. If
	 * preparation fails entirely (unsupported kernel) sb is NULL and we spawn
	 * without a child_setup (fail-open). */
	SirilSandbox *sb = sandbox_prepare(wd, venv_path, allow_network);
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

/* Append (allow file-write* (subpath "PATH")) if path is non-empty. */
static void sb_allow_write_subpath(GString *p, const char *path) {
	if (!path || !*path)
		return;
	g_string_append(p, "(allow file-write* (subpath ");
	sb_append_quoted(p, path);
	g_string_append(p, "))");
}

/* Build a self-contained Seatbelt profile string (passed via sandbox-exec -p).
 * Allow-by-default (so the interpreter reads its libs freely), then confine:
 * deny all writes except the granted roots, and deny IP network egress while
 * leaving AF_UNIX (Siril IPC) untouched. */
static gchar *build_seatbelt_profile(const char *wd, const char *venv_path,
                                     gboolean allow_network) {
	GString *p = g_string_new("(version 1)(allow default)");

	/* Write confinement: global deny, then re-allow the granted roots (SBPL
	 * is last-match-wins, so the allows must follow the deny). */
	g_string_append(p, "(deny file-write*)");
	sb_allow_write_subpath(p, wd);
	sb_allow_write_subpath(p, venv_path);
	sb_allow_write_subpath(p, siril_get_user_data_dir());
	{
		gchar *cfg = g_build_filename(siril_get_config_dir(), PACKAGE, NULL);
		sb_allow_write_subpath(p, cfg);
		g_free(cfg);
	}
	sb_allow_write_subpath(p, g_get_tmp_dir());
	sb_allow_write_subpath(p, "/dev/null");
	sb_allow_write_subpath(p, "/dev/random");
	sb_allow_write_subpath(p, "/dev/urandom");

	/* Network: deny IP egress (exfil). AF_UNIX stays allowed by (allow
	 * default), so the Siril IPC socket is unaffected. */
	if (!allow_network)
		g_string_append(p, "(deny network-outbound (remote ip \"*:*\"))");

	/* Deny reading other processes' state (ptrace/task_for_pid analogue). */
	g_string_append(p, "(deny process-info*)");

	return g_string_free(p, FALSE);
}

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const char *wd, const char *venv_path,
                             gboolean allow_network, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	/* If sandbox-exec is missing (should never happen on macOS), fall back to
	 * an unsandboxed spawn with a warning rather than failing outright. */
	if (!g_file_test(SANDBOX_EXEC_PATH, G_FILE_TEST_IS_EXECUTABLE)) {
		siril_log_message(_("Script sandbox unavailable (%s missing); running "
		                    "Python without confinement.\n"), SANDBOX_EXEC_PATH);
		return g_spawn_async_with_pipes(working_dir, argv, envp,
			G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
			NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);
	}

	gchar *profile = build_seatbelt_profile(wd, venv_path, allow_network);

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
/* Windows: AppContainer confinement — implemented separately.          */
/* Placeholder retained so the tree builds; replaced by the real        */
/* AppContainer CreateProcess in the Windows implementation step.       */
/* ==================================================================== */

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const char *wd, const char *venv_path,
                             gboolean allow_network, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	(void) wd; (void) venv_path; (void) allow_network;
	/* TODO(securescripts): AppContainer CreateProcess. Until then, spawn
	 * unsandboxed so Windows behaviour is unchanged. */
	return g_spawn_async_with_pipes(working_dir, argv, envp,
		G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
		NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);
}

#else /* other POSIX — no confinement available; spawn plainly. */

gboolean siril_sandbox_spawn(const char *working_dir, gchar **argv, gchar **envp,
                             const char *wd, const char *venv_path,
                             gboolean allow_network, GPid *child_pid,
                             gint *stdout_fd, gint *stderr_fd, GError **error) {
	(void) wd; (void) venv_path; (void) allow_network;
	return g_spawn_async_with_pipes(working_dir, argv, envp,
		G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
		NULL, NULL, child_pid, NULL, stdout_fd, stderr_fd, error);
}

#endif /* platform */
