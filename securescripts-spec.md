# securescripts — Linux sandbox for the Python script run phase (implementation spec)

Branch: `securescripts`. This spec is for an implementing sub-agent. Leave all changes in
the working tree (do NOT commit). Build to validate ONLY in `/home/siril-dev/mppbuild`
— never run `meson setup`/`reconfigure` under `/workspace` (it breaks the user's build).

## 1. Goal & threat model

Constrain **Category-B** power of a Siril Python script: the ability to do damage *around*
the sirilpy interface using raw CPython (`open()`, `os.system`, `socket`, `subprocess`,
`ctypes`). Two concrete threats, low-probability / high-impact:

- **Destruction** — `rm -rf ~/*`, overwriting/truncating files outside the project.
- **Exfiltration** — read personal files (`~/.ssh`, `~/Documents`) and POST them out.

Defenses, both enforced by the **kernel** on the spawned child (a script cannot remove
them from inside the interpreter — `ctypes`/native code bypass any Python-level hook):

1. **Landlock** — filesystem access control: **confine writes** to the project working dir
   + venv + temp; allow reads broadly (see §4 for why read-confinement is deliberately out
   of scope for this cut).
2. **seccomp-bpf** — **block internet egress**: deny `socket(2)` for `AF_INET`/`AF_INET6`,
   allow `AF_UNIX` (Siril IPC is an AF_UNIX socket, so this does not touch it).

Scope of THIS cut: **Linux only**. macOS/Windows get no-op stubs so the build and runtime
are unchanged there. Default policy: **network off, writes confined**, applied
unconditionally to every script run (no manifest yet — that is a later layer).

## 2. Why the run-phase spawn only (phase split)

The venv is provisioned by separate `g_spawn_sync(..., env, ...)` calls to `uv`
(around lines 3206, 3266, 4178, 4286, 4834, 4998 in `siril_pythonmodule.c`). Those need
network + cache/venv writes and must stay **unsandboxed**. Only the final script-run spawn,
`g_spawn_async_with_pipes` at **`siril_pythonmodule.c:5770`**, receives the sandbox. So the
phase split requires no restructuring — just hook that one call site.

## 3. Files to add

### `src/io/siril_sandbox.h`
Public API, platform-neutral:

```c
#ifndef SRC_IO_SIRIL_SANDBOX_H
#define SRC_IO_SIRIL_SANDBOX_H
#include <glib.h>

typedef struct _SirilSandbox SirilSandbox;   // opaque

// Build in the PARENT, before fork. `wd`, `venv_path` may be NULL (skipped).
// Returns NULL if sandboxing is unavailable/disabled — caller then spawns
// with no child_setup (fail-open; see §6 policy note).
SirilSandbox *siril_sandbox_prepare(const char *wd,
                                    const char *venv_path,
                                    gboolean allow_network);

// GSpawnChildSetupFunc. Runs in the child AFTER fork, BEFORE exec.
// MUST be async-signal-safe: only raw syscalls, no malloc/glib/logging.
void siril_sandbox_child_setup(gpointer user_data);   // user_data = SirilSandbox*

// Parent-side cleanup after spawn returns (closes the ruleset fd, frees struct).
void siril_sandbox_finish(SirilSandbox *sb);

#endif
```

### `src/io/siril_sandbox.c`
- Non-Linux (`#ifndef __linux__` or `#if !defined(HAVE_LANDLOCK)`): `prepare` returns NULL,
  `child_setup` is a no-op, `finish` frees nothing. Done.
- Linux implementation per §4/§5.

`SirilSandbox` holds only POD that `child_setup` reads post-fork (async-signal-safe to read
parent memory): `int ruleset_fd; int allow_network;` plus the pre-built seccomp
`struct sock_fprog` (or a flag selecting a `static const` filter). No pointers that require
heap access in the child.

## 4. Landlock

> **CORRECTION (found in review, applied in code):** the ruleset must be built in the
> **CHILD** (`child_setup`), NOT the parent. GLib's `GSpawnChildSetupFunc` contract says the
> callback runs *after GLib has performed all its setup* — which includes closing every
> descriptor ≥3 it does not manage. A ruleset fd created in the parent is therefore already
> closed by the time `child_setup` runs, so `landlock_restrict_self(fd)` fails with `EBADF`
> and — because we `_exit(127)` on that — **every script fails to launch** on a
> Landlock-capable kernel. (We must NOT use `G_SPAWN_LEAVE_DESCRIPTORS_OPEN` to keep the fd,
> as that leaks Siril's own fds into the untrusted interpreter.) So: the parent only probes
> the ABI, computes the access mask, and assembles the writable path list (as `g_strdup`'d
> strings, resolving `$TMPDIR` etc. so the child needs no `g_getenv`); the child creates the
> ruleset, adds the rules, and calls `landlock_restrict_self`, using only async-signal-safe
> syscalls (`open`/`fstat`/`close`/`syscall`). seccomp is unaffected — its filter is a
> `static const` array inherited via COW, needing no fd, which is why seccomp-in-child_setup
> is the common, working pattern and Landlock's fd was the wrinkle.
>
> Validated end-to-end against the real `g_spawn_*` path (not a bare `fork()` harness, which
> bypasses GLib's fd-closing and hides this bug): writes confined to wd, `$HOME`/`/etc`
> denied, `AF_INET` blocked, `AF_UNIX` intact.
>
> Build note: `HAVE_LANDLOCK` only lands in `config.h` after a meson **reconfigure**; a bare
> `ninja` may report "no work to do" and silently build the Landlock code out. Force a
> reconfigure (or clean the object) when first wiring the probe.

Original (parent-build) plan below is superseded by the correction above but kept for the
access-mask / path-list details, which still apply.

- Detect ABI: `landlock_create_ruleset(NULL, 0, LANDLOCK_CREATE_RULESET_VERSION)`.
  If it returns `-1` with `ENOSYS`/`EOPNOTSUPP` (kernel <5.13 or Landlock off), Landlock is
  unavailable → see §6.
- Mask the requested access rights down to what the detected ABI supports (older ABI lacks
  `REFER`/`TRUNCATE`/`IOCTL_DEV`); requesting unsupported bits gives `EINVAL`.
- Use `<linux/landlock.h>` for structs; call syscalls via `syscall(__NR_landlock_*)`
  (glibc has no wrappers). Provide the three thin static wrappers.

**Rule set (deliberately confine WRITES, allow READS):**
- Handled write-ish rights: `WRITE_FILE | REMOVE_DIR | REMOVE_FILE | MAKE_CHAR |
  MAKE_DIR | MAKE_REG | MAKE_SOCK | MAKE_FIFO | MAKE_BLOCK | MAKE_SYM |
  TRUNCATE(if ABI≥3) | REFER(if ABI≥2)`.
- **Do NOT handle** the read rights (`READ_FILE`, `READ_DIR`, `EXECUTE`). Under Landlock,
  access rights that are *not handled* are always allowed — so the interpreter can freely
  read its stdlib, `/usr`, `/etc`, the uv-managed python, etc., with no allow-list needed.
  This is why only writes need path rules.
- Grant the handled write rights on these paths only (skip any that are NULL/missing):
  - `wd` (the project working dir — `com.wd`),
  - `venv_path` (python writes `__pycache__`/`.pyc` here; alternatively set
    `PYTHONDONTWRITEBYTECODE=1` in env and you may drop this — but granting is simpler),
  - `TMPDIR` if set, else `/tmp`, plus `/var/tmp`,
  - `/dev/null`, `/dev/zero`, `/dev/full` (needed for normal write() to /dev/null).
- `landlock_restrict_self` requires `PR_SET_NO_NEW_PRIVS` (set in child, §5).

**Why read-confinement is out of scope for this cut:** with egress cut (§5), a script can
read `~/.ssh` but has nowhere to send it. Read-confinement is a later refinement; note it in
a code comment so it isn't mistaken for an oversight.

## 5. seccomp-bpf (build filter in parent, install in child)

Build the `sock_filter` program as a **`static const`** array (or in the struct) in the
parent; the child just calls `prctl` + `seccomp`. In `siril_sandbox_child_setup`, in order:

1. `prctl(PR_SET_NO_NEW_PRIVS, 1, 0, 0, 0)` — required for both unprivileged seccomp and
   Landlock. Abort child (`_exit(127)`) on failure.
2. `landlock_restrict_self(ruleset_fd, 0)` if `ruleset_fd >= 0`. On failure `_exit(127)`.
3. If `!allow_network`: `seccomp(SECCOMP_SET_MODE_FILTER, 0, &prog)` (or the `prctl`
   equivalent). On failure `_exit(127)`.

**Filter logic** (minimal; default ALLOW):
- Load `arch` (`offsetof(struct seccomp_data, arch)`); if not the expected arch, `KILL`
  (guards against arg layout differences). Target `x86_64` and `aarch64` — both have a
  direct `socket(2)` syscall so filtering on `args[0]` (domain) is reliable.
- Load `nr`; if `nr == __NR_socket`, load `args[0]` (domain); if `domain == AF_INET (2)` or
  `AF_INET6 (10)` → `RET_ERRNO(EAFNOSUPPORT)`. Everything else (incl. `AF_UNIX (1)`,
  `AF_NETLINK`) → `ALLOW`.
- All other syscalls → `ALLOW`.
- Handle the two arch NR values (they differ) — emit per-arch, or compile for the build
  arch via `#if defined(__x86_64__) / __aarch64__` and `KILL` on any other arch build so we
  never ship a filter that silently no-ops.

Only `socket()` creation is blocked; with no way to obtain an inet fd, `connect()` etc. need
no filtering. `allow_network=TRUE` (future manifest opt-in) skips step 3 but KEEPS Landlock.

**Async-signal-safety:** `child_setup` must call ONLY `prctl`, `syscall`, `seccomp`,
`_exit`. No `g_*`, no `malloc`, no `siril_log_*`, no `errno`-message formatting. Do not use
`GSpawnChildSetupFunc`'s environment to log; failures surface as the child exiting 127.

## 6. Availability / fail policy

If `siril_sandbox_prepare` cannot build a ruleset (old kernel, Landlock disabled): return
NULL and have the caller spawn with `child_setup = NULL`. i.e. **fail-open with a one-line
`siril_log_message` warning** ("Script sandbox unavailable on this kernel; running Python
without filesystem/network confinement."). Rationale: not breaking users on older kernels
for an MVP. Put the open/closed choice behind a single `#define SANDBOX_FAIL_OPEN 1` so it
can be flipped to fail-closed later. seccomp can still be installed even if Landlock is
absent (they are independent) — attempt each independently; a NULL SirilSandbox means
neither could be prepared.

## 7. Integration at the call site (`siril_pythonmodule.c` ~5768–5786)

- `#include "io/siril_sandbox.h"` at the top.
- Just before the spawn, after `env`/`python_argv` are ready:
  ```c
  SirilSandbox *sb = siril_sandbox_prepare(com.wd, venv_path, FALSE /* allow_network */);
  ```
  Note `venv_path` is still in scope here (freed at line ~5827, after the spawn). `com.wd`
  is the working dir (`core/siril.h:862`).
- Change the spawn call args 5 and 6 from `NULL, NULL` to
  `sb ? siril_sandbox_child_setup : NULL, sb`.
- After `g_spawn_async_with_pipes` returns (success or failure), call
  `siril_sandbox_finish(sb);` (closes the parent's ruleset fd, frees the struct). Make sure
  it runs on BOTH the success and the `!success` cleanup path.
- `GSpawnChildSetupFunc` is ignored by GLib on Windows, and our stub is a no-op, so no
  platform guard is needed at the call site.

## 8. Build system

- Add `'io/siril_sandbox.c',` to `src/meson.build` right after `'io/siril_pythonmodule.c',`
  (line ~204).
- Add a Landlock header probe near the other `cc.has_header` checks:
  ```meson
  if cc.has_header('linux/landlock.h')
    config_h.set('HAVE_LANDLOCK', 1)   # match the project's existing config-define idiom
  endif
  ```
  Find how this project defines config macros (grep `config.h`/`configuration_data`/
  `add_project_arguments('-DHAVE_`) and follow the SAME pattern — do not invent a new one.
- The .c file compiles on all platforms; the Linux body is guarded by `HAVE_LANDLOCK`
  (Landlock) and `__linux__` (seccomp). seccomp needs `<linux/seccomp.h>`,
  `<linux/filter.h>`, `<linux/audit.h>`, `<sys/prctl.h>`, `<sys/syscall.h>`.

## 9. Acceptance checks (build + run in /home/siril-dev/mppbuild only)

1. Builds cleanly on Linux; the sandbox source compiles.
2. A script doing `open('/tmp/../etc/x','w')` or writing anywhere outside `com.wd`/venv/tmp
   raises `PermissionError`; writing a file *inside* the working dir succeeds.
3. `import socket; socket.create_connection(('1.1.1.1',80),2)` raises
   `OSError`/`PermissionError` (EAFNOSUPPORT at socket creation); a normal sirilpy script still
   connects to Siril and runs commands (AF_UNIX unaffected).
4. On a kernel without Landlock, Siril logs the fail-open warning and scripts still run.

Report: files changed, the exact call-site diff, build result, and any of the four checks
you could/couldn't run headless.

## 10. Hardening additions (round 2) — three cheap, high-value vectors

All three land ONLY in `src/io/siril_sandbox.c` (no call-site or meson change). Keep the
existing Landlock design untouched. Everything new that runs in `child_setup` stays
async-signal-safe (raw syscalls / setrlimit / prctl / _exit only).

### 10a. seccomp: deny process-inspection syscalls (a real sandbox ESCAPE otherwise)
Same-uid `ptrace` / `process_vm_readv` lets a script read another process's memory (steal
browser cookies, ssh-agent keys) and — worse — inject code into an UN-sandboxed sibling the
user owns, escaping our per-process confinement entirely. Deny, with `RET_ERRNO(EPERM)` (not
KILL — benign probes shouldn't be fatal): `ptrace`, `process_vm_readv`, `process_vm_writev`,
`process_madvise`, `kcmp`, `pidfd_getfd`. These denials apply in BOTH the network-off and
network-on filters.

### 10b. seccomp: tighten socket() to an AF_UNIX allow-list
Today the filter denies only `AF_INET`/`AF_INET6` and ALLOWS everything else (incl.
`AF_NETLINK`, `AF_PACKET`). Change the default (network-off) filter to allow ONLY
`AF_UNIX` (== 1, needed for Siril IPC) and `RET_ERRNO(EAFNOSUPPORT)` for every other family.
For the future network-on path (`allow_network == TRUE`), install a filter that keeps the
10a ptrace denials but does NOT restrict socket families (inet permitted).

Because seccomp now does more than networking, it must be installed **whenever the arch is
supported, regardless of `allow_network`** (previously it was skipped when network was
allowed). Restructure: store `allow_network` in the struct and select between two
`static const` filters in `child_setup`:
- `sandbox_seccomp_insns` (network off): ptrace denials + AF_UNIX-only.
- `sandbox_seccomp_insns_net` (network on): ptrace denials only.

**BPF assembly guidance (avoid the offset foot-gun):**
- Keep the instruction layout FIXED so `#ifdef`s never shift jump offsets. For syscall
  numbers that may be absent from build headers, define a local fallback that can never
  match a real nr, so the compare instruction is still emitted:
  `#ifdef __NR_process_madvise ... #else #define SBX_NR_PROCESS_MADVISE 0xffffffffu #endif`
  (do this for `process_madvise` and `pidfd_getfd` at least; `ptrace`/`process_vm_*`/`kcmp`
  are old enough to assume present, but guard them the same way for safety).
- Use a shared RET trailer (`RET_KILL`, `RET_DENY_EPERM`, `RET_ERRNO_EAFNOSUPPORT`,
  `RET_ALLOW`) and jump to it, rather than sprinkling RETs — makes the offsets reviewable.
- VALIDATE behaviourally (see below), which catches any offset mistake regardless.

### 10c. setrlimit: DoS hardening (non-fatal, generous)
In `child_setup`, before/after the seccomp install, apply resource limits. These are
hardening, NOT a security boundary — a failed `setrlimit` must be **non-fatal** (ignore the
return; do NOT `_exit`), because breaking a legit script is worse than the DoS risk.
- `RLIMIT_NPROC` — anti-fork-bomb. CAVEAT: this limit is the real-UID's TOTAL live process
  count, not just this subtree, so set it GENEROUS (`#define SANDBOX_MAX_PROCS 4096`) or the
  child may be unable to fork at all on a busy desktop. Comment this prominently.
- `RLIMIT_FSIZE` — cap a single file's size. CAVEAT: astro images/stacks are large and
  exceeding this raises SIGXFSZ (kills the writer), so set it VERY generous
  (`#define SANDBOX_MAX_FSIZE_BYTES (64ULL << 30)` = 64 GiB). It only stops the single-giant-
  file case; total-disk exhaustion needs cgroups/quota (out of scope — note it).
Both via `#define` so they are tunable. Do NOT add memory/CPU limits — they fight legitimate
image processing.

### 10d. Limitations comment (required)
Add a clearly-marked comment block near the top of `siril_sandbox.c` documenting what this
Landlock+seccomp approach CANNOT close, so nobody mistakes it for full containment:
- **AF_UNIX reach to privileged local daemons**: we must allow AF_UNIX for Siril IPC, and
  neither seccomp (can't see the connect path/address) nor Landlock (doesn't mediate unix
  connect) can stop a script connecting to `docker.sock` (→ root), the ssh-agent socket
  (→ the user's keys), or the D-Bus session/system bus. Closing this needs a mount/network
  namespace or an LSM (AppArmor) profile.
- **X11 keylogging/screenshot for GUI (PyQt) scripts**: `$DISPLAY` access = full session
  snooping; cannot be closed while giving the script an X11 GUI. Needs Wayland or a nested X
  server.
- **Trusted computing base**: Siril's own C IPC handlers (`CMD_*`, SHM) run UNsandboxed; a
  memory-safety bug there reachable from a script is a full escape. Hardening those is
  separate/ongoing.
- **Resource limits are generous-by-necessity** and don't contain total disk/CPU exhaustion.
- **No kernel-attack-surface reduction**: the filter is targeted denials, not a syscall
  allow-list, so a kernel LPE via an allowed syscall is not mitigated.

## 11. Windows AppContainer (`#elif defined(_WIN32)` body of siril_sandbox_spawn)

Context: `siril_sandbox_spawn()` is the single cross-platform entry point (§ the
header). Windows AppContainer cannot be applied through `g_spawn` (it needs
`SECURITY_CAPABILITIES` at `CreateProcess` time), so the Windows path owns
process creation and must reproduce g_spawn's contract: fill `*child_pid` (the
process HANDLE cast to GPid), `*stdout_fd`, `*stderr_fd` (CRT fds via
`_open_osfhandle`), return TRUE/FALSE + `g_set_error` on failure.

BUILD ENV: Windows CI builds with **MinGW-w64 / GCC** (MSYS2), NOT MSVC. So:
- `#define _WIN32_WINNT 0x0A00` before `<windows.h>` (AppContainer APIs need ≥8;
  capability-SID helpers want 10).
- Includes: `<windows.h>`, `<userenv.h>` (CreateAppContainerProfile /
  DeriveAppContainerSidFromAppContainerName), `<sddl.h>`, `<aclapi.h>`,
  `<io.h>`, `<fcntl.h>`.
- Link `-luserenv` (advapi32 is auto-linked). Add userenv to the Windows link
  deps in `src/meson.build` (find the existing `host_machine.system() ==
  'windows'` deps block and append; if none, add one).
- Some AppContainer symbols may be absent/partial in the MinGW-w64 headers —
  if so, declare the missing prototypes/constants locally (guarded) and NOTE
  each one in the report so we can confirm on CI.

Implementation sequence:
1. **Container SID**: use a FIXED container name (e.g. `L"org.siril.script_sandbox"`).
   `CreateAppContainerProfile(name, name, desc, NULL, 0, &sid)`; if it returns
   `HRESULT_FROM_WIN32(ERROR_ALREADY_EXISTS)`, call
   `DeriveAppContainerSidFromAppContainerName(name, &sid)`. Zero capabilities
   (no `internetClient` etc.) → network is denied by construction. Free with
   `FreeSid`/`RtlFreeSid` as appropriate.
2. **Filesystem grants**: AppContainer is deny-by-default for user files, so the
   child cannot even load its interpreter without ACEs. Add an inheritable
   access-allowed ACE for the container SID via GetNamedSecurityInfo →
   SetEntriesInAcl (EXPLICIT_ACCESS, `SUB_CONTAINERS_AND_OBJECTS_INHERIT`) →
   SetNamedSecurityInfo, in a small `grant_appcontainer(path, access)` helper:
   - GENERIC_READ|WRITE|EXECUTE (recursive) on: `wd`, `venv_path`,
     `siril_get_user_data_dir()`, `<config>/siril`.
   - GENERIC_READ|EXECUTE on the interpreter tree: the directory of `argv[0]`
     (venv Scripts) AND the base Python install. NB the base-python location is
     not known here; grant the dir of argv[0] and its parent, and CLEARLY log +
     comment that the read-grant set for the interpreter/stdlib will likely need
     dev iteration (this is the most test-dependent part). Windows is thus
     read-confined too (stricter than Linux — acceptable, document it).
   - These ACEs persist on disk (only granting this one container SID). MVP
     leaves them; note it.
3. **IPC pipe (companion change, `src/io/siril_pythonmodule.c` create_connection,
   Windows branch, CreateNamedPipe ~line 1967)**: an AppContainer child cannot
   open the pipe unless its SD grants access. Simplest robust MVP: build a
   SECURITY_ATTRIBUTES whose DACL grants the well-known **ALL APPLICATION
   PACKAGES** SID (`S-1-15-2-1`, via ConvertStringSidToSid) plus the normal
   owner, and pass it to CreateNamedPipe. Guard so behaviour is unchanged when
   sandboxing is off. Flag clearly — without this, IPC breaks under AppContainer.
4. **Spawn**: CreatePipe ×2 (stdout/stderr; write ends inheritable, read ends
   NOT); InitializeProcThreadAttributeList with TWO attrs —
   `PROC_THREAD_ATTRIBUTE_SECURITY_CAPABILITIES` (SECURITY_CAPABILITIES{sid,0,NULL})
   and `PROC_THREAD_ATTRIBUTE_HANDLE_LIST` (the two write handles). Build the
   UTF-16 command line from `argv` with correct quoting (g_utf8_to_utf16 per arg;
   quote args containing spaces/quotes per CommandLineToArgvW rules) and the
   UTF-16 environment block from `envp` (double-NUL terminated) — or pass NULL
   env to inherit if `envp` is NULL. CreateProcessW with
   `EXTENDED_STARTUPINFO_PRESENT | CREATE_UNICODE_ENVIRONMENT`,
   bInheritHandles=TRUE, lpCurrentDirectory from `working_dir`.
5. **Return**: `_open_osfhandle((intptr_t)read_end, _O_RDONLY)` → *stdout_fd /
   *stderr_fd; `*child_pid = (GPid) pi.hProcess`; CloseHandle the write ends and
   pi.hThread in the parent. On any failure: CloseHandle everything, g_set_error
   (G_SPAWN_ERROR/G_SPAWN_ERROR_FAILED with GetLastError), return FALSE.
6. **Fail-open**: if the AppContainer SID cannot be created/derived, log a
   one-line warning and fall back to a plain `g_spawn_async_with_pipes` (matches
   the Linux fail-open policy) rather than blocking scripts.

CANNOT be compiled or run in this environment (no MinGW/Windows). Deliver a
best-effort draft to compile on the MinGW CI; runtime correctness (ACLs, pipe SD,
interpreter read-grants) will be iterated by a dev on Windows. Do NOT claim it is
validated. Report: exact files/functions changed, the meson link change, any
MinGW headers/symbols you had to declare locally, and the command-line/env
UTF-16 conversion approach.

### 10e. Acceptance (build in /home/siril-dev/mppbuild, validate via REAL g_spawn)
A bare `fork()` harness is NOT acceptable for validation (it bypasses GLib's fd-closing and
would hide bugs). Reuse the real-`g_spawn_*` pattern. Confirm:
1. Builds clean (force a rebuild of the sandbox object; a bare ninja may say "no work to do").
2. `socket(AF_INET)` → EAFNOSUPPORT; `socket(AF_NETLINK)` → EAFNOSUPPORT now too;
   `socket(AF_UNIX)` → OK.
3. `ptrace(PTRACE_TRACEME)` (or `os.system("cat /proc/1/maps")` via process_vm) → EPERM;
   a Python `import ctypes; ctypes.CDLL(None).ptrace(...)` attempt fails with EPERM.
4. Landlock write-confinement + AF_UNIX + stdout capture still work (no regression).
Report the filter layout you assembled, the struct/flag change, and each check's outcome.

## 12. Per-script permissions — manifest + remembered consent

User-chosen model (2026-07-14): scripts declare requested permissions in an in-file
`[tool.siril.permissions]` manifest; Siril prompts once and remembers the grant keyed by
script hash. Network is a boolean (per-host is not enforceable by seccomp/Seatbelt/
AppContainer). Three phases:

### Phase A — policy plumbing — DONE (commit 4fa5711eb)
`siril_sandbox_spawn()` takes a `SirilSandboxPolicy { wd, venv_path, allow_network,
extra_write_paths[], extra_read_paths[] }`, honored on all platforms (Windows `allow_network`
now adds the internetClient capability). Call site builds a DEFAULT policy → no elevation yet.

### Phase B — manifest parse — DONE (commit 5c8a54841)
`parse_siril_permissions(source, workdir)` in siril_pythonmodule.c reads
`[tool.siril.permissions]` (`network` bool, `write`/`read` path arrays, tokens
`$WORKDIR/$TMPDIR/$HOME/$USERDATA/$CONFIG`) into a `siril_script_perms` REQUEST. Currently
parsed + logged in execute_python_script(); NOT yet applied.

### Phase C — consent + remembered grants — TODO

Flow in execute_python_script(), after Phase B parsing:
1. If the request does not elevate beyond defaults (no `network`, no extra paths) → run with
   the default policy, no consent needed.
2. Otherwise look up a remembered grant for this script's hash (see store below):
   - grant present & covers the request → build the elevated `SirilSandboxPolicy` and run,
     silently.
   - no grant (script not yet approved) → decide by mode:
     - **GUI**: show a consent dialog listing exactly what is requested (network? which
       paths, read vs write?). Approve → persist grant + apply; Deny → run with DEFAULT policy
       (the script may fail if it truly needs the perm — that is the script's problem, not a
       silent elevation).
     - **Headless / CLI, no `-permissive`**: **default-DENY** — run with the default policy and
       log that elevated permissions were requested but not granted (no interactive prompt is
       possible).

**`-permissive` (new `pyscript` option).** Flips the headless/CLI default from default-deny to
**default-ALLOW + remember**: it approves the script's requested permissions AND persists the
grant to the store, so subsequent runs (even without `-permissive`) are pre-approved. i.e. it
is the CLI way to "approve a script that has not been run before", equivalent to clicking
Approve in the GUI dialog. Semantics: treat `-permissive` as an explicit user approval —
grant the full requested set, persist it, no dialog. (It does NOT grant more than the manifest
requests; it only consents to what was asked.)

- Parsing: `pyscript` currently accepts a leading `-async`. Add `-permissive` as another
  optional leading flag (parse both in a small loop before the script name; order between the
  two flags need not matter). Thread a `gboolean permissive` through
  process_pyscript → pyscript_data → execute_python_script_wrapper → execute_python_script.
- Update `STR_PYSCRIPT` in command_def.h to document `-permissive` (command help lives there,
  NOT in GUI tooltips per project convention).

**Grant store.** `script_perms.json` under the project dir, mirroring script_venvs.json:
versioned, keyed by script hash → { network: bool, write: [paths], read: [paths], approved:
ISO-8601, script_path (informational) }. Reuse the existing script-hash identity
(compute_script_hash / pep723 hash) so a script's grant survives edits the same way its venv
does. A grant only counts if it COVERS the current request (superset); if the manifest later
asks for MORE, re-prompt (or, headless, re-deny unless `-permissive`).

**Safety invariants.** A remembered grant is applied silently, but is never broader than what
the script's manifest requests at run time (defense against a stored over-grant). `-permissive`
and a GUI Approve are the only ways to CREATE a grant. Deny/headless never persists anything.

### `unsandboxed = true` — the all-or-nothing escape hatch (plumbed in Phase A/B; consent = Phase C)

Glue scripts that spawn opaque third-party software (e.g. the RC-Astro tools) have unknowable
needs, so ANY confinement can break the child. `unsandboxed = true` in the manifest requests
NO sandbox at all (plain spawn). Already wired: `siril_script_perms.unsandboxed` (parse) →
`SirilSandboxPolicy.unsandboxed` (every platform's siril_sandbox_spawn early-returns a plain
g_spawn). It OVERRIDES the granular fields.

Phase C consent treatment — deliberately STERNER than granular permissions, because it grants
full user-level access:
- It is the highest tier of request; the consent decision is separate from (and stricter than)
  network/path grants. A grant for network/paths does NOT imply an unsandboxed grant.
- GUI: a distinct, prominent warning dialog ("This script has asked to run with NO sandbox — it
  will have the same full access to your files and network as any program you run. Only approve
  scripts you fully trust."), visually differentiated from the granular-permission dialog.
- Every unsandboxed RUN (even an already-approved one) emits a persistent one-line log
  ("Running <script> with NO sandbox") so it stays visible, not just at approval time.
- `-permissive` question to resolve: does plain `-permissive` approve an unsandboxed request, or
  should that require a more explicit opt-in (e.g. a separate `-permissive-unsandboxed`, or
  `-permissive` + an interactive confirm even in CLI)? Leaning: `-permissive` MAY approve it
  (it is already an explicit user action) but MUST emit the stern log every run; revisit if it
  feels too easy. Store the grant with an explicit `unsandboxed: true` field so it is auditable.
- Windows caveat: an unsandboxed child is a normal-user process, so the IPC-pipe single-ACE DACL
  (ALL APPLICATION PACKAGES only) denies it — the same token-user-ACE fix flagged in §11 is
  required before the unsandboxed (and fail-open) paths work on Windows.
