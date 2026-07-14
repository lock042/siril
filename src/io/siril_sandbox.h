// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SRC_IO_SIRIL_SANDBOX_H
#define SRC_IO_SIRIL_SANDBOX_H

#include <glib.h>

/*
 * Per-script sandbox policy. The caller builds one (typically from the script's
 * [tool.siril.permissions] manifest, gated by user consent) and passes it to
 * siril_sandbox_spawn(). A zeroed struct is the safe default: no network, and
 * writes confined to the built-in roots (working dir, venv, Siril's data/config
 * dirs) with no extra paths.
 *
 *   wd, venv_path        Writable roots always granted (may be NULL).
 *   allow_network        TRUE relaxes the network confinement (boolean: the
 *                        underlying mechanisms cannot filter per-host).
 *   extra_write_paths    NULL-terminated list of additional writable paths, or
 *                        NULL. Granted read+write on every platform.
 *   extra_read_paths     NULL-terminated list of additional read-only paths, or
 *                        NULL. Meaningful where reads are confined (Windows
 *                        AppContainer); a no-op on Linux/macOS, where reads are
 *                        already unrestricted.
 *   unsandboxed          TRUE runs the interpreter with NO confinement at all
 *                        (plain spawn, pre-sandbox behaviour). The escape hatch
 *                        for glue scripts that spawn opaque third-party software
 *                        (e.g. the RC-Astro tools) whose needs are unknown, so
 *                        any sandboxing risks breaking them. Overrides every
 *                        other field. Must only ever be set from an explicitly
 *                        granted `unsandboxed = true` manifest request, gated by
 *                        a sterner consent than the granular permissions.
 */
typedef struct {
	const char *wd;
	const char *venv_path;
	gboolean allow_network;
	const char * const *extra_write_paths;
	const char * const *extra_read_paths;
	gboolean unsandboxed;
} SirilSandboxPolicy;

/*
 * Spawn the Python script interpreter under the platform sandbox, confining
 * its "Category-B" power (filesystem destruction + network exfiltration) that
 * bypasses the sirilpy interface via raw CPython.
 *
 * This mirrors the subset of g_spawn_async_with_pipes() that
 * execute_python_script() relies on, so all platform divergence is hidden here:
 *   - Linux:   g_spawn + a child_setup that installs Landlock + seccomp.
 *   - macOS:   g_spawn with the argv wrapped by /usr/bin/sandbox-exec (Seatbelt).
 *   - Windows: a bespoke AppContainer CreateProcess (g_spawn cannot set the
 *              security capabilities an AppContainer needs).
 *   - other:   plain g_spawn, unsandboxed (with a one-line warning).
 *
 * `working_dir` may be NULL (inherit cwd). `argv`/`envp` follow g_spawn
 * conventions. `wd`/`venv_path` are the writable roots granted to the script
 * (may be NULL). `allow_network` relaxes the network confinement (future
 * manifest opt-in; currently always FALSE).
 *
 * On success returns TRUE and fills `child_pid` (a process HANDLE on Windows),
 * `stdout_fd` and `stderr_fd` (CRT file descriptors, as g_spawn returns). On
 * failure returns FALSE and sets `error`.
 */
gboolean siril_sandbox_spawn(const char *working_dir,
                             gchar **argv,
                             gchar **envp,
                             const SirilSandboxPolicy *policy,
                             GPid *child_pid,
                             gint *stdout_fd,
                             gint *stderr_fd,
                             GError **error);

#endif /* SRC_IO_SIRIL_SANDBOX_H */
