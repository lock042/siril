// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SRC_IO_SIRIL_SANDBOX_H
#define SRC_IO_SIRIL_SANDBOX_H

#include <glib.h>

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
                             const char *wd,
                             const char *venv_path,
                             gboolean allow_network,
                             GPid *child_pid,
                             gint *stdout_fd,
                             gint *stderr_fd,
                             GError **error);

#endif /* SRC_IO_SIRIL_SANDBOX_H */
