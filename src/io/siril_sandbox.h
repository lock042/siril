// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SRC_IO_SIRIL_SANDBOX_H
#define SRC_IO_SIRIL_SANDBOX_H

#include <glib.h>

typedef struct _SirilSandbox SirilSandbox;   // opaque

// Build in the PARENT, before fork. `wd`, `venv_path` may be NULL (skipped).
// Returns NULL if sandboxing is unavailable/disabled — caller then spawns
// with no child_setup (fail-open; see spec §6 policy note).
SirilSandbox *siril_sandbox_prepare(const char *wd,
                                    const char *venv_path,
                                    gboolean allow_network);

// GSpawnChildSetupFunc. Runs in the child AFTER fork, BEFORE exec.
// MUST be async-signal-safe: only raw syscalls, no malloc/glib/logging.
void siril_sandbox_child_setup(gpointer user_data);   // user_data = SirilSandbox*

// Parent-side cleanup after spawn returns (closes the ruleset fd, frees struct).
void siril_sandbox_finish(SirilSandbox *sb);

#endif /* SRC_IO_SIRIL_SANDBOX_H */
