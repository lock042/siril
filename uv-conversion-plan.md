# Plan: convert the Siril Python interface to use `uv`

Branch: `uv`
Author: Adrian Knagg-Baugh
Date: 2026-05-24

## 0. Progress tracking

This document is the source of truth for the conversion. Each
discrete unit of work has a checkbox.

**Rule: every time an item is completed, the corresponding
checkbox in this document MUST be updated from `[ ]` to `[x]` in
the same commit that lands the work** (or in a follow-up commit
on the same branch if completion is mechanically separate from
the code change, e.g. a CI pipeline run). Do not mark an item
complete until it is merged into the `uv` branch and has passed
CI. Partial work stays unchecked.

Summary tables in §10 (rollout phases) and §13 (touch-points) are
the canonical roll-ups; in-section checkboxes mirror these and
are updated in lock-step.

## 1. Goals and non-goals

### Goals

1. Replace `python -m venv` and `python -m pip install …` calls in
   `src/io/siril_pythonmodule.c` and `python_module/sirilpy/` with
   equivalent `uv` invocations, so that Siril benefits from a real
   dependency resolver and a content-addressed wheel cache.
2. Introduce *per-script* venvs (one venv per Python script file)
   that share storage automatically via uv's hardlink cache. This
   structurally solves the "scripts have dependency clashes" problem
   that today silently corrupts the single shared venv.
3. Keep the existing `sirilpy.ensure_installed()` API and the
   `gpuhelper.py` `TorchHelper`/`ONNXHelper`/`JaxHelper` API
   working unchanged so script authors do not have to rewrite their
   scripts. The wins (resolver quality, dedup, speed) accrue
   transparently.
4. Make PEP 723 inline script metadata an opt-in alternative to
   `ensure_installed()` for the deterministic deps of new scripts.

### Non-goals

1. `uv` will not be added as a `dependencies =` entry in
   `python_module/pyproject.toml`. Bootstrapping uv inside the venv
   it is supposed to manage is the wrong layer.
2. We will not break existing scripts that call
   `sirilpy.ensure_installed()`. Script authors decide at their own
   pace whether to migrate to PEP 723.
3. We will not move hardware-conditional torch/onnx/jax selection
   into PEP 723 metadata. That logic stays in
   `sirilpy.gpuhelper`, because it must branch on probed hardware.
4. No replacement of `g_spawn_*` plumbing or the named-pipe/socket
   IPC. Scope is strictly venv management.

## 2. Architecture overview

### Today

```
~/.local/share/siril/
├── venv/                            # one venv, shared by every script
│   └── (sirilpy + every script's runtime install_required pkgs,
│        with the inevitable clashes)
└── .python_module/                  # marker dir tracking installed sirilpy version
```

* `check_or_create_venv()` invokes `python -m venv <path>`.
* `install_module_with_pip()` invokes `python -m pip install <path>`.
* Each script does `sirilpy.ensure_installed("foo")`, which shells
  out to `python -m pip install foo` *inside that shared venv*.
* Hardware-conditional torch/onnx/jax installs are driven by
  `gpuhelper.py` and also shell out to `python -m pip install …`
  with computed `--index-url` / `-f` flags.

### After

```
~/.local/share/siril/
├── uv-cache/                        # UV_CACHE_DIR — on the same volume as venvs
├── python/                          # uv-managed Python toolchain (optional, see §4.1)
├── venv/                            # base venv (sirilpy-only; used when per-script
│                                    #   selector returns base, e.g. no PEP 723 + the
│                                    #   feature gate is off)
├── script_venvs/
│   └── <script-hash>/               # one venv per known script file
│       └── (sirilpy + that script's declared / ensure_installed-added deps)
└── script_venvs.json                # script-hash → venv metadata, last-used time
```

NOTE: an earlier draft of this plan proposed `venvs/_base/` for the
base venv. The implementation kept the base at `venv/` (matching
the legacy layout) and put per-script venvs in a sibling
`script_venvs/` directory. The two-tree layout keeps the legacy
path stable for callers that hardcode `<data>/siril/venv` (e.g.
`get_python_magic_number`) and makes per-script venvs purely
additive.

* A single `uv` binary on `PATH` (shipped or system-provided, see §6)
  is the only new external dependency.
* `check_or_create_venv()` becomes `uv_create_or_check_venv()` which
  shells out to `uv venv` for the `_base` venv and uses
  `uv pip install` for sirilpy. Both fall back to the existing
  pip-based path if `uv` is absent (Stage-1 compatibility).
* On script launch, `execute_python_script()` determines whether the
  script declares PEP 723 metadata. If so, a per-script venv is
  located/created and seeded via `uv pip install --requirements`
  derived from the inline `dependencies = [...]` list. If not, the
  `_base` venv is used (existing behaviour).
* `sirilpy.ensure_installed()` is reimplemented to call
  `uv pip install --python sys.executable <pkg>` when `uv` is on
  `PATH`, otherwise falls back to `pip`. The same change applies to
  every install/uninstall call site in `gpuhelper.py`. Because each
  script now has its own venv, an `ensure_installed()` from one
  script no longer disturbs another.

## 3. On-disk layout, env vars, and identity

### 3.1 Locations

| What                | Path                                                |
| ------------------- | --------------------------------------------------- |
| uv binary           | bundled (Win/macOS), system pkg (Linux/Flatpak)     |
| uv cache            | `$XDG_DATA_HOME/siril/uv-cache`                     |
| uv-managed Python   | `$XDG_DATA_HOME/siril/python` (optional)            |
| Base venv           | `$XDG_DATA_HOME/siril/venv` (kept stable for legacy callers) |
| Per-script venv     | `$XDG_DATA_HOME/siril/script_venvs/<script-hash>`   |
| Venv ledger         | `$XDG_DATA_HOME/siril/script_venvs.json`            |

`UV_CACHE_DIR` must be set explicitly so the cache lives next to the
venvs on a single volume; on Windows the default cache lives under
`%LOCALAPPDATA%` which is frequently on C: while Siril's data dir
may be on another drive, breaking hardlink dedup (see §6.3).

### 3.2 Venv identity

`<script-hash>` = first 16 hex chars of `SHA-256(absolute_script_path)`.

This keeps two scripts with the same basename in different
directories from sharing a venv. The ledger `venvs.json` maps:

```json
{
  "version": 1,
  "entries": {
    "<script-hash>": {
      "script_path": "/home/.../GPU_Manager.py",
      "created": "2026-05-24T10:00:00Z",
      "last_used": "2026-05-24T18:32:00Z",
      "pep723_hash": "ab12…",      // hash of the dependencies block; rebuild if changes
      "sirilpy_version": "1.1.12"  // rebuild if mismatched after sirilpy upgrade
    }
  }
}
```

`rebuild_venv()` in `siril_pythonmodule.c` is split into four
operations (see §4.8) covering surgical, single-script, nuclear,
and routine-maintenance use cases respectively.

## 4. C-side changes

All in `src/io/siril_pythonmodule.c` and `src/io/siril_pythonmodule.h`.

### 4.1 uv discovery (`find_uv_executable`)

- [x] Implemented and merged

New static helper:

```c
static gchar *find_uv_executable(GError **error);
```

Resolution order:

1. `g_getenv("SIRIL_UV")` if set (debug/override).
2. `<bundle>/bin/uv[.exe]` on Windows and macOS (we ship uv).
3. `g_find_program_in_path("uv")` on Linux/Flatpak.
4. NULL → caller falls back to legacy pip path during Stage 1.

A small `validate_uv_version()` companion runs
`uv --version` and ensures it is ≥ a configured floor
(e.g. `0.5.0`). Used by `validate_venv_health()` so a stale bundled
uv triggers the same rebuild path as a broken venv.

### 4.2 Replace `check_or_create_venv` → `uv_create_or_check_base_venv`

- [x] Implemented and merged

Lines 2743–2930 today. The replacement:

1. Resolve uv path via `find_uv_executable()`.
2. If found: invoke
   `uv venv --python <sys_python_or_bundled> --seed <venvs/_base>`.
   `--seed` ensures pip is installed inside the venv so scripts
   that shell out to `pip` continue to work in transition.
3. Wire `UV_CACHE_DIR=<.../siril/uv-cache>` and
   `UV_LINK_MODE=hardlink` (Linux/macOS auto-pick `clone` when
   available — let uv decide unless something goes wrong) into the
   subprocess env.
4. If `uv` not found *and* `SIRIL_REQUIRE_UV` not set: fall through
   to existing `python -m venv` path. Log a one-line debug note
   indicating which backend was used.
5. The existing `validate_system_python()` block (siril_pythonmodule.c:2123) is kept for the
   non-uv fallback path only. With uv we instead rely on `uv venv`
   to surface a clearer error and on uv's ability to download a
   Python toolchain (see §4.4 below).

### 4.3 Replace `install_module_with_pip` → `uv_install_sirilpy`

- [x] Implemented and merged

Lines 2311–2576 today. The replacement:

1. Resolve uv path.
2. `uv pip install --python <venv_python> --no-cache-dir <module_path>`
   — but actually we want to *use* the cache, so just
   `uv pip install --python <venv_python> <module_path>`.
3. The current retry loop (3 attempts × 300 s timeout, with regex
   matching on `stderr`) is removed. `uv` already retries network
   failures internally and exits with distinguishable codes; we
   simply propagate its exit status.
4. The atomic-rename temp-dir dance (`%s.tmp.%d` + `g_rename`,
   lines 2410–2570) goes away. `uv pip install` writes to the venv
   atomically and supports concurrent invocations through its own
   file lock. The lock also closes the race that exists today when
   two Siril instances start simultaneously and both try to
   install sirilpy into the same venv.
5. Version-comparison (`get_installed_module_version()`,
   lines 1976–2058) stays. `uv pip show sirilpy` replaces
   `python -m pip show sirilpy`; output format is the same.

### 4.4 Optional: uv-managed Python

- [x] Env-var gate wired (`SIRIL_UV_MANAGED_PYTHON=1` → `uv venv --python 3.12` + `UV_PYTHON_INSTALL_DIR`)
- [ ] Preference added
- [ ] Default-on for next major release decided

`uv venv --python 3.12` lets uv fetch and install a Python
interpreter under `UV_PYTHON_INSTALL_DIR`. This sidesteps the
long-standing pain in `validate_system_python()` where a distro
ships `python3` without `ensurepip`/`venv` (Linux error message at
siril_pythonmodule.c:2177).

Recommendation: opt-in via a preference
`prefs.python.use_managed_python` (default off in 1.x, default on
for a future major). When enabled and uv is present, we call
`uv python install --quiet 3.12` then point `uv venv` at it. Adds
~30 MB to first-run cost but removes a class of distro packaging
bugs we cannot otherwise solve.

### 4.5 New: per-script venv selection in `execute_python_script`

- [x] Extended `execute_python_script` signature plumbed through
- [x] PEP 723 parser implemented
- [x] `select_venv_for_script` selector implemented
- [x] Per-script venv creation wired up

Replace the single-`venv_path` lookup at lines 3269–3320. The
existing signature is unsuitable because the script editor (see
§4.5.1 below) always passes a tempfile path that differs every
launch, so hashing the *runtime* path would create one orphan venv
per click of the editor's Execute button. We therefore split
runtime path from venv-identity path in the API:

```c
// in siril_pythonmodule.h — extended signature
void execute_python_script(
    gchar       *script_name,         // runtime path (may be a tempfile)
    gboolean     from_file,
    gboolean     sync,
    gchar      **argv_script,
    gboolean     is_temp_file,
    gboolean     from_cli,
    gboolean     debug_mode,
    const gchar *venv_identity_path,  // NEW: canonical path for venv lookup;
                                      // NULL → use buffer-content rules (§4.5.1)
    const gchar *pep723_source        // NEW: text to parse PEP 723 from;
                                      // NULL → read from script_name on disk
);
```

`pep723_source` is needed because in the editor we want to honour
the *buffer's* metadata block, which may differ from whatever is
saved on disk (unsaved edits).

Selection helper:

```c
static gchar *select_venv_for_script(
    const gchar *venv_identity_path, // NULL means "use buffer-derived identity"
    const gchar *pep723_block,       // already-extracted block, or NULL
    gboolean    *out_ephemeral,
    GError     **error);
```

Behaviour:

1. If both `venv_identity_path` and `pep723_block` are NULL, or
   `from_file` is false with no identity (`-c` snippet from the
   CLI), return `venvs/_base`. No per-script venv.
2. If `venv_identity_path` is non-NULL, compute
   `<script-hash> = SHA-256(absolute(venv_identity_path))[:16]`.
3. Otherwise (no identity path but PEP 723 block present —
   editor's unsaved-buffer case), compute
   `<script-hash> = "unsaved-" + SHA-256(canonical_deps_array)[:16]`.
   Two unsaved buffers with identical declared deps share one venv,
   so the editor doesn't accumulate per-click garbage.
4. If `venvs/<script-hash>/` exists and its `venvs.json` entry's
   `pep723_hash` matches the supplied block's hash (or both are
   absent) and `sirilpy_version` matches the installed one →
   return it.
5. Otherwise, create it: `uv venv venvs/<script-hash>`, then
   `uv pip install sirilpy` from the on-disk module path, then if
   the script has PEP 723 metadata, `uv pip install <each dep>`
   from the declared list. Update `venvs.json`.
6. On sirilpy version bump (detected by comparing
   `com.python_version` against the ledger entry), rebuild *only
   that venv*. This is the per-script analogue of
   `rebuild_venv()` and is now cheap because the cache supplies
   most files via hardlink.

The PEP 723 parser is ~40 lines of pure C string handling — a small
GMatchInfo / GRegex pass for the `# /// script … # ///` block plus
a minimal TOML-subset reader for the `dependencies = [...]` array.
We do not need to support the full PEP 723 spec; a strict subset
(comments + a string array) is enough. If parsing fails, log a
warning and fall back to `_base`.

#### 4.5.1 Script editor integration (`src/gui/python_gui.c`)

- [x] PEP 723 extraction from buffer text wired up
- [x] `current_file` path threaded into `execute_python_script`

The editor's run path is `on_action_file_execute` at
python_gui.c:1148. Today it:

1. Dumps the source buffer text to
   `g_file_open_tmp("siril-script-XXXXXX.py", …)` →
   `temp_filename` (a fresh randomised path every click).
2. Calls
   `execute_python_script(temp_filename, TRUE /*from_file*/,
   FALSE /*async*/, script_args, TRUE /*is_temp_file*/, …)`.
3. The tempfile is deleted by `python_process_cleanup()` after the
   child exits.

The editor also maintains a `static GFile *current_file`
(python_gui.c:43) that tracks whatever the buffer was loaded from
or last saved as; it is `NULL` for a fresh unsaved buffer.

Changes:

* Extract the PEP 723 block from the buffer text *before* writing
  the tempfile (we have the text in hand at python_gui.c:1156).
* If `current_file` is non-NULL, call
  `execute_python_script(temp_filename, TRUE, FALSE, script_args,
  TRUE, from_cli, python_debug,
  g_file_get_path(current_file),  /* venv_identity_path */
  text                            /* pep723_source */ );`
  → venv is keyed off the saved file's path, so successive
  Run-from-editor clicks reuse the same env and stay in sync with
  the same venv used when the script is later launched from the
  script menu.
* If `current_file` is NULL, pass `venv_identity_path = NULL` and
  `pep723_source = text`. The selector then either keys off the
  PEP 723 block hash (`unsaved-<hash>` venv, shared across buffers
  that declare the same deps) or falls back to `_base` when the
  buffer declares no deps. Either way the editor never creates a
  fresh venv per click.
* The temp file deletion in `python_process_cleanup` is unchanged
  — the venv lifecycle is now decoupled from the tempfile
  lifecycle.

Stage 0/1 transitional behaviour: in Stage 0 we wire the two new
parameters through but `select_venv_for_script` always returns
`_base`. The behaviour is identical to today. Stage 2 turns on the
real selector. This means the editor change is a single mechanical
plumbing patch in Stage 0/1 and zero further work in Stage 2.

#### 4.5.2 Other callers of `execute_python_script`

- [x] `src/core/command.c:14572` updated
- [x] `src/gui/script_menu.c:208` updated
- [x] `src/io/siril_pythonmodule.c:2957` (`execute_startup_scripts`) updated

`grep -n execute_python_script src/` shows three callers besides
the editor; each needs a one-line update to pass the new
parameters:

| Call site                                                    | `venv_identity_path`               | `pep723_source` |
| ------------------------------------------------------------ | ---------------------------------- | --------------- |
| `src/core/command.c:14572` (`pyscript` CLI command)          | `data->script_name` (canonical)    | NULL (read from file) |
| `src/gui/script_menu.c:208` (script menu / launcher)         | `script_file` (canonical)          | NULL            |
| `src/io/siril_pythonmodule.c:2957` (`execute_startup_scripts`) | `script_path` (canonical)        | NULL            |
| `src/gui/python_gui.c:1189` (editor — see §4.5.1)            | `current_file` path or NULL        | buffer text     |

### 4.6 Ledger management

- [x] Ledger load/lookup/record implemented (yyjson read+write to `script_venvs.json`, atomic via temp+rename)
- [ ] LRU eviction helper implemented (deferred — will land with §4.8 "Clean up unused script environments" preference button)
- [ ] Cross-instance advisory lock implemented (deferred — single-instance assumption holds today; revisit if it bites)

A small `venvs.json` reader/writer using GLib's JSON helpers
(`json-glib`). Functions:

* `ledger_load(GError **error)`
* `ledger_lookup(const gchar *script_hash)`
* `ledger_record(const gchar *script_hash, const gchar *script_path, …)`
* `ledger_evict_lru(gsize max_venvs)` — prune oldest entries when a
  configurable cap is hit (default 50). Used by a settings-level
  "Clean up unused script environments" button.

### 4.7 Cache directory env var

- [x] `UV_CACHE_DIR` / `UV_NO_PROGRESS` env vars set on uv subprocesses (centralised in `inject_uv_env()`; applied at every uv invocation and on the spawned python child too)
- [x] `UV_PYTHON_INSTALL_DIR` set when managed-Python preference is enabled

When spawning any uv subprocess and when spawning the Python script
itself, set:

```
UV_CACHE_DIR=<XDG_DATA_HOME>/siril/uv-cache
UV_PYTHON_INSTALL_DIR=<XDG_DATA_HOME>/siril/python   (if managed Python enabled)
UV_NO_PROGRESS=1                                     (cleaner siril log)
UV_LINK_MODE                                         (do NOT set; let uv pick the best)
```

This goes in the same place that `PYTHONUNBUFFERED`/`PYTHONUTF8` are
already set today (siril_pythonmodule.c:3315–3318).

### 4.8 Venv rebuild and maintenance operations

Today there is a single `rebuild_venv()` function — nuclear-only,
no concept of per-script venvs, always wipes everything. With
per-script venvs in place (§4.5) this conflates three distinct
user needs:

* recovery from corruption (rare, blast-radius minimised by being
  targeted),
* "I changed my system Python and nothing works" (single nuclear
  action),
* "disk is filling up with venvs for scripts I don't use" (routine
  maintenance).

We split into four operations. The first three are recovery; the
fourth is maintenance.

#### 4.8.1 `rebuild_base_venv()` — surgical, base venv only

- [x] Existing `rebuild_venv()` renamed and scoped down
- [x] Preferences confirm-dialog text updated
- [x] Per-script venvs and uv cache explicitly preserved

Scope:

* Delete `<data>/siril/venv/` and `<data>/siril/.python_module/`.
* Leave `<data>/siril/script_venvs/`, `<data>/siril/script_venvs.json`,
  and `<data>/siril/uv-cache/` untouched.
* Call `clear_uv_state()` so the next `initialize_python_venv()`
  re-probes uv (handles "user installed a newer uv on the side").
* Re-run `initialize_python_venv_in_thread()`.

Why this is the right default: rebuilding the base venv is cheap
(one venv + cached wheels) and addresses the common "something is
off with sirilpy" symptom. Per-script venvs survive untouched.

Equivalent to today's `rebuild_venv()` minus the implicit
assumption that there is only one venv to care about.

#### 4.8.2 `rebuild_script_venv(script_hash)` — single per-script venv

- [x] C-side: function in `siril_pythonmodule.c` that wipes
      `script_venvs/<hash>/` and removes the ledger entry under
      a single ledger lock + ledger save (`rebuild_script_venv_by_path`
      — takes path, hashes internally to match the selector)
- [x] GUI integration via `treeview_scripts` context menu (see below)
- [x] Right-click handler refactored from direct-dialog to popup
      menu so the existing "Select Revision" item and the new
      "Rebuild Python environment" item can coexist

Scope:

* Take the ledger mutex.
* Delete `<data>/siril/script_venvs/<hash>/`.
* Remove the ledger entry for that hash.
* Persist the ledger.
* Done. Next launch of the script triggers `ensure_per_script_venv()`
  which rebuilds from the (preserved) uv cache — sub-second for
  most scripts because all wheels are already cached.

GUI integration (`src/gui/uifiles/settings_window.ui` →
`treeview_scripts`, populated by `git_gui.c::fill_script_repo_tree_idle`):

The existing right-click handler at git_gui.c:378–520
(`on_treeview_scripts_button_press`) opens a "Select Revision"
dialog directly. Refactor that to open a `GtkMenu` with two items:

```
Right-click on a Python script row →
┌──────────────────────────────────────┐
│ Open at earlier revision…            │   (existing behaviour)
│ Rebuild Python environment           │   (new, §4.8.2)
└──────────────────────────────────────┘
```

The "Rebuild Python environment" item is sensitive only when:

* the row's `COLUMN_IS_PYTHON` is TRUE (skip SSF scripts), and
* a ledger entry exists for the script's hash (otherwise the
  per-script venv hasn't been built yet — nothing to rebuild).

The handler computes the script hash the same way
`compute_script_hash()` does (SHA-256 of canonicalised
`COLUMN_SCRIPTPATH`), calls `rebuild_script_venv()`, and shows a
confirmation toast or modal — the action is fast enough that no
progress dialog is needed.

Optional follow-up (deferred, but worth noting): a new
`COLUMN_VENV_STATUS` next to `COLUMN_STARTUP` that shows at a
glance whether a per-script venv exists, is up-to-date, or is
stale relative to the script's PEP 723 block. This adds visual
clutter; keep behind a preference if implemented.

#### 4.8.3 `rebuild_all_python_state()` — nuclear option

- [x] C-side function that wipes everything (with optional
      uv-cache-also flag)
- [x] Confirmation dialog rewritten with the warning text and
      "Also clear uv cache" checkbox (default off)

Scope (cache preserved by default):

* `kill_all_python_scripts()`.
* Delete `<data>/siril/venv/`, `<data>/siril/.python_module/`,
  `<data>/siril/script_venvs/`, `<data>/siril/script_venvs.json`.
* If the "Also clear uv cache" checkbox is ticked: delete
  `<data>/siril/uv-cache/` too. Default OFF for the reasons
  discussed below.
* `clear_uv_state()`; re-run `initialize_python_venv_in_thread()`.

Why preserve the cache by default:

* The common trigger for this action is "Python upgraded" — old
  wheels become orphans (keyed by Python version) but cost only
  disk, never break anything. `uv cache prune` removes them on
  demand if disk space matters.
* The cache exists *so that rebuilds are fast*. Wiping it converts
  a sub-second per-venv rebuild into a multi-GB re-download for
  every script that uses heavy deps (numpy / torch / CUDA libs).
* The "cache is itself corrupt" failure mode is rare and better
  served by a separate "Clean uv cache" action (or by deleting the
  directory by hand) than by entangling it with "rebuild venvs."

Confirmation dialog text along the lines of:

> This will delete the base Python environment AND **N**
> per-script environments. They will be rebuilt automatically the
> next time each script runs, using wheels already in the local
> cache.
>
> ☐ Also clear the uv wheel cache (slow; forces re-download of
>   every dependency on next use).
>
> Continue?

The current preference button at `src/gui/preferences.c:1182`
becomes the trigger for this operation, with text updated from
"Rebuild Python venv" to something like "Rebuild all Python
environments…".

#### 4.8.4 `prune_unused_script_venvs(max_age_days)` — routine maintenance

- [x] Ledger eviction pass implemented (orphans pruned regardless;
      LRU only when `max_age_days > 0`, so passing `-days=0` from
      the CLI prunes only orphans)
- [x] Preferences UI: "Clean up unused script environments" button
      with an editable age threshold (default 90 days)

##### Bonus: `prune_uv_cache()` — reclaim disk from wheel cache

- [x] C-side function that shells out to `uv cache prune`
- [x] `pyenv_maint prune_cache` CLI subcommand
- [x] Preferences UI: "Prune uv wheel cache" button in the same grid
- [x] Pairs with `prune_unused_script_venvs`: prune venvs first to
      release wheel references, then prune the cache to collect the
      now-orphaned wheels. Both buttons live next to each other in
      Preferences → Python Interface for that reason.

Scope:

* Walk the ledger.
* For each entry: delete its venv directory + drop its ledger
  entry IF either:
  - the `script_path` no longer exists on disk (orphan venv), OR
  - the `last_used` timestamp is older than the configured
    threshold (default 90 days).
* Persist the ledger.
* Cache untouched.

Why a separate operation from `rebuild_all`: this is "disk getting
full" not "something is broken." Different user intent, different
mental model, different blast radius.

Default age threshold: 90 days. Provide an editable
`GtkSpinButton` in the preference panel alongside the button.

Report the result inline in the prefs panel: "Removed N
environments, freed M MB."

#### 4.8.5 Cross-cutting: the existing `rebuild_venv()` symbol

`rebuild_venv()` is called from `src/gui/preferences.c:1182` and
declared in `src/io/siril_pythonmodule.h:291`. To keep the patch
series easy to review:

* Rename in place: `rebuild_venv` → `rebuild_base_venv`. Update
  callers. This is the minimal change in `siril_pythonmodule.c`'s
  existing logic — same body, scoped name.
* Add the three new functions (`rebuild_script_venv`,
  `rebuild_all_python_state`, `prune_unused_script_venvs`) as
  separate symbols.
* Wire the preferences button at preferences.c:1182 to
  `rebuild_all_python_state()` (with the new confirm-dialog UI),
  not to `rebuild_base_venv()`, because the historical user intent
  of clicking that button is "nuclear option" — preserving that
  user intent is more important than preserving the symbol name.

## 5. Python-side changes (`python_module/sirilpy/`)

### 5.0 Backwards-compatibility constraint (NEW)

This MR ships in Siril 1.5. Siril 1.4 (already released, still in
use) ships an older sirilpy. **All Python-side changes must be
additive: existing function signatures and behaviour are
preserved**, so scripts that work today against Siril 1.4 continue
to work unchanged against Siril 1.5.

Concretely:

* `python_module/pyproject.toml` version bumps from
  **1.1.12 → 1.1.20** (micro-version; the existing
  `sirilpy.utility.check_module_version()` reads this string).
* All new public symbols (`pip_show`, `pip_list`, `pip_uninstall`,
  `declare_dependencies`) are added without modifying existing
  ones.
* Internal rewires (e.g. routing `_install_package()` through uv
  when present) are invisible to script authors — public API
  unchanged.

Any script change in the siril-scripts repository (notably
`GPU_Manager.py`, §13) that *uses* a 1.1.20+ symbol MUST gate the
call so the same file still imports cleanly on Siril 1.4. The
idiomatic gates:

```python
# Version-string gate (preferred — matches existing sirilpy idiom)
import sirilpy as s
if s.utility.check_module_version(">=1.1.20"):
    s.utility.pip_show("torch")
else:
    subprocess.run([sys.executable, "-m", "pip", "show", "torch"], ...)

# Duck-typing gate (also fine; useful when adding many calls)
if hasattr(s.utility, "pip_show"):
    s.utility.pip_show("torch")
else:
    subprocess.run(...)
```

An unconditional call to a 1.1.20+ symbol from a script that is
also expected to run on Siril 1.4 would crash with `AttributeError`
on import. CI for the scripts repository should verify both code
paths.

### 5.1 `utility.py::_install_package` and `ensure_installed`

- [x] `_resolve_installer()` helper added (named `_resolve_uv_path` +
      `_build_install_command` / `_build_uninstall_command` pair, since
      flag names differ between pip and uv — see below)
- [x] `_install_package` switched over
- [x] `uninstall_package` switched over

**Compatibility:** invisible to scripts. The public `ensure_installed`
/ `check_module_version` / `needs_module_version` signatures and
return values are unchanged. The only observable difference on
Siril 1.5 is that the underlying subprocess prefers `uv pip install`
when `uv` is on `PATH` — same semantics, better resolver, faster.
Scripts targeting Siril 1.4 see no behaviour difference because
sirilpy 1.1.12 (and earlier) keeps its existing pip-only path.

Reimplement `_install_package()` (utility.py:433) to prefer uv:

```python
def _resolve_installer() -> list[str]:
    uv = shutil.which("uv") or os.environ.get("SIRIL_UV")
    if uv:
        return [uv, "pip", "install", "--python", sys.executable]
    return [sys.executable, "-m", "pip", "install"]
```

Then build `cmd = _resolve_installer() + [flags…, target]`. The
`from_url` / `index_url` arguments already used by
`gpuhelper.py::_get_onnxruntime_package` (gpuhelper.py:826) map
cleanly: `--index-url` and `--find-links` are both supported by
`uv pip install`.

`uninstall_package()` (utility.py:347) similarly switches to
`uv pip uninstall` when available.

Public API of `ensure_installed`, `check_module_version`,
`needs_module_version` is unchanged.

### 5.2 `gpuhelper.py`

- [x] `TorchHelper.install_torch` routed through resolver
- [x] `TorchHelper.uninstall_torch` routed through resolver
- [x] New `sirilpy.utility.pip_show / pip_list` helpers exposed
      (`pip_uninstall` dropped — redundant with the existing
      `uninstall_package`; `GPU_Manager.py` should use that)
- [x] `ONNXHelper.uninstall_onnxruntime` also rewired (plan §5.2 listed
      only torch sites; ONNX and Jax had analogous direct pip calls,
      caught in implementation)
- [x] `JaxHelper.uninstall_jax` also rewired
- [ ] (scripts repo, separate PR) `GPU_Manager.py` switched to use the new helpers **behind a `check_module_version(">=1.1.20")` gate** with the existing subprocess path as the fallback

**Compatibility:** the `TorchHelper` / `ONNXHelper` / `JaxHelper`
public surface is unchanged. The new `pip_*` helpers are
additive, exported from `sirilpy.utility`, available only on
sirilpy ≥ 1.1.20. Scripts that want to use them must follow the
§5.0 gating pattern. `GPU_Manager.py`'s existing `subprocess.run([
sys.executable, '-m', 'pip', …])` calls stay as the fallback path
so the same file still runs on Siril 1.4.

Three call sites that build `[sys.executable, '-m', 'pip', …]`
arrays directly:

* `TorchHelper.install_torch` (gpuhelper.py:1608) — `pip install` with `--index-url`.
* `TorchHelper.uninstall_torch` (gpuhelper.py:1974) — `pip list` then `pip uninstall -y`.
* The `GPU_Manager.py` worker subprocesses (siril-scripts/core/GPU_Manager.py:74, 111) — these run *outside* sirilpy and call `pip show` directly.

For the helpers: route through the same `_resolve_installer()` (or
its `uninstall`/`list`/`show` siblings). For `GPU_Manager.py`:
because that script lives in the scripts repo and is not part of
sirilpy, the cleanest fix is to expose three new tiny helpers in
sirilpy and have `GPU_Manager.py` call them:

```python
# in sirilpy.utility
def pip_show(package: str) -> dict | None: ...
def pip_list() -> list[tuple[str, str]]: ...
def pip_uninstall(package: str) -> None: ...
```

`GPU_Manager.py` then drops its own `subprocess.run([sys.executable,
'-m', 'pip', …])` lines and uses these instead. Behaviour identical
on pip-only installs; faster and resolver-clean on uv.

### 5.3 New helper: `sirilpy.declare_dependencies` (optional, for transition)

- [x] `declare_dependencies()` implemented (parses PEP 508 requirement
      strings, delegates to `ensure_installed` which already short-circuits
      on installed packages)
- [x] Exported from `sirilpy/__init__.py`

**Compatibility:** `declare_dependencies` is a new public symbol
added in sirilpy 1.1.20. Scripts that call it must guard with
`check_module_version(">=1.1.20")` (or `hasattr`), exactly as for
the `pip_*` helpers. Recommended new-script idiom (works on both
versions):

```python
import sirilpy as s
# Single canonical declaration: works as a no-op when sirilpy 1.1.20+
# is paired with a Siril 1.5 that has already provisioned the venv
# from the PEP 723 block; falls back to runtime ensure_installed on
# older Siril.
if s.utility.check_module_version(">=1.1.20"):
    s.declare_dependencies(["numpy>=1.20", "scipy", "opencv-python"])
else:
    s.ensure_installed("numpy", "scipy", "opencv-python")
```

A no-op runtime helper that scripts can call to keep working both
on old and new Siril releases:

```python
# At top of a new-style script
import sirilpy as s
s.declare_dependencies(["numpy>=1.20", "scipy", "opencv-python"])
```

Internally:

* On a Siril that already created a per-script venv from PEP 723
  metadata: this is a no-op (the deps are already there).
* On an older Siril without PEP 723 support: this is equivalent to
  `s.ensure_installed(*pkgs)`.

This lets the script author write *one* declaration that works
across Siril versions. The recommended idiom for new scripts is:

```python
# /// script
# dependencies = ["numpy>=1.20", "scipy", "opencv-python"]
# ///
import sirilpy as s
s.declare_dependencies(["numpy>=1.20", "scipy", "opencv-python"])  # back-compat shim
```

Optional; not required for the conversion to work.

### 5.4 `pyproject.toml`

- [x] `version = "1.1.12"` → `version = "1.1.20"` (micro bump; new
      surface area is purely additive — see §5.0)

No other change. Deliberately: do *not* add `uv` to `dependencies`.
See non-goal #1.

We could add a build-system marker `[tool.uv]` block to express
authoritative torch-backend pinning for sirilpy itself, but sirilpy
doesn't depend on torch so this is moot for the module proper.

## 6. Build-time dependency changes

### 6.1 Windows installer

- [ ] Fetch+verify step added to `build/windows/native-gitlab-ci/siril-build.sh`
- [ ] `uv.exe` line added to `build/windows/installer/siril64.iss`
- [ ] Checksum file checked in

* Add `uv.exe` (Astral standalone Windows x86_64 build) to the
  bundle layout under `python/bin/uv.exe` or `bin/uv.exe`. Approx
  35 MB.
* Source for the bundled binary: pin a specific uv release tag in
  `build/windows/native-gitlab-ci/siril-build.sh`. Fetch with
  `curl -L https://github.com/astral-sh/uv/releases/download/<tag>/uv-x86_64-pc-windows-msvc.zip`,
  verify checksum against a checked-in `uv.sha256`, unzip into the
  install prefix. Add to `build/windows/installer/siril64.iss`
  similarly to the existing `python/*` line at siril64.iss:55.

### 6.2 macOS bundle

- [ ] `uv-aarch64-apple-darwin` fetched and placed in bundle
- [ ] `uv-x86_64-apple-darwin` fetched and placed in bundle
- [ ] Code-signing/notarisation regression-tested

* Same idea: pre-fetch `uv-aarch64-apple-darwin.tar.gz` *and*
  `uv-x86_64-apple-darwin.tar.gz`, drop the right binary into the
  app bundle's `Contents/Resources/bin/uv`.
* No code-signing complications: the standalone binary is
  Astral-signed and notarised; we just need it inside our own
  signed bundle.

### 6.3 Linux native packages and Flatpak

- [—] deb `Recommends: uv` added (out-of-tree: deb packaging is maintained downstream by Debian/Ubuntu — not in this repo)
- [—] rpm `Recommends: uv` added (out-of-tree; downstream)
- [—] Arch dependency note added (out-of-tree; downstream)
- [x] Flatpak manifest module added (`build/flatpak/org.siril.Siril.json`: new `uv` module with x86_64 + aarch64 sources, sha256-verified)

* Native packages (deb/rpm/Arch): add a `Recommends: uv` (deb) or
  the equivalent. If uv is missing at runtime we silently fall back
  to pip. We do *not* hard-require uv on Linux because the package
  is not yet in every distro.
* Flatpak (`build/flatpak/org.siril.Siril.json`): add a new module
  that downloads the prebuilt `uv-x86_64-unknown-linux-gnu.tar.gz`
  for the runtime architecture and installs the binary into
  `/app/bin/uv`. Cheaper than rebuilding from source; uv has no
  shared-library deps. Pin the release tag in the manifest.

### 6.4 AppImage (`build/appimage/`)

- [x] `wget`+verify+install step added to `build/appimage/generate.sh` (reads version from `build/UV_VERSION`, verifies sha256, installs uv + uvx into `appdir/usr/bin/`)
- [ ] Bundled AppImage sanity-checked (`--appimage-extract` + `uv --version`) — needs to be done by the maintainer running a full CI build of the AppImage. Stand-alone simulation of the download+install step passed locally.

The AppImage is the canonical Linux distribution channel for users
who can't or don't want to install from a distro package. The
recipe at `build/appimage/generate.sh` currently bundles glibc,
gdk-pixbuf loaders, glib schemas, fontconfig, and SSL certs into
`appdir/`. uv must be added so that AppImage users get the same
behaviour as bundled Windows/macOS users; otherwise the AppImage
would silently fall back to the legacy pip path on every machine
and we'd lose the resolver and dedup benefits exactly where the
single shared venv breaks most often (heterogeneous user Pythons).

Add a step to `build/appimage/generate.sh`, after the existing
`apt_bundle` calls (~line 65) and before the `linuxdeployqt`
invocation:

```bash
# Bundle uv (standalone binary; no shared-library deps)
UV_VERSION="0.5.x"   # pinned, kept in sync with §6.1/§6.2
wget -c -nv "https://github.com/astral-sh/uv/releases/download/${UV_VERSION}/uv-x86_64-unknown-linux-gnu.tar.gz"
# Verify the published sha256 (checksum file checked into the repo)
sha256sum -c "${SRCDIR}/build/uv-${UV_VERSION}-x86_64-unknown-linux-gnu.sha256"
mkdir -p tmp-uv && tar -xzf "uv-x86_64-unknown-linux-gnu.tar.gz" -C tmp-uv --strip-components=1
install -m 0755 tmp-uv/uv usr/bin/uv
install -m 0755 tmp-uv/uvx usr/bin/uvx
rm -rf tmp-uv "uv-x86_64-unknown-linux-gnu.tar.gz"
```

No change to `build/appimage/AppRun` is required:
`AppRun:28` already prepends `${HERE}/usr/bin` to `PATH`, so
`g_find_program_in_path("uv")` inside `find_uv_executable()`
locates the bundled uv automatically. We deliberately keep uv
*outside* the `ld-linux-x86-64.so.2 --inhibit-cache --library-path
…` interpreter invocation at AppRun:58 — uv is a statically linked
Rust binary and does not need (and would be confused by) the
patched library path. The kernel runs it directly when Siril spawns
it.

Sanity check for the maintainer: after building the AppImage, run

```
./Siril-x86_64.AppImage --appimage-extract
ls squashfs-root/usr/bin/uv
squashfs-root/usr/bin/uv --version
```

to confirm the binary is present and runs.

### 6.5 Build-image Docker recipes (`build/build-image/`)

- [x] `Dockerfile.debian-oldstable` updated (ARG-based pin, sha256-verified, `uv --version` smoke at build time)
- [x] `Dockerfile.debian-latest` updated (same pattern)
- [—] `Dockerfile.win64-latest` updated (Windows out-of-scope — packaging done separately)
- [x] `build/UV_VERSION` file created (single line: pinned version string). Used by the AppImage script via `cat`; Dockerfile `ARG UV_VERSION` defaults track it, with CI able to override via `--build-arg UV_VERSION=$(cat build/UV_VERSION)`

Three Dockerfiles in this directory are used by the CI pipeline:

* `Dockerfile.debian-oldstable` — builds the AppImage and runs the
  CI integration tests. `python3 / python3-pip / python3-venv`
  are already installed (lines 42–44). Add uv so the test suite
  can exercise the new code paths and produce the bundled AppImage
  (see §6.4):
  ```dockerfile
  # After the existing `RUN pip install meson==1.9.1 ninja!=1.11.1` (line 122)
  RUN UV_VERSION=0.5.x \
      && wget -O /tmp/uv.tar.gz \
         https://github.com/astral-sh/uv/releases/download/${UV_VERSION}/uv-x86_64-unknown-linux-gnu.tar.gz \
      && tar -xzf /tmp/uv.tar.gz -C /tmp \
      && install -m 0755 /tmp/uv-x86_64-unknown-linux-gnu/uv /usr/local/bin/uv \
      && install -m 0755 /tmp/uv-x86_64-unknown-linux-gnu/uvx /usr/local/bin/uvx \
      && rm -rf /tmp/uv.tar.gz /tmp/uv-x86_64-unknown-linux-gnu
  ```
  This is the most important of the three: oldstable verifies that
  uv works against the oldest Python we still claim to support
  (currently Python 3.11 in Debian bookworm — uv supports anything
  ≥ 3.8 so this is comfortably fine).
* `Dockerfile.debian-latest` — symmetrical update; same `RUN`
  block, appended after the `pip install --break-system-packages
  meson … ninja …` step at line 52. Required so the test matrix
  covers both old and new Python/distro lines.
* `Dockerfile.win64-latest` — used for cross-building the Windows
  installer. uv is bundled into the *installer payload* (see §6.1)
  rather than executed inside this image, so no docker-level
  install is strictly required. However it is convenient to also
  install uv here so that fetch-and-verify steps in
  `build/windows/native-gitlab-ci/siril-build.sh` can use the
  *same* uv version as the bundle (consistency between the build
  agent and the artefact). Recommended:
  ```dockerfile
  RUN UV_VERSION=0.5.x \
      && wget -O /tmp/uv.tar.gz \
         https://github.com/astral-sh/uv/releases/download/${UV_VERSION}/uv-x86_64-unknown-linux-gnu.tar.gz \
      && tar -xzf /tmp/uv.tar.gz -C /tmp \
      && install -m 0755 /tmp/uv-x86_64-unknown-linux-gnu/uv /usr/local/bin/uv \
      && rm -rf /tmp/uv.tar.gz /tmp/uv-x86_64-unknown-linux-gnu
  ```

Single source of truth for the uv version: a top-level
`build/UV_VERSION` file containing the pinned tag, sourced by all
three Dockerfiles, the AppImage `generate.sh`, the macOS bundle
recipe (§6.2), and the Windows installer build (§6.1). One bump
updates every artefact.

### 6.6 Source builds

No change required. The build does not need uv at build time. uv is
strictly a *runtime* prerequisite, located via `find_uv_executable`
at first launch. Source builders without uv get the legacy
pip-based path.

### 6.7 Meson

`python_module/meson.build` is unchanged. We are not building uv;
we are bundling a binary. The bundling happens in the per-platform
packaging scripts.

## 7. Runtime dependency changes

| Component                        | Before              | After                                            |
| -------------------------------- | ------------------- | ------------------------------------------------ |
| System Python (≥ 3.9)            | required            | required (still — uv invokes it)                 |
| System `python3-venv`            | required            | not required when uv is present                  |
| System `python3-pip`             | required            | not required when uv is present (uv has its own resolver) |
| `uv` binary                      | n/a                 | required on Win/macOS bundles; recommended on Linux |
| Disk: per script (after dedup)   | varies; growing single venv | one cache copy + small per-venv overhead    |
| First-run network                | several pip rounds  | one uv resolve + download (≈ 5–10× faster)       |

## 8. User-facing changes

1. **First launch is faster.** uv installs sirilpy in seconds
   rather than tens of seconds. The "Installing / updating python
   module in the background" message stays; it just resolves
   sooner.
2. **Scripts no longer interfere with each other.** A script that
   needs `torch+cu128` and a script that needs CPU-only `torch` can
   both work without the second one breaking the first. This is the
   main user-visible reliability improvement.
3. **Disk usage shifts**: each script gets its own venv, but the
   shared cache means the big wheels (torch + CUDA libs) live on
   disk once across N scripts (provided the cache is on the same
   volume — handled by `UV_CACHE_DIR` pinning).
4. **Preferences gains a "Clean up unused script environments"
   button.** Prunes per-script venvs whose scripts no longer exist
   on disk, or whose `last_used` is older than a configurable cut
   (default 90 days).
5. **The "Rebuild Python venv" preference now rebuilds *all* envs.**
   The confirmation dialog text needs updating to match.
6. **No user-facing API breakage.** Existing scripts continue to
   run unchanged.

## 9. Script-writer-facing changes

### Backwards-compatible (Stage 1)

* `sirilpy.ensure_installed()`, `check_module_version()`,
  `needs_module_version()` continue to work identically. Internally
  they now use uv when available; the script does not need to know.
* `sirilpy.gpuhelper.{TorchHelper,ONNXHelper,JaxHelper}` continue
  to work identically. Their install methods are faster and their
  resolver is stricter (an install that would have silently broken
  another script's venv now either succeeds cleanly or fails loudly
  in the script that requested it — the failure is now isolated).
* New optional helpers: `sirilpy.utility.pip_show()`,
  `pip_list()`, `pip_uninstall()` for scripts (like
  `GPU_Manager.py`) that currently spawn `pip` themselves.

### New optional idiom (Stage 2)

```python
# /// script
# dependencies = [
#   "numpy>=1.20",
#   "scipy",
#   "opencv-python",
# ]
# ///

import sirilpy as s
# … no ensure_installed needed for the static deps above …

# Hardware-conditional installs still go through gpuhelper:
from sirilpy.gpuhelper import TorchHelper
TorchHelper().install_torch()      # picks the right backend
```

* The PEP 723 block is parsed by Siril before the script runs.
  Siril seeds the per-script venv with those deps so they are
  installed by the time the script imports them.
* The `gpuhelper` API still owns hardware-conditional installs;
  these cannot be expressed statically in PEP 723.
* Script repository maintainers can adopt PEP 723 gradually;
  scripts using only `ensure_installed()` continue to work.

### Things to *avoid* in new scripts

* Don't put hardware-dependent packages (`torch`, `onnxruntime-*`,
  `jax-*`) in the PEP 723 block. Use `gpuhelper`.
* Don't shell out to `pip` directly. Use `sirilpy.utility.pip_*`
  helpers (added in this work) so uv handles the resolve when
  present.

## 10. Migration / rollout

| Done  | Phase | Scope                                                                                          | Risk   |
| ----- | ----- | ---------------------------------------------------------------------------------------------- | ------ |
| `[x]` | **0** | Wrap uv discovery behind a feature gate (`SIRIL_USE_UV=1`); ship uv in Win/macOS bundles only. | Low    |
| `[x]` | **1** | Default-on uv for venv/install when present. Gate inverted to opt-out (`SIRIL_USE_UV=0`).      | Medium |
| `[x]` | **2** | Default-on per-script venvs with PEP 723 parsing. Opt-out via `SIRIL_PER_SCRIPT_VENVS=0`.      | Medium |
| `[~]` | **3** | Migrate scripts repo: add PEP 723 to the most-installed scripts. Update docs. (GPU_Manager.py done in siril-scripts branch `uv-pip-helpers` — commits `9d580a1` for pip-helper routing, `2088663` for the GPU prefs config. Broader PEP 723 rollout to other scripts deferred to platform maintainers.) | Low    |
| `[ ]` | **4** | (Optional) Default-on uv-managed Python.                                                       | Low    |

Each phase is independently shippable. Phase 0 can land on master
immediately; users opt in with the env var to flush out bundle
problems. Phase 1 requires the bundled binary to be solid and the
ledger code to be well-tested. Phase 2 requires the PEP 723 parser.

The `uv` branch we are on should target Phase 0 + Phase 1 in this
PR series, with Phase 2 as a clear follow-up.

## 11. Testing plan

### 11.1 Unit-level (C)

- [ ] `find_uv_executable()` returns NULL when absent, bundled path when present, env override when set
- [ ] `validate_uv_version()` rejects too-old uv
- [ ] PEP 723 parser golden tests against real scripts in `siril-scripts/`
- [ ] Ledger: round-trip JSON, eviction, hash-based lookup

### 11.2 Integration

- [ ] Test script added to `src/tests/`
- [ ] Clean-init test: `venvs/_base` created and contains `sirilpy`
- [ ] Two-script clash test: each gets its own venv, both succeed
- [ ] Dedup test: `du` confirms single inode across venvs on Linux
- [ ] Cleanup test: deleting a script file → its venv is evicted

### 11.3 Manual cross-platform

- [ ] Linux: deb-installed `uv`
- [ ] Linux: no `uv` (legacy pip fallback)
- [ ] Linux: AppImage with bundled `uv`
- [ ] Linux: Flatpak with bundled `uv`
- [ ] macOS: bundled signed/notarised `uv` from inside the .app
- [ ] Windows: bundled `uv.exe`, cache on same drive as venvs
- [ ] Windows: bundled `uv.exe`, cache forced to a different drive (verify warning + copy fallback)

### 11.4 Regression baseline

- [ ] `AutoBGE.py` timed on pip + uv (first/second run)
- [ ] `AberrationRemover.py` timed on pip + uv
- [ ] `GPU_Manager.py` timed on pip + uv
- [ ] `onnxruntime-tester.py` timed on pip + uv
- [ ] `torch-tester.py` timed on pip + uv
- [ ] Follow-up benchmark note published

Pick five representative scripts from `siril-scripts/`:
`AutoBGE.py` (PyQt6 + scipy), `AberrationRemover.py` (PyQt6 +
scipy), `GPU_Manager.py` (PyQt6 + gpuhelper machinery),
`onnxruntime-tester.py`, `torch-tester.py`. Time end-to-end first
run, second run, both on pip and uv backends.

## 12. Open questions

1. **Should `_base` exist at all?** Alternative: always create a
   per-script venv on first run. Simpler ledger semantics but
   slightly slower first-launch for the very simplest scripts.
   Recommendation: keep `_base` for now; revisit when the ledger
   code lands.
2. **Concurrency between multiple Siril instances.** uv already
   takes a file lock around its cache. We should also take a
   ledger-level advisory lock when mutating `venvs.json`. GLib's
   `GFileIOStream` plus `g_file_replace_contents` provides
   atomic-rename semantics; combine with a `.lock` sentinel.
3. **What is the uv version floor?** Pin to the latest tagged
   release at the time of merge (likely `0.5.x`). PEP 723 has been
   supported since `uv 0.4.x`; `--seed` since `0.3.x`. Anything
   newer than 0.5 should be safe.
4. **Telemetry on legacy script breakage.** If we add a one-shot
   diagnostic that reports the first uv command that failed (without
   exfiltrating script content), we can catch real-world resolver
   issues sooner. Out of scope here; flag for a separate proposal.
5. **`gpuhelper.JaxHelper` is hidden today (per `GPU_Manager.py`
   1.0.2 changelog).** Should we re-enable Jax-via-uv once the
   resolver-induced conflicts with torch's strict CUDnn pin are
   isolated to a per-script venv? Probably yes once Phase 2 lands.

## 13. Summary of code touch-points

This table is the canonical work-item roll-up. Tick each row as
the change merges into the `uv` branch.

| Done  | File                                            | Change                                                |
| ----- | ----------------------------------------------- | ----------------------------------------------------- |
| `[ ]` | `src/io/siril_pythonmodule.c`                   | uv discovery, venv create/install rewrites, PEP 723 parser, ledger, per-script venv selection |
| `[ ]` | `src/io/siril_pythonmodule.h`                   | extended `execute_python_script` signature (adds `venv_identity_path`, `pep723_source`); new ledger/cleanup prototypes |
| `[ ]` | `src/gui/python_gui.c:1148–1192`                | pass `current_file` path + buffer text into the extended `execute_python_script`; PEP 723 extraction from buffer |
| `[ ]` | `src/gui/script_menu.c:208`                     | pass canonical `script_file` for venv identity (one-line update) |
| `[ ]` | `src/core/command.c:14572`                      | pass `data->script_name` for venv identity (one-line update) |
| `[x]` | `src/gui/preferences.c:1175+`                   | rewired existing button to `rebuild_all_python_state()` with cache checkbox; added handlers for "Rebuild base" and "Prune" |
| `[x]` | `src/gui/git_gui.c:378` (`on_treeview_scripts_button_press`) | refactored direct-dialog right-click to popup menu; "Rebuild Python environment" item added (greyed for non-Python rows) |
| `[x]` | `src/gui/uifiles/settings_window.ui`            | expanded "Python Interface" grid from 1×2 to 2×2: added Rebuild base, Rebuild all, Prune buttons + spin button for max-age |
| `[x]` | `src/core/command.c:14662`                      | new `pyenv_maint` CLI subcommand dispatching to all four §4.8 operations (rebuild_base, rebuild_script, rebuild_all [-cache], prune [-days=N]) |
| `[x]` | `src/core/command.h:158`                        | `process_pyenv_maint` declaration |
| `[x]` | `src/core/command_def.h:134`                    | `STR_PYENV_MAINT` definition |
| `[x]` | `src/core/command_list.h:158`                   | `pyenv_maint` dispatch entry |
| `[x]` | `python_module/sirilpy/utility.py`              | `_install_package`, `uninstall_package`, new `pip_show`/`pip_list`/`declare_dependencies` helpers |
| `[x]` | `python_module/sirilpy/gpuhelper.py`            | route all pip subprocesses through `_build_install_command`/`_build_uninstall_command`/`pip_list` |
| `[x]` | `python_module/pyproject.toml`                  | version bump 1.1.12 → 1.1.20 (no `uv` dep — see §5.0/§5.4) |
| `[x]` | `python_module/sirilpy/__init__.py`             | export `pip_show`, `pip_list`, `declare_dependencies` |
| `[—]` | `build/windows/native-gitlab-ci/siril-build.sh` | (Windows packaging done separately by the platform maintainer) |
| `[—]` | `build/windows/installer/siril64.iss`           | (Windows packaging done separately)                   |
| `[—]` | `build/macosx/*`                                | (macOS packaging done separately by the platform maintainer) |
| `[x]` | `build/flatpak/org.siril.Siril.json`            | new module: download+install `uv` into `/app/bin` (x86_64 + aarch64 sources, sha256-verified) |
| `[x]` | `build/appimage/generate.sh`                    | fetch+verify+place `uv` in `appdir/usr/bin/` (AppRun already prepends this to PATH) |
| `[x]` | `build/build-image/Dockerfile.debian-oldstable` | install pinned `uv` (gates the AppImage build and CI tests against the oldest supported Python) |
| `[x]` | `build/build-image/Dockerfile.debian-latest`    | install pinned `uv` (CI parity)                       |
| `[—]` | `build/build-image/Dockerfile.win64-latest`     | (Windows cross-build image — Windows packaging done separately) |
| `[x]` | `build/UV_VERSION` (new file)                   | single pinned uv version tag: read by `generate.sh`, mirrored as the `ARG UV_VERSION=` default in each Dockerfile, embedded in Flatpak source URL+sha |
| `[x]` | `siril-scripts/core/GPU_Manager.py` (scripts repo, branch `uv-pip-helpers`, commit `9d580a1`) | routed `subprocess.run([sys.executable, '-m', 'pip', …])` calls through `sirilpy.pip_show` / `sirilpy.ensure_installed` behind a `check_module_version(">=1.1.20")` gate with the legacy subprocess paths retained for Siril 1.4; ROCm wheel install also picks up uv via `shutil.which` when present. Version bumped to 1.2.0. |

End of plan.
