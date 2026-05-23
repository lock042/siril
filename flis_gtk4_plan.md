# FLIS reimplementation plan — `flis-gtk4` branch

Status: **planning, not started**. No code is to be written until this plan
is reviewed and approved.

---

## 0 · Context and ground rules

### Background

* The original implementation lives on the `flis` branch (1 of `origin/flis`).
  It is a GTK3 implementation built against a pre-separation master where
  `single` carried `nb_layers`/`fit` and GUI/non-GUI code was intermingled.
* Since then `master` has been heavily restructured. In particular:
  * GUI and non-GUI code are completely separated.
    Non-GUI code calls into the GUI only via `core/gui_iface.h`.
  * `single` was simplified to `{filename, fileexist, comment, nb_layers, fit}`.
  * A new GTK4 tree (`src/gui-gtk4/`) is being grown in parallel with
    `src/gui/` (GTK3). The branch we are on (`flis-gtk4`) is the GTK4
    line. Only the GTK4 GUI is required; the GTK3 panel is treated as the
    *behavioural* reference, not as code to port.
  * `src/gui-gtk4/image_display.c` implements a tiled `GtkSnapshot`
    rendering pipeline with per-tile `GdkTexture`s materialised lazily
    from the `view->buf` BGRA buffer, with mip selection and a soft
    eviction cap. This is the opportunity for accelerated compositing.

### Scope, in/out

| In scope                                                                      | Out of scope                                |
| ----------------------------------------------------------------------------- | ------------------------------------------- |
| FLIS v1.0 format read/write, fully spec-compliant                             | New blend modes beyond spec                 |
| Single-window GTK4 layers panel + per-layer property editor                   | Porting the GTK3 panel verbatim             |
| FLIS-aware display compositing, integrated with the tiled snapshot path      | GtkGLArea / custom GL shaders               |
| Layer mask + processing mask plumbing (mask tab integration)                  | New mask types beyond spec                  |
| Layer groups (PASS_THROUGH + NORMAL)                                          | Nested groups (spec forbids them)           |
| Mono tinting (LAYER_COLOR) and LRGB chroma mode                               | Adjustment layers (spec reserved, not v1.0) |
| Undo/redo coverage for every mutating operation                               | Sparse-layer write support beyond what flis branch already has |
| Full geometry awareness (rotate / mirror / crop / resample / binning)        | Sequence-level FLIS (only single-image)     |
| Save dialog FLIS file type, file-browser preview                              |                                             |
| `is_current_image_flis()` guards across all relevant tool/processing paths    |                                             |

### Reuse policy (the user's directive, paraphrased)

* **No new FLIS GUI work in `src/gui/`.** The whole `src/gui/`
  tree is the GTK3 implementation; on this branch it isn't built
  by default (only `src/gui-gtk4/` compiles when meson is configured
  with `gtk=gtk4`, the branch default). All *FLIS-specific* GUI work
  — the layers panel, mask dialogs, registration dialog, layer-match
  dialog, image_display compositing integration, etc. — lives
  exclusively in `src/gui-gtk4/`. Do not port or duplicate any of
  it into `src/gui/`.

  **Exception — "keep-it-building" mirror edits are allowed.** When
  a change to a shared header or core source (`src/core/siril.h`,
  `src/core/processing.{c,h}`, `src/io/…`, etc.) renames or
  restructures a symbol that the GTK3 tree references, mirror the
  minimal change into the affected `src/gui/` file(s) so the GTK3
  build continues to compile. Such edits are **mechanical,
  non-functional, and never introduce FLIS-aware behaviour into the
  GTK3 code path** — they exist purely to keep the GTK3 tree
  consistent with shared headers. Example (stage 1.1): the
  `single.nb_layers` → `single.chans` rename touches
  `src/gui/callbacks.c:1295` because that file reads the field;
  the edit changes only the field name, nothing else. No FLIS
  layer logic, no `flis_*` calls, and no widget changes belong in
  `src/gui/` under any circumstance.
* **Non-GUI code may be copied & adapted** from `flis:src/io/image_format_flis.*`,
  `flis:src/core/undo.{c,h}`, `flis:src/core/processing.{c,h}` (generic_layer_worker),
  parts of `flis:src/core/siril.h` (the `single` extensions + `historic` extensions),
  `flis:src/algos/geometry.c` layer-offset helpers, and similar.
  Wherever the flis branch pulled in `gui/*.h` directly, the GTK4 port
  must instead go through `core/gui_iface.h` or move the call into
  `src/gui-gtk4/`. **No new direct includes of `src/gui*` from non-GUI code.**
* **GUI code must be rewritten from scratch in `src/gui-gtk4/`** in
  idiomatic GTK4: no `GtkDialog`/`gtk_widget_show_all`/`gtk_container_add`,
  no `.glade`/`.ui` files copied verbatim. Reference the GTK3 panel
  in `flis:src/gui/` only for *which* widgets exist, *which* signals
  they wire, and *what* they do — never copy code across.
* A previous unplanned attempt left ~6,300 lines of FLIS code in the
  working tree (`src/io/image_format_flis.{c,h}`, `src/io/flis_compose.c`,
  `src/gui-gtk4/flis_gui.{c,h}`, `src/gui-gtk4/flis_gpu_compose.{c,h}`).
  Those files have been deleted. The plan starts from a clean tree and
  derives all FLIS code from the `flis` branch (non-GUI) or writes it
  fresh (GUI). The scaffold files must not be resurrected — if they are
  useful as reference, consult them via `git show` against an older
  filesystem snapshot, do not check them back in.

### Working principles

* Bottom-up: data model → format I/O → compositing kernel → display
  integration → GUI shell → per-operation FLIS-awareness → undo coverage
  → polish. Each stage produces a build that compiles, runs, and either
  exposes new functionality or holds existing functionality unchanged.
* No deferred work *across* stages. Each stage may have known limitations
  (called out per stage), but the build is shippable at the end of every
  stage. The plan deliberately orders work so the limitations of stage N
  are addressed by stage N+1 or N+2, never indefinitely.
* Tests are written *as part of* each stage, not bolted on. Where the
  existing repo has a `src/tests/` tree we add to it; where we need
  fixtures we add small synthetic FITS/FLIS files.
* Commits land per logical sub-stage, with stage tags
  (`flis-gtk4: stage 1.3 — layer write path`) so reviewers can navigate.

### Definition of done (whole project)

1. Open every example file in the `flis` branch's test corpus and
   round-trip it (open → save → reopen → byte-equivalent metadata,
   pixel-equivalent layers).
2. Every interaction listed in §4 (GUI mapping) is reachable from the
   GTK4 panel and behaves equivalently to the GTK3 panel.
3. Every interaction listed in §4 (GUI mapping) is *also* reachable
   from the Siril command parser (`process_*` table in `command.c`),
   so a `.ssf` script can drive an identical FLIS edit session
   headlessly with no UI. See §C (Command interface) for the master
   list; this is a hard parity requirement, not "nice-to-have".
4. `is_current_image_flis()` guards all operations that would corrupt
   layer state (raw gfit mutations, sequence load, etc.) with either a
   layer-aware code path or a refusal + user-facing message.
5. **In GUI mode**, undo/redo covers every mutating panel operation,
   every processing operation that touches a layer, every geometry op,
   and every group op. **In headless mode**, no undo state is saved
   (Siril design principle, see §C.1a). The same primitives serve both
   modes; the headless suppression is centralised in the undo
   machinery, not duplicated per call site, so command-driven and
   panel-driven operations have identical undo behaviour in whichever
   mode the process is running.
6. Display performance: a 4-layer 24 Mpix FLIS pans and zooms at the
   same frame rate as a single 24 Mpix FITS in the same build (i.e.
   compositing is not on the per-frame critical path).
7. No regressions in plain-FITS open/save/process flows. Validated by
   running through the project's existing scripted test set.
8. No GTK3-period direct GUI includes from non-GUI code.
9. No FLIS GUI work in `src/gui/` (the GTK3 tree); all FLIS GUI
   work lives in `src/gui-gtk4/`. The GTK3 tree may carry
   mechanical "keep-it-building" mirror edits per the §0 exception
   (e.g. shared-field renames) but must contain no FLIS-aware
   logic. Verifiable by inspecting `git diff master...flis-gtk4 --
   src/gui/`: each chunk should be a 1–2 line field-rename or
   signature-update with no `flis_*` symbols introduced.

---

## C · Command interface (parallel track)

This section defines the full FLIS command surface in one place. The
commands are **not** delivered as a separate stage; they ship
progressively, woven into stages 1, 2, 4, and 5, paired with the
underlying primitives so that every command exists from the moment
the feature it controls exists. The rationale:

* **Headless testability from day one.** Stage-1 introspection
  commands (`flis_info`, `flis_layer_list`, …) let us drive
  end-to-end tests through `siril-cli` and `.ssf` scripts without
  needing a single line of GUI code. They are a load-bearing piece
  of the stage-1 test harness.
* **GUI ↔ command parity by construction.** When a panel button is
  added in stage 4, its matching command lands in the same PR. The
  panel handler and the `process_*` function call into the same
  non-GUI primitive (the setters in `image_format_flis.c`), so they
  cannot drift. This is enforced by a per-feature test that runs the
  command in headless mode and asserts the resulting `single` state
  matches the state produced by simulating the equivalent panel
  click in a GUI smoke test.
* **No orphan code paths.** Every command in §C.2 maps to exactly
  one primitive in `image_format_flis.c` and to exactly one panel
  control in §4.2 (or to a CLI-only path explicitly noted as such).

Python coverage is automatic: `sirilpy`'s
`SirilInterface.cmd()` method dispatches strings through the same
Siril command parser used by `.ssf` scripts, so the entire
`flis_*` command surface is immediately callable from Python with
no additional binding work. No sirilpy-specific FLIS sub-project
is required.

### C.1 — Conventions

* **Naming**: `flis_*` prefix throughout. Verb after the prefix.
  Setter form is `flis_set_<property>` to match Siril's existing
  `set*` convention (e.g. `setfindstar`). Group-targeting commands
  use `flis_set_group_<property>` to disambiguate from layer
  setters.
* **Layer identification**: every command that targets a single
  layer accepts `<id|"name">` as the first positional argument.
  A bare integer is interpreted as `ITEM_ID`; anything quoted or
  containing non-digits is interpreted as `LAYER_NAME`. If the
  name is ambiguous (two layers share it) the command errors out
  with a message listing matching IDs; the user resolves by
  passing the ID. Names are matched case-sensitively, exactly
  (no fuzzy / glob).
* **Group identification**: same scheme as layers but `<gid|"name">`,
  scoped to the group namespace.
* **Boolean arguments**: `0`/`1` or `false`/`true`, case-insensitive.
  This matches Siril's existing convention (`set <key> <0|1>`).
* **Floating-point**: `0.0`–`1.0` for opacity and tint; `0`–`100`
  is *not* accepted (avoid the percentage ambiguity that bites the
  panel; commands always use the underlying [0,1] scale).
* **Error reporting**: every command returns `CMD_OK` /
  `CMD_FAILURE_TEST_ERROR` (or the equivalent existing enum) and
  emits a `siril_log_color_message` line on error explaining the
  cause. The same string format is used for both "no FLIS loaded"
  (which is a hard refusal) and "layer not found" (also a hard
  refusal).
* **Refusal when not in FLIS mode**: most `flis_*` commands require
  an active FLIS. Use a new `REQ_CMD_FLIS_IMAGE` flag (modelled on
  `REQ_CMD_SINGLE_IMAGE`) registered in `command_list.h` so the
  parser refuses the command early with a uniform message. The
  exception is `flis_promote`, which works on a plain FITS.
* **Undo behaviour**: every mutating layer operation — whether
  driven from a panel handler or from a `process_flis_*` command —
  is dispatched through `generic_layer_worker` (§1.5). The worker
  is the single owner of the undo-save call, mirroring the existing
  pattern in `generic_image_worker` (see `src/core/processing.c:1621`,
  where `undo_state = args->fit == gfit && !(args->custom_undo ||
  args->for_preview || com.script)` gates the call). Primitives in
  `image_format_flis.c` (`flis_layer_set_opacity` etc.) do **not**
  call `undo_save_flis_*` themselves; the worker calls the
  appropriate `undo_save_flis_*` variant around the hook
  invocation, based on flags in `generic_layer_args`. See §C.1a
  for the headless behaviour that falls out of this.

### C.1a — Undo and headless mode

Siril design principle: **undo state must not be saved when the
process is headless or running a script.** Saving swap files and
history entries in a non-interactive run is pure waste — there's
no UI to drive undo/redo, the process exits when the script
finishes, and the swap files become orphans that have to be
cleaned up. The principle is already implemented for the existing
generic workers: `generic_image_worker` gates on `com.script` at
line 1621 of `src/core/processing.c`, and `generic_mask_worker`
follows the same pattern.

Implementation rule for the FLIS work:

* `generic_layer_worker` (the new framework introduced in §1.5
  and used by every mutating layer op) **gates the undo save on
  exactly the same condition** as `generic_image_worker`: skip
  when `com.script` (script run, includes headless and Python),
  skip when `args->custom_undo` (caller manages own undo), skip
  when a notional `args->for_preview` (reserved for future
  preview-style layer ops). The gate lives in the worker, not in
  the undo primitive, not in `image_format_flis.c` setters, and
  not in the per-operation hooks.
* The primitives in `image_format_flis.c` (setters, layer_add,
  layer_remove, …) are pure state mutators. They never call
  `undo_save_flis_*`. They are safe to call from any context.
* The `undo_save_flis_*` family (added in §1.4) is called only
  from `generic_layer_worker` and from the small number of
  GUI-only special-case paths that cannot run headless by
  construction (the opacity-slider drag-end snapshot, etc. —
  enumerated in §1.4). Those special cases live in `src/gui-gtk4/`
  and therefore never execute headless; they need no explicit
  gate of their own.
* Consequence for the parity contract: command-driven and
  panel-driven operations produce *identical* undo behaviour —
  identical undo entries when running under the GUI (both go
  through `generic_layer_worker`), and identical no-op when
  running headless (both hit the worker's `com.script` gate).
  The parity test in §4.3 verifies state changes (mode-independent);
  undo coverage is verified by a separate GUI-mode test (§4.4).

### C.2 — Master command table

The "Stage" column is the stage in which the command (and its tests)
land. The "Primitive" column is the function in
`image_format_flis.c` (stage 1) the command wraps. The "Panel"
column is the §4.2 widget that drives the same primitive from the
GUI.

| Command                            | Stage | Synopsis                                                                  | Primitive                                | Panel widget                          |
| ---------------------------------- | :---: | ------------------------------------------------------------------------- | ---------------------------------------- | ------------------------------------- |
| **Introspection (read-only)**      |       |                                                                           |                                          |                                       |
| `flis_info`                        | 1.6   | `flis_info` — canvas WxH, layer/group count, ICC, capabilities            | (read-only)                              | implicit (window title, mode label)   |
| `flis_layer_list`                  | 1.6   | `flis_layer_list [-format=text\|csv]` — table of id, order, name, mode, opacity, visible, locked, mono/rgb, tint, has-lmask, has-pmask, group_id, pos_x, pos_y | (read-only)                              | layer list rows                       |
| `flis_group_list`                  | 1.6   | `flis_group_list` — table of group id, name, visible, opacity, blend, member count, collapsed | (read-only)                              | group rows                            |
| `flis_layer_info`                  | 1.6   | `flis_layer_info <id\|"name">` — verbose dump of one layer's props        | (read-only)                              | properties panel display              |
| `flis_group_info`                  | 1.6   | `flis_group_info <gid\|"name">` — verbose dump of one group               | (read-only)                              | group selection display               |
| **File operations**                |       |                                                                           |                                          |                                       |
| `flis_promote`                     | 2     | `flis_promote [-name="base"]` — promote loaded plain FITS to a 1-layer FLIS | `flis_promote_from_gfit`                | toolbar `+` (when in plain-FITS mode) |
| `load`, `save`                     | 2     | (existing commands; FLIS-aware via stage 2 dispatch)                      | `load_flis` / `save_flis`                | open/save dialogs                     |
| **Layer stack**                    |       |                                                                           |                                          |                                       |
| `flis_layer_add`                   | 4.3   | `flis_layer_add "filename" [-name="X"]` — load FITS as new top layer      | `flis_layer_add`                         | `flis_add_btn`                        |
| `flis_layer_remove`                | 4.3   | `flis_layer_remove <id\|"name">`                                          | `flis_layer_remove` + `flis_undo_purge_layer` | `flis_remove_btn`                |
| `flis_layer_duplicate`             | 4.3   | `flis_layer_duplicate <id\|"name">`                                       | `flis_layer_duplicate`                   | `flis_duplicate_btn`                  |
| `flis_layer_move_up`               | 4.3   | `flis_layer_move_up <id\|"name">`                                         | `flis_layer_move_up`                     | `flis_move_up_btn`                    |
| `flis_layer_move_down`             | 4.3   | `flis_layer_move_down <id\|"name">`                                       | `flis_layer_move_down`                   | `flis_move_down_btn`                  |
| `flis_layer_export`                | 4.3   | `flis_layer_export <id\|"name"> "filename.fit"`                           | (uses `savefits` on layer fit)           | "Export current layer as FITS…"       |
| **Layer properties**               |       |                                                                           |                                          |                                       |
| `flis_set_name`                    | 4.3   | `flis_set_name <id\|"name"> "new"`                                        | `flis_layer_set_name`                    | `flis_name_entry`                     |
| `flis_set_blend_mode`              | 4.3   | `flis_set_blend_mode <id\|"name"> <mode>` (case-insensitive; the spec §5 names: normal\|multiply\|screen\|overlay\|soft_light\|hard_light\|color_dodge\|color_burn\|darken\|lighten\|difference\|exclusion\|hue\|saturation\|color\|luminosity\|lrgb_color\|pass_through) | `flis_layer_set_blend_mode`           | `flis_blend_combo`                    |
| `flis_set_opacity`                 | 4.3   | `flis_set_opacity <id\|"name"> <0.0..1.0>`                                | `flis_layer_set_opacity`                 | opacity slider                        |
| `flis_set_visible`                 | 4.3   | `flis_set_visible <id\|"name"> {0\|1}`                                    | `flis_layer_set_visible`                 | row eye toggle                        |
| `flis_set_locked`                  | 4.3   | `flis_set_locked <id\|"name"> {0\|1}`                                     | `flis_layer_set_locked`                  | row lock toggle                       |
| `flis_set_position`                | 4.3   | `flis_set_position <id\|"name"> {-x= -y= \| -dx= -dy=}` — absolute or relative canvas offset | `flis_layer_set_position`        | `flis_drag_toggle_btn` + canvas drag  |
| **Tint (mono layers)**             |       |                                                                           |                                          |                                       |
| `flis_set_tint`                    | 4.3   | `flis_set_tint <id\|"name"> -r= -g= -b=` (each 0.0..1.0)                  | `flis_layer_set_tint`                    | `flis_tint_color_btn`                 |
| `flis_clear_tint`                  | 4.3   | `flis_clear_tint <id\|"name">`                                            | `flis_layer_clear_tint`                  | `flis_tint_check` (uncheck)           |
| **Layer mask**                     |       |                                                                           |                                          |                                       |
| `flis_lmask_add`                   | 4.3   | `flis_lmask_add <id\|"name"> "filename"` (the mask file must match layer dimensions) | `flis_layer_set_lmask`        | `flis_mask_toggle_btn` "Add…"         |
| `flis_lmask_remove`                | 4.3   | `flis_lmask_remove <id\|"name">`                                          | `flis_layer_remove_lmask`                | `flis_mask_toggle_btn` "Remove…"      |
| `flis_lmask_move`                  | 4.3   | `flis_lmask_move <from> <to>`                                             | `flis_layer_move_lmask`                  | `flis_mask_move_btn`                  |
| `flis_lmask_set_active`            | 4.3   | `flis_lmask_set_active <id\|"name"> {0\|1}`                               | (direct field, then invalidate)          | `flis_mask_status_btn`                |
| **Groups**                         |       |                                                                           |                                          |                                       |
| `flis_group_add`                   | 4.3   | `flis_group_add [-name="X"]`                                              | `flis_group_add`                         | `flis_group_btn`                      |
| `flis_group_remove`                | 4.3   | `flis_group_remove <gid\|"name">` (members become ungrouped)              | `flis_group_remove`                      | context menu                          |
| `flis_group_delete_with_layers`    | 4.3   | `flis_group_delete_with_layers <gid\|"name">` (destructive; no undo)      | `flis_group_delete_with_layers`          | context menu (with confirmation)      |
| `flis_layer_set_group`             | 4.3   | `flis_layer_set_group <id\|"name"> <gid\|0>` (0 = ungroup)                | `flis_layer_set_group`                   | "Move layer to group…"                |
| `flis_set_group_name`              | 4.3   | `flis_set_group_name <gid\|"name"> "new"`                                 | `flis_group_set_name`                    | group row name                        |
| `flis_set_group_visible`           | 4.3   | `flis_set_group_visible <gid\|"name"> {0\|1}`                             | `flis_group_set_visible`                 | group row eye toggle                  |
| `flis_set_group_opacity`           | 4.3   | `flis_set_group_opacity <gid\|"name"> <0.0..1.0>`                         | `flis_group_set_opacity`                 | group properties                      |
| `flis_set_group_blend_mode`        | 4.3   | `flis_set_group_blend_mode <gid\|"name"> {normal\|pass_through}`          | `flis_group_set_blend_mode`              | group properties                      |
| **Composite operations**           |       |                                                                           |                                          |                                       |
| `flis_merge_down`                  | 4.3   | `flis_merge_down <id\|"name">` (destructive; purges undo for the two layers; refuses if locked) | `flis_merge_down_layer`     | "Merge Down"                          |
| `flis_flatten`                     | 4.3   | `flis_flatten` (destructive; purges all undo)                             | `flis_flatten_all`                       | "Flatten Image"                       |
| **Integration commands (stage 5)** |       |                                                                           |                                          |                                       |
| `flis_register_layers`             | 5.6   | `flis_register_layers [-ref=<id\|"name">] [-transf=…] [-interp=…]`         | `flis_register_layers` worker            | "Register layers…" dialog             |
| `flis_layers_match`                | 5.7   | `flis_layers_match [-subset=<id\|"name">,<id\|"name">,…]`                  | `flis_background_neutralise_layers`      | "Layers match…" dialog                |
| `mask_from_*`                      | 5.2   | existing mask commands gain `-layermask=<id\|"name">` option              | `generic_mask_args.target_layer_id`      | mask dialogs' layer chooser           |

### C.3 — Per-stage delivery slots

| Stage | Commands shipping this stage                                                                                                                       |
| :---: | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1.6   | All introspection commands (`flis_info`, `flis_layer_list`, `flis_group_list`, `flis_layer_info`, `flis_group_info`). Read-only; no undo concerns. |
| 2     | `flis_promote`. `load` / `save` become FLIS-aware via existing dispatch (stage 2.1 / 2.2).                                                          |
| 3     | None. Display is internal; the introspection commands from §1.6 already let scripts assert the displayed image state via `flis_info`-derived metrics. |
| 4.3   | The bulk: all stack / property / tint / lmask / group / composite commands. Each ships in the same commit as the panel widget that calls the same primitive. |
| 5     | `flis_register_layers`, `flis_layers_match`, the `-layermask=` extension on existing mask commands.                                                 |
| 6     | None new; the `flis_*` commands are user-documented in §6.5 alongside the panel.                                                                    |

### C.4 — Per-command test pattern

Each mutating command lands with a `test_cmd_flis_<name>` headless
test that:

1. Loads a known fixture FLIS into a `siril-cli` instance.
2. Runs the command via the script parser (not via direct C call —
   this is the parity guarantee).
3. Asserts the post-state via the introspection commands' machine-
   readable form (`-format=csv`), comparing against a golden CSV
   committed to the fixture set.

No undo/redo step at this layer — per §C.1a, headless runs save no
undo state, so `undo`/`redo` would be no-ops. Undo coverage is
tested separately in GUI mode (§4.4 `test_panel_<verb>_undo_redo`).

Read-only commands (`flis_info` etc.) get a simpler test that just
asserts output against golden text.

Test fixtures live in `src/tests/fixtures/flis/` and are small
(8×8 or 32×32 pixels) so the test suite stays fast.

### C.5 — Documentation

* `command_def.h` `STR_FLIS_*` entries documenting each command's
  argument list with the same `<b>` / `<i>` markup pattern used by
  the existing commands. These feed `/help <command>` in the script
  console.
* `command_list.h` table entries with proper `REQ_CMD_FLIS_IMAGE` /
  `REQ_CMD_SINGLE_IMAGE` flags.
* `docs/` user guide (stage 6.5) gains a "scripting FLIS" section
  that walks through a worked example: load mono frames, promote
  one to FLIS base, add others as tinted mono layers, set blend
  modes, save. End-to-end script in the user's own `.ssf` file.

---

## 1 · Data model & format I/O (non-GUI)

### 1.1 — `single` and `historic` extensions in `core/siril.h`

Add to `single`:

* `GSList *layers;` ordered ascending by `layer_order`.
* `GSList *groups;` unsorted.
* `gint active_layer;` index into `layers`.
* `gint next_item_id;` initialised to 1.

Replace `nb_layers` with `int chans;` (mirrors active layer's
`fit->naxes[2]`). Keep `fit` as a convenience alias for the active layer.

Add to `historic_struct` the FLIS-aware fields listed in
`flis:src/core/undo.h`'s "IMPORTANT — siril.h change required" block
(`flis_layer_id`, `layer_props`, `lmask_*`, `reorder_*`, `flis_position_*`,
`pmask_only`, `full_layer`, `multi_entries`/`n_multi_entries`).

Forward-declare `flis_layer_t`, `flis_group_t`, `flis_layer_props_t`,
`layermask_t` in `siril.h` (no include of `image_format_flis.h` here).

**Migration**: every existing read of `single.nb_layers` becomes a
read of `single.chans` (or, where the caller wanted the layer count,
of the new `flis_layer_count()`). In practice this is fewer than ten
call sites — well under the original "some dozens" estimate.

`src/gui/` references the field too; per the "keep-it-building"
exception in §0, mirror the rename into `src/gui/` so the GTK3 tree
still compiles. The edit is a literal field-name change with no
FLIS semantics added.

Tests:
* Compile-only. Reorder/rename verification: grep
  `com\.uniq->nb_layers` / `single->nb_layers` across all source
  trees and confirm zero hits.

### 1.2 — `image_format_flis.{c,h}` adopted from `flis` branch

Bring the non-GUI FLIS implementation in:

* Format read (`load_flis`)
* Format write (`save_flis`)
* Layer/group/mask data structures + lifecycle
* All `flis_layer_*` / `flis_group_*` setters/getters
* Canvas geometry helpers (`flis_canvas_rx/ry`, `flis_canvas_to_pixel_index`)
* ICC conversion (`flis_convert_layers_icc`)
* `is_current_image_flis()`
* Sort/lookup helpers
* `flis_promote_from_gfit` (the entry point that turns a freshly-loaded
  FITS into a single-layer FLIS in memory — needed by both the open
  path and "promote to FLIS" panel actions)

Add to `src/meson.build`'s `src_files`:
```
'io/image_format_flis.c',
'io/flis_compose.c',
```

Tests (`src/tests/`):
* `test_flis_roundtrip` — write a synthetic 2-layer FLIS, read it back,
  assert metadata identical and pixel data identical.
* `test_flis_metadata` — verify METADATA column key-value preservation
  (including unknown keys round-tripping per spec §11.2).
* `test_flis_geometry_helpers` — `flis_update_layer_offset_after_*`
  with crop/resize/rotate against known geometry.

**Flatten and merge-down pixel-correctness tests** (depend on §1.3's
kernel being in place; both sets of tests are landed together at the
end of §1.3). These are the user-facing complement to §1.3's
white-box kernel tests: they exercise the same blend arithmetic but
through the actual `flis_flatten_all` and `flis_merge_down_layer`
pipelines a script or panel button would invoke, so they catch bugs
in flatten plumbing (layer deletion, base preservation, mask
consumption, undo purge) as well as bugs in the kernel. Constant-
colour input layers keep the expected output a closed-form value, so
the assertions compare against the W3C blend reference formulas (the
same formulas used by the §1.3 kernel tests; identical reference
means a divergence is unambiguously a bug in flatten plumbing, not
in the test itself).

* `test_flis_flatten_per_blend_mode` — parameterised over every
  blend mode in spec §5 (Normal, Multiply, Screen, Overlay, Soft
  Light, Hard Light, Color Dodge, Color Burn, Darken, Lighten,
  Difference, Exclusion, Hue, Saturation, Color, Luminosity, LRGB
  Chroma). For each mode: build a 2-layer FLIS with constant-colour
  base (br, bg, bb) and top (tr, tg, tb) at opacity 1.0; call
  `flis_flatten_all()`; assert the resulting single layer's pixels
  match `mode_formula(base, top, 1.0)` within ≤1e-5 (float layers)
  or ≤1 ULP (word layers).
* `test_flis_flatten_per_opacity` — same matrix at opacities
  {0.0, 0.25, 0.5, 0.75, 1.0}. 0.0 must produce base unchanged;
  1.0 reduces to the per-blend-mode test above; intermediates
  exercise the Porter-Duff over operator.
* `test_flis_flatten_with_lmask` — three sub-cases per blend mode:
  (a) lmask all-white = same as no lmask; (b) lmask all-mid-grey
  (0.5) = same as opacity 0.5 with no lmask; (c) lmask gradient
  along x = boundary correctness (pixel at x=0 unchanged, pixel at
  x=W-1 fully blended). Asserts the lmask `MASK_ACT` flag is also
  honoured: with `MASK_ACT=F`, the mask is ignored.
* `test_flis_flatten_with_tinted_mono` — mono top layer with
  LAYER_COLOR={1.0, 0.2, 0.1} (Hα red) on an RGB base in Screen
  blend mode; assert per-channel result is
  `screen(base[c], top * tint[c])`.
* `test_flis_flatten_with_group_pass_through` — 2-layer FLIS where
  the top layer is in a PASS_THROUGH group; assert the group's
  visible/opacity scaling is applied and the result equals
  flattening with the layers ungrouped at those effective
  parameters.
* `test_flis_flatten_with_group_normal` — same fixture but group
  blend mode = NORMAL; assert the group's members are composited
  into a temporary buffer first, then blended into the base — and
  that the result differs from the PASS_THROUGH case when the
  group contains layers using non-NORMAL modes (the test that
  catches "I forgot the intermediate buffer" bugs).
* `test_flis_flatten_sparse_layer` — top layer with
  POSITION_X/Y != 0 and extent smaller than the canvas; assert
  pixels outside the layer's extent equal the base unchanged, and
  pixels inside equal the per-mode blend. Catches off-by-one and
  clipping bugs in the sparse-layer code path.
* `test_flis_flatten_invisible_layer_skipped` — top layer with
  `visible=FALSE`; assert result equals base unchanged regardless
  of blend mode.
* `test_flis_flatten_post_state` — after `flis_flatten_all`, assert:
  `flis_layer_count() == 1`; base layer's `blend_mode == NORMAL`,
  `opacity == 1.0`, `visible == TRUE`; base layer's lmask is NULL
  and processing mask is NULL; the entire undo history has been
  purged; all non-base layers' item_ids are gone (no dangling
  references in groups, no orphan undo entries).
* `test_flis_merge_down_per_blend_mode` — same parameterised
  matrix as `test_flis_flatten_per_blend_mode` but using
  `flis_merge_down_layer` on a 2-layer FLIS, asserting the result
  has 1 layer with the expected merged pixels. Also asserts that
  merge-down purges the undo history of the two layers involved
  but leaves the histories of unrelated layers (in a 3+ layer
  fixture) untouched.
* `test_flis_merge_down_consumes_top_lmask` — top layer has an
  lmask; merge-down; assert the lmask is consumed (applied during
  the composite) and the merged layer has no lmask carried over.

### 1.3 — `flis_compose.c` adopted (compositing kernel, non-GUI)

The kernel takes a fully-prepared float[3][H][W] working buffer and
walks the visible layer stack writing into it. Required outputs:

* Display path (called from the GUI to fill a BGRA `view->buf` or
  per-tile buffer — see stage 3).
* Merge-down path (composite one layer onto another and stash the
  result as the new pixels of the lower layer).
* Flatten path (composite the whole stack onto the base layer).
* Thumbnail path (small composite for HDU-0 thumbnail at save time).

API shape (subject to audit refinement):

```c
/* Composite into a caller-owned float[3 * canvas_w * canvas_h] buffer.
 * @subset: NULL = all visible layers; else GSList<flis_layer_t*> filter. */
int flis_compose_to_float_rgb(float *out_rgb,
                              guint canvas_w, guint canvas_h,
                              GSList *subset);

/* Same but writes to a BGRA8 buffer with stretch lookup applied. */
int flis_compose_to_bgra8(uint8_t *out_bgra,
                          guint canvas_w, guint canvas_h,
                          GSList *subset,
                          const struct stretch_lut *lut);

/* Tile-aware variant for the GTK4 tiled-snapshot path (added in stage 3). */
int flis_compose_to_bgra8_region(uint8_t *out_bgra, int stride,
                                 int x0, int y0, int w, int h,
                                 GSList *subset,
                                 const struct stretch_lut *lut);
```

The region variant is the lever for tile-level on-demand compositing
(stage 3.3).

Tests (white-box, kernel direct):
* `test_flis_compose_blend_modes` — for each blend mode in spec §5,
  feed two synthetic constant-colour layers and check the result
  against the W3C reference formulas at a few sample points.
* `test_flis_compose_mono_tint` — mono layer + LAYER_COLOR composited
  over RGB; confirm channel mapping per spec §6.5.
* `test_flis_compose_lmask` — layer mask scales effective alpha
  pixelwise; per-pixel coverage matches the mask values.
* `test_flis_compose_group_pass_through` and `_normal` — group blend
  semantics.
* `test_flis_compose_region` — region matches the full composite when
  re-tiled, to nail down tile-boundary correctness.

The black-box counterpart — same blend arithmetic exercised
through the user-visible `flis_flatten_all` and
`flis_merge_down_layer` pipelines — lives in §1.2's
`test_flis_flatten_*` and `test_flis_merge_down_*` suites. Both
white-box and black-box tests use the same W3C reference
formulas, so a divergence between the kernel and the flatten
pipeline is unambiguous evidence of a plumbing bug rather than a
test-spec disagreement. A third witness — the GPU snapshot path —
is added in §3.5 (`test_display_composite_pixel_equivalence`).

### 1.4 — Undo plumbing in `core/undo.{c,h}`

Bring in from `flis:src/core/undo.c`:
* `flis_undo_purge_layer` (purges all states for a removed item_id)
* `undo_save_flis_layer_props` and `_snapshot` variant
* `undo_save_flis_lmask` and `undo_save_flis_lmask_move`
* `undo_save_flis_layer_reorder`
* `undo_save_flis_multi_layer` (compound: pixels + props for N layers)
* `undo_save_flis_multi_layer_props` (compound: props only, for group drags)
* `undo_save_processing_mask` (pmask-only)
* `undo_save_flis_layer_full` (geometry-changing single-layer state)

Plus the consumer side: `undo_display_data`/`undo_get_data` must
dispatch on the new historic-entry flags (`props_only`, `pmask_only`,
`full_layer`, `multi_entries != NULL`).

**Headless gating — not in these functions.** Per §C.1a, the
`undo_save_flis_*` family does **not** gate on `com.script` /
`com.headless`. They unconditionally push to the stack, mirroring
how the existing `undo_save_state` (`src/core/undo.c:422`) behaves.
The gate lives one level up in `generic_layer_worker` (§1.5) and in
`generic_image_worker` (`src/core/processing.c:1621`), where
`com.script` is already checked. Callers outside the worker
framework — currently a small enumerated set of GUI-only paths
listed in §1.5 — never execute headless by construction and need
no explicit gate.

Tests (all run in GUI mode, with the harness ensuring the
`generic_layer_worker` undo-gating condition is satisfied):
* `test_undo_props_only` — save props, mutate, undo, assert restored.
* `test_undo_lmask_round_trip` — add/remove/move lmask, undo, redo.
* `test_undo_multi_layer` — apply registration-style multi-layer change,
  undo single-step restores all layers.
* `test_undo_purge_on_remove` — remove a layer, confirm its prior
  states are gone but other layers' states remain.
* `test_undo_gated_by_worker_when_script` — separately, see §1.5
  for `test_layer_worker_skips_undo_when_script` covering the gate.

### 1.5 — `generic_layer_worker`: the layer-operation framework

This is the centralisation point for **all** mutating layer
operations. It mirrors `generic_image_worker` and
`generic_mask_worker` and is the single home for:

1. **Undo gating** — exactly the same condition as
   `generic_image_worker` line 1621: skip undo save when
   `args->custom_undo || args->for_preview || com.script`. This is
   the headless gate from §C.1a, in one place.
2. **Undo dispatch** — pick the correct `undo_save_flis_*` variant
   based on flags in `generic_layer_args`:
   * `updates_lmask` → `undo_save_flis_lmask`
   * `geometry_changing` → `undo_save_flis_layer_full`
   * affects multiple layers → `undo_save_flis_multi_layer` (with
     a `layers` GSList field) or `_multi_layer_props` for
     property-only multi-layer ops
   * default → `undo_save_flis_layer_props`
3. **Log hooks** — call `args->log_hook(args->user, DETAILED)` to
   produce the detailed message for the log; call again with
   `SUMMARY` to produce the short label used for the undo entry.
   Identical pattern to `generic_image_worker:1619-1623`.
4. **HISTORY card update** — write the detailed message to the
   active layer's `fit->history` so it's preserved in the FITS
   header on save, regardless of whether undo was saved.
5. **Thread coordination** — runs on the same processing thread as
   the existing workers, respects `start_in_new_thread` /
   `processing_is_job_active`, and posts an end-of-op idle.
6. **Idle posting** — by default `end_generic_layer` posts a refresh
   that calls `flis_display_invalidate(FLIS_INV_STACK)` and
   `flis_gui_update_from_idle()`. Override with
   `args->idle_function` for special cases (the few that need it).
7. **Headless completion** — when `com.headless`, skip the idle
   and free `args` directly, identical to
   `generic_image_worker:1652-1655`.
8. **Error reporting** — uniform progress-bar update and log line
   based on `args->retval`.

Files touched:
* `src/core/processing.{c,h}` — add `struct generic_layer_args` and
  `gpointer generic_layer_worker(gpointer)` plus
  `free_generic_layer_args` and the default `end_generic_layer`
  idle. Adapted from the flis branch but updated to use the current
  `gui_iface` for any GUI feedback (no direct `gui/*.h` includes).
* Additions to `generic_img_args` (the `geometry_changing` flag,
  consumed by the existing `generic_image_worker` so that
  geometry-changing image ops save the right *layer-level* undo
  variant when in FLIS mode) and `generic_mask_args`
  (`target_layer_id`, consumed by `generic_mask_worker` to route
  mask output to a layer's lmask).

Calling convention (the contract for every mutating layer op):

```c
struct generic_layer_args *args = calloc(1, sizeof(*args));
args->layer        = target_layer;       /* primary layer */
args->layer_hook   = my_op_hook;         /* called on worker thread */
args->log_hook     = my_op_log;          /* DETAILED + SUMMARY strings */
args->user         = op_specific_data;
args->description  = g_strdup("My operation");
args->updates_lmask    = FALSE;          /* set TRUE if op writes lmask */
args->geometry_changing = FALSE;         /* set TRUE for rotate/crop/etc */
start_in_new_thread(generic_layer_worker, args);
```

The hook calls the appropriate `flis_layer_set_*` primitive(s) and
returns 0 / non-zero. The worker handles the rest. Both panel
handlers and `process_flis_*` commands use this exact pattern.

Special-case undo-save callers outside the worker (enumerated so
the audit in §4.5 can verify the list is exhaustive):

* `flis_opacity_drag_*` in `src/gui-gtk4/flis_gui.c` — captures a
  property snapshot on drag-press and calls
  `undo_save_flis_layer_props_snapshot` on drag-release, to avoid
  pushing one undo entry per slider tick. GUI-only by construction;
  drag events cannot occur in a script.
* Equivalent drag-end patterns for any other slider/colour-picker
  widgets in §4.2 that batch many setter calls into one undo entry.

The `processing_should_continue` / thread-helper churn shown in the
flis-branch diff was unrelated to FLIS and is dropped from this port.

Tests:
* `test_generic_layer_worker_calls_hook` — minimal hook fires, undo
  state saved (in GUI mode), idle posted.
* `test_layer_worker_skips_undo_when_script` — set `com.script = TRUE`,
  run a property-change layer op, assert no undo entry pushed and
  no swap file created. Then set `com.script = FALSE` and rerun,
  assert undo entry pushed.
* `test_layer_worker_dispatches_undo_variant` — for each combination
  of `updates_lmask` / `geometry_changing` / multi-layer flags,
  assert the correct `undo_save_flis_*` variant was invoked.
* `test_layer_worker_log_hook_messages` — assert DETAILED string
  appears in the log and SUMMARY string appears as the undo label.

### 1.6 — Introspection commands (read-only)

Implement the read-only `flis_*` commands from §C.2 ("Introspection")
in `src/core/command.c`. These have no UI counterpart and exist purely
to enable headless inspection of the FLIS state — they are the test
harness for everything in §1.1–§1.5.

Files touched:
* `src/core/command.c` — new `process_flis_info`, `process_flis_layer_list`,
  `process_flis_group_list`, `process_flis_layer_info`,
  `process_flis_group_info`.
* `src/core/command_list.h` — register the five commands. Use the new
  `REQ_CMD_FLIS_IMAGE` flag (defined here for the first time; rejects
  the command with a uniform message when no FLIS is loaded).
* `src/core/command_def.h` — `STR_FLIS_INFO` etc. for each.
* `src/core/command.h` — prototypes.

Output format: human-readable by default; `-format=csv` for the
table-form commands (`flis_layer_list`, `flis_group_list`) emits a
machine-parseable form suitable for diffing against golden files in
tests. The CSV header columns are fixed by spec so test fixtures
remain stable across reorderings.

Tests:
* `test_cmd_flis_info_golden` — load fixture, run `flis_info`,
  assert against golden text.
* `test_cmd_flis_layer_list_csv` — load fixture, run
  `flis_layer_list -format=csv`, byte-compare to golden CSV.
* `test_cmd_flis_layer_info_golden` — same for `flis_layer_info`.
* `test_cmd_flis_*_refuses_non_flis` — load a plain FITS, run each
  command, assert refusal with `REQ_CMD_FLIS_IMAGE` message.

### 1.7 — Stage 1 build & integration checkoff

- [ ] `meson compile` clean on full tree (headless and GUI configs).
- [ ] Headless smoke: `siril-cli --version` + open + save round-trip
      of a tiny test FLIS produced by `test_flis_roundtrip`.
- [ ] Every blend mode in spec §5 verified correct through three
      independent witnesses: kernel direct (§1.3
      `test_flis_compose_blend_modes`), `flis_flatten_all` pipeline
      (§1.2 `test_flis_flatten_per_blend_mode`), and
      `flis_merge_down_layer` pipeline (§1.2
      `test_flis_merge_down_per_blend_mode`). All three agree to
      within the per-test tolerances.
- [ ] Flatten post-state assertions all pass: 1-layer result, base
      reset to NORMAL/1.0/visible, masks cleared, undo purged
      (§1.2 `test_flis_flatten_post_state`).
- [ ] `siril-cli` can run a `.ssf` script that opens a FLIS, runs
      `flis_info` + `flis_layer_list -format=csv`, and writes the
      expected output. Captured as `test_cmd_flis_introspection_e2e`.
- [ ] No `#include "gui*.h"` from any file outside `src/gui*/`.
- [ ] All `src/tests/test_flis_*` and `test_cmd_flis_*` pass.

At end of stage 1: **no UI integration yet, but full headless
inspection is possible.** The build is fully FLIS-aware in its core
types; a script can open a FLIS, inspect its layer structure, and
save it back unchanged. The next stage wires the loader path through
the open/save dialogs so opening a FLIS through the UI becomes
possible.

---

## 2 · File-open and save plumbing

### 2.1 — Open path (`src/io/image_format_fits.c` + `gui-gtk4/open_dialog.c` + file_browser preview)

In `image_format_fits.c::readfits`, detect `FLIS = T` on the primary HDU
and delegate to `load_flis`. (Pattern copied from flis branch but
adapted: do not include `image_format_flis.h` from a generic header path
chain that bleeds GUI.) Add `.flis` extension to the FITS-like
extension recognisers used by the open dialog / file browser preview.

Update `src/gui-gtk4/file_browser.c` to render thumbnails from the
HDU-0 thumbnail when the file is FLIS (already a small 8-bit RGB),
which is faster than building one from gfit.

Tests:
* Manual: open one of the FLIS sample files; observe the base layer
  loaded into gfit, `is_current_image_flis()` returns TRUE,
  `com.uniq.layers` is non-empty.
* `test_image_format_fits_dispatches_flis` — calling `readfits()` on
  a FLIS file results in `load_flis` being invoked.

### 2.2 — Save path (`src/gui-gtk4/save_dialog.c`)

Add `TYPEFLIS` to `image_type` (in `core/image_type.h` or wherever it
lives), wire it through `get_filetype`/`get_image_type_from_filter`/
`get_type_from_filename` (`.flis`) / `get_type_for_extension`.

Save dispatch: when `type_of_image == TYPEFLIS` (or filename ends in
`.flis`), call `save_flis(filename)` instead of `savefits`. The save
dialog must offer the FLIS filter only when `is_current_image_flis()`
is TRUE (or when a single-layer flatten is acceptable — but for v1
we restrict to the FLIS case to avoid silent flattening).

Behaviour from flis branch to preserve: if the user is in plain-FITS
mode and selects "Save as FLIS…", we promote gfit to a single-layer
FLIS (`flis_promote_from_gfit`) and save. The user is informed via
status bar message.

Tests:
* Manual: open plain FITS, save as `.flis`, reopen — verify a single
  base layer matches the original.
* Manual: open multi-layer FLIS, save as `.fit` — verify save_dialog
  refuses or warns (depending on UX decision at audit).

### 2.3 — `flis_promote` command

`flis_promote [-name="X"]` (see §C.2) — promotes the currently loaded
plain FITS to a single-layer FLIS in memory by calling
`flis_promote_from_gfit`. The base layer takes the given `-name=` or
defaults to `"Background"`. After the command, `is_current_image_flis()`
returns TRUE and subsequent `flis_*` commands work on the image.

This ships in stage 2 (not stage 4.3) because it's the natural pair to
the save-dialog "Save as FLIS…" auto-promote behaviour from §2.2 — both
exist to bridge plain FITS into the FLIS world, and they should land
together so scripts can promote without going through the save dialog.

Files touched:
* `src/core/command.c` — `process_flis_promote`.
* `src/core/command_list.h` / `command_def.h` — registration + help.
  Note: this command uses `REQ_CMD_SINGLE_IMAGE` (not the new
  `REQ_CMD_FLIS_IMAGE`), since it precisely operates on the non-FLIS
  case.

Tests:
* `test_cmd_flis_promote` — load plain FITS, `flis_promote`, assert
  `flis_info` shows a 1-layer FLIS with canvas == FITS dimensions.
* `test_cmd_flis_promote_idempotent_refusal` — calling it twice
  errors (already a FLIS).
* `test_cmd_flis_promote_then_save_then_load` — promote, save as
  `.flis`, close, reopen — assert single base layer, pixel-identical
  to the original FITS.

### 2.4 — Stage 2 checkoff

- [ ] File browser preview shows FLIS thumbnails.
- [ ] Open dialog accepts `.flis` and routes to `load_flis`.
- [ ] Save dialog offers `.flis` and routes to `save_flis`.
- [ ] Plain-FITS → FLIS promotion works on save dialog AND via
      the `flis_promote` command, producing byte-identical output.
- [ ] No FLIS-specific compositing yet — only the *base* layer is
      displayed at this stage. The next stage fixes that.

---

## 3 · Display integration — tile-aware composite

This is the headline GTK4 win. The goal is for composite rebuild to
**not** rebuild a full canvas-sized BGRA buffer every time a single
layer property changes, by hooking into the existing
`siril_image_view_snapshot` tiled pipeline.

### 3.1 — Composite as the source of `view->buf` (correct, slow path)

First, get correctness end-to-end with the simplest possible plumbing:

* `remap_all_vports()` / `remap(vport)` detect FLIS mode and call
  `flis_compose_to_bgra8` to fill `view->buf` from the composite,
  rather than from `gfit`. (`gfit` continues to point at the active
  layer for tool/processing purposes.)
* The existing tile materialisation path is unchanged.

This works but rebuilds the whole composite on every property change.
That's fine as a baseline.

### 3.2 — Per-layer GdkTexture cache (GSK blend path, fast for spec-compliant blends)

Adapt or rewrite (per audit) `flis_gpu_compose.{c,h}`. Per visible layer
that uses a GSK-translatable blend mode (Normal, Multiply, Screen,
Overlay, Darken, Lighten, Color Dodge, Color Burn, Hard Light, Soft
Light, Difference, Exclusion, Hue, Saturation, Color, Luminosity), keep
a cached canvas-sized `GdkTexture` with the stretched pixel data. In
`siril_image_view_snapshot`, push `gtk_snapshot_push_blend(mode)` per
layer and append each texture, instead of painting `view->buf` tiles.

Cache invalidation:
* Layer pixel change → drop that layer's texture.
* Layer property change (opacity, blend, visible, tint, mask): drop
  texture only if the change affects what's *in* the texture (tint
  does; opacity/blend/visible don't — those are pushed at composite
  time).
* Stretch lo/hi change → drop *all* textures.

Fallback to the 4.1 CPU composite for any frame that contains:
* A layer with `has_tint && fit->naxes[2] == 1` (tint requires our
  shader, GSK can't do it),
* A layer with active layer mask (needs alpha mask, GSK push_mask
  exists but pre-CHROMA),
* `FLIS_BLEND_CHROMA` (no GSK equivalent — LRGB),
* Group with non-PASS_THROUGH blend (requires intermediate buffer).

Document the fallback in the cache code's header comment so it's clear
which cases are GPU and which are CPU.

### 3.3 — Tile-aware composite, per-layer textures (the real win)

**Architectural decision (locked in):** integrate the per-layer
`GdkTexture` cache from §3.2 into the existing tiled snapshot tree.
For each visible layer in stack order:

1. Maintain per-layer **per-tile** `GdkTexture`s rather than one
   canvas-sized texture, materialised and evicted by the same lazy
   pipeline that already governs `gui.view[vport].tiles`. Each layer
   gets a parallel `tiles[]` array keyed by `(tx, ty)` covering only
   the canvas region the layer overlaps (sparse-layer aware via
   `position_x/y`).
2. At snapshot time, walk layers in `layer_order`. For each visible
   tile of each layer push a `gtk_snapshot_push_blend(mode)` frame
   for the layer's blend mode and append the per-tile texture inside
   it. Group containers push/pop matching blend frames.
3. Per-layer global opacity is applied via
   `gtk_snapshot_push_opacity` around each layer's tile sequence.
4. Layer masks are pushed as a `gtk_snapshot_push_mask` frame using
   a parallel per-tile mask-texture cache; LMASK_INTENSITY is the
   correct mask mode for alpha-style masks.
5. Mono tint and LRGB chroma — both unrepresentable in GSK natively
   — are handled by baking the contribution at *texture-build time*:
   the per-tile texture for a tinted mono layer already contains
   tint-multiplied RGB pixels; for `FLIS_BLEND_CHROMA`, the
   contribution is computed against the prior composite as a custom
   tile-level kernel using `gtk_snapshot_append_cairo` only for the
   chroma layer (rare; not on the hot path). Tint changes invalidate
   that layer's tile textures (already part of the invalidation
   protocol in §3.4).

Why B, not the per-tile CPU composite alternative:

* Modern GSK lowers most blend frames to GPU shaders. Compositing on
  the GPU avoids round-tripping the whole canvas through CPU memory
  on every property change.
* The existing tile pipeline already handles VRAM pressure: tiles
  evict under the same soft-cap budget. Layered images multiply tile
  count by layer count, but tiles outside the visible rect are not
  materialised, so the working set scales with **visible area ×
  visible layers**, not total canvas area.
* The CPU `flis_compose_to_bgra8` kernel from §1.3 is *not* unused —
  it is the correctness oracle (testing) and the implementation for
  merge-down / flatten / thumbnail / save (offscreen paths). The
  display path uses it only when fallback is triggered (see below).

Testing rigour for this stage is elevated because the per-layer-tile
machinery is the most subtle code in the project:

* **Pixel-equivalence harness** — a headless test renders a fixture
  through the GPU snapshot path (offscreen GSK render) and compares
  pixel-for-pixel against `flis_compose_to_bgra8` for the same
  stretch. Tolerance: ≤1 ULP per channel in BGRA8 (GSK uses
  premultiplied alpha internally; the test accounts for round-trip
  conversion). Driven from `bench/flis_compose_pixel_equivalence`,
  covering every blend mode, every group config, masked and unmasked
  layers, with and without tint.
* **Invalidation correctness** — per-layer mutations (toggle visible,
  change opacity, change blend, change tint, modify lmask, modify
  pixels, reorder) trigger exactly the expected
  `flis_display_invalidate(...)` flags and no others. Asserted via a
  test mode that records invalidation calls and compares against
  hand-authored expectations.
* **VRAM stress** — 8-layer 100 Mpix synthetic FLIS, pan + zoom for
  60 s, assert no unbounded growth in resident texture bytes and no
  visible tile-drop artefacts.

Fallback path (single retained CPU path): when any of the following
hold, drop the per-layer GPU compose for the affected stack and
render the whole composite through the §3.1 CPU path:

* Unknown layer blend mode (forward-compat, spec §11.2).
* GSK renderer version too old (runtime check at startup; falls back
  the entire session and logs once).
* Render budget exceeded — if a snapshot pass takes >50 ms three
  frames in a row, drop the GPU path for that image and re-enable on
  next image. This is a safety net; profiling in §3.5 should keep
  it from triggering in normal use.

### 3.3a — Per-layer tile materialisation details

Concretely, the additions to the rendering state:

```c
/* Parallel to struct image_view, one per FLIS layer. */
struct flis_layer_view {
    gint           item_id;       /* matches flis_layer_t->item_id */
    int            tile_cols;     /* derived from layer extent on canvas */
    int            tile_rows;
    struct image_tile *tiles;     /* same layout as image_view->tiles */
    guint64        lazy_epoch;
    /* Pixel-source cache: BGRA8 layer-image after stretch+tint+lmask
     * bake, in canvas coordinates.  Filled lazily per-tile by the same
     * materialise_tile pattern used for view->buf. */
    uint8_t       *buf;           /* canvas-sized when eager; tile-local when lazy */
    /* Reverse-link for invalidation: which tiles cover sparse-layer pixels. */
};
```

The materialise path for a layer tile mirrors the existing single-image
`materialise_tile`: take a region of the layer's source pixels, apply
the active stretch + tint + lmask, encode to BGRA8, upload as a
`GdkTexture`. The encode is **per-layer per-tile**, so a change to a
single layer rebuilds only that layer's affected tiles.

`siril_image_view_snapshot` becomes, in FLIS mode:

```c
/* Pseudocode — actual implementation per stage 3.3. */
for (layer in com.uniq->layers ordered by layer_order) {
    if (!layer->visible) continue;
    gtk_snapshot_push_opacity(snap, effective_opacity(layer));
    gtk_snapshot_push_blend(snap, gsk_mode_for(layer->blend_mode));
    if (layer->lmask_active && layer->lmask) {
        gtk_snapshot_push_mask(snap, GSK_MASK_MODE_LUMINANCE);
        for_each_visible_tile(layer) append_mask_tile(...);
        gtk_snapshot_pop(snap);  /* mask top → becomes the mask */
    }
    for_each_visible_tile(layer) append_layer_tile(...);
    gtk_snapshot_pop(snap);   /* blend */
    gtk_snapshot_pop(snap);   /* opacity */
}
```

GSK's `push_blend` semantics treat the "top" subtree as the source
and the prior frame's accumulated content as the destination, which
matches the FLIS spec's ordering (`LAYER_ORDER` ascending = bottom up).

### 3.4 — Invalidation hooks

A single `flis_display_invalidate(flis_invalidate_flags_t)` function
that the GUI panel + processing operations call after each mutation.
Flags:
* `FLIS_INV_LAYER_PIXELS(item_id)` — re-build that layer's cache and
  any tiles it covers.
* `FLIS_INV_LAYER_PROPS(item_id)` — no texture rebuild, just queue
  redraw.
* `FLIS_INV_STACK` — order/visibility changed; full re-walk.
* `FLIS_INV_ALL` — drop every cache, full rebuild.

This is the single chokepoint everywhere else hooks into; it makes
display correctness auditable.

### 3.5 — Tests and benchmarks

* `bench/flis_compose_4layer_24mp` — script that opens a fixture FLIS,
  toggles each layer's visibility 10× and reports redraw latency.
  Baseline number captured at start of stage 3 (CPU full-canvas);
  improvement tracked after 4.2 and after 4.3.
* `test_display_composite_pixel_equivalence` — headless: composite a
  fixture and compare against `flis_compose_to_float_rgb` reference.

### 3.5a — FLIS-aware ICC architecture (port from `flis` branch)

**Motivation.** A mono FLIS base in our GTK4 branch was getting a *Gray*
sRGB-TRC profile auto-assigned by `icc_auto_assign_or_convert` on stretch-
tool open. `flis_render_layers` then copied that gray profile onto the
always-RGB composite, and the display pipeline applied a Gray→Monitor
proofing transform to the 3-plane composite — collapsing R/G/B and
erasing every tinted-layer or RGB-layer contribution above the base.

**Architecture (mirrors `flis` branch).** The FLIS invariant is:

* the **base layer** is the only fits* in a FLIS that carries an ICC
  profile, and that profile describes the **RGB composite** — not the
  (possibly mono) base data itself;
* non-base layers always have `icc_profile == NULL`,
  `color_managed == FALSE`;
* the **composite is always 3-channel RGB float**, so its assigned
  profile must also be RGB.

The three helper functions live in `image_format_flis.{c,h}` and are
the seams every ICC call site uses:

* `flis_get_profiled_fit()` — base layer fit (or `gfit` for non-FLIS).
* `flis_composite_naxes2()` — 3 for any FLIS, `gfit->naxes[2]` otherwise.
* `flis_convert_layers_icc(old, new)` — RGB↔RGB conversion across the
  whole stack (mono layers broadcast→transform→collapse via Rec.709),
  plus tint vector conversion.

ICC call sites updated in `core/icc_profile.c`:

* `icc_auto_assign_or_convert` / `icc_auto_assign` redirect from `gfit`
  to `flis_get_profiled_fit()` and choose mono-vs-RGB target profile
  from `flis_composite_naxes2()` (not `fit->naxes[2]`).  This lets a
  mono FLIS base get the RGB working profile — matching the composite.
* `siril_colorspace_transform`: when *assigning* a profile to a FLIS
  base, use `flis_composite_naxes2()` for the channel-compat check.
  Add a re-tag-only safety net for the case where data and assigned
  profile disagree on channel count (mono-base-with-RGB-profile is a
  legitimate FLIS state — running an RGB transform on a mono buffer
  would read OOB).
* `refresh_icc_transforms` / `initialize_proofing_transform` source
  their input profile and channel format from the profiled fit and
  `flis_composite_naxes2()`, not from `gfit` (which is typically the
  active non-base layer with `icc_profile == NULL`).
* `color_manage()` fires the toolbar update when called on the FLIS
  profiled fit, since `gfit != base` is the common FLIS case.

GUI updates: `check_gfit_profile_identical_to_monitor` in
`gui-gtk4/image_display.c` compares the profiled fit against the
monitor profile, not `gfit`.

Load/add invariants (already in `image_format_flis.c`):

* `load_flis`: file-level ICC profile is assigned to the base only;
  any per-layer profile that snuck into a non-base HDU is logged and
  discarded.
* `flis_layer_add`: incoming layer's pixels are converted to the
  base's profile via `siril_colorspace_transform` first, then the
  per-layer profile is stripped.

**Deferred: GTK4-native shader color management.** The current path
runs `cmsCreateProofingTransform` over the composite pixels on the
CPU once per remap (gated by `do_transform`, skipped when source
primaries match the monitor).  For very large FLIS images, the
proofing-transform pass is non-trivial and runs on every redraw that
hits the eager path.  A future stage can move this into the GTK4
snapshot/shader pipeline:

* upload the composite (or the per-layer GdkTextures from §3.2) once;
* express the proofing transform as a 3D LUT texture or a sequence of
  GSK colour-matrix nodes, applied to the textures on the GPU;
* a per-frame swap of the LUT is cheap, no CPU pass needed.

GSK does not currently expose a native colour-managed texture node,
so this would either need a GSK glshader (Linux/Mac) or a custom
render node.  Cross-platform support and HDR considerations make this
a non-trivial undertaking; record as a deferred §3.7 work item and
revisit only if the CPU proofing transform becomes a bottleneck on
realistic FLIS workloads.

### 3.6 — Stage 3 checkoff

- [ ] Multi-layer FLIS displays correctly with all blend modes.
- [ ] Toggling visibility on a hidden layer in a 4-layer 24 Mpix
      FLIS is visually instantaneous (no full-canvas rebuild stall).
- [ ] The composite matches `flis_compose_to_float_rgb` pixel-for-pixel
      after stretch.
- [ ] Selecting a layer in the panel rebinds gfit but does not rebuild
      the composite cache.
- [ ] Opening a stretch tool (asinh / GHS / MTF / …) on a mono base
      FLIS keeps every tinted / RGB layer above visible — i.e. the
      auto-assigned profile lands on the base as RGB and the
      composite renders end-to-end through an RGB proofing transform.
- [ ] No regressions in plain-FITS display.

### 3.7 — Deferred: GTK4-native shader colour management

See §3.5a closing paragraphs.  Not blocking any current Stage 3 goal.

### 3.5b — ICC storage moved to com.uniq (post-§3.5a follow-up)

The §3.5a architecture lived on `fit->icc_profile` / `fit->color_managed`,
with per-layer profile redundantly carried on the FLIS base.  This
follow-up moves the authoritative store to `com.uniq` (the `single`
struct) and **removes** the per-fits ICC fields from the `fits` struct.

Why: ICC undo entries on the per-fits + FLIS-multi-layer architecture
were unnecessarily heavy — undo_save_flis_multi_layer captures every
layer's pixel data via swap files even for ICC operations that never
touch pixels (Assign / Remove).  Moving the profile to com.uniq lets
ICC undo entries snapshot one `cmsHPROFILE` pointer + one boolean,
with no swap-file overhead.

The change in one paragraph: `com.uniq->icc_profile` +
`com.uniq->color_managed` are the single source of truth.  A new set
of accessors in `core/icc_profile.h` (`current_icc_profile`,
`current_image_color_managed`, `current_image_set_icc_profile`,
`current_image_clear_icc_profile`, `current_image_color_manage`)
fronts that storage.  Per-fits helpers (`fit_get_icc_profile`,
`fit_get_color_managed`) return the current-image profile for `gfit`
and the FLIS profiled fit, and NULL/FALSE for intermediate buffers
and sequence frames.  A lightweight `historic->icc_only` undo
flavour with a one-line restore branch handles ICC undo for Assign
and Remove (Convert keeps the multi-layer flavour because it does
rewrite layer pixels).

Migration scope: ~430 call sites across 44 files were updated to
use the accessors, then the two `fits` struct fields were removed.

Knock-on feature regressions (accepted in the migration):
* **Compositing / Pixelmath / Remixer cross-layer profile matching**
  removed.  Input layers are treated as raw pixel data; if a colour-
  managed workflow is needed the user must pre-convert.  These tools
  were not in routine FLIS-era use and the alternative — keeping a
  per-fits profile while everything else moves to com.uniq — was
  rejected in favour of a clean store.
* **Sequence-export ICC handling** removed.  Sequence frames carry
  no profile; this matches the project policy "colour management is
  excluded from sequence operations".
* **flis_layer_add cross-profile conversion** removed.  A newly added
  layer's pixels are assumed to be in the FLIS's colour space; if not,
  the user should convert before adding or via Image → Color Management
  after.

GTK3 sibling files (`src/gui/*.c`) still reference the removed fields
in dead code paths.  They are not built and have not been migrated;
the GTK4 cutover in the broader Siril roadmap will retire them.

Open follow-ups (small):
* Delete `flis_get_profiled_fit` once the few remaining call sites
  (mostly internal to the FLIS-aware ICC accessors and the dialog
  worker dispatch) are migrated to compare against `gfit` directly.
* Consider promoting `fits_initialize_icc` / `check_profile_correct`
  to dedicated load-helpers that take an out-param profile instead
  of conditionally writing to com.uniq based on `fit == gfit`.

---

## 4 · GUI shell — GTK4 layers panel (rewrite)

This stage is the largest single hand-coded chunk. It is **written
from scratch** in GTK4 idioms.

### 4.1 — Panel construction

**Decision (locked in):** the panel is **built procedurally** in
`src/gui-gtk4/flis_gui.c`. Rationale: the layer list and group rows
are dynamic, the property panel's sensitivity flips with selection
type (no layer / layer / group), and the mask sub-frame morphs based
on which mask types are present. Driving all of that from a single
`.ui` file would require either layering signals onto `GtkBuilder`
templates with substantial code-side rewiring on every refresh, or
embedding so many runtime widget swaps that the static description
adds little. Procedural construction keeps construction and refresh
co-located in one file.

The only widget that justifies a `.ui` file is the row template if we
choose to use `GtkBuilderListItemFactory` for `GtkListView` rows; a
small `flis_layer_row.ui` for that one purpose is acceptable. The
panel window, toolbar, property panel, and dialogs are all
procedural.

Top widget:
* `GtkWindow` (no `GtkDialog` — GTK4 deprecates that pattern),
  modeless, hide-on-close, `transient-for` set in code.
* `GtkApplication` action `app.show-layers` (accelerator: Ctrl+L by
  default) toggles panel visibility.
* The panel is allocated lazily on first show and never destroyed
  until app exit; refresh is driven by `flis_gui_update()` /
  `flis_gui_update_from_idle()`.

### 4.2 — Full GUI mapping (every GTK3 widget/menu mapped)

This is the master checklist. Every entry from
`flis:src/gui/uifiles/flis_layers.ui` and the context menu must have
an explicit GTK4 equivalent below. **No widget may be silently
omitted.**

#### Header / mode indicator

| GTK3 widget                 | GTK4 equivalent                     | Notes                                                                 |
| --------------------------- | ----------------------------------- | --------------------------------------------------------------------- |
| `flis_layers_window`        | `GtkWindow`, hide-on-close          | No headerbar; modeless                                                |
| `flis_mode_label` "FITS"    | `GtkLabel`, `.dim-label` CSS class  | Live-updated: "FITS" / "FLIS" / "FLIS·active group"                   |

#### Layer list

| GTK3                                       | GTK4                                                                 | Notes                                                                                  |
| ------------------------------------------ | -------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| `GtkListBox flis_layer_list`               | `GtkListView` with a `GListStore` of layer model objects             | GTK4-idiomatic; supports rubber-band selection if we later want it                     |
| Row visibility toggle (eye icon)           | `GtkToggleButton` with `view-reveal-symbolic` icon                   | Wires `flis_layer_set_visible` + `flis_display_invalidate(FLIS_INV_STACK)`             |
| Row lock toggle (padlock icon)             | `GtkToggleButton` with `changes-prevent-symbolic` icon               | Wires `flis_layer_set_locked`                                                          |
| Row thumbnail                              | `GtkPicture` populated lazily; 32×32                                 | Generated from layer pixels (downsampled) at first use, cached on the layer model      |
| Group row collapse arrow                   | `GtkExpander`                                                        | Hides child rows in the list when collapsed                                            |
| Drag-to-reorder                            | `GtkDragSource`/`GtkDropTarget` on rows                              | Drops between rows recompute `layer_order` to insert (`prev.order + next.order)/2`     |

#### Toolbar

| GTK3 widget id                | GTK4 widget                          | Action                                                                                  |
| ----------------------------- | ------------------------------------ | --------------------------------------------------------------------------------------- |
| `flis_add_btn`                | `GtkButton` "list-add-symbolic"      | Opens FITS chooser, calls `flis_layer_add` (or `flis_promote_from_gfit` if plain FITS)  |
| `flis_remove_btn`             | "list-remove-symbolic"               | `flis_layer_remove` (after undo-purge for that item_id)                                 |
| `flis_duplicate_btn`          | "edit-copy-symbolic"                 | `flis_layer_duplicate`                                                                  |
| `flis_group_btn`              | "folder-new-symbolic"                | `flis_group_add`                                                                        |
| `flis_drag_toggle_btn`        | `GtkToggleButton` "✥"                 | Enters "drag layer in canvas" mode; image_interactions.c handles drag → `position_x/y` updates with undo |
| `flis_move_up_btn`            | "go-up-symbolic"                     | `flis_layer_move_up` with `undo_save_flis_layer_reorder`                                |
| `flis_move_down_btn`          | "go-down-symbolic"                   | `flis_layer_move_down`                                                                  |
| `flis_layer_menu_btn`         | `GtkMenuButton` w/ `GMenu`           | Replaces the GtkMenu context menu (which GTK4 removed)                                  |

#### Property panel (per-layer)

| GTK3                          | GTK4                                                                                | Action                                                                                          |
| ----------------------------- | ----------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `flis_name_entry`             | `GtkEntry` (max-length 32)                                                          | activate / focus-leave → `flis_layer_set_name` with `undo_save_flis_layer_props`                |
| `flis_blend_combo`            | `GtkDropDown` with 18 named modes                                                   | changed → `flis_layer_set_blend_mode` + undo                                                    |
| `flis_opacity_scale`/`spin`   | `GtkScale` + `GtkSpinButton`, shared `GtkAdjustment` 0..100                         | drag-begin records props snapshot; drag-end posts `undo_save_flis_layer_props_snapshot` once (avoids per-tick undo flood) |
| `flis_tint_frame`             | `GtkFrame`                                                                          | Sensitivity follows layer's mono/RGB state                                                      |
| `flis_tint_check`             | `GtkCheckButton`                                                                    | Enables tint colour; calls `flis_layer_set_tint(1,1,1)` initially                               |
| `flis_tint_color_btn`         | `GtkColorDialogButton` (GTK4)                                                       | colour chosen → `flis_layer_set_tint(r,g,b)` + undo                                             |
| `flis_tint_hint_label`        | `GtkLabel` "e.g. red for Ha"                                                        | Pure hint                                                                                       |

#### Layer-mask sub-frame

| GTK3                          | GTK4                                                                | Action                                                                                |
| ----------------------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| `flis_mask_status_btn`        | `GtkButton` with mask-name label                                    | Toggles `lmask_active` (with undo); right-click shows context (remove, properties)    |
| `flis_mask_toggle_btn`        | `GtkButton` "Add…" / "Remove…"                                      | "Add" opens chooser → load mask → `flis_layer_set_lmask`; "Remove" → confirmation → `flis_layer_remove_lmask` |
| `flis_mask_move_btn`          | `GtkButton` "Move…"                                                 | Opens layer chooser → `flis_layer_move_lmask` with compound undo                      |
| `flis_mask_view_row`          | `GtkBox` (visible only when both proc-mask and lmask exist)         |                                                                                       |
| `flis_mask_view_proc_radio` / `flis_mask_view_layer_radio` | Two `GtkCheckButton` in a group                  | Choose which mask the global mask tab + tint overlay shows                            |

#### Context menu items

| GTK3 menu item                            | GTK4 `GMenu` entry                                                                                  |
| ----------------------------------------- | --------------------------------------------------------------------------------------------------- |
| `flis_export_layer_item`                  | "Export current layer as FITS…" — `save_layer_as_fits()`                                            |
| `flis_register_layers_item`               | "Register layers…" — opens the FLIS-aware registration dialog (stage 5)                             |
| `flis_background_neutralise_item`         | "Layers match…" — opens layer-match dialog → `flis_background_neutralise_layers`                    |
| `flis_assign_group_item`                  | "Move layer to group…" — `flis_layer_set_group`                                                     |
| `flis_merge_down_item`                    | "Merge Down" — `flis_merge_down_layer` (purges undo for both layers; confirmation dialog)           |
| `flis_flatten_item`                       | "Flatten Image" — `flis_flatten_all` (purges entire undo; confirmation)                             |

### 4.3 — Operation commands (paired with panel widgets)

For every widget added in §4.2, the matching `flis_*` command from §C.2
("Layer stack" / "Layer properties" / "Tint" / "Layer mask" / "Groups" /
"Composite operations" sections) ships in the same commit. The
parity contract:

* The panel signal handler does *no* mutation directly. It builds an
  argument struct and calls the primitive in `image_format_flis.c`.
* The `process_flis_*` command parses arguments and calls the
  **same** primitive with the **same** argument struct.
* Both the panel handler and the `process_*` function build a
  `generic_layer_args` and call `start_in_new_thread(generic_layer_worker, args)`.
  The worker (§1.5) handles undo-save (gated identically to
  `generic_image_worker` — skipped when `com.script`), invalidation,
  error reporting, logging, and HISTORY update. The primitive
  (`flis_layer_set_opacity` etc.) is a pure mutator called from the
  hook function on the worker thread.
* Therefore both paths produce identical state changes, identical
  log lines, identical HISTORY entries, and — in GUI mode —
  identical undo entries. In script mode neither saves undo; both
  still mutate state and update HISTORY identically.

The full list of commands shipping in this sub-stage is in §C.2 —
all rows tagged "4.3" in the Stage column.

Implementation pattern (per command):
1. Add `process_flis_<verb>` in `src/core/command.c`.
2. Add registration in `command_list.h` with appropriate
   `REQ_CMD_FLIS_IMAGE` flag.
3. Add help string `STR_FLIS_<VERB>` in `command_def.h` documenting
   the argument syntax with the standard markup.
4. Add `test_cmd_flis_<verb>` per §C.4 (golden-CSV diffing for state
   commands; output-text golden for read-only commands). Runs
   headless — no undo step.
5. Wire the panel widget's signal handler to the same primitive.
6. Add `test_panel_<verb>_drives_same_primitive` — a GUI test
   (Wayland headless / Xvfb, with `com.headless = FALSE`) that
   simulates the panel click and asserts post-state matches running
   the command. The same test exercises undo: undo the panel-driven
   change, assert restored; redo, assert reapplied; do the same for
   the command-driven change, then assert the two paths produced
   equivalent undo entries (same `historic_struct` flags, same
   labels) — this is what guarantees undo *parity* in addition to
   state parity.

The parity test in step 6 is the load-bearing piece — it's what
prevents the panel and the command from silently diverging, in both
state and undo behaviour.

### 4.4 — Tests for the panel

All panel tests run in **GUI mode** (`com.headless = FALSE`, Xvfb /
Wayland-headless), so they can exercise undo/redo. This is the
companion track to the headless §C.4 command tests.

* `test_panel_open_closes` — open with no image, with plain FITS,
  with FLIS; close; re-open.
* `test_panel_<verb>_drives_same_primitive` (per §4.3 step 6) —
  for each panel control: simulate the click, capture state;
  rewind; run the equivalent command, capture state; assert state
  equality. Then exercise undo on both, asserting equivalent
  `historic_struct` entries.
* `test_panel_undo_redo_round_trip` — for a representative subset
  of operations (add layer, set opacity, set tint, add lmask,
  reorder, group create, merge down), drive via the panel,
  undo/redo, assert state restored / re-applied each cycle.
* Manual scripted run-through: open multi-layer fixture, drive every
  control, observe display update and undo entry.

### 4.5 — Stage 4 checkoff

- [ ] Every row in §4.2 has a working widget and signal handler.
- [ ] Every mutating widget in §4.2 has a corresponding `flis_*`
      command per §C.2 and a passing `test_cmd_flis_<verb>` test
      (headless, no undo step).
- [ ] Every panel-command pair has a passing
      `test_panel_<verb>_drives_same_primitive` parity test
      (GUI mode, exercises undo equivalence).
- [ ] Mode label correctly reflects FITS / FLIS / group-active.
- [ ] All panel mutations dispatch through `generic_layer_worker`
      (with `args->layer_hook` calling the appropriate
      `flis_layer_set_*` primitive). No panel handler invokes a
      primitive directly outside the worker, and no panel handler
      calls `undo_save_flis_*` directly — except the enumerated
      drag-end snapshot patterns listed in §1.5.
- [ ] Audit: no panel handler bypasses an `image_format_flis.c`
      primitive via direct field assignment (per §7 risk 9).
- [ ] Audit: the special-case undo-save callers (per §1.5) are
      still exactly the enumerated set; no new ones snuck in.
- [ ] All mutations call `flis_display_invalidate(...)` with the
      narrowest applicable flag.
- [ ] Closing the panel does not lose layer state.

---

## 5 · Operations integration — FLIS-awareness across the app

This stage walks the categories of operations that need FLIS hooks.
Each sub-stage is one PR. The categories are derived from the
`flis` branch diff (master..flis showed ~47 files changed); the
non-FLIS-specific refactors in that diff are not ported.

### 5.1 — Geometry ops (crop, rotate, mirror, resample, binning)

`src/algos/geometry.c`, `src/gui-gtk4/menu_gray_geometry.c`:
* Every op that resizes or rotates the active layer must call
  `flis_update_layer_offset_after_*` so non-active layers' positions
  remain correct.
* `generic_img_args.geometry_changing = TRUE` for these ops, so
  `generic_image_worker` saves a `undo_save_flis_layer_full` state.
* Group-level rotation/mirror operates via
  `flis_update_all_layer_offsets_after_rotate` etc.

Tests:
* `test_geometry_crop_updates_layer_offsets`
* `test_geometry_rotate_180_offsets`
* `test_geometry_resize_offsets`

### 5.2 — Masks (`src/gui-gtk4/masks_gui.c`, `mouse_action_functions.c`)

* "Generate mask from chromaticity / image / stars" — port the three
  small dialogs (`mask_from_chromaticity.ui`, `mask_from_image.ui`,
  `mask_from_stars.ui`) into GTK4 `.ui` files in
  `src/gui-gtk4/uifiles/`, route through `generic_mask_worker`.
* Honour `generic_mask_args.target_layer_id`: when set, the mask
  goes to the layer's lmask, not gfit's pmask.
* Mask tab respects the mask-view radio (proc vs layer) from the
  layers panel.

Commands: extend the existing `mask_from_channel`, `mask_from_color`,
`mask_from_lum`, `mask_from_stars` commands with a `-layermask=<id|"name">`
option (mirrors the flis branch's spec). Setting this option routes
the result to the named layer's lmask instead of the processing mask,
identical to setting the GUI mask dialog's "Target: Layer mask of [X]"
combo. Update `command_def.h` STR_MASK_FROM_* help text to mention
the new option (the flis branch already wrote those help strings —
adopt them verbatim).

Tests: `test_cmd_mask_from_channel_layermask_route` — open a FLIS,
run the command with `-layermask=`, assert the named layer's lmask
populated and `gfit->mask` untouched.

### 5.3 — ICC profile (`src/core/icc_profile.c`)

Adapt the flis-branch changes (translated to GTK4 GUI separation):
when the user assigns a profile, call `flis_convert_layers_icc` so
every layer is transformed, not just gfit. The bridge from non-GUI
ICC code to GUI feedback uses `core/gui_iface.h`.

### 5.4 — Single-image / window-title / dialog plumbing

`src/io/single_image.c`, `src/gui-gtk4/callbacks.c`, `dialogs.c`:
* `close_single_image` must free `single.layers`, `single.groups`,
  and call `flis_gpu_compose_free_all`.
* Window title shows `"name.flis — layer 'X' [3/5]"` when in FLIS.
* Modeless layers panel state survives image close (just becomes
  empty); reactivates on next FLIS open.

### 5.5 — Star finder & analysis tools (`src/gui-gtk4/star_finder.c`)

Per flis branch: star finder runs against the **active layer**, not
the composite. Confirm and add a one-line guard.

### 5.6 — Registration ("Register layers…")

Adapt the flis-branch `flis_register_dialog` and the layer-to-layer
registration pipeline. Ports into a GTK4 dialog in
`src/gui-gtk4/uifiles/flis_register_layers.ui`.
Saves a multi-layer compound undo entry.

Command: `flis_register_layers [-ref=<id|"name">] [-transf=…] [-interp=…]`
(see §C.2). Without `-ref=`, the active layer is the reference. The
options mirror the existing `register` command's options where the
semantics carry over. Same primitive as the dialog.

Test: `test_cmd_flis_register_layers_two_layers` — fixture with two
known-misaligned single-star layers, run the command, assert the
non-reference layer's pixels shift to align (correlation peak
detected at offset 0,0 after registration).

### 5.7 — Layer match ("Layers match…" / background neutralise)

Adapt `flis_background_neutralise_layers` and the dialog around it.

Command: `flis_layers_match [-subset=<id|"name">,…]`. With no
subset, all layers participate; with a subset, only those layers
have their scale factors adjusted.

Test: `test_cmd_flis_layers_match_neutral_bg` — fixture with three
tinted mono layers and a known-offset background; run command;
assert composite background is within ε of neutral grey.

### 5.8 — Sequence / livestacking / stacking refusal

`src/io/fits_sequence.c`, `src/livestacking/livestacking.c`,
`src/stacking/stacking.c`:
Refuse to process a sequence whose frames are FLIS; the user must
flatten first. Display a clear error message via the gui_iface.

### 5.9 — Stage 5 checkoff

- [ ] Every geometry op preserves layer relationships and is undoable.
- [ ] Masks targeting a specific layer end up on that layer's lmask
      whether driven from the dialog or via `mask_from_* -layermask=`.
- [ ] ICC re-assign converts all layers.
- [ ] Window title and panel update across image close/open.
- [ ] Star finder and analysis target the active layer.
- [ ] Registration and layer-match dialogs functional and undoable;
      `flis_register_layers` and `flis_layers_match` commands ship
      with passing tests.
- [ ] Sequence operations refuse FLIS frames with a clean message.

---

## 6 · Hardening & polish

### 6.1 — Sparse-layer correctness audit

Spec §6.2 + §12.2 allow layers to extend outside the canvas. Walk
every compositing site, every coordinate transform, and every panel
display to confirm out-of-canvas pixels are clipped, not crash.

### 6.2 — Capability headers on save

`save_flis` writes the §11.1 capability headers reflecting actual
implementation: `FLISCORE=T`, `FLISBLND=131071` (all 17 modes),
`FLISLMSK=T`, `FLISGRP=T`, `FLISICC=T`, `FLISEFF=F`, `FLISEXT='SPARSE'`,
`FLISIMPL='Siril-1.6.0-flis-gtk4'` (or whatever the build sets).

### 6.3 — Forward-compat: preserve-on-resave

When loading, capture every unknown METADATA key, every unknown
ITEM_TYPE row, every unknown HDU header keyword. When saving, emit
them unchanged. This is required by spec §11.2 and §11.4.

Test: `test_flis_unknown_metadata_round_trips` — load a fixture
containing unknown keys, save, load again, assert keys preserved
byte-identical.

### 6.4 — Stress tests

* 16-layer FLIS at 24 Mpix — does it open, display, edit, save in
  acceptable wall time?
* 1-layer 100 Mpix FLIS — tile composite path under VRAM pressure.
* Round-trip every blend mode and every group config.

### 6.5 — Documentation & user-facing

* `docs/` (or wherever Siril docs live): a short FLIS user guide
  with two sections — the GUI walkthrough and a "scripting FLIS"
  worked example (per §C.5) showing an end-to-end `.ssf` that loads
  mono frames, promotes the first to FLIS base, adds the rest as
  tinted layers, sets blend modes, and saves.
* In-app help (status-bar hints in the layers panel) for first-time
  users.

### 6.6 — Final checkoff (whole project DoD)

- [ ] All §0 Definition of Done items satisfied.
- [ ] `git diff master...flis-gtk4 --stat` reviewed; no surprising
      drift in unrelated subsystems.
- [ ] This plan archived under `docs/dev/flis_gtk4_plan.md`.

---

## 7 · Risks & open questions (for review)

1. **GSK pixel-equivalence with the CPU oracle (stage 3.3)** — GSK
   premultiplied-alpha blending and our CPU kernel's straight-alpha
   formulas must produce visually identical results within tolerance.
   If the per-channel ULP test fails for some blend mode, we may
   need to pre-process the per-tile texture (e.g. premultiply
   manually) or bake the blend on the CPU for that specific mode.
   First scope to verify at stage 3 kickoff.
2. **GSK renderer version requirements** — `push_mask` requires
   GTK 4.10+, several enums in `GskBlendMode` require GTK 4.0 but
   their renderer-side implementations differed in early 4.x. The
   plan assumes GTK ≥ 4.12 (Siril's stated GTK4 baseline once the
   port lands). Verify at the start of stage 3.
3. **GTK4 mask group radio** — GTK4 `GtkCheckButton` group API
   differs from GTK3 `GtkRadioButton`; the mask-view radio needs the
   newer `gtk_check_button_set_group` pattern.
4. **GtkColorDialogButton** is GTK 4.10+. Should be fine given the
   baseline above, but flag explicitly at stage 4.
5. **Drag-and-drop to reorder** in `GtkListView` is non-trivial; the
   GTK3 panel used up/down buttons and an "experimental" DnD path.
   v1 of the GTK4 panel ships with up/down only and DnD deferred
   to a tiny follow-up (the one explicit cross-stage deferral).
6. **`save_dialog` filter visibility logic** — should the user be
   able to save a multi-layer FLIS as `.fit` (flattening silently)?
   GTK3 panel did not; recommend same.
7. **Locale handling for blend-mode strings** — the procedural
   panel construction must call `_()` on every user-visible string
   the same way `.ui` `translatable="yes"` does. Pattern to
   establish at start of stage 4.
8. **Per-layer-per-tile VRAM** — 8 layers × full-canvas tiles can
   consume substantial GPU memory at high zoom. The existing tile
   eviction budget must be reviewed to confirm it covers the
   layered case, or extended in stage 3.3 to be layer-count-aware.
9. **Command/GUI parity drift** — the §C parity contract (panel
   handlers and `process_*` functions call the same primitive)
   relies on developer discipline. The §4.3 parity tests will
   catch regressions, but the tests themselves must be maintained.
   At stage 4.5 checkoff, audit that no panel handler bypasses the
   `image_format_flis.c` primitives via direct field assignment.

---

## 8 · Stage delivery summary

| Stage | Deliverable                                                              | Buildable | User-visible change                                                   | Commands shipping (§C)                                              |
| :---: | ------------------------------------------------------------------------ | :-------: | --------------------------------------------------------------------- | ------------------------------------------------------------------- |
| 1     | Core types, FLIS I/O, compositing kernel, undo plumbing + introspection commands |     ✓     | None visually; `siril-cli` can introspect FLIS via 5 read-only commands | `flis_info`, `flis_layer_list`, `flis_group_list`, `flis_layer_info`, `flis_group_info` |
| 2     | Open + save dispatch + promote command                                   |     ✓     | Open FLIS (base only); save plain FITS as FLIS; `flis_promote` script | `flis_promote` (+ `load`/`save` become FLIS-aware)                  |
| 3     | Display composite (CPU oracle → GSK per-layer → tile)                    |     ✓     | Multi-layer FLIS renders correctly                                    | None (display internal)                                             |
| 4     | GTK4 layers panel + per-layer editor + paired operation commands         |     ✓     | Full panel interaction, and every panel operation also scriptable     | ~30 `flis_*` commands (stack / props / tint / lmask / groups / composite) |
| 5     | Geometry/mask/ICC/registration/layer-match integration + integration commands |     ✓     | Operations FLIS-aware; integration dialogs and matching CLI commands  | `flis_register_layers`, `flis_layers_match`, `-layermask=` on mask cmds |
| 6     | Hardening, capability headers, docs (incl. scripting guide)              |     ✓     | Spec-compliant artefacts; user guide covers GUI + scripting           | None new                                                            |

Each stage's checkoff list is the gate for the next stage to begin.
