# Coverity triage — branch covfixes15

Report generated working through `cov.csv` (83 CIDs). Each fixed issue is a separate
commit. This file summarises what was done and what needs your decision.

## TL;DR

- **66 CIDs fixed** across 47 commits (one issue per commit; same-file issues of the same root
  cause are grouped into a single commit that lists every CID). All commits build cleanly and the
  full test suite passes **21/21**; the headless CLI runs.
- **9 CIDs investigated → no change / false positive** (incl. concurrency group (a)).
- **5 CIDs won't-fix** — concurrency group (b): benign lock-free `com.pref.*` config reads
  (per your call).
- **3 CIDs external** (OpenCV system headers) — not actionable.

The concurrency cluster is now resolved: group (a) → false positive, group (b) → won't-fix,
group (c) 647113/647116 → fixed (com.stars locking audit + hardening pass), group (d) 647100 →
fixed (read_single_image no_debayer parameter).

Notable real bugs found (not just analyzer noise): an inverse-FFT channel-loop OOB on RGB+mono
input (647068), an out-of-bounds annotation-visibility array write (647086), a use-after-free in
the file browser (647117), a type-confused wrong free in `mask_from_lum_hook` (641704), an int32
overflow making AVI frame counts negative (530794), and a Python-protocol double-response bug
(530688). Several previously "Fix Submitted" issues were genuinely unresolved and are now fixed.

To re-mark statuses in Coverity: everything under "Fixed" → *Fixed*; "Investigated" → *False
Positive*/*Intentional*; the concurrency group (a) entries → *False Positive*.

## Excluded from work (no action)

These rows are not actionable in this pass:

- **External code (OpenCV headers, not ours)**: CIDs 432129, 585112 (`opencv2/core/matx.inl.hpp`),
  637647 (`opencv2/flann/any.h`). Cannot/should not patch vendored system headers.

Note: rows marked `Triaged` / `Fix Submitted` (641702, 641689, 530672, 530766, 641704,
530688, 641563) are **re-examined**, not skipped — their reappearance in a fresh master
scan means the earlier fix did not fully resolve the defect.

---

## Fixed (one commit each)

- **647072** (newdeconv.c) — removed dead `if (!ret)` error-dialog block (async launcher
  has no sync return). Commit `98dcc87aa`.
- **647119** (documentation.c) — same dead-code pattern removed. Commit `4110f78b2`.
- **647068** (fft.c) — copy-paste error: inverse-FFT channel loop iterated `fit->naxes[2]`
  instead of operand `tmp->naxes[2]`; RGB gfit + mono modulus → OOB channel access. `9363840dd`.
- **647086** (annotation_catalogues.c) — OOB write/read: visibility accessors clamped to
  `CAT_AN_USER_TEMP` (index 11) on an 11-element array; clamped to `CAT_AN_SSO_VECTORS`. `112f5acbc`.
- **647117** (file_browser.c) — UAF: `first_selected_info` returned a pointer it had already
  unref'd; made it transfer-full and unref in all 4 callers. `1e88f1b66`.
- **647109** (background_extraction.c) — OOB read: caller bbox could overrun the luminance
  buffer; clamp selection to image bounds up front. `13252648b`.
- **647073/647089/647093/647104/647115** (icc_profile.c) — leaked `siril_file_chooser_get_filename`
  results (redundant `g_strdup`, discarded NULL-test results). One commit `a66ef8e07`.
- **647061/647070/647122** (plot.c) — AAVSO dialog leaked owned obstype/filter strings on early
  returns and the whole struct on normal exit. `7f8053d6e`.
- **647074** (nina_light_curve.c) — leaked `nina_file` on the no-WCS early return; also removed a
  latent double-free of it on the thread-start failure path. `445d9c17e`.
- **647081/647084/647098** (preferences.c) — three leaked `siril_file_chooser_get_filename`
  results (interface language, swap dir, six catalogue paths). `e8b049ee9`.
- **647077** (annotate.c) — leaked combo text in `get_cat_index_from_combo` on both returns. `0e9b57a41`.

- **647107** (script_menu.c) — `open_script_in_editor` leaked the buffer when
  `g_file_get_contents` succeeded with length 0. `1eada0b7e`.
- **647087/647106** (livestacking.c) — leaked `replace_ext` path and unchecked `any_to_fits`
  return in the TYPERAW branch. `6427b757e`.
- **647060/647114/647076** (command.c) — leaked a parsed SirilWorldCS in `process_eqcrop` and a
  freshly-loaded sequence in `process_seq_ghs`'s COL_SAT error path. `c2dcf18c5`.
- **647085** (photometry.c) — leaked `siril_catalogue` in `catmag_mono_worker` on every exit
  path; added NULL-safe `siril_catalog_free` to each. `bca5acf20`.
- **647083** (healpix.cpp) — leaked HTTP `buffer`/`data_buffer` on two error returns. `cc8e5be8b`.
- **641726/642868** (deconvolution.c) — leaked `ss_string` (log hook) and `seqEntry`
  (sequence command alloc-failure path). `61ab60988`.
- **647062/647066/647094** (message_dialog.c) — `strip_last_ret_char` deref'd its argument before
  the callers' null checks; made it NULL-safe (and fixed an empty-string `str[-1]` OOB). `3a27a281c`.
- **647101** (mpp.cpp) — NULL deref of optional `run->included` in the quality-sort loop; treat
  NULL as "all frames included". `da9c14954`.
- **647088** (healpix.cpp) — unguarded `*error_status = READ_ERROR` though error_status is
  optional. `bf8053e98`.
- **647120** (siril-window.c) — `data->result` written in the `!data` branch. `3f76abab1`.
- **647090** (siril_pythonmodule.c) — `g_strstrip` before the `!stdout_data` check. `be387c338`.
- **647110** (script_menu.c) — dropped an impossible null check on `capitalized` that Coverity
  paired with a later use. `3315d9d18`.
- **647069** (mpp_sidecar.c) — removed two dead version-dispatch branches (version pinned to 9);
  `mpp_sidecar_test` still 4/4. `39d6cc25a`.
- **647091/647108/647097** (siril_pythonmodule.c) — dead `spawn_error->message` fallbacks (g_spawn_sync
  was passed NULL for its GError**), dead `if (python_exe)` cleanup branch, and unchecked
  `g_spawn_sync` return. `958b88d19`.
- **647096/647125** (dialog_preview.c) — unchecked `ser_read_frame` (NULL `pdata` deref) and a
  div-by-zero when the row-accumulation loop breaks immediately (`m` stays 0). `e20c34d5b`.
- **647065** (DA3D.cpp) — `abs(det)` (could bind to integer `::abs`) → `std::fabs` so the
  singularity guard protecting the `/det` division is sound. `96979b302`.
- **647078** (mpp_stack.cpp) — guard the `nsum` division (provably >0, but Coverity can't see it);
  `mpp_stack_test` still 19/19. `668a05306`.
- **647071** (stackminmax.c) — INT_MIN-shift overflow (`x - INT_MIN`); skip invalid-shift frames
  instead of the insufficient `+= 1` nudge. `779e4a036`.
- **530794** (pipp_avi_write.cpp) — `0xFFFFFFFF` into an int32_t field overflowed to -1; compute the
  RIFF byte budget in int64_t. `6d067c9c5`.
- **647058** (applyreg.c) — `(void)` on the intentionally-ignored `read_fits_metadata_from_path_first_HDU`. (committed)
- **647079** (utils.c) — `(void)` on the best-effort `g_file_move` in `delete_directory`. (committed)
- **647092** (undo.c) — `(void)` on best-effort `g_unlink`/`SetFileInformationByHandle`. (committed)
- **647102** (conversion.c) — `(void)` on `g_strlcpy` of a short fixed Bayer string. (committed)
- **647123** (fits_keywords.c) — `(void)` on `g_strlcpy` of a FITS keyword value. (committed)
- **647118** (call_nlbayes.cpp) — `std::move(orig)` into `input` (last use of `orig`). `ad1458300`.

### Re-examined `Triaged`/`Fix Submitted` rows (earlier fix did not fully resolve)

- **641702/641704** (masks.c) — `fits_to_mask` NULL-deref when a float image has BYTE/USHORT
  `orig_bitpix` (data==NULL); `mask_from_lum_hook` called `destroy_mtf_data(data)` on the wrong
  type, leaking `gi_data`. `94d741f90`.
- **641689** (deconvolution.c) — `deconvolution_finalize_hook` deref'd `data` before the
  `else if (data && …)` null check; collapsed the duplicate branches. (committed)
- **530672** (sum.c) — same INT_MIN-overflow as 647071; the `+= 1` nudge was ineffective. Skip the
  frame. (committed)
- **530766** (comparison_stars.c) — `var_stars_cat` leaked when `sort_compstars` returned early;
  free it in `free_compstars_arg` and NULL it after the in-place free. (committed)
- **530688** (siril_pythoncommands.c) — out-of-range index sent STATUS_ERROR then fell through to a
  second STATUS_OK (overwriting the error result); skip the OK path on error. (committed)
- **641563** (opencv.cpp) — `contours = { pts }` copied; `std::move` pts into the vector. (committed)

### Other
- **647127** (dialog_preview.c / .h) — removed the disused GtkFileChooser preview pipeline (no
  callers; the leaking `siril_file_chooser_add_preview` vbox tree, plus `update_preview`,
  `end_update_preview_cb`, `new_preview_object`, `siril_preview_free`, `is_callback_called`, the
  preview structs/globals, and the `siril_get_file_info`/jpeg-probe helpers only it used). Previews
  are handled by the custom file browser's `siril_file_browser_default_preview`. Verified the live
  preview path is correct; live extractors untouched. `28eb456fc`.
- **647121** (mpp_sidecar.c) — range-checked untrusted header fields in the binary sidecar reader:
  `frame_rows`/`frame_cols` in `[1, MPP_SIDECAR_MAX_DIM]`, `num_layers` in `[1, 3]`, `bitdepth`
  restricted to 8/16/32. Valid sidecars unaffected (`mpp_sidecar_test` 4/4). `5e962088d`.
- **647063** (processing.c) — `generic_image_worker` flagged for a possible NULL `orig` at the mask
  blend; the invariant guarantees it's non-NULL, but Coverity can't see `hook_fit == argfit`. Added
  an explicit `else if (orig)` guard to harden the invariant. `e4eb85d21`.

## Investigated — no change needed (false positives / already correct)

- **647080** (file_browser.c `on_new_folder_create`) — flagged leak of `name` at line 958, but
  `name` is freed on all three return paths and the `child` GFile is balanced
  (`apply_current_folder` g_object_ref's it, function unrefs at the end). False positive in
  current source (function last changed 2026-06-13, before the scan).
- **647095** (gui_iface_impl.c `impl_update_sequence_overlay_async`) — `g_thread_join()` consumes
  the reference returned by `g_thread_new()`, so the GThread is not leaked. Coverity does not
  model `g_thread_join`'s ref consumption. False positive.
- **647105** (siril_pythonmodule.c `teardown_connection_and_worker`) — same class as 647095: the
  worker GThread (created once at line ~3329) is consumed exactly once by `g_thread_join` (sync
  path) or by `python_teardown_worker`'s join (async path); `conn` is freed by
  `teardown_connection` on both paths. False positive (g_thread_join ref consumption not modelled).
- **647082/647099** (file_browser.c `reset_browser_state`) — `fb` is null-checked with an early
  `return` at the top, and the line-3085 access is itself inside `if (fb->search_entry)`. The
  flagged deref-after-null-check is the known Coverity pattern where the null test lives in a
  callee (`install_folder_monitor(fb, …)`). False positive.
- **647059** (mpp.cpp:1485) — dereference of `seq`, which is guarded by the early
  `return MPP_EINVAL` at the function top. False positive (the real `run->included` deref was
  647101, fixed above).

---

## Concurrency cluster (11 CIDs) — resolved

Originally flagged for a locking-model decision; dispositioned with your guidance.

**(a) False positive / intentional — mark as such in Coverity** (no code change)
- **647064** (processing_thread.c `processing_system_shutdown:285`) — the unlocked read of
  `worker_running` is *required*: `queue_mutex` is only initialised after `system_initialized`,
  so locking before the early-out would touch an uninitialised mutex if shutdown runs before init.
  Only the main thread writes it; the worker only reads it. No real race.
- **647112** (stacking.c `init_stacking_args`) — field initialisation of a freshly-created args
  struct before any worker thread exists.
- **647067** (geometry.c `rotation_image_hook`) — per-frame workers only *read* the shared
  `rotation_args`; read-only sharing, no writer.

**(b) Won't-fix — benign lock-free `com.pref.*` config access** (per your call)
- **647075 / 647103 / 647111** (callbacks.c `initialize_preprocessing`), **647126**
  (OS_utils.c `get_max_memory_in_MB`), **647124** (preferences.c `update_preferences_from_model`).
  Siril reads `com.pref.*` lock-free in hundreds of places; these are aligned scalar reads/writes
  whose worst case is a momentarily-stale config value. A real fix would be a cross-cutting
  `com.pref` locking/snapshot strategy — explicitly deferred.

**(c) `com.stars` locking — FIXED** (commit `5a8e4ad90`)
- **647113 / 647116** (deconvolution.c). Did a full audit of every `com.stars` /
  `com.star_is_seqdata` access against the lock contract. The two flagged deconvolution sites are
  already correctly locked in the current tree (the scan predates the wrapping). The audit found
  the real residual gaps — all the same root cause: snapshot the `com.stars` *pointer* under the
  reader lock, release, then dereference the array on a worker thread (use-after-free if a
  concurrent clear/replace frees it). Fixed by:
  - `clear_stars_list()` now reads/resets `com.star_is_seqdata` inside the writer section (and
    closed a latent empty-array leak);
  - new `snapshot_com_stars()` returns a private deep copy under the reader lock (via
    `duplicate_psf`), used by `mask_create_from_stars`, `generate_synthstars`,
    `generate_binary_starmask`, and the Python `CMD_GET_PSFSTARS` handler.
  - Follow-up (commit `d2b0ba32e`): fixed the deconvolution `com.stars = detected` overwrite. It
    leaked any prior list, split the check from the assignment across two lock acquisitions
    (TOCTOU), and on a 0-star (empty non-NULL) detection left the array in `com.stars` unfreed.
    Added `replace_com_stars()` (atomic free-old + install-new under the writer lock) and now check
    the detection result before installing. The other two direct `com.stars =` sites
    (`update_stars_idle` clears first, `registration.c` only assigns when NULL) were verified clean.

**(d) Non-atomic global mutation — FIXED** (commit `2fbe7fcbe`)
- **647100** (siril_pythoncommands.c). The Python load command no longer toggles the global
  `com.pref.debayer.open_debayer` around `read_single_image`; it passes an explicit `no_debayer`
  argument instead (other callers pass `FALSE`, unchanged behaviour).

---

# Second pass — `Thorough.csv` (re-scan with dismissed / FP / intentional filtered in)

This scan (90 CIDs) re-included previously-dismissed/FP/intentional issues plus 4 brand-new
(06/27) CIDs. Summary of this pass:

- **14 CIDs fixed** across 8 commits: 4 regressions I introduced + 4 genuine leaks wrongly
  dismissed + 6 (638293/511593/515929/358911/515030/358842) actioned from the deferred list at your
  direction — plus 1 real use-after-free found in passing. All build clean, 21/21 tests pass.
- **5 of my first-pass fixes** re-appeared only because the scan baseline predates this branch —
  verified all still present in the code (already fixed).
- **Nothing left deferred** — every actionable CID in this scan is fixed, already-fixed, or
  re-confirmed correctly dismissed (FP / intentional / vendored).

## Fixed this pass

### Regressions from my first-pass `snapshot_com_stars()` refactor (NEW, 06/27)
The refactor made the borrowed `com.stars` an *owned* deep copy, so every early exit after the
snapshot must free it. Commit `c75a21dcf`.
- **647137** (unpurple.c `generate_binary_starmask`) — **real leak**: the no-stars and
  `new_fit_image`-failure returns are reachable with a non-NULL owned `stars`. Now freed.
  **Correction (dip-check vs the Coverity trace):** the trace's `overwrite_var` step points
  at line 92 `stars = peaker(...)`, *not* the early returns. `snapshot_com_stars` can return a
  non-NULL but empty array (first `duplicate_psf` OOM → `copy` allocated, `*nb_out == 0`), which
  enters the `nb_stars < 1` branch and overwrites `stars` without freeing the snapshot. The
  early-return frees in `c75a21dcf` did not cover this; added a free of the snapshot before the
  `peaker()` reassignment. Commit `8e3eb6eac`.
- **647139 / 647138 / 647140** (masks.c / synthstar.c / siril_pythoncommands.c) — the flagged
  `sf_data`-alloc-failure exits sit in the "no com.stars" branch where `stars` is NULL in practice,
  but the owned-resource contract should hold on every path; free it there too (NULL-safe).

#### Audit follow-up — same overwrite leak in the other three snapshot callers (not Coverity CIDs)
Triggered by the 647137 dip-check: the same non-NULL-empty return (`snapshot_com_stars` first
`duplicate_psf` OOM) reaches `findstar_worker`, whose `*args->stars = stars` (star_finder.c:1455)
overwrites the borrowed pointer **without freeing it**. The three callers that pass `&stars` into
`findstar_worker` (`generate_synthstar` synthstar.c, `mask_create_from_stars` masks.c, the
`CMD_GET_PSFSTARS` handler in siril_pythoncommands.c) therefore leak the snapshot array on that
OOM edge. Fixed at the call sites (free the snapshot at the top of the `nb_stars < 1` branch,
mirroring the unpurple fix) rather than in `findstar_worker` — it's an output-param producer whose
other callers pass uninitialised/borrowed pointers. Coverity did not flag these (it never traced
the non-NULL-empty return into them). Commit `db0daa1c7`.

### Genuine leaks wrongly dismissed — now fixed (commit `c63c7bd3b`, `862dea01d`)
- **509293 / 509291** (extraction.c) — `extractHaOIII_image_hook` / `split_cfa_image_hook` freed
  only the per-channel fits on the error path, leaking the `images[]` array + `multi_data` struct
  **every failed frame** (per-frame sequence hooks). `c63c7bd3b`.
- **641576** (astrometry_solver.c `add_object_in_tree_view`) — leaked `args` on the "No object
  found" early return. `c63c7bd3b`.
- **452432** (cut.c `on_cut_sequence_apply_from_gui`) — leaked `arg` (+3 strdup'd strings) when the
  cut struct was invalid; freed via `free_cut_args`. `c63c7bd3b`.
- **433622 / 494308** (siril_update.c) — `siril_check_updates` / `siril_check_notifications`
  discarded the `g_thread_new()` handle; now `g_thread_unref`'d (matching the SPCC sibling). `c63c7bd3b`.
- **583967** (PSF_list.c `update_star_list`) — leaked `args` on the headless path; freed (matches
  the no-op gui_iface stub's no-ownership contract). `c63c7bd3b`.
- **530799** (siril_pythonmodule.c `execute_python_script`) — `script_name` (caller transfers
  ownership) leaked on the success path; freed after the `temp_filename` copy. `862dea01d`.
- **637631 / 583996** (siril_pythonmodule.c `execute_python_script`) — `stdout_data`/`stderr_data`
  GDataInputStreams leaked their creation ref (cleanup unref'd the underlying stream, not the data
  stream), leaking the object + fd per run. Drop the creation ref. `862dea01d`.

### Real bug found while reviewing (bonus, not a listed CID) — `c63c7bd3b`
- **python_gui.c `on_action_file_execute`** (LANG_SSF path) — use-after-free: `text` was handed to
  a worker-thread input stream with a NULL destroy-notify and then `g_free`'d immediately. Gave the
  stream ownership (`g_free` destroy-notify); `execute_script` unrefs the stream when done.

## Already fixed in this branch (first pass) — verified present, no action
- **647121** (mpp_sidecar range checks), **647109** (background_extraction clamp), **647087**
  (livestacking), **647076** (process_seq_ghs), **647065** (DA3D `std::fabs`), **530766**
  (comparison_stars `var_stars_cat`). They appear only because the scan baseline predates the branch.
- All the 06/26 "New" rows (647126/647124/647112/647111/647103/647075/647105/647099/647095/647082/
  647080/647067/647064/647059) are the FP / won't-fix / intentional items already documented above
  (concurrency groups a/b, the g_thread_join FPs, the file_browser/mpp FPs).

## Deferred items — now actioned per your direction

1. **638293 + 511593** (siril_pythonmodule.c) — **FIXED** `ff54e99a0`. Cap untrusted width/height to
   100000 and compute the pixel-count / pdata-offset products in 64-bit (`size_t`) so they can't wrap.
2. **515929** (siril_pythoncommands.c `unpack_plot_data`) — **FIXED** `6698467bb`. Bounds-check every
   read against `buffer_size` (overflow-safe `NEED()`), copy strings only up to a NUL found inside the
   buffer, verify all point bytes before allocating, lenient `num_series <= buffer_size` cap, single
   cleanup path.
3. **358911** (image_formats_internal.c `import_pnm_to_fits`) — **FIXED** `0c319f254`. Reject zero and
   cap PNM dimensions to 100000/axis before any allocation.
4. **515030** (pixelmath.c `pixel_math_evaluate`) — **FIXED** `d860e7cec`. `handed_off` flag +
   full cleanup epilogue (free_pm_var + varname + expressions + fit + args) on every error path,
   skipped only once the worker takes ownership (no double-free on success).

### 5. 358842 (seqwriter.c `seqwriter_set_max_active_blocks`) — FIXED `7231ba33c` (option a)

Per your direction, documented the "called once at single-threaded setup" contract (which is why
the unlocked writes are safe) and deleted the dead, racy, buggy up-scaling branch. Details below
for reference.

**What the function does.** It sets `configured_max_active_blocks` (the cap on how many decoded
image blocks may sit in the write queue at once) and resets `nb_blocks_active` (how many are queued
right now). Both are shared with the writer pool: `seqwriter_wait_for_memory()` /
`seqwriter_release_memory()` / `notify_data_freed()` read-modify them **under `pool_mutex`**;
`seqwriter_set_max_active_blocks()` writes them **without** the lock — that inconsistency is what
Coverity flags.

**What "dynamic up-scaling" means.** Lines 268-273: if the function is called a *second* time with a
*larger* `max` while a cap is already in force (`configured_max_active_blocks > 0`), it computes
`more = max - old` and calls `seqwriter_release_memory()` that many times. The intent is to *raise
the queue limit mid-run* and immediately wake that many blocked producer threads so they can use the
newly-granted slots — i.e. grow the memory pool on the fly. That's the only path that would run
concurrently with active workers, so it's the only genuinely racy one.

**Why it's benign today.** Every caller (`seq_prepare_hook` → processing.c:590, sequence_export,
conversion, command.c) calls it **exactly once, during the single-threaded prepare phase**
(generic_sequence_worker runs `prepare_hook` at processing.c:140, *before* the parallel frame loop at
line 210 and before producers call `seqwriter_wait_for_memory`). So nothing reads the fields while
it writes them, and the up-scaling branch is never taken. There's also a latent logic oddity: even
when the branch runs, line 276 unconditionally sets `nb_blocks_active = 0` afterwards, which would
discard both the `more` decrements *and* any genuinely in-flight count.

**Options / recommendation.**
- (a) **Document the contract** (call only once, at setup, before the writer runs) and leave the code
  — lowest risk, reflects reality. I'd pair it with deleting the dead up-scaling branch (it's unused
  and its interaction with the `= 0` reset is buggy).
- (b) **Make it lock-safe** if you intend dynamic scaling to be a real feature: take `pool_mutex`
  around the updates, and inline the slot releases (`nb_blocks_active--; g_cond_signal`) instead of
  calling `seqwriter_release_memory()` (which re-takes `pool_mutex` → self-deadlock), and drop the
  unconditional `nb_blocks_active = 0`.

My lean is **(a)** — there's no current need for runtime up-scaling. Tell me which way and I'll apply it.

## Re-confirmed correctly dismissed (no action)

- **Already correctly handled / FP** (ownership passes to a worker, freed in finalize/idle, value
  can't be NULL/untrusted, etc.): 641698, 641677, 641629, 641594, 433620, 433619, 487611 (taint
  bounded), 358748, 358903, 358840, 358839, 358805, 530703, 532536, 583995, 494309, 433625, 358923,
  358772, 583981, 530745, 531973, 515199, 637615/583962/583959 (process_pyscript — ownership passed
  to the worker thread).
- **Intentional / acceptable**: 455007, 466086, 358965, 358855, 358879 (best-effort `g_unlink` /
  safe `g_strlcpy` truncation), 494095 (non-OpenMP default), 358958 (reachable abort check), 488375
  (cond-wait, not infinite), 358937 (deliberate in-place expansion), 358930 (mkfifo EEXIST handled),
  499074 (degraded RNG-seed fallback).
- **Vendored third-party** (don't patch): 462101, 462099 (Little-CMS `lcms_acceleration`); and the
  `/usr/include` C++/OpenCV headers (637647, 585112, 432129, 583310, 583253, 530034, 453412, 352844,
  346185, 126278, 75640).

---

# Bonus fix found during headless testing (not a Coverity CID)

While exercising the fixes via siril-cli, `stack <seq> sum|min|max` **segfaulted** on a
heterogeneous/unregistered sequence. Root cause: `sum_stacking_image_hook` /
`minmax_stacking_image_hook` read `fit->fpdata[layer]` (when `input_32bits`) or
`fit->pdata[layer]` for every output layer; a frame whose bit depth differs from the accumulator
(e.g. a 32-bit float frame in a 16-bit stack) or with fewer channels than the output leaves that
buffer NULL → NULL deref. Added an up-front compatibility check to both hooks that aborts
(`ST_SEQUENCE_ERROR`) with a clear message naming the offending frame, instead of crashing.
Verified: the offending sequence now aborts cleanly ("Image #4 … is incompatible …"), a valid
homogeneous sequence still stacks. Commit `98e4ab09e`. (Pre-existing bug, independent of the
Coverity work; my 530672/647071 changes only touched the unrelated INT_MIN branch.)

---

# Third pass — verification against saved Coverity traces (`123456.txt` in project root)

The user exported the detailed trace for each issue still showing open. Cross-checking the **actual
current code** against each trace (not trusting the earlier log) — prompted by the 647137 miss,
where my "fix" addressed a different defect than the trace flagged. Several traces here likewise
describe a *different* defect than my earlier log attributed to that CID number (Coverity appears to
have re-pointed some CIDs between scans), so each was re-derived from source.

**Two real bugs found and fixed (647087 leak, 647059 null-deref); the other 12 confirmed genuinely
FP / intentional / residual.** (647059 was initially mis-called an FP here and corrected after the
user supplied its full transitive trace — see its entry below.)

## Real leak the earlier fix missed — now fixed

- **647087** (livestacking) — the trace flags `str` at livestacking.c:223 (`livestacking_display(str,
  TRUE)` then `return`), **not** the `replace_ext` path my `6427b757e` touched. `livestacking_display`'s
  second arg `free_after_display=TRUE` means the callee owns `str`, and the interactive GUI path frees
  it (via `set_label` → `label_update_idle`). But **both headless paths break the contract**:
  `livestacking_gui.c` logs and `return`s without freeing when `com.headless`, and the
  `headless_stubs.c` stub ignored `str` entirely. Every `TRUE` caller (livestacking.c:224, 774, 806)
  therefore leaks in headless mode. Fixed at the source — both implementations now `g_free(str)` when
  `free_after_display`. Compiles clean (gui + headless stub objects). Commit `46d9724a2`.

## Confirmed false positives (verified against current source)

- **530766** (compstars_worker, `query_args->item`) & **641698** (process_catsearch, `query_args->name`)
  — `free_sky_object_query` (siril_catalogues.c:1636) **does** free `->name`, `->item` (and the item's
  contents). Coverity doesn't model the field frees. FP. (My separate 530766 `var_stars_cat` fix is a
  different leak and stands.)
- **530799 / 637631** (execute_python_script) — `python_process_cleanup` frees `cleanup->temp_filename`
  (3255) and `cleanup` (3256); `script_name` is freed at 3533; the async path hands `cleanup` to the
  child-watch. Nothing leaks at the flagged points — Coverity can't trace the callee frees / async
  ownership handoff. FP (residual after `862dea01d`).
- **647059** (mpp_stack_apply_impl) — **REAL null-deref, not an FP** (corrected after the user supplied
  the full trace). `mpp_classify_sequence_input(seq)` at mpp.cpp:1485 dereferences `seq->imgparam`
  *transitively*, three calls deep: → `fits_seq_is_cfa` → `mpp_seq_read_frame` (SEQ_REGULAR) →
  `fit_sequence_get_image_filename_checkext`, which reads `seq->imgparam[index].filenum`
  (catalogues path, line 1375). The function itself treats `seq->imgparam` as possibly-NULL at line
  1397, so the null state is reachable; the classify call would then crash. Every downstream consumer
  needs `imgparam` (all read frames by filename), so a NULL can't yield a valid stack — added
  `!seq->imgparam` to the precondition guard (fail fast with `MPP_EINVAL`) and dropped the now-redundant
  `seq->imgparam &&` at 1397. My earlier "classify never derefs imgparam" was wrong — I stopped at
  `fits_seq_is_cfa`/`mpp_seq_read_frame` and didn't follow into the filename helper. Commit `<pending>`.
- **647076** (process_seq_ghs) — the trace flags the **success** path, but `apply_ght_to_sequence`
  transfers `seqdata` to the async worker via `args->user`, freed in `ght_finalize_hook` (frees `data`
  + `params`). Ownership handoff Coverity can't see. FP. (My `c2dcf18c5` COL_SAT error-path fix stands.)
- **647080** (on_new_folder_create) — `name` is freed on all three returns (file_browser.c:960/968/986). FP.
- **647082 / 647099** (reset_browser_state → update_nav_sensitivity) — `back_history`/`forward_history`
  are allocated **unconditionally** in the constructor (file_browser.c:2497-2498) and never nulled;
  `update_nav_sensitivity` derefs them safely. The defensive `if (fb->…history)` guards in
  reset_browser_state are what mislead Coverity. FP. *(Earlier "null-test-in-callee" reasoning was
  imprecise — corrected.)*
- **647109** (generate_samples_random) — `tmp` is built mono (1 channel, bg_extraction.c:701), and
  `find_unlinked_midtones_balance` only writes `results[0..naxes[2]-1]` = `results[0]`, so the singleton
  `&params` is correct. FP (ARRAY_VS_SINGLETON heuristic). *(Unrelated to the bbox-clamp `13252648b`,
  which the CID number formerly pointed at.)*

## Confirmed intentional / residual (no change)

- **647064** (processing_system_shutdown) & **647067** (rotation_image_hook) — LOCK_EVASION. 647064's
  unlocked `worker_running` read is required (queue_mutex is initialised after `system_initialized`).
  647067's unlocked `com.selection.w` check-then-write is benign: parallel per-frame rotate hooks all
  expand the (idempotent) global selection to the full image under `com.mutex`; worst case is a
  redundant identical write. *(647067's earlier "workers only read rotation_args" reasoning was wrong —
  the flagged field is the global `com.selection`; corrected.)*
- **647065** (ComputeRegressionPlane) — the `std::fabs(det) < epsilon` guard from `96979b302` is
  present; the `/det` in the `else` runs only when `|det| ≥ epsilon`, so it can't divide by zero.
  Residual DIVIDE_BY_ZERO FP.
- **647121** (mpp_sidecar_read) — the range checks from `5e962088d` are present; the shifts-section
  `snf`/`sna` are rejected unless they equal the already-validated `num_frames`/`aps_count`, so
  `success_count = snf*sna` is bounded. Residual TAINTED_SCALAR FP (Coverity doesn't propagate the
  equality). Could add an explicit cap purely to silence it, but no real defect.
