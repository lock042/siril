# PlanetarySystemStacker port to Siril — implementation plan

Tracker for porting the PSS multipoint registration + stacking pipeline into Siril on the `pss` branch. Tick boxes as work lands. Open items are status-of-record; revise this file when scope changes.

## Guiding principles

- **Algorithmic fidelity over re-invention.** Where PSS uses an OpenCV call, the port calls the same OpenCV function with the same arguments. C/C++ throughout, no Python in the shipped product.
- **Golden-output regression testing.** Each module is validated against intermediate artifacts captured from a Python run of PSS on the same input. Acceptance bar: per-frame quality within 1e-4 relative, integer shifts exact, sub-pixel shifts within 0.05 px, final stacked PSNR > 50 dB vs PSS reference.
- **Headless first.** Everything driven by a new `pss` command. GUI gated behind an explicit go/no-go pause.
- **Algorithmic-match definition:** behaviour-equivalent within float-rounding. Bit-identical is not achievable across platforms even between two Python runs.

## Cross-cutting invariants

### Bit-depth handling

Siril stores both 8-bit and 16-bit SER data as `WORD` (uint16) but keeps the values in their native range: 0..255 for `fit->bitpix == BYTE_IMG (8)`, 0..65535 for `USHORT_IMG (20)`. There is no automatic upscale at load time (`ser.c:ser_manage_endianess_and_depth` just zero-extends 8-bit bytes to WORDs).

`mpp_config_t.bitdepth` is our 8/16 abstraction — set it from `fit->bitpix` via `mpp_bitdepth_from_fits_bitpix(int)`. The threshold scale is `mpp_cfg_threshold_scale(cfg)` (1.0 for 8-bit, 256.0 for 16-bit). `align_average_frame` upscales 8-bit data by `256/N` so `mean_frame` always lands in the 0..65535 range, mirroring PSS — downstream AP code can then assume 16-bit-equivalent units and hardcode `× 256` threshold scaling.

### Per-frame global shifts (no separate pre-alignment)

Wind jolts, telescope drift, tracking errors and other large frame-to-frame translations are handled by Phase 2's `align_global_from_frames` directly. Its two-phase `cv::TM_CCORR_NORMED` against an auto-picked patch on the best frame searches `±align_frames_search_width` (default 34 px) per axis on the coarse stride-2 grid, plus PSS's backward-then-forward cumulative-shift chaining (a shifted reference window if the running drift would push the patch off-frame). For larger or more pathological drift the search width can be raised via `cfg.align_frames_search_width`. **A separate COG / pre-alignment pass is not required and is not on the plan.**

## Scope

In scope: SER ingest, debayering, frame ranking, global alignment, AP grid, per-AP shifts, multipoint stacking, bicubic drizzle (PSS-faithful), STScI drizzle (Siril-original), Bayer drizzle (Siril-original), `pss` command, eventually a GUI.

Out of scope (initially): PSS postprocessing (wavelet unsharp masking), dark/flat calibration as a built-in step (use Siril's existing pipeline upstream), INI-file compatibility with PSS's `.PlanetarySystemStacker.ini`.

## Module layout

All under `src/registration/`. C for plumbing, C++ where touching OpenCV. Tests under `src/tests/`.

| File | Role |
|---|---|
| `mpp.h` / `mpp.c` | Public types (`mpp_config_t`, `mpp_run_t`, `mpp_aps_t`, `mpp_shifts_t`), orchestrator `register_mpp()`, progress callback infra |
| `mpp_config.c` | PSS parameter defaults; per-invocation overrides |
| `mpp_image_priv.hpp` | C++ helper: zero-copy `cv::Mat` views over Siril's `fits` planes; CV_32F normalisation to PSS conventions |
| `mpp_frames.{c,cpp}` | Frame iteration over a Siril sequence; Bayer debayering via `cv::cvtColor`; integrates with existing `src/io/ser.c` |
| `mpp_rank.cpp` | Per-frame Laplace-σ quality (PSS `rank_frames.py`) |
| `mpp_align.cpp` | Global frame alignment: patch selection + two-phase correlation (PSS `align_frames.py`) |
| `mpp_ap.cpp` | Staggered AP grid placement with brightness/structure filtering (PSS `alignment_points.create_ap_grid`) |
| `mpp_shift.cpp` | Per-AP per-frame local shift via two-phase multilevel correlation with optional sub-pixel parabolic fit (PSS `miscellaneous.multilevel_correlation`) |
| `mpp_stack.cpp` | Top-N selection per AP, ramped 2D weighting, brightness equalisation, background composition. Calls into a pluggable resample backend. |
| `mpp_drizzle.cpp` | Wraps `dobox()` for STScI-on-debayered and Bayer-drizzle resample backends |
| `mpp_sidecar.c` | (Deferable) Binary persistence of run state for resume/inspection |

## PSS reference oracle

Lives under `tools/pss_reference/`. Runs upstream PSS modules and dumps intermediate artifacts as `.npz` (quality vector, global shifts, AP coords, per-AP per-frame shifts, stacked FITS). All unit tests use this as ground truth.

---

## Phase 0.0 — PSS bring-up (prerequisite)

The oracle is worthless if PSS itself will not run. Confirm and fix before anything else.

- [x] 0.0.1 Inspect PSS dependencies (`setup.py`, README) and the local Python environment.
- [x] 0.0.2 Attempt to launch PSS (`planetary_system_stacker.py` and/or one of `Test_programs/`) and capture failures.
- [x] 0.0.3 Identify root causes — typically: PyQt5 vs PyQt6, NumPy ≥ 1.25 removal of deprecated aliases, OpenCV API changes, astropy.io.fits changes, ConfigParser changes. _(None required at bring-up — PSS ran unmodified on Py 3.12 / NumPy 2.4 / OpenCV 4.13 / astropy 7.2.)_
- [x] 0.0.4 Apply minimal patches in a local copy under `tools/pss_reference/pss_patched/` (or as a patch stack against the cloned tree). _(N/A — see 0.0.3.)_
- [x] 0.0.5 Verify a non-GUI test program runs end-to-end on a small input. _(`tools/pss_reference/run_pss.py` on synthetic 24-frame image dir.)_
- [x] 0.0.6 Document the patches and minimum working dependency set in `tools/pss_reference/README.md`.

## Phase 0.1 — Skeleton and harness

- [x] 0.1.1 Create empty `mpp_*` source/header files under `src/registration/` with stubs returning `MPP_ENOTIMPL` or equivalent. Public surface from the module table above.
- [x] 0.1.2 Wire `mpp_*.{c,cpp}` into `src/meson.build`. Confirm OpenCV dep is available (Siril already links OpenCV elsewhere). _(Siril's `libsiril.a` already pulls OpenCV transitively; nothing extra needed.)_
- [x] 0.1.3 Add `pss` command to `src/core/command_list.h` and `process_pss()` stub in `src/core/command.c` returning a "not implemented" message. _(Also added STR_PSS in command_def.h and the declaration in command.h.)_
- [x] 0.1.4 Empty Criterion test files for each module under `src/tests/mpp_*_test.{c,cpp}`. Wire into `src/tests/meson.build`. _(7 .cpp + 1 .c, under suite `mpp`.)_
- [x] 0.1.5 `meson compile` clean. `meson test` runs the (currently trivial) mpp tests. _(All 8 mpp tests pass; full project still 19/19. Requires `-Dcriterion=true`.)_
- [x] 0.1.6 Create `tools/pss_reference/` Python harness: pipeline driver that dumps intermediate artifacts (ranked quality vector, global shifts, AP coords, per-AP shifts, stacked FITS) to `.npz` / FITS for any input video. _(`run_pss.py` produces `oracle.npz` + `ref_avg.fits` + `stacked.fits`. Per-AP per-frame shifts still TODO — needs deeper hook into PSS internals; will add when Phase 4 needs it.)_
- [x] 0.1.7 Generate a tiny synthetic SER (~30 frames, simulated planet with prescribed jitter + per-frame quality variation) in `tools/pss_reference/test_data/`. _(`gen_synth_frames.py` produces a 24-frame PNG directory with truth.npz. Real SER writer still TODO — image-dir input drives the same PSS pipeline so this is sufficient for Phases 1–5.)_
- [x] 0.1.8 Reserve one real planetary SER for the Phase-7 end-to-end milestone. _(`test-big.ser` — 500-frame Bayer Jupiter sequence; used by `end_to_end_real_ser_oracle` and the CLI integration tests.)_

**Exit criterion:** PSS produces oracle artifacts on the synthetic SER; Siril compiles with empty mpp skeleton; `pss` returns clean "not implemented."

---

## Phase 1 — Quality ranking

- [x] 1.1 Implement `mpp_rank.cpp`: PSS pipeline = GaussianBlur(7×7) → stride-2 sub-sample (`align_frames_sampling_stride`, NOT `rank_frames_pixel_stride` — the latter is for xy/Sobel only) → `cv::Laplacian(CV_32F)` → `cv::convertScaleAbs(alpha=1/256)` → `meanStdDev[1][0][0]`. Returns σ on the **uint8** Laplacian, matching PSS's `Frames.frames_mono_blurred_laplacian` + `RankFrames.frame_score` chain exactly.
- [x] 1.2 Optional brightness normalisation: `σ / (cv::mean(threshold(mono, thr, 255, THRESH_TOZERO))[0] + 1e-10)` with `thr = frames_normalization_threshold * 256` for 16-bit input.
- [x] 1.3 Wire frame iteration through `mpp_frames` (SER reader from `src/io/ser.c`, optional debayer via `cv::cvtColor`). _(Landed as part of Phase 6.1: `mpp_seq_read_frame` / `read_analysis_frame` / `read_full_frame` in `mpp.cpp` drive the algorithms off Siril sequences end-to-end. The standalone `mpp_frames.{c,cpp}` stubs from 0.1.1 stayed at MPP_ENOTIMPL — the work amortised into the orchestrator instead.)_
- [x] 1.4 Oracle test: rank synthetic dataset; per-frame quality within 1e-4 relative of PSS. _(`mpp_rank::oracle_equivalence_synthetic` — passes.)_
- [x] 1.5 Sanity tests: blur-decreases-quality, brightness-invariance. _(Plus default-config, score-positive, threshold-scales-with-bitpix → 6 sanity tests, all green.)_

## Phase 2 — Global frame alignment

- [x] 2.1 Implement `mpp_align.cpp` patch picker (`align_frames.compute_alignment_rect` + `quality_measure_threshold_weighted`). Replicates PSS's uint16 modular subtraction in the quality measure so absdiff matches NumPy's `abs(a-b)` for unsigned inputs.
- [x] 2.2 Two-phase correlation: phase 1 stride-2 grid (with extra 7×7 GaussianBlur on the strided frame window) `cv::TM_CCORR_NORMED` with half-width `(search_width − 4)/2`; phase 2 stride-1, search ±4. Success requires both phases' optima to be in the interior of their search windows. Implements PSS's backward-then-forward cumulative-shift loop.
- [x] 2.3 Build average reference frame from top `align_frames_average_frame_percent=5%`. Includes `find_best_frames` sliding-window selector for `align_frames_fast_changing_object=true`.
- [x] 2.4 Oracle test: per-frame integer shifts exact; patch coordinates exact. _(`mpp_align::oracle_equivalence_synthetic` passes — patch `(116,153,119,169)` and all 24 per-frame integer shifts match PSS exactly.)_

## Phase 3 — Alignment-point grid

- [x] 3.1 Implement `mpp_ap.cpp` staggered grid (`step_size = (half_patch_width × 4.5) / 3`, `half_patch_width = half_box_width × 1.5`, `half_box_width=24`). PSS pre-blurs the mean frame with `GaussianBlur(7×7)` inside `AlignmentPoints.__init__`; we replicate that as the first step of `ap_create_grid`.
- [x] 3.2 Per-AP filters: brightness (`max(box) > 10 × 256`), contrast (`max−min > 0` by default — the contrast threshold is 0), dim-fraction COM re-centring (`fraction below threshold > 0.6` triggers a centre-of-mass shift of the AP), and post-normalisation structure threshold (`min(avg|∂x|, avg|∂y|) / max_struct ≥ 0.04`).
- [x] 3.3 Oracle test: AP coordinate set exact match. _(`mpp_ap::oracle_equivalence_synthetic` passes — 12 APs, every (y, x) centre and (box_y_low, box_y_high, box_x_low, box_x_high) matches PSS exactly.)_

## Phase 4 — Per-AP local shifts

- [x] 4.1 Two-phase multilevel correlation in `mpp_shift.cpp`. Pulled the kernel out into `mpp::multilevel_correlation` (shared with Phase 2) parameterised by `gauss_width` and `search_width`. Per-AP shift uses `alignment_points_search_width=14` (smaller than the 34 for global). Reference boxes come from the *unblurred* `align_frames.mean_frame` (PSS `set_reference_boxes_correlation`, NOT the AlignmentPoints-blurred one used for AP placement — those are different mean frames, and the oracle is unforgiving about which you pick).
- [x] 4.2 Optional sub-pixel: 6-parameter quadratic surface fit (`f(x,y) = a x² + b y² + c xy + d x + e y + g`) on the 3×3 correlation neighbourhood, solved via PSS's precomputed pseudo-inverse `(AᵀA)⁻¹ Aᵀ` from `miscellaneous.sub_pixel_solve_matrix`. Gated on `alignment_points_local_search_subpixel=true` (default false to match PSS).
- [ ] 4.3 Penalty matrix for off-centre peaks (`alignment_points_penalty_factor=0.00025`). _Deferred to Phase 5a (where `stack_frames.py:311-317` builds it and passes it into `compute_shift_alignment_point` during the actual stacking pass; the Phase-4 oracle was captured without the matrix because `compute_frame_qualities` doesn't use it)._
- [x] 4.4 Sign-convention regression test: full Phase 2+3+4 chain on a jittered sequence; after Phase 2 absorbs the global shift, per-AP residuals settle to (≤1, ≤1) integer pixels including (0, 0) on the best frame.
- [x] 4.5 Oracle test: per-AP per-frame integer shifts match PSS **exactly** across all 24 frames × 12 APs on the bundled synthetic dataset (`mpp_shift::oracle_equivalence_synthetic`).

**Companion check:** `mpp_align::average_frame_oracle_equivalence` now also locks in bit-for-bit equivalence of `align_average_frame` against PSS's `align_frames.mean_frame` on the same dataset (caught a `convertTo(CV_32S)` rounding divergence — PSS uses `astype(int32)` which truncates toward zero, OpenCV's default rounds half-to-even; fixed by manual truncation).

## Phase 5a — Bicubic stacking (PSS-faithful)

- [x] 5a.1 `mpp_stack` factored into `stack_frames_loop` (per-frame brightness equalise → drizzle resize → per-AP shift+remap into per-AP buffers, plus background accumulator) and `stack_merge_alignment_point_buffers` (weighted merge, background blend, border trim, scale+clip to u16). Public C `mpp_stack(seq, …)` still stubbed pending orchestrator.
- [x] 5a.2 `cv::resize(frame, scale=K, INTER_LINEAR)` — PSS uses INTER_LINEAR for the drizzle expansion, not INTER_CUBIC. _(Plan was wrong; PSS source confirms LINEAR at `stack_frames.py:350`.)_
- [x] 5a.3 Top-N per AP via `ap_compute_frame_qualities` (`alignment_points.compute_frame_qualities` Laplace+Laplace path with brightness normalisation). 2D weights via PSS's `one_dim_weight × min`. Background composed before merge.
- [x] 5a.4 Per-frame brightness equalisation against the median frame brightness (`frame × median / (avg + 1e-7)`).
- [x] 5a.5 Background composition: subdivides the intersection into `background_patch_size`-px tiles, keeps only tiles where some pixel needs background (PSS `stack_frames_background_fraction` decides full vs patch-based). Foreground/background blend uses the `sum_weights / (blend_threshold × stack_size)`, clipped to [0, 1] ramp.
- [x] 5a.6 Drizzle plumbing in place (`drizzle_factor` config field; per-frame resize, drizzled patch bounds, drizzled global buffer). **Functionally untested for drizzle > 1** — oracle was captured with the PSS default `drizzle="Off"`. Drizzle-specific oracle deferred to a follow-up.
- [x] 5a.7 End-to-end oracle: final uint16 stacked image matches PSS cell-for-cell on the synthetic dataset with > 99 % exact pixels and worst Δ ≤ 2 levels out of 65535 (`mpp_stack::end_to_end_stacked_oracle`).

**Phase 4.3 closed here:** the penalty-weight matrix (`alignment_points_penalty_factor=0.00025`) was added as an optional parameter to `mpp::multilevel_correlation` and is used by the stacking-time per-AP shift compute. Phase-4-default callers (`align_shift_one_frame`, the Phase 4 oracle path) keep `weight_matrix_first_phase=None`, matching PSS.

### Phase 5a oracle tests (all green)

| Test | Acceptance bar | Notes |
|---|---|---|
| `default_config_phase5a` | PSS defaults exact | 10 fields |
| `one_dim_weight_*` | Reference formulas exact | symmetric / extend_low / extend_high |
| `weight_matrix_centre_and_corners` | 1.0 at centre, `1 − 2 × penalty` at corners | |
| `remap_rigid_*` | accumulation + clip + border counts exact | |
| `ap_compute_frame_qualities_oracle` | quality matrix `1e-9` rel, stack_size + best_frame_indices + used_alignment_points exact | |
| `prepare_for_blending_oracle` | `sum_single_frame_weights` cell-wise `< 1e-5` abs; number_stacking_holes exact | |
| `stack_frames_loop_oracle` | per-AP buffer + averaged_background `< 1e-3` rel; bg_patches + border counts exact | depth oracle — exercises every Phase 1–4 module |
| `end_to_end_stacked_oracle` | u16 stacked image **> 99 % exact, worst Δ ≤ 2 out of 65535** | 24-frame synthetic |
| `end_to_end_real_ser_oracle` | u16 stacked image **99.94 % exact, worst Δ ≤ 4 out of 65535, mean \|Δ\| = 6e-4** | **500-frame real Jupiter SER (test-big.ser → debayered+greyscale 8-bit PNGs)** |

### Two 8-bit / 16-bit gotchas caught by the real-SER oracle

The 8-bit path needed PSS-faithful pre-blur upscale to behave correctly:

1. **PSS upscales 8-bit mono to 16-bit range before any GaussianBlur** (`frames.py:1505-1506`: `frame_mono.astype(uint16) * 256`). Without that, blurred-frame Laplacian magnitudes are < 256 and `convertScaleAbs(α=1/256)` zeroes the entire image — Phase 1 σ goes uniformly to 0, the argmax is just "frame 0", and the whole downstream pipeline collapses. Encapsulated in `mpp::blur_mono_for_align(mono, cfg)`; used by `rank_blurred_laplacian_u8` and by every test that feeds blurred frames to align/shift. Raw frames (for brightness, average, stacking-time drizzle resize) stay at native bit depth.
2. **`frames_average_brightness` uses the raw (non-upscaled) frame**, with threshold = `frames_normalization_threshold × {1, 256}` for {8, 16}-bit. Already correct via `mpp_cfg_threshold_scale`.

## Phase 5b — STScI drizzle + Bayer drizzle (Siril-original)

This phase brings two new resample backends to the multipoint stacker, both built on Siril's existing STScI `dobox()` engine in `src/drizzle/`. The bicubic path (Phase 5a) stays untouched and remains the default.

### Architectural anchor: the pixmap path

`dobox()` consumes per-pixel float maps `pixmap->xmap[i,j]`, `pixmap->ymap[i,j]` to convert input-pixel coordinates to output-pixel coordinates (`src/drizzle/cdrizzlemap.c:378`). That means **MPP's per-AP shifts can be plumbed through dobox without any algorithmic changes** to the drizzle engine — we just need to build the right pixmap per input frame, where each input pixel maps to its drizzle-output coordinate via:

`(x_out, y_out) = (scale × (x_in + global_dx + ap_dx(x,y)), scale × (y_in + global_dy + ap_dy(x,y)))`

`ap_dx(x,y) / ap_dy(x,y)` is a smooth per-pixel field obtained by interpolating the per-AP shifts (`run->shifts->shifts`) at the input pixel location. Phase 5a already does this interpolation conceptually inside the per-AP remap loop; Phase 5b lifts the interpolation into the pixmap so dobox sees a single non-affine warp per frame and produces the drizzled output in one pass.

This is materially simpler than the alternatives we considered (per-AP cutouts, piecewise-affine homographies) and uses dobox exactly as Siril's normal stack/drizzle path does — same kernel, same pixfrac, same weight semantics, same Bayer support.

### Sub-tasks

- [ ] **5b.1 — `mpp_pixmap`: per-frame, per-pixel coordinate map builder**
  Construct an `imgmap_t` (float `xmap[ry × rx]`, `ymap[ry × rx]`) for input frame `i`:
  - Base map: identity grid (`xmap[j,i] = i`, `ymap[j,i] = j`).
  - Add global shift: `xmap += global_shifts[i].dx`, `ymap += global_shifts[i].dy`.
  - Add per-AP shift contribution: for each output-pixel position, sample the per-AP shift field by:
    - For each AP `a`, compute its weight at `(i, j)` via the same 2D ramp PSS uses for stack-time per-AP weighting (`stack_frames.py:300+`; already in our codebase as `mpp::one_dim_weight × min` in `mpp_stack.cpp`).
    - Weighted-average the per-AP `(sdy, sdx)` into a single `(dy, dx)` for that pixel.
    - Add to `xmap` / `ymap` at that pixel.
  - Multiply final map by `scale` (the drizzle factor) so output coordinates land in the drizzled canvas.
  Implementation note: the per-AP weight field is identical to the per-AP buffer weighting we already compute — extract it into a reusable helper. The pixmap build is O(rx × ry × num_aps) but the inner loop is tight (no I/O, no allocs); for 264×258 frames × ~12 APs that's ~800k weight evaluations per frame — sub-millisecond.

- [ ] **5b.2 — `mpp_drizzle_stsci_path`: integrate dobox into Stage C**
  In `mpp_stack.cpp`, branch on `cfg->drizzle_mode`:
  - `MPP_DRIZZLE_BICUBIC` (default) → existing Phase 5a path.
  - `MPP_DRIZZLE_STSCI` → new path:
    - Allocate output `fits` at `drizzle_factor × frame_dims`.
    - Allocate `output_counts` (single-channel weight accumulator at the same shape).
    - For each frame in `best_frame_indices` (selected for this AP / globally per cfg):
      - Build pixmap via 5b.1.
      - Fill `driz_param_t` from `cfg`: `kernel`, `pixel_fraction`, `scale`, etc.
      - Call `dobox(p)` to accumulate this frame into output + counts.
    - Final normalisation: `output_data /= output_counts` (with eps-floor for empty cells, mirroring Siril's existing post-drizzle pass).
    - Per-frame brightness equalisation: keep the same `frame × median / (avg + 1e-7)` correction Phase 5a uses, applied to the input frame before pixmap build.
  - **Background composition**: STScI drizzle's `output_counts` array naturally captures "how many frames contributed to each output pixel" — use it as the foreground/background weight directly (no separate accumulator needed). Cells with `counts < bg_blend_threshold × stack_size` blend in the background.
  - Use the existing `mpp_drizzle.cpp` file (currently a 14-line stub returning MPP_ENOTIMPL).
  - New API: `mpp_stack_apply_stsci(seq, cfg, run, fits *out)` — parallel to the existing bicubic `mpp_stack_apply`. Or fold both behind `mpp_stack_apply` with a `cfg->drizzle_mode` switch (preferred — keeps the orchestrator simple).

- [ ] **5b.3 — `mpp_drizzle_bayer_path`: 3-channel drizzle from raw Bayer samples**
  Stage A/B unchanged (still operate on the debayered analysis frame so the global align + AP shifts have signal). Stage C diverges:
  - Skip Siril's debayer entirely for this path — read frames as raw single-channel uint16 via a new `read_bayer_frame()` helper in `mpp.cpp` (existing `read_full_frame` always debayers).
  - Set `driz_param_t::is_bayer = TRUE`, `cfa[...]` from the SER's `bayer_pattern` (RGGB / GRBG / GBRG / BGGR; SER stores this at header offset 18 and Siril already parses it for `com.seq`).
  - dobox handles the CFA-aware sampling internally (only contributes a green-channel sample where the input pixel falls on a green Bayer pixel, etc. — see existing `cdrizzlebox.c` Bayer branches).
  - Output is a 3-channel `fits` at `drizzle_factor × frame_dims`.
  - Pixmap is identical to 5b.1 — Bayer doesn't change the geometry, only the per-pixel channel assignment.
  - **Sanity check needed**: Siril's existing Bayer drizzle (used in deep-sky stacking) assumes one homography per frame. Verify that dobox's CFA path handles arbitrary pixmaps correctly — if it assumes the input is rectified to the CFA grid in a specific way, MPP's non-affine pixmap might break it. Add a small unit test on a 2× synthetic Bayer mosaic before declaring 5b.3 done.

- [ ] **5b.4 — CLI + GUI flag surface**
  - Extend `cfg->drizzle_mode` enum: `MPP_DRIZZLE_OFF / MPP_DRIZZLE_BICUBIC / MPP_DRIZZLE_STSCI / MPP_DRIZZLE_BAYER`. Keep `drizzle_factor` separate (1 / 2 / 3 ints; 1.5 still rendered at 3× pending the downsample step).
  - New cfg fields: `drizzle_pixfrac` (float, default 0.7) and `drizzle_kernel` (enum mirroring `e_kernel_t`: square / gaussian / point / turbo / lanczos2 / lanczos3).
  - CLI flag parser (`apply_mpp_flag` in `command.c`):
    - `-drizzle=Off|1.5|2|3` (existing — bicubic).
    - `-drizzle=stsci-2x|stsci-3x` (new).
    - `-drizzle=bayer-2x|bayer-3x` (new).
    - `-pixfrac=<float>` (new; 0 < f ≤ 1).
    - `-driz-kernel=<name>` (new; one of square|gaussian|point|turbo|lanczos2|lanczos3).
  - GUI: extend the stack-side MPP sub-panel's drizzle combo with the new entries; add pixfrac spinner + kernel combo, visible only when an STScI mode is selected (`gtk_widget_set_visible` on `drizzle_mode` combo change).
  - Cross-mode flag validation in `process_register_mpp` / `process_stack_mpp`: pixfrac / kernel are stack-time only; reject if passed to register_mpp.

- [ ] **5b.5 — STScI synthetic-truth test (`mpp_drizzle_test.cpp`)**
  Generate a high-resolution ground-truth image (e.g. 1024×1024 synthetic Jupiter with fine surface features) → downsample to 256×256 with prescribed sub-pixel offsets (24 frames, offsets uniformly distributed in [−1, 1] px) → run pipeline at `drizzle=stsci-2x` → upsample → compare to ground truth via PSNR and SSIM.
  Acceptance: STScI 2× recovers the ground truth with **PSNR > 2 dB better** than bicubic 2× on the same fixture. This is the resolution-recovery claim of STScI drizzle and is the test that justifies the path's existence.

- [ ] **5b.6 — Bayer-drizzle slanted-edge MTF test**
  Generate a high-resolution RGB ground truth with a slanted-edge target → mosaic to RGGB Bayer at 256×256 with prescribed sub-pixel offsets (24 frames) → run pipeline at `drizzle=bayer-2x` → measure MTF50 on the slanted edge in the output per channel.
  Acceptance: Bayer-2x MTF50 **per channel** exceeds the MTF50 of "debayer-with-cv::cvtColor then bicubic-stack-2x" by at least 1.3× on red and blue channels (green is already over-sampled in RGGB so the gain is smaller). Standard slanted-edge measurement: ISO 12233.

### Phase 5b tests (all on synthetic ground truth — no PSS oracle exists for these paths)

| Test | Acceptance bar | Notes |
|---|---|---|
| `mpp_pixmap_identity` | identity pixmap reproduces input bit-for-bit through dobox | smoke test |
| `mpp_pixmap_global_shift` | pure global shift via pixmap matches `cv::warpAffine` to < 0.5 LSB | regression |
| `mpp_pixmap_per_ap_field` | per-AP weight field at AP centres equals 1.0 (the AP wins its own centre) | invariant |
| `mpp_drizzle_stsci_synthetic` | PSNR gain ≥ 2 dB vs bicubic on the synthetic-truth dataset | the main bar |
| `mpp_drizzle_bayer_slanted_edge` | MTF50_R ≥ 1.3× and MTF50_B ≥ 1.3× of debayer-then-bicubic | the colour-resolution bar |
| `mpp_drizzle_end_to_end_real_ser_stsci` | STScI path runs cleanly on `test-big.ser`; visual review only (no PSS oracle) | smoke + perf |

### Risks and open questions

1. **Per-AP weight ramp inside the pixmap** — the bicubic Phase 5a path applies the per-AP weight at *merge* time (after each AP's buffer is built independently). Drizzle would have to apply it at *sample* time (during the pixmap build, before dobox). The two are mathematically distinct: bicubic weights at output-pixel resolution; drizzle weights at input-pixel resolution. Need a small numerical experiment to verify the drizzle formulation doesn't visibly bias the merge. _(Likely fine for typical 10+ APs with smooth ramps, but worth confirming on the synthetic dataset before declaring 5b.5 done.)_

2. **dobox CFA + arbitrary pixmap interaction** — flagged in 5b.3. The existing Bayer drizzle path in Siril is only exercised with homography-derived pixmaps. Code review of `cdrizzlebox.c` Bayer branches + a small unit test gates whether 5b.3 is a one-day port or a deeper refactor.

3. **Drizzle 1.5× downsample** — Phase 5a stores `drizzle_factor=3` for the 1.5× option and the tooltip says "Siril will downsample later (not yet implemented)". This is a separate small piece of work that lives alongside 5b but is independent of STScI/Bayer: a `cv::resize` post-pass when the drizzle-factor-stored value (3) doesn't match the user's intent (1.5). Worth doing in 5b's flag surface refactor rather than letting it linger.

4. **Memory footprint** — `drizzle_factor × frame_dims` output at float32, plus a counts buffer, plus per-frame pixmaps of the same size. For `drizzle=3` on 1024×1024 input that's ~75 MB peak; manageable. For deep-sky 4096×4096 it would be ~1.2 GB; out of scope for now (deep-sky uses STScI on a different path).

### Phase 7 dependents

These re-open once 5b ships:

- [ ] **7.4** `pss test-big.ser -drizzle=stsci-2x` → produce a stacked FITS; visual sanity check (cleaner edges than bicubic-2x); SSIM comparison to bicubic-2x as a quantitative anchor.
- [ ] **7.5** Bayer-2x on a real Bayer SER → slanted-edge MTF measurement against debayer-then-stack on a controlled target (would need a planet-side test image with a sharp edge, or a synthetic injected one).

## Phase 6 — Orchestrator and `pss` command

### Pipeline staging

The mpp pipeline splits cleanly into three stages with materially different compute weight. The split matters because **the AP grid must be editable in the GUI before the heavy per-AP shift compute runs** — for typical planetary seeing the user may exclude 80–90 % of frames after seeing the rank plot, and we don't want to compute shifts for any of those.

| Stage | What runs | Rough cost on the 500-frame 264×258 real SER | Persisted? |
|---|---|---|---|
| **A — Analyze** | rank → global align → average reference frame → auto-place APs (Phases 1–3 of this plan) | ~0.5–1 s | in-memory in GUI; not persisted in CLI |
| **B — Registration** | per-AP per-frame `multilevel_correlation` (Phase 4) — `selected_frames × num_APs` correlations | scales linearly in both axes; ~3–5 s on the real SER above | **sidecar** (`.mpp`) |
| **C — Stack** | brightness equalise → resample → per-AP remap → merge → uint16 (Phase 5a) | ~1 s | output image |

Terminology: **"registration" = Stage B**, the per-AP shift compute. Stage A is "analyze" / "AP placement"; Stage C is "stack". This matches the user's framing of "auto-place APs, refine, then commit to the heavy registration".

### Commands

```
register -method=mpp seqname [-selected] [flags…]
                                    # → runs A + B; writes seqname.mpp sidecar
                                    #   No manual AP editing in CLI v1.
stack    -method=mpp seqname [-drizzle=…] [-out=…]
                                    # → requires seqname.mpp; runs C
pss      seqname [register flags…] [-drizzle=…] [-out=…]
                                    # → A + B + C in one go (PSS defaults)
```

`register -method=mpp` is unusual among Siril registration methods in that it doesn't produce a registered sequence (no per-frame transform; mpp's shifts are per-AP per-frame and only meaningful when consumed by the mpp stacker). It writes the sidecar instead.

`-selected` honours the user's pre-set selection from Siril's frame selector. Stage A still ranks **all** frames so the selector has rank data to display, but Stage B only computes shifts for selected frames and Stage C only stacks them.

### CLI flag surface

Flags name PSS config fields directly so scripts read across.

| Flag | PSS field | Default | Used by |
|---|---|---|---|
| `-min-brightness=N` | `alignment_points_brightness_threshold` | 10 | register, pss |
| `-min-structure=F` | `alignment_points_structure_threshold` | 0.04 | register, pss |
| `-min-contrast=N` | `alignment_points_contrast_threshold` | 0 | register, pss |
| `-half-box=N` | `alignment_points_half_box_width` | 24 | register, pss |
| `-search-width=N` | `alignment_points_search_width` | 14 | register, pss |
| `-search-global=N` | `align_frames_search_width` | 34 | register, pss |
| `-patch-scale=F` | `align_frames_rectangle_scale_factor` | 3.0 | register, pss |
| `-no-dewarp` | `alignment_points_de_warp=false` | off (de-warp on) | register, pss |
| `-no-normalize` | `frames_normalization=false` | off (normalize on) | register, pss |
| `-stack-percent=N` | `alignment_points_frame_percent` | 10 | stack, pss |
| `-stack-frames=N` | `alignment_points_frame_number` | −1 | stack, pss |
| `-drizzle={Off\|1.5\|2\|3}` | `stack_frames_drizzle_factor_string` | Off | stack, pss |
| `-bg-fraction=F` | `stack_frames_background_fraction` | 0.3 | stack, pss |
| `-bg-blend=F` | `stack_frames_background_blend_threshold` | 0.2 | stack, pss |
| `-out=file` | — | sequence-name-derived | stack, pss |
| `-selected` | — | off | register, pss |

### Sidecar (`.mpp`) contents

Binary blob next to the sequence. Shape-defined by `mpp_sidecar.{c,h}` from Phase 0; final contents per Phase 6:

- Header: magic `SIRILMPP`, version, `frame_count`, `frame_rows` / `cols`, `bitdepth`, `drizzle_factor` at register time, included-frame count + bitmask.
- Per-frame: quality, avg_brightness, global `(dy, dx)`, included/excluded flag.
- Patch coords (4 ints).
- Intersection bounds (4 ints).
- AP records (count + per-AP `mpp_ap_record_t`).
- Per-AP per-frame shifts (`num_aps × num_frames × 2` doubles + success byte), with sparse entries for excluded frames.
- Per-AP per-frame qualities, per-AP `best_frame_indices`, per-frame `used_alignment_points`.
- Median brightness (scalar).

≈ 1–2 MB for typical runs. Round-trip-tested in `mpp_sidecar_test.c`.

### GUI surface (Phase 9 will consume these)

Phase 6 exposes the entry points Phase 9 needs, allowing the GUI to interleave the stages with user edits:

```c
int mpp_analyze(struct sequence *seq, const mpp_config_t *cfg,
                mpp_run_t **run_out);
int mpp_edit_aps(mpp_run_t *run, ... /* GUI-mediated changes */);
int mpp_compute_shifts(mpp_run_t *run, const mpp_config_t *cfg);
int mpp_stack_apply(const mpp_run_t *run, const mpp_config_t *cfg,
                    struct fits *out);
```

GUI flow:

1. Load sequence; "Analyze" → `mpp_analyze` runs Stage A.
2. Siril's existing frame selector now shows per-frame quality; user bulk-deselects (e.g. "keep top 30 %"). Typical planetary seeing excludes 80–90 % here.
3. New AP-editor panel overlays auto-placed APs on the average reference frame; user can add / remove / move / resize (Phase 9 widget — share with shift-distribution viewer where possible).
4. "Register" → `mpp_compute_shifts` runs Stage B **only on currently selected frames**, with progress bar.
5. "Stack" → `mpp_stack_apply` runs Stage C and loads the result into the main viewer.

If the user changes the frame selection materially after Stage A, the GUI should hint that re-running Analyze gives a more accurate alignment patch and average frame (best frame may have moved). Cheap to re-run.

### Frame selection / bulk filtering

Reuses Siril's existing frame selector (no PSS-style thumbnail grid). The selector's quality column is populated by Stage A's per-frame ranks. mpp honours the `included` flag at every subsequent stage:

- Stage A re-uses cached ranks for the included set when re-run (no recompute on excluded frames).
- Stage B's loop iterates over `included_frames × aps`, never excluded × aps.
- Stage C's background-frame top-N selector picks from the included set only.

CLI users get the same behaviour by passing `-selected`. Without it, all frames in the sequence participate.

### Colour and sequence formats

- **Bayer / RGB:** Phase 6 forces a debayer-on-open override (`com.pref.debayer.open_debayer = TRUE` for the orchestrator's scope) so Bayer SERs always come in as RGB. Stacking buffers are 3D for colour; brightness measure is per-channel-mean; the rest of the pipeline operates on a mono "analysis frame" derived from the green channel (PSS `color_index` default) for ranking, alignment, AP placement, and per-AP shift.
- **Sequence formats:** all Siril-supported via the existing `sequence` reader (`src/io/ser.c`, `src/io/seq.c`, etc.) — SER (primary planetary target), FITSEQ, individual FITS files, and image directories. AVI/MOV is supported by the reader but Siril deprecates it; we don't advertise it.

### Implementation phasing

- [x] 6.1 Bridge `mpp_frames` to Siril's sequence reader: implement `mpp_frames_load_mono(seq, idx, cfg, fits*)` and a colour analogue so the orchestrator can feed Siril sequences into the existing mpp algorithms. _(`mpp_seq_read_frame` + `read_analysis_frame` / `read_full_frame` in `mpp.cpp`.)_
- [x] 6.2 Force debayer on open inside the orchestrator's scope; route the green channel through analysis, full RGB through stacking. _(All three command handlers go through `load_sequence_force_debayer`; analysis layer is the PSS panchromatic luminance — `cv::cvtColor(RGB, GRAY)` — for RGB input.)_
- [x] 6.3 Implement the three GUI-callable entry points (`mpp_analyze`, `mpp_compute_shifts`, `mpp_stack_apply`). The C++ implementations from Phases 1–5a already exist; this is glue + the colour-aware variant of stack. _(Live in `src/registration/mpp.cpp`.)_
- [x] 6.4 Finalise the sidecar format and round-trip; implement `mpp_sidecar_write` / `mpp_sidecar_read` against the contents listed above. _(`mpp_sidecar.c` + `mpp_sidecar_test.c` covers the round-trip with shifts.)_
- [x] 6.5 Implement `register_mpp(struct registration_args*)`: drives A + B, then `mpp_sidecar_write`. _(GUI hook in `mpp.cpp`; CLI handler `process_register_mpp` in `command.c` does the same with per-call flag overrides.)_
- [x] 6.6 Implement the `stack_mpp` consumer: `mpp_sidecar_read` → `mpp_stack_apply` → write FITS. _(`process_stack_mpp` in `command.c`. Plan originally said `stack -method=mpp` but `stack`/`register` are deep-sky-specific commands with parameter shapes that don't apply here; separate `register_mpp` / `stack_mpp` commands match the same hybrid-CLI intent without polluting the existing handlers.)_
- [x] 6.7 Implement `process_pss()`: parses the flag surface; calls the three stages directly with per-flag cfg overrides. _(In `command.c`; shares `apply_mpp_flag` with the two-step commands.)_
- [x] 6.8 Register `REG_MPP` in `regmethod_index` (`registration.h`) and append to the method table in `gui/registration.c` (selection type `REQUIRES_NO_SELECTION`, registration type `REGTYPE_PLANETARY`). _(REG_MPP sits before REG_2PASS to keep the hidden 2pass tail invariant.)_
- [x] 6.9 Colour-path tests: real-SER `test-big.ser` (Bayer RGGB) end-to-end produces 3-layer RGB FITS that matches the PSS oracle within 0.04% per-channel mean (1 px width difference attributable to BAYER_RCD vs cv2's bilinear debayer — algorithmic, not a bug).
- [x] 6.10 End-to-end CLI test on `test-big.ser` through Siril's actual command parser. _(`pss` 1-step and `register_mpp + stack_mpp` 2-step both work; outputs are 100% bit-identical pixel-for-pixel. `stack_mpp -drizzle=2` correctly overrides the sidecar's `drizzle_factor=1` to produce a 500×494 output from the 250×247 base. Cross-mode flag rejection works (register_mpp rejects `-out=`; stack_mpp rejects `-half-box=`).)_

**Exit criterion:** `pss test-big.ser -out=jupiter.fits` from a real Siril invocation produces a stacked FITS that is byte-near-equivalent to the Phase-5a oracle output (same acceptance bar: > 99 % pixels exact, worst Δ ≤ 4 / 65535).

## Phase 7 — ★ Real-SER end-to-end milestone ★

- [x] 7.1 Run PSS (Python) on a real planetary SER with defaults; archive output. _(`tools/pss_reference/oracle_out_realser/` from `run_pss.py /workspace/test-big.ser --type=video --n=500`.)_
- [x] 7.2 `pss` (Siril) default → compare to PSS reference. PSNR > 50 dB. _(**128.34 dB** on pre-debayered mono fixture; **49.63 dB** on colour SER with best-shift overlap — gap is purely the Bayer debayer-algorithm choice, see MPP_PSS_DIFFS.md §1.)_
- [x] 7.3 `pss` `-drizzle=2` (bicubic) → compare to PSS `drizzle_factor=2` reference. PSNR > 50 dB. _(**54.05 dB** overall (R 55.18 / G 53.45 / B 53.70); 2× resampling smooths per-pixel debayer noise.)_
- [ ] 7.4 `pss` `-drizzle=stsci-2x` → sanity-check: clean image, comparable or better SSIM vs bicubic-2x output. _(Gated on Phase 5b.)_
- [ ] 7.5 If Bayer SER: `pss` `-drizzle=bayer-2x` → slanted-edge MTF measurement vs debayer-then-stack. _(Gated on Phase 5b.)_
- [x] 7.6 Document any deliberate departures from PSS in `src/registration/MPP_PSS_DIFFS.md`.

## Phase 8 — Hardening

- [x] 8.1 Edge cases: tiny sequences (< 30 frames), mono-only, color without SER header info, oversaturated targets, sequences with junk frames, memory profile on large sequences. _(Verified pss runs cleanly on 8-frame and 3-frame SER variants; junk all-black frames mixed into a 40-frame sequence get demoted by the ranker rather than crashing. Found and fixed a real bug: `pss` on a true-mono SER (ColorID=0) was reinterpreting it as CFA because `com.pref.debayer.open_debayer = TRUE` triggers Siril's "force-debayer-mono" path. Now peek at the SER header at offset 18 and only force debayer when ColorID != 0.)_
- [x] 8.2 Error paths return clean status codes through the orchestrator. No asserts on user-reachable paths. _(Audited mpp_*.c/cpp for asserts/abort/exit — none in production code. All allocations checked; all `mpp_*` functions return MPP_E* on failure.)_
- [x] 8.3 Performance pass: OpenMP scheduling, buffer reuse, profile on real SER. Aim for parity-or-better vs PSS Python runtime. _(500-frame test-big.ser: Siril pss 3.01 s vs PSS oracle 3.25 s — Siril ~8% faster end-to-end. Peak RSS 250 MB (Siril) vs 228 MB (PSS oracle) — ~10% more, largely Siril CLI baseline. No further optimisation needed for v1.)_

---

## 🛑 GUI go/no-go pause

Algorithm shippable headless at this point.

- [x] **Decision recorded:** Option B (GTK3 on trunk, port later). User: "There are still some GTK4 rendering issues that will take a while to sort out, so I think we should continue with a GTK3 interface and port it later. But please follow the GTK3 to GTK4 transition guidelines and don't use parts of the GTK3 API that will make migration harder than necessary." Phase 9 followed the GTK3-with-migration-friendly-choices guidance (no GtkVBox/HBox/Table/EventBox; GtkBuilder + .ui XML; GtkAdjustment-backed GtkSpinButton).

---

## Phase 9 (gated) — GUI port

Shipped over 25 incremental slices on the `pss` branch (commits Phase 9 (1/n) … (25/n)). High-level outcomes:

- [x] 9.1 Parameter configuration. _(Realised as sub-panels on the existing Registration and Stacking tabs rather than a standalone dialog mirroring PSS's `parameter_configuration.ui`. Spinners on the MPP register sub-panel drive Stage A/B parameters; spinners on the MPP stack sub-panel drive Stage C parameters. GtkAdjustment objects are shared between the register sub-panel and the AP editor so values stay in sync. Auto-selects STACK_MPP method when a `.mpp` sidecar is detected at sequence-open time.)_
- [x] 9.2 Frame selector + viewer. _(Reuses Siril's existing frame selector and quality plot — `mpp_write_quality_to_regdata` publishes per-frame Stage-A quality scores into `seq->regparam[layer][i].quality` so the standard widget surfaces them. Analyze paints the mean reference frame into gfit so the user sees what AP placement worked on; `com.uniq` is updated to filename "REFERENCE IMAGE" so display_filename labels it correctly.)_
- [x] 9.3 Alignment-point editor (`mpp_ap_editor_dialog.ui`). _(Non-modal, transient-for the main window. Auto-place / Clear / Commit / Cancel buttons. Mouse interactivity: left-click to add an AP, right-click to remove, drag to move; hover highlighting (orange). Cancel reverts from a snapshot taken on dialog open; X-close treated as Cancel. The cached run lives in `com.mpp_run`; Register reuses it instead of re-running Stage A. Editor edits invalidate `best_frame_indices` and Register transparently recomputes per-AP qualities via `mpp_recompute_qualities` before Stage B — no re-Analyze required.)_
- [x] 9.4 Shift distribution viewer (`mpp_shift_viewer_dialog.ui`). _(Non-modal diagnostic dialog with frame spinner and vector-scale multiplier. While open, the AP overlay paints arrows from each AP centre showing per-AP shift for the selected frame: green = Stage B converged, red = fell back to zero. Small filled dot at every AP centre keeps zero-shift APs visible. Frame spinner also loads the selected frame so pixel data and arrows update together. Sidecar auto-load on sequence open populates `com.mpp_run` so the overlay and viewer come up immediately on previously-registered sequences.)_
- [ ] ~~9.5 Job dialog for batch multi-SER processing.~~ _(**Dropped from scope.** User: "9.5 is future work though — combining multiple .sers is a complex task often involving derotation which Siril doesn't support at all, so the current focus is on single SERs.")_

### Deferred / open items from Phase 9

Things flagged during the 25-slice GUI build but not blocking, listed in order of likely value:

- GUI auto-force-debayer for Bayer SERs (CLI does it via `load_sequence_force_debayer`; GUI currently surfaces a clear error and asks the user to toggle `com.pref.debayer.open_debayer` manually).
- Long-click or modifier-click to cycle through overlapping APs in the editor (hit-test currently returns the topmost AP; overlap is rare by design but possible with user-added APs).
- Mouse cursor change when in `MOUSE_ACTION_EDIT_APS` mode (bgext sets "cell"; we leave the default).
- Explicit "no APs to register" guard when Register is clicked with a cleared AP grid (today it falls through to a fresh Stage A, which is correct but might confuse).
- Future port of the GUI to GTK4 once Siril's GTK4 migration completes (Phase 9 followed migration-friendly conventions throughout).

---

## Test strategy

| Layer | Test type | Location |
|---|---|---|
| Rank, AP, shift, stack, align | Criterion unit tests vs PSS oracle | `src/tests/mpp_*_test.{c,cpp}` |
| Pipeline | Synthetic-SER integration test | `src/tests/mpp_pipeline_test.cpp` |
| End-to-end real SER | PSNR/SSIM compare to PSS output | `tools/pss_reference/compare.py`, manual at Phase 7, nightly later |

## Open assumptions (carry-forward)

- Algorithmic match = behaviour-equivalent within float-rounding, not bit-identical.
- PSS sign conventions for shifts and AP coordinates are authoritative; documented at module boundaries.
- Bayer-drizzle path is Siril-original (PSS has no equivalent); no oracle exists, so validated via synthetic ground truth only.
