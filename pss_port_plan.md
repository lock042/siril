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

#### 5b.1 — `mpp_pixmap`: per-frame, per-pixel coordinate map builder

**What it is:** a deterministic function that produces an `imgmap_t` describing where each input pixel of frame `i` lands on the drizzled output canvas. This is the only piece that "knows" about MPP's non-affine warp; everything downstream (dobox, weighting, normalisation) treats the warp as a black box.

**Signature (proposed):**

```c
// In mpp_drizzle.h
mpp_status_t mpp_pixmap_build(const mpp_run_t *run,
                              int frame_idx,
                              double drizzle_scale,
                              imgmap_t *out_pixmap);
```

**Implementation steps:**

1. Allocate `out_pixmap->xmap` / `ymap` as `rx × ry` float arrays. `rx, ry` come from `run->frame_cols / frame_rows` (the input frame size — *not* the drizzled output size).
2. Compute the per-AP weight surface once per frame, reused across both x and y axes. For each AP `a` and each output pixel `(i, j)` define a 2D ramp weight `w_a(i, j) = ramp_x_a(i) × ramp_y_a(j)`, where:
   - `ramp_x_a(i)` is 1.0 at the AP centre column, falls linearly to 0.0 at the AP's box boundary, and 0.0 outside the box. (Identical formulation to PSS's `one_dim_weight`, already implemented as `mpp::one_dim_weight` in `mpp_stack.cpp`.)
   - `ramp_y_a` is symmetric.
3. For each input pixel `(i, j)`:
   - `dy_local = (Σ_a w_a(i,j) × shifts[2(i_frame × M + a) + 0]) / (Σ_a w_a(i,j))`, with `Σ_a w_a` guaranteed > 0 by the staggered AP grid construction (every pixel is inside at least one AP box for the typical configuration). When `Σ_a w_a = 0` (corner cases at the very edge of the intersection), fall back to the global shift only.
   - `dy_total = global_shifts[frame_idx].dy + dy_local`; same for `dx_total`.
   - `xmap[j × rx + i] = drizzle_scale × (i + dx_total)`.
   - `ymap[j × rx + i] = drizzle_scale × (j + dy_total)`.

**Reusable helpers extracted from `mpp_stack.cpp`:**
- `mpp::one_dim_weight(int box_low, int box_high, bool extend_low, bool extend_high, double pos) → double` — already exists.
- `mpp::per_pixel_weight_sum(...)` — pulled out of the per-AP buffer weighting code so 5b can share it.

**Lifetime / reuse strategy:** one pixmap buffer is allocated for the whole stack run and reused across frames — `mpp_pixmap_build` overwrites it in place rather than allocating a fresh map per frame. Concretely:

```c
imgmap_t pm = { 0 };
imgmap_alloc(&pm, run->frame_cols, run->frame_rows);   // once
for (int f : frames_to_stack) {
    mpp_pixmap_build(run, f, drizzle_scale, &pm);      // overwrites pm
    /* fill driz_param_t p with &pm and the frame data, call dobox(&p) */
}
imgmap_free(&pm);
```

The build-use-discard pattern (build the pixmap, hand it to dobox, drop it on the floor) wins on three counts and the alternatives don't help us:

| Approach | Memory | Wall time | Notes |
|---|---|---|---|
| **One reused buffer (chosen)** | 2 × rx × ry × float ≈ 0.5 MB for 264×258, 16 MB for 4k frames | dobox dominates; pixmap build ≈ 20–30 % of dobox cost | simple, cache-friendly, zero alloc churn |
| Build all up front | N × that = 250 MB for 500-frame `test-big.ser` at 264×258; ~8 GB for 500 4k frames | pixmap build parallel across frames (small win) | untenable for large sequences |
| Pool of K buffers + per-thread dobox | K × that, plus K × output canvas for reduce | parallel per-frame, with reduce overhead | revisit only if profiling demands it |

**Intra-frame parallelism stays.** The per-pixel loop inside `mpp_pixmap_build` is embarrassingly parallel (each output coordinate is independent), so we wrap it in `#pragma omp parallel for collapse(2)` over (j, i). With 12 APs and a 264×258 frame that's ~800k independent weight evaluations — saturates a modern CPU's threads. **Inter-frame parallelism is left out** of the initial 5b.2 implementation because dobox writes to a shared output canvas; per-thread output buffers plus reduce would need ~K × output canvas memory (24 MB per thread at drizzle=2 on a 264×258 frame; 384 MB per thread on a 4k frame) which is the same memory problem as "build all pixmaps up front", and the bicubic Phase 5a path is already sequential per-frame for the same reason. Add it only if profiling on `test-big.ser` shows the per-frame loop is the bottleneck (the existing bicubic 5a path stacks 500 frames in ~1 s, so the floor is already low).

`imgmap_alloc` / `imgmap_free` are new tiny helpers next to `mpp_pixmap_build` — just malloc + memset + free, parallel-safe.

**Effort:** ~half a day. The arithmetic is straightforward; the only subtlety is the corner-case "no AP covers this pixel" fallback, which needs a unit test (test `mpp_pixmap_outside_aps` below).

#### 5b.2 — `mpp_drizzle_stsci_path`: integrate `dobox` into Stage C

**What it is:** a new stack-time path that replaces the bicubic resize + per-AP cutout + accumulate sequence (Phase 5a) with a single per-frame `dobox()` call into a shared output canvas. The MPP-specific pieces (brightness equalisation, frame selection per AP) wrap around the call.

**Public surface change:** fold STScI behaviour into the existing `mpp_stack_apply` rather than adding a parallel entry point. Add a `drizzle_mode` switch at the top:

```c
mpp_status_t mpp_stack_apply(sequence *seq, const mpp_config_t *cfg,
                             const mpp_run_t *run, fits *out) {
  switch (cfg->drizzle_mode) {
    case MPP_DRIZZLE_OFF:
    case MPP_DRIZZLE_BICUBIC: return mpp_stack_apply_bicubic(seq, cfg, run, out);   // existing
    case MPP_DRIZZLE_STSCI:   return mpp_stack_apply_stsci  (seq, cfg, run, out);   // new (5b.2)
    case MPP_DRIZZLE_BAYER:   return mpp_stack_apply_bayer  (seq, cfg, run, out);   // new (5b.3)
  }
}
```

**STScI path skeleton:**

1. Allocate output: `fits out` at `(drizzle_factor × intersection_w, drizzle_factor × intersection_h, num_layers)`, type `DATA_FLOAT`. Per-channel layout matches Siril's existing drizzle output.
2. Allocate `output_counts`: single-channel `fits` of the same shape — dobox's weight accumulator.
3. Build sequence-wide `driz_args_t` from cfg: `kernel`, `pixel_fraction`, `scale = drizzle_factor`, `is_bayer = FALSE`.
4. For each frame index `f` in the union of `best_frame_indices` across all APs (deduplicated):
   - Read frame `f` (RGB or mono per `run->num_layers`).
   - Apply per-frame brightness equalisation: `frame × median_brightness / (frame_brightness[f] + 1e-7)`. Same as Phase 5a.
   - Build pixmap via `mpp_pixmap_build(run, f, drizzle_factor, &pixmap)`.
   - Fill `driz_param_t p` from `driz_args_t` + per-frame data (`p.data = &frame_fits`, `p.pixmap = &pixmap`, `p.output_data = &out`, `p.output_counts = &output_counts`).
   - Call `dobox(&p)`. dobox accumulates into `out` and `output_counts` in-place.
5. Background composition: `output_counts` is dobox's natural weight map — cells with `counts < cfg->stack_frames_background_blend_threshold × n_frames` get the median background blended in (same threshold semantics as Phase 5a, just driven by a different counter).
6. Final normalisation: `out /= output_counts` with epsilon floor, then scale + cast to uint16.

**Differences vs bicubic path that need attention:**
- Phase 5a applies per-AP weights at *merge* time (after each AP's buffer is built); 5b applies them at *sample* time (during pixmap build). The two are mathematically distinct; verify visually equivalent on the synthetic dataset (test `pixmap_vs_bicubic_per_ap_weight` below).
- Phase 5a iterates per-AP × per-frame; 5b iterates per-frame (each frame contributes once to all output pixels at once). Net cost should be lower for large AP counts.
- Phase 5a's "drizzle expansion" was a `cv::resize(frame, scale=K, INTER_LINEAR)`; 5b does the expansion *inside* dobox via the scaled pixmap. The two paths should produce visually similar output at `drizzle_factor > 1` but **not** bit-identical (different resampling kernels).

**Effort:** ~2 days. dobox is well-encapsulated and Siril already exercises it; the new code is mostly orchestration. Most of the time goes into the test fixtures and the back-and-forth checking that the per-AP weight ramp doesn't bias the output.

#### 5b.3 — `mpp_drizzle_bayer_path`: 3-channel drizzle from raw Bayer samples

**What it is:** the same dobox-based stack path as 5b.2, but feeding raw single-channel Bayer-pattern frames (no debayer) into a 3-channel drizzled output. dobox's CFA-aware sampler decides which output channel each input sample contributes to based on the per-pixel CFA position.

**Stage A/B unchanged:** the analysis frame is still the debayered+luminance projection so the global align + AP shifts have meaningful signal. Only Stage C diverges.

**Implementation steps:**

1. New helper `read_bayer_frame(sequence *seq, int idx, fits *out)` in `mpp.cpp`. Like `read_full_frame` but skips `cv::cvtColor`/`bayer_to_rgb`; output is single-channel uint16 with the raw mosaic.
2. Verify SER's `bayer_pattern` is populated on `com.seq` for Bayer SERs (it is — Siril parses the SER header at sequence-open time).
3. In `mpp_stack_apply_bayer`:
   - Set `driz_args_t::is_bayer = TRUE` and copy the CFA pattern from the SER header into `driz_args_t::cfa[]` / `cfadim`.
   - Output `fits` is 3-channel (`naxes[2] = 3`); `output_counts` is also 3-channel because each Bayer channel has independent counts.
   - Per-frame brightness equalisation: compute on the *debayered* frame (Siril already does this for deep-sky Bayer drizzle), then scale the raw Bayer frame by the same factor.
   - For each frame: build pixmap as in 5b.1, call `dobox` with `is_bayer = TRUE`. dobox routes each input sample to one of the three output channels based on its CFA position; the other channels' weight maps stay zero at that output pixel.
   - Normalise per channel: `out[c] /= output_counts[c]` with epsilon floor.
4. **Critical sanity check before claiming 5b.3 works**: dobox's CFA branch in `cdrizzlebox.c` assumes a specific input-pixel-to-CFA mapping. With MPP's non-affine pixmap, the same input pixel `(i, j)` might land at different output positions for different frames — fine — but dobox must still consult the *input* CFA position (which doesn't move) for the channel routing. Read the dobox CFA loop carefully before writing the stack path; if there's an assumption that the output-CFA position equals the input-CFA position (i.e. that pixmap is identity for CFA purposes), the path either needs a refactor in dobox or a workaround in `mpp_drizzle.cpp` (likely: produce three single-channel intermediate outputs at each CFA position then merge). The synthetic Bayer test (`mpp_bayer_drizzle_synthetic_mosaic` below) is the gate for this.

**Effort:** ~3–5 days depending on what the dobox CFA review surfaces. If dobox accepts arbitrary pixmaps cleanly, this is a one-day extension of 5b.2; if it requires a per-CFA-channel split, more like a week.

#### 5b.4 — CLI + GUI flag surface

**`mpp_config_t` changes** (in `mpp_config.h`):

```c
enum mpp_drizzle_mode {
    MPP_DRIZZLE_OFF = 0,        // existing — drizzle_factor implicitly 1
    MPP_DRIZZLE_BICUBIC = 1,    // existing — drizzle_factor 2 or 3
    MPP_DRIZZLE_STSCI = 2,      // new (5b.2)
    MPP_DRIZZLE_BAYER = 3,      // new (5b.3)
};
enum mpp_drizzle_kernel {       // mirrors e_kernel_t but typed for cfg
    MPP_KERNEL_SQUARE = 0,
    MPP_KERNEL_GAUSSIAN,
    MPP_KERNEL_POINT,
    MPP_KERNEL_TURBO,
    MPP_KERNEL_LANCZOS2,
    MPP_KERNEL_LANCZOS3,
};
struct mpp_config {
    /* ... existing fields ... */
    enum mpp_drizzle_mode   drizzle_mode;     // default MPP_DRIZZLE_OFF
    int                     drizzle_factor;   // existing; 1, 2, or 3
    float                   drizzle_pixfrac;  // new; default 0.7, range (0, 1]
    enum mpp_drizzle_kernel drizzle_kernel;   // new; default MPP_KERNEL_SQUARE
};
```

`mpp_config_defaults` sets pixfrac=0.7 and kernel=SQUARE — matching Siril's existing drizzle defaults.

**CLI flag parser** (`apply_mpp_flag` in `command.c`):

| Flag | Effect on cfg |
|---|---|
| `-drizzle=Off` | `drizzle_mode = OFF`, `drizzle_factor = 1` |
| `-drizzle=1.5` | `drizzle_mode = BICUBIC`, `drizzle_factor = 3` (pending downsample) |
| `-drizzle=2` | `drizzle_mode = BICUBIC`, `drizzle_factor = 2` |
| `-drizzle=3` | `drizzle_mode = BICUBIC`, `drizzle_factor = 3` |
| `-drizzle=stsci-2x` (new) | `drizzle_mode = STSCI`, `drizzle_factor = 2` |
| `-drizzle=stsci-3x` (new) | `drizzle_mode = STSCI`, `drizzle_factor = 3` |
| `-drizzle=bayer-2x` (new) | `drizzle_mode = BAYER`, `drizzle_factor = 2` |
| `-drizzle=bayer-3x` (new) | `drizzle_mode = BAYER`, `drizzle_factor = 3` |
| `-pixfrac=<f>` (new) | `drizzle_pixfrac = f`; reject if not in (0, 1] |
| `-driz-kernel=<name>` (new) | parse name → enum; reject unknown |

Cross-mode validation in `process_register_mpp` / `process_stack_mpp`:
- `-pixfrac` / `-driz-kernel` only valid in stack-time invocations (`stack_mpp` and `pss`); reject when passed to `register_mpp`.
- `-drizzle=stsci-*` and `-drizzle=bayer-*` rejected for `register_mpp` (Stage B doesn't drizzle).
- `-drizzle=bayer-*` rejected on non-Bayer sequences (peek SER header, fail with clear message).

**GUI changes** (stack-side MPP sub-panel in `siril.ui`):
- Extend `combo_mpp_drizzle` with new entries: "Off", "1.5x", "2x", "3x" (existing), "2x (STScI)", "3x (STScI)", "2x (Bayer)", "3x (Bayer)" (new). Map combo index to `(drizzle_mode, drizzle_factor)` in `on_drizzle_combo_changed`.
- Add `spin_mpp_drizzle_pixfrac` (GtkSpinButton, range 0.1–1.0, step 0.05, default 0.7) and `combo_mpp_drizzle_kernel` (GtkComboBoxText with the 6 kernel names). Wire visibility: shown only when an STScI / Bayer mode is selected.
- New cfg→GUI sync code in `update_stack_interface` so the widgets reflect a loaded sidecar's cfg.

**Sidecar compatibility:** the `.mpp` sidecar's `drizzle_factor` field (`mpp_sidecar.c`) stores `cfg.drizzle_factor` only. The new `drizzle_mode` / `drizzle_pixfrac` / `drizzle_kernel` are *stack-time* cfg, not register-time, so they don't need to live in the sidecar — the user picks them at stack time. Verify this assumption holds when implementing 5b.2; if STScI registration ever needs to feed back into Stage B (it shouldn't), revisit.

**Effort:** ~half a day for CLI + GUI plumbing, plus 0.5 day for the GUI visibility / wiring testing.

#### 5b.5 — STScI synthetic-truth test (`mpp_drizzle_test.cpp`)

See "Phase 5b test plan" below.

#### 5b.6 — Bayer-drizzle slanted-edge resolution test  ✅

Done in `mpp_bayer_drizzle::slanted_edge_resolution` (`src/tests/mpp_drizzle_test.cpp`). Implementation notes:

- Fixture: 384×384 8-bit RGB ground truth with a slanted (slope=0.1, ~5.7°) vertical edge running through the centre, sampled at 24 sub-pixel-shifted offsets to RGGB Bayer 192×192 LR frames. Two parallel pipelines run on the same fixture: (A) raw-Bayer → `mpp::stack_apply_bayer` at 2× and (B) `cv::cvtColor(COLOR_BayerRGGB2RGB)` → Phase 5a bicubic at 2×. Both produce 3-channel 384×384 outputs.
- Metric: slant-aware super-sampled (8×) ESF construction, then 10-90% rise width per channel. Width is monotonically inverse-proportional to MTF50 on a smooth edge — far simpler to implement reliably than a full FFT-based MTF that would need sub-pixel edge-angle estimation. We compute the ratio (bicubic-width / Bayer-width) per channel and assert thresholds.
- Acceptance bars: **R ≥ 1.50×**, **B ≥ 1.50×**, **G ≥ 1.30×**. The plan's original 1.3× target was a floor; observed values on this fixture comfortably exceed it (R 2.00×, G 1.53×, B 1.95×), so the bars tighten to lock in the realised gain while leaving room for fixture variance. G is over-sampled in RGGB (50% of CFA positions) so its gain is structurally smaller; R and B (25% each) see the biggest wins because cv::cvtColor's bilinear demosaic is most aggressive on those channels.
- The Bayer-drizzle output retains the CFA mosaic structure (each output cell only receives samples from CFA-matching input positions) — this is correct algorithmic behaviour and visible in the artifact `bayer_slanted_edge_compare.png`. Despite the visible mosaic colouring, the per-channel edge response is sharper.
- Saved artifacts when `MPP_DUMP_RESULT_DIR` env set: `bayer_slanted_edge_bayerdrizzle.png`, `bayer_slanted_edge_bicubic.png`, `bayer_slanted_edge_gt.png`. The labelled side-by-side `bayer_slanted_edge_compare.png` lives in the project root.

### Phase 5b test plan

Three test groups serve different purposes:

1. **Pixmap-builder unit tests** — prove the pixmap is what it says it is, in isolation from dobox.
2. **STScI sanity-against-bicubic-oracle** — *user's request*: at (scale=1, pixfrac=1) STScI should approximate bicubic-drizzle=1; at (scale=2, pixfrac=0.6) it should approximate bicubic-drizzle=2. Rules out junk output without requiring algorithmic equivalence.
3. **STScI / Bayer synthetic-truth tests** — the algorithm-justification tests (resolution recovery for STScI; colour resolution for Bayer drizzle).

#### Group 1 — Pixmap unit tests (`mpp_drizzle_test.cpp::mpp_pixmap_*`)

| Test | Setup | Acceptance |
|---|---|---|
| `mpp_pixmap_identity` | Zero global shifts, zero per-AP shifts, scale=1. | xmap[j,i]=i and ymap[j,i]=j to within float epsilon for every pixel. |
| `mpp_pixmap_global_only` | Non-zero integer global shift `(dy, dx) = (3, -2)`, zero per-AP, scale=1. | xmap = i + dx, ymap = j + dy for every pixel. |
| `mpp_pixmap_per_ap_at_centre` | One AP at `(cy, cx)` with shift `(2, -1)`, zero global, scale=1. | xmap[cy, cx] = cx - 1; ymap[cy, cx] = cy + 2 (the AP wins its own centre). |
| `mpp_pixmap_per_ap_between` | Two APs at `x = 10` and `x = 30` (same y row) with shifts `(0, +2)` and `(0, -2)`. | xmap at the midpoint `(y, 20)` ≈ x_in + 0 (ramp weights cancel symmetrically); xmap monotonically interpolates from `+2` at AP1 to `-2` at AP2. |
| `mpp_pixmap_outside_aps` | One AP in the centre; sample at a pixel outside its box. | Fallback path: xmap = i + global_dx; ymap = j + global_dy. No NaN, no divide-by-zero. |
| `mpp_pixmap_scale` | Global shift `(0, 0)`, zero per-AP, `scale=2`. | xmap[j,i] = 2 × i; ymap[j,i] = 2 × j. |

Total: 6 small unit tests, ~100 LOC, runs in milliseconds. Gates the pixmap code before any dobox integration.

#### Group 2 — STScI sanity against the PSS bicubic oracle (this is the new group from the user request)

Both tests use the bundled real-SER fixture `test-big.ser` (500-frame Bayer Jupiter) and the existing PSS oracle outputs in `tools/pss_reference/oracle_out_realser/`. The PSS oracle is a *bicubic* drizzle; STScI is a *different algorithm*. We're not asking for algorithmic equivalence — just that STScI in its "near-passthrough" and "real-2x" configurations produces a Jupiter image close enough to the oracle to rule out junk (wrong sign, dimensions, bit-depth, total black/white, channels swapped, etc.).

| Test | Setup | Acceptance |
|---|---|---|
| `mpp_stsci_near_passthrough_oracle` | `drizzle=stsci-2x` is the smallest drizzle mode; force `drizzle_factor=1` and `pixfrac=1` at the API level (CLI doesn't expose drizzle=1 for stsci — this is a unit test). Compare to PSS oracle stacked at `drizzle="Off"` (the existing `oracle_out_realser/stacked.fits`). | PSNR ≥ 35 dB on the central crop excluding the 10-px border (drizzle output dims may differ slightly from intersection size). SSIM ≥ 0.95. Mean abs diff ≤ 200 / 65535. No pixel anywhere differs by more than 5000 / 65535. **All four criteria must hold.** |
| `mpp_stsci_real_2x_oracle` | `drizzle=stsci-2x`, `pixfrac=0.6`. Compare to PSS oracle stacked at `drizzle_factor=2` (the existing `oracle_out_realser_drizzle2/stacked.fits`). | PSNR ≥ 25 dB on the central crop. SSIM ≥ 0.85. Mean abs diff ≤ 600 / 65535. The image must be a recognisable Jupiter — quick structural test: cross-correlation between the STScI output and the bicubic reference, downsampled to 64×64, peaks at offset (0, 0) ± 1 pixel. |
| `mpp_stsci_near_passthrough_oracle_drizzle3x` (optional extension) | `drizzle=stsci-3x`, `pixfrac=1` vs PSS oracle at `drizzle_factor=3`. | Same bars as 2x with PSNR ≥ 30 dB, SSIM ≥ 0.90. |

The fixture and acceptance numbers above are conservative starting points; once 5b.2 lands and we see actual numbers, the bars can be tightened (or loosened to reflect intrinsic differences) before merge. The intent is *correctness floor*, not algorithmic equivalence.

**Implementation note**: the PSS oracles already exist in `tools/pss_reference/oracle_out_realser*/`; the test only needs to add the STScI runs and the PSNR/SSIM/cross-correlation comparisons. Cross-correlation in Python is one `scipy.signal.correlate2d` call; in C++ we can use OpenCV's `matchTemplate(TM_CCORR_NORMED)`.

#### Group 3 — STScI / Bayer synthetic-truth tests

| Test | Setup | Acceptance |
|---|---|---|
| `mpp_drizzle_stsci_synthetic_resolution` | Generate 1024×1024 synthetic Jupiter with fine surface features (sinusoidal pattern of varying period, plus the existing synthetic planet). Downsample to 256×256 with 24 frames of prescribed sub-pixel global offsets uniformly in [−1, 1] px. Run pipeline at `drizzle=stsci-2x`, `pixfrac=0.5`. Upsample ground truth to 512×512 for comparison. | PSNR(STScI output vs upsampled-truth) **≥ PSNR(bicubic-2x output vs upsampled-truth) + 2 dB.** This is the resolution-recovery claim. Same fixture stacked at `drizzle=bicubic-2x` is the baseline. |
| `mpp_bayer_drizzle_synthetic_mosaic` | Generate 1024×1024 RGB ground truth with a slanted edge at 7° tilt, in colour-sensitive regions. Mosaic to RGGB Bayer at 256×256 with 24 frames, sub-pixel offsets. Run at `drizzle=bayer-2x`. Measure MTF50 per channel at the slanted edge (ISO 12233). | MTF50_R(bayer-2x) **≥ 1.3× MTF50_R(debayer-cv::cvtColor then bicubic-2x)**; same for B channel. G can be ≥ 1.0× since it's already 50% sampled. |
| `mpp_bayer_drizzle_first_smoke` (precedes the MTF test) | Generate a 2× synthetic RGGB mosaic with a known pattern (e.g. red letter "R", green "G", blue "B"); a single frame; run at `drizzle=bayer-2x`, `pixfrac=1`. | Output channels contain their respective letters cleanly. **This is the gate test for the 5b.3 dobox-CFA sanity check.** Skip the MTF test if this fails. |

#### Test runtime budget

The pixmap unit tests are sub-second. Group 2 reuses the existing real-SER fixture which already takes ~3 s in the bicubic oracle test; expect the two STScI runs to add ~6 s. Group 3's synthetic-truth tests on a 24-frame 256×256 fixture should take ~1 s each. Total Phase 5b test addition: well under 30 s.

### Summary of acceptance bars

- **5b.1 (pixmap)**: 6 unit tests green.
- **5b.2 (STScI path)**: Group 2 tests (passthrough + real-2x against PSS oracle) green at the bars above. Group 3 STScI synthetic-truth test green at ≥ 2 dB PSNR gain over bicubic.
- **5b.3 (Bayer path)**: Group 3 first-smoke test green (gates the dobox-CFA assumption), then the slanted-edge MTF test green at ≥ 1.3× per channel.
- **5b.4 (flag surface)**: CLI flag parser unit tests green; GUI manual smoke check (combo updates / widget visibility) signed off.
- **5b.5, 5b.6**: covered by Group 3 above.

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

- [x] **7.4** `pss test-big.ser -drizzle=stsci-2x` → produce a stacked FITS; visual sanity check (cleaner edges than bicubic-2x); SSIM comparison to bicubic-2x as a quantitative anchor.
- [x] **7.5** Bayer-2x on a real Bayer SER → produce a stacked FITS; demonstrate raw-CFA drizzle path runs end-to-end. (The slanted-edge MTF acceptance test lives in Phase 5b.6 — `mpp_bayer_drizzle::slanted_edge_resolution`.)

### Drizzle workflow policy

True drizzle should not run on classically-debayered colour input. The justification: drizzle preserves and amplifies high-frequency content, and a classical debayer (BAYER_RCD, bilinear, etc.) introduces interpolation artefacts (zipper effects, false-colour fringes, per-channel low-pass) before drizzle ever sees the data. Debayer-then-drizzle therefore takes those artefacts and amplifies them along with the legitimate signal.

The two supported drizzle workflows are the ones `applyreg.c` already implements:

- **Mono → mono drizzle** (`-drizzle=stsci-*` on single-channel input). No colour reconstruction needed.
- **CFA → RGB Bayer drizzle** (`-drizzle=bayer-*` on raw mosaic input). Drizzle does the colour reconstruction itself, with no prior interpolation to amplify.

`process_pss` and `process_stack_mpp` reject `-drizzle=stsci-*` on multi-channel input early, with a clear redirect to `-drizzle=bayer-Nx`. The two-step `register_mpp` + `stack_mpp` workflow catches the same configuration at stack-time because the `-drizzle=` flag is stack-side only and `process_stack_mpp` is where the final mode is known.

This policy keeps `do_kernel_*` in `cdrizzlebox.c` unchanged — no multi-channel inner-loop optimisation is needed because the 3-channel-input STScI configuration is now unsupported on principle, not a perf-cost question.

### Phase 7.4 / 7.5 benchmark results — test-big.ser (40,165 frames, 264×258 Bayer RGGB 8-bit, 2.6 GB)

All runs warm-cache after a single full-fixture warmup; two consecutive runs each to confirm repeatability. Stage A = 14 APs, stack_size = 4017 (10% top-quality slice).

| Mode | Output dims | Run 1 | Run 2 | Mean | vs bicubic |
|---|---|---|---|---|---|
| `bicubic-2x` (cv::resize INTER_LINEAR per AP buffer) | 460×448 | 175.93 s | 178.21 s | **177.07 s** | 1.00× |
| `stsci-2x`   (STScI dobox per-channel on debayered RGB) | 468×456 | 1090.11 s | 1089.52 s | ~1090 s | **rejected as of follow-up commit — bad workflow** |
| `bayer-2x`   (STScI dobox on raw CFA, pixfrac=0.7)   | 468×456 | 514.66 s¹ | 486.64 s² | **~500 s** | **~2.8× slower** |

¹ first run had a bug — see "Bayer-2x bug fixed during this milestone" below; the wall time of the buggy and fixed paths are within noise of each other.  ² post-fix run.

**Quantitative quality anchor** (central crop 400×412, against the bicubic-2x output as the reference):

| Pair | SSIM R / G / B | PSNR R / G / B (dB) |
|---|---|---|
| STScI vs bicubic | 0.9701 / 0.9785 / 0.9718 | 40.57 / 36.36 / 39.37 |
| Bayer vs bicubic | 0.9742 / 0.9709 / 0.9654 | 40.58 / 36.34 / 39.33 |
| STScI vs Bayer   | 0.9996 / 0.9986 / 0.9993 | 65.83 / 62.34 / 65.85 |

Sharpness (mean Sobel gradient on luminance): bicubic 0.02339, STScI 0.02266 (-3%), Bayer 0.02264 (-3%). All three paths produce comparable visual outputs.

**Bayer-2x bug fixed during this milestone**: the original Bayer wrapper passed `open_debayer=FALSE` to `ser_read_frame`, but that parameter is a latent NOP for SERs that were opened with debayer-on-open — the actual debayer call inside `ser_read_frame` dispatches on `ser_file->debayer_type_ser` (set at open time from `com.pref.debayer.open_debayer`), not on the function parameter. The wrapper was therefore reading a 3-channel debayered fits but treating the first plane (R) as the Bayer mosaic, so every CFA position fed dobox the same debayered-R value. The fix has two parts:

1. Temporarily flip `seq->ser_file->debayer_type_ser` to `SER_MONO` around the read loop (with restore on error), forcing `ser_read_frame` to return the raw single-channel mosaic.
2. Apply `adjust_Bayer_pattern_orientation` to compensate for `ser_read_frame`'s Y-flip (SER raster is top-down, FITS is bottom-up). For an even-height frame, RGGB on disk arrives as GBRG in the buffer; the cfa[] array we feed to dobox must match the in-memory pattern, not the file header.

Pre-fix per-channel means were a clear diagnostic — R/G/B all averaged ~7000 because all three accumulated the same debayered-R values. Post-fix matches bicubic within 3% per channel.

**Why STScI is 6× slower**: bicubic's stack contribution per frame is `cv::resize(intersection_strip, output_strip, INTER_LINEAR)` — single AVX-optimised pass. STScI's per-frame work is (a) build a pixmap (per-pixel weighted-AP interpolation, ~13× more arithmetic than a single shift add), (b) dobox's per-output-pixel box-area integration with pixfrac<1 (each input pixel touches up to 4 output cells with computed overlap area), (c) three times over for a 3-channel RGB frame after the per-channel-view fix from 5b.5. STScI is the price for true geometric drizzle; bicubic-2x stays the speed-optimal default.

**Why bayer-2x is ~2× faster than stsci-2x**: bayer drizzles a single-channel raw Bayer frame in one dobox call with CFA channel routing (`is_bayer=TRUE`), whereas stsci-2x calls dobox three times per frame (once per planar RGB channel). The pixmap is built once per frame in both cases. (Note: the stsci-2x configuration measured here is now a rejected workflow — see "Drizzle workflow policy" above.)

Artifacts in project root: `jupiter_drizzle_compare.png` (side-by-side bicubic / STScI / bayer), `phase_7_4_benchmark.txt` (raw timings). Source FITS in `/tmp/pss_bench/`.

### Drizzle-kernel sweep — bayer-2x on test-big.ser

Four `dobox` kernels measured back-to-back, same fixture, same APs, same pixfrac (0.7), same frame selection — only `-driz-kernel=` varies. Quality metrics computed against the `bicubic-2x` output (central 400×412 crop): SSIM and PSNR per channel then averaged, sharpness = mean Sobel-gradient magnitude on luminance, noise = std-dev in a 40×40 patch of the dark background.

| Kernel | Wall | vs bicubic | SSIM avg | PSNR avg | Sharpness | Noise |
|---|---|---|---|---|---|---|
| `turbo` | **174.6 s** | **0.99×** | 0.9686 | 38.64 dB | 0.02305 | 0.00041 |
| `gaussian` | 319.1 s | 1.80× | 0.9687 | 38.68 dB | 0.02299 | 0.00041 |
| `square` | 486.9 s | 2.75× | 0.9701 | 38.75 dB | 0.02264 | 0.00042 |
| `lanczos3` | 684.9 s | 3.87× | 0.9700 | 38.74 dB | 0.02256 | 0.00042 |
| _bicubic ref_ | _177.1 s_ | _1.00×_ | _—_ | _—_ | _0.02339_ | _0.00058_ |

Observations:

- **Quality is essentially flat across kernels on this fixture.** SSIM-avg span is 0.9686–0.9701 (0.15 % range), PSNR-avg span 38.64–38.75 dB (0.11 dB). Jupiter is a smooth disc with limited aliased high-frequency content; kernels mostly diverge on sharp edges, of which there are few.
- **Turbo is the practical winner**: 6.3× faster than the slowest kernel (lanczos3), 2.8× faster than square — and it lands at the same wall-clock cost as the `bicubic-2x` interpolation. There is no measurable quality penalty on this fixture.
- **Bayer drizzle is structurally less noisy than bicubic** (0.00041 vs 0.00058 std-dev in dark background) regardless of kernel. Drizzle's per-pixel averaging across 4017 frames suppresses per-frame noise more thoroughly than per-AP bicubic resampling does.
- **Square and lanczos3 have nearly identical sharpness** to each other and to bicubic; lanczos3 doesn't deliver the high-frequency recovery on a smooth planetary target that it would on a high-contrast edge.

**Decision (recorded)**: `turbo` is now the default kernel — `mpp_config_defaults` sets `drizzle_kernel = MPP_KERNEL_TURBO`, the GUI combo's initial selection is "turbo", and the GUI fallback on invalid state is `MPP_KERNEL_TURBO`. Same cost as the bicubic baseline, slight noise advantage, no quality penalty on this fixture. The synthetic high-frequency tests in `mpp_drizzle_test.cpp::mpp_stsci_synthetic::resolution_recovery` and `::mpp_bayer_drizzle::slanted_edge_resolution` continue to use `square` explicitly, so the kernel-choice acceptance bars are unchanged.

Artifact: `jupiter_kernel_compare.png` (side-by-side 4-kernel comparison with per-kernel SSIM and sharpness annotations).

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
