# PlanetarySystemStacker port to Siril — implementation plan

Tracker for porting the PSS multipoint registration + stacking pipeline into Siril on the `pss` branch. Tick boxes as work lands. Open items are status-of-record; revise this file when scope changes.

## Guiding principles

- **Algorithmic fidelity over re-invention.** Where PSS uses an OpenCV call, the port calls the same OpenCV function with the same arguments. C/C++ throughout, no Python in the shipped product.
- **Golden-output regression testing.** Each module is validated against intermediate artifacts captured from a Python run of PSS on the same input. Acceptance bar: per-frame quality within 1e-4 relative, integer shifts exact, sub-pixel shifts within 0.05 px, final stacked PSNR > 50 dB vs PSS reference.
- **Headless first.** Everything driven by a new `pss` command. GUI gated behind an explicit go/no-go pause.
- **Algorithmic-match definition:** behaviour-equivalent within float-rounding. Bit-identical is not achievable across platforms even between two Python runs.

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
- [ ] 0.1.8 Reserve one real planetary SER for the Phase-7 end-to-end milestone. _(User to nominate file or maintainer to pick well-known sample.)_

**Exit criterion:** PSS produces oracle artifacts on the synthetic SER; Siril compiles with empty mpp skeleton; `pss` returns clean "not implemented."

---

## Phase 1 — Quality ranking

- [x] 1.1 Implement `mpp_rank.cpp`: PSS pipeline = GaussianBlur(7×7) → stride-2 sub-sample (`align_frames_sampling_stride`, NOT `rank_frames_pixel_stride` — the latter is for xy/Sobel only) → `cv::Laplacian(CV_32F)` → `cv::convertScaleAbs(alpha=1/256)` → `meanStdDev[1][0][0]`. Returns σ on the **uint8** Laplacian, matching PSS's `Frames.frames_mono_blurred_laplacian` + `RankFrames.frame_score` chain exactly.
- [x] 1.2 Optional brightness normalisation: `σ / (cv::mean(threshold(mono, thr, 255, THRESH_TOZERO))[0] + 1e-10)` with `thr = frames_normalization_threshold * 256` for 16-bit input.
- [ ] 1.3 Wire frame iteration through `mpp_frames` (SER reader from `src/io/ser.c`, optional debayer via `cv::cvtColor`). _(Deferred: stub still returns ENOTIMPL; the algorithm is validated, plumbing through Siril sequence types can ride alongside Phase 2/3 to amortize the work.)_
- [x] 1.4 Oracle test: rank synthetic dataset; per-frame quality within 1e-4 relative of PSS. _(`mpp_rank::oracle_equivalence_synthetic` — passes.)_
- [x] 1.5 Sanity tests: blur-decreases-quality, brightness-invariance. _(Plus default-config, score-positive, threshold-scales-with-bitpix → 6 sanity tests, all green.)_

## Phase 2 — Global frame alignment

- [ ] 2.1 Implement `mpp_align.cpp` patch picker (`align_frames.compute_alignment_rect` + `quality_measure_threshold_weighted`). Bounds 0.2–0.7 of frame, border 10 px, black threshold 10240.
- [ ] 2.2 Two-phase correlation: phase 1 stride-2 grid, `cv::TM_CCORR_NORMED`, search half-width `(align_frames_search_width − 4)/2`; phase 2 stride-1, search ±4.
- [ ] 2.3 Build average reference frame from top `align_frames_average_frame_percent=5%`.
- [ ] 2.4 Oracle test: per-frame integer shifts exact; patch coordinates exact.

## Phase 3 — Alignment-point grid

- [ ] 3.1 Implement `mpp_ap.cpp` staggered grid (`step_size = (half_patch_width × 4.5) / 3`, `half_patch_width = half_box_width × 1.5`, `half_box_width=24`).
- [ ] 3.2 Per-AP filters: brightness threshold (`10 × 256` in 16-bit), structure threshold (`0.04 × max_σ_in_frame`).
- [ ] 3.3 Oracle test: AP coordinate set exact match.

## Phase 4 — Per-AP local shifts

The single highest-risk module. Most numerical landmines live here.

- [ ] 4.1 Two-phase multilevel correlation in `mpp_shift.cpp`:
    - Phase 1: stride-2 grid, Gaussian-blurred reference and frame windows, `cv::TM_CCORR_NORMED`, integer shift = `2 × ((search_width_phase1, search_width_phase1) − maxLoc)`.
    - Phase 2: stride-1, search ±4, no blur, `TM_CCORR_NORMED`.
- [ ] 4.2 Optional sub-pixel: fit `f(x,y) = a + bx + cy + dx² + ey² + fxy` over 3×3 correlation neighbourhood; solve 2×2 linear system. Gated on `alignment_points_local_search_subpixel=true`.
- [ ] 4.3 Penalty matrix for off-centre peaks per `alignment_points_penalty_factor=0.00025`.
- [ ] 4.4 Sign-convention regression tests: known synthetic shift recovered.
- [ ] 4.5 Oracle test: per-AP per-frame shifts match within 0.05 px (sub-pixel) / exact (integer).

## Phase 5a — Bicubic stacking (PSS-faithful)

- [ ] 5a.1 Factor `mpp_stack.cpp` into per-frame resample stage and AP composition stage. Pluggable resample backend interface.
- [ ] 5a.2 Implement `mpp_resample_bicubic()`: `cv::resize(frame, scale=N, INTER_CUBIC)`.
- [ ] 5a.3 Implement `mpp_compose_aps()`: top-N per AP via local quality, sub-pixel patch extraction via `cv::getRectSubPix`, ramped 2D weight `min(y_ramp, x_ramp)`, accumulate Σ patch × weight and Σ weight, normalise at end.
- [ ] 5a.4 Per-frame brightness equalisation (mean match).
- [ ] 5a.5 Background composition: if AP coverage > `1 − background_fraction` (0.3) use only APs; else compose from non-AP regions and blend at AP edges with ramp width `background_blend_threshold` (0.2).
- [ ] 5a.6 Drizzle 1.5/2/3 via output buffer dimensions; PSS-faithful path.
- [ ] 5a.7 Oracle test: synthetic SER end-to-end; per-channel PSNR > 50 dB vs PSS reference.

## Phase 5b — STScI drizzle + Bayer drizzle (Siril-original)

- [ ] 5b.1 Implement `mpp_drizzle.cpp`: build `driz_param_t` per frame from `mpp_shift_t`. Decide whether per-AP residuals go through `dobox()` as a piecewise-affine map or are absorbed in the composition stage — resolve once `dobox()` non-affine handling is verified.
- [ ] 5b.2 STScI-on-debayered path: `mpp_resample_stsci()`. Default pixfrac 0.7, kernel "square"; expose `-pixfrac` and `-driz-kernel`.
- [ ] 5b.3 Bayer drizzle path: analysis stages still operate on debayered frame, but stack-time `mpp_resample_bayer_stsci()` feeds raw Bayer samples into a 3-channel output. Per-channel weight tracking and normalisation.
- [ ] 5b.4 `pss` flag surface: `-drizzle=Off|1.5|2|3` (bicubic), `-drizzle=stsci-2x|3x`, `-drizzle=bayer-2x|3x`; `-pixfrac=<f>`, `-driz-kernel=<name>`.
- [ ] 5b.5 Synthetic-truth test: downsample known high-res ground truth with prescribed sub-pixel offsets, assert STScI path recovers it with higher PSNR than bicubic path.
- [ ] 5b.6 Bayer-drizzle test: mosaic a high-res RGB ground truth, run bayer-drizzle; assert color resolution exceeds debayer-then-stack on a slanted-edge target.

## Phase 6 — Orchestrator and `pss` command

- [ ] 6.1 Implement `register_mpp()` in `mpp.c`: chain load → rank → global-align → place-APs → per-AP shifts → stack → write FITS. Progress callbacks via Siril's existing infra.
- [ ] 6.2 Implement `process_pss()` in `command.c` with the flag surface above. PSS defaults; sequence is the currently-loaded one.
- [ ] 6.3 Register `REG_PSS` in `regmethod_index` (`registration.h`) so `register -method=pss` also works.
- [ ] 6.4 Command help text.

**Exit criterion:** `pss` runs the synthetic SER end-to-end from the CLI and produces a stacked FITS comparable to the oracle.

## Phase 7 — ★ Real-SER end-to-end milestone ★

- [ ] 7.1 Run PSS (Python) on a real planetary SER with defaults; archive output.
- [ ] 7.2 `pss` (Siril) default → compare to PSS reference. PSNR > 50 dB.
- [ ] 7.3 `pss` `-drizzle=2` (bicubic) → compare to PSS `drizzle_factor=2` reference. PSNR > 50 dB.
- [ ] 7.4 `pss` `-drizzle=stsci-2x` → sanity-check: clean image, comparable or better SSIM vs bicubic-2x output.
- [ ] 7.5 If Bayer SER: `pss` `-drizzle=bayer-2x` → slanted-edge MTF measurement vs debayer-then-stack.
- [ ] 7.6 Document any deliberate departures from PSS in `src/registration/MPP_PSS_DIFFS.md`.

## Phase 8 — Hardening

- [ ] 8.1 Edge cases: tiny sequences (< 30 frames), mono-only, color without SER header info, oversaturated targets, sequences with junk frames, memory profile on large sequences.
- [ ] 8.2 Error paths return clean status codes through the orchestrator. No asserts on user-reachable paths.
- [ ] 8.3 Performance pass: OpenMP scheduling, buffer reuse, profile on real SER. Aim for parity-or-better vs PSS Python runtime.

---

## 🛑 GUI go/no-go pause

Algorithm shippable headless at this point.

- [ ] **Decision recorded:** Option A (wait for GTK4 merge) / Option B (GTK3 on trunk, port later).

---

## Phase 9 (gated) — GUI port

- [ ] 9.1 Parameter configuration dialog (mirrors `parameter_configuration.ui`).
- [ ] 9.2 Frame selector + viewer (`frame_selector_gui.ui`, `frame_viewer_gui.ui`).
- [ ] 9.3 Alignment-point editor (`alignment_point_editor_gui.ui`) — most complex panel.
- [ ] 9.4 Shift distribution viewer (per-AP shift vector field).
- [ ] 9.5 Job dialog for batch multi-SER processing.

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
