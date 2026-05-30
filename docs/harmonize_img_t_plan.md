# Plan: Harmonize da3d & nlbayes onto `img_t` (branch `harmonize_img_t`)

## Goal
Replace da3d's `Image`/`DftPatch` and nlbayes's `std::vector<float>`+`ImageSize` with the
deconvolution `img_t<T>` template, extending `img_t`/`img_expr_t` where they lack an operation
these filters need — **without changing denoising output (within tolerance) or regressing
performance**. `LibMatrix` is kept as a standalone linear-algebra unit, but made generic, faster,
and connected copy-free to `img_t`.

---

## Key findings that shape the plan

| Module | Image type | Layout | FFT | Special structures |
|---|---|---|---|---|
| **deconvolution** | `img_t<T>` | **planar** (channel-major, `data[x + w*(y + h*dd)]`) | batched multi-ch `fftwf_plan_many`, **plan cache keyed by (h,w,d)** | `img_expr_t` expression templates, `process_in_slices` (symmetric reflect), `rgb2ycbcr` |
| **nlbayes** | `std::vector<float>` + `ImageSize` | **planar** (`c*wh` offsets, LibImages.cpp:562) | **none** (`fftw3.h` is a dead include) | `LibMatrix` = patch covariance/inverse/GEMM (algorithm, *not* image ops) |
| **da3d** | `da3d::Image` | **interleaved** (`row*cols*ch + col*ch + chan`, Image.hpp:79) | `DftPatch` = per-patch batched r2c/c2r | `WeightMap` pyramid for patch selection |

Decisive consequences:
- **nlbayes is already planar → lower-risk migration.** Its `LibMatrix` patch algebra stays out of `img_t` scope.
- **da3d is interleaved → migration changes channel layout.** This mirrors the deconvolution planar fix and
  likely *exposes/fixes a latent colour bug* (note the contradictory `call_nlbayes.cpp:222` "convert bgrbgr back
  to planar" vs the flat planar copy at `:98`). Resolving the exact layout contract is a prerequisite.
- **`DftPatch` is the hard part.** `img_t` can do batched multi-channel FFT with a plan cache, but DA3D transforms
  thousands of small patches; a naïve `img_t<complex>` per patch with `FFTW_MEASURE` would be catastrophically slow.
  A thin patch-DFT abstraction reusing cached plans must be spiked and benchmarked first.
- **`WeightMap`** is algorithm-specific; keep it as a class, back it with `img_t`.

---

## Risks / decisions to confirm before coding
- [x] **Layout contract at the da3d boundary** — RESOLVED. `do_nlbayes` (call_nlbayes.cpp:92-101) fills `bgr_f`
      by a flat copy of the planar `fit->fdata`/`fit->data`, so the buffer is **planar** `[R…R,G…G,B…B]`. nlbayes
      consumes it as planar (correct). da3d builds `Image input(bgr_f,…)` (line 199) but `da3d::Image` indexes
      **interleaved** (Image.hpp:79), and `ColorTransform` (DA3D.cpp:58-63) reads `val(col,row,0/1/2)` as a pixel's
      R/G/B. **For colour images the DA3D stage misinterprets planar data as interleaved → a latent colour bug** (same
      class as the deconvolution one; the line-222 "convert bgrbgr back to planar" comment is stale — no de-interleave
      happens). Mono is unaffected. **Consequence:** the da3d migration is a *bug-fix* for colour, so da3d colour golden
      output will intentionally differ; capture before/after.
- [ ] **Exact RGB↔YUV coefficients** — nlbayes `transformColorSpace` (LibImages.cpp:465) uses a specific opponent
      transform; `img_t::rgb2ycbcr` may differ. Must match nlbayes math exactly, not "a" YUV.
- [ ] **Tolerance policy** — bit-exact is unrealistic (FFT/threading reorderings). Propose `max abs err ≤ 1e-5` on
      [0,1] data for golden regression.
- [ ] **Scope of `img_t`/`img_expr_t` extension** — agree the new API surface before implementing.

---

## Phase 0 — Setup & de-risking
- [x] Create branch `harmonize_img_t` from `master`.
- [ ] **Golden baselines:** capture current da3d and nlbayes output on a fixed mono + a fixed colour test FITS
      (small, committed fixtures), at fixed sigma/params. These reference arrays are the regression oracle for the
      whole effort. *(Deferred to the start of the migration phases, where there is a refactored version to compare
      against; needs real FITS fixtures + a built Siril.)*
- [x] Add a Criterion test target to `src/tests/meson.build` (model on `imoper_test.c`), gated on
      `criterion_dep.found()`. *(Done: `harmonize_img_t_test.cpp`, suite `harmonize`. NOTE: clang is absent in the dev
      container; configured a gcc meson build dir `_build_gcc` with `-Dcriterion=true` for building/running tests.)*
- [x] **Resolve the layout-contract question** (see Risks above; da3d colour path is a latent bug).
- [x] **Spike: patch-DFT performance.** RESOLVED by code analysis: `img_t`'s `fft::r2c` runs a *full complex* FFT
      (fft.hpp:117-133) and heap-allocates a new buffer per call, while DA3D needs the real-optimized half-spectrum
      transform with persistent buffers. **Decision: dedicated helper.** Implemented `imgops::dft_patch`
      (`image_dft.hpp`) backed by planar `img_t` storage using `fftwf_plan_many_dft_r2c/c2r` with reused plans.

**Unit tests:** harness smoke test; baseline-capture reproducibility (same input → same bytes twice).

---

## Phase 1 — Extend `img_t` / `img_expr_t` to cover the gaps (TDD: test first)
Write the Criterion test *before* each implementation, asserting against a hand-computed or current-code reference.

- [x] **Symmetric (half-sample) boundary pad / unpad** — `imgops::pad_symmetric`/`unpad` + `symmetric_coordinate`
      (image_ops.hpp). Confirmed da3d `SymmetricCoordinate` and nlbayes `symetrizeImage` use the *same* half-sample
      reflection, so one primitive serves both. Tests: reflection cases, interior preserved, pad→unpad identity.
      *(Byte-for-byte match vs nlbayes `symetrizeImage` to be asserted during Phase 2.)*
- [ ] **RGB↔YUV transform** — DEFERRED to the per-filter migrations: nlbayes `transformColorSpace` and da3d
      `ColorTransform` (orthonormal opponent: `(r+g+b)/√3, (r-b)/√2, (r-2g+b)/√6`) are *different* transforms, so each
      becomes a thin img_t-based function in its own module rather than a shared Phase-1 op.
- [x] **Tile split / merge** — `imgops::compute_tiling`/`split_tiles`/`merge_tiles` (image_ops.hpp), faithful port of
      da3d `ComputeTiling`/`SplitTiles`/`MergeTiles`. Tests: tiling cases, split→merge round-trip reconstructs original.
- [x] **Weighted aggregation helper** — folded into `merge_tiles` (overlap-add accumulate + divide by summed weight).
      Test: direct weighted-average (value/weight) check.
- [x] **Monochrome detect / collapse** — `imgops::is_monochrome` + `imgops::make_monochrome` (double-accumulating mean,
      matching da3d). Tests: equal-channel detection, colour mismatch, channel-average value.
- [x] **Patch-DFT helper** — `imgops::dft_patch` (image_dft.hpp); see Phase 0 spike. Tests: DC component, forward→inverse
      identity. *(Equality vs da3d `DftPatch` on a fixed patch to be asserted during Phase 3.)*
- [x] **Axis-reduce + axis-broadcast for `img_expr_t`.** Added `reduce_axis(expr, axis, reductor, init)` and
      `broadcast_axis(expr, axis, n)` (image_expr.hpp) over `AXIS_X/AXIS_Y/AXIS_D`, plus `imgops::mean_axis`. NL-Bayes
      centering is now `group - broadcast_axis(mean_axis(group, AXIS_X), AXIS_X, group.w)`; tested incl. the centering
      composition. **Side fix:** constrained the `DEFINE_EXPR_1/2` arithmetic operators with a new `is_imglike_v` trait
      — they were unconstrained and matched `std::vector` iterator arithmetic, breaking any TU that mixes
      `image_expr.hpp` with `std::vector<img_t<T>>` (required by the tiling ops).
- [ ] Confirm/parameterise threading (`com.max_thread`/`fftw_max_thread`) for the new ops.

### Tiling review — `img_t::process_in_slices` vs `imgops` tiling
The deconvolution `img_t::process_in_slices` (image.hpp) and the new `imgops::split_tiles`/`merge_tiles` look similar
but are different abstractions, differing on three axes:
1. **Sizing:** process_in_slices picks slice sizes from a *memory budget* (smallest-number / FFTW-optimal /
   best-compromise, `com.pref.max_slice_size`); imgops takes a *fixed tile count* (`compute_tiling`).
2. **Reflection:** process_in_slices uses *whole-sample* reflection; imgops uses *half-sample* (to match da3d/nlbayes).
3. **Stitch:** process_in_slices discards the overlap and hard-cuts the non-overlap core into the output; imgops does a
   *weighted overlap-add blend*.

Decision: a forced merge would be a flag-heavy god-function (3 mode axes) that is *more* complex, not less — so it is
**not** done. The one genuinely-duplicated piece, the boundary reflection math, IS now shared via `image_boundary.hpp`
(`symmetric_coordinate` + `reflect_whole_sample`), used by both. Full unification of the tiling *framework* is deferred
to Phase 2/3 (rule of three): once nlbayes (`subDivide`/`subBuild`) and da3d (`SplitTiles`/`MergeTiles`) are migrated
onto `imgops`, there will be three concrete consumers and any genuinely-common shape beyond the reflection can be
extracted with evidence rather than speculation.

**Exit criterion:** all Phase-1 ops have green Criterion tests; no migration started yet.

---

## Phase 2 — Migrate nlbayes (lower risk, already planar)
- [ ] Replace `std::vector<float>`+`ImageSize` with `img_t<float>` through `NlBayes.cpp`/`LibImages.cpp`, routing
      boundary/tiling/colour through Phase-1 ops. Keep `LibMatrix` algebra untouched.
- [ ] Update entry point `do_nlbayes` (call_nlbayes.cpp) to build an `img_t` directly from `fit->fdata` (planar, no
      copy/reorder).
- [ ] Delete `ImageSize` and now-dead `LibImages` functions superseded by `img_t` ops.
- [ ] **Golden regression:** nlbayes output matches Phase-0 baseline within tolerance (mono + colour).
- [ ] Per-function Criterion tests for any nlbayes logic that changed shape (patch extraction indices, aggregation).

### Section 2.5 — LibMatrix integration, genericization & OMP optimization
LibMatrix stays a standalone linear-algebra unit (GEMM / Cholesky-inverse / covariance / transpose). We improve how it
*connects* to `img_t`/`img_expr_t` and how it *performs*, without folding its algebra into the expression templates
(wrong access pattern: contractions + sequential factorization, and the mandated patch-contiguous layout is orthogonal
to the one axis `img_expr_t` can reduce).

**Central constraint:** NL-Bayes already parallelizes at the **sub-image level** — `#pragma omp parallel for` over
`nbParts = 2*nbThreads` tiles (NlBayes.cpp:252, 317), one thread per sub-image. Every LibMatrix call therefore runs
**nested inside an already-parallel region.** This dictates the OMP strategy below.

**2.5.1 — Container-agnostic API (`std::span`)**
- [ ] Replace every `std::vector<float>&` parameter with `std::span<T>` (dims already passed explicitly).
  - **Why:** an `img_t`-owned buffer (`std::vector<T, fftw_alloc<T>>`) can't bind to `std::vector<float>&`, forcing a
    copy on every group→LibMatrix handoff. `std::span` makes it zero-copy.
  - Tests: call each function via a `span` over an `img_t` buffer and over a raw array; identical results to the
    current `vector` API on fixed inputs.

**2.5.2 — Template on `T` for precision parity with `img_t<T>`**
- [ ] Template all four functions on `T`. Allow `covarianceMatrix`/`inverseMatrix` to accumulate in `double` even when
      storage is `float` (covariance is where NL-Bayes is most precision-sensitive).
  - Tests: `float` instantiation matches current output within ε; `double`-accumulation reduces error on an
    ill-conditioned covariance fixture; `A·inverse(A) ≈ I` for SPD `A` (both `T`).

**2.5.3 — Represent the patch group as `img_t<T>(nSimP, sP2)`**
- [ ] Store `io_group3d[c]` as `img_t<T>(nSimP, sP2)` (patch index = x/contiguous, pixel = y), mapping **bit-for-bit**
      to LibMatrix's `i*nSimP + k` layout (LibMatrix.cpp:143). Same fftw-aligned buffer feeds `span`-based LibMatrix with
      no copy, and exposes the elementwise estimator steps to `img_expr_t`.
- [ ] Route the elementwise update `group -= valDiag · groupT` (NlBayes.cpp:793–795) through `img_expr_t`:
      `group.map(group - valDiag * groupT)`.
- [ ] Leave `centerData` + baricenter add/sub (a reduce-over-x then broadcast-over-x) as hand loops **unless** the
      optional Phase-1 axis primitives are built.
  - Tests: golden regression — refactored estimator matches current NL-Bayes output within tolerance (mono + colour).

**2.5.4 — OMP optimization (no intrinsics; multi-platform)**
- [ ] **Primary — `#pragma omp simd reduction(+:...)` on the inner dot-products** (compose with outer threading; pure
      SIMD within the current thread):
  - `covarianceMatrix` inner `k` loop (LibMatrix.cpp:142–144)
  - `productMatrix` inner `k` loop (179–181)
  - `inverseMatrix` three inner reduction loops (59–61, 74–76, 84–86)
- [ ] **Canonicalize pointer-walk loops to index form first** (e.g. inverse's dual-incrementing `r`/`s` →
      `for k: z += mat[base1+k]*mat[base2+k]`) so the vectorizer accepts the `simd` pragmas.
- [ ] **No unconditional thread-level `parallel for` inside LibMatrix** — the caller already saturates threads
      (`nbParts = 2·nbThreads`). Where a serial-context caller could benefit (e.g. a future da3d use), gate thread-level
      parallelism on spare capacity using `img_t`'s pattern: `available = com.max_thread - omp_get_num_threads(); if
      (available > 1)` or `!omp_in_parallel()`. Applies to `covarianceMatrix`'s outer triangular `i` loop
      (`schedule(guided)` for the triangular imbalance) and `productMatrix`'s outer `i` loop (with `q0` made
      **private/thread-local** — currently shared, would race).
- [ ] **`inverseMatrix` stays scalar-outer** (Cholesky columns are sequential); only inner reductions get `simd`. For
      small `sP2` the gain is modest — measure before keeping.
  - Tests: results bit-stable vs pre-optimization (within FP-reassociation ε for the reductions); determinism across
    runs/thread counts; **micro-benchmark** at representative sizes (`sP2 ≈ 25–49`, `nSimP ≈ 30–100`), recorded
    before/after — fail on regression.

**2.5.5 — Sequencing**
- [ ] 2.5.1 → 2.5.2 → 2.5.4 land independently of the image migration (pure LibMatrix, fully unit-testable). 2.5.3
      depends on Phase-2 `img_t` adoption. Do the LibMatrix-internal work (span + template + OMP) **first** as a
      self-contained, golden-tested PR, then wire it to `img_t`.

---

## Phase 3 — Migrate da3d (interleaved→planar, DftPatch→patch-DFT)
- [ ] Replace `da3d::Image` with `img_t<float>`; convert all `val(col,row,chan)` to planar `operator()` — removing the
      interleaved assumption.
- [ ] Replace `DftPatch` with the Phase-0/1 patch-DFT approach; port freq-domain ops (DA3D.cpp:331–342: modulus via
      `std::norm`, scalar multiply, DC special-case).
- [ ] Re-back `WeightMap` (WeightMap.cpp pyramid) with `img_t`; keep its min-pyramid algorithm.
- [ ] Route tiling through Phase-1 split/merge; port `ColorTransform`, `ComputeRegressionPlane`, `BilateralWeight`,
      `ModifyPatch` to planar/per-channel.
- [ ] Update entry point (call_nlbayes.cpp:199) to construct `img_t` from the (layout-confirmed) buffer.
- [ ] **Golden regression** within tolerance (mono + colour); if colour now differs because the old interleaved path
      was buggy, document it as an intended fix with before/after images.
- [ ] **Performance benchmark** vs baseline; investigate any >5% regression.

---

## Phase 4 — Cleanup & validation
- [ ] Delete dead files/types: da3d `Image.hpp`, `DftPatch.hpp`; nlbayes `LibImages.*` remnants, `ImageSize`; dead
      `fftw3.h` include.
- [ ] Update `meson.build` source lists for both filters.
- [ ] Full build under clang (warnings-clean) + gcc.
- [ ] Run the full Criterion suite + both golden regressions + perf benchmarks; record results in the PR.
- [ ] Manual end-to-end: run NL-Bayes and DA3D from the Siril UI on a real colour image; confirm sane output.

---

## Cross-cutting testing strategy
- **Unit (Criterion, `src/tests/`):** one test file per new `img_t`/`img_expr_t` op and per migrated sub-routine,
  modelled on `imoper_test.c`/`soper_test.c`; assert against current-code reference or analytic values.
- **Golden regression (the real safety net):** fixed input FITS → output arrays captured on `master` before any change;
  every phase re-asserts ≤ tolerance. This is what guarantees "denoising results don't change."
- **Performance:** committed micro-benchmark for patch-DFT (DA3D hot path) and LibMatrix ops, plus whole-image timing for
  both filters; fail the phase on unexplained regression.

---

## What stays out of scope
- `LibMatrix`'s core algebra is **not** folded into `img_expr_t` (contractions + sequential factorization don't fit the
  elementwise/stencil model; the patch-contiguous layout is orthogonal to `reduce_d`'s axis).
- `LibImages` matrix-adjacent helpers that are pure linear algebra stay with LibMatrix.
- No CPU intrinsics anywhere (multi-platform) — optimization is `-O3` + OpenMP (`simd`, and threading only where not
  nested).
