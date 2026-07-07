# mpp_improve — closing the detail gap vs AutoStakkert

## Problem

On the 2026-07-07 Saturn L data, AutoStakkert (ap198) preserves globe banding
that the Siril mpp stack (half-box 24, 50 %) cannot deliver no matter how it
is sharpened. Both outputs have Saturn at the same pixel scale, so the gap is
algorithmic, not sampling. The mpp chain is bit-near-equivalent to PSS
(128 dB gauge), so the gap is PSS's ceiling, and it concentrates in **per-AP
shift measurement accuracy on low-contrast structure**: every tenth of a
pixel of per-frame per-AP shift error convolves the stack with that blur
kernel. Saturn's belts are low-amplitude high-frequency signal on a smooth
limb-darkened gradient — exactly what a ~0.2–0.5 px RMS shift-noise floor
erases before post-processing ever sees it.

Root causes, ranked (from the July 2026 code review):

1. **Blurry correlation reference.** Stage B correlates each (blurred) frame
   box against the *seeing-averaged* mean frame. Banding is washed out of
   that reference, so it cannot contribute to the peak; correlation peaks on
   globe APs are broad and the sub-pixel fit rides noise.
2. **`TM_CCORR_NORMED` without mean subtraction.** Plain NCC is
   gain-invariant but not offset-invariant: brightness offset mismatches
   between a frame and the averaged reference (transparency drift, haze —
   Stage B applies no brightness normalisation to the correlation inputs)
   bias the peak by pixels. *(Empirical correction during implementation:
   a pedestal or gradient common to both windows cancels in the
   normalisation and does NOT bias plain NCC — the demonstrated failure
   regime is the gain+offset mismatch, locked in by the
   `zero_mean_recovers_shift_under_gradient` fixture.)*
3. Heavy pre-blur (ksize 7, twice in phase 1) + 3×3 parabolic sub-pixel fit
   → pixel-locking bias, ~0.1–0.2 px accuracy floor even at good APs.
4. AP pitch hard-coupled to box size (step = 2.25 × half-box → 54 px at
   half-box 24, ~60 APs where AS! used 198): warp field sampled ~3× coarser,
   intra-patch differential warp averaged into blur.
5. Per-AP frame ranking quantized to uint8 (`convertScaleAbs` α=1/256) on
   stride-2 blurred data → per-AP "best 50 %" selection is weakly informative
   exactly where the detail is.

## Plan

### Phase 1 (this branch, first commits): items (1) + (2)

**(2) Zero-mean Stage B correlation** — config-gated switch of the per-AP
`matchTemplate` method from `TM_CCORR_NORMED` to `TM_CCOEFF_NORMED`
(zero-mean NCC = Pearson correlation, invariant to frame-vs-reference
brightness gain and offset mismatch), applied in both correlation phases of
Stage B only. Global alignment and the Phase-4 `shift_compute_all` path keep
`TM_CCORR_NORMED` (PSS-faithful; used by the equivalence oracle tests).

- `cfg.alignment_points_zero_mean`, default **true** (deliberate divergence,
  same precedent as the Planet-mode default; documented in
  `MPP_PSS_DIFFS.md`). CLI opt-out `-no-zero-mean` on `mpp`/`register_mpp`.

**(1) Refined-reference second registration pass** — after the normal Stage B
pass, build a small, *locally registered* stack and re-measure all per-AP
shifts against it, replacing the seeing-blurred mean frame as the correlation
reference:

1. Pass 1: Stage B as today (shifts vs mean frame).
2. Reference build: stack the top **K** frames *per AP* (existing per-AP
   ranking, hard cut, weight 1) with the classical engine at scale 1.0,
   mono analysis channel, using the pass-1 shifts. Merge, then paste into a
   mean-frame-sized canvas (mean frame fills the trimmed border ring) so AP
   box coordinates stay valid.
3. Pass 2: rerun `stack_compute_shifts_streamed` with reference boxes cut
   from that stacked canvas; the refined shifts replace pass-1 in
   `run->shifts` and flow into the sidecar and Stage C unchanged.

**K (first-stage stack size): as small as gives a useful reference.**
Reference noise adds in quadrature to frame noise in the NCC, so K ≥ ~8
makes the reference contribution negligible (≤ 1/8 of the variance), while
*small* K keeps the reference sharp (only the per-AP best moments of seeing,
and less residual pass-1 shift jitter averaged in). Auto value:
`K = clamp(ceil(5 % of included frames), 8, 32)`, capped by the baked per-AP
selection; explicit override `-refine-frames=N`
(`cfg.alignment_points_reference_frames`, 0 = auto).

- `cfg.alignment_points_refine_reference`, default **true**; CLI opt-out
  `-no-refine`.
- Only frames appearing in some AP's top-K (plus the global top-K for the
  background composite) are read for the reference build, so the extra cost
  is roughly one partial read pass + one extra Stage B pass (analysis-cache
  hits in the common workflow).
- Cancellation propagates (MPP_EINTR); reference-build failure (OOM etc.) is
  non-fatal — pass-1 shifts are kept with a log message.

Also in phase 1: sidecar version bump (raw `mpp_config_t` snapshot changes
layout), CLI flag plumbing + command help, `MPP_PSS_DIFFS.md` sections, unit
test for the zero-mean flag, defaults pinned to the PSS-faithful values in
any test asserting oracle equivalence.

Validation: A/B on the Saturn SER and on `test-big.ser` — compare
`-no-refine -no-zero-mean` (baseline = current behaviour, must stay
bit-identical) vs defaults; inspect banding with the à trous level-1 ×40
boost standard. Watch Stage B failure counters and the shift-viewer spread
(per-AP shift dispersion should tighten in pass 2).

### Phase 2 (next): measurement accuracy floor

- Upsampled-correlation sub-pixel refinement (Guizar-Sicairos-style local
  DFT ×8–16, or 5×5 centroid) replacing the 3×3 quadratic fit, ideally on
  less-blurred data once phase 2 has localized. Config-gated.
- Float-precision per-AP ranking (σ on CV_32F |Laplacian|, skip the uint8
  quantization) behind the same PSS-equivalence gate.

### Phase 3 (later): warp-field stacking

- Decouple AP pitch from box size (overlap factor), or better: interpolate a
  smooth per-pixel warp from the AP shift lattice and remap each frame once
  (Lanczos), replacing rigid-patch pastes + Hann blending. Builds on the
  pss-derot warp-field engine. Removes intra-patch warp averaging and blend
  softening; enables AS!-like effective AP density.

### GUI (follows once CLI behaviour is validated)

- Registration tab: "Zero-mean correlation" and "Refine reference (second
  pass)" checkbuttons + a reference-frames spin (0 = auto). Defaults on.
  (No "PSS" in user-facing strings.)
