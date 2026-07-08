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

### Phase 1 results (Saturn A/B, sat.ser 512×320×17934 8-bit, 1.5×, 50%)

- hb 24 (26 APs): refined+zero-mean vs baseline is a wash on the globe —
  pixel-scale ratio 0.995–1.016, band-scale (6–32 px) +0.7–1.4 %, visually
  marginal, no artefacts, +18 s. Stage B failures 3.5 % → 5.3 % of pairs.
- Density probe via `-half-box=12` (76 APs): **net loss** — band-scale
  −0.3…−2.3 %, failures 10 %/19 %, visible AP-lattice imprints on the globe.
  Shrinking the half-box starves the 24 px correlation template; density
  gained that way is worthless. AS! used 198 APs on this data (7.6× denser
  than our 26) with *large overlapping* measurement boxes — density and
  measurement support must be decoupled.

### Phase 1b (implemented): decoupled AP grid pitch

`cfg.alignment_points_step` / `-ap-step=N` (0 = legacy PSS geometry,
bit-identical). Explicit step keeps the correlation box at 2×half-box,
places APs every `step` px, and sizes the paste patch to
step + 0.75×half-box (legacy blend margin). Box may exceed patch (it only
feeds measurement + placement filters); per-AP frame ranking region is the
patch ∪ box union so a shrunken patch doesn't degrade ranking. Sidecar v11.
Validation: `-ap-step=27` (≈4× APs at unchanged box SNR) vs auto on the
Saturn data.

### Phase 1c (implemented): per-frame shift-field smoothing

Sidecar shift-field analysis (96 APs): the adjacent-AP disagreement structure
function is flat 20–200 px (RMS ≈ 1.1 px dy / 2.5 px dx) ⇒ per-(frame,AP)
measurement noise dominates true differential warp; dx ≫ dy is the aperture
problem on banded content. `cfg.alignment_points_smooth_radius`
(`-shift-smooth=F`, grid-step units, default 2.5, 0 = off): robust local
plane fit per frame per AP (tricube × Tukey bisquare, 3 IRLS rounds) over the
successful neighbours; failed pairs get the fit prediction. Applied after
each Stage B pass (the refined reference is built from the smoothed pass-1
field). Sidecar v12.

**Validation & the key negative result:** smoothing demonstrably denoises the
field (structure function 0.94/2.18 → 0.40/0.79 px at 20–35 px, now growing
with separation = physical warp), yet the stacked output is unchanged
(band-scale ratios 0.999). MTF arithmetic explains the whole plateau: after
Hann blending averages 2–4 APs, effective paste noise is ~0.3–0.5 px, and a
0.3 px Gaussian costs only ~3 % contrast at λ = 8 px, <1 % at λ = 16 px.
**Registration accuracy is not the binding constraint for band-scale detail
on this data** — phases 1/1b/1c all land within ±2 % of each other because
they were all pushing on a non-binding constraint. (Smoothing stays default
ON: measurably cleaner warp field, better failed-pair handling, no downside;
re-validate on lunar/solar where warp is real.)

### Phase 1d: zero-mean default flipped OFF (newsat.ser ablation)

The true problem capture (newsat.ser, 752×394×9064, 8-bit RGGB — much
harder: 12.5 % baseline failure rate vs 3.5 % on sat.ser) exposed zero-mean
as a net harm on low-SNR data: a 4-run ablation traced a 3–4 % band-power
loss, a doubling of the Stage B failure rate (12.5→23.7 %), and a
catastrophic 159-row border trim at legacy pitch directly to that one flag.
Mean-subtracting weak boxes leaves noise → honest rails → fallbacks that
cost more than plain NCC's DC-anchored near-zero answers. Default now OFF
(`-zero-mean` to opt in); refine + smoothing + dense pitch survive the
ablation clean. Best newsat config: `-ap-step=27 -stack-percent=50`
(199 APs, 11.7 % failures, band ratios ≈ 1.00 vs baseline, clean render).

**Where the AS! gap actually lives (next):** a heavily mid-band-boosted
render of the current stack shows banding approaching the AS! reference, so
part of the historical gap was the post-processing chain. Remaining
candidates for the residual fine-scale gap, in test order: (1) per-frame
RCD debayer of noisy 8-bit CFA vs AS!'s Bayer handling — try pre-debayered
input / green-extracted mono; (2) frame-quality ranking discrimination
(uint8-quantised Laplace σ) — float ranking, phase 2; (3) sub-pixel peak
bias (3×3 parabolic) — matters at fine scale only; upsampled-DFT fit.

### Phase 2a (mpp_improve_experimental): debayer A/B — DONE

`-debayer={rcd|bilinear|vng|ahd|amaze|dcb|hphd|igv|lmmse}` (stack-side;
`cfg.debayer_method`, sidecar v13) replaces the hard-coded BAYER_RCD in
`ser_read_frame`'s CFA path for the duration of Stage C. 5-way A/B on
newsat.ser (same sidecar, same shifts): **luminance band-scale detail is
independent of the demosaicing algorithm** (all variants within ±0.6 % of
RCD at every scale) — the debayer is ruled out as the band-detail gap.
Material finding: **LMMSE cuts stacked chroma noise 19 %** vs RCD (AMaZE/
VNG ~12 %; bilinear +4 % worse) at ~+14 s. **LMMSE is now the mpp stack
default** (`cfg.debayer_method`; `-debayer=rcd` restores the application-
wide SER choice). The `ser_set_debayer_method` override is scoped to
Stage C via an RAII guard, so nothing outside mpp sees it.

Remaining band-detail suspects: per-AP frame-quality ranking precision
(uint8-quantised Laplace σ → float; cheap test: stack-percent sweep
10/25/50 from the same sidecar — if percent barely matters, ranking is
weakly discriminating), 3×3 parabolic sub-pixel bias (fine-scale only),
and the AS! reference's third-hand processing chain (unquantifiable —
no unprocessed AS! stack available).

### Phase 2b: stack-percent sweep — the ranking DOES discriminate

Sweep 5/10/25/50/100 % from one sidecar (newsat, 199 APs, LMMSE): globe
band power rises **monotonically as the percent drops** — vs the 50 %
operating point, 10 % gives +13.8 %/+8.1 %/+5.8 % (3-8/6-16/12-32 px) and
5 % gives +21.9 %/+12.2 %/+10.7 %, while 100 % LOSES ~11 %. Sky-noise
correction confirms the gains are sharpness, not noise (even 5 % = 453
frames/AP is not noise-limited on this 9064-frame capture; lumHP rises
only ~19 % at 5 %). Two conclusions: (1) the uint8-quantised per-AP
ranking is genuinely informative — the float-ranking upgrade is a
refinement, not a fix; (2) **the 50 % operating point was costing more
band detail than any pipeline deficiency we found** — the AS!-convention
best-5–25 % is the right regime for long captures. User-side: prefer
10 % (or 5 % if the extra noise tolerates the denoise step).

### Phase 2c: float-precision ranking + lunar/solar regression suite

`cfg.rank_float_precision` (default true, `-no-float-rank`, sidecar v14):
Laplace-σ computed on CV_32F |lap|·α instead of PSS's uint8
convertScaleAbs image (quantised to steps of 256 16-bit ADU, saturating).
Same α scale so score magnitudes stay comparable; one shared quality
image feeds the global ranking and the per-AP selection. Unit test locks
the below-one-quantisation-step ordering case.

**Regression suite (branch defaults vs legacy flags):**
- Moon.ser (1500×1440×979 8-bit RGGB, surface mode, 557 APs): **SAFE** —
  band ratios 0.984/0.992/1.000, no AP-lattice imprint in the high-pass
  seam check, failure rate 0.6 % both ways, +43 % wall time.
- 11_29_46.ser (2180² ×1221 16-bit mono Hα full disc, 1218 APs): **SAFE,
  large win** — refined reference cuts Stage B failures 43.5 %→25.9 %
  and fine-detail band power rises **+21.7 %** (+8.7 % mid). Limb
  AP-pitch power spectra identical to 3 decimals at top and left arcs
  (no beading/scallop regression); prominences render identically; DC
  corrections drop 84→47 ADU. The granulation-scale texture that made
  per-AP correlation fail on the blurry mean reference is exactly what
  the refined reference fixes — solar is where phase 1 pays off.

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

### GUI — DONE (cog button + advanced settings window)

`mpp_advanced_dialog.ui` (hide-on-close GtkWindow, transient over the
control window) opened by a cog button (`emblem-system-symbolic`) on the
registration tab's multipoint page. Holds: refine reference (+ frames
spin, greyed when refine off), shift-smooth radius, AP grid pitch,
zero-mean opt-in, full-precision ranking, and the stacking demosaicing
drop-down (indices = interpolation_method 0..8), plus Reset-to-defaults.
Widget state is read by `fill_registration_structure_from_GUI` whether or
not the window is open; the demosaicing choice travels to the stack via
the sidecar cfg snapshot (the stack tab's overrides deliberately do not
touch it). Headless smoke test: builder loads with no missing-handler
warnings; UI files validate; dialog render verified. (No "PSS" in
user-facing strings.)
