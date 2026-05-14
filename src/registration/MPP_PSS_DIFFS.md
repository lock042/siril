# Deliberate departures from PlanetarySystemStacker

This document records every place where Siril's `mpp` pipeline diverges from
upstream [PlanetarySystemStacker](https://github.com/Rolf-Hempel/PlanetarySystemStacker)
(PSS) by intent. Inadvertent divergences are bugs; everything below is a
considered decision with the reasoning attached.

The yardstick: on a pre-debayered mono fixture (where neither side's choice
of debayer algorithm is in play) the Siril stacked output and the PSS oracle
match at **128 dB PSNR** on the bundled real-SER 500-frame test fixture,
worst Δ = 1 out of 65535 (Phase 7.2). The mpp algorithmic chain is therefore
bit-near-equivalent to PSS for every PSS-supplied input regime; the
divergences listed here are everything *outside* that chain.

## 1. Bayer debayer algorithm

| Aspect | PSS | Siril mpp |
|---|---|---|
| Library | OpenCV `cv2.cvtColor` with `COLOR_BayerRG2BGR` etc. | Siril's `ser.c` reader |
| Algorithm | Bilinear (cv2 default) | RCD (hard-coded; see `ser.c:897`) |
| Output dimensions on `test-big.ser` | 264 × 258 | 263 × 258 (1 px narrower) |

**Why:** Siril's SER reader pipeline is shared across the whole application
and is the authoritative path for every other tool in Siril. Re-implementing
a Bayer→RGB pass inside the mpp module would split that authority and create
two ways for a Bayer SER to be interpreted within one Siril run. Keeping the
common reader is the right call for consistency and maintenance even though
it costs ~80 dB of headline PSNR vs the PSS oracle in the colour case.

**Observable impact:** On the 500-frame colour comparison, Siril produces a
250 × 247 × 3 stacked image and PSS produces 251 × 247 × 3. After a best-shift
overlap search (`(dy, dx) = (1, -1)`) and crop to the 246 × 250 common region:

| Drizzle | Overall PSNR | R | G | B |
|---|---|---|---|---|
| Off (1×) | 49.63 dB | 51.15 dB | 48.32 dB | 49.87 dB |
| 2× | 54.05 dB | 55.18 dB | 53.45 dB | 53.70 dB |

The green channel is worst, consistent with RGGB carrying two green samples
per Bayer cell that interpolate differently between bilinear and RCD. The
drizzle=2 numbers are higher because the 2× resampling smooths per-pixel
debayer-algorithm noise.

**If you need PSS-exact colour output**, pre-debayer the SER outside Siril
with the same `cv2.cvtColor` call PSS uses and feed Siril a 3-layer FITS
sequence — Siril will then leave the frames alone (no Bayer interpretation
in `mpp_seq_read_frame` for non-CFA input).

## 2. Float-rounding ordering

PSS chains NumPy and OpenCV operations in a specific order. Each operation
has well-defined IEEE 754 rounding; the order in which results compose
determines the last-bit pattern of intermediate results. Three rounding
behaviours sit at boundaries between the two libraries and are matched in
Siril mpp **by design** rather than by accident:

| Place | Behaviour | Source |
|---|---|---|
| Average reference frame cast back to integer | NumPy `astype(int32)` truncation toward zero | `align_frames.compute_average_frame` |
| Per-AP correlation peak coords (no subpixel) | OpenCV `cvRound` half-to-even | `cv::minMaxLoc` |
| Final stacked output uint16 cast | skimage `img_as_uint` round-half-away-from-zero | `stack_frames._u16_stack` |

`mpp_align.cpp::align_average_frame` and
`mpp_stack.cpp::stack_merge_alignment_point_buffers` explicitly implement
the truncate-toward-zero and round-half-away-from-zero behaviours
respectively, even though they're not the libraries' defaults at those call
sites — both are bit-equivalent to PSS on every fixture in the test suite.

## 3. Frame analysis layer for colour input

| Aspect | PSS | Siril mpp |
|---|---|---|
| Analysis-layer source | `Frames.frames_mono_channel = 'panchromatic'` (`color_index = 3`) | Same |
| Conversion | `cv2.cvtColor(rgb, COLOR_RGB2GRAY)` (Rec. 601 luma) | `cv::cvtColor(rgb, COLOR_RGB2GRAY)` (same coefficients) |

**No departure** — explicitly listed here because earlier Siril prototypes
incorrectly used the green channel alone (`color_index = 1`), which produces
subtly different AP placement on bright-saturated targets like Jupiter. The
current `mpp.cpp::read_analysis_frame` matches PSS exactly.

## 4. Penalty matrix for off-centre correlation peaks

`alignment_points_penalty_factor = 0.00025` is the small radial penalty PSS
applies to discourage shift solutions that drift far from `(0, 0)` in the
correlation surface. Siril mpp applies it in **Stage B only** (per-AP per-frame
shift compute used during stacking) and **not in Stage A** (where Phase 4's
oracle was captured). This matches PSS exactly:
`alignment_points.compute_frame_qualities` does not pass the weight matrix
into the correlation, but `stack_frames.stack_frames` does.

## 5. Drizzle interpolation

| Aspect | PSS | Siril mpp |
|---|---|---|
| Resize for drizzle | `cv2.resize(frame, scale=K, INTER_LINEAR)` | Same (`cv::resize` with `cv::INTER_LINEAR`) |

**No departure.** The plan originally specified `INTER_CUBIC` for bicubic
drizzle, which would have produced sharper results. PSS source actually
uses `INTER_LINEAR` at `stack_frames.py:350`; the plan was wrong, and the
implementation matches PSS rather than the plan to keep the bit-equivalence
gauge meaningful.

## 6. Drizzle factor encoding

| Aspect | PSS | Siril mpp |
|---|---|---|
| Surface | `stack_frames_drizzle_factor_string` ∈ {`Off`, `1.5x`, `2x`, `3x`} | `cfg.drizzle_factor` ∈ {1, 2, 3} |
| `1.5x` representation | `drizzle_factor = 3` + flag `drizzle_factor_is_1_5 = True` | `cfg.drizzle_factor = 3` (no separate flag yet) |

**Why:** Siril doesn't yet implement the post-stack 0.5× downsample that PSS
applies when `drizzle_factor_is_1_5` is set, so the user-facing
`-drizzle=1.5` flag currently behaves as `-drizzle=3`. Track at Phase 5b.

## 7. STScI drizzle, Bayer-drizzle (Phase 5b)

Deferred per the plan's scope decisions. PSS does not offer STScI-style
drizzle or true raw-Bayer drizzle, so neither feature can be PSS-compared
on equal footing in Phase 7. They will gain their own synthetic-truth tests
when Phase 5b lands.

## 8. Where to look in code

| Topic | File |
|---|---|
| Reference frame integer cast | `src/registration/mpp_align.cpp::align_average_frame` |
| Final uint16 cast | `src/registration/mpp_stack.cpp::stack_merge_alignment_point_buffers` |
| Analysis layer for RGB input | `src/registration/mpp.cpp::read_analysis_frame` |
| Drizzle resize | `src/registration/mpp_stack.cpp::stack_frames_loop` |
| Debayer choice (out of mpp's hands) | `src/io/ser.c::ser_read_frame` (`BAYER_RCD` literal) |
