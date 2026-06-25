# Multi-sequence planetary derotation (Option B — zero-copy)

Branch: `pss-derot-multi` (off `pss-derot`).

## Goal
Stack several sequences of the same planet **together, as one accumulation**,
without an intermediate merged file, by:

1. derotating every sequence to a **common reference epoch** (each keeps its own
   disk fit), then
2. MPP-registering every sequence's frames to **one user-designated reference
   sequence's reference frame**, then
3. MPP-stacking all of them into that common canvas.

For **mono** inputs, sequences may be tagged **R / G / B** (or luminance): each
channel runs its own multi-source pipeline against the *same* shared reference,
so the channel stacks are co-registered and compose directly into RGB.

Single interpolation throughout (derotation + global + AP residual fold into one
`cv::remap` per frame, exactly as the single-sequence warp engine does today).

## Why this over "merge"
`merge` already gives a working same-filter workflow, but copies data and forces
one in-image pole angle for the whole set (loses per-segment field rotation on
alt-az, can't do per-sequence fits). Option B is zero-copy, keeps a fit per
sequence, and extends cleanly to RGB.

## Key design decisions (settled with maintainer)
- **Scale invariant** — same camera assumed; no `r_eq` normalisation. The warp
  scale stays the drizzle factor.
- **Rotation from each sequence's own fit** (measurable for Jupiter/Saturn; for
  Mars/ice giants all channels share the ephemeris pole-PA at the epoch).
- **One global reference sequence**, user-selected, shared across all channels
  (so R/G/B co-register). The global alignment is disc-CoG (filter-invariant);
  AP residuals may be slightly weaker cross-filter but are small.
- Common epoch = midpoint of the union of all sequences' spans
  (`mpp_derot_union_epoch`, already on `pss-derot`).

## Components
- **C1 — session model** (`registration/mpp/mpp_session.{c,h}`). Pure data +
  logic: the list of sequences, each with a channel tag and capture span; the
  designated reference; body/system; shared-epoch computation; channel grouping.
  Unit-tested. *(this branch's first commit)*
- **C2 — per-sequence `.derot` at the shared epoch.** Reuses `derotate` /
  `mpp_derot_build` with the session epoch and each sequence's own fit. Mostly
  exists (`-epoch-from`); wire it to the session.
- **C3 — derot-map generalisation.** `mpp_derot_frame_map` currently uses one
  disk fit for both the epoch (output) and frame (source) sides. Split them: the
  **output/epoch** side = the reference sequence's geometry+canvas (common), the
  **frame/source** side = each sequence's own fit. This is what lands every
  sequence in the one reference canvas.
- **C4 — multi-source frame provider.** A global frame index → (sequence, local
  index) mapping over the union of the channel's sequences; reads the right SER
  and applies the right per-frame derotation plan.
- **C5 — multi-source analysis.** Reference mean + AP grid from the designated
  reference sequence (common canvas); option to build the mean from all frames
  later for S/N.
- **C6 — multi-source shift measurement.** Per-AP shifts for the union of frames
  against the common reference.
- **C7 — multi-source stacking.** Warp every frame (derot + global + AP) once
  into the common canvas; accumulate per channel.
- **C8 — orchestration.** Per-channel pipelines sharing the one reference; CLI
  entry first (headless-testable), then GUI.
- **C9 — GUI.** Add sequences (multi-select File Open), tag R/G/B/mono/lum,
  designate the reference, fit each disk, run. Produces per-channel stacks.
- **C10 — RGB compose.** The per-channel stacks are co-registered; compose via
  the existing Compositing channel-assembly (or a thin built-in step).

## Testing
Synthetic multi-"sequence" recovery: render one rotating disk into N short
"sequences" at different epochs/positions, derotate to the shared epoch,
multi-source register+stack → assert the single combined stack recovers the disk
and beats each sequence stacked alone; and that two channels co-register to
sub-pixel.

## Status
- C1: in progress.
