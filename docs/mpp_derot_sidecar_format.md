# `.derot` derotation sidecar — format spec

## Purpose

A `.derot` sidecar holds a per-sequence **planetary derotation plan**: for every
frame, the disk geometry (sub-observer latitude, central-meridian longitude, and
pole position angle) at that frame's time, plus the reference epoch and the disk
fit. It lets the multipoint pipeline derotate a long capture — or several
captures taken minutes apart — to a common instant, folding the rotation into
the stack's single resample so planetary rotation no longer smears the result.

It is written next to the sequence as `<seqname>.derot` by:

- the `derotate` command (`process_derotate`, `src/core/command.c`)
- the derotation tool window (`src/gui-gtk4/derotation.c`)

both via the shared builder `mpp_derot_build` (`src/registration/mpp/mpp_derot_build.c`).

It is read by:

- `register_mpp` / `mpp_analyze` / `mpp_compute_shifts` — derotation-aware
  registration (measures per-AP shifts in the epoch frame)
- `stack_mpp_handler` / `process_stack_mpp` — Stage C, which forces the
  warp-field stacking engine when a `.derot` is present

It is deliberately **independent of the `.mpp` sidecar** so that re-running
registration (which rewrites `.mpp`) never clobbers the derotation plan.

## Endianness, alignment, types

- **Byte order**: little-endian (Siril builds for LE hosts), written field by
  field with no compiler padding between blocks.
- **Scalar types**: `int32_t`, `uint32_t`, IEEE-754 `double` (8 bytes).
- **No checksum / no trailer.** Truncation is detected only by short reads;
  `mpp_derot_read` returns `MPP_EIO` on any short read or magic/version mismatch.

## Top-level layout

```
+----------------------------------------------------------+
| 8  bytes : magic       "SIRILDRT"                        |
| 4  bytes : uint32 version    (currently 2)               |
| 4  bytes : uint32 flags      (reserved; 0)               |
+----------------------------------------------------------+
| 6 × int32  : body, rot_system, ephem_version,            |
|              num_frames (N), frame_rows, frame_cols      |
+----------------------------------------------------------+
| 14 × double : header scalars (see below)                 |
+----------------------------------------------------------+
| N × double : jd[N]            (per-frame UTC Julian date) |
| N × double : sub_obs_lat[N]   (B, degrees)               |
| N × double : cm[N]            (central meridian, degrees) |
| N × double : pole_pa[N]       (P, degrees)               |
+----------------------------------------------------------+
```

### Header int32 block

| Index | Field          | Notes                                              |
|------:|----------------|----------------------------------------------------|
| 0     | `body`         | `planet_body_t` (0 = Jupiter, 1 = Saturn, 2 = Mars)|
| 1     | `rot_system`   | `rot_system_t` (0/1/2 = System I/II/III)           |
| 2     | `ephem_version`| ephemeris-model version, provenance only           |
| 3     | `num_frames` N | must equal the sequence frame count to apply       |
| 4     | `frame_rows`   | sequence frame height the disk fit is in            |
| 5     | `frame_cols`   | sequence frame width the disk fit is in             |

### Header double block (14 values)

```
epoch_jd            reference epoch, UTC Julian date
obs_lat, obs_lon    observer geodetic latitude / longitude (deg); NaN = geocentre
obs_elev            observer elevation (m); NaN = geocentre
cx, cy              disk centre, ORIGINAL full-frame pixels
r_eq                apparent equatorial radius, pixels
flattening          geometric flattening f = (Req-Rpol)/Req (from the ephemeris)
parity              +1 normal, -1 mirrored image
pole_angle_epoch    in-image position angle of the projected north pole, RADIANS
epoch_sub_obs_lat   B at the epoch, degrees
epoch_cm            central-meridian longitude at the epoch, degrees
epoch_pole_pa       P at the epoch, degrees
(reserved)          0.0
```

## Per-frame arrays

Four `double[N]` arrays, written in order: `jd`, `sub_obs_lat`, `cm`,
`pole_pa`. They are **resolved at compute time** from the built-in ephemeris and
stored, so stacking is reproducible regardless of later ephemeris-model changes
or timestamp re-reads. Angles are degrees; `cm` is in the chosen `rot_system`.

## Conventions

- **Disk fit** (`cx`, `cy`, `r_eq`) is in the sequence's original full-frame
  pixel coordinates. Stage C scales it into the drizzled intersection canvas;
  Stage A/B use it at full resolution (`mpp_derot_frame_map`).
- **`cm` / `pole_angle`**: derotation uses *differences* of these between a
  frame and the epoch, so any constant offset (delta-T, J2000-vs-of-date) cancels
  — see `planet_ephem.h`.

## Relationship to the `.mpp` sidecar — the consistency guard

Applying derotation on top of per-AP shifts that were measured *without*
derotation (against a rotationally-smeared mean) would double-count the
rotation. To prevent that:

- `mpp_derot_fingerprint()` is a 64-bit FNV-1a hash over the whole plan (body,
  system, epoch, observer, disk fit, and the per-frame arrays).
- When derotation-aware registration runs, that fingerprint is stamped into the
  `.mpp` (`mpp_run_t.derot_fingerprint`, sidecar v12).
- Stage C **refuses to stack** unless the `.mpp`'s stored fingerprint equals the
  fingerprint of the `.derot` being applied — i.e. the sequence must be
  re-registered after any change to (or addition of) the derotation plan.

## Compatibility & evolution

- **Magic + version are checked strictly.** Mismatch → `MPP_EIO`.
- The current writer always emits version 2.
- **v1 → v2** added `frame_rows`/`frame_cols` so a plan from a different-
  resolution capture (with the same frame count) is rejected rather than
  mis-placing the disk. v1 files are not readable.

## File lifecycle

- **Written**: `mpp_derot_write` (`src/registration/mpp/mpp_derot_sidecar.c`) —
  overwrite-on-write; no atomic rename.
- **Read**: `mpp_derot_read` returns a freshly-allocated `mpp_derot_t*` owned by
  the caller; free via `mpp_derot_free`.
- **Location**: `<seqname>.derot`, alongside the sequence file.

## Size estimate

Dominated by the four per-frame `double[N]` arrays: `32 × N` bytes (plus a
~150-byte header). For N = 40 000 frames that is ~1.3 MB — small relative to the
sequence.
