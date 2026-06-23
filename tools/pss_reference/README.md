# PSS reference oracle

Drives upstream PlanetarySystemStacker headlessly and dumps intermediate
artifacts to disk, for use as ground truth by Siril's `pss` unit tests.

## Environment

PSS targets Python 3.5/3.6 with `numpy < 1.23` per its `setup.py`. As of bring-up
on 2026-05-14 it nonetheless **runs unmodified** on Python 3.12 / NumPy 2.4 /
OpenCV 4.13 / astropy 7.2 / scipy 1.17 / scikit-image 0.26 / matplotlib 3.10 /
PyQt5. No upstream PSS source patches have been required so far. If runtime
errors appear once richer inputs are exercised, log them in `PSS_PATCHES.md` (to
be created at that point) and patch the upstream clone in place — it is
untracked in this repo and lives at `../../PlanetarySystemStacker/`.

## Setup

    cd tools/pss_reference
    python3 -m venv .venv
    .venv/bin/pip install -U pip wheel
    .venv/bin/pip install numpy opencv-python astropy scipy scikit-image \
        matplotlib psutil PyQt5 Pillow

## Generating synthetic frames

    .venv/bin/python gen_synth_frames.py

Writes 24 16-bit-PNG frames + a `truth.npz` (per-frame ground-truth jitter and
blur sigma) to `test_data/synth_planet/`.

## Running the oracle

    QT_QPA_PLATFORM=offscreen \
        .venv/bin/python run_pss.py test_data/synth_planet --type image \
            --out-dir oracle_out

`QT_QPA_PLATFORM=offscreen` is required on headless boxes because PSS imports
PyQt5 unconditionally and instantiates a QApplication at module load time.

For video input use `--type video` and pass the path to an `.avi`/`.mov`/`.mp4`
/`.ser` file.

Outputs:

| File | Contents |
|---|---|
| `oracle.npz` | `quality` (per-frame), `global_shifts` (n×2, dy/dx), `ap_yx` (m×2), `ap_box` (m×4, y_lo/y_hi/x_lo/x_hi), `align_patch` (y_lo,y_hi,x_lo,x_hi) |
| `ref_avg.fits` | Average reference frame computed from the top N% best frames |
| `stacked.fits` | Final stacked image |

These are the artifacts the Siril unit tests will compare against.

## Sign conventions (verified at bring-up)

PSS reports `global_shifts[i] = [dy, dx]` such that the *frame must be shifted
by `(dy, dx)` to align with the reference*. This is the negative of the
generator's "frame i was shifted by `(dy, dx)` away from truth" convention used
in `truth.npz`. The Siril port must document and conform to PSS's convention.
