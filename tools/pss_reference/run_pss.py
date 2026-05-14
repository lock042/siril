"""
Headless PSS pipeline driver.

Based on PSS's own Test_programs/main_program.py::workflow but stripped of GUI
dependencies and instrumented to dump intermediate artifacts to disk for use
as a unit-test oracle.

Outputs (in --out-dir):
  oracle.npz   per-frame quality, global shifts, AP coordinates, per-AP shifts
  ref_avg.fits average reference frame from top N% of frames
  stacked.fits final stacked image
"""

import argparse
import os
import sys
import traceback
import warnings
from pathlib import Path

import numpy as np

# Suppress harmless Windows-path SyntaxWarnings from PSS frames*.py
warnings.filterwarnings("ignore", category=SyntaxWarning)

PSS_ROOT = Path(__file__).resolve().parent.parent.parent / "PlanetarySystemStacker" / "planetary_system_stacker"
sys.path.insert(0, str(PSS_ROOT))

# Headless matplotlib backend before any PSS import.
import matplotlib
matplotlib.use("Agg")

# Avoid PyQt5 import side-effects where possible; create a QApplication only if PSS demands it.
try:
    from PyQt5 import QtWidgets
    _qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
except Exception:
    _qapp = None

from configuration import Configuration  # noqa: E402
from frames import Frames  # noqa: E402
from rank_frames import RankFrames  # noqa: E402
from align_frames import AlignFrames  # noqa: E402
from alignment_points import AlignmentPoints  # noqa: E402
from stack_frames import StackFrames  # noqa: E402
from timer import timer  # noqa: E402


def run(input_path, input_type, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    my_timer = timer()
    my_timer.create("Execution over all")

    if input_type == "video":
        names = input_path
    else:
        names = sorted(
            os.path.join(input_path, name)
            for name in os.listdir(input_path)
            if not name.endswith(".npz") and not name.startswith(".")
        )

    print(f"+++ Driving PSS oracle on {input_path} ({input_type})")
    config = Configuration()
    config.initialize_configuration()

    print("+++ Reading frames")
    frames = Frames(config, names, type=input_type)
    print(f"  number={frames.number} shape={frames.shape} color={frames.color}")

    print("+++ Ranking frames")
    rank = RankFrames(frames, config)
    rank.frame_score()
    quality = np.asarray(rank.frame_ranks, dtype=np.float64)
    print(f"  best={rank.frame_ranks_max_index} q[best]={quality.max():.6g}")

    print("+++ Global alignment")
    align = AlignFrames(frames, rank, config)
    patch_yxyx = None
    if config.align_frames_mode == "Surface":
        y_low, y_high, x_low, x_high = align.compute_alignment_rect(
            config.align_frames_rectangle_scale_factor
        )
        patch_yxyx = (y_low, y_high, x_low, x_high)
        print(f"  patch y=({y_low},{y_high}) x=({x_low},{x_high})")
    align.align_frames()
    shifts = np.asarray(align.frame_shifts, dtype=np.int64)  # shape (n, 2) — (dy, dx)

    print("+++ Average reference frame")
    average = align.average_frame()

    print("+++ Alignment-point grid")
    aps_obj = AlignmentPoints(config, frames, rank, align)
    aps_obj.create_ap_grid()
    ap_list = aps_obj.alignment_points
    print(f"  aps_kept={len(ap_list)} "
          f"dim_dropped={aps_obj.alignment_points_dropped_dim} "
          f"struct_dropped={aps_obj.alignment_points_dropped_structure}")
    ap_yx = np.array([(ap["y"], ap["x"]) for ap in ap_list], dtype=np.int64)
    ap_box = np.array(
        [(ap["box_y_low"], ap["box_y_high"], ap["box_x_low"], ap["box_x_high"]) for ap in ap_list],
        dtype=np.int64,
    )

    print("+++ Per-AP frame qualities")
    aps_obj.compute_frame_qualities()

    print("+++ Stacking")
    stack = StackFrames(config, frames, rank, align, aps_obj, my_timer)
    stack.stack_frames()
    stacked = stack.merge_alignment_point_buffers()
    if config.drizzle_factor_is_1_5:
        stack.half_stacked_image_buffer_resolution()

    print("+++ Saving oracle artifacts")
    np.savez(
        out_dir / "oracle.npz",
        quality=quality,
        global_shifts=shifts,
        ap_yx=ap_yx,
        ap_box=ap_box,
        align_patch=np.array(patch_yxyx if patch_yxyx is not None else (-1, -1, -1, -1), dtype=np.int64),
    )
    # Plain-text mirrors of the artifacts that C++ unit tests consume directly,
    # one value per line (or whitespace-separated for multi-column).
    np.savetxt(out_dir / "quality.csv", quality, fmt="%.17g")
    np.savetxt(out_dir / "global_shifts.csv", shifts, fmt="%d")
    np.savetxt(out_dir / "ap_yx.csv", ap_yx, fmt="%d")
    np.savetxt(out_dir / "ap_box.csv", ap_box, fmt="%d")

    # Save reference + stacked as FITS via PSS's own writer (keeps any bit-depth conventions).
    Frames.save_image(str(out_dir / "ref_avg.fits"), average, color=frames.color,
                      header=config.global_parameters_version)
    Frames.save_image(str(out_dir / "stacked.fits"), stacked, color=frames.color,
                      header=config.global_parameters_version)

    my_timer.stop("Execution over all")
    my_timer.print()
    print(f"+++ Done. Artifacts in {out_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Path to video file or image directory")
    ap.add_argument("--type", choices=["video", "image"], default="image")
    ap.add_argument("--out-dir", default="oracle_out")
    args = ap.parse_args()
    try:
        run(args.input, args.type, args.out_dir)
    except Exception:
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
