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


def run(input_path, input_type, out_dir, max_frames=0):
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
    if max_frames > 0 and frames.number > max_frames:
        frames.number = max_frames
        frames.number_original = max_frames
        # Truncate every per-frame array PSS allocated so subsequent code
        # respects the new count.
        if hasattr(frames, 'frames_monochrome'):
            frames.frames_monochrome = frames.frames_monochrome[:max_frames]
        if hasattr(frames, 'frames_monochrome_blurred'):
            frames.frames_monochrome_blurred = frames.frames_monochrome_blurred[:max_frames]
        if hasattr(frames, 'frames_monochrome_blurred_laplacian'):
            frames.frames_monochrome_blurred_laplacian = frames.frames_monochrome_blurred_laplacian[:max_frames]
        if hasattr(frames, 'frames_average_brightness'):
            frames.frames_average_brightness = frames.frames_average_brightness[:max_frames]
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
    # Phase 5a oracle artifacts.
    stack_size = aps_obj.stack_size
    n_aps = len(aps_obj.alignment_points)
    n_frames = frames.number
    ap_qualities = np.zeros((n_aps, n_frames), dtype=np.float64)
    best_frame_indices = np.zeros((n_aps, stack_size), dtype=np.int64)
    for aidx, ap in enumerate(aps_obj.alignment_points):
        ap_qualities[aidx] = ap['frame_qualities']
        best_frame_indices[aidx] = ap['best_frame_indices']
    # Per-frame avg brightness (used by stacking for median-equalisation).
    frame_brightness = np.array(
        [frames.average_brightness(i) for i in range(n_frames)],
        dtype=np.float64)
    print(f"  stack_size={stack_size}")

    print("+++ Per-AP per-frame shifts (Phase 4 oracle)")
    # PSS computes per-AP per-frame shifts inside stack_frames; we drive
    # compute_shift_alignment_point directly here so we don't run the full
    # stacker.  Output is a (num_frames, num_aps, 2) tensor of integer
    # shifts (subpixel disabled to match the default config).
    aps_obj.set_reference_boxes_correlation()
    n_frames = frames.number
    n_aps = len(aps_obj.alignment_points)
    ap_shifts = np.zeros((n_frames, n_aps, 2), dtype=np.int64)
    ap_shift_success = np.zeros((n_frames, n_aps), dtype=np.uint8)
    for fidx in range(n_frames):
        mono_blurred = frames.frames_mono_blurred(fidx)
        for aidx in range(n_aps):
            shift, ok = aps_obj.compute_shift_alignment_point(
                mono_blurred, fidx, aidx,
                de_warp=True, subpixel_solve=False)
            ap_shifts[fidx, aidx, 0] = int(shift[0])
            ap_shifts[fidx, aidx, 1] = int(shift[1])
            ap_shift_success[fidx, aidx] = 1 if ok else 0

    print("+++ Stacking")
    stack = StackFrames(config, frames, rank, align, aps_obj, my_timer)
    stack.stack_frames()
    # stack_frames internally calls prepare_for_stack_blending and doesn't
    # mutate sum_single_frame_weights further, so capturing here is correct.
    sum_weights_buf = stack.sum_single_frame_weights.astype(np.float32).copy()
    number_holes = int(stack.number_stacking_holes)
    # Phase 5a oracle artifacts captured BEFORE merge (the merge step
    # consumes them).
    bg_patches = list(stack.background_patches) if stack.background_patches else []
    ap_stacking_bufs = [ap['stacking_buffer'].astype(np.float32).copy()
                        for ap in aps_obj.alignment_points]
    border_counts = (int(stack.border_y_low), int(stack.border_y_high),
                     int(stack.border_x_low), int(stack.border_x_high))
    avg_bg = (stack.averaged_background.astype(np.float32).copy()
              if stack.averaged_background is not None else None)
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
    if patch_yxyx is not None:
        np.savetxt(out_dir / "align_patch.csv",
                   np.array(patch_yxyx, dtype=np.int64).reshape(1, 4),
                   fmt="%d")
    # Flatten ap_shifts to a 2D CSV of shape (num_frames, 2 * num_aps) so the
    # whitespace-separated reader in the C++ tests can ingest it row by row.
    np.savetxt(out_dir / "ap_shifts.csv",
               ap_shifts.reshape(n_frames, 2 * n_aps), fmt="%d")
    np.savetxt(out_dir / "ap_shift_success.csv",
               ap_shift_success, fmt="%d")
    # Also persist the intersection bounds so the C++ test can build matching
    # per-frame offsets without re-deriving them.
    intersect = np.array([align.intersection_shape[0][0],
                          align.intersection_shape[0][1],
                          align.intersection_shape[1][0],
                          align.intersection_shape[1][1]], dtype=np.int64)
    np.savetxt(out_dir / "intersection.csv", intersect.reshape(1, 4), fmt="%d")
    # Dump PSS's mean_frame (int32) as a raw binary so C++ tests can diff
    # against it without depending on a FITS reader. Layout: little-endian
    # int32, rows * cols values. Companion dims file gives the shape.
    mf = align.mean_frame.astype(np.int32)
    mf.tofile(out_dir / "mean_frame_raw.bin")
    np.savetxt(out_dir / "mean_frame_dims.csv",
               np.array([mf.shape[0], mf.shape[1]], dtype=np.int64).reshape(1, 2),
               fmt="%d")
    # Also the post-AP blurred mean_frame (same convention) — this is what
    # set_reference_boxes_correlation samples from.
    mfb = aps_obj.mean_frame.astype(np.int32)
    mfb.tofile(out_dir / "mean_frame_blurred_raw.bin")
    # Phase 5a oracles.
    np.savetxt(out_dir / "ap_qualities.csv", ap_qualities, fmt="%.17g")
    np.savetxt(out_dir / "best_frame_indices.csv", best_frame_indices, fmt="%d")
    np.savetxt(out_dir / "frame_brightness.csv", frame_brightness, fmt="%.17g")
    with open(out_dir / "stack_size.txt", "w") as f:
        f.write(f"{stack_size}\n")
    sum_weights_buf.tofile(out_dir / "sum_single_frame_weights.bin")
    np.savetxt(out_dir / "sum_single_frame_weights_dims.csv",
               np.array([sum_weights_buf.shape[0], sum_weights_buf.shape[1]],
                        dtype=np.int64).reshape(1, 2), fmt="%d")
    with open(out_dir / "number_stacking_holes.txt", "w") as f:
        f.write(f"{number_holes}\n")
    # Background patches: one row per patch, (y_lo y_hi x_lo x_hi).
    if bg_patches:
        np.savetxt(out_dir / "background_patches.csv",
                   np.array([(p['patch_y_low'], p['patch_y_high'],
                              p['patch_x_low'], p['patch_x_high']) for p in bg_patches],
                            dtype=np.int64),
                   fmt="%d")
    else:
        # Touch an empty file to make the test detect "no patches" deterministically.
        (out_dir / "background_patches.csv").write_text("")
    # Border counts: single row of 4 ints (y_lo y_hi x_lo x_hi).
    np.savetxt(out_dir / "stack_borders.csv",
               np.array([border_counts], dtype=np.int64), fmt="%d")
    # Per-AP stacking buffers: concatenated binary + a meta CSV with one row
    # per AP (rows cols byte_offset).
    meta = []
    off = 0
    with open(out_dir / "ap_stacking_buffers.bin", "wb") as f:
        for buf in ap_stacking_bufs:
            assert buf.dtype == np.float32
            f.write(buf.tobytes(order="C"))
            meta.append((buf.shape[0], buf.shape[1], off))
            off += buf.nbytes
    np.savetxt(out_dir / "ap_stacking_buffers_meta.csv",
               np.array(meta, dtype=np.int64), fmt="%d")
    if avg_bg is not None:
        avg_bg.tofile(out_dir / "averaged_background.bin")
        np.savetxt(out_dir / "averaged_background_dims.csv",
                   np.array([avg_bg.shape[0], avg_bg.shape[1]],
                            dtype=np.int64).reshape(1, 2), fmt="%d")
    # Final uint16 stacked image (post-merge, post-trim, post-clip).
    stacked_u16 = np.asarray(stacked, dtype=np.uint16)
    stacked_u16.tofile(out_dir / "stacked_u16.bin")
    np.savetxt(out_dir / "stacked_u16_dims.csv",
               np.array([stacked_u16.shape[0], stacked_u16.shape[1]],
                        dtype=np.int64).reshape(1, 2), fmt="%d")
    # used_alignment_points is jagged; flatten with row lengths.
    used_aps = frames.used_alignment_points
    with open(out_dir / "used_alignment_points.csv", "w") as f:
        for fidx in range(n_frames):
            lst = used_aps[fidx] if fidx < len(used_aps) else []
            f.write(" ".join(str(a) for a in lst) + "\n")

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
    ap.add_argument("--n", type=int, default=0, help="Truncate to first N frames (0 = all).")
    args = ap.parse_args()
    try:
        run(args.input, args.type, args.out_dir, max_frames=args.n)
    except Exception:
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
