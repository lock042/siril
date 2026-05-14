"""
Extract N frames from a SER file as 8-bit mono PNGs.

Reads raw frames via PSS's ser_parser, debayers Bayer SERs to grayscale via
OpenCV, then writes 8-bit mono PNGs to an output directory. Designed for
generating realistic real-data oracle inputs for the mpp port — the
resulting PNGs feed both PSS (--type image) and the Siril unit tests with
identical content.

For non-Bayer / non-color SERs the frames are saved as-is.
"""

import argparse
import sys
from pathlib import Path

import cv2
import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent
                       / "PlanetarySystemStacker" / "planetary_system_stacker"))
from ser_parser import SERParser  # noqa: E402


BAYER_CODE_TO_CV = {
    # SER ColorID → OpenCV BayerXX2GRAY enum. PSS ColorIDs 8-11 are Bayer.
    8:  cv2.COLOR_BayerRG2GRAY,   # RGGB
    9:  cv2.COLOR_BayerGR2GRAY,   # GRBG
    10: cv2.COLOR_BayerGB2GRAY,   # GBRG
    11: cv2.COLOR_BayerBG2GRAY,   # BGGR
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ser_path")
    ap.add_argument("out_dir")
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--start", type=int, default=0)
    args = ap.parse_args()

    p = SERParser(args.ser_path)
    n_total = p.header["FrameCount"]
    color_id = p.header["ColorID"]
    depth = p.header["PixelDepthPerPlane"]
    print(f"SER: {n_total} frames, {p.header['ImageWidth']}×{p.header['ImageHeight']}, "
          f"ColorID={color_id}, depth={depth}")
    n = min(args.n, n_total - args.start)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    cv_bayer = BAYER_CODE_TO_CV.get(color_id)
    for i in range(n):
        raw = p.read_frame(args.start + i)
        if cv_bayer is not None and raw.ndim == 2:
            mono = cv2.cvtColor(raw, cv_bayer)
        elif raw.ndim == 3:
            mono = cv2.cvtColor(raw, cv2.COLOR_RGB2GRAY)
        else:
            mono = raw
        if mono.dtype != np.uint8:
            # 16-bit input: clip to uint8 by rounding (we keep test inputs
            # as 8-bit mono so the PSS pipeline uses the 8-bit path).
            mono = np.clip(mono, 0, 255).astype(np.uint8) if depth <= 8 \
                else (mono >> 8).astype(np.uint8)
        Image.fromarray(mono, mode="L").save(out / f"frame_{i:04d}.png")
        if (i + 1) % 50 == 0:
            print(f"  wrote {i + 1}/{n}")
    print(f"Wrote {n} frames to {out}")


if __name__ == "__main__":
    main()
