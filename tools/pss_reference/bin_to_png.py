"""
Convert raw uint16 .bin outputs (paired with a `*_dims.csv`) to viewable
16-bit PNGs. Run it on one of the oracle_out directories after the mpp
stack oracle test has been driven with MPP_DUMP_RESULT_DIR set so both
PSS's and Siril's u16 outputs are present on disk.

Also writes a `stacked_u16_diff_x256.png` showing |PSS − Siril| amplified
by 256× for visibility (raw diffs are 0..few levels and would be invisible).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image


def load(bin_path: Path, dims_csv: Path) -> np.ndarray:
    H, W = (int(v) for v in dims_csv.read_text().split())
    return np.fromfile(bin_path, dtype=np.uint16).reshape(H, W)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("oracle_dir", help="e.g. tools/pss_reference/oracle_out_real")
    args = ap.parse_args()
    d = Path(args.oracle_dir)

    arrays = {}
    for stem in ("stacked_u16", "stacked_u16_siril"):
        bin_p = d / f"{stem}.bin"
        dims_p = d / f"{stem}_dims.csv"
        if not (bin_p.exists() and dims_p.exists()):
            print(f"  skip (missing): {stem}")
            continue
        img = load(bin_p, dims_p)
        png = d / f"{stem}.png"
        Image.fromarray(img, mode="I;16").save(png)
        arrays[stem] = img
        print(f"  wrote {png}  shape={img.shape}, min={img.min()}, max={img.max()}")

    if "stacked_u16" in arrays and "stacked_u16_siril" in arrays:
        pss, siril = arrays["stacked_u16"], arrays["stacked_u16_siril"]
        if pss.shape != siril.shape:
            print(f"  shape mismatch: PSS {pss.shape} vs Siril {siril.shape}")
            return 1
        absd = np.abs(pss.astype(np.int32) - siril.astype(np.int32))
        absd_vis = np.clip(absd * 256, 0, 65535).astype(np.uint16)
        Image.fromarray(absd_vis, mode="I;16").save(d / "stacked_u16_diff_x256.png")
        n = absd.size
        exact = int((absd == 0).sum())
        print(f"  wrote {d / 'stacked_u16_diff_x256.png'}  (|Δ| × 256)")
        print(f"  exact={100 * exact / n:.4f}%  worst Δ={absd.max()}  "
              f"mean |Δ|={absd.mean():.6f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
