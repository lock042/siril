"""
Generate a tiny synthetic image-directory dataset for PSS bring-up tests.

Produces N frames of a simulated planetary disk with prescribed per-frame
sub-pixel jitter and quality variation (Gaussian blur sigma drifting with
frame index). Output is a directory of 16-bit PNGs ready to feed into PSS
with `input_type='image'`. This is *not* yet a SER; SER generation is a
later Phase-0 task.
"""

import argparse
import os
from pathlib import Path

import numpy as np
from PIL import Image
import cv2


def make_truth_frame(h, w):
    """Static high-contrast scene: bright disk + small bright spots + textured
    background. Returns float32 in [0, 1]."""
    yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)
    cy, cx = h / 2.0, w / 2.0
    r = np.hypot(yy - cy, xx - cx)
    disk_radius = min(h, w) * 0.30
    edge = 2.0
    # Soft-edged bright disk.
    disk = 0.75 * 0.5 * (1.0 - np.tanh((r - disk_radius) / edge))
    # Sparse bright spots on the disk surface.
    rng = np.random.default_rng(seed=42)
    spots = np.zeros_like(disk)
    for _ in range(30):
        sy = cy + rng.uniform(-disk_radius * 0.7, disk_radius * 0.7)
        sx = cx + rng.uniform(-disk_radius * 0.7, disk_radius * 0.7)
        sr = rng.uniform(1.5, 3.5)
        amp = rng.uniform(0.10, 0.20)
        spots += amp * np.exp(-((yy - sy) ** 2 + (xx - sx) ** 2) / (2 * sr * sr))
    # Faint texture across the disk.
    tex = 0.04 * rng.standard_normal((h, w)).astype(np.float32)
    tex = cv2.GaussianBlur(tex, (5, 5), 0)
    mask = (r < disk_radius + edge).astype(np.float32)
    img = disk + (spots + tex) * mask
    return np.clip(img, 0.0, 1.0)


def jitter_frame(truth, dy, dx, blur_sigma, noise_sigma, rng):
    """Apply sub-pixel shift, blur (seeing proxy), and additive read noise."""
    M = np.array([[1.0, 0.0, dx], [0.0, 1.0, dy]], dtype=np.float32)
    h, w = truth.shape
    shifted = cv2.warpAffine(truth, M, (w, h), flags=cv2.INTER_CUBIC,
                             borderMode=cv2.BORDER_REFLECT_101)
    if blur_sigma > 0:
        k = max(3, int(2 * round(3 * blur_sigma) + 1))
        shifted = cv2.GaussianBlur(shifted, (k, k), blur_sigma)
    shifted += rng.normal(0.0, noise_sigma, shifted.shape).astype(np.float32)
    return np.clip(shifted, 0.0, 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="test_data/synth_planet",
                    help="Output directory (relative to script's directory)")
    ap.add_argument("--n", type=int, default=24, help="Number of frames")
    ap.add_argument("--h", type=int, default=200, help="Frame height")
    ap.add_argument("--w", type=int, default=240, help="Frame width")
    ap.add_argument("--max-jitter", type=float, default=3.0,
                    help="Max per-frame translation in pixels")
    ap.add_argument("--noise", type=float, default=0.005)
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    truth = make_truth_frame(args.h, args.w)

    out_root = Path(__file__).resolve().parent / args.out
    out_root.mkdir(parents=True, exist_ok=True)

    # Quality profile: a handful of best frames, rest progressively blurred.
    # Linearly varying blur sigma in [0.2, 1.8] but shuffled so 'best' is not
    # at index 0.
    sigmas = np.linspace(0.2, 1.8, args.n)
    perm = rng.permutation(args.n)
    sigmas = sigmas[perm]

    truth_meta = []
    for i in range(args.n):
        dy = rng.uniform(-args.max_jitter, args.max_jitter)
        dx = rng.uniform(-args.max_jitter, args.max_jitter)
        f = jitter_frame(truth, dy, dx, sigmas[i], args.noise, rng)
        u16 = np.round(f * 65535.0).astype(np.uint16)
        path = out_root / f"frame_{i:03d}.png"
        Image.fromarray(u16, mode="I;16").save(path)
        truth_meta.append((i, dy, dx, float(sigmas[i])))

    # Persist the ground-truth jitter table for unit-test assertions.
    np.savez(out_root / "truth.npz",
             dy=np.array([m[1] for m in truth_meta], dtype=np.float64),
             dx=np.array([m[2] for m in truth_meta], dtype=np.float64),
             blur_sigma=np.array([m[3] for m in truth_meta], dtype=np.float64))

    print(f"Wrote {args.n} frames to {out_root}")


if __name__ == "__main__":
    main()
