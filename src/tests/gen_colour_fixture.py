#!/usr/bin/env python3
# Regenerates the small synthetic 3-channel FITS used as a colour golden-regression
# fixture for the nlbayes harmonization (branch harmonize_img_t). Deterministic.
# Usage: python3 gen_colour_fixture.py [out.fit]
import sys
import numpy as np

W, H, C = 48, 48, 3
yy, xx = np.mgrid[0:H, 0:W].astype(np.float64)
img = np.zeros((C, H, W), dtype=np.float64)
img[0] = 0.3 + 0.3 * np.sin(xx / 5.0) + 0.1 * np.cos(yy / 3.0)
img[1] = 0.5 + 0.2 * np.cos(xx / 4.0) + 0.1 * np.sin(yy / 6.0)
img[2] = 0.4 + 0.25 * np.sin((xx + yy) / 7.0)
img += 0.05 * np.sin(13.0 * xx + 7.0 * yy)  # high-frequency "noise"
img = np.clip(img, 0.0, 1.0).astype('>f4')  # big-endian float32, planar (C,H,W)


def card(k, v):
    if isinstance(v, bool):
        v = 'T' if v else 'F'
    s = f"{k:<8}= {str(v):>20}"
    return (s + ' ' * (80 - len(s)))[:80]


cards = [card('SIMPLE', True), card('BITPIX', -32), card('NAXIS', 3),
         card('NAXIS1', W), card('NAXIS2', H), card('NAXIS3', C),
         ('END' + ' ' * 77)[:80]]
hdr = ''.join(cards).encode('ascii')
hdr += b' ' * ((2880 - len(hdr) % 2880) % 2880)
data = img.tobytes()
data += b'\x00' * ((2880 - len(data) % 2880) % 2880)

out = sys.argv[1] if len(sys.argv) > 1 else 'test_colour.fit'
open(out, 'wb').write(hdr + data)
print('wrote', out, W, H, C)
