"""Pack the archive_xp_continuous.csv into a binary blob of SourceEntryXPcts
records that xpcts_validate.cpp can consume.

Pads BP/RP coefficient arrays with zeros if the source has fewer than 55
parameters. Stuffs source_id into the ra_scaled/dec_scaled slots (big-endian
upper/lower halves) purely so the C harness can echo it back for matching.
"""
import argparse
import struct
import sys
import numpy as np
import pandas as pd
from os.path import join, dirname

sys.path.insert(0, join(dirname(__file__), ".."))
from numpy_repro import parse_paren_array

XPCTS_NBASES = 55
RECORD_SIZE = 4 + 4 + 2 + 2 + 2 + 1 + 1 + 4 * XPCTS_NBASES + 4 * XPCTS_NBASES
assert RECORD_SIZE == 456, RECORD_SIZE  # struct invariant
MAGIC = 0x58504354  # "XPCT"


def pack_record(source_id, bp, rp, bp_nrel, rp_nrel):
    bp = np.zeros(XPCTS_NBASES, dtype=np.float32) if bp is None else bp.astype(np.float32)
    rp = np.zeros(XPCTS_NBASES, dtype=np.float32) if rp is None else rp.astype(np.float32)
    if bp.size < XPCTS_NBASES:
        bp = np.pad(bp, (0, XPCTS_NBASES - bp.size))
    if rp.size < XPCTS_NBASES:
        rp = np.pad(rp, (0, XPCTS_NBASES - rp.size))
    bp = bp[:XPCTS_NBASES]
    rp = rp[:XPCTS_NBASES]
    sid_hi = (source_id >> 32) & 0xFFFFFFFF
    sid_lo = source_id & 0xFFFFFFFF
    # Stuff source_id into ra/dec int32 slots (signed view doesn't matter — we
    # cast back to unsigned in the harness).
    ra = struct.unpack("<i", struct.pack("<I", sid_hi))[0]
    dec = struct.unpack("<i", struct.pack("<I", sid_lo))[0]
    header = struct.pack("<iihhhBB", ra, dec, 0, 0, 0,
                         min(bp_nrel, 255), min(rp_nrel, 255))
    return header + bp.tobytes() + rp.tobytes()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("out_bin")
    ap.add_argument("--truncation", type=int, default=0,
                    help="0 = full 55, -1 = use_hint, N = forced")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    blobs = []
    for _, row in df.iterrows():
        bp = parse_paren_array(row["bp_coefficients"])
        rp = parse_paren_array(row["rp_coefficients"])
        bp_nr = int(row["bp_n_relevant_bases"])
        rp_nr = int(row["rp_n_relevant_bases"])
        blobs.append(pack_record(int(row["source_id"]), bp, rp, bp_nr, rp_nr))

    with open(args.out_bin, "wb") as f:
        f.write(struct.pack("<IIIi", MAGIC, RECORD_SIZE, len(blobs), args.truncation))
        for b in blobs:
            assert len(b) == RECORD_SIZE
            f.write(b)
    print(f"wrote {len(blobs)} records ({len(blobs) * RECORD_SIZE} bytes data)"
          f" to {args.out_bin}, truncation={args.truncation}")


if __name__ == "__main__":
    main()
