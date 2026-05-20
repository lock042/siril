"""Compare xpcts_validate's CSV output against numpy_repro reconstruction.

Three reconstruction passes are checked:
  truncation=0  : full 55 bases per band (the production default)
  truncation=-1 : honour bp/rp_n_relevant_bases from struct
  truncation=30 : forced low N (stress-test)

Reports per-source band-integrated relative error matching the metric used
in the dtype sweep, plus per-sample max in the SPCC range.
"""
import sys
import numpy as np
import pandas as pd
from os.path import join, dirname

sys.path.insert(0, join(dirname(__file__), ".."))
from numpy_repro import (SAMPLING, build_design_matrix, merge_weights,
                         load_bases, load_dispersion, load_response,
                         parse_paren_array, CONFIG_PATH)
from os.path import join as pjoin

BANDS = [("B",430,490),("V",500,580),("R",600,680),("I",750,870),
         ("Ha",655,660),("OIII",499,502),("BPw",380,640),("RPw",640,950),("Gw",400,850)]

def integ(flux, lo, hi):
    m = (SAMPLING >= lo) & (SAMPLING <= hi)
    return np.trapezoid(flux[m], SAMPLING[m])


def numpy_reference(csv_path, design):
    bp_d, rp_d, bp_w, rp_w = design
    df = pd.read_csv(csv_path)
    out = {}
    nrel = {}
    for _, row in df.iterrows():
        sid = int(row["source_id"])
        # Match what pack_xpcts.py does: cast coefficients to float32 (and
        # back to float64 for accumulation) so the reference matches what
        # the C++ runtime sees from the struct.
        bp = parse_paren_array(row["bp_coefficients"]).astype(np.float32).astype(np.float64)
        rp = parse_paren_array(row["rp_coefficients"]).astype(np.float32).astype(np.float64)
        bp = np.pad(bp, (0, max(0, 55 - bp.size)))[:55]
        rp = np.pad(rp, (0, max(0, 55 - rp.size)))[:55]
        out[sid] = (bp, rp)
        nrel[sid] = (int(row["bp_n_relevant_bases"]), int(row["rp_n_relevant_bases"]))
    return out, nrel, design


def reconstruct(bp, rp, bp_n, rp_n, design):
    bp_d, rp_d, bp_w, rp_w = design
    bp = bp.copy(); rp = rp.copy()
    if bp_n is not None: bp[bp_n:] = 0
    if rp_n is not None: rp[rp_n:] = 0
    return (bp @ bp_d) * bp_w + (rp @ rp_d) * rp_w


def main():
    if len(sys.argv) != 3:
        print("usage: diff_against_numpy.py <archive.csv> <c_output.csv>")
        sys.exit(2)
    archive_csv, c_csv = sys.argv[1], sys.argv[2]

    bp_bases = load_bases(pjoin(CONFIG_PATH, "bpC03_v375wi_bases.csv"))
    rp_bases = load_bases(pjoin(CONFIG_PATH, "rpC03_v142r_bases.csv"))
    bp_disp = load_dispersion(pjoin(CONFIG_PATH, "bpC03_v375wi_dispersion.csv"))
    rp_disp = load_dispersion(pjoin(CONFIG_PATH, "rpC03_v142r_dispersion.csv"))
    bp_resp = load_response(pjoin(CONFIG_PATH, "bpC03_v375wi_response.csv"))
    rp_resp = load_response(pjoin(CONFIG_PATH, "rpC03_v142r_response.csv"))
    bp_w, rp_w = merge_weights(SAMPLING)
    bp_d = np.nan_to_num(build_design_matrix(SAMPLING, bp_bases, *bp_disp, *bp_resp, bp_w))
    rp_d = np.nan_to_num(build_design_matrix(SAMPLING, rp_bases, *rp_disp, *rp_resp, rp_w))
    design = (bp_d, rp_d, bp_w, rp_w)

    coeffs, nrel, _ = numpy_reference(archive_csv, design)

    # Read C output.
    c_df = pd.read_csv(c_csv)
    c_df["source_id"] = (c_df["source_id_hi"].astype(np.uint64).to_numpy() << np.uint64(32)) | c_df["source_id_lo"].astype(np.uint64).to_numpy()
    flux_cols = [f"f{i}" for i in range(343)]

    # Determine truncation regime from CSV filename (lazy but cheap).
    if "trunc0" in c_csv:
        regime = ("truncation=0 (full)", lambda sid: (None, None))
    elif "truncm1" in c_csv:
        regime = ("truncation=-1 (use_hint)", lambda sid: nrel[sid])
    elif "trunc30" in c_csv:
        regime = ("truncation=30 (forced)", lambda sid: (30, 30))
    else:
        regime = ("unknown", lambda sid: (None, None))
    print(f"\n=== {regime[0]} — {len(c_df)} sources ===")

    persamp_max, persamp_p99, persamp_med = [], [], []
    band_max_list = []
    spcc_mask = (SAMPLING >= 380) & (SAMPLING <= 950)
    for _, row in c_df.iterrows():
        sid = int(row["source_id"])
        if sid not in coeffs:
            continue
        bp, rp = coeffs[sid]
        bp_n, rp_n = regime[1](sid)
        ref = reconstruct(bp, rp, bp_n, rp_n, design)
        c_flux = row[flux_cols].astype(np.float64).values
        # per-sample
        m = spcc_mask & (np.abs(ref) > 0)
        rel = np.abs((c_flux - ref)[m] / ref[m])
        persamp_max.append(rel.max()); persamp_p99.append(np.percentile(rel, 99)); persamp_med.append(np.median(rel))
        # band-integrated
        ref_int = np.array([integ(ref, lo, hi) for _, lo, hi in BANDS])
        c_int = np.array([integ(c_flux, lo, hi) for _, lo, hi in BANDS])
        bre = np.max(np.abs((c_int - ref_int) / ref_int))
        band_max_list.append(bre)

    persamp_max = np.array(persamp_max); band_max_arr = np.array(band_max_list)
    print(f"  per-sample max |rel err|: max-over-sources = {persamp_max.max():.3e},"
          f" p99-over-sources = {np.percentile(persamp_max,99):.3e}")
    print(f"  band-integrated:          max = {band_max_arr.max():.3e},"
          f" p99 = {np.percentile(band_max_arr,99):.3e},"
          f" med = {np.median(band_max_arr):.3e}")
    if band_max_arr.max() < 1e-5:
        print("  PASS: C++ runtime matches numpy reference to <1e-5 band error")
    elif band_max_arr.max() < 1e-3:
        print("  OK: <0.1% band error (within fp32 expectations for full bases)")
    else:
        print("  WARN: band error >0.1% — check truncation regime, this is expected for forced=30")


if __name__ == "__main__":
    main()
