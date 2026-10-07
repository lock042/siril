"""
Independent numpy reproduction of GaiaXPy's external-calibration flux path,
built from first principles (Hermite recurrence, dispersion spline, inverse-
bases & transformation matrices, response-based normalisation, BP/RP linear
merge). Validated bit-faithfully against GaiaXPy's calibrate() output to
machine precision (~1e-13 max relative error), and used as both:

  - the reference implementation that tools/bake_xp_design.py invokes to
    generate src/io/healpix/gaia_xp_design.c, and
  - the truth source that tools/xpcts_validate/diff_against_numpy.py compares
    Siril's C++ runtime against.

Re-run the bake/diff scripts after editing this file. Path defaults below
assume a sibling GaiaXPy checkout; override via env vars or function args.
"""
import os
import numpy as np
import pandas as pd
from scipy import interpolate
from os.path import join

# Constants (mirror gaiaxpy.core.nature, gaiaxpy.core.satellite)
C = 2.99792458e8
PLANCK = 6.62607004e-34
TELESCOPE_PUPIL_AREA = 0.7278
HC_NM = 1.0e9 * C * PLANCK
BP_WL_LOW, BP_WL_HIGH = 330.0, 643.0
RP_WL_LOW, RP_WL_HIGH = 635.0, 1020.0

CONFIG_PATH = os.environ.get("GAIAXPY_CONFIG", "../GaiaXPy/src/gaiaxpy/config")
FIXTURE = os.environ.get("GAIAXPY_FIXTURE",
                         "../GaiaXPy/tests/files/xp_continuous/XP_CONTINUOUS_RAW.csv")
SAMPLING = np.arange(336.0, 1021.0, 2.0)
assert len(SAMPLING) == 343


# ---------- Hermite functions (orthonormal "physicist's" form) -------------
SQRT_4_PI = np.pi ** -0.25  # = 1/pi^(1/4)


def hermite_function_vec(n_max, x):
    """Return shape (n_max, len(x)) with rows psi_0..psi_{n_max-1} evaluated at x.

    Uses the stable two-term recurrence from GaiaXPy's _hermite_function:
        psi_0 = pi^(-1/4) * exp(-x^2/2)
        psi_1 = pi^(-1/4) * sqrt(2) * x * exp(-x^2/2)
        psi_n = sqrt(2/n) * x * psi_{n-1} - sqrt((n-1)/n) * psi_{n-2}
    """
    x = np.asarray(x, dtype=np.float64)
    out = np.empty((n_max, x.size), dtype=np.float64)
    out[0] = SQRT_4_PI * np.exp(-0.5 * x * x)
    if n_max > 1:
        out[1] = SQRT_4_PI * np.sqrt(2.0) * x * np.exp(-0.5 * x * x)
    for n in range(2, n_max):
        c1 = np.sqrt(2.0 / n)
        c2 = -np.sqrt((n - 1.0) / n)
        out[n] = c1 * x * out[n - 1] + c2 * out[n - 2]
    return out


# ---------- External instrument model loaders ------------------------------
def load_dispersion(path):
    arr = np.genfromtxt(path, delimiter=",")
    return arr[0], arr[1]  # wavelength, pseudo-wavelength


def load_response(path):
    arr = np.genfromtxt(path, delimiter=",")
    return arr[0], arr[1]


def parse_paren_array(s):
    s = s.strip()
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    return np.fromstring(s, sep=",", dtype=np.float64)


def load_bases(path):
    """Parse single-row CSV with parenthesised arrays.
    Columns: nBases, pwlRangeMin, pwlRangeMax, normRangeMin, normRangeMax,
             nInverseBasesCoefficients, inverseBasesCoefficients,
             nTransformedBases, transformationMatrix
    """
    df = pd.read_csv(path)
    row = df.iloc[0]
    n_bases = int(row["nBases"])
    n_inv = int(row["nInverseBasesCoefficients"])
    n_tr = int(row["nTransformedBases"])
    inv = parse_paren_array(row["inverseBasesCoefficients"]).reshape(n_bases, n_inv)
    tr = parse_paren_array(row["transformationMatrix"]).reshape(n_bases, n_tr)
    return {
        "nBases": n_bases,
        "nInv": n_inv,
        "nTr": n_tr,
        "pwlRangeMin": float(row["pwlRangeMin"]),
        "pwlRangeMax": float(row["pwlRangeMax"]),
        "normRangeMin": float(row["normRangeMin"]),
        "normRangeMax": float(row["normRangeMax"]),
        "inverseBases": inv,
        "transformation": tr,
    }


# ---------- Build per-band design matrix on user grid ----------------------
def build_design_matrix(sampling, bases, disp_wl, disp_pwl, resp_wl, resp_r, weights):
    """Reproduce SampledBasisFunctions.from_external_instrument_model.

    Returns design_matrix of shape (nBases, len(sampling)).
    `weights` (length len(sampling)) zeroes out the Hermite evaluation where 0.
    """
    sampling = np.asarray(sampling, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)

    # 1. Linear rescale of pwl into the normalised Hermite domain.
    scale = ((bases["normRangeMax"] - bases["normRangeMin"]) /
             (bases["pwlRangeMax"] - bases["pwlRangeMin"]))
    offset = bases["normRangeMin"] - bases["pwlRangeMin"] * scale

    # 2. wl -> pwl via the dispersion spline (scipy splrep, default cubic).
    tck_disp = interpolate.splrep(disp_wl, disp_pwl, s=0)
    sampling_pwl = interpolate.splev(sampling, tck_disp, der=0)
    rescaled_pwl = sampling_pwl * scale + offset

    # 3. Hermite functions psi_0..psi_{nInv-1} at every rescaled position.
    H = hermite_function_vec(bases["nInv"], rescaled_pwl)  # (nInv, n_samples)
    H = H * (weights > 0)  # GaiaXPy zeroes contributions where weight==0

    # 4. design_matrix = transformation @ inverseBases @ H, shape (nBases, n_samples).
    interim = bases["inverseBases"] @ H                         # (nBases, n_samples)
    transformed = bases["transformation"] @ interim             # (nBases, n_samples)

    # 5. Response-based normalisation per output wavelength.
    tck_resp = interpolate.splrep(resp_wl, resp_r, s=0)
    r = interpolate.splev(sampling, tck_resp, der=0)
    norm = np.where(r > 0, HC_NM / (TELESCOPE_PUPIL_AREA * r * sampling), 0.0)
    design = transformed * norm  # broadcasts norm over rows

    return design


def merge_weights(sampling):
    """Reproduce calibrator.__create_merge for both bands."""
    s = np.asarray(sampling, dtype=np.float64)
    bp = np.where(s < RP_WL_LOW, 1.0,
         np.where(s > BP_WL_HIGH, 0.0,
                  1.0 - (s - RP_WL_LOW) / (BP_WL_HIGH - RP_WL_LOW)))
    rp = np.where(s < RP_WL_LOW, 0.0,
         np.where(s > BP_WL_HIGH, 1.0,
                  (s - RP_WL_LOW) / (BP_WL_HIGH - RP_WL_LOW)))
    return bp, rp


# ---------- Fixture loader (parses the GaiaXPy CSV format) -----------------
def load_fixture_coefficients(path):
    df = pd.read_csv(path)
    out = []
    for _, row in df.iterrows():
        bp_c = parse_paren_array(row["bp_coefficients"])
        rp_c = parse_paren_array(row["rp_coefficients"])
        out.append({
            "source_id": int(row["source_id"]),
            "bp": bp_c,
            "rp": rp_c,
            "bp_n_relevant": int(row["bp_n_relevant_bases"]),
            "rp_n_relevant": int(row["rp_n_relevant_bases"]),
        })
    return out


# ---------- Driver ---------------------------------------------------------
def main():
    bp_bases = load_bases(join(CONFIG_PATH, "bpC03_v375wi_bases.csv"))
    rp_bases = load_bases(join(CONFIG_PATH, "rpC03_v142r_bases.csv"))
    bp_disp = load_dispersion(join(CONFIG_PATH, "bpC03_v375wi_dispersion.csv"))
    rp_disp = load_dispersion(join(CONFIG_PATH, "rpC03_v142r_dispersion.csv"))
    bp_resp = load_response(join(CONFIG_PATH, "bpC03_v375wi_response.csv"))
    rp_resp = load_response(join(CONFIG_PATH, "rpC03_v142r_response.csv"))

    bp_w, rp_w = merge_weights(SAMPLING)

    bp_design = build_design_matrix(SAMPLING, bp_bases, *bp_disp, *bp_resp, bp_w)
    rp_design = build_design_matrix(SAMPLING, rp_bases, *rp_disp, *rp_resp, rp_w)

    sources = load_fixture_coefficients(FIXTURE)
    flux_test = np.zeros((len(sources), SAMPLING.size), dtype=np.float64)
    for i, s in enumerate(sources):
        bp_flux = s["bp"] @ bp_design
        rp_flux = s["rp"] @ rp_design
        flux_test[i] = bp_flux * bp_w + rp_flux * rp_w

    # Compare to reference.
    ref = np.load("/workspace/xpcts_work/gaiaxpy_reference.npz")
    flux_ref = ref["flux_full"]
    pos_ref = ref["pos"]
    assert np.array_equal(pos_ref, SAMPLING), "Reference grid mismatch"

    print(f"sources: {[s['source_id'] for s in sources]}")
    print(f"reference shape: {flux_ref.shape}, test shape: {flux_test.shape}")
    diff = flux_test - flux_ref
    rel = np.where(flux_ref != 0, diff / flux_ref, np.nan)
    abs_rel = np.abs(rel)
    print(f"\nPer-source results (full reconstruction, no truncation):")
    for i, s in enumerate(sources):
        finite = np.isfinite(abs_rel[i])
        print(f"  source {s['source_id']}:"
              f"  max |rel err| = {np.nanmax(abs_rel[i]):.3e}"
              f"  median |rel err| = {np.nanmedian(abs_rel[i]):.3e}"
              f"  max |abs err| = {np.max(np.abs(diff[i])):.3e}"
              f"  ref range = [{flux_ref[i].min():.3e}, {flux_ref[i].max():.3e}]")
        # Report worst sample location
        if finite.any():
            wmax = np.nanargmax(np.where(finite, abs_rel[i], -np.inf))
            print(f"    worst-rel sample idx {wmax} (lambda={SAMPLING[wmax]:.1f}nm),"
                  f" ref={flux_ref[i, wmax]:.3e}, test={flux_test[i, wmax]:.3e}")

    print(f"\nOverall: max |rel err| = {np.nanmax(abs_rel):.3e},"
          f" max |abs err| = {np.max(np.abs(diff)):.3e}")
    np.savez("/workspace/xpcts_work/numpy_repro.npz",
             source_id=np.array([s["source_id"] for s in sources]),
             flux=flux_test, pos=SAMPLING)


if __name__ == "__main__":
    main()
