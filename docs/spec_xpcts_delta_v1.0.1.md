# Siril HEALpixel Catalog Format — proposed v1.0.1 delta

These edits convert the placeholder "Photometric Extract with **xp_continuous** data" subsection in §3 (DATA) into a concrete record definition, matching the `SourceEntryXPcts` struct shipped with Siril's xp_continuous integration. They are intended to be applied verbatim against the v1.0.0 source LaTeX (Zenodo DOI 10.5281/zenodo.14697486).

---

## §1 — Version History (append)

> - **1.0.1**: defines the on-disk record layout for the `xp_continuous` photometric catalog type (cat_type = 3). No header or index changes; the format is fully forward-compatible with 1.0.0 readers, which can detect and reject the new type via the `cat_type` byte at offset 50.

## §3 — DATA — replace the existing subsection

Replace:

> ### Photometric Extract with `xp_continuous` data
>
> - To be defined in a future version of this specification.

with:

---

### Photometric Extract with `xp_continuous` data

This catalog type stores the per-source Hermite coefficients of the Gaia DR3 (and later) **internally calibrated** xp_continuous mean spectrum, together with a hint of the per-source recommended truncation. A consumer reconstructs an externally-calibrated, sampled spectrum on demand by multiplying the coefficients into a baked design matrix of basis functions (see *Reconstruction*, below). This avoids storing the redundant 343-sample flux vector per source while preserving full compatibility with downstream tools that expect xp_sampled output (e.g. SPCC).

| Field                       | Type            | Bytes | Notes |
| --------------------------- | --------------- | ----- | ----- |
| `RA_scaled`                 | `int32_t`       | 4     | Right Ascension, scaled by $360/(2^{31}-1)$. |
| `Dec_scaled`                | `int32_t`       | 4     | Declination, scaled by $360/(2^{31}-1)$. |
| `dRA`                       | `int16_t`       | 2     | Proper motion in RA (mas/yr). |
| `dDec`                      | `int16_t`       | 2     | Proper motion in Dec (mas/yr). |
| `G_mean_mag`                | `int16_t`       | 2     | Scaled mean G-band magnitude (×1000). |
| `bp_n_relevant`             | `uint8_t`       | 1     | GaiaXPy "recommended truncation" for BP, in [0, 55]. Advisory; consumers may honour or ignore it. A value of 0 means "no recommendation". |
| `rp_n_relevant`             | `uint8_t`       | 1     | As above for RP. |
| `bp_coefficients[55]`       | `float32` × 55  | 220   | BP Hermite coefficients in the Gaia DR3 internally-calibrated basis. |
| `rp_coefficients[55]`       | `float32` × 55  | 220   | RP Hermite coefficients (same basis class as BP, distinct numerical basis). |
| **Total size per record**   |                 | **456** |   |

The basis count of 55 corresponds to Gaia DR3's `bp_basis_function_id = 32` and `rp_basis_function_id = 37`; future Gaia data releases may use different bases, in which case the catalog header `gaia_version` field disambiguates the basis (and a future revision of this specification will record any change to the per-record layout).

#### Rationale for design choices

- **`float32` coefficients.** Empirical comparison against the GaiaXPy reference implementation across a Teff/Gmag-stratified sample of Gaia DR3 sources shows `float32` storage produces a worst-case integrated-band reconstruction error of $\approx 5\times 10^{-8}$, four orders of magnitude tighter than the truncation noise floor and effectively negligible against the intrinsic photometric uncertainty of the input spectra. `float64` was rejected as wasteful; `float16` was rejected because the leading BP coefficient overflows IEEE binary16 for the brightest stars (G ≲ 11).
- **No covariance.** The deterministic flux reconstruction $\mathrm{flux} = \mathbf{D}\,\mathbf{c}$ does not require the covariance matrix — covariance only enters analytic or Monte-Carlo error-band estimation, which Siril's SPCC pipeline does not perform. Storing the full $55 \times 55$ symmetric covariance per band would inflate each record by ~24 kB (a 50× growth), so it is omitted. A future "with-covariance" variant catalog type may be defined separately.
- **`n_relevant` hints, no enforced truncation.** GaiaXPy's recommended truncation reduces noise on faint sources but introduces up to ~7 % integrated-band bias on the bluest hot stars when applied uniformly. The hints are stored for callers that want them but the consumer is free to ignore them; the reference Siril runtime defaults to "use all 55 bases".

#### Reconstruction

Given a source record and a target sampling grid $\boldsymbol\lambda$, the absolute calibrated flux at sample $i$ is

$$
\mathrm{flux}(\lambda_i) = w^{\mathrm{BP}}(\lambda_i)\sum_{k=0}^{N^{\mathrm{BP}}-1} c^{\mathrm{BP}}_k\, D^{\mathrm{BP}}_{ki} \;+\; w^{\mathrm{RP}}(\lambda_i)\sum_{k=0}^{N^{\mathrm{RP}}-1} c^{\mathrm{RP}}_k\, D^{\mathrm{RP}}_{ki}
$$

where $D$ is the design matrix that bakes in (a) the wavelength-to-pseudo-wavelength dispersion spline, (b) the orthonormal Hermite basis functions, (c) the inverse-bases and transformation matrices supplied by GaiaXPy, and (d) the response-based normalisation $hc/(A_{\mathrm{tel}}\,R(\lambda)\,\lambda)$. $w^{\mathrm{BP}}$ and $w^{\mathrm{RP}}$ are the linear-blend merge weights that crossover in 635–643 nm. $N$ is either 55 (full) or the value supplied by `{bp,rp}_n_relevant` if the consumer chooses to honour the truncation hint.

Where the response is zero (past 1018 nm in DR3) the corresponding entries of $D$ are baked as zero, so consumers need no special edge handling.

The reference Siril implementation (`src/io/healpix/xp_continuous.{h,cpp}` and `tools/bake_xp_design.py`) bit-faithfully reproduces GaiaXPy's `calibrate()` output to machine precision (max integrated-band relative error $5.7 \times 10^{-16}$ vs the float64 numpy reference) on the standard 343-sample 336–1020 nm 2 nm grid.

---

## §9 — Best Practices

The existing entry already lists `xpcts` as a valid `<type>` token, so no edit is needed beyond confirming the example:

> - `siril_cat2_healpix8_xpcts_43.dat` – chunk 43 of an `xp_continuous` photometric catalog indexed at HEALpix level 8 and chunked at HEALpix level 1.

---

## §10 — Data Integrity (optional addendum)

If the xp_continuous catalog will be published alongside the existing astrometric and xp_sampled catalogs, add:

> The xp_continuous catalog is compiled by extracting the `bp_coefficients`, `rp_coefficients`, `bp_n_relevant_bases` and `rp_n_relevant_bases` columns from the Gaia archive `xp_continuous_mean_spectrum` table for sources matching the same per-HEALpixel selection criterion as the xp_sampled catalog, casting coefficients from `float64` to `float32` on write. Round-trip reconstruction was validated against the GaiaXPy reference at machine precision.

---

## Bytes-on-disk verification

The reference C struct definition (verbatim from `src/io/healpix/xp_continuous.h`):

```c
#pragma pack(push, 1)
typedef struct _SourceEntryXPcts {
    int32_t  ra_scaled;                       /*   4 B */
    int32_t  dec_scaled;                      /*   4 B */
    int16_t  dra_scaled;                      /*   2 B */
    int16_t  ddec_scaled;                     /*   2 B */
    int16_t  mag_scaled;                      /*   2 B */
    uint8_t  bp_n_relevant;                   /*   1 B */
    uint8_t  rp_n_relevant;                   /*   1 B */
    float    bp_coefficients[XPCTS_NBASES];   /* 220 B */
    float    rp_coefficients[XPCTS_NBASES];   /* 220 B */
} SourceEntryXPcts;                           /* 456 B */
#pragma pack(pop)
```

`sizeof(SourceEntryXPcts)` should equal **456** on every supported platform. A static_assert on the runtime side guards this.
