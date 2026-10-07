/*
 * xp_continuous.h — Gaia DR3 XP continuous spectra: catalogue struct and
 * external-calibration runtime.
 *
 * SourceEntryXPcts holds the per-source BP/RP Hermite coefficients (float32)
 * plus the GaiaXPy-recommended truncation hints. xpcts_to_xpsampled() is the
 * runtime that reconstructs the absolute, sampled spectrum on Siril's
 * 343-point 2 nm grid (336..1020 nm) — drop-in replacement for the
 * SourceEntryXPsamp.flux[] array consumed by SPCC.
 *
 * Calibration uses the design matrices baked in gaia_xp_design.c.
 */
#ifndef SIRIL_XP_CONTINUOUS_H
#define SIRIL_XP_CONTINUOUS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define XPCTS_NBASES 55

#pragma pack(push, 1)
typedef struct _SourceEntryXPcts {
    int32_t  ra_scaled;                       /* hours * 1e6           */
    int32_t  dec_scaled;                      /* degrees * 1e5         */
    int16_t  dra_scaled;                      /* mas/yr                */
    int16_t  ddec_scaled;
    int16_t  mag_scaled;                      /* mag * 1000 (G band)   */
    uint8_t  bp_n_relevant;                   /* GaiaXPy truncation    */
    uint8_t  rp_n_relevant;                   /* hint, advisory only   */
    float    bp_coefficients[XPCTS_NBASES];   /* 220 B                 */
    float    rp_coefficients[XPCTS_NBASES];   /* 220 B                 */
} SourceEntryXPcts;                           /* total: 456 B          */
#pragma pack(pop)

/*
 * Reconstruct the absolute calibrated sampled spectrum for one source.
 *
 * `truncation` controls how many coefficients per band are used:
 *   0                    -> use all XPCTS_NBASES (recommended default)
 *   XPCTS_USE_HINT       -> honour src->bp_n_relevant / rp_n_relevant
 *   1..XPCTS_NBASES      -> force the same N for both bands
 *
 * Output `out` is exactly XPSAMPLED_LEN (=343) doubles in W*nm^-1*m^-2,
 * suitable as a drop-in for SourceEntryXPsamp.flux after the standard
 * fexpo decode.
 *
 * Past lambda=1018 nm the Gaia external calibration response is zero, so
 * those samples are baked as 0.0 and need no special handling.
 */
#define XPCTS_USE_HINT (-1)

void xpcts_to_xpsampled(const SourceEntryXPcts *src,
                        int truncation,
                        double *out /* [343] */);

#ifdef __cplusplus
}
#endif

#endif
