/*
 * xp_continuous.cpp — runtime for Gaia DR3 XP-continuous external calibration.
 *
 * Pure GEMV: per band, flux_band[i] = sum_k (double)c[k] * D[k*N + i].
 * The design matrices D bake in everything else (Hermite evaluation,
 * dispersion spline, response normalisation), so the runtime has no
 * dependency beyond the static rodata blob in gaia_xp_design.c.
 */
#include "xp_continuous.h"
#include "gaia_xp_design.h"

#include <algorithm>

extern "C" void xpcts_to_xpsampled(const SourceEntryXPcts *src,
                                   int truncation,
                                   double *out)
{
    int bp_n;
    int rp_n;
    if (truncation == XPCTS_USE_HINT) {
        bp_n = src->bp_n_relevant ? src->bp_n_relevant : XPCTS_NBASES;
        rp_n = src->rp_n_relevant ? src->rp_n_relevant : XPCTS_NBASES;
    } else if (truncation == 0) {
        bp_n = XPCTS_NBASES;
        rp_n = XPCTS_NBASES;
    } else {
        bp_n = std::min(truncation, XPCTS_NBASES);
        rp_n = std::min(truncation, XPCTS_NBASES);
    }
    bp_n = std::min(bp_n, (int)XPCTS_NBASES);
    rp_n = std::min(rp_n, (int)XPCTS_NBASES);

    constexpr int N = GAIA_XP_NSAMPLES;
    static_assert(N == 343, "GAIA_XP_NSAMPLES must match XPSAMPLED_LEN");
    static_assert(GAIA_XP_NBASES == XPCTS_NBASES, "basis count mismatch");
    static_assert(sizeof(SourceEntryXPcts) == 456,
                  "SourceEntryXPcts on-disk size must match catalog format spec");

    /* Accumulate in double; coefficients are float32 in the struct. */
    for (int i = 0; i < N; ++i) {
        double bp_flux = 0.0;
        double rp_flux = 0.0;
        for (int k = 0; k < bp_n; ++k)
            bp_flux += (double)src->bp_coefficients[k] * gaia_xp_bp_design[k * N + i];
        for (int k = 0; k < rp_n; ++k)
            rp_flux += (double)src->rp_coefficients[k] * gaia_xp_rp_design[k * N + i];
        out[i] = bp_flux * gaia_xp_bp_merge[i] + rp_flux * gaia_xp_rp_merge[i];
    }
}
