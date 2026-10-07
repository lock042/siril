/*
 * xpcts_validate — round-trips a binary file of SourceEntryXPcts records
 * through xpcts_to_xpsampled() and writes the reconstructed flux to CSV.
 *
 * Input format:
 *   uint32_t  magic       (= 0x58504354 "XPCT")
 *   uint32_t  record_size (= sizeof(SourceEntryXPcts))
 *   uint32_t  count       (number of records)
 *   int32_t   truncation  (0 = full, -1 = USE_HINT, N = forced)
 *   <count> packed SourceEntryXPcts records back-to-back
 *
 * Output: CSV with one row per source: source_id_hi,source_id_lo,flux_343...
 * source_id is reconstructed from the 4-byte synthetic value we pack into
 * a spare slot — see pack_xpcts.py.
 */
#include "xp_continuous.h"
#include "gaia_xp_design.h"

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

static constexpr uint32_t MAGIC = 0x58504354u; /* "XPCT" little-endian */

int main(int argc, char **argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <input.bin> <output.csv>\n", argv[0]);
        return 2;
    }

    FILE *fin = std::fopen(argv[1], "rb");
    if (!fin) { std::perror("fopen input"); return 1; }

    uint32_t magic, record_size, count;
    int32_t truncation;
    if (std::fread(&magic, 4, 1, fin) != 1
        || std::fread(&record_size, 4, 1, fin) != 1
        || std::fread(&count, 4, 1, fin) != 1
        || std::fread(&truncation, 4, 1, fin) != 1) {
        std::fprintf(stderr, "short read on header\n"); return 1;
    }
    if (magic != MAGIC) {
        std::fprintf(stderr, "bad magic 0x%08x (expected 0x%08x)\n", magic, MAGIC);
        return 1;
    }
    if (record_size != sizeof(SourceEntryXPcts)) {
        std::fprintf(stderr, "record_size mismatch: file=%u runtime=%zu\n",
                     record_size, sizeof(SourceEntryXPcts));
        return 1;
    }

    /* We piggyback the source_id into the otherwise-unused ra_scaled+dec_scaled
     * fields purely for round-trip identification — the calibration math
     * only reads the coefficient arrays. */
    std::vector<SourceEntryXPcts> recs(count);
    if (std::fread(recs.data(), record_size, count, fin) != count) {
        std::fprintf(stderr, "short read on records\n"); return 1;
    }
    std::fclose(fin);

    FILE *fout = std::fopen(argv[2], "w");
    if (!fout) { std::perror("fopen output"); return 1; }
    std::fprintf(fout, "source_id_hi,source_id_lo");
    for (int i = 0; i < GAIA_XP_NSAMPLES; ++i) std::fprintf(fout, ",f%d", i);
    std::fputc('\n', fout);

    double flux[GAIA_XP_NSAMPLES];
    for (uint32_t r = 0; r < count; ++r) {
        xpcts_to_xpsampled(&recs[r], truncation, flux);
        /* Reconstruct 64-bit source_id from ra/dec slots. */
        uint64_t hi = (uint32_t)recs[r].ra_scaled;
        uint64_t lo = (uint32_t)recs[r].dec_scaled;
        std::fprintf(fout, "%llu,%llu",
                     (unsigned long long)hi, (unsigned long long)lo);
        for (int i = 0; i < GAIA_XP_NSAMPLES; ++i)
            std::fprintf(fout, ",%.17e", flux[i]);
        std::fputc('\n', fout);
    }
    std::fclose(fout);
    std::fprintf(stderr, "wrote %u sources to %s (truncation=%d)\n",
                 count, argv[2], truncation);
    return 0;
}
