#ifndef _CATALOGUES_H
#define _CATALOGUES_H

#include <glib.h>
#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "algos/photometry.h"
#include "io/siril_catalogues.h"

// the 16-byte struct
// this is the struct we effectively use, units are the correct ones
typedef struct {
	int32_t RA;	// hours times 1000000
	int32_t Dec;	// degrees times 100000
	int16_t dRA;	// in mas per year
	int16_t dDec;
	int16_t B;	// mag times 1000
	int16_t V;
} deepStarData;

// Gaia source data structure with packed attributes to match file format
#pragma pack(push, 1)  // Ensure no padding between members
struct SourceEntryAstro {
	uint64_t source_id;  // 8 bytes
	int32_t ra_scaled;   // 4 bytes
	int32_t dec_scaled;  // 4 bytes
	int16_t dra_scaled;  // 2 bytes
	int16_t ddec_scaled; // 2 bytes
	int16_t mag_scaled;  // 2 bytes
};

struct SourceEntryPhoto {
	uint64_t source_id;  // 8 bytes
	int32_t ra_scaled;   // 4 bytes
	int32_t dec_scaled;  // 4 bytes
	int16_t dra_scaled;  // 2 bytes
	int16_t ddec_scaled; // 2 bytes
	int16_t mag_scaled;  // 2 bytes
	// The remaining fields are only read for SPCC
	uint8_t fexpo;       // 1 byte
	int16_t flux[343];   // 686 bytes: xp_sampled flux values
};
#pragma pack(pop)

void initialize_local_catalogues_paths();
gboolean local_catalogues_available();

int siril_catalog_get_stars_from_local_catalogues(siril_catalogue *siril_cat);
gpointer write_trixels(gpointer p);
gpointer list_trixels(gpointer p);

#endif
