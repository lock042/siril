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
	int16_t B;	// B mag times 1000 (abused for Teff in offline_gaia catalogue)
	int16_t V;	// mag times 1000
} deepStarData;

// Gaia source data structure with packed attributes to match file format
#pragma pack(push, 1)  // Ensure no padding between members
typedef struct _SourceEntryAstro {
	int32_t ra_scaled;   // 4 bytes
	int32_t dec_scaled;  // 4 bytes
	int16_t dra_scaled;  // 2 bytes, mas per year
	int16_t ddec_scaled; // 2 bytes, mas per year
	uint16_t teff;       // 2 bytes, K between 0-65535 (clipped)
	int16_t mag_scaled;  // 2 bytes
} SourceEntryAstro;

typedef struct _SourceEntryXPsamp {
	int32_t ra_scaled;   // 4 bytes
	int32_t dec_scaled;  // 4 bytes
	int16_t dra_scaled;  // 2 bytes, mas per year
	int16_t ddec_scaled; // 2 bytes, mas per year
	int16_t mag_scaled;  // 2 bytes
	// The remaining fields are only read for SPCC
	uint8_t fexpo;       // 1 byte
	int16_t flux[XPSAMPLED_LEN];   // 686 bytes: xp_sampled flux values
} SourceEntryXPsamp;
#pragma pack(pop)

void initialize_local_catalogues_paths();
gboolean local_kstars_available();
gboolean local_gaia_available();
gboolean local_catalogues_available();
siril_cat_index get_local_catalogue_index();

int siril_catalog_get_stars_from_local_catalogues(siril_catalogue *siril_cat);
gpointer write_trixels(gpointer p);
gpointer list_trixels(gpointer p);

#endif
