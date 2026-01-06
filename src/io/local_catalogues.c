/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/* This file serves as an interface to KStars binary catalogue format. It is a
 * C implementation derived from the KStars tool nomadbinfiletester.c, simpler
 * than the actual C++ implementation of KStars for catalogue access.
 * The structure of the file is documented at
 * https://api.kde.org/kstars/html/Stars.html#BinaryFormat
 * See also ../../subprojects/htmesh/README for complete information.
 */

#include "io/kstars/binfile.h"
#include "io/kstars/htmesh_wrapper.h"
#include "core/siril.h"
#include "core/arithm.h" // for half_to_float()
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "algos/photometry.h"
#include "algos/siril_wcs.h"
#include "registration/matching/degtorad.h"
#include "siril_catalogues.h"
#include "local_catalogues.h"
#include "io/healpix/healpix_cat.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef fseek64
#ifdef _WIN32
#define fseek64 _fseeki64
#else
#define fseek64 fseeko
#endif
#endif

/* Used for sanitizing catalogue input. 1<<20 seems reasonable as a maximum limit
 * both for the number of trixels in a catalogue and for the number of stars in a
 * single trixel. These values are mostly used to sanitize the data to prevent
 * asking malloc() to allocate nonsensical amounts of memory.
 */
#define MAX_NUM_TRIXELS 1<<20
#define MAX_STARS_PER_TRIXEL 1<<20

#define NAMEDSTARS_DAT "~/.local/share/kstars/namedstars.dat"
#define UNNAMEDSTARS_DAT "~/.local/share/kstars/unnamedstars.dat"
#define TYCHOSTARS_DAT "~/.local/share/kstars/deepstars.dat"
#define NOMAD_DAT "~/.local/share/kstars/USNO-NOMAD-1e8.dat"
#define GAIA_ASTRO_DAT "~/.local/share/siril/gaia_astrometric.dat"
#define GAIA_PHOTO_DAT "~/.local/share/siril/gaia_photometric.dat"

// All these are searched in the legacy local catalogue search, so we don't add GAIA_ASTRO_DAT to this
// array as it is handled separately
const char *default_catalogues_paths[] = { NAMEDSTARS_DAT, UNNAMEDSTARS_DAT, TYCHOSTARS_DAT, NOMAD_DAT, GAIA_ASTRO_DAT, GAIA_PHOTO_DAT };

#define DEBUG_LOCALCAT 0 // set to 1 to print out trixel star search verbose
#define NTRIXELS 2
#define NMAG 19

#define cat_debug_print(fmt, ...) \
	do { if (DEBUG_LOCALCAT) fprintf(stdout, fmt, ##__VA_ARGS__); } while (0)

struct catalogue_index {
	uint32_t trixelID;
	uint32_t offset;	// file size is < 4G, it's fine
	uint32_t nrecs;
};

struct catalogue_file {
	FILE *f;
	dataElement *de;
	uint16_t nfields;
	long index_offset, data_offset;
	char byteswap;
	uint32_t ntrixels;
	int8_t HTM_Level;
	struct catalogue_index *indices;
};

struct expansion_field {
	int16_t faint_mag;	// magnitude times 1000
	uint8_t HTM_level;
	uint16_t max_stars;
};

// alternative 16-byte struct for plate solving, replaces proper motion by distance
typedef struct {
	int32_t RA;	// hours times 1000000
	int32_t Dec;	// degrees times 100000
	float distance;	// for catalogue queries, distance from target, in arcmin
	int16_t B;	// mag times 1000
	int16_t V;
} deepStarData_dist;

// the 32-byte struct
typedef struct {
	int32_t RA;
	int32_t Dec;
	int32_t dRA;
	int32_t dDec;
	int32_t parallax;
	uint32_t HD;	// unsigned 32-bit Henry Draper Index.
	int16_t mag;
	int16_t bv_index;
	char spec_type[2];
	char flags;
	char unused;
} shallowStarData;

struct top_header {
	char description[124];
	int16_t endian_id;
	uint8_t version_no;
};

static struct catalogue_file *catalogue_read_header(FILE *f);
static int read_trixels_of_target(double ra, double dec, double radius, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);
static int read_trixels_by_ID(int ID, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);
static int read_trixel(int trixel, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);

static void bswap_stardata(deepStarData *stardata) {
    stardata->RA = bswap_32(stardata->RA);
    stardata->Dec = bswap_32(stardata->Dec);
    stardata->dRA = bswap_16(stardata->dRA);
    stardata->dDec = bswap_16(stardata->dDec);
    stardata->B = bswap_16(stardata->B);
    stardata->V = bswap_16(stardata->V);
}

/* returns the complete list of stars for a catalogue's list of trixels */
static int read_trixels_from_catalogue(const char *path, double ra, double dec, double radius, deepStarData **trixel_stars, uint32_t *trixel_nb_stars) {
	cat_debug_print("reading data from catalogue %s\n", path);
	FILE *f = g_fopen(path, "rb");
	if (!f) {
		siril_log_message(_("Could not open local catalogue %s\n"), path);
		return 1;
	}

	struct catalogue_file *cat = catalogue_read_header(f);
	if (!cat) {
		siril_log_message(_("Failed to read the local catalogue %s\n"), path);
		fclose(f);
		return 1;
	}

	//if (read_trixel_of_target(ra, dec, cat, &trixel_stars, &trixel_nb_stars)) {
	if (read_trixels_of_target(ra, dec, radius, cat, trixel_stars, trixel_nb_stars)) {
		free(cat);
		fclose(f);
		return 1;
	}
	fclose(f);
	free(cat->indices);
	free(cat->de);
	free(cat);
	return 0;
}

/* returns the complete list of stars for a catalogue's list of trixels */
static int read_trixelID_from_catalogue(const char *path, int ID, deepStarData **trixel_stars, uint32_t *trixel_nb_stars) {
	cat_debug_print("reading data from catalogue %s\n", path);
	if (ID < 0 || ID > 512) {
		siril_log_message(_("Wrong trixel ID\n"));
		return 1;
	}
	FILE *f = g_fopen(path, "rb");
	if (!f) {
		siril_log_message(_("Could not open local catalogue%s\n"), path);
		return 1;
	}

	struct catalogue_file *cat = catalogue_read_header(f);
	if (!cat) {
		siril_log_message(_("Failed to read the local catalogue %s\n"), path);
		fclose(f);
		return 1;
	}

	if (read_trixels_by_ID(ID, cat, trixel_stars, trixel_nb_stars)) {
		free(cat);
		fclose(f);
		return 1;
	}
	fclose(f);
	free(cat->indices);
	free(cat->de);
	free(cat);
	return 0;
}

static struct catalogue_file *catalogue_read_header(FILE *f) {
	struct top_header top;
	if (!fread(&top, 127, 1, f)) {
		siril_debug_print("error reading top header\n");
		return NULL;
	}
	top.description[123] = '\0'; // we may lose a character, but it shouldn't be bad

	//siril_debug_print("Opened catalogue '%s'\n", top.description);

	struct catalogue_file *cat = calloc(1, sizeof(struct catalogue_file));
	if (!cat) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	if (top.endian_id != 0x4B53) {
		if (top.endian_id != 0x534B) {
			siril_debug_print("invalid endian ID in header\n");
			free(cat);
			return NULL;
		}
		siril_debug_print("Byteswapping required\n");
		cat->byteswap = 1;
	}
	//siril_debug_print("Version number: %d\n", top.version_no);

	/* reading the fields */
	if (!fread(&cat->nfields, 2, 1, f)) {
		siril_debug_print("error reading fields\n");
		free(cat);
		return NULL;
	}
	if (cat->nfields != 6 && cat->nfields != 11) {
		siril_debug_print("the number of field is not recognized as one of the main catalogues\n");
		free(cat);
		return NULL;
	}

	if (cat->byteswap)
		cat->nfields = bswap_16(cat->nfields);
	//siril_debug_print("%d fields reported:\n", cat->nfields);
	cat->de = malloc(cat->nfields * sizeof(dataElement));
	if (!cat->de) {
		PRINT_ALLOC_ERR;
		free(cat);
		return NULL;
	}

	for (int i = 0; i < cat->nfields; ++i) {
		if (!fread(&(cat->de[i]), sizeof(dataElement), 1, f)) {
			siril_debug_print("error reading field descriptor %d\n", i);
			free(cat->de);
			free(cat);
			return NULL;
		}
		if (cat->byteswap)
			cat->de[i].scale = bswap_32(cat->de[i].scale);
		//displayDataElementDescription(&(cat->de[i]));
	}

	/* reading the trixel index table */
	if (!fread(&cat->ntrixels, 4, 1, f)) {
		siril_debug_print("error reading number of trixels\n");
		free(cat->de);
		free(cat);
		return NULL;
	}
	if (cat->byteswap)
		cat->ntrixels = bswap_32(cat->ntrixels);

	uint32_t ntrixels = cat->ntrixels;
	cat->HTM_Level = -1;
	while (ntrixels >>= 2) ++cat->HTM_Level;
	/* 512 for level 3, 32768 for NOMAD level 6 */
	//siril_debug_print("Number of trixels reported = %d (levels: %d)\n",
	//		cat->ntrixels, cat->HTM_Level);

	// ntrixels is tainted (read from external file): sanitize values
	// Maximum of 2^20 (MAX_NUM_TRIXELS) currently seems reasonable (and allows plenty
	// of headroom)
	if (cat->ntrixels < 1 || cat->ntrixels > MAX_NUM_TRIXELS) {
		siril_log_color_message(_("Error: number of trixels reported by file is out of limits.\n"), "red");
		free(cat->de);
		free(cat);
		return NULL;
	}

	cat->index_offset = ftell(f);

	cat->indices = malloc(cat->ntrixels * sizeof(struct catalogue_index));
	if (!cat->indices) {
		PRINT_ALLOC_ERR;
		free(cat->de);
		free(cat);
		return NULL;
	}

	if (fread(cat->indices, sizeof(struct catalogue_index), cat->ntrixels, f) < cat->ntrixels) {
		siril_debug_print("unexpected read failure in index table\n");
		free(cat->indices);
		free(cat->de);
		free(cat);
		return NULL;
	}
	for (uint32_t i = 0; i < cat->ntrixels; ++i) {
		struct catalogue_index *current = cat->indices + i;
		if (cat->byteswap) {
			current->trixelID = bswap_32(current->trixelID);
			current->offset = bswap_32(current->offset);
			current->nrecs = bswap_32(current->nrecs);
		}
		if (current->trixelID != i) {
			siril_debug_print("expected trixel ID of %d, got %u\n", i, current->trixelID);
			// if this is not right, we won't be able to use the indexing of the index
			free(cat->indices);
			free(cat->de);
			free(cat);
			return NULL;
		}
		if (current->nrecs > MAX_STARS_PER_TRIXEL) {
			siril_debug_print("catalogue claims excessive number (%d) of stars in trixel ID %u. Potential data error or corruption\n", current->nrecs, current->trixelID);
			// memory safety as this value is used for memory allocation
			free(cat->indices);
			free(cat->de);
			free(cat);
			return NULL;
		}
	}

	cat->data_offset = ftell(f);
	cat->f = f;

	struct expansion_field header;
	if (!fread(&header, sizeof(struct expansion_field), 1, cat->f)) {
		siril_debug_print("failed to read the trixel header\n");
	}
	else {
		if (cat->byteswap) {
			header.faint_mag = bswap_16(header.faint_mag);
			header.max_stars = bswap_16(header.max_stars);
		}
		//siril_debug_print("faint magnitude: %.2f\n", header.faint_mag / 1000.0);
	}

	//siril_debug_print("read the trixel index table, header read complete\n");
	return cat;
}

/* radius is in degrees */
static int read_trixels_of_target(double ra, double dec, double radius, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars) {
	int *trixels;
	int nb_trixels;
	if (get_htm_indices_around_target(ra, dec, radius, cat->HTM_Level, &trixels, &nb_trixels)) {
		*stars = NULL;
		*nb_stars = 0;
		return 0;
	}
	cat_debug_print("trixel search found %d trixels\n", nb_trixels);
	deepStarData **stars_list;
	uint32_t *nb_stars_list;
	stars_list = malloc(nb_trixels * sizeof(deepStarData *));
	if (!stars_list) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	nb_stars_list = malloc(nb_trixels * sizeof(uint32_t));
	if (!nb_stars_list) {
		free(stars_list);
		PRINT_ALLOC_ERR;
		return -1;
	}
	int retval = 0;
	uint32_t total_star_count = 0;
	for (int i = 0; i < nb_trixels; i++) {
		retval = read_trixel(trixels[i], cat, stars_list + i, nb_stars_list + i);
		if (retval) break;
		cat_debug_print("trixel %d (%d) contained %u stars\n", i, trixels[i], nb_stars_list[i]);
		total_star_count += nb_stars_list[i];
	}
	free(trixels);

	if (!retval) {
		// aggregate
		*stars = malloc(total_star_count * sizeof(deepStarData));
		if (!*stars) {
			free(stars_list);
			free(nb_stars_list);
			PRINT_ALLOC_ERR;
			return -1;
		}
		*nb_stars = total_star_count;
		uint32_t offset = 0;
		for (int i = 0; i < nb_trixels; i++) {
			memcpy(*stars + offset, stars_list[i], nb_stars_list[i] * sizeof(deepStarData));
			offset += nb_stars_list[i];
			free(stars_list[i]);
		}
	}
	free(stars_list);
	free(nb_stars_list);
	return 0;
}

static int read_trixels_by_ID(int ID, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars) {
	int *trixels;
	int nb_trixels;

	if (cat->HTM_Level == 3) {
		trixels = malloc(1 * sizeof(int));
		trixels[0] = ID;
		nb_trixels = 1;
	} else if (cat->HTM_Level == 6) {
		trixels = malloc(64 * sizeof(int));
		nb_trixels = 64;
		int base = 64 * ID;
		for (int i = 0; i < 64; i++) {
			trixels[i] = base + i;
		}
	} else {
		*stars = NULL;
		*nb_stars = 0;
		return 1;
	}
	cat_debug_print("trixel search found %d trixels\n", nb_trixels);
	deepStarData **stars_list;
	uint32_t *nb_stars_list;
	stars_list = malloc(nb_trixels * sizeof(deepStarData *));
	if (!stars_list) {
		free(trixels);
		PRINT_ALLOC_ERR;
		return -1;
	}
	nb_stars_list = malloc(nb_trixels * sizeof(uint32_t));
	if (!nb_stars_list) {
		free(trixels);
		free(stars_list);
		PRINT_ALLOC_ERR;
		return -1;
	}
	int retval = 0;
	uint32_t total_star_count = 0;
	for (int i = 0; i < nb_trixels; i++) {
		// double ra1, dec1, ra2, dec2, ra3, dec3;
		// get_vertices_for_index(trixels[i], cat->HTM_Level, &ra1, &dec1, &ra2, &dec2, &ra3, &dec3);
		retval = read_trixel(trixels[i], cat, stars_list + i, nb_stars_list + i);
		if (retval) break;
		cat_debug_print("trixel %d (%d) contained %u stars\n", i, trixels[i], nb_stars_list[i]);
		total_star_count += nb_stars_list[i];
	}
	free(trixels);

	if (!retval) {
		// aggregate
		*stars = malloc(total_star_count * sizeof(deepStarData));
		if (!*stars) {
			free(stars_list);
			free(nb_stars_list);
			PRINT_ALLOC_ERR;
			return -1;
		}
		*nb_stars = total_star_count;
		uint32_t offset = 0;
		for (int i = 0; i < nb_trixels; i++) {
			memcpy(*stars + offset, stars_list[i], nb_stars_list[i] * sizeof(deepStarData));
			offset += nb_stars_list[i];
			free(stars_list[i]);
		}
	}
	free(stars_list);
	free(nb_stars_list);
	return 0;
}

/* this function reads all stars of a trixel and returns them as deep star data
 * (the 16-byte struct) even if this is a 32-byte catalog */
static int read_trixel(int trixel, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars) {
	struct catalogue_index *index = cat->indices + trixel;
	if (index->trixelID != trixel) {
		siril_debug_print("INDEX IS WRONG, trixel ID did not match\n");
		return 1;
	}

	//siril_debug_print("offset for trixel %u: %u\n", trixel, index->offset);
	if (fseek64(cat->f, index->offset, SEEK_SET)) {
		siril_debug_print("failed to seek to the trixel offset %u\n", index->offset);
		return 1;
	}

	deepStarData *trix_stars = malloc(index->nrecs * sizeof(deepStarData));
	if (!trix_stars) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	*stars = trix_stars;
	*nb_stars = index->nrecs;

	if (cat->nfields == 6) {
		if (fread(trix_stars, sizeof(deepStarData), index->nrecs, cat->f) < index->nrecs) {
			siril_debug_print("error reading trixel data\n");
			free(trix_stars);
			return 1;
		}
	} else {
		shallowStarData *read_stars = malloc(index->nrecs * sizeof(shallowStarData));
		if (!read_stars) {
			PRINT_ALLOC_ERR;
			free(trix_stars);
			return 1;
		}
		if (fread(read_stars, sizeof(shallowStarData), index->nrecs, cat->f) < index->nrecs) {
			siril_debug_print("error reading trixel data\n");
			free(read_stars);
			free(trix_stars);
			return 1;
		}
		for (uint32_t i = 0; i < index->nrecs; ++i) {
			trix_stars[i].RA = read_stars[i].RA;
			trix_stars[i].Dec = read_stars[i].Dec;
			trix_stars[i].dRA = read_stars[i].dRA / 10;
			trix_stars[i].dDec = read_stars[i].dDec / 10;
			trix_stars[i].V = read_stars[i].mag * 10;
			trix_stars[i].B = (read_stars[i].bv_index + read_stars[i].mag) * 10;
			/* the BV - mag trick may not be correct, but we only use the V and
			 * B fields to compute B-V anyway */
			/* the * 10 trick is because the scaling is 100 for this catalogue,
			 * while it is 1000 for the deep star catalogue, but we don't know
			 * where they come from after it has returned */
		}
		free(read_stars);
	}

	if (cat->byteswap) {
		for (uint32_t i = 0; i < index->nrecs; ++i) {
			bswap_stardata(trix_stars + i);
		}
	}

	return 0;
}

/* Fetches the stars from local catalogs in the form of deepStarData stars
   Filters out stars w/o photometric data if photometric is TRUE
*/
static int get_raw_stars_from_local_catalogues(double target_ra, double target_dec, double radius,
		float max_mag, gboolean photometric, deepStarData **stars, uint32_t *nb_stars) {
	int nb_catalogues = 4;
	deepStarData **catalogue_stars = malloc(nb_catalogues * sizeof(deepStarData *));
	uint32_t *catalogue_nb_stars = malloc(nb_catalogues * sizeof(uint32_t));
	uint32_t total_nb_stars = 0;
	int retval = 0, catalogue = 0;
	radius /= 60.; // converting to degrees
	double radius_h = pow(sin(0.5 * radius * DEGTORAD), 2);
	int16_t mag_threshold = (int16_t)roundf(max_mag * 1000.f);

	siril_debug_print("looking for stars in local catalogues for target %f, %f, radius %f, magnitude %.2f, photometric: %d\n", target_ra, target_dec, radius, max_mag, photometric);
	for (; catalogue < nb_catalogues; catalogue++) {
		if (catalogue == 3 && max_mag < 12.0) {
			cat_debug_print("not querying NOMAD for this limit magnitude\n");
			catalogue_stars[catalogue] = NULL;
			catalogue_nb_stars[catalogue] = 0;
			continue;
		}

		retval = read_trixels_from_catalogue(com.pref.catalogue_paths[catalogue],
				target_ra, target_dec, radius,
				catalogue_stars + catalogue, catalogue_nb_stars + catalogue);
		if (retval)
			break;
		total_nb_stars += catalogue_nb_stars[catalogue];
		cat_debug_print("%d raw stars from catalogue %d\n", catalogue_nb_stars[catalogue], catalogue);
	}

	if (catalogue == nb_catalogues) {
		// aggregate and filter
		*stars = malloc(total_nb_stars * sizeof(deepStarData));
		if (!*stars) {
			PRINT_ALLOC_ERR;
			retval = 1;
		} else {
			uint32_t j = 0;
			for (catalogue = 0; catalogue < nb_catalogues; catalogue++) {
				deepStarData *cat_stars = catalogue_stars[catalogue];
				uint32_t cat_nb_stars = catalogue_nb_stars[catalogue];

				for (uint32_t i = 0; i < cat_nb_stars ; i++) {
					if (cat_stars[i].V > mag_threshold)
						continue;
					if (cat_stars[i].B >= 30000 && photometric)
						continue;
					// catalogue has RA in hours, hence the x15
					double ra = cat_stars[i].RA * .000015;
					double dec = cat_stars[i].Dec * .00001;
					double dist_h = compute_coords_distance_h(ra, dec, target_ra, target_dec);
					if (dist_h > radius_h)
						continue;
					deepStarData *sdd = &cat_stars[i];
					(*stars)[j] = *sdd;
					j++;
				}
				free(catalogue_stars[catalogue]);
			}

			*nb_stars = j;
			*stars = realloc(*stars, j * sizeof(deepStarData));
		}
	}

	free(catalogue_nb_stars);
	free(catalogue_stars);
	return retval;
}

static int get_raw_stars_from_local_catalogues_byID(int ID, float max_mag, gboolean photometric, deepStarData **stars, uint32_t *nb_stars) {
	int nb_catalogues = sizeof(default_catalogues_paths) / sizeof(const char *);
	deepStarData **catalogue_stars = malloc(nb_catalogues * sizeof(deepStarData *));
	uint32_t *catalogue_nb_stars = malloc(nb_catalogues * sizeof(uint32_t));
	uint32_t total_nb_stars = 0;
	int retval = 0, catalogue = 0;
	int16_t mag_threshold = (int16_t)roundf(max_mag * 1000.f);

	siril_debug_print("looking for stars in local catalogues for trixel %4d\n", ID);
	for (; catalogue < nb_catalogues; catalogue++) {
		retval = read_trixelID_from_catalogue(com.pref.catalogue_paths[catalogue],
				ID,
				catalogue_stars + catalogue, catalogue_nb_stars + catalogue);
		if (retval)
			break;
		total_nb_stars += catalogue_nb_stars[catalogue];
		cat_debug_print("%d raw stars from catalogue %d\n", catalogue_nb_stars[catalogue], catalogue);

	}

	if (catalogue == nb_catalogues) {
		// aggregate and filter
		*stars = malloc(total_nb_stars * sizeof(deepStarData));
		if (!*stars) {
			PRINT_ALLOC_ERR;
			retval = 1;
		} else {
			uint32_t j = 0;
			for (catalogue = 0; catalogue < nb_catalogues; catalogue++) {
				deepStarData *cat_stars = catalogue_stars[catalogue];
				uint32_t cat_nb_stars = catalogue_nb_stars[catalogue];

				for (uint32_t i = 0; i < cat_nb_stars ; i++) {
					if (cat_stars[i].V > mag_threshold)
						continue;
					if (cat_stars[i].B >= 30000 && photometric)
						continue;
					deepStarData *sdd = &cat_stars[i];
					(*stars)[j] = *sdd;
					j++;
				}
				free(catalogue_stars[catalogue]);
			}

			*nb_stars = j;
			*stars = realloc(*stars, j * sizeof(deepStarData));
		}
	}

	free(catalogue_nb_stars);
	free(catalogue_stars);
	return retval;
}

void initialize_local_catalogues_paths() {
	int nb_catalogues = sizeof(default_catalogues_paths) / sizeof(const char *);
	int maxpath = get_pathmax();
	for (int catalogue = 0; catalogue < nb_catalogues; catalogue++) {
		if (com.pref.catalogue_paths[catalogue] &&
				com.pref.catalogue_paths[catalogue][0] != '\0')
			continue;
		char path[maxpath];
		strncpy(path, default_catalogues_paths[catalogue], maxpath - 1);
		expand_home_in_filename(path, maxpath);
		com.pref.catalogue_paths[catalogue] = g_strdup(path);
	}
}

gboolean local_catalogues_available() {
	return local_kstars_available() || local_gaia_available();
}

siril_cat_index get_local_catalogue_index() {
	if (local_gaia_available())
		return CAT_LOCAL_GAIA_ASTRO;
	if (local_kstars_available())
		return CAT_LOCAL_KSTARS;
	return CAT_UNDEF;
}

gboolean local_kstars_available() {
	int nb_catalogues = 4;
	for (int catalogue = 0; catalogue < nb_catalogues; catalogue++) {
		if (!is_readable_file(com.pref.catalogue_paths[catalogue]))
			return FALSE;
	}
	return TRUE;
}

gboolean local_gaia_available() {
	if (!is_readable_file(com.pref.catalogue_paths[4]))
		return FALSE;
	return TRUE;
}

/* This function is the main interface to collect a local catalogue
   It sends a conesearch around given center, within given radius and for stars below limit_mag
   Internally, it uses get_raw_stars_from_local_catalogues which returns a list of deepStarData stars
   It fills the siril_catalogue given in input
   Returns the number of stars fetched, -1 if none found, 0 if error
*/
int siril_catalog_get_stars_from_local_catalogues(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return 0;
	if (siril_cat->cat_index != CAT_LOCAL_KSTARS &&
			siril_cat->cat_index != CAT_LOCAL_GAIA_ASTRO &&
			siril_cat->cat_index != CAT_LOCAL_GAIA_XPSAMP &&
			siril_cat->cat_index != CAT_LOCAL_TRIX) {
		siril_debug_print("Local cat query - Should not happen\n");
		return 0;
	}
	deepStarData *stars = NULL;
	uint32_t nb_stars;

	double ra_mult = siril_catalog_ra_multiplier(siril_cat->cat_index);
	double dec_mult = siril_catalog_dec_multiplier(siril_cat->cat_index);

	if (siril_cat->cat_index == CAT_LOCAL_GAIA_XPSAMP) {
		SourceEntryXPsamp *stars = NULL;
		if (get_raw_stars_from_local_gaia_xpsampled_catalogue(siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius, siril_cat->limitmag, &stars, &nb_stars))
			return 0;
		siril_cat->nbitems = (int)nb_stars;
		siril_cat->nbincluded = (int)nb_stars;
		siril_cat->cat_items = calloc(siril_cat->nbitems, sizeof(cat_item));
		for (int i = 0; i < siril_cat->nbitems; i++) {
			siril_cat->cat_items[i].xp_sampled = malloc(XPSAMPLED_LEN * sizeof(double));
			siril_cat->cat_items[i].ra = (double)stars[i].ra_scaled * ra_mult;
			siril_cat->cat_items[i].dec = (double)stars[i].dec_scaled * dec_mult;
			siril_cat->cat_items[i].pmra = (double)stars[i].dra_scaled;
			siril_cat->cat_items[i].pmdec = (double)stars[i].ddec_scaled;
			siril_cat->cat_items[i].mag = (float)stars[i].mag_scaled * 0.001;
			float powexp = pow(10.f, stars[i].fexpo);
			for (int j = 0 ; j < XPSAMPLED_LEN ; j++) {
				float d = half_to_float(stars[i].flux[j]);
				siril_cat->cat_items[i].xp_sampled[j] = d / powexp;
			}
			siril_cat->cat_items[i].included = TRUE;
		}
		free(stars);
		if (!siril_cat->nbitems)
			return -1;
		return siril_cat->nbitems;
	}

	if (siril_cat->cat_index == CAT_LOCAL_GAIA_ASTRO && get_raw_stars_from_local_gaia_astro_catalogue(siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius, siril_cat->limitmag, siril_cat->phot, &stars, &nb_stars))
		return 0;

	if (siril_cat->cat_index == CAT_LOCAL_KSTARS && get_raw_stars_from_local_catalogues(siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius, siril_cat->limitmag,
				siril_cat->phot, &stars, &nb_stars))
		return 0;
	if (siril_cat->cat_index == CAT_LOCAL_TRIX && get_raw_stars_from_local_catalogues_byID(siril_cat->trixel, siril_cat->limitmag, siril_cat->phot, &stars, &nb_stars))
		return 0;

	// compute the trixel centroid position from its 3 vertices
	if (siril_cat->cat_index == CAT_LOCAL_TRIX) {
		// double ra1, dec1, ra2, dec2, ra3, dec3;
		// get_vertices_for_index(siril_cat->trixel, 3, &ra1, &dec1, &ra2, &dec2, &ra3, &dec3);
		double ra[3], dec[3], x = 0., y = 0., z = 0., N;
		get_vertices_for_index(siril_cat->trixel, 3, ra, dec, ra + 1, dec + 1, ra + 2, dec + 2);
		for (int i = 0; i < 3; i++) {
			double a = ra[i] * DEGTORAD;
			double d = dec[i] * DEGTORAD;
			x += cos(d) * cos(a);
			y += cos(d) * sin(a);
			z += sin(d);
		}
		x /= 3.; y /= 3.; z /= 3.;
		N = sqrt( x * x + y * y + z * z);
		x /= N; y /= N; z /= N;
		siril_cat->center_dec = asin(z) * RADTODEG;
		siril_cat->center_ra = atan2(y, x) * RADTODEG;
		if (siril_cat->center_ra < 0.)
			siril_cat->center_ra +=360.;
	}

	siril_cat->nbitems = (int)nb_stars;
	siril_cat->nbincluded = (int)nb_stars;
	siril_cat->cat_items = calloc(siril_cat->nbitems, sizeof(cat_item));
	// Gaia catalogs present RA in deg, others present it in hours
	double CAT_RA_TO_DEG = siril_cat->cat_index == CAT_LOCAL_GAIA_ASTRO ? ra_mult : ra_mult * 15.0;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		siril_cat->cat_items[i].ra = (double)stars[i].RA * CAT_RA_TO_DEG;
		siril_cat->cat_items[i].dec = (double)stars[i].Dec * dec_mult;
		siril_cat->cat_items[i].pmra = (double)stars[i].dRA;
		siril_cat->cat_items[i].pmdec = (double)stars[i].dDec;
		siril_cat->cat_items[i].mag = (float)stars[i].V * .001;
		if (siril_cat->cat_index == CAT_LOCAL_GAIA_ASTRO)
			siril_cat->cat_items[i].teff = (float)(uint16_t)stars[i].B;
		else
			siril_cat->cat_items[i].bmag = (float)stars[i].B * .001;
		siril_cat->cat_items[i].included = TRUE;
	}
	free(stars);
	if (!siril_cat->nbitems)
		return -1;
	return siril_cat->nbitems;
}

gpointer write_trixels(gpointer p) {
	int nb[NTRIXELS];
	double ra[NTRIXELS], dec[NTRIXELS];
	int table[NTRIXELS][NMAG];
	double rav[NTRIXELS][3];
	double decv[NTRIXELS][3];
	for (int i = 0; i < NTRIXELS; i++) {
		siril_catalogue *siril_cat = siril_catalog_new(CAT_LOCAL_TRIX);
		siril_cat->trixel = i;
		siril_cat->limitmag = 25;
		if (!siril_catalog_conesearch(siril_cat)) {
			siril_catalog_free(siril_cat);
			siril_add_idle(end_generic, NULL);
			return GINT_TO_POINTER(1);
		}
		sort_cat_items_by_mag(siril_cat);
		siril_log_message("trixel #%4d - nbstars: %8d - ra: %5.1f, dec: %+5.1f\n", i, siril_cat->nbitems, siril_cat->center_ra, siril_cat->center_dec);
		nb[i] = siril_cat->nbitems;
		ra[i] = siril_cat->center_ra;
		dec[i] = siril_cat->center_dec;
		int n = 0;
		for (int j = 0; j < NMAG; j++) {
			double m = (double)j;
			while (n < siril_cat->nbitems && siril_cat->cat_items[n].mag < m)
				n++;
			// siril_log_message("%4.1f:%7d\n", m, n);
			table[i][j] = n;
		}
		siril_catalog_free(siril_cat);
		get_vertices_for_index(i, 3, rav[i], decv[i],  rav[i] + 1,decv[i] + 1,  rav[i] + 2, decv[i] + 2);
	}

	GError *error = NULL;
	GFile *file = g_file_new_for_path("trixels.csv");
	if (g_file_test("trixels.csv", G_FILE_TEST_EXISTS)) { // file already exists, we removed it
		if (!g_file_delete(file, NULL, &error)) {
			siril_log_color_message(_("File %s cannot be deleted (%s), aborting\n"), "red", "trixels.csv", (error) ? error->message : "unknown error");
			g_clear_error(&error);
			g_object_unref(file);
			siril_add_idle(end_generic, NULL);
			return GINT_TO_POINTER(1);
		}
	}
	// and we open the stream to write to
	GOutputStream *output_stream = (GOutputStream*) g_file_create(file, G_FILE_CREATE_NONE, NULL, &error);
	if (!output_stream) {
		siril_log_color_message(_("Cannot create catalogue file %s (%s)\n"), "red", "trixels.csv", (error) ? error->message : "unknown error");
		g_clear_error(&error);
		g_object_unref(file);
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}
	gsize s;
	g_output_stream_printf(output_stream, &s, NULL, &error, "%s", "trixel,nb,ra,dec,ra1,dec1,ra2,dec2,ra3,dec3");
	for (int i = 0; i < NMAG; i++) {
		g_output_stream_printf(output_stream, &s, NULL, &error, ",%d", i);
	}
	g_output_stream_printf(output_stream, &s, NULL, &error, "\n");
	for (int i = 0; i < NTRIXELS; i++) {
		g_output_stream_printf(output_stream, &s, NULL, &error, "%d,%d,%.4f,%+.4f,", i, nb[i], ra[i], dec[i]);
		g_output_stream_printf(output_stream, &s, NULL, &error, "%.4f,%+.4f,%.4f,%+.4f,%.4f,%+.4f", rav[i][0], decv[i][0], rav[i][1], decv[i][1], rav[i][2], decv[i][2]);
		for (int j = 0; j < NMAG; j++)
			g_output_stream_printf(output_stream, &s, NULL, &error, ",%d", table[i][j]);
		g_output_stream_printf(output_stream, &s, NULL, &error, "\n");
	}
	g_object_unref(file);
	g_object_unref(output_stream);
	siril_add_idle(end_generic, NULL);
	siril_log_message(_("Trixels list saved to trixels.csv\n"));
	return GINT_TO_POINTER(0);
}

gpointer list_trixels(gpointer p) {
	if (!has_wcs(gfit))
		return GINT_TO_POINTER(1);
	int nb_trixels;
	int *trixels;
	double ra[4], dec[4];
	double rx = (double)gfit->rx;
	double ry = (double)gfit->ry;
	pix2wcs(gfit, 0., 0., ra, dec);
	pix2wcs(gfit, rx, 0., ra + 1, dec + 1);
	pix2wcs(gfit, rx, ry, ra + 2, dec + 2);
	pix2wcs(gfit, 0., ry, ra + 3, dec + 3);
	int status = get_htm_indices_around_rectangle(ra, dec, 3, &trixels, &nb_trixels);
	if (status) {
		siril_log_color_message(_("Failed\n"), "red");
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(0);
	}
	siril_log_message(_("Found %d trixels:\n"), nb_trixels);
	for (int i = 0; i < nb_trixels; i++) {
		siril_log_message("%d\n", trixels[i]);
	}
	free(trixels);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(0);
}
