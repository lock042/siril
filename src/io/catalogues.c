/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/proto.h"
#include "algos/photometry.h"
#include "algos/siril_wcs.h"
#include "algos/astrometry_solver.h"
#include "registration/matching/degtorad.h"
#include "registration/matching/misc.h"
#include "catalogues.h"
#include <stdio.h>
#include <stdlib.h>

#define NAMEDSTARS_DAT "~/.local/share/kstars/namedstars.dat"
#define UNNAMEDSTARS_DAT "~/.local/share/kstars/unnamedstars.dat"
#define TYCHOSTARS_DAT "~/.local/share/kstars/deepstars.dat"
#define NOMAD_DAT "~/.local/share/kstars/USNO-NOMAD-1e8.dat"
const char *catalogues_paths[] = { NAMEDSTARS_DAT, UNNAMEDSTARS_DAT, TYCHOSTARS_DAT, NOMAD_DAT };

struct catalogue_index {
	uint32_t trixelID;
	uint32_t offset;
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

static struct catalogue_file *catalogue_read_header(FILE *f);
//static int read_trixel_of_target(double ra, double dec, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);
static int read_trixels_of_target(double ra, double dec, double radius, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);
static int read_trixel(int trixel, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);
static int update_coords_with_proper_motion(double *ra, double *dec, double dRA, double dDec, double jmillenia);

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
	siril_debug_print("reading data from catalogue %s\n", path);
	FILE *f = fopen(path, "rb");
	if (!f) {
		siril_log_message(_("Could not open local NOMAD catalogue\n"));
		return 1;
	}

	struct catalogue_file *cat = catalogue_read_header(f);
	if (!cat) {
		siril_log_message(_("Failed to read the local NOMAD catalogue\n"));
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
	return 0;
}

static int get_projected_stars_from_local_catalogue(const char *path, double ra, double dec, double radius, gboolean use_proper_motion, fits *fit, float max_mag, pcc_star **stars, int *nb_stars) {
	deepStarData *trixel_stars;
	uint32_t trixel_nb_stars;
	if (read_trixels_from_catalogue(path, ra, dec, radius, &trixel_stars, &trixel_nb_stars))
		return 1;

	// project to image
	*stars = malloc(sizeof(pcc_star) * trixel_nb_stars);
	if (!*stars) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	GDateTime *dt = g_date_time_new_now_utc();
	gdouble jd = date_time_to_Julian(dt);
	g_date_time_unref(dt);
	double J2000 = 2451545.0;
	double jmillenia = (jd - J2000) / 365250.;
	int16_t mag_threshold = (int16_t)roundf(max_mag * 1000.f);

	int j = 0;
	for (int i = 0; i < trixel_nb_stars; i++) {
		// filter out stars with no valid photometry
		if (trixel_stars[i].B >= 30000 || trixel_stars[i].V >= 30000)
			continue;
		if (trixel_stars[i].V >= mag_threshold)
			continue;
		// catalogue has RA in hours, hence the x15
		double ra = trixel_stars[i].RA * .000015;
		double dec = trixel_stars[i].Dec * .00001;
		double Bmag = trixel_stars[i].B * 0.001;
		double Vmag = trixel_stars[i].V * 0.001;
		if (use_proper_motion) {
			update_coords_with_proper_motion(&ra, &dec,
					/* pass dRA and dDec per millenia */
					trixel_stars[i].dRA * 0.001, trixel_stars[i].dDec * 0.001,
					jmillenia);
		}

		double x, y;
		if (wcs2pix(fit, ra, dec, &x, &y))
			continue;

		(*stars)[j].x = x;
		(*stars)[j].y = fit->ry - y - 1;
		(*stars)[j].mag = Vmag;
		(*stars)[j].BV = Bmag - Vmag;
		j++;
		//siril_debug_print("star at %f,\t%f,\tV=%.2f B=%.2f\timage coords (%.1f, %.1f)\tpm (%hd, %hd) mas/yr\n", ra, dec, Vmag, Bmag, (*stars)[j].x, (*stars)[j].y, trixel_stars[i].dRA, trixel_stars[i].dDec);
	}
	*nb_stars = j;

	return 0;
}

/* ra and dec in degrees, dRA and dDec in mas/millenia, jmillenia in thousand years
 * ra and dec are input and output
 * return 1 if coords were updated, 0 else */
static int update_coords_with_proper_motion(double *ra, double *dec, double dRA, double dDec, double jmillenia) {
	double pmms = dRA * dRA + dDec * dDec;
	if (pmms * jmillenia * jmillenia < .01)
		return 0;
	siril_debug_print("PM passed the threshold\n");
	double ra_rad = *ra * DEGTORAD, dec_rad = *dec * DEGTORAD;
	double cosRa = cos(ra_rad), sinRa = sin(ra_rad);
	double cosDec = cos(dec_rad), sinDec = sin(dec_rad);
	double scale = jmillenia * (M_PI / (180.0 * 3600.0));

	// Note: Below assumes that dRA is already pre-scaled by cos(delta), as it is for Hipparcos
	double net_pmRA = dRA * scale, net_pmDec = dDec * scale;
	// net_pm are now in radian

	double x0 = cosDec * cosRa, y0 = cosDec * sinRa, z0 = sinDec;
	double dX = - net_pmRA * sinRa - net_pmDec * sinDec * cosRa;
	double dY = net_pmRA * cosRa - net_pmDec * sinDec * sinRa;
	double dZ = net_pmDec * cosDec;
	double x = x0 + dX, y = y0 + dY, z = z0 + dZ;

	double alpha = atan2(y, x) * RADTODEG;
	double delta = atan2(z, sqrt(x * x + y * y)) * RADTODEG;
	alpha = alpha - 360.0 * floor(alpha / 360.0);

	if (fabs(*ra - alpha) > .0000277 || fabs(*dec - delta) > .0000277)
		siril_debug_print("star moved from more than a tenth of arcsec on an axis (%f, %f)\n", *ra, *dec);
	*ra = alpha;
	*dec = delta;
	return 1;
}

/* get stars with coordinates projected on image (pcc_star), only those that
 * have B-V photometric information, up to the max_mag magnitude.
 * fit must have WCS inforation */
int get_stars_from_local_catalogues(double ra, double dec, double radius, fits *fit, float max_mag, pcc_star **stars, int *nb_stars) {
	if (!has_wcs(fit))
		return 1;
	int nb_catalogues = sizeof(catalogues_paths) / sizeof(const char *);
	pcc_star **catalogue_stars = malloc(nb_catalogues * sizeof(pcc_star *));
	int *catalogue_nb_stars = malloc(nb_catalogues * sizeof(int));
	int total_nb_stars = 0, retval = 0, catalogue = 0;

	for (; catalogue < nb_catalogues; catalogue++) {
		char path[1024];
		strcpy(path, catalogues_paths[catalogue]);
		expand_home_in_filename(path, 1024);
		retval = get_projected_stars_from_local_catalogue(path, ra, dec, radius,
				// Tycho-2 proper motions seem to be garbage, disabling PM computation for it
				/* catalogue != 2 */ FALSE, fit, max_mag,
				catalogue_stars + catalogue, catalogue_nb_stars + catalogue);
		if (retval)
			break;
		total_nb_stars += catalogue_nb_stars[catalogue];
		siril_debug_print("%d stars from catalogue %d\n", catalogue_nb_stars[catalogue], catalogue);
	}

	if (!retval) {
		// aggregate
		*stars = malloc(total_nb_stars * sizeof(pcc_star));
		if (!*stars) {
			PRINT_ALLOC_ERR;
			retval = 1;
		} else {
			uint32_t offset = 0;
			*nb_stars = total_nb_stars;
			for (int catalogue = 0; catalogue < nb_catalogues; catalogue++) {
				memcpy(*stars + offset, catalogue_stars[catalogue],
						catalogue_nb_stars[catalogue] * sizeof(pcc_star));
				offset += catalogue_nb_stars[catalogue];
			}
		}
	}

	free(catalogue_nb_stars);
	free(catalogue_stars);
	return retval;
}


struct top_header {
	char description[124];
	int16_t endian_id;
	uint8_t version_no;
};

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
	if (cat->byteswap)
		cat->nfields = bswap_16(cat->nfields);
	//siril_debug_print("%d fields reported:\n", cat->nfields);
	cat->de = malloc(cat->nfields * sizeof(dataElement));
	if (!cat->de) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	for (int i = 0; i < cat->nfields; ++i) {
		if (!fread(&(cat->de[i]), sizeof(dataElement), 1, f)) {
			siril_debug_print("error reading field descriptor %d\n", i);
			free(cat);
			return NULL;
		}
		if (cat->byteswap)
			cat->de[i].scale = bswap_32(cat->de[i].scale);
		//displayDataElementDescription(&(cat->de[i]));
	}

	if (cat->nfields != 6 && cat->nfields != 11) {
		siril_debug_print("the number of field is not recognized as one of the main catalogues\n");
		free(cat);
		return NULL;
	}

	/* reading the trixel index table */
	if (!fread(&cat->ntrixels, 4, 1, f)) {
		siril_debug_print("error reading number of trixels\n");
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

	cat->index_offset = ftell(f);

	cat->indices = malloc(cat->ntrixels * sizeof(struct catalogue_index));
	if (!cat->indices) {
		PRINT_ALLOC_ERR;
		free(cat);
		return NULL;
	}

	if (fread(cat->indices, sizeof(struct catalogue_index), cat->ntrixels, f) < cat->ntrixels) {
		siril_debug_print("unexpected read failure in index table\n");
		free(cat->indices);
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

/*static int read_trixel_of_target(double ra, double dec, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars) {
	int trixelID = get_htm_index_for_coords(ra, dec, cat->HTM_Level);

	if (trixelID > cat->ntrixels) {
		siril_debug_print("invalid trixel number %u for %u available\n", trixelID, cat->ntrixels);
		return 1;
	}

	return read_trixel(trixelID, cat, stars, nb_stars);
}*/

/* radius is in degrees */
static int read_trixels_of_target(double ra, double dec, double radius, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars) {
	int *trixels;
	int nb_trixels;
	if (get_htm_indices_around_target(ra, dec, radius, cat->HTM_Level, &trixels, &nb_trixels)) {
		*stars = NULL;
		*nb_stars = 0;
		return 0;
	}

	siril_debug_print("trixel search found %d trixels\n", nb_trixels);
	deepStarData **stars_list;
	uint32_t *nb_stars_list;
	stars_list = malloc(nb_trixels * sizeof(deepStarData *));
	if (!stars_list) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	nb_stars_list = malloc(nb_trixels * sizeof(uint32_t *));
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
		siril_debug_print("trixel %d (%d) contained %u stars\n", i, trixels[i], nb_stars_list[i]);
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
	int fseek_retval;
	/* If offset > 2^31 - 1, do the fseek in two steps */
	if (index->offset > INT_MAX) {
		if ((fseek_retval = fseek(cat->f, INT_MAX, SEEK_SET)))
			fseek_retval = fseek(cat->f, index->offset - ((uint32_t)1 << 31) + 1, SEEK_CUR);
	}
	else fseek_retval = fseek(cat->f, index->offset, SEEK_SET);
        if (fseek_retval) {
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

gboolean local_catalogues_available() {
	/* static path for now */
	int nb_catalogues = sizeof(catalogues_paths) / sizeof(const char *);
	for (int catalogue = 0; catalogue < nb_catalogues; catalogue++) {
		char path[1024];
		strcpy(path, catalogues_paths[catalogue]);
		expand_home_in_filename(path, 1024);
		if (!is_readable_file(path))
			return FALSE;
	}
	return TRUE;
}

// in degrees
static double compute_coords_distance(double ra1, double dec1, double ra2, double dec2) {
	double ra1minra2 = ra1 - ra2;
	double ra_diff;
	if (ra1 > ra2) {
		double ra1o = ra1 - 360.0;
		if (ra1minra2 < fabs(ra1o - ra2))
			ra_diff = ra1minra2;
		else ra_diff = ra1o - ra2;
	} else {
		double ra2o = ra2 - 360.0;
		if (fabs(ra1minra2) < fabs(ra1 - ra2o))
			ra_diff = ra1minra2;
		else ra_diff = ra1 - ra2o;
	}

	double dec_diff = dec1 - dec2;
	return sqrt(ra_diff * ra_diff + dec_diff * dec_diff);
}

/* similar to get_stars_from_local_catalogues except it doesn't project the
 * stars and filters them based on the radius instead of image bounds */
static int get_raw_stars_from_local_catalogues(double target_ra, double target_dec, double radius,
		float max_mag, gboolean photometric, deepStarData_dist **stars, uint32_t *nb_stars) {
	int nb_catalogues = sizeof(catalogues_paths) / sizeof(const char *);
	deepStarData **catalogue_stars = malloc(nb_catalogues * sizeof(deepStarData *));
	uint32_t *catalogue_nb_stars = malloc(nb_catalogues * sizeof(uint32_t));
	uint32_t total_nb_stars = 0;
	int retval = 0, catalogue = 0;

	siril_debug_print("looking for stars in local catalogues for target %f, %f, radius %f, magnitude %.2f, photometric: %d\n", target_ra, target_dec, radius, max_mag, photometric);
	for (; catalogue < nb_catalogues; catalogue++) {
		if (catalogue == 3 && max_mag < 12.0) {
			siril_debug_print("not querying NOMAD for this limit magnitude\n");
			catalogue_stars[catalogue] = NULL;
			catalogue_nb_stars[catalogue] = 0;
			continue;
		}

		char path[1024];
		strcpy(path, catalogues_paths[catalogue]);
		expand_home_in_filename(path, 1024);
		// Tycho-2 proper motions seem to be garbage, disabling PM computation for it
		retval = read_trixels_from_catalogue(path, target_ra, target_dec, radius,
				catalogue_stars + catalogue, catalogue_nb_stars + catalogue);
		if (retval)
			break;
		total_nb_stars += catalogue_nb_stars[catalogue];
		siril_debug_print("%d raw stars from catalogue %d\n", catalogue_nb_stars[catalogue], catalogue);
	}

	if (catalogue == nb_catalogues) {
		// aggregate and filter
		*stars = malloc(total_nb_stars * sizeof(deepStarData_dist));
		if (!*stars) {
			PRINT_ALLOC_ERR;
			retval = 1;
		} else {
			uint32_t j = 0;
			int16_t mag_threshold = (int16_t)roundf(max_mag * 1000.f);
			for (int catalogue = 0; catalogue < nb_catalogues; catalogue++) {
				deepStarData *cat_stars  = catalogue_stars[catalogue];
				uint32_t cat_nb_stars = catalogue_nb_stars[catalogue];

				for (uint32_t i = 0; i < cat_nb_stars ; i++) {
					int16_t mag;
					if (cat_stars[i].B >= 30000) {
						if (photometric) continue;
						mag = cat_stars[i].V;
					}
					else mag = cat_stars[i].B;
					if (cat_stars[i].V >= 30000) {
						if (photometric) continue;
						mag = cat_stars[i].B;
					}
					else mag = cat_stars[i].V;
					if (mag >= 30000) continue;
					if (mag >= mag_threshold)
						continue;

					// catalogue has RA in hours, hence the x15
					double ra = cat_stars[i].RA * .000015;
					double dec = cat_stars[i].Dec * .00001;
					double dist = compute_coords_distance(ra, dec, target_ra, target_dec);
					if (dist > radius)
						continue;
					// a bit of fun, storing the distance in the struct
					deepStarData_dist *sdd = (deepStarData_dist*)&cat_stars[i];
					sdd->distance = (float)(dist * 60.0);
					(*stars)[j] = *sdd;
					j++;
				}
			}

			*nb_stars = j;
			*stars = realloc(*stars, j * sizeof(deepStarData_dist));
		}
	}

	free(catalogue_nb_stars);
	free(catalogue_stars);
	return retval;
}

/* project a list of stars for plate solving
 * see registration/matching/project_coords.c proc_star_file() for original
 * code, based on a downloaded text file catalog */
static int project_local_catalog(deepStarData_dist *stars, uint32_t nb_stars, double center_ra, double center_dec, GFile *file_out, int doASEC) {
	double xx, yy, xi, eta;

	GError *error = NULL;
	GOutputStream *output_stream = (GOutputStream *)g_file_append_to(file_out, G_FILE_CREATE_NONE, NULL, &error);
	if (!output_stream) {
		if (error != NULL) {
			siril_debug_print("project_local_catalog: can't open file %s for output. [%s]", g_file_peek_path(file_out), error->message);
			g_clear_error(&error);
			return 1;
		}
	}

	for (uint32_t i = 0; i < nb_stars; i++) {
		double ra = stars[i].RA * .000015;
		double dec = stars[i].Dec * .00001;
		double Bmag = stars[i].B * 0.001;
		double Vmag = stars[i].V * 0.001;
		double delta_ra;
		if (ra < 10 && center_ra > 350) {
			delta_ra = (ra + 360) - center_ra;
		} else if (ra > 350 && center_ra < 10) {
			delta_ra = (ra - 360) - center_ra;
		} else {
			delta_ra = ra - center_ra;
		}
		delta_ra *= DEGTORAD;

		/*
		 * let's transform from (delta_RA, delta_Dec) to (xi, eta),
		 */
		double dec_rad = dec * DEGTORAD;
		double cent_dec_rad = center_dec * DEGTORAD;
		xx = cos(dec_rad) * sin(delta_ra);
		yy = sin(cent_dec_rad) * sin(dec_rad)
				+ cos(cent_dec_rad) * cos(dec_rad) * cos(delta_ra);
		xi = (xx / yy);

		xx = cos(cent_dec_rad) * sin(dec_rad)
				- sin(cent_dec_rad) * cos(dec_rad) * cos(delta_ra);
		eta = (xx / yy);

		/* if desired, convert xi and eta from radians to arcsec */
		if (doASEC > 0) {
			xi *= RADtoASEC;
			eta *= RADtoASEC;
		}

		gchar line[1000];
		if (doASEC)
			sprintf(line, "%f %12.5f %12.5f %f %f\n", stars[i].distance, xi, eta, Vmag, Bmag);
		else 	sprintf(line, "%f %13.6e %13.6e %f %f\n", stars[i].distance, xi, eta, Vmag, Bmag);

		g_output_stream_write_all(output_stream, line, strlen(line), NULL, NULL, NULL);
	}

	return 0;
}

/* project a list of stars into /tmp/catalog.proj for plate solving (function managing the file)
 * copied from project_catalog() in algos/astrometry_solver.c
 * radius is in degrees */
gchar *get_and_project_local_catalog(SirilWorldCS *catalog_center, double radius, double max_mag,
		gboolean for_photometry) {
	GError *error = NULL;
	gchar *foutput = NULL;
	double center_ra = siril_world_cs_get_alpha(catalog_center);
	double center_dec = siril_world_cs_get_delta(catalog_center);

	siril_debug_print("using local catalogues for plate solving\n");
	GFile *fproj = g_file_new_build_filename(g_get_tmp_dir(), "catalog.proj", NULL);
	/* We want to remove the file if already exisit */
	if (!g_file_delete(fproj, NULL, &error)
			&& !g_error_matches(error, G_IO_ERROR, G_IO_ERROR_NOT_FOUND)) {
		// deletion failed for some reason other than the file not existing:
		// so report the error
		g_warning("Failed to delete %s: %s", g_file_peek_path(fproj),
				error->message);
	}

	deepStarData_dist *stars = NULL;
	uint32_t nb_stars;
	if (get_raw_stars_from_local_catalogues(center_ra, center_dec, radius, max_mag,
				for_photometry, &stars, &nb_stars))
		return NULL;
	siril_debug_print("got %u stars from local catalogues\n", nb_stars);
	if (project_local_catalog(stars, nb_stars, center_ra, center_dec, fproj, 1))
		return NULL;
	foutput = g_file_get_path(fproj);
	g_object_unref(fproj);

	return foutput;
}

#if 0
void filter_deepStars_for_photometry(deepStarData_dist *stars, uint32_t nb_stars) {
	uint32_t j = 0;
	for (uint32_t i = 0; i < nb_stars; i++) {
		if (stars[i].B >= 30000 || stars[i].V >= 30000)
			continue;
		if (i == j)
			continue;
		stars[j] = stars[i];
		j++;
	}
	// realloc?
}
#endif
