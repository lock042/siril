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
 * See also ../subprojects/htmesh/README for complete information.
 */

#include "io/kstars/binfile.h"
#include "io/kstars/htmesh_wrapper.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "algos/photometry.h"
#include "algos/siril_wcs.h"
#include "catalogues.h"
#include <stdio.h>
#include <stdlib.h>

#define NOMAD_DAT "~/.local/share/kstars/USNO-NOMAD-1e8.dat"

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

typedef struct
{
    int32_t RA;
    int32_t Dec;
    int16_t dRA;
    int16_t dDec;
    int16_t B;
    int16_t V;
} deepStarData;

static struct catalogue_file *catalogue_read_header(FILE *f);
//static int read_trixel_of_target(double ra, double dec, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);
static int read_trixels_of_target(double ra, double dec, double radius, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);
static int read_trixel(int trixel, struct catalogue_file *cat, deepStarData **stars, uint32_t *nb_stars);

static void bswap_stardata(deepStarData *stardata)
{
    stardata->RA   = bswap_32(stardata->RA);
    stardata->Dec  = bswap_32(stardata->Dec);
    stardata->dRA  = bswap_16(stardata->dRA);
    stardata->dDec = bswap_16(stardata->dDec);
    stardata->B    = bswap_16(stardata->B);
    stardata->V    = bswap_16(stardata->V);
}

int get_stars_from_local_nomad(double ra, double dec, double radius, fits *fit, float max_mag, pcc_star **stars, int *nb_stars) {
	char path[1024];
	strcpy(path, NOMAD_DAT);
	expand_home_in_filename(path, 1024);
	FILE *f = fopen(path, "r");
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

	deepStarData *trixel_stars;
	uint32_t trixel_nb_stars;
	if (read_trixels_of_target(ra, dec, radius, cat, &trixel_stars, &trixel_nb_stars)) {
		free(cat);
		fclose(f);
		return 1;
	}
	fclose(f);

	// project to image
	*stars = malloc(sizeof(pcc_star) * trixel_nb_stars);
	if (!*stars) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	int j = 0;
	for (int i = 0; i < trixel_nb_stars; i++) {
		if (trixel_stars[i].B >= 30000)
			continue;
		if (trixel_stars[i].V >= (int16_t)roundf(max_mag * 1000.f))
			continue;
		// catalogue has RA in hours, hence the x15
		double ra = trixel_stars[i].RA / 1000000.0 * 15.0;
		double dec= trixel_stars[i].Dec / 100000.0;
		double Bmag = trixel_stars[i].B / 1000.0;
		double Vmag = trixel_stars[i].V / 1000.0;
		// TODO: manage dRA and dDec

		double x, y;
		wcs2pix(fit, ra, dec, &x, &y);

		(*stars)[j].x = x;
		(*stars)[j].y = fit->ry - y - 1;
		(*stars)[j].mag = Vmag;
		(*stars)[j].BV = Bmag - Vmag;
		j++;
		//siril_debug_print("star at %f,\t%f,\tV=%.2f B=%.2f\timage coords (%.1f, %.1f)\n", ra, dec, Vmag, Bmag, (*stars)[j].x, (*stars)[j].y);
	}
	*nb_stars = j;

	return 0;
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

	siril_debug_print("Opened catalogue '%s'\n", top.description);

	struct catalogue_file *cat = calloc(1, sizeof(struct catalogue_file));
	if (!cat) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	if (top.endian_id != 0x4B53) {
		siril_debug_print("Byteswapping required\n");
		cat->byteswap = 1;
	}
	siril_debug_print("Version number: %d\n", top.version_no);

	/* reading the fields */
	if (!fread(&cat->nfields, 2, 1, f)) {
		siril_debug_print("error reading fields\n");
		free(cat);
		return NULL;
	}
	if (cat->byteswap)
		cat->nfields = bswap_16(cat->nfields);
	siril_debug_print("%d fields reported\n", cat->nfields);
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
		displayDataElementDescription(&(cat->de[i]));
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
	/* 512 should be for level 3, 8192 for 6, but it's not quite the case (32768 for NOMAD level 6) */
	siril_debug_print("Number of trixels reported = %d (levels: %d)\n",
			cat->ntrixels, cat->HTM_Level);

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
		siril_debug_print("faint magnitude: %.2f\n", header.faint_mag / 1000.0);
	}

	siril_debug_print("read the trixel index table, header read complete\n");
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
	if (!get_htm_indices_around_target(ra, dec, radius, cat->HTM_Level, &trixels, &nb_trixels)) {
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
	}
	else {
		*stars = NULL;
		*nb_stars = 0;
	}
	return 0;
}

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
			fseek_retval = fseek(cat->f, index->offset - ((u_int32_t)1 << 31) + 1, SEEK_CUR);
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

	if (fread(trix_stars, sizeof(deepStarData), index->nrecs, cat->f) < index->nrecs) {
		siril_debug_print("error reading trixel data\n");
		free(trix_stars);
		return 1;
	}

	if (cat->byteswap) {
		for (uint32_t i = 0; i < index->nrecs; ++i) {
			bswap_stardata(trix_stars + i);
			// TODO: manage dRA and dDec
		}
	}

	return 0;
}
