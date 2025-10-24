/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/fits_keywords.h"
#include "algos/demosaicing.h"
#include "algos/extraction.h"


const char *filter_pattern[] = {
	"RGGB",
	"BGGR",
	"GBRG",
	"GRBG",

/* XTRANS */
	"GGRGGB" // ----> XTRANS_1
	"GGBGGR"
	"BRGRBG"
	"GGBGGR"
	"GGRGGB"
	"RBGBRG",

	"RBGBRG" // ----> XTRANS_2
	"GGRGGB"
	"GGBGGR"
	"BRGRBG"
	"GGBGGR"
	"GGRGGB",

	"GRGGBG"
	"BGBRGR"
	"GRGGBG"
	"GBGGRG"
	"RGRBGB"
	"GBGGRG",

	"GBGGRG"
	"RGRBGB"
	"GBGGRG"
	"GRGGBG"
	"BGBRGR"
	"GRGGBG"
};
const size_t num_filter_patterns = G_N_ELEMENTS(filter_pattern);

////////////////////////////////
// Functions to get pattern info
////////////////////////////////
static gboolean get_debayer_orientation(fits *fit, gboolean *forced, gboolean *header) {
	gboolean top_down = FALSE;
	if (header)
		*header = FALSE;
	if (forced)
		*forced = com.pref.debayer.orientation == ROW_ORDER_FORCE_BOTTOMUP || com.pref.debayer.orientation == ROW_ORDER_FORCE_TOPDOWN;
	if (com.pref.debayer.orientation == ROW_ORDER_FORCE_BOTTOMUP)
		return FALSE;
	if (com.pref.debayer.orientation == ROW_ORDER_FORCE_TOPDOWN)
		return TRUE;
	if (com.pref.debayer.orientation == ROW_ORDER_HEADER_TOPDOWN || com.pref.debayer.orientation == ROW_ORDER_HEADER_BOTTOMUP) {
		if (!g_strcmp0(fit->keywords.row_order, "TOP-DOWN")) {
			if (header)
				*header = TRUE;
			top_down = TRUE;
		} else if (!g_strcmp0(fit->keywords.row_order, "BOTTOM-UP")) {
			if (header)
				*header = TRUE;
			top_down = FALSE;
		} else {
			top_down = com.pref.debayer.orientation == ROW_ORDER_HEADER_TOPDOWN;
			siril_debug_print("No row_order keyword found, falling back to preference\n");
		}
	}
	return top_down;
}

static void adjust_Bayer_pattern_offset(sensor_pattern *pattern, int xbayeroff, int ybayeroff) {
	if (!pattern)
		return;
	if (xbayeroff % 2 != 0) {
		switch (*pattern) {
		case BAYER_FILTER_RGGB:
			*pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_BGGR:
			*pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_GBRG:
			*pattern = BAYER_FILTER_BGGR;
			break;
		case BAYER_FILTER_GRBG:
			*pattern = BAYER_FILTER_RGGB;
			break;
		default:
			return;
		}
	}
	if (ybayeroff % 2 != 0) {
		switch (*pattern) {
		case BAYER_FILTER_RGGB:
			*pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_BGGR:
			*pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_GBRG:
			*pattern = BAYER_FILTER_RGGB;
			break;
		case BAYER_FILTER_GRBG:
			*pattern = BAYER_FILTER_BGGR;
			break;
		default:
			return;
		}
	}
}

static int adjust_Bayer_pattern(fits *fit, sensor_pattern *pattern, gboolean flip, int xbayeroff, int ybayeroff) {
	if (!fit || !pattern) {
		siril_log_color_message(_("Invalid FITS file or pattern for debayering\n"), "red");
		return 1;
	}
	if (*pattern < BAYER_FILTER_MIN || *pattern > BAYER_FILTER_MAX) {
		siril_log_color_message(_("Invalid Bayer pattern for debayering: %d\n"), "red", *pattern);
		return 1;
	}
	unsigned int ry = (!fit->orig_ry) ? fit->ry : fit->orig_ry;
	adjust_Bayer_pattern_orientation(pattern, ry, flip);
	adjust_Bayer_pattern_offset(pattern, xbayeroff, ybayeroff);
	return 0;
}

/* convert the string-described X-Trans pattern into an int array with value corresponding to filter */
static int compile_XTrans_pattern(const char *bayer, unsigned int xtrans[6][6], gboolean flip, int xbayeroff, int ybayeroff, int flip_offset) {
	int i = 0;

	if (strlen(bayer) != 36) {
		siril_log_color_message(_("FITS header does not contain a proper XTRANS pattern, demosaicing cannot be done\n"), "red");
		return 1;
	}
	if (flip) {
		unsigned int orig[6][6];
		memcpy(orig, xtrans, 36 * sizeof(unsigned int));
		for (int i = 0; i < 6; i++) {
			int y = (5 - i + flip_offset) % 6;
			for (int j = 0; j < 6; j++) {
				xtrans[i][j] = orig[y][j];
			}
		}
	}
	for (int y = 0; y < 6; y++) {
		int yoff = (y + ybayeroff) % 6; // apply y offset
		for (int x = 0; x < 6; x++) {
			int xoff = (x + xbayeroff) % 6; // apply x offset
			switch (bayer[i]) {
			case 'R':
				xtrans[yoff][xoff] = 0;
				break;
			case 'G':
				xtrans[yoff][xoff] = 1;
				break;
			case 'B':
				xtrans[yoff][xoff] = 2;
				break;
			default:
				siril_log_color_message(_("Invalid character in X-Trans filter pattern: %c\n"), "red", bayer[i]);
				return 1;
			}
			i++;
		}
	}
	return 0;
}

void adjust_Bayer_pattern_orientation(sensor_pattern *pattern, unsigned int ry, gboolean flip) {
	if (!pattern || !flip || *pattern == BAYER_FILTER_NONE)
		return;
	// we perform the pattern flip accounting for offset in case 
	// the image does not have an even number of rows
	const char *bayer = filter_pattern[*pattern];
	char bayer_flipped[5];
	memset(bayer_flipped, 0, sizeof(bayer_flipped));
	int offset = ry % 2;
	for (int i = 0; i < 2; i++) {
		int y = (1 - i + offset) % 2;
		for (int j = 0; j < 2; j++) {
			int index_in = y * 2 + j;
			int index_out = i * 2 + j;
			bayer_flipped[index_out] = bayer[index_in];
		}
	}
	*pattern = get_cfa_pattern_index_from_string(bayer_flipped);
}

// gets the index in filter_pattern
sensor_pattern get_cfa_pattern_index_from_string(const char *bayer) {
	if (bayer[0] == '\0') // so we avoid to read and compare against the whole list
		return BAYER_FILTER_NONE;
	for (int i = 0; i < G_N_ELEMENTS(filter_pattern); i++) {
		if (g_ascii_strcasecmp(bayer, filter_pattern[i]) == 0) {
			return i;
		}
	}
	return BAYER_FILTER_NONE;
}

// Returns the sensor_pattern based on header info and/or preferences.
// For real bayer pattern, the returned pattern also accounts for orientation
// as no matter what flipping or offset we have, we always get back to one of the 4 patterns
// For Xtrans patterns, flips and offsets are handled later with get_compiled_pattern
static sensor_pattern get_bayer_pattern(fits *fit, gboolean force_debayer, gboolean verbose) {
	/* Get Bayer informations from header if available,
		according to settings otherwise
	*/
	if (!fit)
		return BAYER_FILTER_NONE;

	sensor_pattern tmp_pattern = com.pref.debayer.bayer_pattern;
	gboolean from_header = FALSE;
	if (com.pref.debayer.use_bayer_header) {
		sensor_pattern bayer = get_cfa_pattern_index_from_string(fit->keywords.bayer_pattern);
		if (bayer == BAYER_FILTER_NONE) {
			siril_debug_print("No Bayer pattern found in the header file.\n");
			if (!force_debayer)
				return BAYER_FILTER_NONE;
		} else {
			from_header = TRUE;
			tmp_pattern = bayer;
		}
	}
	gboolean forced = FALSE, header = FALSE;
	gboolean top_down = get_debayer_orientation(fit, &forced, &header);
	if (verbose && tmp_pattern >= BAYER_FILTER_MIN) {
		siril_log_color_message(_("Filter Pattern: %s %s, Orientation: %s %s %s\n"), "blue",
		filter_pattern[tmp_pattern],
		from_header ? _("from header") : _("from settings"),
		top_down ? _("top-down") : _("bottom-up"),
		header ? _("from header") : _("from settings"),
		forced ? _("forced") : "");
	}

	if (tmp_pattern >= BAYER_FILTER_MIN && tmp_pattern <= BAYER_FILTER_MAX) {
		int xbayeroff = 0, ybayeroff = 0;
		if (!com.pref.debayer.use_bayer_header) {
			xbayeroff = com.pref.debayer.xbayeroff;
			ybayeroff = com.pref.debayer.ybayeroff;
		} else {
			xbayeroff = (fit->keywords.bayer_xoffset == DEFAULT_INT_VALUE) ? 0: fit->keywords.bayer_xoffset;
			ybayeroff = (fit->keywords.bayer_yoffset == DEFAULT_INT_VALUE) ? 0: fit->keywords.bayer_yoffset;
		}
		// x_offset amd y_offset are non null only when we do a partial read
		xbayeroff += fit->x_offset;
		ybayeroff += fit->y_offset;
	 	if (adjust_Bayer_pattern(fit, &tmp_pattern, !top_down, xbayeroff, ybayeroff)) {
			return BAYER_FILTER_NONE;
		}
	}
	// For X-trans, we don't update the pattern here,
	// it will be handled in get_compiled_pattern
	// just in time when needed. This enables to handle bottom-up images
	// and images with a height not multiple of 6 (or both).
	return tmp_pattern;
}

// This function is used to get the validated CFA pattern for debayering:
// If fit->debayer_checked is TRUE, it means the CFA pattern has already been validated
// it only returns the pattern from the header (this is TRUE for SER).
// Otherwise, it calls get_bayer_pattern to get the pattern from the FITS file 
// and according to preferences. It is much more chatty in this case.
// force-debayer is passed only when debayer has not been checked
sensor_pattern get_validated_cfa_pattern(fits *fit, gboolean force_debayer, gboolean verbose) {
	if (!fit) {
		siril_debug_print("No FITS file provided to get_validated_cfa_pattern\n");
		return BAYER_FILTER_NONE;
	}
	sensor_pattern pattern = BAYER_FILTER_NONE;
	if (fit->debayer_checked) {// TRUE for SER images
		pattern = get_cfa_pattern_index_from_string(fit->keywords.bayer_pattern);
	} else
		pattern = get_bayer_pattern(fit, force_debayer, verbose);
	siril_debug_print("Pattern to debayer: %s (%d)\n", filter_pattern[pattern], fit->debayer_checked);
	return pattern;
}


/* from the header description of the color filter array, we create a BYTE mask that
 * contains filter values */
int get_compiled_pattern(fits *fit, BYTE pattern[36], int *pattern_size, gboolean verbose) {
	sensor_pattern idx = get_validated_cfa_pattern(fit, FALSE, verbose);

	if (idx >= BAYER_FILTER_MIN && idx <= BAYER_FILTER_MAX) {
		/* 2x2 Bayer matrix */
		const gchar *cfa_str = filter_pattern[idx];
		for (int i = 0; i < 4; i++)
			pattern[i] = (BYTE)(cfa_str[i] == 'R' ? RLAYER :
								cfa_str[i] == 'G' ? GLAYER :
								cfa_str[i] == 'B' ? BLAYER : 3);
		*pattern_size = 2;
		return 0;
	}
	else if (idx >= XTRANS_FILTER_1 && idx <= XTRANS_FILTER_4) {
		/* 6x6 X-Trans matrix */
		unsigned int xtrans[6][6];
		const gchar *xtrans_str = filter_pattern[idx];
		gboolean top_down = get_debayer_orientation(fit, NULL, NULL);
		int flipoffset = fit->ry % 6;
		if (flipoffset)
			siril_debug_print("Image with an X-Trans sensor doesn't have a height multiple of 6\n");
		int xbayeroff = fit->keywords.bayer_xoffset;
		int ybayeroff = fit->keywords.bayer_yoffset;
		compile_XTrans_pattern(xtrans_str, xtrans, !top_down, xbayeroff, ybayeroff, flipoffset);
		for (int i = 0; i < 36; i++)
			pattern[i] = (BYTE)(((unsigned int *)xtrans)[i]);
		*pattern_size = 6;
		return 0;
	}
	return 1;
}

/** Calculate the filter color from a pattern array(in BYTE) using the row, column and pattern size
 * For size == 2, it uses bitwise operations that perform faster
 * Other wise it uses modulo which is slower
 **/
int FC_array(int row, int col, BYTE* bpattern, int size) {
	if (size == 2) {
		return bpattern[(row & 1) << 1 | (col & 1)];
	} else {
		// Direct array access with modulo
		return bpattern[((row % size) * size + (col % size))];
	}
}

gboolean compare_compiled_pattern(BYTE *refpattern, BYTE *pattern, int pattern_size) {
	int length  = pattern_size * pattern_size;
	for (int i = 0; i < length; i++) {
		if (refpattern[i] != pattern[i])
			return FALSE;
	}
	return TRUE;
}

// tests fits bayer_pattern to see if the first character is one of RGB
gboolean fit_is_cfa(fits *fit) {
	if (!fit)
		return FALSE;
	const char t = fit->keywords.bayer_pattern[0];
	return fit->naxes[2] == 1 && (t == 'R' || t == 'G' || t == 'B');
}

/* From an area, get the area corresponding to the debayer data for all colors,
 * the dashed area below.
 * 0 1 2 3 4 5
 * - - - - - -
 * - - - - - -
 * - - G R - -
 * - - B G - -
 * - - - - - -
 * - - - - - -
 *
 * area is the requested area of an image (simplified as GRBG above)
 * debayer_area is the result of this function, the area with enough pixels to
 *	have a valid debayer
 * image_area is the size of the image, to avoid going out of bounds
 * debayer_offset_x and y are the offset that need to be applied to the debayer
 *	data to find the original area (between 0 and 3).
 */
void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y) {
	int right, bottom;	// temp debayer negative offsets

	/* left side */
	if (area->x & 1)
		*debayer_offset_x = 11;
	else
		*debayer_offset_x = 10;
	if (area->x - *debayer_offset_x < 0) {
		debayer_area->x = 0;
		*debayer_offset_x = area->x;
	} else {
		debayer_area->x = area->x - *debayer_offset_x;
	}

	/* right side */
	int xend = area->x + area->w - 1;
	if (xend & 1)
		right = 10;
	else
		right = 11;
	if (xend + right >= image_area->w) {
		right = image_area->w - xend - 1;
	}
	debayer_area->w = area->w + (area->x - debayer_area->x) + right;

	/* top */
	if (area->y & 1)
		*debayer_offset_y = 11;
	else
		*debayer_offset_y = 10;
	if (area->y - *debayer_offset_y < 0) {
		debayer_area->y = 0;
		*debayer_offset_y = area->y;
	} else {
		debayer_area->y = area->y - *debayer_offset_y;
	}

	/* bottom */
	int yend = area->y + area->h - 1;
	if (yend & 1)
		bottom = 10;
	else
		bottom = 11;
	if (yend + bottom >= image_area->h) {
		bottom = image_area->h - yend - 1;
	}
	debayer_area->h = area->h + (area->y - debayer_area->y) + bottom;

	assert(debayer_area->x < image_area->w);
	assert(debayer_area->y < image_area->h);
	assert(debayer_area->h > 2);
	assert(debayer_area->w > 2);
}

//////////////////////////////////
// Functions to perform debayering
//////////////////////////////////

static int debayer_ushort(fits *fit, interpolation_method interpolation, sensor_pattern pattern) {
	int width = fit->rx;
	int height = fit->ry;
	WORD *buf = fit->data;

	unsigned int xtrans[6][6];
	if (interpolation == XTRANS) {
		const gchar *xtrans_str = filter_pattern[pattern];
		gboolean top_down = get_debayer_orientation(fit, NULL, NULL);
		int flipoffset = fit->ry % 6;
		if (flipoffset)
			siril_debug_print("Image with an X-Trans sensor doesn't have a height multiple of 6\n");
		int xbayeroff = fit->keywords.bayer_xoffset;
		int ybayeroff = fit->keywords.bayer_yoffset;
		compile_XTrans_pattern(xtrans_str, xtrans, !top_down, xbayeroff, ybayeroff, flipoffset);
	}
	// use librtprocess debayer
	WORD *newbuf = debayer_buffer_new_ushort(buf, &width, &height, interpolation, pattern, xtrans, fit->bitpix);
	if (!newbuf)
		return 1;

	fit_debayer_buffer(fit, newbuf);
	/* we remove Bayer header because not needed now */
	clear_Bayer_information(fit);

	/* The data is no longer mono. It's almost certainly linear, but we will
	 * not assign a color profile to it at this stage, it is likely one of
	 * many subs to be sequence processed along with others and it is more
	 * efficient to ignore color managemet until the final stacked image is
	 * available, and the user can then assign a profile as they choose.
	 */
	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = NULL;
	color_manage(fit, FALSE);
	return 0;
}

static int debayer_float(fits* fit, interpolation_method interpolation, sensor_pattern pattern) {
	int width = fit->rx;
	int height = fit->ry;
	float *buf = fit->fdata;

	unsigned int xtrans[6][6];
	if (interpolation == XTRANS) {
		const gchar *xtrans_str = filter_pattern[pattern];
		gboolean top_down = get_debayer_orientation(fit, NULL, NULL);
		int flipoffset = fit->ry % 6;
		if (flipoffset)
			siril_debug_print("Image with an X-Trans sensor doesn't have a height multiple of 6\n");
		int xbayeroff = fit->keywords.bayer_xoffset;
		int ybayeroff = fit->keywords.bayer_yoffset;
		compile_XTrans_pattern(xtrans_str, xtrans, !top_down, xbayeroff, ybayeroff, flipoffset);
	}

	float *newbuf = debayer_buffer_new_float(buf, &width, &height, interpolation, pattern, xtrans);
	if (!newbuf)
		return 1;

	fit_debayer_buffer(fit, newbuf);
	/* we remove Bayer header because not needed now */
	clear_Bayer_information(fit);

	/* The data is no longer mono. It's almost certainly linear, but we will
	 * not assign a color profile to it at this stage, it is likely one of
	 * many subs to be sequence processed along with others and it is more
	 * efficient to ignore color managemet until the final stacked image is
	 * available, and the user can then assign a profile as they choose.
	 */
	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = NULL;
	color_manage(fit, FALSE);

	return 0;
}

int debayer(fits *fit, interpolation_method interpolation, sensor_pattern pattern) {

	if (fit->type == DATA_USHORT)
		return debayer_ushort(fit, interpolation, pattern);
	else if (fit->type == DATA_FLOAT)
		return debayer_float(fit, interpolation, pattern);
	else return -1;
}

// debayers the image if it's a FITS image and if debayer is activated globally
// or if the force argument is passed
int debayer_if_needed(image_type imagetype, fits *fit, gboolean force_debayer) {
	if ((imagetype != TYPEFITS && imagetype != TYPETIFF && imagetype != TYPEXISF && imagetype != TYPERAW) || (!com.pref.debayer.open_debayer && !force_debayer))
		return 0;

	if (fit->naxes[2] != 1) {
		siril_log_message(_("Cannot perform debayering on image with more than one channel\n"));
		return 0;
	}
	sensor_pattern pattern = get_validated_cfa_pattern(fit, TRUE, TRUE);
	if (pattern < BAYER_FILTER_MIN) {
		siril_log_color_message(_("No pattern found to debayer\n"), "salmon");
		return 1;
	}

	/* Get Bayer informations from header if available */
	interpolation_method algo = com.pref.debayer.bayer_inter;
	if (pattern >= XTRANS_FILTER_1) {
		algo = XTRANS;
	}

	int retval = debayer(fit, algo, pattern);
	if (retval) {
		siril_log_message(_("Cannot perform debayering\n"));
	}
	return retval;
}

//////////////////////////////////////////
// Functions related to merging CFA images
//////////////////////////////////////////

static int mergecfa_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	required = 8 * MB_per_image;
	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
				thread_limit = com.max_thread;
		limit = thread_limit;
		if (for_writer)
			return 1;
	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
#endif
	}
	return limit;
}

void update_bayer_pattern_information(fits *fit, sensor_pattern pattern) {
	switch (pattern) {
		case BAYER_FILTER_RGGB:;
			sprintf(fit->keywords.bayer_pattern, "RGGB");
			break;
		case BAYER_FILTER_BGGR:;
			sprintf(fit->keywords.bayer_pattern, "BGGR");
			break;
		case BAYER_FILTER_GBRG:;
			sprintf(fit->keywords.bayer_pattern, "GBRG");
			break;
		case BAYER_FILTER_GRBG:;
			sprintf(fit->keywords.bayer_pattern, "GRBG");
			break;
		case XTRANS_FILTER_1:;
			sprintf(fit->keywords.bayer_pattern, "GGRGGBGGBGGRBRGRBGGGBGGRGGRGGBRBGBRG");
			break;
		case XTRANS_FILTER_2:;
			sprintf(fit->keywords.bayer_pattern, "RBGBRGGGRGGBGGBGGRBRGRBGGGBGGRGGRGGB");
			break;
		case XTRANS_FILTER_3:;
			sprintf(fit->keywords.bayer_pattern, "GRGGBGBGBRGRGRGGBGGBGGRGRGRBGBGBGGRG");
			break;
		case XTRANS_FILTER_4:;
			sprintf(fit->keywords.bayer_pattern, "GBGGRGRGRBGBGBGGRGGRGGBGBGBRGRGRGGBG");
			break;
		default:;
			break;
	}
}

gint64 mergecfa_compute_size_hook(struct generic_seq_args *args, int nb_frames) {
	double ratio = 4.;
	double fullseqsize = seq_compute_size(args->seq, nb_frames, args->output_type);
	return (gint64)(fullseqsize * ratio);
}

int mergecfa_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	fits metadata = { 0 };
	if (seq_read_frame_metadata(args->seq, out_index, &metadata)) {
		siril_log_message(_("Could not load metadata\n"));
		return 1;
	}
	int retval = 0;
	struct merge_cfa_data *merge_cfa_args = (struct merge_cfa_data*) args->user;
	fits cfa1 = { 0 };
	fits cfa2 = { 0 };
	fits cfa3 = { 0 };
	fits *out = { 0 };
	retval = seq_read_frame(merge_cfa_args->seq1, out_index, &cfa1, args->seq->bitpix == 16, -1);
	if(retval != 0) {
		siril_log_message(_("Image %d: error opening CFA1 file\n"), args->seq->current);
		goto CLEANUP_MERGECFA;
	}
	retval = seq_read_frame(merge_cfa_args->seq2, out_index, &cfa2, args->seq->bitpix == 16, -1);
	if(retval != 0) {
		siril_log_message(_("Image %d: error opening CFA2 file\n"), args->seq->current);
		goto CLEANUP_MERGECFA;
	}
	retval = seq_read_frame(merge_cfa_args->seq3, out_index, &cfa3, args->seq->bitpix == 16, -1);
	if(retval != 0) {
		siril_log_message(_("Image %d: error opening CFA3 file\n"), args->seq->current);
		goto CLEANUP_MERGECFA;
	}
	out = merge_cfa(fit, &cfa1, &cfa2, &cfa3, merge_cfa_args->pattern);
	if (out != NULL) {
		clearfits(fit);
		copyfits(out, fit, (CP_ALLOC | CP_COPYA | CP_FORMAT), -1);
		clearfits(out);
		siril_free(out);
	}
	copy_fits_metadata(&metadata, fit);
	update_sampling_information(fit, 0.5f);
	update_bayer_pattern_information(fit, merge_cfa_args->pattern);
	clearfits(&metadata);

CLEANUP_MERGECFA:
	clearfits(&cfa1);
	clearfits(&cfa2);
	clearfits(&cfa3);

	return retval;
}

int mergecfa_finalize_hook(struct generic_seq_args *args) {
	struct merge_cfa_data *data = (struct merge_cfa_data *) args->user;
	int retval = seq_finalize_hook(args);
	if (data->seq0 != args->seq && !check_seq_is_comseq(data->seq0)) {
		free_sequence(data->seq0, TRUE);
		data->seq0 = NULL;
	}
	if (data->seq1 != args->seq && !check_seq_is_comseq(data->seq1)) {
		free_sequence(data->seq1, TRUE);
		data->seq1 = NULL;
	}
	if (data->seq2 != args->seq && !check_seq_is_comseq(data->seq2)) {
		free_sequence(data->seq2, TRUE);
		data->seq2 = NULL;
	}
	if (data->seq3 != args->seq && !check_seq_is_comseq(data->seq3)) {
		free_sequence(data->seq3, TRUE);
		data->seq3 = NULL;
	}
	siril_free(data);
	return retval;
}

void apply_mergecfa_to_sequence(struct merge_cfa_data *merge_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(merge_cfa_args->seq0);
	args->seq = merge_cfa_args->seq0;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = merge_cfa_args->seq0->selnum;
	args->compute_mem_limits_hook = mergecfa_compute_mem_limits;
	args->compute_size_hook = mergecfa_compute_size_hook;
	args->prepare_hook = seq_prepare_hook;
	args->image_hook = mergecfa_image_hook;
	args->finalize_hook = mergecfa_finalize_hook;
	args->description = _("Merge CFA");
	args->has_output = TRUE;
	args->new_seq_prefix = strdup(merge_cfa_args->seqEntryOut);
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = merge_cfa_args;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		if (!check_seq_is_comseq(merge_cfa_args->seq0))
			free_sequence(merge_cfa_args->seq0, TRUE);
		if (!check_seq_is_comseq(merge_cfa_args->seq1))
			free_sequence(merge_cfa_args->seq1, TRUE);
		if (!check_seq_is_comseq(merge_cfa_args->seq2))
			free_sequence(merge_cfa_args->seq2, TRUE);
		if (!check_seq_is_comseq(merge_cfa_args->seq3))
			free_sequence(merge_cfa_args->seq3, TRUE);
		siril_free(merge_cfa_args);
		free_generic_seq_args(args, TRUE);
	}
}

//
// Re-mosaic 4 subpattern images previously separated with split_cfa
// This routine is Bayer pattern agnostic: split_cfa generates files
// CFA0, CFA1, CFA2 and CFA3 which can be RGGB, BGGR etc. This routine
// will assemble them in the same pattern as long as the files are
// provided in the correct order CFA0 to CFA3.
//
// (returns 0 on success, -1 on failure)
//
fits* merge_cfa (fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3, sensor_pattern pattern) {
	fits *out = NULL;

	// Check input files are compatible
	gboolean x_compat = (cfa0->naxes[0] == cfa1->naxes[0] && cfa1->naxes[0] == cfa2->naxes[0] && cfa2->naxes[0] == cfa3->naxes[0]);
	gboolean y_compat = (cfa0->naxes[1] == cfa1->naxes[1] && cfa1->naxes[1] == cfa2->naxes[1] && cfa2->naxes[1] == cfa3->naxes[1]);
	gboolean c_compat = (cfa0->naxes[2] == cfa1->naxes[2] && cfa1->naxes[2] == cfa2->naxes[2] && cfa2->naxes[2] == cfa3->naxes[2] && cfa3->naxes[2] == 1);
	gboolean t_compat = (cfa0->type == cfa1->type && cfa1->type == cfa2->type && cfa2->type == cfa3->type);
	if (!(x_compat && y_compat && c_compat && t_compat)) {
		siril_log_color_message(_("Input files are incompatible (all must be mono with the same size and bit depth). Aborting...\n"), "red");
		if(!x_compat)
			siril_log_message(_("X dimensions incompatible\n"));
		if(!y_compat)
			siril_log_message(_("Y dimensions incompatible\n"));
		if(!c_compat)
			siril_log_message(_("Channels not all mono\n"));
		if(!t_compat)
			siril_log_message(_("Input files not all the same bit depth\n"));
		return NULL;
	}
	int datatype = cfa0->type;

	// Create output fits twice the width and height of the cfa fits files
	if (new_fit_image(&out, cfa0->rx << 1, cfa0->ry << 1, 1, datatype)) {
		siril_log_color_message(_("Error creating output image\n"), "red");
		return NULL;
	}

//	out->header = copy_header(cfa0);
	copy_fits_metadata(cfa0, out);
	if (out->header)
		fprintf(stdout, "header ok\n");
	for (size_t outx = 0 ; outx < out->rx; outx += 2) {
		for(size_t outy = 0 ; outy < out->ry ; outy += 2) {
			size_t cfax = outx >> 1;
			size_t cfay = outy >> 1;
			size_t indexcfa = cfax + cfay * (size_t) cfa0->rx;
			size_t indexout0 = outx + outy * out->rx;
			size_t indexout1 = (outx + 1) + outy * out->rx;
			size_t indexout2 = outx + (outy + 1) * out->rx;
			size_t indexout3 = (outx + 1) + (outy + 1) * out->rx;
			switch (datatype) {
				case DATA_FLOAT:
					out->fdata[indexout0] = cfa0->fdata[indexcfa];
					out->fdata[indexout1] = cfa1->fdata[indexcfa];
					out->fdata[indexout2] = cfa2->fdata[indexcfa];
					out->fdata[indexout3] = cfa3->fdata[indexcfa];
					break;
				case DATA_USHORT:
					out->data[indexout0] = cfa0->data[indexcfa];
					out->data[indexout1] = cfa1->data[indexcfa];
					out->data[indexout2] = cfa2->data[indexcfa];
					out->data[indexout3] = cfa3->data[indexcfa];
					break;
			}
		}
	}
	// Set Bayer pattern in FITS header
	switch (pattern) {
		case BAYER_FILTER_RGGB:
			strcpy(out->keywords.bayer_pattern, "RGGB");
			break;
		case BAYER_FILTER_BGGR:
			strcpy(out->keywords.bayer_pattern, "BGGR");
			break;
		case BAYER_FILTER_GRBG:
			strcpy(out->keywords.bayer_pattern, "GRBG");
			break;
		case BAYER_FILTER_GBRG:
			strcpy(out->keywords.bayer_pattern, "GBRG");
			break;
		default:
			break;
	}
	clearfits(cfa0);
	clearfits(cfa1);
	clearfits(cfa2);
	clearfits(cfa3);
	siril_debug_print("Merge CFA complete\n");
	return out;
}

//////////////////////////////////////////////
// Functions related to extracting CFA buffers
//////////////////////////////////////////////

WORD *extract_CFA_buffer_ushort(fits *fit, int layer, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, pattern, &pattern_size, FALSE))
		return NULL;

	// alloc buffer
	size_t npixels = fit->naxes[0] * fit->naxes[1];
	WORD *buf = siril_malloc(npixels * sizeof(WORD));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	WORD *input = fit->data;

	// and we copy the pixels top-down
	size_t i = 0, j = 0;
	int x, y, pattern_y = 0;
	for (y = 0; y < fit->ry; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		for (x = 0; x < fit->rx; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void * tmp = siril_realloc(buf, j * sizeof(WORD));
	if (!tmp) {
		PRINT_ALLOC_ERR;
		siril_free(buf);
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}

WORD *extract_CFA_buffer_area_ushort(fits *fit, int layer, rectangle *bounds, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, pattern, &pattern_size, FALSE))
		return NULL;
	siril_debug_print("CFA buffer extraction with area\n");

	// alloc buffer
	size_t npixels = bounds->w * bounds->h;
	WORD *buf = siril_malloc(npixels * sizeof(WORD));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	WORD *input = fit->data;

	// copy the pixels top-down
	size_t j = 0;
	int start_y = fit->ry - bounds->y - bounds-> h;
	int pattern_y = start_y % pattern_size;
	for (int y = start_y; y < fit->ry - bounds->y; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		size_t i = y * fit->rx + bounds->x;
		for (int x = bounds->x; x < bounds->x + bounds->w; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void *tmp = siril_realloc(buf, j * sizeof(WORD));
	if (!tmp) {
		siril_free(buf);
		PRINT_ALLOC_ERR;
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}

float *extract_CFA_buffer_float(fits *fit, int layer, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, pattern, &pattern_size, FALSE))
		return NULL;

	// alloc buffer
	size_t npixels = fit->naxes[0] * fit->naxes[1];
	float *buf = siril_malloc(npixels * sizeof(float));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	float *input = fit->fdata;

	// copy the pixels top-down
	size_t i = 0, j = 0;
	int x, y, pattern_y = 0;
	for (y = 0; y < fit->ry; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		for (x = 0; x < fit->rx; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void *tmp = siril_realloc(buf, j * sizeof(float));
	if (!tmp) {
		siril_free(buf);
		PRINT_ALLOC_ERR;
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}

float *extract_CFA_buffer_area_float(fits *fit, int layer, rectangle *bounds, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, pattern, &pattern_size, FALSE))
		return NULL;
	siril_debug_print("CFA buffer extraction with area\n");

	// alloc buffer
	size_t npixels = bounds->w * bounds->h;
	float *buf = siril_malloc(npixels * sizeof(float));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	float *input = fit->fdata;

	// copy the pixels top-down
	size_t j = 0;
	int start_y = fit->ry - bounds->y - bounds-> h;
	int pattern_y = start_y % pattern_size;
	for (int y = start_y; y < fit->ry - bounds->y; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		size_t i = y * fit->rx + bounds->x;
		for (int x = bounds->x; x < bounds->x + bounds->w; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void *tmp = siril_realloc(buf, j * sizeof(float));
	if (!tmp) {
		siril_free(buf);
		PRINT_ALLOC_ERR;
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}
