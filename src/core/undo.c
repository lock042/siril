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

/* TODO: is there an opportunity to make this more efficient so that
 * either an image-only undo state, a mask-only undo state or a "both"
 * undo state can be saved, to reduce the size of undo states in storage?
 * The "both" state would apply to geometry changing operations. */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "gui/siril_preview.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/annotation_catalogues.h"
#include "gui/flis_gui.h"
#include "core/undo.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"

#ifndef O_BINARY
#define O_BINARY 0
#endif

/* *filename must be freed */
static int undo_build_swapfile(fits *fit, char **filename) {
	gchar *nameBuff;
	char name[] = "siril_swp-XXXXXX";
	gchar *tmpdir;
	int fd;

	tmpdir = com.pref.swap_dir;
	nameBuff = g_build_filename(tmpdir, name, NULL);
	fd = g_mkstemp(nameBuff);
	if (fd < 1) {
		siril_log_message(_("File I/O Error: Unable to create swap file in %s: [%s]\n"),
				tmpdir, strerror(errno));
		g_free(nameBuff);
		return 1;
	}

	size_t size = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

	errno = 0;
	// Write some data to the temporary file
	if (fit->type == DATA_USHORT) {
		if (-1 == write(fd, fit->data, size * sizeof(WORD))) {
			siril_log_message(_("File I/O Error: Unable to write swap file in %s: [%s]\n"),
					tmpdir, strerror(errno));
			g_free(nameBuff);
			g_close(fd, NULL);
			return 1;
		}
	} else if (fit->type == DATA_FLOAT) {
		if (-1 == write(fd, fit->fdata, size * sizeof(float))) {
			siril_log_message(_("File I/O Error: Unable to write swap file in %s: [%s]\n"),
					tmpdir, strerror(errno));
			g_free(nameBuff);
			g_close(fd, NULL);
			return 1;
		}
	}
	*filename = nameBuff;
	g_close(fd, NULL);

	return 0;
}

/* *filename must be freed */
static int undo_build_mask_swapfile(fits *fit, char **filename) {
	gchar *nameBuff;
	char name[] = "siril_msk-XXXXXX";
	gchar *tmpdir;
	int fd;

	if (!fit->mask || !fit->mask->data) {
		*filename = NULL;
		return 0;
	}

	tmpdir = com.pref.swap_dir;
	nameBuff = g_build_filename(tmpdir, name, NULL);
	fd = g_mkstemp(nameBuff);
	if (fd < 1) {
		siril_log_message(_("File I/O Error: Unable to create mask swap file in %s: [%s]\n"),
				tmpdir, strerror(errno));
		g_free(nameBuff);
		return 1;
	}

	size_t n_pixels = fit->naxes[0] * fit->naxes[1];
	size_t elem_size;

	switch (fit->mask->bitpix) {
		case 8:
			elem_size = sizeof(uint8_t);
			break;
		case 16:
			elem_size = sizeof(uint16_t);
			break;
		case 32:
			elem_size = sizeof(float);
			break;
		default:
			siril_log_message(_("Error: Invalid mask bitpix value: %d\n"), fit->mask->bitpix);
			g_free(nameBuff);
			g_close(fd, NULL);
			return 1;
	}

	size_t size = n_pixels * elem_size;

	errno = 0;
	if (-1 == write(fd, fit->mask->data, size)) {
		siril_log_message(_("File I/O Error: Unable to write mask swap file in %s: [%s]\n"),
				tmpdir, strerror(errno));
		g_free(nameBuff);
		g_close(fd, NULL);
		return 1;
	}

	*filename = nameBuff;
	g_close(fd, NULL);

	return 0;
}

/* Write a layermask_t to a new temporary swap file.
 * On success sets *filename, *out_w, *out_h, *out_bitpix and returns 0.
 * On failure *filename is NULL and returns 1. */
static int undo_build_lmask_swapfile(layermask_t *lm,
                                     gchar **filename,
                                     size_t *out_w, size_t *out_h,
                                     guint8 *out_bitpix) {
	*filename  = NULL;
	*out_w     = 0;
	*out_h     = 0;
	*out_bitpix = 0;

	if (!lm || !lm->data)
		return 0; /* no mask — caller stores NULL filename */

	size_t elem_size;
	switch (lm->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			siril_log_message(_("undo_build_lmask_swapfile: unsupported bitpix %d\n"),
			                  (int)lm->bitpix);
			return 1;
	}

	gchar *nameBuff = g_build_filename(com.pref.swap_dir, "siril_lmsk-XXXXXX", NULL);
	int fd = g_mkstemp(nameBuff);
	if (fd < 1) {
		siril_log_message(_("File I/O Error: Unable to create lmask swap file in %s: [%s]\n"),
		                  com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return 1;
	}

	size_t total = lm->w * lm->h * elem_size;
	errno = 0;
	if (write(fd, lm->data, total) == -1) {
		siril_log_message(_("File I/O Error: Unable to write lmask swap file: [%s]\n"),
		                  strerror(errno));
		g_close(fd, NULL);
		g_unlink(nameBuff);
		g_free(nameBuff);
		return 1;
	}
	g_close(fd, NULL);

	*filename   = nameBuff;
	*out_w      = lm->w;
	*out_h      = lm->h;
	*out_bitpix = lm->bitpix;
	return 0;
}

/* Free all heap resources owned by a historic_layer_entry_t.
 * Does NOT free the entry struct itself (it is typically part of an array). */
static void free_layer_entry(historic_layer_entry_t *e) {
	if (!e) return;
	if (e->filename) {
		if (g_unlink(e->filename))
			siril_debug_print("g_unlink() of multi-layer swap failed\n");
		g_free(e->filename);
		e->filename = NULL;
	}
	if (e->mask_filename) {
		if (g_unlink(e->mask_filename))
			siril_debug_print("g_unlink() of multi-layer mask swap failed\n");
		g_free(e->mask_filename);
		e->mask_filename = NULL;
	}
	if (e->lmask_filename) {
		if (g_unlink(e->lmask_filename))
			siril_debug_print("g_unlink() of multi-layer lmask swap failed\n");
		g_free(e->lmask_filename);
		e->lmask_filename = NULL;
	}
	if (e->wcslib) {
		wcsfree(e->wcslib);
		free(e->wcslib);
		e->wcslib = NULL;
	}
	if (e->icc_profile) {
		cmsCloseProfile(e->icc_profile);
		e->icc_profile = NULL;
	}
	g_free(e->layer_props);
	e->layer_props = NULL;
	memset(&e->wcsdata, 0, sizeof(wcs_info));
}

static int undo_remove_item(historic *histo, int index) {
	if (histo[index].filename) {
		if (g_unlink(histo[index].filename))
			siril_debug_print("g_unlink() failed\n");
		if (histo[index].icc_profile)
			cmsCloseProfile(histo[index].icc_profile);
		g_free(histo[index].filename);
		histo[index].filename = NULL;
		memset(&histo[index].wcsdata, 0, sizeof(wcs_info));
	}
	if (histo[index].mask_filename) {
		if (g_unlink(histo[index].mask_filename))
			siril_debug_print("g_unlink() of mask failed\n");
		g_free(histo[index].mask_filename);
		histo[index].mask_filename = NULL;
	}
	memset(histo[index].history, 0, FLEN_VALUE);
	histo[index].mask_bitpix = 0;
	histo[index].flis_layer_id = FLIS_UNDO_LAYER_NONE;
	histo[index].flis_position_x = 0;
	histo[index].flis_position_y = 0;
	g_free(histo[index].layer_props);
	histo[index].layer_props = NULL;
	if (histo[index].lmask_filename) {
		if (g_unlink(histo[index].lmask_filename))
			siril_debug_print("g_unlink() of lmask swap failed\n");
		g_free(histo[index].lmask_filename);
		histo[index].lmask_filename = NULL;
	}
	histo[index].lmask_layer_id  = FLIS_UNDO_LAYER_NONE;
	histo[index].lmask_dest_layer_id = FLIS_UNDO_LAYER_NONE;
	histo[index].lmask_w      = 0;
	histo[index].lmask_h      = 0;
	histo[index].lmask_bitpix = 0;
	histo[index].reorder_layer_a_id    = FLIS_UNDO_LAYER_NONE;
	histo[index].reorder_layer_a_order = 0;
	histo[index].reorder_layer_b_id    = FLIS_UNDO_LAYER_NONE;
	histo[index].reorder_layer_b_order = 0;
	histo[index].pmask_only  = FALSE;
	histo[index].full_layer  = FALSE;
	/* Free compound multi-layer entries */
	if (histo[index].multi_entries) {
		for (guint k = 0; k < histo[index].n_multi_entries; k++)
			free_layer_entry(&histo[index].multi_entries[k]);
		g_free(histo[index].multi_entries);
		histo[index].multi_entries   = NULL;
		histo[index].n_multi_entries = 0;
	}
	return 0;
}

static void undo_add_item(fits *fit, char *filename, char *mask_filename,
                          const char *histo, gint flis_layer_id,
                          gint flis_pre_pos_x, gint flis_pre_pos_y) {

	if (!com.history) {
		com.hist_size = HISTORY_SIZE;
		com.history = calloc(com.hist_size, sizeof(historic));
		com.hist_current = 0;
		com.hist_display = 0;
	}
	/* when undo, we remove all further items being after */
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}
	int status = -1;
	com.history[com.hist_current].filename = filename;
	com.history[com.hist_current].mask_filename = mask_filename;
	com.history[com.hist_current].mask_bitpix = (fit->mask && fit->mask->data) ? fit->mask->bitpix : 0;
	com.history[com.hist_current].rx = fit->rx;
	com.history[com.hist_current].ry = fit->ry;
	com.history[com.hist_current].nchans = fit->naxes[2];
	com.history[com.hist_current].type = fit->type;
	com.history[com.hist_current].wcsdata = fit->keywords.wcsdata;
	com.history[com.hist_current].wcslib = wcs_deepcopy(fit->keywords.wcslib, &status);
	if (status)
		siril_debug_print("could not copy wcslib struct\n");
	com.history[com.hist_current].focal_length = fit->keywords.focal_length;
	com.history[com.hist_current].icc_profile = copyICCProfile(fit->icc_profile);
	com.history[com.hist_current].flis_layer_id  = flis_layer_id;
	/* Position was snapshotted in undo_save_state before the operation ran */
	com.history[com.hist_current].flis_position_x = flis_pre_pos_x;
	com.history[com.hist_current].flis_position_y = flis_pre_pos_y;
	snprintf(com.history[com.hist_current].history, FLEN_VALUE, "%s", histo);

	if (com.hist_current == com.hist_size - 1) {
		/* we must shift all elements except 0 that must always match with the original file
		 * 0  1  2  3  4  5  6  7  8  9 10 become
		 * 0  2  3  4  5  6  7  8  9 10 11 and
		 * 0  3  4  5  6  7  8  9 10 11 12 and so on
		 */
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
				(com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}
	com.hist_current++;
	com.hist_display = com.hist_current;
}

/* Build a historic_layer_entry_t snapshot from a live FLIS layer.
 * Writes pixel and processing-mask swap files.
 * Returns 0 on success; on failure all already-written files are cleaned up
 * and 1 is returned. */
static int build_layer_entry(flis_layer_t *lay, historic_layer_entry_t *e) {
	memset(e, 0, sizeof(*e));
	e->flis_layer_id = lay->item_id;
	e->position_x    = lay->position_x;
	e->position_y    = lay->position_y;

	fits *fit = lay->fit;
	if (!fit) return 0;

	e->rx     = fit->rx;
	e->ry     = fit->ry;
	e->nchans = (int)fit->naxes[2];
	e->type   = fit->type;

	if (undo_build_swapfile(fit, &e->filename))
		return 1;

	if (undo_build_mask_swapfile(fit, &e->mask_filename)) {
		if (e->filename) { g_unlink(e->filename); g_free(e->filename); e->filename = NULL; }
		return 1;
	}
	e->mask_bitpix = (fit->mask && fit->mask->data) ? fit->mask->bitpix : 0;

	/* Layer mask snapshot */
	if (undo_build_lmask_swapfile(lay->lmask,
	                              &e->lmask_filename,
	                              &e->lmask_w, &e->lmask_h, &e->lmask_bitpix)) {
		if (e->mask_filename) { g_unlink(e->mask_filename); g_free(e->mask_filename); e->mask_filename = NULL; }
		if (e->filename) { g_unlink(e->filename); g_free(e->filename); e->filename = NULL; }
		return 1;
	}

	e->wcsdata = fit->keywords.wcsdata;
	if (fit->keywords.wcslib) {
		int status = -1;
		e->wcslib = wcs_deepcopy(fit->keywords.wcslib, &status);
		if (status)
			siril_debug_print("build_layer_entry: could not copy wcslib\n");
	}
	e->focal_length = fit->keywords.focal_length;
	e->icc_profile  = copyICCProfile(fit->icc_profile);

	e->layer_props = g_new(flis_layer_props_t, 1);
	e->layer_props->blend_mode   = lay->blend_mode;
	e->layer_props->opacity      = lay->opacity;
	e->layer_props->visible      = lay->visible;
	e->layer_props->locked       = lay->locked;
	e->layer_props->has_tint     = lay->has_tint;
	e->layer_props->tint         = lay->layer_tint;
	e->layer_props->lmask_active = lay->lmask_active;
	e->layer_props->position_x   = lay->position_x;
	e->layer_props->position_y   = lay->position_y;
	g_strlcpy(e->layer_props->name,
	          lay->layer_name ? lay->layer_name : "",
	          sizeof(e->layer_props->name));
	return 0;
}

static int undo_get_data_ushort(fits *fit, historic *hist) {
	int fd;

	if ((fd = g_open(hist->filename, O_RDONLY | O_BINARY, 0)) == -1) {
		printf("Error opening swap file : %s\n", hist->filename);
		return 1;
	}

	errno = 0;
	fit->rx = fit->naxes[0] = hist->rx;
	fit->ry = fit->naxes[1] = hist->ry;

	size_t n = fit->naxes[0] * fit->naxes[1];
	size_t size = n * fit->naxes[2] * sizeof(WORD);
	WORD *buf = calloc(1, size);
	// read the data from temporary file
	if ((read(fd, buf, size)) < size) {
		printf("Undo Read of [%s], failed with error [%s]\n", hist->filename, strerror(errno));
		free(buf);
		g_close(fd, NULL);
		return 1;
	}
	/* need to reallocate data as size may have changed */
	WORD *newdata = (WORD*) realloc(fit->data, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		free(newdata);
		free(buf);
		g_close(fd, NULL);
		return 1;
	}
	fit->data = newdata;
	memcpy(fit->data, buf, size);
	fit->pdata[RLAYER] = fit->data;
	if (fit->naxes[2] > 1) {
		fit->pdata[GLAYER] = fit->data + n;
		fit->pdata[BLAYER] = fit->data + n * 2;
	} else {
		fit->pdata[GLAYER] = fit->pdata[BLAYER] = fit->pdata[RLAYER];
	}
	memcpy(&fit->keywords.wcsdata, &hist->wcsdata, sizeof(wcs_info));
	if (hist->wcslib) {
		int status = -1;
		fit->keywords.wcslib = wcs_deepcopy(hist->wcslib, &status);
		if (status)
			siril_debug_print("could not copy wcslib struct\n");
	} else {
		free_wcs(fit);
		reset_wcsdata(fit);
	}
	fit->keywords.focal_length = hist->focal_length;

	full_stats_invalidation_from_fit(fit);
	free(buf);
	g_close(fd, NULL);
	return 0;
}

static int undo_get_data_float(fits *fit, historic *hist) {
	int fd;

	if ((fd = g_open(hist->filename, O_RDONLY | O_BINARY, 0)) == -1) {
		printf("Error opening swap file : %s\n", hist->filename);
		return 1;
	}

	errno = 0;
	fit->rx = fit->naxes[0] = hist->rx;
	fit->ry = fit->naxes[1] = hist->ry;

	size_t n = fit->naxes[0] * fit->naxes[1];
	size_t size = n * fit->naxes[2] * sizeof(float);
	float *buf = calloc(1, size);
	// read the data from temporary file
	if ((read(fd, buf, size) < size)) {
		printf("Undo Read of [%s], failed with error [%s]\n", hist->filename, strerror(errno));
		free(buf);
		g_close(fd, NULL);
		return 1;
	}
	/* need to reallocate data as size may have changed */
	float *newdata = (float*) realloc(fit->fdata, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		free(newdata);
		free(buf);
		g_close(fd, NULL);
		return 1;
	}
	fit->fdata = newdata;
	memcpy(fit->fdata, buf, size);
	fit->fpdata[RLAYER] = fit->fdata;
	if (fit->naxes[2] > 1) {
		fit->fpdata[GLAYER] = fit->fdata + n;
		fit->fpdata[BLAYER] = fit->fdata + n * 2;
	} else {
		fit->fpdata[GLAYER] = fit->fpdata[BLAYER] = fit->fpdata[RLAYER];
	}
	memcpy(&fit->keywords.wcsdata, &hist->wcsdata, sizeof(wcs_info));
	if (hist->wcslib) {
		int status = -1;
		fit->keywords.wcslib = wcs_deepcopy(hist->wcslib, &status);
		if (status)
			siril_debug_print("could not copy wcslib struct\n");
	} else {
		free_wcs(fit);
		reset_wcsdata(fit);
	}
	fit->keywords.focal_length = hist->focal_length;

	full_stats_invalidation_from_fit(fit);
	free(buf);
	g_close(fd, NULL);
	return 0;
}

static int undo_get_mask_data(fits *fit, historic *hist) {
	int fd;

	/* If there's no mask file saved, free any existing mask */
	if (!hist->mask_filename) {
		if (fit->mask) {
			if (fit->mask->data) {
				free(fit->mask->data);
				fit->mask->data = NULL;
			}
			free(fit->mask);
			fit->mask = NULL;
		}
		return 0;
	}

	/* Allocate or reallocate mask structure */
	if (!fit->mask) {
		fit->mask = calloc(1, sizeof(mask_t));
		if (!fit->mask) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	fit->mask->bitpix = hist->mask_bitpix;

	if ((fd = g_open(hist->mask_filename, O_RDONLY | O_BINARY, 0)) == -1) {
		printf("Error opening mask swap file : %s\n", hist->mask_filename);
		return 1;
	}

	size_t n_pixels = hist->rx * hist->ry;
	size_t elem_size;

	switch (hist->mask_bitpix) {
		case 8:
			elem_size = sizeof(uint8_t);
			break;
		case 16:
			elem_size = sizeof(uint16_t);
			break;
		case 32:
			elem_size = sizeof(float);
			break;
		default:
			siril_log_message(_("Error: Invalid mask bitpix value in history: %d\n"), hist->mask_bitpix);
			g_close(fd, NULL);
			return 1;
	}

	size_t size = n_pixels * elem_size;
	void *buf = calloc(1, size);
	if (!buf) {
		PRINT_ALLOC_ERR;
		g_close(fd, NULL);
		return 1;
	}

	errno = 0;
	if ((read(fd, buf, size)) < size) {
		printf("Undo Read of mask [%s], failed with error [%s]\n", hist->mask_filename, strerror(errno));
		free(buf);
		g_close(fd, NULL);
		return 1;
	}

	/* Reallocate mask data if needed */
	void *newdata = realloc(fit->mask->data, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		free(buf);
		g_close(fd, NULL);
		return 1;
	}

	fit->mask->data = newdata;
	memcpy(fit->mask->data, buf, size);

	free(buf);
	g_close(fd, NULL);
	return 0;
}

/* Restore pixel data, processing mask, WCS, position and props from a
 * historic_layer_entry_t into a live FLIS layer.  The entry must remain
 * valid for the duration of the call (its swap files are read but not freed). */
static int restore_layer_entry(flis_layer_t *lay, const historic_layer_entry_t *e) {
	/* Props-only entry: restore layer_props and nothing else */
	if (e->props_only) {
		if (e->layer_props) {
			const flis_layer_props_t *p = e->layer_props;
			lay->blend_mode   = p->blend_mode;
			lay->opacity      = p->opacity;
			lay->visible      = p->visible;
			lay->locked       = p->locked;
			lay->has_tint     = p->has_tint;
			lay->layer_tint   = p->tint;
			lay->lmask_active = p->lmask_active;
			lay->position_x   = p->position_x;
			lay->position_y   = p->position_y;
			g_free(lay->layer_name);
			lay->layer_name = g_strdup(p->name);
		}
		return 0;
	}

	fits *fit = lay->fit;
	if (!fit) return 0;

	/* Pixel data */
	if (e->filename) {
		/* Build a minimal historic that the data helpers understand.
		 * We only set the fields they actually read; everything else is 0. */
		historic tmp;
		memset(&tmp, 0, sizeof(tmp));
		tmp.filename     = e->filename;
		tmp.rx           = e->rx;
		tmp.ry           = e->ry;
		tmp.nchans       = e->nchans;
		tmp.type         = e->type;
		tmp.wcsdata      = e->wcsdata;
		tmp.wcslib       = e->wcslib;   /* read-only; not freed by helpers */
		tmp.focal_length = e->focal_length;

		if (fit->icc_profile)
			cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = copyICCProfile(e->icc_profile);
		color_manage(fit, (fit->icc_profile != NULL));
		fits_change_depth(fit, e->nchans);

		int retval;
		if (e->type == DATA_USHORT) {
			if (fit->type != DATA_USHORT) {
				size_t ndata = (size_t)fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
				fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
			}
			retval = undo_get_data_ushort(fit, &tmp);
		} else if (e->type == DATA_FLOAT) {
			if (fit->type != DATA_FLOAT) {
				size_t ndata = (size_t)fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
				fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
			}
			retval = undo_get_data_float(fit, &tmp);
		} else {
			retval = 1;
		}
		if (retval) return retval;
	}

	/* Processing mask — always call; undo_get_mask_data handles NULL (removes mask) */
	{
		historic tmp_mask;
		memset(&tmp_mask, 0, sizeof(tmp_mask));
		tmp_mask.mask_filename = e->mask_filename;
		tmp_mask.mask_bitpix   = e->mask_bitpix;
		tmp_mask.rx            = e->rx;
		tmp_mask.ry            = e->ry;
		if (undo_get_mask_data(fit, &tmp_mask))
			return 1;
	}

	/* Layer mask */
	if (!e->lmask_filename) {
		flis_layer_remove_lmask(lay);
	} else {
		size_t npix = e->lmask_w * e->lmask_h;
		size_t elem_size;
		switch (e->lmask_bitpix) {
			case 8:  elem_size = sizeof(uint8_t);  break;
			case 16: elem_size = sizeof(uint16_t); break;
			case 32: elem_size = sizeof(float);    break;
			default:
				siril_log_color_message(
					_("Undo: invalid lmask bitpix %d in compound entry\n"),
					"red", (int)e->lmask_bitpix);
				return 1;
		}
		int lfd = g_open(e->lmask_filename, O_RDONLY | O_BINARY, 0);
		if (lfd == -1) {
			siril_log_color_message(
				_("Undo: cannot open lmask swap %s\n"), "red", e->lmask_filename);
			return 1;
		}
		size_t total = npix * elem_size;
		void *buf = malloc(total);
		if (!buf) { PRINT_ALLOC_ERR; g_close(lfd, NULL); return 1; }
		if ((size_t)read(lfd, buf, total) < total) {
			siril_log_color_message(_("Undo: short read on lmask swap\n"), "red");
			free(buf); g_close(lfd, NULL); return 1;
		}
		g_close(lfd, NULL);
		layermask_t *lm = calloc(1, sizeof(layermask_t));
		if (!lm) { PRINT_ALLOC_ERR; free(buf); return 1; }
		lm->w      = e->lmask_w;
		lm->h      = e->lmask_h;
		lm->bitpix = e->lmask_bitpix;
		lm->data   = buf;
		if (flis_layer_set_lmask(lay, lm)) {
			layermask_free(lm);
			return 1;
		}
	}

	/* Layer offset */
	lay->position_x = e->position_x;
	lay->position_y = e->position_y;

	/* Layer properties */
	if (e->layer_props) {
		const flis_layer_props_t *p = e->layer_props;
		lay->blend_mode   = p->blend_mode;
		lay->opacity      = p->opacity;
		lay->visible      = p->visible;
		lay->locked       = p->locked;
		lay->has_tint     = p->has_tint;
		lay->layer_tint   = p->tint;
		lay->lmask_active = p->lmask_active;
		g_free(lay->layer_name);
		lay->layer_name = g_strdup(p->name);
	}

	full_stats_invalidation_from_fit(fit);
	return 0;
}

static int undo_get_data(fits *fit, historic *hist) {
	/* Compound multi-layer state: restore each layer in turn */
	if (hist->n_multi_entries > 0 && hist->multi_entries) {
		if (!is_current_image_flis()) return 1;
		for (guint k = 0; k < hist->n_multi_entries; k++) {
			const historic_layer_entry_t *e = &hist->multi_entries[k];
			flis_layer_t *lay = flis_layer_get_by_id(e->flis_layer_id);
			if (!lay) {
				siril_log_color_message(
					_("Undo: target layer (id %d) no longer exists — skipping\n"),
					"salmon", e->flis_layer_id);
				continue;
			}
			if (restore_layer_entry(lay, e)) {
				siril_log_color_message(
					_("Undo: failed to restore layer (id %d)\n"),
					"red", e->flis_layer_id);
				return 1;
			}
		}
		return 0;
	}

	/* Processing-mask-only state: restore fit->mask without touching pixels. */
	if (hist->pmask_only) {
		return undo_get_mask_data(fit, hist);
	}

	/* Full single-layer FLIS geometry state: pixels + pmask + lmask + props. */
	if (hist->full_layer) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *lay = flis_layer_get_by_id(hist->flis_layer_id);
		if (!lay) {
			siril_log_color_message(
				_("Undo: target layer (id %d) no longer exists\n"),
				"salmon", hist->flis_layer_id);
			return 1;
		}
		/* ICC profile */
		if (fit->icc_profile) cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = copyICCProfile(hist->icc_profile);
		color_manage(fit, (fit->icc_profile != NULL));
		fits_change_depth(fit, hist->nchans);
		/* Pixels */
		int retval;
		if (hist->type == DATA_USHORT) {
			if (fit->type != DATA_USHORT) {
				size_t ndata = (size_t)fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
				fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
			}
			retval = undo_get_data_ushort(fit, hist);
		} else if (hist->type == DATA_FLOAT) {
			if (fit->type != DATA_FLOAT) {
				size_t ndata = (size_t)fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
				fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
			}
			retval = undo_get_data_float(fit, hist);
		} else {
			retval = 1;
		}
		if (retval) return retval;
		/* Processing mask */
		retval = undo_get_mask_data(fit, hist);
		if (retval) return retval;
		/* Layer mask */
		if (!hist->lmask_filename) {
			flis_layer_remove_lmask(lay);
		} else {
			size_t npix = hist->lmask_w * hist->lmask_h;
			size_t elem_size;
			switch (hist->lmask_bitpix) {
				case 8:  elem_size = sizeof(uint8_t);  break;
				case 16: elem_size = sizeof(uint16_t); break;
				case 32: elem_size = sizeof(float);    break;
				default:
					siril_log_color_message(
						_("Undo: invalid lmask bitpix %d in full_layer state\n"),
						"salmon", (int)hist->lmask_bitpix);
					return 1;
			}
			int lfd = g_open(hist->lmask_filename, O_RDONLY | O_BINARY, 0);
			if (lfd == -1) {
				siril_log_color_message(
					_("Undo: cannot open lmask swap %s\n"), "red", hist->lmask_filename);
				return 1;
			}
			size_t total = npix * elem_size;
			void *buf = malloc(total);
			if (!buf) { PRINT_ALLOC_ERR; g_close(lfd, NULL); return 1; }
			if ((size_t)read(lfd, buf, total) < total) {
				siril_log_color_message(_("Undo: short read on lmask swap\n"), "red");
				free(buf); g_close(lfd, NULL); return 1;
			}
			g_close(lfd, NULL);
			layermask_t *lm = calloc(1, sizeof(layermask_t));
			if (!lm) { PRINT_ALLOC_ERR; free(buf); return 1; }
			lm->w      = hist->lmask_w;
			lm->h      = hist->lmask_h;
			lm->bitpix = hist->lmask_bitpix;
			lm->data   = buf;
			if (flis_layer_set_lmask(lay, lm)) { layermask_free(lm); return 1; }
		}
		/* Layer properties (blend mode, opacity, position, …) */
		if (hist->layer_props) {
			const flis_layer_props_t *p = hist->layer_props;
			lay->blend_mode   = p->blend_mode;
			lay->opacity      = p->opacity;
			lay->visible      = p->visible;
			lay->locked       = p->locked;
			lay->has_tint     = p->has_tint;
			lay->layer_tint   = p->tint;
			lay->lmask_active = p->lmask_active;
			lay->position_x   = p->position_x;
			lay->position_y   = p->position_y;
			g_free(lay->layer_name);
			lay->layer_name = g_strdup(p->name);
		}
		full_stats_invalidation_from_fit(fit);
		return 0;
	}

	/* Property-only undo state: no pixel swap file, just layer properties */
	if (!hist->filename && hist->layer_props) {
		if (hist->flis_layer_id == FLIS_UNDO_LAYER_NONE || !is_current_image_flis())
			return 1;
		flis_layer_t *layer = flis_layer_get_by_id(hist->flis_layer_id);
		if (!layer) {
			siril_log_color_message(
				_("Undo: target layer (id %d) no longer exists\n"),
				"salmon", hist->flis_layer_id);
			return 1;
		}
		const flis_layer_props_t *p = hist->layer_props;
		/* Apply properties directly, bypassing the mutating functions that
		 * would call flis_layer_touch_modified() and save a new undo state */
		layer->blend_mode   = p->blend_mode;
		layer->opacity      = p->opacity;
		layer->visible      = p->visible;
		layer->locked       = p->locked;
		layer->has_tint     = p->has_tint;
		layer->layer_tint   = p->tint;
		layer->lmask_active = p->lmask_active;
		layer->position_x   = p->position_x;
		layer->position_y   = p->position_y;
		g_free(layer->layer_name);
		layer->layer_name = g_strdup(p->name);
		return 0;
	}

	/* Atomic layer-mask move state: lmask_layer_id = source, lmask_dest_layer_id = dest.
	 * Restore by moving the mask in the reverse direction (dest → source). */
	if (!hist->filename && !hist->layer_props &&
	    hist->lmask_layer_id  != FLIS_UNDO_LAYER_NONE &&
	    hist->lmask_dest_layer_id != FLIS_UNDO_LAYER_NONE) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *src  = flis_layer_get_by_id(hist->lmask_layer_id);
		flis_layer_t *dest = flis_layer_get_by_id(hist->lmask_dest_layer_id);
		if (!src || !dest) {
			siril_log_color_message(
				_("Undo: layer for mask move no longer exists\n"), "salmon");
			return 1;
		}
		/* dest currently holds the mask; move it back to src */
		flis_layer_move_lmask(dest, src);
		return 0;
	}

	/* Layer reorder state: swap the saved layer_order values back and re-sort. */
	if (!hist->filename && !hist->layer_props &&
	    hist->reorder_layer_a_id != FLIS_UNDO_LAYER_NONE) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *a = flis_layer_get_by_id(hist->reorder_layer_a_id);
		flis_layer_t *b = flis_layer_get_by_id(hist->reorder_layer_b_id);
		if (!a || !b) {
			siril_log_color_message(
				_("Undo: layer for reorder no longer exists\n"), "salmon");
			return 1;
		}
		a->layer_order = hist->reorder_layer_a_order;
		b->layer_order = hist->reorder_layer_b_order;
		flis_sort_layer_stack();
		return 0;
	}

	/* Layer mask (lmask) add/remove state */
	if (!hist->filename && !hist->layer_props &&
	    hist->lmask_layer_id != FLIS_UNDO_LAYER_NONE) {
		if (!is_current_image_flis()) return 1;

		flis_layer_t *layer = flis_layer_get_by_id(hist->lmask_layer_id);
		if (!layer) {
			siril_log_color_message(
				_("Undo: target layer (id %d) no longer exists\n"),
				"salmon", hist->lmask_layer_id);
			return 1;
		}

		if (!hist->lmask_filename) {
			/* Saved state had no mask — remove current mask */
			flis_layer_remove_lmask(layer);
			return 0;
		}

		/* Saved state had a mask — restore it from the swap file */
		size_t npix = hist->lmask_w * hist->lmask_h;
		size_t elem_size;
		switch (hist->lmask_bitpix) {
			case 8:  elem_size = sizeof(uint8_t); break;
			case 16: elem_size = sizeof(uint16_t); break;
			case 32: elem_size = sizeof(float); break;
			default:
				siril_log_color_message(
					_("Undo: invalid lmask bitpix %d\n"),
					"salmon", (int)hist->lmask_bitpix);
				return 1;
		}

		int fd = g_open(hist->lmask_filename, O_RDONLY | O_BINARY, 0);
		if (fd == -1) {
			siril_log_color_message(
				_("Undo: cannot open lmask swap file %s\n"),
				"red", hist->lmask_filename);
			return 1;
		}

		size_t total = npix * elem_size;
		void *buf = malloc(total);
		if (!buf) { PRINT_ALLOC_ERR; g_close(fd, NULL); return 1; }

		if ((size_t)read(fd, buf, total) < total) {
			siril_log_color_message(
				_("Undo: short read on lmask swap file\n"), "red");
			free(buf); g_close(fd, NULL); return 1;
		}
		g_close(fd, NULL);

		layermask_t *lm = calloc(1, sizeof(layermask_t));
		if (!lm) { PRINT_ALLOC_ERR; free(buf); return 1; }
		lm->w      = hist->lmask_w;
		lm->h      = hist->lmask_h;
		lm->bitpix = hist->lmask_bitpix;
		lm->data   = buf;

		if (flis_layer_set_lmask(layer, lm)) {
			layermask_free(lm);
			return 1;
		}
		return 0;
	}

	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = copyICCProfile(hist->icc_profile);
	color_manage(fit, (fit->icc_profile != NULL));
	fits_change_depth(fit, hist->nchans);

	int retval = 0;

	if (hist->type == DATA_USHORT) {
		if (gfit->type != DATA_USHORT) {
			size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
			fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
			gui_function(set_precision_switch, NULL);
		}
		retval = undo_get_data_ushort(fit, hist);
	} else if (hist->type == DATA_FLOAT) {
		if (gfit->type != DATA_FLOAT) {
			size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
			fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
			gui_function(set_precision_switch, NULL);
		}
		retval = undo_get_data_float(fit, hist);
	} else {
		retval = 1;
	}

	/* Restore mask data regardless of image data success/failure */
	if (retval == 0) {
		retval = undo_get_mask_data(fit, hist);
	}

	/* Restore the layer offset if this state belongs to a FLIS layer.
	 * This is needed for geometry operations (e.g. crop) that modify
	 * both pixel dimensions and position_x/y. */
	if (retval == 0 && hist->flis_layer_id != FLIS_UNDO_LAYER_NONE
	    && is_current_image_flis()) {
		flis_layer_t *lay = flis_layer_get_by_id(hist->flis_layer_id);
		if (lay) {
			lay->position_x = hist->flis_position_x;
			lay->position_y = hist->flis_position_y;
		}
	}

	return retval;
}

gboolean is_undo_available() {
    return (com.history && com.hist_display > 0);
}

gboolean is_redo_available() {
    return (com.history && (com.hist_display < com.hist_current - 1));
}

int undo_save_state(fits *fit, const char *message, ...) {
	gchar *filename;
	gchar *mask_filename;
	va_list args;
	va_start(args, message);

	if (single_image_is_loaded()) {
		char histo[FLEN_VALUE] = { 0 };

		if (message != NULL) {
			vsnprintf(histo, FLEN_VALUE, message, args);
		}

		/* For a FLIS image, record which layer this state belongs to using
		 * the layer's stable item_id.  This allows undo to find and switch
		 * to the correct layer even after the stack has been reordered.
		 * For plain FITS images, flis_layer_id is set to FLIS_UNDO_LAYER_NONE.
		 *
		 * Snapshot position_x/y here — before the swap file is built and
		 * before any offset update (e.g. flis_update_layer_offset_after_crop)
		 * can have run on the caller's side. undo_add_item uses these values
		 * directly rather than re-reading from the live layer. */
		gint flis_layer_id = FLIS_UNDO_LAYER_NONE;
		gint flis_pre_pos_x = 0, flis_pre_pos_y = 0;
		if (is_current_image_flis()) {
			flis_layer_t *active = flis_active_layer();
			if (active) {
				flis_layer_id   = active->item_id;
				flis_pre_pos_x  = active->position_x;
				flis_pre_pos_y  = active->position_y;
			}
		}

		if (undo_build_swapfile(fit, &filename)) {
			va_end(args);
			return 1;
		}

		if (undo_build_mask_swapfile(fit, &mask_filename)) {
			g_free(filename);
			va_end(args);
			return 1;
		}

		undo_add_item(fit, filename, mask_filename, histo,
		              flis_layer_id, flis_pre_pos_x, flis_pre_pos_y);

		/* update menus */
		gui_function(update_MenuItem, NULL);
	}
	va_end(args);
	return 0;
}

/* Shared inner implementation — takes an already-built props snapshot and
 * adds it to the ring buffer.  Used by both undo_save_flis_layer_props()
 * (which builds the snapshot from the live layer) and the opacity drag-end
 * path in flis_gui.c (which builds the snapshot at drag-start). */
static int undo_push_flis_layer_props(gint item_id,
                                      const flis_layer_props_t *props,
                                      const char *histo) {
	flis_layer_props_t *copy = g_new(flis_layer_props_t, 1);
	*copy = *props;

	if (!com.history) {
		com.hist_size    = HISTORY_SIZE;
		com.history      = calloc(com.hist_size, sizeof(historic));
		com.hist_current = 0;
		com.hist_display = 0;
	}

	/* Discard any redo states above hist_display */
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}

	/* Handle ring-buffer wrap */
	if (com.hist_current == com.hist_size - 1) {
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
		        (com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}

	historic *h = &com.history[com.hist_current];
	memset(h, 0, sizeof(*h));
	h->filename              = NULL;
	h->layer_props           = copy;
	h->flis_layer_id         = item_id;
	h->lmask_layer_id        = FLIS_UNDO_LAYER_NONE;
	h->lmask_dest_layer_id   = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_a_id    = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_b_id    = FLIS_UNDO_LAYER_NONE;
	snprintf(h->history, FLEN_VALUE, "%s", histo ? histo : "");

	com.hist_current++;
	com.hist_display = com.hist_current;

	gui_function(update_MenuItem, NULL);
	return 0;
}

int undo_save_flis_layer_props(flis_layer_t *layer, const char *message, ...) {
	if (!layer || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	char histo[FLEN_VALUE] = { 0 };
	if (message) {
		va_list args;
		va_start(args, message);
		vsnprintf(histo, FLEN_VALUE, message, args);
		va_end(args);
	}

	flis_layer_props_t props;
	props.blend_mode   = layer->blend_mode;
	props.opacity      = layer->opacity;
	props.visible      = layer->visible;
	props.locked       = layer->locked;
	props.has_tint     = layer->has_tint;
	props.tint         = layer->layer_tint;
	props.lmask_active = layer->lmask_active;
	props.position_x   = layer->position_x;
	props.position_y   = layer->position_y;
	g_strlcpy(props.name,
	          layer->layer_name ? layer->layer_name : "", sizeof(props.name));

	return undo_push_flis_layer_props(layer->item_id, &props, histo);
}

int undo_save_flis_layer_props_snapshot(gint item_id,
                                        const flis_layer_props_t *props,
                                        const char *message) {
	if (!is_current_image_flis() || !single_image_is_loaded()) return 1;
	return undo_push_flis_layer_props(item_id, props, message);
}

/* undo_save_flis_lmask:
 *
 * Saves the current layer-mask state of @layer as an undo entry.
 * Call this BEFORE making any change to layer->lmask.
 *
 * If @layer currently has no mask, an "empty mask" marker is stored —
 * on undo, the current mask (whatever it is at that point) will be removed.
 * If @layer currently has a mask, its pixel data is written to a swap file —
 * on undo, the mask is restored from that file.
 *
 * This handles three operations:
 *   add mask:    call before setting the new mask — saves "no mask" marker
 *   remove mask: call before removing — saves current mask pixels
 *   move mask:   call on BOTH source and destination layers before the move
 */
int undo_save_flis_lmask(flis_layer_t *layer, const char *message) {
	if (!layer || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	if (!com.history) {
		com.hist_size    = HISTORY_SIZE;
		com.history      = calloc(com.hist_size, sizeof(historic));
		com.hist_current = 0;
		com.hist_display = 0;
	}

	/* Discard any redo states above hist_display */
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}

	/* Handle ring-buffer wrap */
	if (com.hist_current == com.hist_size - 1) {
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
		        (com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}

	historic *h = &com.history[com.hist_current];
	memset(h, 0, sizeof(*h));
	h->filename              = NULL;
	h->layer_props           = NULL;
	h->flis_layer_id         = FLIS_UNDO_LAYER_NONE;
	h->lmask_layer_id        = layer->item_id;
	h->lmask_dest_layer_id   = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_a_id    = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_b_id    = FLIS_UNDO_LAYER_NONE;
	snprintf(h->history, FLEN_VALUE, "%s", message ? message : "");

	if (!layer->lmask || !layer->lmask->data) {
		/* Layer currently has no mask — store the "no mask" marker */
		h->lmask_filename = NULL;
		h->lmask_w        = 0;
		h->lmask_h        = 0;
		h->lmask_bitpix   = 0;
	} else {
		/* Write the current mask pixels to a swap file */
		layermask_t *lm = layer->lmask;

		size_t elem_size;
		switch (lm->bitpix) {
			case 8:  elem_size = sizeof(uint8_t);  break;
			case 16: elem_size = sizeof(uint16_t); break;
			case 32: elem_size = sizeof(float);    break;
			default:
				siril_log_color_message(
					_("FLIS undo: unsupported lmask bitpix %d\n"),
					"red", (int)lm->bitpix);
				memset(h, 0, sizeof(*h));
				h->flis_layer_id  = FLIS_UNDO_LAYER_NONE;
				h->lmask_layer_id = FLIS_UNDO_LAYER_NONE;
				return 1;
		}

		gchar *nameBuff = g_build_filename(
		    com.pref.swap_dir, "siril_lmsk-XXXXXX", NULL);
		int fd = g_mkstemp(nameBuff);
		if (fd < 1) {
			siril_log_message(
				_("File I/O Error: Unable to create lmask swap file in %s: [%s]\n"),
				com.pref.swap_dir, strerror(errno));
			g_free(nameBuff);
			memset(h, 0, sizeof(*h));
			h->flis_layer_id  = FLIS_UNDO_LAYER_NONE;
			h->lmask_layer_id = FLIS_UNDO_LAYER_NONE;
			return 1;
		}

		size_t total = lm->w * lm->h * elem_size;
		errno = 0;
		if (write(fd, lm->data, total) == -1) {
			siril_log_message(
				_("File I/O Error: Unable to write lmask swap file: [%s]\n"),
				strerror(errno));
			g_close(fd, NULL);
			g_unlink(nameBuff);
			g_free(nameBuff);
			memset(h, 0, sizeof(*h));
			h->flis_layer_id  = FLIS_UNDO_LAYER_NONE;
			h->lmask_layer_id = FLIS_UNDO_LAYER_NONE;
			return 1;
		}
		g_close(fd, NULL);

		h->lmask_filename = nameBuff;
		h->lmask_w        = lm->w;
		h->lmask_h        = lm->h;
		h->lmask_bitpix   = lm->bitpix;
	}

	com.hist_current++;
	com.hist_display = com.hist_current;

	gui_function(update_MenuItem, NULL);
	return 0;
}

/* undo_save_flis_lmask_move:
 *
 * Saves an atomic undo state for a mask-move operation.  A single entry
 * records both the source and destination layer IDs so undo restores the
 * mask to its original layer in one step, with no broken intermediate state.
 * No swap file is written — the mask data stays alive in memory.
 *
 * Call this BEFORE flis_layer_move_lmask(source, dest).
 */
int undo_save_flis_lmask_move(flis_layer_t *source, flis_layer_t *dest,
                               const char *message) {
	if (!source || !dest || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	if (!com.history) {
		com.hist_size    = HISTORY_SIZE;
		com.history      = calloc(com.hist_size, sizeof(historic));
		com.hist_current = 0;
		com.hist_display = 0;
	}
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}
	if (com.hist_current == com.hist_size - 1) {
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
		        (com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}

	historic *h = &com.history[com.hist_current];
	memset(h, 0, sizeof(*h));
	h->flis_layer_id         = FLIS_UNDO_LAYER_NONE;
	h->lmask_layer_id        = source->item_id;
	h->lmask_dest_layer_id   = dest->item_id;
	h->lmask_filename        = NULL;   /* no swap file needed */
	h->reorder_layer_a_id    = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_b_id    = FLIS_UNDO_LAYER_NONE;
	snprintf(h->history, FLEN_VALUE, "%s", message ? message : "");

	com.hist_current++;
	com.hist_display = com.hist_current;
	gui_function(update_MenuItem, NULL);
	return 0;
}

/* undo_save_flis_layer_reorder:
 *
 * Saves an atomic undo state for a layer reorder (move-up / move-down).
 * Records the current layer_order of both layers before the swap so that
 * undo can restore both values atomically without a broken intermediate.
 *
 * Call this BEFORE flis_layer_move_up() or flis_layer_move_down(), passing
 * the layer being moved as layer_a and its swap partner as layer_b.
 */
int undo_save_flis_layer_reorder(flis_layer_t *layer_a, flis_layer_t *layer_b,
                                  const char *message) {
	if (!layer_a || !layer_b || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	if (!com.history) {
		com.hist_size    = HISTORY_SIZE;
		com.history      = calloc(com.hist_size, sizeof(historic));
		com.hist_current = 0;
		com.hist_display = 0;
	}
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}
	if (com.hist_current == com.hist_size - 1) {
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
		        (com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}

	historic *h = &com.history[com.hist_current];
	memset(h, 0, sizeof(*h));
	h->flis_layer_id         = FLIS_UNDO_LAYER_NONE;
	h->lmask_layer_id        = FLIS_UNDO_LAYER_NONE;
	h->lmask_dest_layer_id   = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_a_id    = layer_a->item_id;
	h->reorder_layer_a_order = layer_a->layer_order;
	h->reorder_layer_b_id    = layer_b->item_id;
	h->reorder_layer_b_order = layer_b->layer_order;
	snprintf(h->history, FLEN_VALUE, "%s", message ? message : "");

	com.hist_current++;
	com.hist_display = com.hist_current;
	gui_function(update_MenuItem, NULL);
	return 0;
}

/* Returns TRUE if the currently active layer (gfit) has either a processing
 * mask or a layer mask.  Used to decide whether to switch away from the mask
 * tab after a pure FLIS structural undo/redo. */
static gboolean flis_active_layer_has_any_mask(void) {
	if (gfit && gfit->mask) return TRUE;
	if (!is_current_image_flis() || !com.uniq) return FALSE;
	for (GSList *l = com.uniq->layers; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (lay && lay->fit == gfit)
			return (lay->lmask != NULL);
	}
	return FALSE;
}

/* Switch the active FLIS layer to match a historic state, so that
 * undo_get_data() restores pixel data into the correct layer's fits*.
 *
 * Returns TRUE if the switch succeeded or was not needed, FALSE if the
 * target layer no longer exists (e.g. it was deleted after the state was
 * saved) — the caller should skip this history entry in that case. */
static gboolean flis_undo_switch_layer(const historic *hist) {
	/* Compound multi-layer state: undo_get_data() iterates all entries by ID;
	 * there is no single target layer to switch gfit to. */
	if (hist->n_multi_entries > 0)
		return TRUE;

	/* Atomic lmask-move, lmask add/remove, and reorder states all address
	 * layers directly by ID in undo_get_data(); no gfit switch needed.
	 * pmask_only states DO need a layer switch (flis_layer_id is set). */
	if (!hist->filename && !hist->layer_props && !hist->pmask_only)
		return TRUE;

	if (hist->flis_layer_id == FLIS_UNDO_LAYER_NONE)
		return TRUE; /* plain FITS state: nothing to do */

	if (!is_current_image_flis() || !com.uniq) {
		siril_log_color_message(
			_("Undo: state was saved from a FLIS layer but no FLIS is loaded — skipping\n"),
			"salmon");
		return FALSE;
	}

	flis_layer_t *target = flis_layer_get_by_id(hist->flis_layer_id);
	if (!target) {
		siril_log_color_message(
			_("Undo: target layer (id %d) no longer exists — state cannot be restored\n"),
			"salmon", hist->flis_layer_id);
		return FALSE;
	}

	gint idx = flis_layer_get_index(target);
	if (idx < 0) {
		siril_log_color_message(
			_("Undo: could not locate target layer in stack — skipping\n"),
			"salmon");
		return FALSE;
	}

	/* Switch the active layer so gfit points at the target layer's fits* */
	uniq_set_active_layer(com.uniq, idx);
	return TRUE;
}

/* For FLIS images, init_right_tab() uses isrgb(gfit) to choose between
 * RGB_VPORT and RED_VPORT, but gfit is always the active layer (mono).
 * This wrapper uses com.uniq->chans, which we set to reflect the composite
 * colour model before calling close_tab(), so the correct tab is activated. */
static gboolean flis_init_right_tab(gpointer user_data) {
	(void)user_data;
	if (is_current_image_flis() && com.uniq)
		activate_tab(com.uniq->chans >= 3 ? RGB_VPORT : RED_VPORT);
	else
		init_right_tab(NULL);
	return FALSE;
}

int undo_display_data(int dir) {
	if (!com.history) {
		return 1;
	}
	switch (dir) {
	case UNDO:
		if (is_undo_available()) {
			// Avoid any issues with ROI or preview
			gboolean preview_was_active = is_preview_active();
			/* Writer lock: covers the ROI metadata reads (rx/ry/naxes[2]),
			 * undo_save_state (reads pixels), and undo_get_data (writes pixels).
			 * The entire save+restore must be atomic against the Python thread. */
			g_rw_lock_writer_lock(&gfit->rwlock);
			// Can't reactivate the ROI if the size has changed
			gboolean roi_was_active = (gui.roi.active && gfit->rx == com.history[com.hist_display - 1].rx
					&& gfit->ry == com.history[com.hist_display - 1].ry
					&& gfit->naxes[2] == com.history[com.hist_display - 1].nchans);
			rectangle roi_rect;
			memcpy(&roi_rect, &gui.roi.selection, sizeof(rectangle));
			siril_preview_hide();
			on_clear_roi();
			if (com.hist_current == com.hist_display) {
				/* We are at the tip of the history — save a matching "redo target"
				 * state before stepping back.  The saved type must mirror the state
				 * being undone so redo can restore it correctly. */
				historic *next = &com.history[com.hist_display - 1];
				gboolean is_flis = is_current_image_flis();

				if (next->n_multi_entries > 0 && next->multi_entries && is_flis) {
					/* Compound multi-layer state: rebuild layer list and snapshot.
					 * Use props-only variant if every entry is props-only. */
					gboolean all_props_only = TRUE;
					GSList *layers = NULL;
					for (guint k = 0; k < next->n_multi_entries; k++) {
						flis_layer_t *lay = flis_layer_get_by_id(
						        next->multi_entries[k].flis_layer_id);
						if (lay) layers = g_slist_prepend(layers, lay);
						if (!next->multi_entries[k].props_only)
							all_props_only = FALSE;
					}
					if (layers) {
						layers = g_slist_reverse(layers);
						if (all_props_only)
							undo_save_flis_multi_layer_props(layers, NULL);
						else
							undo_save_flis_multi_layer(layers, NULL);
						g_slist_free(layers);
					} else {
						undo_save_state(gfit, NULL);
					}
				} else if (!next->filename && is_flis) {
					if (next->layer_props) {
						/* Props-only state */
						flis_layer_t *active = flis_layer_get_by_id(next->flis_layer_id);
						if (active)
							undo_save_flis_layer_props(active, NULL);
						else
							undo_save_state(gfit, NULL);
					} else if (next->lmask_dest_layer_id != FLIS_UNDO_LAYER_NONE) {
						/* Atomic mask-move state: save reverse direction */
						flis_layer_t *src  = flis_layer_get_by_id(next->lmask_layer_id);
						flis_layer_t *dest = flis_layer_get_by_id(next->lmask_dest_layer_id);
						if (src && dest)
							undo_save_flis_lmask_move(dest, src, NULL);
						else
							undo_save_state(gfit, NULL);
					} else if (next->reorder_layer_a_id != FLIS_UNDO_LAYER_NONE) {
						/* Reorder state: save current positions of both layers */
						flis_layer_t *a = flis_layer_get_by_id(next->reorder_layer_a_id);
						flis_layer_t *b = flis_layer_get_by_id(next->reorder_layer_b_id);
						if (a && b)
							undo_save_flis_layer_reorder(a, b, NULL);
						else
							undo_save_state(gfit, NULL);
					} else if (next->lmask_layer_id != FLIS_UNDO_LAYER_NONE) {
						/* Lmask add/remove state */
						flis_layer_t *layer = flis_layer_get_by_id(next->lmask_layer_id);
						if (layer)
							undo_save_flis_lmask(layer, NULL);
						else
							undo_save_state(gfit, NULL);
					} else if (next->pmask_only) {
						/* Processing-mask-only state */
						undo_save_processing_mask(gfit, NULL);
					} else {
						undo_save_state(gfit, NULL);
					}
				} else if (next->full_layer && is_flis) {
					/* Full FLIS geometry state: snapshot current layer state */
					flis_layer_t *lay = flis_layer_get_by_id(next->flis_layer_id);
					if (lay) {
						flis_layer_props_t cur_props;
						cur_props.blend_mode   = lay->blend_mode;
						cur_props.opacity      = lay->opacity;
						cur_props.visible      = lay->visible;
						cur_props.locked       = lay->locked;
						cur_props.has_tint     = lay->has_tint;
						cur_props.tint         = lay->layer_tint;
						cur_props.lmask_active = lay->lmask_active;
						cur_props.position_x   = lay->position_x;
						cur_props.position_y   = lay->position_y;
						g_strlcpy(cur_props.name,
						          lay->layer_name ? lay->layer_name : "",
						          sizeof(cur_props.name));
						undo_save_flis_layer_full(gfit, lay, lay->lmask, &cur_props, NULL);
					} else {
						undo_save_state(gfit, NULL);
					}
				} else {
					undo_save_state(gfit, NULL);
				}
				com.hist_display--;
			}
			com.hist_display--;
			siril_log_message(_("Undo: %s\n"), com.history[com.hist_display].history);

			if (!flis_undo_switch_layer(&com.history[com.hist_display])) {
				/* Target layer gone; skip this state and step back one more */
				siril_log_color_message(
					_("Undo: skipped unreachable state\n"), "salmon");
				break;
			}

			undo_get_data(gfit, &com.history[com.hist_display]);

			/* Rebuild the FLIS composite since a layer's pixels changed */
			if (is_current_image_flis())
				flis_invalidate_composite();

			invalidate_gfit_histogram();
			invalidate_stats_from_fit(gfit);
			g_rw_lock_writer_unlock(&gfit->rwlock); // Finished with writer lock
			g_rw_lock_reader_lock(&gfit->rwlock);   // But still need reader lock
			update_gfit_histogram_if_needed();
			gui_function(close_tab, NULL); // These 2 lines account for possible change from mono to RGB
			g_rw_lock_reader_unlock(&gfit->rwlock);
			gui_function(update_MenuItem, NULL);
			lock_display_transform();
			if (gui.icc.proofing_transform)
				cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
			unlock_display_transform();
			refresh_annotations(TRUE);
			/* Update com.uniq->chans to reflect the composite colour model
			 * before close_tab() reads it — same correction as in
			 * open_single_image_from_gfit(). */
			if (is_current_image_flis() && com.uniq) {
				gboolean composite_rgb = FALSE;
				for (GSList *l = com.uniq->layers; l && !composite_rgb; l = l->next) {
					flis_layer_t *lay = (flis_layer_t *)l->data;
					if (!lay || !lay->fit) continue;
					if (lay->fit->naxes[2] >= 3) composite_rgb = TRUE;
					if (lay->has_tint)           composite_rgb = TRUE;
				}
				com.uniq->chans = composite_rgb ? 3 : 1;
			}
			/* Reinitialise channel tabs when pixel data changed.  For pure
			 * FLIS structural states (no swap file), only switch away from
			 * the mask tab if no mask remains on the active layer — e.g.
			 * when undoing mask creation.  If a mask still exists (e.g.
			 * undoing a blur), stay on the mask tab.  Always update mask
			 * tab visibility explicitly since close_tab() may be skipped. */
			if (com.history[com.hist_display].filename) {
				gui_function(close_tab, NULL);
				gui_function(flis_init_right_tab, NULL);
			} else {
				show_or_hide_mask_tab();
				if (!flis_active_layer_has_any_mask())
					gui_function(flis_init_right_tab, NULL);
			}
			flis_gui_update();
			/* redraw_mask_idle posts an idle — must be called outside any gfit lock */
			redraw_mask_idle(NULL);
			if (!com.pref.gui.mask_tints_vports) {// redraw() is called in redraw_mask_idle if this is TRUE
				g_rw_lock_reader_lock(&gfit->rwlock);
				notify_gfit_data_modified();
				g_rw_lock_reader_unlock(&gfit->rwlock);
				redraw(REMAP_ALL);
			}
			if (preview_was_active) {
				g_rw_lock_reader_lock(&gfit->rwlock);
				copy_gfit_to_backup();
				g_rw_lock_reader_unlock(&gfit->rwlock);
				siril_log_message(_("Following undo / redo with a preview active you may need "
						"to toggle the preview off and on again to reactivate the preview effect\n"));
				// TODO: To be perfect, we would need a register of preview functions
				// look up the correct one for the open dialog and re-apply the preview
			}
			if (roi_was_active) {
				memcpy(&com.selection, &roi_rect, sizeof(rectangle));
				on_set_roi();
			}
			g_rw_lock_reader_lock(&gfit->rwlock);
			update_fits_header(gfit);
			g_rw_lock_reader_unlock(&gfit->rwlock);
		}
		break;
	case REDO:
		if (is_redo_available()) {
			// Avoid any issues with ROI or preview
			gboolean preview_was_active = is_preview_active();
			/* Writer lock: covers the ROI metadata reads and undo_get_data (writes pixels). */
			g_rw_lock_writer_lock(&gfit->rwlock);
			// Can't reactivate the ROI if the size has changed
			gboolean roi_was_active = (gui.roi.active && gfit->rx == com.history[com.hist_display + 1].rx
					&& gfit->ry == com.history[com.hist_display + 1].ry
					&& gfit->naxes[2] == com.history[com.hist_display + 1].nchans);
			rectangle roi_rect;
			memcpy(&roi_rect, &gui.roi.selection, sizeof(rectangle));
			on_clear_roi();
			siril_preview_hide();
			siril_log_message(_("Redo: %s\n"), com.history[com.hist_display].history);
			com.hist_display++;

			if (!flis_undo_switch_layer(&com.history[com.hist_display])) {
				siril_log_color_message(
					_("Redo: skipped unreachable state\n"), "salmon");
				break;
			}

			undo_get_data(gfit, &com.history[com.hist_display]);

			/* Rebuild the FLIS composite since a layer's pixels changed */
			if (is_current_image_flis())
				flis_invalidate_composite();

			invalidate_gfit_histogram();
			invalidate_stats_from_fit(gfit);
			g_rw_lock_writer_unlock(&gfit->rwlock); // Finished with writer lock
			g_rw_lock_reader_lock(&gfit->rwlock);   // But still need reader lock
			update_gfit_histogram_if_needed();
			g_rw_lock_reader_unlock(&gfit->rwlock);
			gui_function(update_MenuItem, NULL);
			refresh_annotations(TRUE);
			lock_display_transform();
			if (gui.icc.proofing_transform)
				cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
			unlock_display_transform();
			/* Update com.uniq->chans to reflect the composite colour model
			 * before close_tab() reads it. */
			if (is_current_image_flis() && com.uniq) {
				gboolean composite_rgb = FALSE;
				for (GSList *l = com.uniq->layers; l && !composite_rgb; l = l->next) {
					flis_layer_t *lay = (flis_layer_t *)l->data;
					if (!lay || !lay->fit) continue;
					if (lay->fit->naxes[2] >= 3) composite_rgb = TRUE;
					if (lay->has_tint)           composite_rgb = TRUE;
				}
				com.uniq->chans = composite_rgb ? 3 : 1;
			}
			if (com.history[com.hist_display].filename) {
				gui_function(close_tab, NULL);
				gui_function(flis_init_right_tab, NULL);
			} else {
				show_or_hide_mask_tab();
				if (!flis_active_layer_has_any_mask())
					gui_function(flis_init_right_tab, NULL);
			}
			flis_gui_update();
			/* redraw_mask_idle posts an idle — must be called outside any gfit lock */
			redraw_mask_idle(NULL);
			if (!com.pref.gui.mask_tints_vports) { // redraw() is called in redraw_mask_idle if this is TRUE
				g_rw_lock_reader_lock(&gfit->rwlock);
				notify_gfit_data_modified();
				g_rw_lock_reader_unlock(&gfit->rwlock);
				redraw(REMAP_ALL);
			}
			if (preview_was_active) {
				g_rw_lock_reader_lock(&gfit->rwlock);
				copy_gfit_to_backup();
				g_rw_lock_reader_unlock(&gfit->rwlock);
			}
			if (roi_was_active) {
				memcpy(&gui.roi.selection, &roi_rect, sizeof(rectangle));
				on_set_roi();
			}
			g_rw_lock_reader_lock(&gfit->rwlock);
			update_fits_header(gfit);
			g_rw_lock_reader_unlock(&gfit->rwlock);
		}
		break;
	default:
		printf("ERROR\n");
		return -1;
	}
	return 0;
}

gboolean undo_in_thread(gpointer user_data) {
	undo_display_data(UNDO);
	return FALSE;
}

gboolean redo_in_thread(gpointer user_data) {
	undo_display_data(REDO);
	return FALSE;
}

int undo_flush() {
	if (!com.history) {
		return 1;
	}
	for (int i = 0; i < com.hist_current; i++) {
		undo_remove_item(com.history, i);
	}
	free(com.history);
	com.history = NULL;
	com.hist_current = 0;
	com.hist_display = 0;
	return 0;
}

/* flis_undo_purge_layer:
 *
 * Removes all undo/redo states that belong to the specified FLIS layer from
 * the ring buffer and frees their swap files.  Called before a layer is
 * permanently removed from the stack so that stale states referencing a
 * now-nonexistent item_id do not clutter the history or trigger warnings.
 *
 * States belonging to all other layers are preserved unchanged.
 *
 * When states are purged from the middle of the ring buffer the entries are
 * compacted: everything above the removed slot is shifted down by one so the
 * buffer remains contiguous.  hist_current and hist_display are adjusted to
 * match.
 *
 * This replaces the previous flis_undo_notify_structural_change() which
 * flushed the entire history.  Add, duplicate, and reorder operations do
 * not need to touch the history at all because:
 *   - Add / duplicate: the new layer has no prior states (new item_id).
 *   - Reorder: flis_layer_get_by_id() locates layers by ID not position,
 *     so reordering is transparent to the undo mechanism.
 */
void flis_undo_purge_layer(gint item_id) {
	if (!com.history || item_id == FLIS_UNDO_LAYER_NONE) return;

	int purged = 0;
	int i = 0;
	while (i < com.hist_current) {
		/* Check whether this compound entry references item_id */
		gboolean multi_match = FALSE;
		if (com.history[i].n_multi_entries > 0 && com.history[i].multi_entries) {
			for (guint k = 0; k < com.history[i].n_multi_entries; k++) {
				if (com.history[i].multi_entries[k].flis_layer_id == item_id) {
					multi_match = TRUE;
					break;
				}
			}
		}
		if (multi_match ||
		    com.history[i].flis_layer_id == item_id ||
		    com.history[i].lmask_layer_id == item_id ||
		    com.history[i].lmask_dest_layer_id == item_id ||
		    com.history[i].reorder_layer_a_id == item_id ||
		    com.history[i].reorder_layer_b_id == item_id) {
			undo_remove_item(com.history, i);
			/* Compact: shift everything above this slot down */
			if (i < com.hist_current - 1)
				memmove(&com.history[i], &com.history[i + 1],
				        (com.hist_current - i - 1) * sizeof(*com.history));
			/* Zero the vacated top slot */
			memset(&com.history[com.hist_current - 1], 0,
			       sizeof(*com.history));
			com.history[com.hist_current - 1].flis_layer_id =
			        FLIS_UNDO_LAYER_NONE;
			com.hist_current--;
			/* Clamp the display pointer in case it pointed into the
			 * vacated region */
			if (com.hist_display > com.hist_current)
				com.hist_display = com.hist_current;
			purged++;
			/* Do not advance i: the slot we just vacated now holds the
			 * next entry that needs checking */
		} else {
			i++;
		}
	}

	if (purged) {
		siril_debug_print("FLIS undo: purged %d state(s) for layer id %d\n",
		                  purged, item_id);
		gui_function(update_MenuItem, NULL);
	}
}

/* undo_save_flis_multi_layer:
 *
 * Saves an atomic compound undo state covering every layer in @layers.
 * For each layer one pixel swap file and one processing-mask swap file
 * (if the layer has a mask) are written.  A flis_layer_props_t snapshot
 * is also stored.  A single undo step reverts all layers at once.
 *
 * Call BEFORE the operation runs so the saved state holds pre-operation values.
 * Returns 0 on success; on failure all swap files already written are cleaned
 * up and no state is pushed.
 */
int undo_save_flis_multi_layer(GSList *layers, const char *message, ...) {
	if (!layers || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	guint n = g_slist_length(layers);
	if (n == 0) return 1;

	char histo[FLEN_VALUE] = { 0 };
	if (message) {
		va_list args;
		va_start(args, message);
		vsnprintf(histo, FLEN_VALUE, message, args);
		va_end(args);
	}

	historic_layer_entry_t *entries = g_new0(historic_layer_entry_t, n);
	guint built = 0;

	for (GSList *l = layers; l; l = l->next, built++) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay) {
			siril_log_color_message(
				_("undo_save_flis_multi_layer: NULL layer in list — aborting\n"), "red");
			goto fail;
		}
		if (build_layer_entry(lay, &entries[built])) {
			siril_log_color_message(
				_("undo_save_flis_multi_layer: failed to build entry for layer %d\n"),
				"red", lay->item_id);
			goto fail;
		}
	}

	/* All entries built successfully — push onto the ring buffer */
	if (!com.history) {
		com.hist_size    = HISTORY_SIZE;
		com.history      = calloc(com.hist_size, sizeof(historic));
		com.hist_current = 0;
		com.hist_display = 0;
	}
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}
	if (com.hist_current == com.hist_size - 1) {
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
		        (com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}

	{
		historic *h = &com.history[com.hist_current];
		memset(h, 0, sizeof(*h));
		h->flis_layer_id        = FLIS_UNDO_LAYER_NONE;
		h->lmask_layer_id       = FLIS_UNDO_LAYER_NONE;
		h->lmask_dest_layer_id  = FLIS_UNDO_LAYER_NONE;
		h->reorder_layer_a_id   = FLIS_UNDO_LAYER_NONE;
		h->reorder_layer_b_id   = FLIS_UNDO_LAYER_NONE;
		h->multi_entries        = entries;
		h->n_multi_entries      = n;
		snprintf(h->history, FLEN_VALUE, "%s", histo);
	}

	com.hist_current++;
	com.hist_display = com.hist_current;
	gui_function(update_MenuItem, NULL);
	return 0;

fail:
	for (guint k = 0; k < built; k++)
		free_layer_entry(&entries[k]);
	g_free(entries);
	return 1;
}

/* Shared ring-buffer push boilerplate.  Allocates the ring if needed, trims
 * any redo states above hist_display, and handles wrap-around.
 * Returns a zeroed historic* ready to be filled in, or NULL on OOM. */
static historic *undo_ring_prepare(void) {
	if (!com.history) {
		com.hist_size    = HISTORY_SIZE;
		com.history      = calloc(com.hist_size, sizeof(historic));
		if (!com.history) { PRINT_ALLOC_ERR; return NULL; }
		com.hist_current = 0;
		com.hist_display = 0;
	}
	while (com.hist_display < com.hist_current) {
		com.hist_current--;
		undo_remove_item(com.history, com.hist_current);
	}
	if (com.hist_current == com.hist_size - 1) {
		undo_remove_item(com.history, 1);
		memmove(&com.history[1], &com.history[2],
		        (com.hist_size - 2) * sizeof(*com.history));
		com.hist_current = com.hist_size - 2;
	}
	historic *h = &com.history[com.hist_current];
	memset(h, 0, sizeof(*h));
	h->flis_layer_id       = FLIS_UNDO_LAYER_NONE;
	h->lmask_layer_id      = FLIS_UNDO_LAYER_NONE;
	h->lmask_dest_layer_id = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_a_id  = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_b_id  = FLIS_UNDO_LAYER_NONE;
	return h;
}

int undo_save_flis_multi_layer_props(GSList *layers, const char *message, ...) {
	if (!layers || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	guint n = g_slist_length(layers);
	if (n == 0) return 1;

	char histo[FLEN_VALUE] = { 0 };
	if (message) {
		va_list args;
		va_start(args, message);
		vsnprintf(histo, FLEN_VALUE, message, args);
		va_end(args);
	}

	historic_layer_entry_t *entries = g_new0(historic_layer_entry_t, n);
	guint built = 0;

	for (GSList *l = layers; l; l = l->next, built++) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay) {
			siril_log_color_message(
				_("undo_save_flis_multi_layer_props: NULL layer — aborting\n"), "red");
			g_free(entries);
			return 1;
		}
		historic_layer_entry_t *e = &entries[built];
		e->flis_layer_id = lay->item_id;
		e->props_only    = TRUE;
		e->layer_props   = g_new(flis_layer_props_t, 1);
		e->layer_props->blend_mode   = lay->blend_mode;
		e->layer_props->opacity      = lay->opacity;
		e->layer_props->visible      = lay->visible;
		e->layer_props->locked       = lay->locked;
		e->layer_props->has_tint     = lay->has_tint;
		e->layer_props->tint         = lay->layer_tint;
		e->layer_props->lmask_active = lay->lmask_active;
		e->layer_props->position_x   = lay->position_x;
		e->layer_props->position_y   = lay->position_y;
		g_strlcpy(e->layer_props->name,
		          lay->layer_name ? lay->layer_name : "",
		          sizeof(e->layer_props->name));
	}

	historic *h = undo_ring_prepare();
	if (!h) { g_free(entries); return 1; }

	h->flis_layer_id        = FLIS_UNDO_LAYER_NONE;
	h->lmask_layer_id       = FLIS_UNDO_LAYER_NONE;
	h->lmask_dest_layer_id  = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_a_id   = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_b_id   = FLIS_UNDO_LAYER_NONE;
	h->multi_entries        = entries;
	h->n_multi_entries      = n;
	snprintf(h->history, FLEN_VALUE, "%s", histo);

	com.hist_current++;
	com.hist_display = com.hist_current;
	gui_function(update_MenuItem, NULL);
	return 0;
}

int undo_save_processing_mask(fits *fit, const char *message, ...) {
	if (!fit || !single_image_is_loaded()) return 1;

	char histo[FLEN_VALUE] = { 0 };
	if (message) {
		va_list args;
		va_start(args, message);
		vsnprintf(histo, FLEN_VALUE, message, args);
		va_end(args);
	}

	historic *h = undo_ring_prepare();
	if (!h) return 1;

	h->pmask_only  = TRUE;
	h->rx          = fit->rx;
	h->ry          = fit->ry;
	h->mask_bitpix = (fit->mask && fit->mask->data) ? fit->mask->bitpix : 0;
	snprintf(h->history, FLEN_VALUE, "%s", histo);

	if (undo_build_mask_swapfile(fit, &h->mask_filename)) {
		memset(h, 0, sizeof(*h));
		h->flis_layer_id = FLIS_UNDO_LAYER_NONE;
		return 1;
	}

	if (is_current_image_flis()) {
		flis_layer_t *active = flis_active_layer();
		if (active)
			h->flis_layer_id = active->item_id;
	}

	com.hist_current++;
	com.hist_display = com.hist_current;
	gui_function(update_MenuItem, NULL);
	return 0;
}

int undo_save_flis_layer_full(fits *fit_snapshot,
                               flis_layer_t *lay,
                               layermask_t *lmask_snapshot,
                               const flis_layer_props_t *props_snapshot,
                               const char *message, ...) {
	if (!fit_snapshot || !lay || !props_snapshot || !single_image_is_loaded())
		return 1;
	if (!is_current_image_flis()) return 1;

	char histo[FLEN_VALUE] = { 0 };
	if (message) {
		va_list args;
		va_start(args, message);
		vsnprintf(histo, FLEN_VALUE, message, args);
		va_end(args);
	}

	historic *h = undo_ring_prepare();
	if (!h) return 1;

	/* Pixel data */
	if (undo_build_swapfile(fit_snapshot, &h->filename)) {
		memset(h, 0, sizeof(*h)); h->flis_layer_id = FLIS_UNDO_LAYER_NONE; return 1;
	}
	/* Processing mask */
	if (undo_build_mask_swapfile(fit_snapshot, &h->mask_filename)) {
		if (h->filename) { g_unlink(h->filename); g_free(h->filename); h->filename = NULL; }
		memset(h, 0, sizeof(*h)); h->flis_layer_id = FLIS_UNDO_LAYER_NONE; return 1;
	}
	/* Layer mask */
	if (undo_build_lmask_swapfile(lmask_snapshot,
	                               &h->lmask_filename,
	                               &h->lmask_w, &h->lmask_h, &h->lmask_bitpix)) {
		if (h->mask_filename) { g_unlink(h->mask_filename); g_free(h->mask_filename); h->mask_filename = NULL; }
		if (h->filename) { g_unlink(h->filename); g_free(h->filename); h->filename = NULL; }
		memset(h, 0, sizeof(*h)); h->flis_layer_id = FLIS_UNDO_LAYER_NONE; return 1;
	}

	h->full_layer     = TRUE;
	h->flis_layer_id  = lay->item_id;
	h->rx             = fit_snapshot->rx;
	h->ry             = fit_snapshot->ry;
	h->nchans         = (int)fit_snapshot->naxes[2];
	h->type           = fit_snapshot->type;
	h->mask_bitpix    = (fit_snapshot->mask && fit_snapshot->mask->data)
	                      ? fit_snapshot->mask->bitpix : 0;
	h->wcsdata        = fit_snapshot->keywords.wcsdata;
	if (fit_snapshot->keywords.wcslib) {
		int status = -1;
		h->wcslib = wcs_deepcopy(fit_snapshot->keywords.wcslib, &status);
		if (status) siril_debug_print("undo_save_flis_layer_full: wcs copy failed\n");
	}
	h->focal_length   = fit_snapshot->keywords.focal_length;
	h->icc_profile    = copyICCProfile(fit_snapshot->icc_profile);
	/* position_x/y captured in layer_props; also store in flis_position for
	 * the position-restore path that runs at the end of undo_get_data. */
	h->flis_position_x = props_snapshot->position_x;
	h->flis_position_y = props_snapshot->position_y;

	h->layer_props = g_new(flis_layer_props_t, 1);
	*h->layer_props = *props_snapshot;

	snprintf(h->history, FLEN_VALUE, "%s", histo);

	com.hist_current++;
	com.hist_display = com.hist_current;
	gui_function(update_MenuItem, NULL);
	return 0;
}

void set_undo_redo_tooltip() {
	if (is_undo_available()) {
		gchar *str = g_strdup_printf(_("Undo: \"%s\""), com.history[com.hist_display - 1].history);
		gtk_widget_set_tooltip_text(lookup_widget("header_undo_button"), str);
		g_free(str);
	}
	else gtk_widget_set_tooltip_text(lookup_widget("header_undo_button"), _("Nothing to undo"));
	if (is_redo_available()) {
		gchar *str = g_strdup_printf(_("Redo: \"%s\""), com.history[com.hist_display].history);
		gtk_widget_set_tooltip_text(lookup_widget("header_redo_button"), str);
		g_free(str);
	}
	else gtk_widget_set_tooltip_text(lookup_widget("header_redo_button"), _("Nothing to redo"));
}
