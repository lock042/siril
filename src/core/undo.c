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
#include "io/annotation_catalogues.h"
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
	memset(histo[index].history, 0, FLEN_VALUE);
	return 0;
}

static void undo_add_item(fits *fit, char *filename, const char *histo) {

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

static int undo_get_data(fits *fit, historic *hist) {
	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = copyICCProfile(hist->icc_profile);
	color_manage(fit, (fit->icc_profile != NULL));
	fits_change_depth(fit, hist->nchans);

	if (hist->type == DATA_USHORT) {
		if (gfit->type != DATA_USHORT) {
			size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
			fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
			gui_function(set_precision_switch, NULL);
		}
		return undo_get_data_ushort(fit, hist);
	} else if (hist->type == DATA_FLOAT) {
		if (gfit->type != DATA_FLOAT) {
			size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
			fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
			gui_function(set_precision_switch, NULL);
		}
		return undo_get_data_float(fit, hist);
	}
	return 1;
}

gboolean is_undo_available() {
    return (com.history && com.hist_display > 0);
}

gboolean is_redo_available() {
    return (com.history && (com.hist_display < com.hist_current - 1));
}

int undo_save_state(fits *fit, const char *message, ...) {
	gchar *filename;
	va_list args;
	va_start(args, message);

	if (single_image_is_loaded()) {
		char histo[FLEN_VALUE] = { 0 };

		if (message != NULL) {
			vsnprintf(histo, FLEN_VALUE, message, args);
		}

		if (undo_build_swapfile(fit, &filename)) {
			va_end(args);
			return 1;
		}

		undo_add_item(fit, filename, histo);

		/* update menus */
		gui_function(update_MenuItem, NULL);
	}
	va_end(args);
	return 0;
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
			// Can't reactivate the ROI if the size has changed
			gboolean roi_was_active = (gui.roi.active && gfit->rx == com.history[com.hist_display - 1].rx
					&& gfit->ry == com.history[com.hist_display - 1].ry
					&& gfit->naxes[2] == com.history[com.hist_display - 1].nchans);
			rectangle roi_rect;
			memcpy(&roi_rect, &gui.roi.selection, sizeof(rectangle));
			siril_preview_hide();
			on_clear_roi();
			if (com.hist_current == com.hist_display) {
				undo_save_state(gfit, NULL);
				com.hist_display--;
			}
			com.hist_display--;
			siril_log_message(_("Undo: %s\n"), com.history[com.hist_display].history);
			undo_get_data(gfit, &com.history[com.hist_display]);
			invalidate_gfit_histogram();
			invalidate_stats_from_fit(gfit);
			update_gfit_histogram_if_needed();
			gui_function(update_MenuItem, NULL);
			lock_display_transform();
			if (gui.icc.proofing_transform)
				cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
			unlock_display_transform();
			refresh_annotations(TRUE);
			gui_function(close_tab, NULL); // These 2 lines account for possible change from mono to RGB
			redraw(REMAP_ALL);
			if (preview_was_active) {
				copy_gfit_to_backup();
				siril_log_message(_("Following undo / redo with a preview active you may need "
						"to toggle the preview off and on again to reactivate the preview effect\n"));
				// TODO: To be perfect, we would need a register of preview functions
				// look up the correct one for the open dialog and re-apply the preview
			}
			if (roi_was_active) {
				memcpy(&com.selection, &roi_rect, sizeof(rectangle));
				on_set_roi();
			}
			update_fits_header(gfit);
		}
		break;
	case REDO:
		if (is_redo_available()) {
			// Avoid any issues with ROI or preview
			gboolean preview_was_active = is_preview_active();
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
			undo_get_data(gfit, &com.history[com.hist_display]);
			invalidate_gfit_histogram();
			invalidate_stats_from_fit(gfit);
			update_gfit_histogram_if_needed();
			gui_function(update_MenuItem, NULL);
			refresh_annotations(TRUE);
			lock_display_transform();
			if (gui.icc.proofing_transform)
				cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
			unlock_display_transform();
			gui_function(close_tab, NULL); // These 2 lines account for possible change from mono to RGB
			redraw(REMAP_ALL);
			if (preview_was_active)
				copy_gfit_to_backup();
			if (roi_was_active) {
				memcpy(&gui.roi.selection, &roi_rect, sizeof(rectangle));
				on_set_roi();
			}
			update_fits_header(gfit);
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
