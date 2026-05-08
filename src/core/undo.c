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
#include <wcslib.h>

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

	// Write some data to the temporary file
	errno = 0;
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

	errno = 0;
	if (-1 == write(fd, fit->mask->data, n_pixels * elem_size)) {
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

static void undo_free_item(historic *h) {
	if (h->filename) {
		if (g_unlink(h->filename))
			siril_debug_print("g_unlink() failed\n");
		g_free(h->filename);
	}
	if (h->mask_filename) {
		if (g_unlink(h->mask_filename))
			siril_debug_print("g_unlink() of mask failed\n");
		g_free(h->mask_filename);
	}
	if (h->wcslib) {
		wcsfree(h->wcslib);
		free(h->wcslib);
	}
	if (h->icc_profile)
		cmsCloseProfile(h->icc_profile);
	g_free(h);
}

/* Save current gfit pixels to a new swap file and push the resulting historic
 * entry to *stack. label is the operation name (may be empty, not NULL). */
static int undo_push_to(GList **stack, fits *fit, const char *label) {
	gchar *filename = NULL, *mask_filename = NULL;

	if (undo_build_swapfile(fit, &filename))
		return 1;
	if (undo_build_mask_swapfile(fit, &mask_filename)) {
		g_free(filename);
		return 1;
	}

	historic *h = g_new0(historic, 1);
	h->filename = filename;
	h->mask_filename = mask_filename;
	h->mask_bitpix = (fit->mask && fit->mask->data) ? fit->mask->bitpix : 0;
	h->rx = fit->rx;
	h->ry = fit->ry;
	h->nchans = fit->naxes[2];
	h->type = fit->type;
	h->wcsdata = fit->keywords.wcsdata;
	int status = -1;
	h->wcslib = wcs_deepcopy(fit->keywords.wcslib, &status);
	if (status)
		siril_debug_print("could not copy wcslib struct\n");
	h->focal_length = fit->keywords.focal_length;
	h->icc_profile = copyICCProfile(fit->icc_profile);
	snprintf(h->history, FLEN_VALUE, "%s", label ? label : "");

	*stack = g_list_prepend(*stack, h);
	return 0;
}

static int undo_get_data_ushort(fits *fit, historic *hist) {
	int fd;

	// read the data from temporary file
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
	/* read() returns ssize_t; comparing it directly against size_t silently
	 * promotes -1 to SIZE_MAX and hides errors. Capture the signed return. */
	ssize_t bytes_read = read(fd, buf, size);
	if (bytes_read < 0 || (size_t)bytes_read < size) {
		printf("Undo Read of [%s], failed with error [%s]\n", hist->filename, strerror(errno));
		free(buf);
		g_close(fd, NULL);
		return 1;
	}
	/* need to reallocate data as size may have changed */
	WORD *newdata = (WORD*) realloc(fit->data, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
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
	/* Always free the existing wcslib first: assigning the deepcopy result
	 * directly to fit->keywords.wcslib would leak a previously-set struct. */
	free_wcs(fit);
	if (hist->wcslib) {
		int status = -1;
		fit->keywords.wcslib = wcs_deepcopy(hist->wcslib, &status);
		if (status)
			siril_debug_print("could not copy wcslib struct\n");
	} else {
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

	// read the data from temporary file
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
	ssize_t bytes_read = read(fd, buf, size);
	if (bytes_read < 0 || (size_t)bytes_read < size) {
		printf("Undo Read of [%s], failed with error [%s]\n", hist->filename, strerror(errno));
		free(buf);
		g_close(fd, NULL);
		return 1;
	}
	/* need to reallocate data as size may have changed */
	float *newdata = (float*) realloc(fit->fdata, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
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
	free_wcs(fit);
	if (hist->wcslib) {
		int status = -1;
		fit->keywords.wcslib = wcs_deepcopy(hist->wcslib, &status);
		if (status)
			siril_debug_print("could not copy wcslib struct\n");
	} else {
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
	ssize_t bytes_read = read(fd, buf, size);
	if (bytes_read < 0 || (size_t)bytes_read < size) {
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

static int undo_restore(fits *fit, historic *hist) {
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
	if (retval == 0)
		retval = undo_get_mask_data(fit, hist);

	return retval;
}

gboolean is_undo_available() {
	return com.undo_stack != NULL;
}

gboolean is_redo_available() {
	return com.redo_stack != NULL;
}

int undo_save_state(fits *fit, const char *message, ...) {
	if (!single_image_is_loaded())
		return 0;

	char histo[FLEN_VALUE] = { 0 };
	if (message != NULL) {
		va_list args;
		va_start(args, message);
		vsnprintf(histo, FLEN_VALUE, message, args);
		va_end(args);
	}

	/* discard redo stack: a new operation invalidates the redo branch */
	g_list_free_full(com.redo_stack, (GDestroyNotify) undo_free_item);
	com.redo_stack = NULL;

	if (undo_push_to(&com.undo_stack, fit, histo))
		return 1;

	gui_function(update_MenuItem, NULL);
	return 0;
}

int undo_display_data(int dir) {
	switch (dir) {
	case UNDO:
		if (is_undo_available()) {
			// Avoid any issues with ROI or preview
			gboolean preview_was_active = is_preview_active();
			historic *top = (historic *) com.undo_stack->data;
			// Can't reactivate the ROI if the size has changed
			gboolean roi_was_active = (gui.roi.active
					&& gfit->rx == top->rx
					&& gfit->ry == top->ry
					&& gfit->naxes[2] == top->nchans);
			rectangle roi_rect;
			memcpy(&roi_rect, &gui.roi.selection, sizeof(rectangle));

			/* siril_preview_hide() and on_clear_roi() can transitively reach
			 * notify_gfit_data_modified() → copy_roi_into_gfit(), which acquires
			 * gfit->rwlock as a writer. They MUST run with no gfit lock held —
			 * otherwise the same thread re-entering the writer lock self-deadlocks
			 * (silently, since GLib's GRWLock makes recursive lock UB and Linux
			 * pthread_rwlock blocks indefinitely). */
			siril_preview_hide();
			on_clear_roi();

			/* Writer lock: undo_push_to (reads pixels) and undo_restore
			 * (writes pixels) must be atomic against the Python thread. */
			g_rw_lock_writer_lock(&gfit->rwlock);

			/* save current state to redo stack before restoring */
			if (undo_push_to(&com.redo_stack, gfit, top->history)) {
				g_rw_lock_writer_unlock(&gfit->rwlock);
				return 1;
			}

			siril_log_message(_("Undo: %s\n"), top->history);

			/* pop and restore */
			com.undo_stack = g_list_remove_link(com.undo_stack, com.undo_stack);
			undo_restore(gfit, top);
			undo_free_item(top);

			invalidate_gfit_histogram();
			invalidate_stats_from_fit(gfit);
			g_rw_lock_writer_unlock(&gfit->rwlock); // Finished with writer lock
			g_rw_lock_reader_lock(&gfit->rwlock);   // But still need reader lock
			update_gfit_histogram_if_needed();
			gui_function(close_tab, NULL); // These 2 lines account for possible change from mono to RGB
			g_rw_lock_reader_unlock(&gfit->rwlock);
			/* update menus */
			gui_function(update_MenuItem, NULL);
			lock_display_transform();
			if (gui.icc.proofing_transform)
				cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
			unlock_display_transform();
			refresh_annotations(TRUE);
			/* redraw_mask_idle posts an idle — must be called outside any gfit lock */
			redraw_mask_idle(NULL);
			if (!com.pref.gui.mask_tints_vports) { // redraw() is called in redraw_mask_idle if this is TRUE
				/* No reader lock here: notify_gfit_data_modified() may call
				 * copy_roi_into_gfit() which acquires the writer lock —
				 * mirrors the contract redraw_mask_idle relies on. */
				notify_gfit_data_modified();
				redraw(REMAP_ALL);
			}
			if (preview_was_active) {
				g_rw_lock_reader_lock(&gfit->rwlock);
				copy_gfit_to_backup();
				g_rw_lock_reader_unlock(&gfit->rwlock);
				// TODO: To be perfect, we would need a register of preview functions
				// look up the correct one for the open dialog and re-apply the preview
				siril_log_message(_("Following undo / redo with a preview active you may need "
						"to toggle the preview off and on again to reactivate the preview effect\n"));
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
			historic *top = (historic *) com.redo_stack->data;
			// Can't reactivate the ROI if the size has changed
			gboolean roi_was_active = (gui.roi.active
					&& gfit->rx == top->rx
					&& gfit->ry == top->ry
					&& gfit->naxes[2] == top->nchans);
			rectangle roi_rect;
			memcpy(&roi_rect, &gui.roi.selection, sizeof(rectangle));

			/* See UNDO case: these can transitively try to take the writer
			 * lock and must run with no gfit lock held. */
			on_clear_roi();
			siril_preview_hide();

			/* Writer lock: undo_push_to (reads pixels) and undo_restore
			 * (writes pixels) must be atomic against the Python thread. */
			g_rw_lock_writer_lock(&gfit->rwlock);

			/* save current state to undo stack before restoring */
			if (undo_push_to(&com.undo_stack, gfit, top->history)) {
				g_rw_lock_writer_unlock(&gfit->rwlock);
				return 1;
			}

			siril_log_message(_("Redo: %s\n"), top->history);

			/* pop and restore */
			com.redo_stack = g_list_remove_link(com.redo_stack, com.redo_stack);
			undo_restore(gfit, top);
			undo_free_item(top);

			invalidate_gfit_histogram();
			invalidate_stats_from_fit(gfit);
			g_rw_lock_writer_unlock(&gfit->rwlock); // Finished with writer lock
			g_rw_lock_reader_lock(&gfit->rwlock);   // But still need reader lock
			update_gfit_histogram_if_needed();
			g_rw_lock_reader_unlock(&gfit->rwlock);
			/* update menus */
			gui_function(update_MenuItem, NULL);
			refresh_annotations(TRUE);
			lock_display_transform();
			if (gui.icc.proofing_transform)
				cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
			unlock_display_transform();
			gui_function(close_tab, NULL); // These 2 lines account for possible change from mono to RGB
			/* redraw_mask_idle posts an idle — must be called outside any gfit lock */
			redraw_mask_idle(NULL);
			if (!com.pref.gui.mask_tints_vports) { // redraw() is called in redraw_mask_idle if this is TRUE
				/* No reader lock: notify_gfit_data_modified() may call
				 * copy_roi_into_gfit() which acquires the writer lock. */
				notify_gfit_data_modified();
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
	g_list_free_full(com.undo_stack, (GDestroyNotify) undo_free_item);
	com.undo_stack = NULL;
	g_list_free_full(com.redo_stack, (GDestroyNotify) undo_free_item);
	com.redo_stack = NULL;
	return 0;
}

void set_undo_redo_tooltip() {
	if (is_undo_available()) {
		historic *h = (historic *) com.undo_stack->data;
		gchar *str = g_strdup_printf(_("Undo: \"%s\""), h->history);
		gtk_widget_set_tooltip_text(lookup_widget("header_undo_button"), str);
		g_free(str);
	} else {
		gtk_widget_set_tooltip_text(lookup_widget("header_undo_button"), _("Nothing to undo"));
	}
	if (is_redo_available()) {
		historic *h = (historic *) com.redo_stack->data;
		gchar *str = g_strdup_printf(_("Redo: \"%s\""), h->history);
		gtk_widget_set_tooltip_text(lookup_widget("header_redo_button"), str);
		g_free(str);
	} else {
		gtk_widget_set_tooltip_text(lookup_widget("header_redo_button"), _("Nothing to redo"));
	}
}

/* ---- Long-press popover for undo/redo history navigation ---- */

static gboolean destroy_widget_idle(gpointer data) {
	gtk_widget_destroy(GTK_WIDGET(data));
	g_object_unref(G_OBJECT(data));
	return G_SOURCE_REMOVE;
}

static void on_undo_popover_closed(GtkPopover *popover, gpointer user_data) {
	/* Schedule destruction via idle to avoid re-entrancy: gtk_popover_popdown
	 * emits "closed" while still holding internal GTK state on the widget. Calling
	 * gtk_widget_destroy synchronously from here crashes in g_type_check_instance
	 * because popdown continues to access the widget after the signal returns. */
	g_object_ref(popover);  /* keep alive until the idle fires */
	g_idle_add(destroy_widget_idle, popover);
}

static void on_undo_popover_row_activated(GtkListBox *box, GtkListBoxRow *row, gpointer user_data) {
	GtkWidget *popover = GTK_WIDGET(g_object_get_data(G_OBJECT(box), "popover"));
	int dir   = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(row), "dir"));
	int level = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(row), "level"));

	/* Hide the popover (fires "closed" → on_undo_popover_closed → idle destroy).
	 * We use gtk_widget_hide rather than gtk_popover_popdown because popdown
	 * keeps internal GTK references alive across the "closed" emission, making
	 * synchronous destruction unsafe. dir/level are already local copies. */
	gtk_widget_hide(GTK_WIDGET(popover));

	set_cursor_waiting(TRUE);
	for (int i = 0; i < level; i++) {
		if (dir == UNDO && !is_undo_available()) break;
		if (dir == REDO && !is_redo_available()) break;
		undo_display_data(dir);
	}
	set_cursor_waiting(FALSE);
}

static void show_undo_history_popover(GtkWidget *button, int dir) {
	GList *stack = (dir == UNDO) ? com.undo_stack : com.redo_stack;
	if (!stack) return;

	GtkWidget *popover = gtk_popover_new(button);
	g_signal_connect(popover, "closed", G_CALLBACK(on_undo_popover_closed), NULL);

	GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_widget_set_margin_start(vbox, 4);
	gtk_widget_set_margin_end(vbox, 4);
	gtk_widget_set_margin_top(vbox, 6);
	gtk_widget_set_margin_bottom(vbox, 4);

	/* heading */
	const gchar *title = (dir == UNDO) ? _("Undo history") : _("Redo history");
	GtkWidget *heading = gtk_label_new(title);
	PangoAttrList *attrs = pango_attr_list_new();
	pango_attr_list_insert(attrs, pango_attr_weight_new(PANGO_WEIGHT_BOLD));
	gtk_label_set_attributes(GTK_LABEL(heading), attrs);
	pango_attr_list_unref(attrs);
	gtk_widget_set_margin_bottom(heading, 4);
	gtk_box_pack_start(GTK_BOX(vbox), heading, FALSE, FALSE, 0);

	GtkWidget *sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_widget_set_margin_bottom(sep, 2);
	gtk_box_pack_start(GTK_BOX(vbox), sep, FALSE, FALSE, 0);

	GtkWidget *scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
	                               GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_max_content_height(GTK_SCROLLED_WINDOW(scroll), 300);
	gtk_scrolled_window_set_propagate_natural_height(GTK_SCROLLED_WINDOW(scroll), TRUE);

	GtkWidget *listbox = gtk_list_box_new();
	gtk_list_box_set_selection_mode(GTK_LIST_BOX(listbox), GTK_SELECTION_NONE);
	g_object_set_data(G_OBJECT(listbox), "dir",     GINT_TO_POINTER(dir));
	g_object_set_data(G_OBJECT(listbox), "popover", popover);
	g_signal_connect(listbox, "row-activated", G_CALLBACK(on_undo_popover_row_activated), NULL);

	int n = 0;
	for (GList *l = stack; l; l = l->next, n++) {
		historic *h = (historic *)l->data;
		const gchar *label = (h->history[0] != '\0') ? h->history : _("(unnamed)");

		GtkWidget *row_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		gtk_widget_set_margin_start(row_box, 6);
		gtk_widget_set_margin_end(row_box, 6);
		gtk_widget_set_margin_top(row_box, 3);
		gtk_widget_set_margin_bottom(row_box, 3);

		GtkWidget *lbl = gtk_label_new(label);
		gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
		gtk_box_pack_start(GTK_BOX(row_box), lbl, TRUE, TRUE, 0);

		GtkWidget *list_row = gtk_list_box_row_new();
		gtk_container_add(GTK_CONTAINER(list_row), row_box);
		g_object_set_data(G_OBJECT(list_row), "dir",   GINT_TO_POINTER(dir));
		g_object_set_data(G_OBJECT(list_row), "level", GINT_TO_POINTER(n + 1));
		gtk_list_box_insert(GTK_LIST_BOX(listbox), list_row, -1);
	}

	gtk_container_add(GTK_CONTAINER(scroll), listbox);
	gtk_box_pack_start(GTK_BOX(vbox), scroll, TRUE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(popover), vbox);

	gtk_widget_show_all(popover);
	gtk_popover_popup(GTK_POPOVER(popover));
}

static void on_long_press_undo(GtkGestureLongPress *gesture, gdouble x, gdouble y, gpointer user_data) {
	gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
	show_undo_history_popover(GTK_WIDGET(user_data), UNDO);
}

static void on_long_press_redo(GtkGestureLongPress *gesture, gdouble x, gdouble y, gpointer user_data) {
	gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);
	show_undo_history_popover(GTK_WIDGET(user_data), REDO);
}

void setup_undo_redo_long_press(void) {
	GtkWidget *undo_btn = lookup_widget("header_undo_button");
	GtkWidget *redo_btn = lookup_widget("header_redo_button");

	/* Intentionally not unreffed: the widget holds no strong ref in GTK 3.24,
	 * so we keep the reference alive for the lifetime of the application. */
	GtkGesture *undo_gesture = gtk_gesture_long_press_new(undo_btn);
	gtk_gesture_single_set_touch_only(GTK_GESTURE_SINGLE(undo_gesture), FALSE);
	gtk_event_controller_set_propagation_phase(GTK_EVENT_CONTROLLER(undo_gesture), GTK_PHASE_CAPTURE);
	g_signal_connect(undo_gesture, "pressed", G_CALLBACK(on_long_press_undo), undo_btn);

	GtkGesture *redo_gesture = gtk_gesture_long_press_new(redo_btn);
	gtk_gesture_single_set_touch_only(GTK_GESTURE_SINGLE(redo_gesture), FALSE);
	gtk_event_controller_set_propagation_phase(GTK_EVENT_CONTROLLER(redo_gesture), GTK_PHASE_CAPTURE);
	g_signal_connect(redo_gesture, "pressed", G_CALLBACK(on_long_press_redo), redo_btn);
}
