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
#ifdef _WIN32
#include <windows.h>
#include <io.h>  /* _get_osfhandle */
#endif

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "core/gui_iface.h"
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

/* Mark fd so the file is deleted automatically when the last handle closes.
 * On POSIX: unlink removes the directory entry immediately; the data lives
 * until the fd is closed (even on SIGKILL).
 * On Windows: SetFileInformationByHandle with FileDispositionInfo achieves the
 * same effect - the OS deletes the file when the process exits or the fd is
 * closed, including on abnormal termination. */
static void swap_mark_delete_on_close(int fd, const gchar *path) {
#ifndef _WIN32
	(void)fd;
	g_unlink(path);
#else
	(void)path;
	HANDLE h = (HANDLE)_get_osfhandle(fd);
	if (h != INVALID_HANDLE_VALUE) {
		FILE_DISPOSITION_INFO fdi = { .DeleteFile = TRUE };
		SetFileInformationByHandle(h, FileDispositionInfo, &fdi, sizeof(fdi));
	}
#endif
}

/* Returns an open fd on success, -1 on error.
 * The file is marked delete-on-close on both POSIX and Windows: it vanishes
 * automatically when the fd is closed, even if Siril is killed. */
static int undo_build_swapfile(fits *fit) {
	char name[] = "siril_swp-XXXXXX";
	gchar *nameBuff = g_build_filename(com.pref.swap_dir, name, NULL);
#ifndef _WIN32
	int fd = g_mkstemp(nameBuff);
#else
	int fd = g_mkstemp_full(nameBuff, _O_RDWR | _O_CREAT | _O_EXCL | _O_BINARY | _O_TEMPORARY, 
    _S_IREAD | _S_IWRITE);
#endif
	if (fd < 0) {
		siril_log_error(_("File I/O Error: Unable to create swap file in %s: [%s]\n"),
				com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return -1;
	}
	swap_mark_delete_on_close(fd, nameBuff);
	g_free(nameBuff);

	size_t size = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	errno = 0;
	ssize_t written = -1;
	if (fit->type == DATA_USHORT)
		written = write(fd, fit->data, size * sizeof(WORD));
	else if (fit->type == DATA_FLOAT)
		written = write(fd, fit->fdata, size * sizeof(float));

	if (written == -1) {
		siril_log_error(_("File I/O Error: Unable to write swap file: [%s]\n"), strerror(errno));
		g_close(fd, NULL);
		return -1;
	}
	return fd;
}

/* Returns an open fd on success, -1 if there is no mask (not an error), -2 on error. */
static int undo_build_mask_swapfile(fits *fit) {
	if (!fit->mask || !fit->mask->data)
		return -1;  /* no mask - not an error */

	char name[] = "siril_msk-XXXXXX";
	gchar *nameBuff = g_build_filename(com.pref.swap_dir, name, NULL);
#ifndef _WIN32
	int fd = g_mkstemp(nameBuff);
#else
	int fd = g_mkstemp_full(nameBuff, _O_RDWR | _O_CREAT | _O_EXCL | _O_BINARY | _O_TEMPORARY, 
    _S_IREAD | _S_IWRITE);
#endif
	if (fd < 0) {
		siril_log_error(_("File I/O Error: Unable to create mask swap file in %s: [%s]\n"),
				com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return -2;
	}
	swap_mark_delete_on_close(fd, nameBuff);
	g_free(nameBuff);

	size_t n_pixels = fit->naxes[0] * fit->naxes[1];
	size_t elem_size;
	switch (fit->mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			siril_log_error(_("Error: Invalid mask bitpix value: %d\n"), fit->mask->bitpix);
			g_close(fd, NULL);
			return -2;
	}

	errno = 0;
	if (-1 == write(fd, fit->mask->data, n_pixels * elem_size)) {
		siril_log_error(_("File I/O Error: Unable to write mask swap file: [%s]\n"), strerror(errno));
		g_close(fd, NULL);
		return -2;
	}
	return fd;
}

static void undo_free_item(historic *h) {
	/* Closing the fd is all that is needed on both platforms: the file was
	 * marked delete-on-close at creation (POSIX unlink / Windows
	 * FileDispositionInfo) so the OS reclaims it here automatically. */
	if (h->fd >= 0)
		g_close(h->fd, NULL);
	if (h->mask_fd >= 0)
		g_close(h->mask_fd, NULL);
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
	int fd = undo_build_swapfile(fit);
	if (fd < 0)
		return 1;

	int mask_fd = undo_build_mask_swapfile(fit);
	if (mask_fd == -2) {  /* -1 means "no mask", -2 means error */
		g_close(fd, NULL);  /* delete-on-close cleans up the swap file */
		return 1;
	}

	historic *h = g_new0(historic, 1);
	h->fd = fd;
	h->mask_fd = mask_fd;  /* -1 if no mask, >= 0 if mask exists */
	h->mask_bitpix = (fit->mask && fit->mask->data) ? fit->mask->bitpix : 0;
	h->rx = fit->rx;
	h->ry = fit->ry;
	h->nchans = fit->naxes[2];
	h->type = fit->type;
	h->wcsdata = fit->keywords.wcsdata;
	int status = -1;
	h->wcslib = wcs_deepcopy(fit->keywords.wcslib, &status);
	if (status)
		siril_log_debug("could not copy wcslib struct\n");
	h->focal_length = fit->keywords.focal_length;
	h->icc_profile = copyICCProfile(fit->icc_profile);
	snprintf(h->history, FLEN_VALUE, "%s", label ? label : "");

	*stack = g_list_prepend(*stack, h);
	return 0;
}

static int undo_get_data_ushort(fits *fit, historic *hist) {
	if (lseek(hist->fd, 0, SEEK_SET) == (off_t)-1) {
		printf("Error seeking swap file: [%s]\n", strerror(errno));
		return 1;
	}

	errno = 0;
	fit->rx = fit->naxes[0] = hist->rx;
	fit->ry = fit->naxes[1] = hist->ry;

	size_t n = fit->naxes[0] * fit->naxes[1];
	size_t size = n * fit->naxes[2] * sizeof(WORD);
	WORD *buf = calloc(1, size);
	/* Cast size to ssize_t so a -1 return from read() compares as negative,
	 * not as SIZE_MAX after silent promotion (which would mask errors). */
	if (read(hist->fd, buf, size) < (ssize_t)size) {
		printf("Undo Read failed with error [%s]\n", strerror(errno));
		free(buf);
		return 1;
	}
	/* need to reallocate data as size may have changed */
	WORD *newdata = (WORD*) realloc(fit->data, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		free(buf);
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
			siril_log_debug("could not copy wcslib struct\n");
	} else {
		reset_wcsdata(fit);
	}
	fit->keywords.focal_length = hist->focal_length;

	full_stats_invalidation_from_fit(fit);
	free(buf);
	return 0;
}

static int undo_get_data_float(fits *fit, historic *hist) {
	if (lseek(hist->fd, 0, SEEK_SET) == (off_t)-1) {
		printf("Error seeking swap file: [%s]\n", strerror(errno));
		return 1;
	}

	errno = 0;
	fit->rx = fit->naxes[0] = hist->rx;
	fit->ry = fit->naxes[1] = hist->ry;

	size_t n = fit->naxes[0] * fit->naxes[1];
	size_t size = n * fit->naxes[2] * sizeof(float);
	float *buf = calloc(1, size);
	if (read(hist->fd, buf, size) < (ssize_t)size) {
		printf("Undo Read failed with error [%s]\n", strerror(errno));
		free(buf);
		return 1;
	}
	/* need to reallocate data as size may have changed */
	float *newdata = (float*) realloc(fit->fdata, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		free(buf);
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
			siril_log_debug("could not copy wcslib struct\n");
	} else {
		reset_wcsdata(fit);
	}
	fit->keywords.focal_length = hist->focal_length;

	full_stats_invalidation_from_fit(fit);
	free(buf);
	return 0;
}

static int undo_get_mask_data(fits *fit, historic *hist) {
	/* mask_fd == -1 means no mask was saved */
	if (hist->mask_fd < 0) {
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

	if (lseek(hist->mask_fd, 0, SEEK_SET) == (off_t)-1) {
		printf("Error seeking mask swap file: [%s]\n", strerror(errno));
		return 1;
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

	size_t n_pixels = hist->rx * hist->ry;
	size_t elem_size;
	switch (hist->mask_bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			siril_log_error(_("Error: Invalid mask bitpix value in history: %d\n"), hist->mask_bitpix);
			return 1;
	}

	size_t size = n_pixels * elem_size;
	void *buf = calloc(1, size);
	if (!buf) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	errno = 0;
	if (read(hist->mask_fd, buf, size) < (ssize_t)size) {
		printf("Undo Read of mask failed with error [%s]\n", strerror(errno));
		free(buf);
		return 1;
	}

	/* Reallocate mask data if needed */
	void *newdata = realloc(fit->mask->data, size);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		free(buf);
		return 1;
	}

	fit->mask->data = newdata;
	memcpy(fit->mask->data, buf, size);

	free(buf);
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
			gui_iface.on_precision_changed();
		}
		retval = undo_get_data_ushort(fit, hist);
	} else if (hist->type == DATA_FLOAT) {
		if (gfit->type != DATA_FLOAT) {
			size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
			fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
			gui_iface.on_precision_changed();
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

	gui_iface.update_menu_state();
	return 0;
}

int undo_display_data(int dir) {
	switch (dir) {
	case UNDO:
		if (is_undo_available()) {
			// Avoid any issues with ROI or preview
			gboolean preview_was_active = gui_iface.is_preview_active();
			historic *top = (historic *) com.undo_stack->data;
			// Can't reactivate the ROI if the size has changed
			gboolean roi_was_active = (gui_iface.roi_is_active()
					&& gfit->rx == top->rx
					&& gfit->ry == top->ry
					&& gfit->naxes[2] == top->nchans);
			rectangle roi_rect;
			gui_iface.get_roi_selection(&roi_rect);

			/* hide_preview() and clear_roi() can transitively reach
			 * notify_gfit_data_modified() → copy_roi_into_gfit(), which acquires
			 * gfit->rwlock as a writer. They MUST run with no gfit lock held —
			 * otherwise the same thread re-entering the writer lock self-deadlocks
			 * (silently, since GLib's GRWLock makes recursive lock UB and Linux
			 * pthread_rwlock blocks indefinitely). */
			gui_iface.hide_preview();
			gui_iface.clear_roi();

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

			gui_iface.invalidate_histogram();
			invalidate_stats_from_fit(gfit);
			g_rw_lock_writer_unlock(&gfit->rwlock); // Finished with writer lock
			g_rw_lock_reader_lock(&gfit->rwlock);   // But still need reader lock
			gui_iface.update_histogram();
			gui_iface.on_channel_count_changed(); // These 2 lines account for possible change from mono to RGB
			g_rw_lock_reader_unlock(&gfit->rwlock);
			gui_iface.update_menu_state();
			gui_iface.reset_display_transform();
			refresh_annotations(TRUE);
			/* No gfit lock held here: notify_gfit_data_modified() may call
			 * copy_roi_into_gfit() which acquires the writer lock, and
			 * redraw_mask_idle takes the reader lock itself. */
			notify_gfit_data_modified();
			gui_iface.redraw_image(REDRAW_ALL);
			/* The remap above already applied any mask tint, so only the
			 * mask vport buffer needs refreshing. */
			gui_iface.redraw_mask_idle(FALSE);
			if (preview_was_active) {
				g_rw_lock_reader_lock(&gfit->rwlock);
				gui_iface.copy_gfit_to_backup();
				g_rw_lock_reader_unlock(&gfit->rwlock);
				// TODO: To be perfect, we would need a register of preview functions
				// look up the correct one for the open dialog and re-apply the preview
				siril_log_message(_("Following undo / redo with a preview active you may need "
						"to toggle the preview off and on again to reactivate the preview effect\n"));
			}
			if (roi_was_active) {
				gui_iface.restore_roi(&roi_rect);
			}
			g_rw_lock_reader_lock(&gfit->rwlock);
			update_fits_header(gfit);
			g_rw_lock_reader_unlock(&gfit->rwlock);
		}
		break;

	case REDO:
		if (is_redo_available()) {
			// Avoid any issues with ROI or preview
			gboolean preview_was_active = gui_iface.is_preview_active();
			historic *top = (historic *) com.redo_stack->data;
			// Can't reactivate the ROI if the size has changed
			gboolean roi_was_active = (gui_iface.roi_is_active()
					&& gfit->rx == top->rx
					&& gfit->ry == top->ry
					&& gfit->naxes[2] == top->nchans);
			rectangle roi_rect;
			gui_iface.get_roi_selection(&roi_rect);

			/* See UNDO case: these can transitively try to take the writer
			 * lock via notify_gfit_data_modified() → copy_roi_into_gfit()
			 * and must run with no gfit lock held. */
			gui_iface.clear_roi();
			gui_iface.hide_preview();

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

			gui_iface.invalidate_histogram();
			g_rw_lock_writer_unlock(&gfit->rwlock); // Finished with writer lock
			g_rw_lock_reader_lock(&gfit->rwlock);   // But still need reader lock
			gui_iface.update_histogram();
			g_rw_lock_reader_unlock(&gfit->rwlock);
			gui_iface.update_menu_state();
			refresh_annotations(TRUE);
			gui_iface.reset_display_transform();
			gui_iface.on_channel_count_changed(); // These 2 lines account for possible change from mono to RGB
			/* No gfit lock held here: notify_gfit_data_modified() may call
			 * copy_roi_into_gfit() which acquires the writer lock, and
			 * redraw_mask_idle takes the reader lock itself. */
			notify_gfit_data_modified();
			gui_iface.redraw_image(REDRAW_ALL);
			/* The remap above already applied any mask tint, so only the
			 * mask vport buffer needs refreshing. */
			gui_iface.redraw_mask_idle(FALSE);
			if (preview_was_active) {
				g_rw_lock_reader_lock(&gfit->rwlock);
				gui_iface.copy_gfit_to_backup();
				g_rw_lock_reader_unlock(&gfit->rwlock);
			}
			if (roi_was_active) {
				gui_iface.restore_roi(&roi_rect);
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

