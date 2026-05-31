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
#include "io/image_format_flis.h"
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

/* =====================================================================
 * FLIS swap-file builders — filename-mode (stage 1.4)
 * =====================================================================
 *
 * Filename-mode is required by multi-layer FLIS undo because keeping N
 * open fds per undo entry (one per layer) exhausts the OS fd limit at
 * modest layer counts.  Filename-mode trades the delete-on-close
 * orphan-safety of fd-mode for unbounded multi-file capacity: the
 * caller stores the returned path and must call g_unlink + g_free when
 * the entry is destroyed (undo_free_item handles this uniformly).
 *
 * The fd-mode helpers (undo_build_swapfile / undo_build_mask_swapfile)
 * stay in place above and continue to back undo_save_state.  These
 * _named variants exist purely for FLIS callers.
 * ===================================================================== */

/* Returns 0 on success, 1 on error.  On success *out_filename owns a
 * heap-allocated path that the caller is responsible for. */
static int undo_build_swapfile_named(fits *fit, char **out_filename) {
	*out_filename = NULL;
	gchar *nameBuff = g_build_filename(com.pref.swap_dir, "siril_swp-XXXXXX", NULL);
	int fd = g_mkstemp(nameBuff);
	if (fd < 0) {
		siril_log_error(_("File I/O Error: Unable to create swap file in %s: [%s]\n"),
				com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return 1;
	}

	size_t size = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	errno = 0;
	ssize_t written = -1;
	if (fit->type == DATA_USHORT)
		written = write(fd, fit->data, size * sizeof(WORD));
	else if (fit->type == DATA_FLOAT)
		written = write(fd, fit->fdata, size * sizeof(float));

	g_close(fd, NULL);
	if (written < 0) {
		siril_log_error(_("File I/O Error: Unable to write swap file: [%s]\n"), strerror(errno));
		g_unlink(nameBuff);
		g_free(nameBuff);
		return 1;
	}
	*out_filename = nameBuff;
	return 0;
}

/* Returns 0 on success (whether or not a mask was present); on success
 * *out_filename is non-NULL if a mask was written, NULL if the fit had
 * no mask.  Returns 1 on error. */
static int undo_build_mask_swapfile_named(fits *fit, char **out_filename) {
	*out_filename = NULL;
	if (!fit->mask || !fit->mask->data)
		return 0; /* no mask — success with NULL filename */

	gchar *nameBuff = g_build_filename(com.pref.swap_dir, "siril_msk-XXXXXX", NULL);
	int fd = g_mkstemp(nameBuff);
	if (fd < 0) {
		siril_log_error(_("File I/O Error: Unable to create mask swap file in %s: [%s]\n"),
				com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return 1;
	}

	size_t n_pixels = fit->naxes[0] * fit->naxes[1];
	size_t elem_size;
	switch (fit->mask->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			siril_log_error(_("Error: Invalid mask bitpix value: %d\n"), fit->mask->bitpix);
			g_close(fd, NULL);
			g_unlink(nameBuff);
			g_free(nameBuff);
			return 1;
	}

	errno = 0;
	ssize_t written = write(fd, fit->mask->data, n_pixels * elem_size);
	g_close(fd, NULL);
	if (written < 0) {
		siril_log_error(_("File I/O Error: Unable to write mask swap file: [%s]\n"), strerror(errno));
		g_unlink(nameBuff);
		g_free(nameBuff);
		return 1;
	}
	*out_filename = nameBuff;
	return 0;
}

/* Write a layermask_t to a new temporary swap file.  On success sets
 * *out_filename, *out_w, *out_h, *out_bitpix and returns 0.
 * On a NULL or empty mask, *out_filename is set to NULL and 0 is
 * returned (not an error).  On real error returns 1. */
static int undo_build_lmask_swapfile(layermask_t *lm,
                                     char **out_filename,
                                     size_t *out_w, size_t *out_h,
                                     guint8 *out_bitpix) {
	*out_filename = NULL;
	*out_w = 0;
	*out_h = 0;
	*out_bitpix = 0;

	if (!lm || !lm->data)
		return 0;

	size_t elem_size;
	switch (lm->bitpix) {
		case 8:  elem_size = sizeof(uint8_t);  break;
		case 16: elem_size = sizeof(uint16_t); break;
		case 32: elem_size = sizeof(float);    break;
		default:
			siril_log_warning(_("undo_build_lmask_swapfile: unsupported bitpix %d\n"),
			                  (int)lm->bitpix);
			return 1;
	}

	gchar *nameBuff = g_build_filename(com.pref.swap_dir, "siril_lmsk-XXXXXX", NULL);
	int fd = g_mkstemp(nameBuff);
	if (fd < 0) {
		siril_log_error(_("File I/O Error: Unable to create lmask swap file in %s: [%s]\n"),
		                com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return 1;
	}

	size_t total = lm->w * lm->h * elem_size;
	errno = 0;
	ssize_t written = write(fd, lm->data, total);
	g_close(fd, NULL);
	if (written < 0) {
		siril_log_error(_("File I/O Error: Unable to write lmask swap file: [%s]\n"),
		                strerror(errno));
		g_unlink(nameBuff);
		g_free(nameBuff);
		return 1;
	}

	*out_filename = nameBuff;
	*out_w        = lm->w;
	*out_h        = lm->h;
	*out_bitpix   = lm->bitpix;
	return 0;
}

/* Read an entire swap file into a freshly allocated buffer.  Returns the
 * buffer (caller frees) or NULL on error. */
static void *read_swap_file(const char *filename, size_t expected) {
	int fd = g_open(filename, O_RDONLY | O_BINARY, 0);
	if (fd == -1) {
		siril_log_error(_("Undo: cannot open swap file %s: [%s]\n"),
		                filename, strerror(errno));
		return NULL;
	}
	void *buf = malloc(expected);
	if (!buf) { PRINT_ALLOC_ERR; g_close(fd, NULL); return NULL; }
	ssize_t got = read(fd, buf, expected);
	g_close(fd, NULL);
	if (got < (ssize_t)expected) {
		siril_log_error(_("Undo: short read on swap file %s\n"), filename);
		free(buf);
		return NULL;
	}
	return buf;
}

/* =====================================================================
 * Compound-undo per-layer entry helpers
 * ===================================================================== */

static void free_layer_entry(historic_layer_entry_t *e) {
	if (!e) return;
	if (e->filename)       { g_unlink(e->filename);       g_free(e->filename); }
	if (e->mask_filename)  { g_unlink(e->mask_filename);  g_free(e->mask_filename); }
	if (e->lmask_filename) { g_unlink(e->lmask_filename); g_free(e->lmask_filename); }
	g_free(e->layer_props);
	if (e->wcslib) { wcsfree(e->wcslib); free(e->wcslib); }
	if (e->icc_profile) cmsCloseProfile(e->icc_profile);
	/* note: callers free the array of entries; we only free contents. */
}

/* Build a per-layer snapshot inside a compound entry: pixels + pmask +
 * lmask + all metadata + props.  Returns 0 on success, 1 on failure
 * (in which case any partially-built fields are cleaned up). */
static int build_layer_entry(flis_layer_t *lay, historic_layer_entry_t *e) {
	memset(e, 0, sizeof(*e));
	if (!lay || !lay->fit) return 1;
	fits *fit = lay->fit;

	e->flis_layer_id = lay->item_id;
	e->rx     = fit->rx;
	e->ry     = fit->ry;
	e->nchans = fit->naxes[2];
	e->type   = fit->type;
	e->wcsdata = fit->keywords.wcsdata;
	int s = -1;
	e->wcslib = wcs_deepcopy(fit->keywords.wcslib, &s);
	e->focal_length = fit->keywords.focal_length;
	/* ICC profile is image-level (com.uniq); the per-layer entry no longer
	 * needs to snapshot it.  Captured by the lightweight icc-only undo
	 * flavour for ICC-altering operations. */
	e->icc_profile  = NULL;
	e->position_x   = lay->position_x;
	e->position_y   = lay->position_y;

	/* Props snapshot */
	e->layer_props = g_new0(flis_layer_props_t, 1);
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

	/* Pixel swap */
	if (undo_build_swapfile_named(fit, &e->filename)) {
		free_layer_entry(e);
		return 1;
	}
	/* Processing mask swap (NULL filename if fit had no mask) */
	if (undo_build_mask_swapfile_named(fit, &e->mask_filename)) {
		free_layer_entry(e);
		return 1;
	}
	if (fit->mask) e->mask_bitpix = fit->mask->bitpix;

	/* Layer mask swap */
	if (undo_build_lmask_swapfile(lay->lmask,
	                              &e->lmask_filename,
	                              &e->lmask_w, &e->lmask_h, &e->lmask_bitpix)) {
		free_layer_entry(e);
		return 1;
	}
	return 0;
}

/* Restore a layer from a compound-entry snapshot.  Returns 0 on success,
 * 1 on failure (partial restore possible). */
static int restore_layer_entry(flis_layer_t *lay, const historic_layer_entry_t *e) {
	if (!lay || !lay->fit || !e) return 1;
	fits *fit = lay->fit;

	/* Props-only entry: skip pixel/mask restore */
	if (e->props_only && e->layer_props) {
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
		return 0;
	}

	/* ICC profile no longer lives on the fits struct.  For the current
	 * image, profile restoration happens via the lightweight icc-only
	 * undo flavour.  Non-current fits never had a profile to restore. */
	fits_change_depth(fit, e->nchans);

	/* Adjust precision if needed */
	if (e->type == DATA_USHORT && fit->type != DATA_USHORT) {
		size_t ndata = (size_t)fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
		fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
	} else if (e->type == DATA_FLOAT && fit->type != DATA_FLOAT) {
		size_t ndata = (size_t)fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
		fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
	}

	/* Pixels */
	fit->rx = fit->naxes[0] = e->rx;
	fit->ry = fit->naxes[1] = e->ry;
	size_t n = (size_t)e->rx * e->ry;
	size_t pix_size = n * e->nchans * (e->type == DATA_USHORT ? sizeof(WORD) : sizeof(float));
	void *buf = read_swap_file(e->filename, pix_size);
	if (!buf) return 1;
	if (e->type == DATA_USHORT) {
		WORD *newdata = realloc(fit->data, pix_size);
		if (!newdata) { PRINT_ALLOC_ERR; free(buf); return 1; }
		fit->data = newdata;
		memcpy(fit->data, buf, pix_size);
		fit->pdata[RLAYER] = fit->data;
		if (fit->naxes[2] > 1) {
			fit->pdata[GLAYER] = fit->data + n;
			fit->pdata[BLAYER] = fit->data + n * 2;
		} else {
			fit->pdata[GLAYER] = fit->pdata[BLAYER] = fit->pdata[RLAYER];
		}
	} else {
		float *newdata = realloc(fit->fdata, pix_size);
		if (!newdata) { PRINT_ALLOC_ERR; free(buf); return 1; }
		fit->fdata = newdata;
		memcpy(fit->fdata, buf, pix_size);
		fit->fpdata[RLAYER] = fit->fdata;
		if (fit->naxes[2] > 1) {
			fit->fpdata[GLAYER] = fit->fdata + n;
			fit->fpdata[BLAYER] = fit->fdata + n * 2;
		} else {
			fit->fpdata[GLAYER] = fit->fpdata[BLAYER] = fit->fpdata[RLAYER];
		}
	}
	free(buf);

	/* WCS */
	memcpy(&fit->keywords.wcsdata, &e->wcsdata, sizeof(wcs_info));
	free_wcs(fit);
	if (e->wcslib) {
		int s = -1;
		fit->keywords.wcslib = wcs_deepcopy(e->wcslib, &s);
	} else {
		reset_wcsdata(fit);
	}
	fit->keywords.focal_length = e->focal_length;

	/* Processing mask */
	if (e->mask_filename) {
		size_t mask_elem;
		switch (e->mask_bitpix) {
			case 8:  mask_elem = sizeof(uint8_t);  break;
			case 16: mask_elem = sizeof(uint16_t); break;
			case 32: mask_elem = sizeof(float);    break;
			default: mask_elem = 0;
		}
		if (mask_elem) {
			size_t mask_size = n * mask_elem;
			void *mbuf = read_swap_file(e->mask_filename, mask_size);
			if (mbuf) {
				if (!fit->mask) fit->mask = calloc(1, sizeof(mask_t));
				fit->mask->bitpix = e->mask_bitpix;
				void *newmask = realloc(fit->mask->data, mask_size);
				if (newmask) {
					fit->mask->data = newmask;
					memcpy(fit->mask->data, mbuf, mask_size);
				}
				free(mbuf);
			}
		}
	} else if (fit->mask) {
		/* saved state had no mask — drop current */
		free(fit->mask->data);
		free(fit->mask);
		fit->mask = NULL;
	}

	/* Layer mask */
	if (e->lmask_filename) {
		size_t lmask_elem;
		switch (e->lmask_bitpix) {
			case 8:  lmask_elem = sizeof(uint8_t);  break;
			case 16: lmask_elem = sizeof(uint16_t); break;
			case 32: lmask_elem = sizeof(float);    break;
			default: lmask_elem = 0;
		}
		if (lmask_elem) {
			size_t lmask_size = e->lmask_w * e->lmask_h * lmask_elem;
			void *lbuf = read_swap_file(e->lmask_filename, lmask_size);
			if (lbuf) {
				layermask_t *lm = calloc(1, sizeof(layermask_t));
				lm->w = e->lmask_w;
				lm->h = e->lmask_h;
				lm->bitpix = e->lmask_bitpix;
				lm->data = lbuf;
				flis_layer_set_lmask(lay, lm);
			}
		}
	} else {
		flis_layer_remove_lmask(lay);
	}

	/* Props */
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

	full_stats_invalidation_from_fit(fit);
	return 0;
}

static void undo_free_item(historic *h) {
	/* fd-mode: closing the fd is all that is needed on both platforms — the
	 * file was marked delete-on-close at creation (POSIX unlink / Windows
	 * FileDispositionInfo) so the OS reclaims it here automatically. */
	if (h->fd >= 0)
		g_close(h->fd, NULL);
	if (h->mask_fd >= 0)
		g_close(h->mask_fd, NULL);
	/* filename-mode (FLIS): explicitly unlink + free the path strings. */
	if (h->filename) {
		g_unlink(h->filename);
		g_free(h->filename);
	}
	if (h->mask_filename) {
		g_unlink(h->mask_filename);
		g_free(h->mask_filename);
	}
	if (h->wcslib) {
		wcsfree(h->wcslib);
		free(h->wcslib);
	}
	if (h->icc_profile)
		cmsCloseProfile(h->icc_profile);
	/* FLIS-aware fields */
	g_free(h->layer_props);
	if (h->lmask_filename) {
		g_unlink(h->lmask_filename);
		g_free(h->lmask_filename);
	}
	if (h->multi_entries) {
		for (guint i = 0; i < h->n_multi_entries; i++)
			free_layer_entry(&h->multi_entries[i]);
		g_free(h->multi_entries);
	}
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
	/* Plain undo entries capture pixel state of the current image.  ICC
	 * state is restored separately via the icc-only undo flavour. */
	h->icc_profile = (fit == gfit) ? copyICCProfile(current_icc_profile()) : NULL;
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

/* Restore a "regular" pixel + mask state — the original fd-backed flavour
 * that backs undo_save_state.  Split out of undo_restore so the dispatcher
 * (added below) can call it for plain-FITS entries. */
static int undo_restore_plain(fits *fit, historic *hist) {
	/* Restore the image-level ICC profile if the entry captured one and
	 * we're restoring to the current image. */
	if (fit == gfit && hist->icc_profile) {
		current_image_set_icc_profile(copyICCProfile(hist->icc_profile));
		current_image_color_manage(TRUE);
	} else if (fit == gfit) {
		current_image_clear_icc_profile();
	}
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

/* Apply a flis_layer_props_t snapshot directly to a layer, bypassing the
 * mutating setters that would call flis_layer_touch_modified() and push a
 * new undo state.  Shared between the props-only and full-layer restore
 * branches. */
static void apply_layer_props(flis_layer_t *layer, const flis_layer_props_t *p) {
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
}

/* Restore an lmask from a swap file into a layer (or remove it if the
 * saved state had no mask).  Returns 0 on success, 1 on error. */
static int restore_lmask_for_layer(flis_layer_t *layer,
                                   const char *filename,
                                   size_t w, size_t h, guint8 bitpix) {
	if (!filename) {
		flis_layer_remove_lmask(layer);
		return 0;
	}
	size_t elem;
	switch (bitpix) {
		case 8:  elem = sizeof(uint8_t);  break;
		case 16: elem = sizeof(uint16_t); break;
		case 32: elem = sizeof(float);    break;
		default:
			siril_log_warning(_("Undo: invalid lmask bitpix %d\n"), (int)bitpix);
			return 1;
	}
	size_t total = w * h * elem;
	void *buf = read_swap_file(filename, total);
	if (!buf) return 1;
	layermask_t *lm = calloc(1, sizeof(layermask_t));
	if (!lm) { PRINT_ALLOC_ERR; free(buf); return 1; }
	lm->w = w; lm->h = h; lm->bitpix = bitpix; lm->data = buf;
	if (flis_layer_set_lmask(layer, lm)) {
		layermask_free(lm);
		return 1;
	}
	return 0;
}

/* Dispatch a historic entry to the right restore path based on its
 * flavour flags.  All FLIS entries route through here; plain-FITS entries
 * fall through to undo_restore_plain. */
static int undo_restore(fits *fit, historic *hist) {
	/* Compound multi-layer: restore each layer in turn. */
	if (hist->n_multi_entries > 0 && hist->multi_entries) {
		if (!is_current_image_flis()) return 1;
		for (guint k = 0; k < hist->n_multi_entries; k++) {
			const historic_layer_entry_t *e = &hist->multi_entries[k];
			flis_layer_t *lay = flis_layer_get_by_id(e->flis_layer_id);
			if (!lay) {
				siril_log_warning(_("Undo: target layer (id %d) no longer exists — skipping\n"),
				                  e->flis_layer_id);
				continue;
			}
			if (restore_layer_entry(lay, e)) {
				siril_log_error(_("Undo: failed to restore layer (id %d)\n"),
				                e->flis_layer_id);
				return 1;
			}
		}
		return 0;
	}

	/* Processing-mask-only: restore fit->mask, leave pixels alone. */
	if (hist->pmask_only) {
		return undo_get_mask_data(fit, hist);
	}

	/* ICC-profile-only: install the snapshot via the accessor so the
	 * gfit / FLIS-base mirrors stay consistent.  Then refresh the
	 * proofing transform so the display picks up the restored profile. */
	if (hist->icc_only) {
		current_image_set_icc_profile(hist->icc_profile
			? copyICCProfile(hist->icc_profile) : NULL);
		current_image_color_manage(hist->icc_was_managed);
		refresh_icc_transforms();
		return 0;
	}

	/* Full single-layer geometry state: pixels + pmask + lmask + props. */
	if (hist->full_layer) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *lay = flis_layer_get_by_id(hist->flis_layer_id);
		if (!lay) {
			siril_log_warning(_("Undo: target layer (id %d) no longer exists\n"),
			                  hist->flis_layer_id);
			return 1;
		}
		/* Pixels + pmask via plain restore (uses hist->fd path) */
		int rv = undo_restore_plain(fit, hist);
		if (rv) return rv;
		/* Layer mask */
		if (restore_lmask_for_layer(lay, hist->lmask_filename,
		                            hist->lmask_w, hist->lmask_h,
		                            hist->lmask_bitpix))
			return 1;
		/* Props */
		if (hist->layer_props) apply_layer_props(lay, hist->layer_props);
		return 0;
	}

	/* Property-only: no pixel swap, just apply props to the layer. */
	if (!hist->filename && hist->fd < 0 && hist->layer_props
	    && hist->flis_layer_id != FLIS_UNDO_LAYER_NONE) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *layer = flis_layer_get_by_id(hist->flis_layer_id);
		if (!layer) {
			siril_log_warning(_("Undo: target layer (id %d) no longer exists\n"),
			                  hist->flis_layer_id);
			return 1;
		}
		apply_layer_props(layer, hist->layer_props);
		return 0;
	}

	/* Atomic lmask-move: lmask_layer_id = source, lmask_dest_layer_id = dest. */
	if (!hist->filename && hist->fd < 0 && !hist->layer_props
	    && hist->lmask_layer_id != FLIS_UNDO_LAYER_NONE
	    && hist->lmask_dest_layer_id != FLIS_UNDO_LAYER_NONE) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *src  = flis_layer_get_by_id(hist->lmask_layer_id);
		flis_layer_t *dest = flis_layer_get_by_id(hist->lmask_dest_layer_id);
		if (!src || !dest) {
			siril_log_warning(_("Undo: layer for mask move no longer exists\n"));
			return 1;
		}
		/* dest currently holds the mask; move it back to src */
		flis_layer_move_lmask(dest, src);
		return 0;
	}

	/* Layer reorder: swap saved layer_order values back. */
	if (!hist->filename && hist->fd < 0 && !hist->layer_props
	    && hist->reorder_layer_a_id != FLIS_UNDO_LAYER_NONE) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *a = flis_layer_get_by_id(hist->reorder_layer_a_id);
		flis_layer_t *b = flis_layer_get_by_id(hist->reorder_layer_b_id);
		if (!a || !b) {
			siril_log_warning(_("Undo: layer for reorder no longer exists\n"));
			return 1;
		}
		a->layer_order = hist->reorder_layer_a_order;
		b->layer_order = hist->reorder_layer_b_order;
		flis_sort_layer_stack();
		return 0;
	}

	/* Lmask add/remove (no move, no reorder, no compound). */
	if (!hist->filename && hist->fd < 0 && !hist->layer_props
	    && hist->lmask_layer_id != FLIS_UNDO_LAYER_NONE) {
		if (!is_current_image_flis()) return 1;
		flis_layer_t *layer = flis_layer_get_by_id(hist->lmask_layer_id);
		if (!layer) {
			siril_log_warning(_("Undo: target layer (id %d) no longer exists\n"),
			                  hist->lmask_layer_id);
			return 1;
		}
		return restore_lmask_for_layer(layer, hist->lmask_filename,
		                               hist->lmask_w, hist->lmask_h,
		                               hist->lmask_bitpix);
	}

	/* Default: plain pixel + mask restore. */
	int rv = undo_restore_plain(fit, hist);

	/* If this state belongs to a FLIS layer, also restore the layer offset
	 * — geometry operations (crop, rotate) save position_x/y alongside the
	 * pixels and we need to put it back. */
	if (rv == 0 && hist->flis_layer_id != FLIS_UNDO_LAYER_NONE
	    && is_current_image_flis()) {
		flis_layer_t *lay = flis_layer_get_by_id(hist->flis_layer_id);
		if (lay) {
			lay->position_x = hist->flis_position_x;
			lay->position_y = hist->flis_position_y;
		}
	}
	return rv;
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
			/* redraw_mask_idle posts an idle - must be called outside any gfit lock */
			gui_iface.redraw_mask_idle();
			if (!com.pref.gui.mask_tints_vports) { // redraw() is called in redraw_mask_idle if this is TRUE
				/* No reader lock here: notify_gfit_data_modified() may call
				 * copy_roi_into_gfit() which acquires the writer lock —
				 * mirrors the contract redraw_mask_idle relies on. */
				notify_gfit_data_modified();
				gui_iface.redraw_image(REDRAW_ALL);
			}
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
			/* FLIS panel + canvas-properties dialog snapshot layer
			 * geometry on draw; an undo that restored layer positions
			 * or dims needs them to repaint. */
			gui_iface.flis_gui_update();
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
			/* redraw_mask_idle posts an idle - must be called outside any gfit lock */
			gui_iface.redraw_mask_idle();
			if (!com.pref.gui.mask_tints_vports) { // redraw() is called in redraw_mask_idle if this is TRUE
				/* No reader lock: notify_gfit_data_modified() may call
				 * copy_roi_into_gfit() which acquires the writer lock. */
				notify_gfit_data_modified();
				gui_iface.redraw_image(REDRAW_ALL);
			}
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
			/* See UNDO case — mirror the notification so the FLIS panel
			 * and canvas dialog repaint after a redo. */
			gui_iface.flis_gui_update();
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


/* =====================================================================
 * FLIS undo API — public save functions (stage 1.4)
 * =====================================================================
 *
 * All of these:
 *   • short-circuit (return 0) when com.script is set — headless runs save
 *     no undo state per the §C.1a principle.  This is in addition to the
 *     gating in generic_layer_worker (§1.5); duplicated here so callers
 *     outside the worker framework (e.g. opacity drag-end) are also safe.
 *   • short-circuit (return 1) when no single image is loaded.
 *   • discard the redo stack (new operation invalidates the redo branch).
 *   • g_list_prepend the new historic entry to com.undo_stack.
 *   • call gui_iface.update_menu_state() so the UI's undo button enables.
 *
 * The historic entry's flavour is encoded by which fields are populated:
 *   • props-only:  layer_props != NULL,  fd < 0,  filename == NULL
 *   • lmask:       lmask_layer_id != NONE
 *   • lmask-move:  lmask_layer_id != NONE && lmask_dest_layer_id != NONE
 *   • reorder:     reorder_layer_a_id != NONE
 *   • pmask-only:  pmask_only flag set
 *   • full-layer:  full_layer flag set
 *   • compound:    n_multi_entries > 0
 * The undo_restore dispatcher above keys off these to choose the right
 * restore path.
 * ===================================================================== */

/* Internal helper: initialise a freshly-allocated historic to the
 * "no flavour set" baseline.  Fields not touched by the caller stay zero,
 * so the dispatcher's first matching branch wins. */
static historic *flis_alloc_historic(const char *message) {
	historic *h = g_new0(historic, 1);
	h->fd = -1;
	h->mask_fd = -1;
	h->flis_layer_id        = FLIS_UNDO_LAYER_NONE;
	h->lmask_layer_id       = FLIS_UNDO_LAYER_NONE;
	h->lmask_dest_layer_id  = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_a_id   = FLIS_UNDO_LAYER_NONE;
	h->reorder_layer_b_id   = FLIS_UNDO_LAYER_NONE;
	snprintf(h->history, FLEN_VALUE, "%s", message ? message : "");
	return h;
}

/* Push an entry to the undo stack and discard any pending redo. */
static void flis_push_historic(historic *h) {
	g_list_free_full(com.redo_stack, (GDestroyNotify) undo_free_item);
	com.redo_stack = NULL;
	com.undo_stack = g_list_prepend(com.undo_stack, h);
	gui_iface.update_menu_state();
}

/* Format a printf-style varargs message into msg_var, identically to
 * undo_save_state.  The caller declares `char msg_var[FLEN_VALUE] = {0};`
 * before invoking — declaring inside a macro do-while block would create
 * a new scope and the buffer would be invisible to the caller. */
#define FLIS_FORMAT_MESSAGE_INTO(msg_var, fmt_var) do {                   \
	if ((fmt_var) != NULL) {                                              \
		va_list _args;                                                    \
		va_start(_args, (fmt_var));                                       \
		vsnprintf((msg_var), FLEN_VALUE, (fmt_var), _args);               \
		va_end(_args);                                                    \
	}                                                                     \
} while (0)

/* Internal: push a props-only entry for the given (item_id, props).  Used
 * by both undo_save_flis_layer_props and …_snapshot. */
static int undo_push_flis_layer_props(gint item_id,
                                      const flis_layer_props_t *props,
                                      const char *histo) {
	historic *h = flis_alloc_historic(histo);
	h->layer_props   = g_memdup2(props, sizeof(*props));
	h->flis_layer_id = item_id;
	flis_push_historic(h);
	return 0;
}

int undo_save_flis_layer_props(flis_layer_t *layer, const char *message, ...) {
	if (com.script) return 0;
	if (!layer || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	char msg_buf[FLEN_VALUE] = { 0 };
	FLIS_FORMAT_MESSAGE_INTO(msg_buf, message);

	flis_layer_props_t props = {
		.blend_mode   = layer->blend_mode,
		.opacity      = layer->opacity,
		.visible      = layer->visible,
		.locked       = layer->locked,
		.has_tint     = layer->has_tint,
		.tint         = layer->layer_tint,
		.lmask_active = layer->lmask_active,
		.position_x   = layer->position_x,
		.position_y   = layer->position_y,
	};
	g_strlcpy(props.name, layer->layer_name ? layer->layer_name : "", sizeof(props.name));
	return undo_push_flis_layer_props(layer->item_id, &props, msg_buf);
}

int undo_save_flis_layer_props_snapshot(gint item_id,
                                        const flis_layer_props_t *props,
                                        const char *message) {
	if (com.script) return 0;
	if (!is_current_image_flis() || !single_image_is_loaded() || !props)
		return 1;
	return undo_push_flis_layer_props(item_id, props, message);
}

int undo_save_flis_lmask(flis_layer_t *layer, const char *message) {
	if (com.script) return 0;
	if (!layer || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	historic *h = flis_alloc_historic(message);
	h->lmask_layer_id = layer->item_id;

	if (layer->lmask) {
		if (undo_build_lmask_swapfile(layer->lmask, &h->lmask_filename,
		                              &h->lmask_w, &h->lmask_h, &h->lmask_bitpix)) {
			undo_free_item(h);
			return 1;
		}
	}
	/* else: empty mask state — h->lmask_filename stays NULL, sentinel for
	 * "no mask was present at save time". */

	flis_push_historic(h);
	return 0;
}

int undo_save_flis_lmask_move(flis_layer_t *source, flis_layer_t *dest,
                              const char *message) {
	if (com.script) return 0;
	if (!source || !dest || !is_current_image_flis() || !single_image_is_loaded())
		return 1;
	historic *h = flis_alloc_historic(message);
	h->lmask_layer_id      = source->item_id;
	h->lmask_dest_layer_id = dest->item_id;
	flis_push_historic(h);
	return 0;
}

int undo_save_flis_layer_reorder(flis_layer_t *layer_a, flis_layer_t *layer_b,
                                 const char *message) {
	if (com.script) return 0;
	if (!layer_a || !layer_b || !is_current_image_flis() || !single_image_is_loaded())
		return 1;
	historic *h = flis_alloc_historic(message);
	h->reorder_layer_a_id    = layer_a->item_id;
	h->reorder_layer_a_order = layer_a->layer_order;
	h->reorder_layer_b_id    = layer_b->item_id;
	h->reorder_layer_b_order = layer_b->layer_order;
	flis_push_historic(h);
	return 0;
}

int undo_save_flis_multi_layer(GSList *layers, const char *message, ...) {
	if (com.script) return 0;
	if (!layers || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	char msg_buf[FLEN_VALUE] = { 0 };
	FLIS_FORMAT_MESSAGE_INTO(msg_buf, message);

	guint n = g_slist_length(layers);
	historic *h = flis_alloc_historic(msg_buf);
	h->multi_entries   = g_new0(historic_layer_entry_t, n);
	h->n_multi_entries = n;

	guint i = 0;
	for (GSList *l = layers; l; l = l->next, i++) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (build_layer_entry(lay, &h->multi_entries[i])) {
			/* On failure, clean up everything we've built so far */
			undo_free_item(h);
			return 1;
		}
	}
	flis_push_historic(h);
	return 0;
}

int undo_save_flis_multi_layer_props(GSList *layers, const char *message, ...) {
	if (com.script) return 0;
	if (!layers || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	char msg_buf[FLEN_VALUE] = { 0 };
	FLIS_FORMAT_MESSAGE_INTO(msg_buf, message);

	guint n = g_slist_length(layers);
	historic *h = flis_alloc_historic(msg_buf);
	h->multi_entries   = g_new0(historic_layer_entry_t, n);
	h->n_multi_entries = n;

	guint i = 0;
	for (GSList *l = layers; l; l = l->next, i++) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		historic_layer_entry_t *e = &h->multi_entries[i];
		e->flis_layer_id = lay->item_id;
		e->props_only    = TRUE;
		e->layer_props   = g_new0(flis_layer_props_t, 1);
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
		/* pixel/mask filenames remain NULL */
	}
	flis_push_historic(h);
	return 0;
}

/* Lightweight ICC-only undo entry.  Captures the current image's
 * profile + color_managed flag from com.uniq — no swap files, no
 * pixel data.  Used by the Assign / Convert / Remove ICC dialog
 * paths and by icc_auto_assign_or_convert.  Restore branch in
 * undo_restore re-installs the snapshot via the current_image_*
 * accessors. */
int undo_save_icc_state(const char *message, ...) {
	if (com.script) return 0;
	if (!single_image_is_loaded())
		return 1;
	char msg_buf[FLEN_VALUE] = { 0 };
	FLIS_FORMAT_MESSAGE_INTO(msg_buf, message);

	historic *h = flis_alloc_historic(msg_buf);
	h->icc_only        = TRUE;
	h->icc_was_managed = com.uniq ? com.uniq->color_managed : FALSE;
	h->icc_profile     = (com.uniq && com.uniq->icc_profile)
	                   ? copyICCProfile(com.uniq->icc_profile) : NULL;
	flis_push_historic(h);
	return 0;
}

int undo_save_processing_mask(fits *fit, const char *message, ...) {
	if (com.script) return 0;
	if (!fit || !single_image_is_loaded())
		return 1;

	char msg_buf[FLEN_VALUE] = { 0 };
	FLIS_FORMAT_MESSAGE_INTO(msg_buf, message);

	historic *h = flis_alloc_historic(msg_buf);
	h->pmask_only  = TRUE;
	h->rx = fit->rx;
	h->ry = fit->ry;
	h->mask_bitpix = (fit->mask && fit->mask->data) ? fit->mask->bitpix : 0;
	/* Build the mask swap file via the fd path (same as undo_save_state).
	 * Mask-only entries are short-lived and rare, so the fd cost is fine. */
	int mfd = undo_build_mask_swapfile(fit);
	if (mfd == -2) {  /* error (mfd == -1 means "no mask", which is fine) */
		undo_free_item(h);
		return 1;
	}
	h->mask_fd = mfd;
	flis_push_historic(h);
	return 0;
}

int undo_save_flis_layer_full(fits *fit_snapshot,
                              flis_layer_t *lay,
                              layermask_t *lmask_snapshot,
                              const flis_layer_props_t *props_snapshot,
                              const char *message, ...) {
	if (com.script) return 0;
	if (!fit_snapshot || !lay || !props_snapshot
	    || !is_current_image_flis() || !single_image_is_loaded())
		return 1;

	char msg_buf[FLEN_VALUE] = { 0 };
	FLIS_FORMAT_MESSAGE_INTO(msg_buf, message);

	historic *h = flis_alloc_historic(msg_buf);
	h->full_layer       = TRUE;
	h->flis_layer_id    = lay->item_id;
	h->flis_position_x  = lay->position_x;
	h->flis_position_y  = lay->position_y;
	h->rx = fit_snapshot->rx;
	h->ry = fit_snapshot->ry;
	h->nchans = fit_snapshot->naxes[2];
	h->type   = fit_snapshot->type;
	h->wcsdata = fit_snapshot->keywords.wcsdata;
	int s = -1;
	h->wcslib = wcs_deepcopy(fit_snapshot->keywords.wcslib, &s);
	h->focal_length = fit_snapshot->keywords.focal_length;
	/* fit_snapshot is a layer snapshot — never carries a profile.  ICC
	 * state for the FLIS lives on com.uniq and is preserved unchanged
	 * across geometry operations that this entry type backs. */
	h->icc_profile  = NULL;
	h->mask_bitpix  = (fit_snapshot->mask && fit_snapshot->mask->data)
	                  ? fit_snapshot->mask->bitpix : 0;

	/* Pixel + pmask: use the fd path so undo_restore_plain can read them. */
	int fd = undo_build_swapfile(fit_snapshot);
	if (fd < 0) { undo_free_item(h); return 1; }
	h->fd = fd;
	int mfd = undo_build_mask_swapfile(fit_snapshot);
	if (mfd == -2) { undo_free_item(h); return 1; }
	h->mask_fd = mfd;

	/* Layer mask snapshot (NULL filename means saved state had no mask) */
	if (lmask_snapshot) {
		if (undo_build_lmask_swapfile(lmask_snapshot, &h->lmask_filename,
		                              &h->lmask_w, &h->lmask_h, &h->lmask_bitpix)) {
			undo_free_item(h);
			return 1;
		}
	}

	/* Props snapshot */
	h->layer_props = g_memdup2(props_snapshot, sizeof(*props_snapshot));

	flis_push_historic(h);
	return 0;
}

/* =====================================================================
 * flis_undo_purge_layer
 * =====================================================================
 *
 * Walk both stacks and remove any entry that addresses the layer being
 * deleted, so we never try to restore pixels into a layer that no longer
 * exists.  Compound entries are handled specially: if exactly one
 * sub-entry matches the deleted layer, we shrink the entry's array; if
 * every sub-entry matches, we remove the whole entry.
 *
 * Call BEFORE the layer is freed so item_id is still valid.
 * ===================================================================== */

static gboolean entry_addresses_layer(const historic *h, gint id) {
	if (h->flis_layer_id == id) return TRUE;
	if (h->lmask_layer_id == id || h->lmask_dest_layer_id == id) return TRUE;
	if (h->reorder_layer_a_id == id || h->reorder_layer_b_id == id) return TRUE;
	if (h->n_multi_entries > 0 && h->multi_entries) {
		for (guint i = 0; i < h->n_multi_entries; i++)
			if (h->multi_entries[i].flis_layer_id == id) return TRUE;
	}
	return FALSE;
}

/* For a compound entry that addresses @id, remove the matching sub-entries
 * in place.  Returns TRUE if the compound is now empty (caller should
 * remove the whole historic), FALSE if it still has surviving entries. */
static gboolean compact_compound(historic *h, gint id) {
	if (h->n_multi_entries == 0 || !h->multi_entries) return FALSE;
	guint dst = 0;
	for (guint src = 0; src < h->n_multi_entries; src++) {
		historic_layer_entry_t *e = &h->multi_entries[src];
		if (e->flis_layer_id == id) {
			free_layer_entry(e);
		} else {
			if (dst != src) h->multi_entries[dst] = *e;
			dst++;
		}
	}
	h->n_multi_entries = dst;
	return (dst == 0);
}

static void purge_layer_from_stack(GList **stack, gint id) {
	GList *cur = *stack;
	while (cur) {
		GList *next = cur->next;
		historic *h = (historic *)cur->data;
		if (entry_addresses_layer(h, id)) {
			gboolean remove_whole = TRUE;
			if (h->n_multi_entries > 0)
				remove_whole = compact_compound(h, id);
			if (remove_whole) {
				*stack = g_list_delete_link(*stack, cur);
				undo_free_item(h);
			}
		}
		cur = next;
	}
}

void flis_undo_purge_layer(gint item_id) {
	if (item_id == FLIS_UNDO_LAYER_NONE) return;
	purge_layer_from_stack(&com.undo_stack, item_id);
	purge_layer_from_stack(&com.redo_stack, item_id);
	gui_iface.update_menu_state();
}
