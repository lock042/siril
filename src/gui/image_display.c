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

#include <math.h>
#include <ctype.h>

#include "git-version.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "algos/astrometry_solver.h"
#include "algos/colors.h"
#include "algos/ccd-inspector.h"
#include "gui/ccd-inspector.h"
#include "algos/background_extraction.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "algos/sorting.h"
#include "io/annotation_catalogues.h"
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_shift.h"
#include "gui/mpp_ap_editor.h"
#include "gui/mpp_shift_viewer.h"
#include "filters/mtf.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/image_interactions.h"
#include "gui/registration_preview.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/siril_preview.h"
#include "gui/user_polygons.h"
#include "livestacking/livestacking.h"
#include "histogram.h"
#include "registration/matching/degtorad.h"
#include "registration/registration.h"
#include "opencv/opencv.h"

#include <wcslib.h>
#include <wcsfix.h>

#include "image_display.h"

/* is gfit->icc_profile identical to the monitor profile, if so we can avoid the
 * transform */
static cmsBool identical = FALSE;

#define ANGLE_TOP 315. * DEGTORAD
#define ANGLE_BOT 45. * DEGTORAD

/* remap index data, an index for each layer */
static float last_pente;
static display_mode last_mode;

/* STF (auto-stretch) data */
static gboolean stf_computed = FALSE; // Flag to know if STF parameters are available
static struct mtf_params stf[3];

/* widgets for draw_reg_data*/
GtkComboBox *seqcombo;
GtkToggleButton *drawframe;
static GtkWidget *rotation_dlg = NULL;

/* widgets for cut tool*/
static GtkWidget *cut_dialog = NULL, *cut_cdialog = NULL, *cut_sdialog = NULL;
static GtkToggleButton *tri_cut_toggle = NULL;
static GtkSpinButton *tri_cut_spin_step = NULL;

static GtkApplicationWindow *imgdisp_app_win = NULL;
static GtkCheckMenuItem *imgdisp_autohd_item = NULL;
static GtkWidget *imgdisp_drawing_rgb = NULL;
static GtkWidget *imgdisp_drawing_r = NULL;

static void image_display_init_statics(void) {
	if (imgdisp_app_win) return;
	imgdisp_app_win = GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	imgdisp_autohd_item = GTK_CHECK_MENU_ITEM(gtk_builder_get_object(gui.builder, "autohd_item"));
	imgdisp_drawing_rgb = GTK_WIDGET(gtk_builder_get_object(gui.builder, "drawingareargb"));
	imgdisp_drawing_r = GTK_WIDGET(gtk_builder_get_object(gui.builder, "drawingarear"));
	seqcombo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "seqlist_dialog_combo"));
	drawframe = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "drawframe_check"));
	rotation_dlg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rotation_dialog"));
	cut_dialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cut_dialog"));
	cut_cdialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cut_coords_dialog"));
	cut_sdialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cut_spectroscopy_dialog"));
	tri_cut_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "cut_tri_cut"));
	tri_cut_spin_step = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "cut_tricut_step"));
}

static void invalidate_image_render_cache(int vport);

static int allocate_full_surface(struct image_view *view) {
	g_mutex_lock(&gui.cairo_mutex);

	/* Cairo surfaces are limited to 32767 px per side (int16_t internally).
	 * For images that exceed this limit, scale down to at most 4096 px per
	 * side so the display buffer stays within reasonable memory bounds.
	 * (Step 2 - GtkGLArea - will provide full-resolution tile rendering.) */
#define SIRIL_MAX_SURFACE_DIM 4096
	double scale = 1.0;
	if (gfit->rx > 32767 || gfit->ry > 32767)
		scale = MIN((double)SIRIL_MAX_SURFACE_DIM / gfit->rx,
		            (double)SIRIL_MAX_SURFACE_DIM / gfit->ry);
#undef SIRIL_MAX_SURFACE_DIM
	gui.surface_scale = scale;

	int sw = (int)(gfit->rx * scale);
	int sh = (int)(gfit->ry * scale);
	if (sw < 1) sw = 1;
	if (sh < 1) sh = 1;

	int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, sw);

	// allocate the image surface if not already done or the same size
	if (stride != view->full_surface_stride
				|| sh != view->full_surface_height
				|| !view->full_surface || !view->buf) {
		siril_log_debug("display buffers and full_surface (re-)allocation %p\n", view);

		guchar *tmp = realloc(view->buf, (size_t)stride * sh * sizeof(guchar));
		if (!tmp) {
			PRINT_ALLOC_ERR;
			free(view->buf);
			view->buf = NULL;
			if (view->full_surface) {
				cairo_surface_destroy(view->full_surface);
				view->full_surface = NULL;
			}
			view->full_surface_stride = 0;
			view->full_surface_height = 0;
			g_mutex_unlock(&gui.cairo_mutex);
			return 1;
		}
		view->buf = tmp;
		// Only update dimension fields once the allocation succeeded
		view->full_surface_stride = stride;
		view->full_surface_height = sh;

		if (view->full_surface)
			cairo_surface_destroy(view->full_surface);
		view->full_surface = cairo_image_surface_create_for_data(view->buf,
					CAIRO_FORMAT_RGB24, sw, sh, stride);
		if (cairo_surface_status(view->full_surface) != CAIRO_STATUS_SUCCESS) {
			siril_log_debug("Error creating the cairo image full_surface for the RGB image\n");
			cairo_surface_destroy(view->full_surface);
			view->full_surface = NULL;
			free(view->buf);
			view->buf = NULL;
			view->full_surface_stride = 0;
			view->full_surface_height = 0;
			g_mutex_unlock(&gui.cairo_mutex);
			return 1;
		}
	}
	g_mutex_unlock(&gui.cairo_mutex);
	return 0;
}

void check_gfit_profile_identical_to_monitor() {
	if (!com.headless && gfit->icc_profile && gfit->color_managed)
		identical = profiles_identical(gfit->icc_profile, com.gui_icc.monitor);
	siril_log_debug("gfit profile identical to monitor profile: %d\n", identical);
}

static void remaprgb(void) {
	guint32 *dst;
	const guint32 *bufr, *bufg, *bufb;
	gint i;
	int nbdata;

	siril_log_debug("remaprgb\n");
	if (!isrgb(gfit))
		return;

	struct image_view *rgbview = &gui.view[RGB_VPORT];
	if (allocate_full_surface(rgbview))
		return;

	// Source pointers are captured inside the lock to avoid a race with a
	// concurrent allocate_full_surface() realloc between its internal unlock
	// and the lock acquisition below.
	g_mutex_lock(&gui.cairo_mutex);
	bufr = (const guint32*) gui.view[RED_VPORT].buf;
	bufg = (const guint32*) gui.view[GREEN_VPORT].buf;
	bufb = (const guint32*) gui.view[BLUE_VPORT].buf;
	if (bufr == NULL || bufg == NULL || bufb == NULL) {
		siril_log_debug("remaprgb: gray buffers not allocated for display\n");
		g_mutex_unlock(&gui.cairo_mutex);
		return;
	}
	dst = (guint32*) rgbview->buf;	// index is j
	/* Channel bufs are at surface dimensions (may be smaller than gfit when image
	 * exceeds Cairo's 32767-px limit); use the actual surface size here. */
	nbdata = (rgbview->full_surface_stride / 4) * rgbview->full_surface_height;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (i = 0; i < nbdata; ++i) {
		dst[i] = (bufr[i] & 0xFF0000) | (bufg[i] & 0xFF00) | (bufb[i] & 0xFF);
	}

	// flush to ensure all writing to the image was done and redraw the surface
	cairo_surface_flush(rgbview->full_surface);
	cairo_surface_mark_dirty(rgbview->full_surface);
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; must be called after our unlock
	invalidate_image_render_cache(RGB_VPORT);
}

void allocate_hd_remap_indices() {
	gui.hd_remap_max = 1 << (guint) com.pref.hd_bitdepth;
	for (unsigned i=0; i < 3; i++) {
		if (gui.hd_remap_index[i] != NULL)
			free(gui.hd_remap_index[i]);
		gui.hd_remap_index[i] = (BYTE*) calloc(gui.hd_remap_max + 1, sizeof(BYTE));
		if (gui.hd_remap_index[i] == NULL) {
			siril_log_error(_("Error: memory allocaton failure when instantiating HD LUTs. Reverting to standard 16 bit LUTs.\n"));
			gui.use_hd_remap = FALSE;
			image_display_init_statics();
			gtk_check_menu_item_set_active(imgdisp_autohd_item, FALSE);
			hd_remap_indices_cleanup();
			return;
		}
	}
}

void hd_remap_indices_cleanup() {
	for (unsigned i=0 ; i < 3; i++) {
		if (gui.hd_remap_index[i] != NULL) {
			free(gui.hd_remap_index[i]);
			gui.hd_remap_index[i] = NULL;
		}
	}
}

static int make_index_for_current_display(int vport);

static int make_hd_index_for_current_display(int vport);

static int make_index_for_rainbow(BYTE index[][3]);

static void remap_mask(mask_t *mask) {
	siril_log_debug("mask remap\n");

	int vport = MASK_VPORT;
	struct image_view *view = &gui.view[vport];
	if (allocate_full_surface(view))
		return;

	g_mutex_lock(&gui.cairo_mutex);
	// Cairo setup
	unsigned char *dst_base = cairo_image_surface_get_data(view->full_surface);
	int dst_stride = cairo_image_surface_get_stride(view->full_surface);

	// Mask setup
	void *src_data = mask->data;
	guint rx = gfit->rx;
	guint ry = gfit->ry;
	int bitpix = mask->bitpix;

	/* When the image exceeds Cairo's 32767-px limit, the surface is smaller than
	 * gfit; nearest-neighbour downsample into the surface dimensions. */
	const gboolean downscaled = (gui.surface_scale < 1.0);
	const int sh = view->full_surface_height;
	const int sw = (int)(rx * gui.surface_scale);
	const double inv_scale = downscaled ? 1.0 / gui.surface_scale : 1.0;

	if (downscaled) {
		switch (bitpix) {
			case 8: {
				uint8_t *s8 = (uint8_t *)src_data;
				#ifdef _OPENMP
				#pragma omp parallel for num_threads(com.max_thread) schedule(static)
				#endif
				for (int dy = 0; dy < sh; dy++) {
					int sy = (int)(dy * inv_scale);
					uint8_t *src_row = s8 + (guint)sy * rx;
					uint32_t *dst_row = (uint32_t *)(dst_base + (sh - 1 - dy) * dst_stride);
					for (int dx = 0; dx < sw; dx++) {
						int sx = (int)(dx * inv_scale);
						uint8_t val = src_row[sx];
						dst_row[dx] = (val << 16) | (val << 8) | val;
					}
				}
				break;
			}
			case 16: {
				uint16_t *s16 = (uint16_t *)src_data;
				#ifdef _OPENMP
				#pragma omp parallel for num_threads(com.max_thread) schedule(static)
				#endif
				for (int dy = 0; dy < sh; dy++) {
					int sy = (int)(dy * inv_scale);
					uint16_t *src_row = s16 + (guint)sy * rx;
					uint32_t *dst_row = (uint32_t *)(dst_base + (sh - 1 - dy) * dst_stride);
					for (int dx = 0; dx < sw; dx++) {
						int sx = (int)(dx * inv_scale);
						uint32_t val = ((uint32_t)src_row[sx] * 255 + 32895) >> 16;
						dst_row[dx] = (val << 16) | (val << 8) | val;
					}
				}
				break;
			}
			case 32: {
				float *sf = (float *)src_data;
				#ifdef _OPENMP
				#pragma omp parallel for num_threads(com.max_thread) schedule(static)
				#endif
				for (int dy = 0; dy < sh; dy++) {
					int sy = (int)(dy * inv_scale);
					float *src_row = sf + (guint)sy * rx;
					uint32_t *dst_row = (uint32_t *)(dst_base + (sh - 1 - dy) * dst_stride);
					for (int dx = 0; dx < sw; dx++) {
						int sx = (int)(dx * inv_scale);
						uint8_t val = roundf_to_BYTE(src_row[sx] * UCHAR_MAX_SINGLE);
						dst_row[dx] = (val << 16) | (val << 8) | val;
					}
				}
				break;
			}
			default:
				siril_log_debug("Error: invalid mask bitpix\n");
				break;
		}
	} else {
	// We switch ONCE. This requires duplicating the loop code,
	// but ensures the tightest possible inner loops for the compiler to vectorize.
	switch (bitpix) {
		case 8: {
			uint8_t *s8 = (uint8_t *)src_data;
			#ifdef _OPENMP
			#pragma omp parallel for num_threads(com.max_thread) schedule(static)
			#endif
			for (guint y = 0; y < ry; y++) {
				// Pre-calculate pointers for this row
				uint8_t *src_row = s8 + (y * rx);
				// Vertical flip logic
				uint32_t *dst_row = (uint32_t *)(dst_base + ((ry - 1 - y) * dst_stride));

				for (guint x = 0; x < rx; x++) {
					uint8_t val = src_row[x];
					// Alpha is unused/ignored in Cairo RGB24 (upper 8 bits)
					dst_row[x] = (val << 16) | (val << 8) | val;
				}
			}
			break;
		}

		case 16: {
			uint16_t *s16 = (uint16_t *)src_data;
			#ifdef _OPENMP
			#pragma omp parallel for num_threads(com.max_thread) schedule(static)
			#endif
			for (guint y = 0; y < ry; y++) {
				uint16_t *src_row = s16 + (y * rx);
				uint32_t *dst_row = (uint32_t *)(dst_base + ((ry - 1 - y) * dst_stride));

				for (guint x = 0; x < rx; x++) {
					uint32_t val = ((uint32_t)src_row[x] * 255 + 32895) >> 16;
					dst_row[x] = (val << 16) | (val << 8) | val;
				}
			}
			break;
		}

		case 32: {
			float *sf = (float *)src_data;

			#ifdef _OPENMP
			#pragma omp parallel for num_threads(com.max_thread) schedule(static)
			#endif
			for (guint y = 0; y < ry; y++) {
				float *src_row = sf + (y * rx);

				// Destination pointer calculation (Vertical Flip)
				uint32_t *dst_row = (uint32_t *)(dst_base + ((ry - 1 - y) * dst_stride));

				for (guint x = 0; x < rx; x++) {
					// Scale
					float v = src_row[x] * UCHAR_MAX_SINGLE;

					// Round & Clamp
					uint8_t val = roundf_to_BYTE(v);

					// Pack to ARGB (Alpha unused)
					dst_row[x] = (val << 16) | (val << 8) | val;
				}
			}
			break;
		}

		default: {
			siril_log_debug("Error: invalid mask bitpix\n");
			break;
		}
	}
	} // end !downscaled

	cairo_surface_flush(view->full_surface);
	cairo_surface_mark_dirty(view->full_surface);
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; must be called after our unlock
	invalidate_image_render_cache(vport);
	test_and_allocate_reference_image(vport);
}

/*
 * Idle wrapper for set_viewer_mode_widgets_sensitive().
 *
 * remap() and remap_all_vports() may be called from a non-GTK thread
 * (via notify_gfit_data_modified()).  Widget-sensitivity changes must only
 * happen on the GTK main thread, so we dispatch them as idle callbacks.
 * The gint value (0 or 1) is passed via GINT_TO_POINTER / GPOINTER_TO_INT.
 */
static gboolean viewer_mode_sensitive_idle(gpointer data) {
	set_viewer_mode_widgets_sensitive(GPOINTER_TO_INT(data));
	return G_SOURCE_REMOVE;
}

// remapping one vport at a time is used for DISPLAY_STF and DISPLAY_HISTEQ
static void remap(int vport) {
	// This function maps fit data with a linear LUT between lo and hi levels
	// to the buffer to be displayed; display only is modified
	BYTE *dst, *index, rainbow_index[UCHAR_MAX + 1][3];
	WORD *src;
	float *fsrc;
	gboolean inverted;
	siril_log_debug("HISTEQ / STF remap %d\n", vport);
	if (vport == RGB_VPORT) {
		remaprgb();
		return;
	}
	if (gfit->type == DATA_UNSUPPORTED) {
		siril_log_debug("data is not loaded yet\n");
		return;
	}

	struct image_view *view = &gui.view[vport];
	if (allocate_full_surface(view))
		return;

	/* Cache the two GAction pointers on first call.
	 * g_action_map_lookup_action() is not thread-safe; caching here (initialised
	 * from the GTK main thread on first image display) avoids calling it from
	 * worker threads.  g_action_get_state() on a GSimpleAction is thread-safe
	 * (atomic pointer read) and is therefore safe to call from any thread. */
	static GAction *action_neg_cached = NULL;
	static GAction *action_color_cached = NULL;
	if (action_neg_cached == NULL) {
		image_display_init_statics();
		action_neg_cached   = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
		action_color_cached = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "color-map");
	}
	GVariant *neg_state = g_action_get_state(action_neg_cached);
	inverted = g_variant_get_boolean(neg_state);
	g_variant_unref(neg_state);
	neg_state = NULL;

	if (gui.rendering_mode == HISTEQ_DISPLAY) {
		double hist_sum, nb_pixels;
		size_t i, hist_nb_bins;
		gsl_histogram *histo;
		compute_histo_for_fit(gfit);
		histo = com.layers_hist[vport];
		hist_nb_bins = gsl_histogram_bins(histo);
		nb_pixels = (double)(gfit->rx * gfit->ry);
		// build the remap_index
		index = gui.remap_index[0];
		index[0] = 0;
		hist_sum = gsl_histogram_get(histo, 0);
		for (i = 1; i < hist_nb_bins; i++) {
			hist_sum += gsl_histogram_get(histo, i);
			index[i] = round_to_BYTE((hist_sum / nb_pixels) * UCHAR_MAX_DOUBLE);
		}

		last_mode = gui.rendering_mode;
		histo = com.layers_hist[vport];
		/* Widget-sensitivity changes must happen on the GTK main thread. */
		siril_add_idle(viewer_mode_sensitive_idle, GINT_TO_POINTER(FALSE));
	} else {
		if (gui.rendering_mode == STF_DISPLAY && !stf_computed) {
			if (gui.unlink_channels)
				find_unlinked_midtones_balance_default(gfit, stf);
			else find_linked_midtones_balance_default(gfit, stf);
			stf_computed = TRUE;
		}
		if (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit->type == DATA_FLOAT) {
			make_hd_index_for_current_display(vport);
		}
		else
			make_index_for_current_display(vport);
		siril_add_idle(viewer_mode_sensitive_idle,
		               GINT_TO_POINTER(gui.rendering_mode != STF_DISPLAY));
	}

	src = gfit->pdata[vport];
	fsrc = gfit->fpdata[vport];

	g_mutex_lock(&gui.cairo_mutex);
	dst = view->buf;

	GVariant *rainbow_state = g_action_get_state(action_color_cached);
	color_map color = g_variant_get_boolean(rainbow_state);
	g_variant_unref(rainbow_state);
	rainbow_state = NULL;

	if (color == RAINBOW_COLOR)
		make_index_for_rainbow(rainbow_index);
	int target_index = gui.rendering_mode == STF_DISPLAY && gui.unlink_channels ? vport : 0;

	gboolean hd_mode = (gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gfit->type == DATA_FLOAT);
	if (hd_mode) {
		index = gui.hd_remap_index[target_index];
	}
	else
		index = gui.remap_index[target_index];

	// Check if mask overlay is active
	gboolean apply_mask = (gfit->mask != NULL && gfit->mask_active && com.pref.gui.mask_tints_vports);
	uint8_t *mask_u8 = NULL;
	uint16_t *mask_u16 = NULL;
	float *mask_f32 = NULL;
	int mask_bitpix = 0;

	if (apply_mask) {
		mask_bitpix = gfit->mask->bitpix;
		switch (mask_bitpix) {
			case 8:
				mask_u8 = (uint8_t *)gfit->mask->data;
				break;
			case 16:
				mask_u16 = (uint16_t *)gfit->mask->data;
				break;
			case 32:
				mask_f32 = (float *)gfit->mask->data;
				break;
			default:
				apply_mask = FALSE; // Invalid bitpix, disable mask
				break;
		}
	}

	const guint width = gfit->rx;
	const guint height = gfit->ry;

/* Macro to compute a pixel from source index src_i and write it to dst+dst_i */
#define REMAP_WRITE_PIXEL(src_i, dst_i) \
	do { \
		BYTE dpv; \
		if (gfit->type == DATA_USHORT) { \
			const WORD sv = src[(src_i)]; \
			dpv = hd_mode ? index[sv * gui.hd_remap_max / USHRT_MAX] : index[sv]; \
		} else { \
			dpv = hd_mode ? index[float_to_max_range(fsrc[(src_i)], gui.hd_remap_max)] \
			              : index[roundf_to_WORD(fsrc[(src_i)] * USHRT_MAX_SINGLE)]; \
		} \
		if (inverted) dpv = UCHAR_MAX - dpv; \
		if (apply_mask) { \
			BYTE mv; \
			if (mask_bitpix == 8) mv = mask_u8[(src_i)]; \
			else if (mask_bitpix == 16) mv = mask_u16[(src_i)] >> 8; \
			else mv = (BYTE)(mask_f32[(src_i)] * UCHAR_MAX); \
			const uint16_t imv = UCHAR_MAX - mv; \
			if (color == RAINBOW_COLOR) { \
				const uint16_t rs = imv + rainbow_index[dpv][0]; \
				*(guint32*)(dst + (dst_i)) = \
					((rs > UCHAR_MAX ? UCHAR_MAX : (BYTE)rs) << 16) | \
					((rainbow_index[dpv][1] * mv) >> 8) << 8 | \
					((rainbow_index[dpv][2] * mv) >> 8); \
			} else { \
				const uint16_t rs = imv + dpv; \
				*(guint32*)(dst + (dst_i)) = \
					((rs > UCHAR_MAX ? UCHAR_MAX : (BYTE)rs) << 16) | \
					((dpv * mv) >> 8) << 8 | ((dpv * mv) >> 8); \
			} \
		} else { \
			if (color == RAINBOW_COLOR) \
				*(guint32*)(dst + (dst_i)) = rainbow_index[dpv][0] << 16 | \
				                             rainbow_index[dpv][1] << 8 | \
				                             rainbow_index[dpv][2]; \
			else \
				*(guint32*)(dst + (dst_i)) = dpv << 16 | dpv << 8 | dpv; \
		} \
	} while (0)

	if (gui.surface_scale < 1.0) {
		/* Downscaled path: buf/full_surface smaller than gfit; nearest-neighbour. */
		const int sw = (int)(width * gui.surface_scale);
		const int sh = view->full_surface_height;
		const int dstride = view->full_surface_stride;
		const double inv_scale = 1.0 / gui.surface_scale;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (int dy = 0; dy < sh; dy++) {
			const int sy = (int)(dy * inv_scale);
			const guint src_row_start = (guint)sy * width;
			const guint dst_row_start = (guint)(sh - 1 - dy) * (guint)dstride;
			for (int dx = 0; dx < sw; dx++) {
				const int sx = (int)(dx * inv_scale);
				REMAP_WRITE_PIXEL(src_row_start + (guint)sx, dst_row_start + (guint)dx * 4);
			}
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (guint y = 0; y < height; y++) {
			const guint src_row_start = y * width;
			const guint dst_row_start = (height - 1 - y) * width * 4;
			for (guint x = 0; x < width; ++x) {
				REMAP_WRITE_PIXEL(src_row_start + x, dst_row_start + x * 4);
			}
		}
	}
#undef REMAP_WRITE_PIXEL

	// flush to ensure all writing to the image was done and redraw the surface
	cairo_surface_flush(view->full_surface);
	cairo_surface_mark_dirty(view->full_surface);
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; must be called after our unlock
	invalidate_image_render_cache(vport);
	test_and_allocate_reference_image(vport);
}

static void remap_all_vports() {
	gboolean inverted;
	/* Snapshot gui.hi/gui.lo: init_layers_hi_and_lo_values() may write them
	 * from the worker thread concurrently while this path runs from
	 * remap_all() called directly on the worker (via notify_gfit_data_modified). */
	g_mutex_lock(&com.mutex);
	WORD remap_hi = gui.hi;
	WORD remap_lo = gui.lo;
	g_mutex_unlock(&com.mutex);

	/* Cache GAction pointers on first call — same reasoning as in remap() above. */
	static GAction *action_neg_cached = NULL;
	static GAction *action_color_cached = NULL;
	if (action_neg_cached == NULL) {
		image_display_init_statics();
		action_neg_cached   = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
		action_color_cached = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "color-map");
	}
	struct image_view *view[3] = { &gui.view[0], &gui.view[1], &gui.view[2] };
	GVariant *state_neg = g_action_get_state(action_neg_cached);
	inverted = g_variant_get_boolean(state_neg);
	g_variant_unref(state_neg);
	state_neg = NULL;
	// We are now dealing with a 3-channel image

	// Check if we need a rainbow color map
	BYTE rainbow_index[UCHAR_MAX + 1][3];
	GVariant* rainbow_state = g_action_get_state(action_color_cached);
	color_map color = g_variant_get_boolean(rainbow_state);
	g_variant_unref(rainbow_state);
	rainbow_state = NULL;
	if (color == RAINBOW_COLOR)
		make_index_for_rainbow(rainbow_index);

	// Check if mask overlay is active
	gboolean apply_mask = (gfit->mask != NULL && gfit->mask_active && com.pref.gui.mask_tints_vports);
	uint8_t *mask_u8 = NULL;
	uint16_t *mask_u16 = NULL;
	float *mask_f32 = NULL;
	int mask_bitpix = 0;

	if (apply_mask) {
		mask_bitpix = gfit->mask->bitpix;
		switch (mask_bitpix) {
			case 8:
				mask_u8 = (uint8_t *)gfit->mask->data;
				break;
			case 16:
				mask_u16 = (uint16_t *)gfit->mask->data;
				break;
			case 32:
				mask_f32 = (float *)gfit->mask->data;
				break;
			default:
				apply_mask = FALSE; // Invalid bitpix, disable mask
				break;
		}
	}

	// This function maps fit data with a linear LUT between lo and hi levels
	// to the buffer to be displayed; display only is modified
	BYTE *dst[3], *index[3];
	WORD *src[3];
	float *fsrc[3];

	if (gfit->type == DATA_UNSUPPORTED) {
		siril_log_debug("data is not loaded yet\n");
		return;
	}

	lock_display_transform();
	if (gfit->color_managed) {
		// Set the transform in case it is missing
		if (!com.gui_icc.proofing_transform) {
			com.gui_icc.proofing_transform = initialize_proofing_transform();
			com.gui_icc.profile_changed = TRUE;
		}
		if (com.gui_icc.profile_changed) {
			com.gui_icc.same_primaries = same_primaries(gfit->icc_profile, com.gui_icc.monitor, (com.gui_icc.soft_proof && com.pref.icc.soft_proofing_profile_active) ? com.gui_icc.soft_proof : NULL);
//			com.gui_icc.same_primaries = FALSE;
			check_gfit_profile_identical_to_monitor();
			// Calling color_manage() like this updates the color management button tooltip
			color_manage(gfit, gfit->color_managed);
			if (is_preview_active())
				copy_gfit_icc_to_backup();
		}
	}

	make_index_for_current_display(0);
	index[0] = gui.remap_index[0];
	if (gfit->color_managed) {
		for (int i = 1 ; i < 3 ; i++) {
			make_index_for_current_display(i);
			index[i] = gui.remap_index[i];
		}
	}
	unlock_display_transform();

	com.gui_icc.profile_changed = FALSE;
	/* Widget-sensitivity changes must happen on the GTK main thread. */
	siril_add_idle(viewer_mode_sensitive_idle, GINT_TO_POINTER(TRUE));

	last_mode = gui.rendering_mode;

	// Allocate surfaces and capture src/fsrc outside the lock.
	// dst[] pointers are captured below, inside the lock, to prevent a stale-
	// pointer race: a concurrent allocate_full_surface() could realloc view->buf
	// between its internal unlock and the cairo_mutex acquisition below.
	for (int i = 0 ; i < 3 ; i++) {
		src[i] = gfit->pdata[i];
		fsrc[i] = gfit->fpdata[i];
		if (allocate_full_surface(view[i]))
			return;
	}
	g_mutex_lock(&gui.cairo_mutex);
	for (int i = 0 ; i < 3 ; i++)
		dst[i] = view[i]->buf;

	int norm = (int) get_normalized_value(gfit);
	const guint width = gfit->rx;
	const guint height = gfit->ry;

	/* When image exceeds Cairo's 32767-px limit, the surface is downscaled. */
	const gboolean downscaled = (gui.surface_scale < 1.0);
	const int loop_max = downscaled ? view[0]->full_surface_height : (int)height;
	const int sw = downscaled ? (int)(width * gui.surface_scale) : (int)width;
	const int dstride = downscaled ? view[0]->full_surface_stride : (int)(width * 4);
	const double inv_scale = downscaled ? 1.0 / gui.surface_scale : 1.0;

	// Shared flag for per-thread allocation failures; checked after the parallel block.
	gboolean alloc_error = FALSE;

	{
		siril_log_debug((com.gui_icc.proofing_transform && !identical && (!com.gui_icc.same_primaries || com.gui_icc.profile_changed)) ? "Non-identical primaries: doing expensive color transform\n" : "");
		const gboolean do_transform = (com.gui_icc.proofing_transform && !identical && (!com.gui_icc.same_primaries || com.gui_icc.profile_changed));

		if (do_transform)
			lock_display_transform();

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) shared(alloc_error)
{
		// Thread-local buffers - allocated once per thread
		WORD *pixelbuf = malloc(width * 3 * sizeof(WORD));
		WORD *linebuf[3] = { pixelbuf, (pixelbuf + width), (pixelbuf + 2 * width) };
		BYTE *pixelbuf_byte = malloc(width * 3);
		BYTE *linebuf_byte[3] = { pixelbuf_byte, (pixelbuf_byte + width), (pixelbuf_byte + 2 * width) };
		BYTE *mask_row = apply_mask ? malloc(width * sizeof(BYTE)) : NULL;

		if (!pixelbuf || !pixelbuf_byte || (apply_mask && !mask_row)) {
			PRINT_ALLOC_ERR;
#pragma omp atomic write
			alloc_error = TRUE;
		} else {
#pragma omp for schedule(static)
#else
		// Single-threaded: allocate once
		WORD *pixelbuf = malloc(width * 3 * sizeof(WORD));
		WORD *linebuf[3] = { pixelbuf, (pixelbuf + width), (pixelbuf + 2 * width) };
		BYTE *pixelbuf_byte = malloc(width * 3);
		BYTE *linebuf_byte[3] = { pixelbuf_byte, (pixelbuf_byte + width), (pixelbuf_byte + 2 * width) };
		BYTE *mask_row = apply_mask ? malloc(width * sizeof(BYTE)) : NULL;

		if (!pixelbuf || !pixelbuf_byte || (apply_mask && !mask_row)) {
			PRINT_ALLOC_ERR;
			free(pixelbuf);
			free(pixelbuf_byte);
			free(mask_row);
			alloc_error = TRUE;
		} else {
#endif
		for (int row = 0; row < loop_max; row++) {
			guint x;
			/* row == destination row; sy == source row */
			const guint sy = downscaled ? (guint)(row * inv_scale) : (guint)row;
			const guint src_i = sy * width;

			// Precompute mask values for entire source row if mask is active
			if (apply_mask) {
				if (mask_bitpix == 8) {
#pragma omp simd
					for (x = 0; x < width; ++x) {
						mask_row[x] = mask_u8[src_i + x];
					}
				} else if (mask_bitpix == 16) {
#pragma omp simd
					for (x = 0; x < width; ++x) {
						mask_row[x] = mask_u16[src_i + x] >> 8;
					}
				} else { // 32
#pragma omp simd
					for (x = 0; x < width; ++x) {
						mask_row[x] = (BYTE)(mask_f32[src_i + x] * UCHAR_MAX);
					}
				}
			}

			if (gfit->type == DATA_FLOAT) {
				for (int c = 0 ; c < 3 ; c++) {
					WORD *line = linebuf[c];
					float *source = fsrc[c];
#pragma omp simd
					for (x = 0 ; x < width ; x++)
						line[x] = roundf_to_WORD(source[src_i + x] * USHRT_MAX_SINGLE);
				}
			} else if (norm == UCHAR_MAX) {
				for (int c = 0 ; c < 3 ; c++) {
					WORD *line = linebuf[c];
					const WORD *source = src[c];
#pragma omp simd
					for (x = 0 ; x < width ; x++)
						line[x] = source[src_i + x] << 8;
				}
			} else {
				for (int c = 0 ; c < 3 ; c++)
// No omp simd here as memcpy should already be highly optimized
					memcpy(linebuf[c], src[c] + src_i, width * sizeof(WORD));
			}
			if (gfit->type == DATA_USHORT && norm == UCHAR_MAX) {
				for (int c = 0 ; c < 3 ; c++) {
					WORD *line = linebuf[c];
#pragma omp simd
					for (x = 0 ; x < width ; x++)
						line[x] = line[x] >> 8;
				}
			}
			for (int c = 0 ; c < 3 ; c++) {
				const int cc = gfit->color_managed ? c : 0;
				for (x = 0; x < width; ++x) {
					WORD val = linebuf[c][x];
					if (gui.cut_over && val > remap_hi) {	// cut
						linebuf_byte[c][x] = 0;
					} else {
						linebuf_byte[c][x] = index[cc][val - remap_lo < 0 ? 0 : val - remap_lo];
					}
					if (inverted)
						linebuf_byte[c][x] = UCHAR_MAX - linebuf_byte[c][x];
				}
			}
			if (do_transform) {
				cmsDoTransformLineStride(com.gui_icc.proofing_transform, pixelbuf_byte, pixelbuf_byte, width, 1, width * 3, width * 3, width, width);
			}

			const guint dst_row_start = downscaled
				? (guint)(loop_max - 1 - row) * (guint)dstride
				: (height - 1 - sy) * (guint)width * 4;
			const int out_w = sw;  /* number of destination pixels to write */

			if (color == RAINBOW_COLOR) {
				if (apply_mask) {
					// Rainbow with mask
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
						for (int dx = 0 ; dx < out_w ; dx++) {
							const int sx = downscaled ? (int)(dx * inv_scale) : dx;
							const guint dst_index = dst_row_start + (guint)dx * 4;
							const BYTE pixel_val = line_byte[sx];
							const BYTE mask_val = mask_row[sx];

							const BYTE rainbow_r = rainbow_index[pixel_val][0];
							const BYTE rainbow_g = rainbow_index[pixel_val][1];
							const BYTE rainbow_b = rainbow_index[pixel_val][2];

							const uint16_t inv_mask = UCHAR_MAX - mask_val;
							const uint16_t r_sum = inv_mask + rainbow_r;
							const BYTE r_val = (r_sum > UCHAR_MAX) ? UCHAR_MAX : (BYTE)r_sum;
							const BYTE g_val = (rainbow_g * mask_val) >> 8;
							const BYTE b_val = (rainbow_b * mask_val) >> 8;

							*(guint32*)(dst[c] + dst_index) = r_val << 16 | g_val << 8 | b_val;
						}
					}
				} else {
					// Rainbow without mask
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
						for (int dx = 0 ; dx < out_w ; dx++) {
							const int sx = downscaled ? (int)(dx * inv_scale) : dx;
							const guint dst_index = dst_row_start + (guint)dx * 4;
							const BYTE pixel_val = line_byte[sx];
							*(guint32*)(dst[c] + dst_index) = rainbow_index[pixel_val][0] << 16 |
							                                   rainbow_index[pixel_val][1] << 8 |
							                                   rainbow_index[pixel_val][2];
						}
					}
				}
			} else {
				// NORMAL_COLOR
				if (apply_mask) {
					// Normal with mask
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
						for (int dx = 0 ; dx < out_w ; dx++) {
							const int sx = downscaled ? (int)(dx * inv_scale) : dx;
							const guint dst_index = dst_row_start + (guint)dx * 4;
							const BYTE pixel_val = line_byte[sx];
							const BYTE mask_val = mask_row[sx];

							const uint16_t inv_mask = UCHAR_MAX - mask_val;
							const uint16_t r_sum = inv_mask + pixel_val;
							const BYTE r_val = (r_sum > UCHAR_MAX) ? UCHAR_MAX : (BYTE)r_sum;
							const BYTE g_val = (pixel_val * mask_val) >> 8;
							const BYTE b_val = (pixel_val * mask_val) >> 8;

							*(guint32*)(dst[c] + dst_index) = r_val << 16 | g_val << 8 | b_val;
						}
					}
				} else {
					// Normal without mask
					for (int c = 0 ; c < 3 ; c++) {
						const BYTE *line_byte = linebuf_byte[c];
#pragma omp simd
						for (x = 0 ; x < (guint)out_w ; x++) {
							const int sx = downscaled ? (int)(x * inv_scale) : (int)x;
							const guint dst_index = dst_row_start + x * 4;
							const BYTE pixel_val = line_byte[sx];
							*(guint32*)(dst[c] + dst_index) = pixel_val << 16 | pixel_val << 8 | pixel_val;
						}
					}
				}
			}
		}
		// Close the else-branch of the allocation check, then free buffers.
		// In the OpenMP path free() is called by every thread for its own
		// allocations; free(NULL) is safe for the error path where some
		// allocations may have failed.
		}
#ifdef _OPENMP
		free(pixelbuf);
		free(pixelbuf_byte);
		free(mask_row);
}
#else
		free(pixelbuf);
		free(pixelbuf_byte);
		free(mask_row);
#endif
		if (do_transform)
			unlock_display_transform();
	}

	// Bail out cleanly if any thread's allocation failed
	if (alloc_error) {
		g_mutex_unlock(&gui.cairo_mutex);
		return;
	}

	// flush to ensure all writing to the image was done and redraw the surfaces
	for (int vport = 0 ; vport < 3 ; vport++) {
		cairo_surface_flush(view[vport]->full_surface);
		cairo_surface_mark_dirty(view[vport]->full_surface);
	}
	g_mutex_unlock(&gui.cairo_mutex);

	// invalidate_image_render_cache is self-locking; both it and
	// test_and_allocate_reference_image must be called after our unlock
	for (int vport = 0 ; vport < 3 ; vport++) {
		invalidate_image_render_cache(vport);
		test_and_allocate_reference_image(vport);
	}
}

static int make_hd_index_for_current_display(int vport) {
	float slope;
	int i;
	BYTE *index;
	float pxl;
	// Check if the bit depth matches the LUT size, if not we need to realloc
	if (gui.hd_remap_max != 1 << com.pref.hd_bitdepth) {
		gui.hd_remap_max = 1 << com.pref.hd_bitdepth;
		allocate_hd_remap_indices();
	}

	/* initialization of data required to build the remap_index
	 *
	 * The HD remap curve is only used with STF rendering mode
	 * as it is the only mode that results in a slope steep enough
	 * to cause noticeable quantization of levels */
	slope = UCHAR_MAX_SINGLE;
	/************* Building the HD remap_index **************/
	siril_log_debug("Rebuilding HD remap_index\n");
	int target_index = gui.rendering_mode == STF_DISPLAY && gui.use_hd_remap && gui.unlink_channels ? vport : 0;
	index = gui.hd_remap_index[target_index];

	for (i = 0; i <= gui.hd_remap_max; i++) {
		pxl = (float) i / gui.hd_remap_max;
		index[i] = roundf_to_BYTE((MTFp(pxl, stf[target_index])) * slope);
		// check for maximum overflow, given that df/di > 0. Should not happen with round_to_BYTE
		if (index[i] == UCHAR_MAX)
			break;
	}
	if (i != gui.hd_remap_max + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= gui.hd_remap_max; i++)
			index[i] = UCHAR_MAX;
	}
	return 0;
}

static int make_index_for_current_display(int vport) {
	g_mutex_lock(&com.mutex);
	WORD lo = gui.lo;
	WORD hi = gui.hi;
	g_mutex_unlock(&com.mutex);
	float slope, delta = hi - lo;
	int i;
	BYTE *index;
	float pxl;
	/* initialization of data required to build the remap_index */
	switch (gui.rendering_mode) {
		case LINEAR_DISPLAY:
			slope = UCHAR_MAX_SINGLE / delta;
			break;
		case LOG_DISPLAY:
			slope = fabsf(UCHAR_MAX_SINGLE / logf(delta * 0.1f));
			break;
		case SQRT_DISPLAY:
			slope = UCHAR_MAX_SINGLE / sqrtf(delta);
			break;
		case SQUARED_DISPLAY:
			slope = UCHAR_MAX_SINGLE / SQR(delta);
			break;
		case ASINH_DISPLAY:
			slope = UCHAR_MAX_SINGLE / asinhf(delta * 0.001f);
			break;
		case STF_DISPLAY:
			slope = UCHAR_MAX_SINGLE;
			break;
		default:
			return 1;
	}
	if(!(slope == last_pente && gui.rendering_mode == last_mode))
		com.gui_icc.profile_changed = TRUE;

	if ((gui.rendering_mode != HISTEQ_DISPLAY && gui.rendering_mode != STF_DISPLAY) &&
			slope == last_pente && gui.rendering_mode == last_mode && !com.gui_icc.profile_changed) {
		siril_log_debug("Re-using previous gui.remap_index\n");
		return 0;
	}

	/************* Building the remap_index **************/
	siril_log_debug("Rebuilding gui.remap_index %d\n", vport);
	// target_index only used for STF mode
	int target_index = gui.rendering_mode == STF_DISPLAY && gui.unlink_channels ? vport : 0;
	index = gui.remap_index[vport];
	for (i = 0; i <= USHRT_MAX; i++) {
		switch (gui.rendering_mode) {
			case LOG_DISPLAY:
				// ln(5.56*10^110) = 255
				if (i < 10)
					index[i] = 0; /* avoid null and negative values */
				else
					index[i] = roundf_to_BYTE(logf((float) i / 10.f) * slope); //10.f is arbitrary: good matching with ds9
				break;
			case SQRT_DISPLAY:
				// sqrt(2^16) = 2^8
				index[i] = roundf_to_BYTE(sqrtf((float) i) * slope);
				break;
			case SQUARED_DISPLAY:
				// pow(2^4,2) = 2^8
				index[i] = roundf_to_BYTE(SQR((float)i) * slope);
				break;
			case ASINH_DISPLAY:
				// asinh(2.78*10^110) = 255
				index[i] = roundf_to_BYTE(asinhf((float) i / 1000.f) * slope); //1000.f is arbitrary: good matching with ds9, could be asinhf(a*Q*i)/Q
				break;
			case LINEAR_DISPLAY:
				index[i] = roundf_to_BYTE((float) i * slope);
				break;
			case STF_DISPLAY:
				pxl = (gfit->orig_bitpix == BYTE_IMG ?
						(float) i / UCHAR_MAX_SINGLE :
						(float) i / USHRT_MAX_SINGLE);
				index[i] = roundf_to_BYTE((MTFp(pxl, stf[target_index])) * slope);
				break;
			default:
				return 1;
		}
		// check for maximum overflow, given that df/di > 0. Should not happen with round_to_BYTE
		if (index[i] == UCHAR_MAX)
			break;
	}
	if (i != USHRT_MAX + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= USHRT_MAX; i++) {
			index[i] = UCHAR_MAX;
		}
	}
	if (gfit->color_managed && com.gui_icc.same_primaries && com.gui_icc.proofing_transform && gui.rendering_mode != STF_DISPLAY)
		display_index_transform(index, vport);

	last_pente = slope;
	last_mode = gui.rendering_mode;
	return 0;
}

static int make_index_for_rainbow(BYTE index[][3]) {
	int i;
	double h, s, v, r, g, b;

	for (i = 0; i < UCHAR_MAX + 1; i++) {
		r = g = b = (double) i / UCHAR_MAX_DOUBLE;
		rgb_to_hsv(r, g, b, &h, &s, &v);
		double off = 300.0 / 360.0; /* Arbitrary: we want h from 300 to 0 deg */
		h = (off - (double) i * (off / UCHAR_MAX_DOUBLE));
		s = 1.;
		v = 1.; /* Saturation and Value are set to 100% */
		hsv_to_rgb(h, s, v, &r, &g, &b);
		index[i][0] = round_to_BYTE(r * UCHAR_MAX_DOUBLE);
		index[i][1] = round_to_BYTE(g * UCHAR_MAX_DOUBLE);
		index[i][2] = round_to_BYTE(b * UCHAR_MAX_DOUBLE);
	}
	return 0;
}

/*****************************************************************************
 * ^ ^ ^ above:     R E M A P P I N G     I M A G E     D A T A        ^ ^ ^ *
 *                                                                           *
 * v v v below:          R E D R A W I N G     W I D G E T S           v v v *
 *****************************************************************************/

typedef struct label_point_struct {
	double x, y, ra, dec, angle;
	gboolean isRA;
	int border;
} label_point;

static void request_gtk_redraw_of_cvport() {
	//siril_log_debug("image redraw requested (vport %d)\n", gui.cvport);
	GtkWidget *widget = gui.view[gui.cvport].drawarea;
	gtk_widget_queue_draw(widget);
}

/* forward declaration - redraw_drawingarea is defined later in this file */
gboolean redraw_drawingarea(GtkWidget *widget, cairo_t *cr, gpointer data);

/* Block / unblock the draw signal handler on all viewports.  Must be called
 * from the GTK main thread.  Used by generic_image_worker (via
 * execute_idle_and_wait_for_it) to prevent redraws with stale Cairo buffers
 * while the worker is running. */
void block_drawarea_handlers(void) {
	for (int i = 0; i < MAXVPORT; i++)
		if (gui.view[i].drawarea)
			g_signal_handlers_block_by_func(gui.view[i].drawarea,
					redraw_drawingarea, NULL);
}

void unblock_drawarea_handlers(void) {
	for (int i = 0; i < MAXVPORT; i++)
		if (gui.view[i].drawarea)
			g_signal_handlers_unblock_by_func(gui.view[i].drawarea,
					redraw_drawingarea, NULL);
}

static void draw_empty_image(const draw_data_t* dd) {
	static GdkPixbuf *siril_pix = NULL;
	cairo_t *cr = dd->cr;
	guint width = dd->window_width;
	guint height = dd->window_height;
	guint pix_size = height / 3;
	guint offset;
#ifdef SIRIL_UNSTABLE
	offset = 32;
#else
	offset = 2;
#endif

	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	/* Create pixbuf from siril.svg file */
	if (siril_pix == NULL) {
		siril_pix = gdk_pixbuf_new_from_resource_at_scale("/org/siril/ui/pixmaps/siril.svg", 256, 256, FALSE, NULL);
	}

	GdkPixbuf *pixbuf = gdk_pixbuf_scale_simple(siril_pix, pix_size, pix_size, GDK_INTERP_BILINEAR);

	gdk_cairo_set_source_pixbuf(cr, pixbuf, (width - pix_size) / 2, (height - pix_size) / offset);
	cairo_paint(cr);
	cairo_fill(cr);

	g_object_unref(pixbuf);


	image_display_init_statics();
	GtkWidget *widget = imgdisp_drawing_rgb;
	GtkStyleContext *context = gtk_widget_get_style_context(widget);
	GtkStateFlags state = gtk_widget_get_state_flags(widget);
	PangoLayout *layout;
	gchar *msg;
	GtkAllocation allocation;
	gdouble scale;
	GdkRGBA color;
	gint w, h;

	layout = gtk_widget_create_pango_layout(widget, NULL);

#ifdef SIRIL_UNSTABLE

	msg = g_strdup_printf(_("<big>Unstable Development Version</big>\n\n"
			    "<small>%c%s</small>\n"
				"<small>commit <tt>%s</tt></small>\n"
				"<small>Please test bugs against "
				"latest git master branch\n"
				"before reporting them.</small>"),
			toupper(PACKAGE_STRING[0]),
			(char *)PACKAGE_STRING + 1,
			SIRIL_GIT_VERSION_ABBREV);
	pango_layout_set_markup(layout, msg, -1);
	g_free(msg);
	pango_layout_set_alignment(layout, PANGO_ALIGN_CENTER);

	pango_layout_get_pixel_size(layout, &w, &h);
	gtk_widget_get_allocation(widget, &allocation);

	scale = MIN(((gdouble ) allocation.width / 2.0) / (gdouble ) w,
			((gdouble ) allocation.height / 2.0) / (gdouble ) h / 2);

	gtk_style_context_get_color(context, state, &color);
	gdk_cairo_set_source_rgba(cr, &color);

	cairo_move_to(cr, (allocation.width - (w * scale)) / 2,
			(allocation.height - (h * scale)) / 2);

#else
	msg = g_strdup_printf("%c%s", toupper(PACKAGE_STRING[0]), (char *)PACKAGE_STRING + 1);

	pango_layout_set_markup(layout, msg, -1);
	g_free(msg);
	pango_layout_set_alignment(layout, PANGO_ALIGN_CENTER);

	pango_layout_get_pixel_size(layout, &w, &h);
	gtk_widget_get_allocation(widget, &allocation);

	scale = MIN(((gdouble ) allocation.width / 4.0) / (gdouble ) w,
			((gdouble ) allocation.height / 4.0) / (gdouble ) h / 4);

	gtk_style_context_get_color(context, state, &color);
	gdk_cairo_set_source_rgba(cr, &color);

	cairo_move_to(cr, (allocation.width - (w * scale)) / 2,
			3 * (allocation.height - (h * scale)) / 4);

#endif /* SIRIL_UNSTABLE */
	cairo_scale(cr, scale, scale);

	pango_cairo_show_layout(cr, layout);

	g_object_unref(layout);
}

/* Fill disp_surface at full image resolution when zoom >= 1.0 on images that
 * exceed Cairo's 32767-pixel limit (gui.surface_scale < 1.0).  Reads directly
 * from gfit pixel data and applies the same LUT as the regular remap path,
 * avoiding the blurry upscale of the downsampled full_surface.
 * Must be called with gui.cairo_mutex held.
 * Returns TRUE on success, FALSE on allocation failure. */
static gboolean fill_hires_disp(const draw_data_t *dd, struct image_view *view,
                                 WORD remap_hi, WORD remap_lo) {
	unsigned char *disp_data = cairo_image_surface_get_data(view->disp_surface);
	const int disp_stride = cairo_image_surface_get_stride(view->disp_surface);
	if (!disp_data) return FALSE;

	const int win_w = dd->window_width;
	const int win_h = dd->window_height;
	const int img_w = (int)gfit->rx;
	const int img_h = (int)gfit->ry;
	const double zoom  = dd->zoom;
	const double off_x = gui.display_offset.x;
	const double off_y = gui.display_offset.y;

	const gboolean inverted = dd->neg_view;
	const gboolean cut_over = gui.cut_over;
	const gboolean hd_mode  = (gui.rendering_mode == STF_DISPLAY &&
	                            gui.use_hd_remap && gfit->type == DATA_FLOAT);
	const gboolean is_rgb   = (dd->vport == RGB_VPORT);

	/* x range of screen columns that overlap the image */
	const int dx_first = MAX(0, (int)floor(off_x));
	const int dx_last  = MIN(win_w - 1, (int)ceil(off_x + (double)img_w * zoom) - 1);
	const int visible_w = (dx_last >= dx_first) ? (dx_last - dx_first + 1) : 0;

	/* Clear entire surface to opaque black */
	memset(disp_data, 0, (size_t)win_h * disp_stride);
	if (visible_w <= 0) return TRUE;

	/* RGB path: per-row channel buffers + optional ICC transform */
	BYTE *lb_rgb = NULL;
	gboolean do_transform = FALSE;
	cmsHTRANSFORM transform = NULL;
	if (is_rgb) {
		lb_rgb = malloc((size_t)visible_w * 3);
		if (!lb_rgb) return FALSE;
		lock_display_transform();
		if (gfit->color_managed && com.gui_icc.proofing_transform && !identical &&
				!com.gui_icc.same_primaries) {
			do_transform = TRUE;
			transform = com.gui_icc.proofing_transform;
		}
		unlock_display_transform();
	}

	for (int dy = 0; dy < win_h; dy++) {
		/* Screen row → image buf row → gfit row.
		 * buf stores rows bottom-up (gfit row 0 = buf last row), so we flip. */
		const double fy  = (dy - off_y) / zoom;
		const int buf_row = (int)fy;
		if (buf_row < 0 || buf_row >= img_h) continue;
		const int iy = img_h - 1 - buf_row;

		uint32_t *row_out = (uint32_t *)(disp_data + dy * disp_stride);

		if (!is_rgb) {
			/* Single-channel vport: mirrors REMAP_WRITE_PIXEL logic */
			const int vp = (dd->vport < 3) ? dd->vport : 0;
			const int ti = (gui.rendering_mode == STF_DISPLAY && gui.unlink_channels) ? vp : 0;
			const BYTE *lut = hd_mode ? gui.hd_remap_index[ti] : gui.remap_index[ti];

			for (int dx = dx_first; dx <= dx_last; dx++) {
				const int ix = (int)((dx - off_x) / zoom);
				if ((unsigned)ix >= (unsigned)img_w) continue;
				const size_t si = (size_t)iy * img_w + ix;
				BYTE val;
				if (gfit->type == DATA_FLOAT) {
					val = hd_mode
						? lut[float_to_max_range(gfit->fpdata[vp][si], gui.hd_remap_max)]
						: lut[roundf_to_WORD(gfit->fpdata[vp][si] * USHRT_MAX_SINGLE)];
				} else {
					const WORD sv = gfit->pdata[vp][si];
					val = hd_mode ? lut[(guint)sv * gui.hd_remap_max / USHRT_MAX] : lut[sv];
				}
				if (inverted) val = UCHAR_MAX - val;
				row_out[dx] = 0xFF000000u | ((uint32_t)val << 16) | ((uint32_t)val << 8) | val;
			}
		} else {
			/* RGB composite vport: mirrors remap_all_vports logic */
			BYTE *lb_r = lb_rgb;
			BYTE *lb_g = lb_rgb + visible_w;
			BYTE *lb_b = lb_rgb + 2 * visible_w;

			for (int dx = dx_first; dx <= dx_last; dx++) {
				const int ix = (int)((dx - off_x) / zoom);
				const int li = dx - dx_first;
				if ((unsigned)ix >= (unsigned)img_w) {
					lb_r[li] = lb_g[li] = lb_b[li] = 0;
					continue;
				}
				const size_t si = (size_t)iy * img_w + ix;
				BYTE *channels[3] = { lb_r, lb_g, lb_b };
				for (int c = 0; c < 3; c++) {
					const int cc = gfit->color_managed ? c : 0;
					WORD w;
					if (gfit->type == DATA_FLOAT)
						w = roundf_to_WORD(gfit->fpdata[c][si] * USHRT_MAX_SINGLE);
					else
						w = gfit->pdata[c][si];
					const int sh = (int)w - remap_lo;
					BYTE val = (cut_over && w > remap_hi) ? 0 : gui.remap_index[cc][sh < 0 ? 0 : sh];
					if (inverted) val = UCHAR_MAX - val;
					channels[c][li] = val;
				}
			}

			if (do_transform && transform)
				cmsDoTransformLineStride(transform, lb_rgb, lb_rgb, visible_w, 1,
					visible_w * 3, visible_w * 3, visible_w, visible_w);

			for (int dx = dx_first; dx <= dx_last; dx++) {
				const int li = dx - dx_first;
				row_out[dx] = 0xFF000000u
					| ((uint32_t)lb_r[li] << 16)
					| ((uint32_t)lb_g[li] <<  8)
					|  (uint32_t)lb_b[li];
			}
		}
	}

	free(lb_rgb);
	cairo_surface_flush(view->disp_surface);
	cairo_surface_mark_dirty(view->disp_surface);
	return TRUE;
}

static void draw_vport(const draw_data_t* dd) {
	/* Snapshot lo/hi before cairo_mutex to preserve lock ordering
	 * (remap takes com.mutex then cairo_mutex; never the reverse). */
	g_mutex_lock(&com.mutex);
	WORD remap_hi = gui.hi;
	WORD remap_lo = gui.lo;
	g_mutex_unlock(&com.mutex);

	g_mutex_lock(&gui.cairo_mutex);
	struct image_view *view = &gui.view[dd->vport];
	if (!view->disp_surface) {
		cairo_surface_t *target = cairo_get_target(dd->cr);
		view->disp_surface = cairo_surface_create_similar_image(target, CAIRO_FORMAT_ARGB32,
					dd->window_width, dd->window_height);
		if (cairo_surface_status(view->disp_surface) != CAIRO_STATUS_SUCCESS) {
			siril_log_debug("Error creating the cairo image disp_surface for vport %d\n", dd->vport);
			cairo_surface_destroy(view->disp_surface);
			view->disp_surface = NULL;
			g_mutex_unlock(&gui.cairo_mutex);
			return;
		}
		view->view_width = dd->window_width;
		view->view_height = dd->window_height;

		/* For large images at zoom >= 1.0 (surface_scale < 1.0), paint at full
		 * resolution directly from gfit data rather than upscaling the small
		 * surface.
		 * Skip when generic_image_worker is running: gfit->pdata may be
		 * mid-modification; fall back to view->buf via full_surface instead. */
		const gboolean do_hires = (gui.surface_scale < 1.0 && dd->zoom >= 1.0
		                           && !g_atomic_int_get(&gui.suppress_drawarea_redraw));
		if (do_hires && fill_hires_disp(dd, view, remap_hi, remap_lo)) {
			/* hires surface filled - skip standard surface paint */
		} else {
		cairo_t *cached_cr = cairo_create(view->disp_surface);
		cairo_matrix_t y_reflection_matrix, flipped_matrix;
		cairo_matrix_init_identity(&y_reflection_matrix);
		if (livestacking_is_started() && !g_strcmp0(gfit->keywords.row_order, "TOP-DOWN")) {
			y_reflection_matrix.yy = -1.0;
			/* y0 must be the surface height, not the image height, so the
			 * flip stays in surface coordinate space. */
			y_reflection_matrix.y0 = (double)view->full_surface_height;
		}
		/* When gui.surface_scale < 1, full_surface is smaller than gfit.
		 * Compensate in the matrix so the surface renders at the correct
		 * zoom level (image-space coordinates remain unchanged for overlays). */
		cairo_matrix_t surf_disp_matrix;
		double ss = (view->full_surface_height > 0 && gfit->ry > 0)
			? (double)view->full_surface_height / (double)gfit->ry : 1.0;
		double surf_zoom = gui.display_matrix.xx / ss;
		cairo_matrix_init(&surf_disp_matrix,
				surf_zoom, 0, 0, surf_zoom,
				gui.display_matrix.x0, gui.display_matrix.y0);
		cairo_matrix_multiply(&flipped_matrix, &y_reflection_matrix, &surf_disp_matrix);
		cairo_transform(cached_cr, &flipped_matrix);
		cairo_set_source_surface(cached_cr, view->full_surface, 0, 0);
		cairo_pattern_set_filter(cairo_get_source(cached_cr), dd->filter);
		cairo_paint(cached_cr);
		cairo_destroy(cached_cr);
		}

//		siril_log_debug("@@@\t\t\tcache surface created (%d x %d)\t\t\t@@@\n",
//				view->view_width, view->view_height);
	}
	cairo_set_source_surface(dd->cr, view->disp_surface, 0, 0);
	cairo_paint(dd->cr);

	// prepare the display matrix for remaining drawing (selection, stars, ...)
	cairo_transform(dd->cr, &gui.display_matrix);
	g_mutex_unlock(&gui.cairo_mutex);

}

static void draw_main_image(const draw_data_t* dd) {
	g_mutex_lock(&gui.cairo_mutex);
	gboolean has_buf = (gui.view[dd->vport].buf != NULL);
	g_mutex_unlock(&gui.cairo_mutex);

	if (has_buf) {
		draw_vport(dd);
	} else {
		draw_empty_image(dd);
	}
}

gboolean get_context_rotation_matrix(double rotation, cairo_matrix_t *transform, gboolean invert) {
	if (rotation == 0.) return FALSE;
	double dx = (double)com.selection.x + (double)com.selection.w * 0.5;
	double dy = (double)com.selection.y + (double)com.selection.h * 0.5;
	cairo_matrix_init_translate(transform, dx, dy);
	cairo_matrix_rotate(transform, rotation * DEGTORAD);
	cairo_matrix_translate(transform, -dx, -dy);
	if (invert) return (cairo_matrix_invert(transform) == CAIRO_STATUS_SUCCESS);
	return TRUE;
}

static void rotate_context(cairo_t *cr, double rotation) {
	cairo_matrix_t transform;
	if (!get_context_rotation_matrix(rotation, &transform, FALSE)) return;
	cairo_transform(cr, &transform);
}

static void draw_roi(const draw_data_t *dd) {
	double r, g, b;
	if (gui.roi.operation_supports_roi) {
		r = 0.3; g = 1.0; b = 0.3;
	} else {
		r = 1.0; g = 0.0; b = 0.0;
	}
	if (gui.roi.selection.w > 0 && gui.roi.selection.h > 0 && gui.roi.active) {
		cairo_t *cr = dd->cr;
		static double dash_format[] = { 4.0, 2.0 };
		cairo_set_line_width(cr, 1.5 / dd->zoom);
		cairo_set_dash(cr, dash_format, 2, 0);
		cairo_set_source_rgb(cr, r, g, b);
		cairo_save(cr); // save the original transform
		cairo_rectangle(cr, (double) gui.roi.selection.x, (double) gui.roi.selection.y,
						(double) gui.roi.selection.w, (double) gui.roi.selection.h);
		cairo_stroke(cr);
		cairo_restore(cr);
	}
}

static void draw_selection(const draw_data_t* dd) {
	if (com.selection.w > 0 && com.selection.h > 0) {
		if ((com.selection.x + com.selection.w > gfit->rx) ||
		(com.selection.y + com.selection.h > gfit->ry)) {
			rectangle area = {0, 0, gfit->rx, gfit->ry};
			memcpy(&com.selection, &area, sizeof(rectangle));
		}
		if (!rotation_dlg) image_display_init_statics();
		cairo_t *cr = dd->cr;
		static double dash_format[] = { 4.0, 2.0 };
		cairo_set_line_width(cr, 1.5 / dd->zoom);
		cairo_set_dash(cr, dash_format, 2, 0);
		cairo_set_source_rgb(cr, 0.8, 1.0, 0.8);
		cairo_save(cr); // save the original transform
		if (gtk_widget_is_visible(rotation_dlg)) {
			double dashes2[]={5.0, 5.0};
			cairo_set_dash(cr, dashes2, 2, 0);
			cairo_set_line_width(cr, 0.5 / dd->zoom);
			cairo_rectangle(cr, (double) com.selection.x, (double) com.selection.y,
						(double) com.selection.w, (double) com.selection.h);
			cairo_stroke(cr);
			cairo_set_line_width(cr, 3. / dd->zoom);
			cairo_set_source_rgb(cr, 0.8, 0.0, 0.0);
			rotate_context(cr, -gui.rotation); // cairo is positive CW while opencv is positive CCW

			// draw a circle at top left corner to visualize rots larger than 90
			double size = 10. / dd->zoom;
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_arc(cr, com.selection.x, com.selection.y, size * 0.5, 0., 2. * M_PI);
			cairo_stroke_preserve(cr);
			cairo_fill(cr);
			cairo_set_dash(cr, dash_format, 2, 0);
		}
		cairo_rectangle(cr, (double) com.selection.x, (double) com.selection.y,
						(double) com.selection.w, (double) com.selection.h);
		cairo_stroke(cr);

		// display a grid when the selection is being made / modified, when it is big enough
		if (com.pref.gui.selection_guides > 1 && gui.drawing && com.selection.w > 40 / dd->zoom && com.selection.h > 40 / dd->zoom) {
			cairo_set_line_width(cr, 0.4 / dd->zoom);
			cairo_set_dash(cr, NULL, 0, 0);
			for (int i = 1; i < com.pref.gui.selection_guides; i++) {
				int x = com.selection.x + com.selection.w * i / com.pref.gui.selection_guides;
				int y = com.selection.y + com.selection.h * i / com.pref.gui.selection_guides;
				cairo_move_to(cr, x, com.selection.y);
				cairo_line_to(cr, x, com.selection.y + com.selection.h);
				cairo_move_to(cr, com.selection.x, y);
				cairo_line_to(cr, com.selection.x + com.selection.w, y);
			}
			cairo_stroke(cr);
		}

		// display a mini cross when the selection is being dragged
		if ((gui.freezeX && gui.freezeY) || gtk_widget_is_visible(rotation_dlg)) {
			cairo_set_line_width(cr, 1.0 / dd->zoom);
			point selection_center = { com.selection.x + (double)com.selection.w / 2.,
				com.selection.y + (double)com.selection.h / 2. };
			cairo_move_to(cr, selection_center.x, selection_center.y - 5 / dd->zoom);
			cairo_line_to(cr, selection_center.x, selection_center.y + 5 / dd->zoom);
			cairo_move_to(cr, selection_center.x - 5 / dd->zoom, selection_center.y);
			cairo_line_to(cr, selection_center.x + 5 / dd->zoom, selection_center.y);
			cairo_stroke(cr);
		}
		cairo_restore(cr); // restore the original transform
	}
}

static void draw_cut_line(const draw_data_t* dd) {
	if (!cut_dialog)
		image_display_init_statics();
	if (!(gtk_widget_get_visible(cut_dialog) || gtk_widget_get_visible(cut_cdialog) || gtk_widget_get_visible(cut_sdialog)))
		return;
	if (gui.cut.cut_end.x == -1 || gui.cut.cut_end.y == -1 || gui.cut.seq)
		return;
	gboolean tri = gtk_toggle_button_get_active(tri_cut_toggle);
	double offstartx, offstarty, offendx, offendy, step;

	cairo_t *cr = dd->cr;
	static double dash_format[] = { 4.0, 2.0 };
	static double solid_format[] = { 1.0, 0.0 };
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	cairo_set_dash(cr, dash_format, 2, 0);

	if (tri) {
		point delta;
		delta.x = gui.cut.cut_end.x - gui.cut.cut_start.x;
		delta.y = gui.cut.cut_end.y - gui.cut.cut_start.y;
		double length = sqrt(delta.x * delta.x + delta.y * delta.y);
		if (length < 1.) return;
		int nbr_points = (int) length;
		double point_spacing_x = delta.x / nbr_points;
		double point_spacing_y = delta.y / nbr_points;
		step = gtk_spin_button_get_value(tri_cut_spin_step);
		double line_r[3] = { 0.58, 0.0, 0.34 }; // These colours match the 3 lines plotted by siril plot
		double line_g[3] = { 0.0, 0.62, 0.70 };
		double line_b[3] = { 0.83, 0.45, 0.91 };
		double arrow_length = 10 / dd->zoom;
		double arrow_angle = 0.5;
		double angle = atan2(gui.cut.cut_end.y - gui.cut.cut_start.y, gui.cut.cut_end.x - gui.cut.cut_start.x);
		for (int offset = -1 ; offset < 2 ; offset++) {
			cairo_set_dash(cr, dash_format, 2, 0);
			offstartx = gui.cut.cut_start.x + (offset * point_spacing_y * step);
			offstarty = gui.cut.cut_start.y - (offset * point_spacing_x * step);
			offendx = gui.cut.cut_end.x + (offset * point_spacing_y * step);
			offendy = gui.cut.cut_end.y - (offset * point_spacing_x * step);
			cairo_set_source_rgb(cr, line_r[offset+1], line_g[offset+1], line_b[offset+1]);
			cairo_save(cr);
			cairo_move_to(cr, offstartx + 0.5, offstarty + 0.5);
			cairo_line_to(cr, offendx + 0.5, offendy + 0.5);
			cairo_stroke(cr);
			// Draw arrowheads at the end
			cairo_set_dash(cr, solid_format, 0, 0); // Draw the arrow heads solid
			point pt1 = { offendx + 0.5 - arrow_length * cos(angle - arrow_angle), offendy + 0.5 - arrow_length * sin(angle - arrow_angle) };
			point pt2 = { offendx + 0.5 - arrow_length * cos(angle + arrow_angle), offendy + 0.5 - arrow_length * sin(angle + arrow_angle) };
			cairo_line_to(cr, offendx + 0.5, offendy + 0.5);
			cairo_line_to(cr, pt1.x, pt1.y);
			cairo_move_to(cr, offendx + 0.5, offendy + 0.5);
			cairo_line_to(cr, pt2.x, pt2.y);
			cairo_stroke(cr);
			cairo_restore(cr);
		}
	} else {
		cairo_set_source_rgb(cr, 0.0, 0.62, 0.70); // This matches the single line plotted by siril plot
		cairo_save(cr);
		cairo_move_to(cr, gui.cut.cut_start.x + 0.5, gui.cut.cut_start.y + 0.5);
		cairo_line_to(cr, gui.cut.cut_end.x + 0.5, gui.cut.cut_end.y + 0.5);
		// Draw an arrowhead at the end
		double arrow_length = 10 / dd->zoom;
		double arrow_angle = 0.5;
		double angle = atan2(gui.cut.cut_end.y - gui.cut.cut_start.y, gui.cut.cut_end.x - gui.cut.cut_start.x);
		cairo_stroke(cr);
		cairo_set_dash(cr, solid_format, 0, 0); // Draw the arrow heads solid
		point pt1 = { gui.cut.cut_end.x + 0.5 - arrow_length * cos(angle - arrow_angle), gui.cut.cut_end.y + 0.5 - arrow_length * sin(angle - arrow_angle) };
		point pt2 = { gui.cut.cut_end.x + 0.5 - arrow_length * cos(angle + arrow_angle), gui.cut.cut_end.y + 0.5 - arrow_length * sin(angle + arrow_angle) };
		cairo_line_to(cr, gui.cut.cut_end.x + 0.5, gui.cut.cut_end.y + 0.5);
		cairo_line_to(cr, pt1.x, pt1.y);
		cairo_move_to(cr, gui.cut.cut_end.x + 0.5, gui.cut.cut_end.y + 0.5);
		cairo_line_to(cr, pt2.x, pt2.y);
		cairo_stroke(cr);
		cairo_restore(cr);
	}
}

static void draw_measurement_line(const draw_data_t* dd) {
	if (gui.measure_start.x == -1)
		return;
	cairo_t *cr = dd->cr;
	static double dash_format[] = { 4.0, 2.0 };
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	cairo_set_dash(cr, dash_format, 2, 0);
	cairo_set_source_rgb(cr, 0.8, 1.0, 0.8);
	cairo_save(cr);
	cairo_move_to(cr, gui.measure_start.x + 0.5, gui.measure_start.y + 0.5);
	cairo_line_to(cr, gui.measure_end.x + 0.5, gui.measure_end.y + 0.5);
	cairo_stroke(cr);
	cairo_restore(cr);
}

static void draw_user_polygons(const draw_data_t *dd) {
	if (gui.user_polygons == NULL)
		return;
	cairo_t *cr = dd->cr;
	static double dash_format[] = { 4.0, 2.0 };
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	cairo_set_dash(cr, dash_format, 2, 0);
	GSList *l;
	for (l = gui.user_polygons; l != NULL; l = l->next) {
		UserPolygon *polygon = (UserPolygon *)l->data;
		if (polygon->n_points < 2)
			continue;

		// Calculate center of polygon for legend placement
		double center_x = 0, center_y = 0;
		for (int i = 0; i < polygon->n_points; i++) {
			center_x += polygon->points[i].x;
			center_y += polygon->points[i].y;
		}
		center_x /= polygon->n_points;
		center_y /= polygon->n_points;

		if (polygon->fill) {
			cairo_save(cr);
			// Set color for filling
			cairo_set_source_rgba(cr,
								  polygon->color[0],
						 polygon->color[1],
						 polygon->color[2],
						 polygon->color[3]);
			cairo_move_to(cr, polygon->points[0].x + 0.5, polygon->points[0].y + 0.5);
			for (int i = 1; i < polygon->n_points; i++) {
				cairo_line_to(cr, polygon->points[i].x + 0.5, polygon->points[i].y + 0.5);
			}
			cairo_close_path(cr);
			cairo_fill_preserve(cr);
			/* slightly darker outline when filled */
			cairo_set_source_rgba(cr,
								  polygon->color[0] * 0.8,
						 polygon->color[1] * 0.8,
						 polygon->color[2] * 0.8,
						 polygon->color[3]);
			cairo_stroke(cr);
			cairo_restore(cr);
		} else {
			cairo_save(cr);
			cairo_set_source_rgba(cr,
								  polygon->color[0],
						 polygon->color[1],
						 polygon->color[2],
						 polygon->color[3]);
			cairo_move_to(cr, polygon->points[0].x + 0.5, polygon->points[0].y + 0.5);
			for (int i = 1; i < polygon->n_points; i++) {
				cairo_line_to(cr, polygon->points[i].x + 0.5, polygon->points[i].y + 0.5);
			}
			cairo_close_path(cr);
			cairo_stroke(cr);
			cairo_restore(cr);
		}

		// Draw legend text if it exists
		if (polygon->legend != NULL && polygon->legend[0] != '\0') {
			cairo_save(cr);

			// Set up text properties
			cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
			cairo_set_font_size(cr, 12.0 / dd->zoom);

			// Get text dimensions to center properly
			cairo_text_extents_t text_extents;
			cairo_text_extents(cr, polygon->legend, &text_extents);

			// Create a background for the text for better visibility
			double padding = 4.0 / dd->zoom;
			cairo_set_source_rgba(cr, 0, 0, 0, 0.5);  // Darken the background
			cairo_rectangle(cr,
							center_x - text_extents.width/2 - padding,
				   center_y - text_extents.height/2 - padding,
				   text_extents.width + 2*padding,
				   text_extents.height + 2*padding);
			cairo_fill(cr);

			// Draw the text
			cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);  // White text
			cairo_move_to(cr,
						  center_x - text_extents.width/2,
				 center_y + text_extents.height/2);
			cairo_show_text(cr, polygon->legend);
			cairo_restore(cr);
		}
	}
}

static void draw_stars(const draw_data_t* dd) {
	cairo_t *cr = dd->cr;
	int i = 0;

	if (!com.script && (single_image_is_loaded() || sequence_is_loaded()) &&
			g_rw_lock_reader_trylock(&com.stars_lock)) {
		/* com.stars is a NULL-terminated array */
		if (com.stars) {
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
		cairo_set_line_width(cr, 1.5 / dd->zoom);

		while (com.stars[i]) {
			double size = com.stars[i]->fwhmx * 2.0;
			if (size <= 0.0) size = com.pref.phot_set.aperture;
			if (i == gui.selected_star) {
				// We draw horizontal and vertical lines to show the star
				cairo_set_line_width(cr, 2.0 / dd->zoom);
				cairo_set_source_rgba(cr, 0.0, 0.4, 1.0, 0.6);

				cairo_move_to(cr, com.stars[i]->xpos, 0);
				cairo_line_to(cr, com.stars[i]->xpos, dd->image_height);
				cairo_stroke(cr);
				cairo_move_to(cr, 0, com.stars[i]->ypos);
				cairo_line_to(cr, dd->image_width, com.stars[i]->ypos);
				cairo_stroke(cr);

				cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
				cairo_set_line_width(cr, 1.5 / dd->zoom);
			}
			if (com.stars[i]->has_saturated) {
				cairo_set_source_rgba(cr, 0.75, 0.22, 1.0, 0.9);
				cairo_set_line_width(cr, 3.0 / dd->zoom);
			}
			cairo_save(cr); // save the original transform
			cairo_translate(cr, com.stars[i]->xpos, com.stars[i]->ypos);
			cairo_rotate(cr, M_PI * 0.5 + com.stars[i]->angle * M_PI / 180.);
			double r = com.stars[i]->fwhmx > 0.0 ? com.stars[i]->fwhmy / com.stars[i]->fwhmx : 1.0;
			cairo_scale(cr, r, 1);
			cairo_arc(cr, 0., 0., size, 0., 2 * M_PI);
			cairo_restore(cr); // restore the original transform
			cairo_stroke(cr);
			/* to keep  for debugging boxes adjustements */
			// if (com.stars[i]->R > 0)
			// 	cairo_rectangle(cr, com.stars[i]->xpos - (double)com.stars[i]->R, com.stars[i]->ypos - (double)com.stars[i]->R, (double)com.stars[i]->R * 2 + 1, (double)com.stars[i]->R * 2 + 1);
			// cairo_stroke(cr);
			if (com.stars[i]->has_saturated) {
				cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
				cairo_set_line_width(cr, 1.5 / dd->zoom);
			}

			i++;
		}
		} // if (com.stars)
		g_rw_lock_reader_unlock(&com.stars_lock);
	}

	/* quick photometry */
	if (!com.script && gui.qphot && mouse_status == MOUSE_ACTION_PHOTOMETRY) {
		double size = (!com.pref.phot_set.force_radius) ? 0.5 * gui.qphot->fwhmx * com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
		if (size <= 0.0) size = com.pref.phot_set.aperture;

		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
		cairo_set_line_width(cr, 1.5 / dd->zoom);

		/* fwhm * 2: first circle */
		cairo_arc(cr, gui.qphot->xpos, gui.qphot->ypos, size, 0., 2. * M_PI);
		cairo_stroke(cr);

		/* sky annulus */
		if (dd->neg_view) {
			cairo_set_source_rgba(cr, 0.5, 0.0, 0.7, 0.9);
		} else {
			cairo_set_source_rgba(cr, 0.5, 1.0, 0.3, 0.9);
		}

		cairo_arc(cr, gui.qphot->xpos, gui.qphot->ypos, com.pref.phot_set.inner, 0., 2. * M_PI);
		cairo_stroke(cr);
		cairo_arc(cr, gui.qphot->xpos, gui.qphot->ypos, com.pref.phot_set.outer, 0., 2. * M_PI);
		cairo_stroke(cr);
		cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
		cairo_set_font_size(cr, 40.0 / dd->zoom);
		cairo_move_to(cr, gui.qphot->xpos + com.pref.phot_set.outer + 5, gui.qphot->ypos);
		cairo_show_text(cr, "V");  // was missing - stroke on empty path was a no-op
		cairo_stroke(cr);
	}

	/* draw seqpsf stars */
	if (sequence_is_loaded() && com.seq.current >= 0) {
		for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++) {
			psf_star *the_psf = com.seq.photometry[i][com.seq.current];
			if (the_psf) {
				double size = (!com.pref.phot_set.force_radius && the_psf->fwhmx > 0.0) ?
					0.5 * the_psf->fwhmx * com.pref.phot_set.auto_aperture_factor : com.pref.phot_set.aperture;
				cairo_set_dash(cr, NULL, 0, 0);
				// make the aperture slightly brighter
				cairo_set_source_rgba(cr, min(com.seq.photometry_colors[i][0] + 0.2, 1.0),
						min(com.seq.photometry_colors[i][1] + 0.2, 1.0),
						min(com.seq.photometry_colors[i][2] + 0.2, 1.0), 1.0);
				cairo_set_line_width(cr, 2.0 / dd->zoom);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, size, 0., 2. * M_PI);
				cairo_stroke(cr);

				cairo_set_source_rgba(cr, com.seq.photometry_colors[i][0],
						com.seq.photometry_colors[i][1],
						com.seq.photometry_colors[i][2], 1.0);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.pref.phot_set.inner, 0.,
						2. * M_PI);
				cairo_stroke(cr);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.pref.phot_set.outer, 0.,
						2. * M_PI);
				cairo_stroke(cr);
				cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
				cairo_set_font_size(cr, 40.0 / dd->zoom);
				cairo_move_to(cr, the_psf->xpos + com.pref.phot_set.outer + 5, the_psf->ypos);
				if (i == 0) {
					cairo_show_text(cr, "V");
				} else {
					char tmp[16];
					sprintf(tmp, "%d", i);
					cairo_show_text(cr, tmp);
				}
				cairo_stroke(cr);
			}
		}

		/* draw a cross on excluded images */
		if (com.seq.imgparam && com.seq.current >= 0 &&
				!com.seq.imgparam[com.seq.current].incl) {
			int w = dd->image_width > gfit->rx ? gfit->rx : dd->image_width;
			int h = dd->image_height > gfit->ry ? gfit->ry : dd->image_height;
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_set_source_rgb(cr, 1.0, 0.8, 0.7);
			cairo_set_line_width(cr, 2.0 / dd->zoom);
			cairo_move_to(cr, 0, 0);
			cairo_line_to(cr, w, h);
			cairo_move_to(cr, 0, h);
			cairo_line_to(cr, w, 0.0);
			cairo_stroke(cr);
		}

		/* draw preview rectangles for the manual registration */
		for (i = 0; i < PREVIEW_NB; i++) {
			if (com.seq.previewX[i] >= 0) {
				int textX, textY;
				gchar *text;
				cairo_set_line_width(cr, 1.0 / dd->zoom);
				cairo_set_source_rgb(cr, 0.1, 0.6, 0.0);
				cairo_rectangle(cr,
						com.seq.previewX[i] - com.seq.previewW[i] / 2,
						com.seq.previewY[i] - com.seq.previewH[i] / 2,
						com.seq.previewW[i], com.seq.previewH[i]);
				cairo_stroke(cr);

				textX = com.seq.previewX[i] - com.seq.previewW[i] / 2;
				textX += 0.1 * com.seq.previewW[i];

				textY = com.seq.previewY[i] - com.seq.previewH[i] / 2;
				textY += 0.1 * com.seq.previewH[i];

				text = g_strdup_printf("%d", i + 1);

				cairo_set_font_size(cr, 12.0 / dd->zoom);
				cairo_move_to(cr, textX, textY);
				cairo_show_text(cr, text);
				cairo_stroke(cr);
				g_free(text);
			}
		}
	}
}

/* Draw alignment-point boxes from the cached MPP run (com.mpp_run).
 * Called from redraw_drawingarea right after draw_brg_boxes so APs sit
 * in the same overlay stratum as bgext samples.
 *
 * AP records hold mean-frame coordinates. When gfit shows the mean ref
 * frame (== analyze just ran or user re-clicked it), no shift compensation
 * is needed; otherwise, the user is viewing a source frame and we apply
 * the per-frame offset (intersection - global_shift) so AP boxes track
 * the same physical features as they slew through the sequence.
 *
 * The "showing ref frame" test compares gfit's dims to the run's mean
 * frame dims. If global shifts are zero across the sequence, intersection
 * == frame dims and the test gives a false positive — but then dx,dy = 0
 * anyway, so the rendered result is still correct.
 *
 * Hovered AP (mpp_ap_editor_get_hover_idx) is drawn in orange with a
 * thicker line so the user knows which AP a click would affect. */
static void draw_mpp_aps(const draw_data_t* dd) {
	mpp_run_t *run = mpp_get_cached_run();
	if (!run || !run->aps || run->aps->count <= 0 || !run->cfg) return;
	/* AP coordinates live in the run's mean-frame space. Show them
	 * when gfit's dimensions match a known coordinate system in the
	 * run — either the mean frame itself (the Analyse-painted ref
	 * image, where ap.{x,y} maps directly), or a sequence frame at
	 * the run's frame_cols/_rows (where per-frame shift compensates).
	 * Hide otherwise (stacking result at non-mean-frame dims, an
	 * unrelated single image the user has loaded, etc.). The check
	 * is dim-based rather than com.seq.current-based so APs reappear
	 * after re-Analyse even when com.seq.current is still RESULT_IMAGE
	 * from a prior stack run. */
	const gboolean showing_ref = (gfit
	    && (int) gfit->rx == run->mean_frame_cols
	    && (int) gfit->ry == run->mean_frame_rows);
	const gboolean showing_seq_frame = (gfit
	    && com.seq.current >= 0
	    && com.seq.current < run->num_frames
	    && (int) gfit->rx == run->frame_cols
	    && (int) gfit->ry == run->frame_rows
	    && run->global_shifts
	    && sequence_is_loaded());
	if (!showing_ref && !showing_seq_frame) return;

	const int hb = run->cfg->alignment_points_half_box_width;
	if (hb <= 0) return;
	const int side = 2 * hb;

	int dx = 0, dy = 0;
	if (showing_seq_frame) {
		const int i = com.seq.current;
		/* mean_frame row r maps to frame row
		 *   r + intersection[0] - global_shifts[i]
		 * (matches mpp::offsets_from_run). Adding `dy` to ap->y
		 * gives the AP's pdata-row position on frame i. */
		dy = run->intersection[0] - run->global_shifts[2 * i + 0];
		dx = run->intersection[2] - run->global_shifts[2 * i + 1];
	}

	/* ap->y is a pdata-row index on the mean_frame (Siril's pdata
	 * convention has row 0 at the bottom of the displayed scene for
	 * both native FITS and SER-after-load-time-flip). The Cairo
	 * overlay context has y=0 at the TOP of the surface, so the row
	 * index has to be flipped to draw at the correct vertical position.
	 * remap_all_vports already does the same flip when filling
	 * view->buf from pdata — that's what makes the image itself display
	 * right-side-up; overlays in pdata-row coords have to redo it. */
	const int H = (int) gfit->ry;
	const int hover = mpp_ap_editor_get_hover_idx();
	for (int i = 0; i < run->aps->count; ++i) {
		const mpp_ap_record_t *ap = &run->aps->records[i];
		if (i == hover) {
			cairo_set_line_width(dd->cr, 2.0 / dd->zoom);
			cairo_set_source_rgba(dd->cr, 1.0, 0.5, 0.0, 1.0);   /* orange */
		} else {
			cairo_set_line_width(dd->cr, 1.0 / dd->zoom);
			cairo_set_source_rgba(dd->cr, 1.0, 1.0, 0.0, 0.7);   /* yellow */
		}
		const int box_y_top = (H - 1) - (ap->y + dy) - hb;
		cairo_rectangle(dd->cr, ap->x - hb + dx, box_y_top, side, side);
		cairo_stroke(dd->cr);
	}

	/* Diagnostic: only when the shift viewer is open do we overlay
	 * per-AP shift arrows. Drawing them at all times when shift data
	 * is cached would clutter the normal review workflow, especially
	 * for users without sub-pixel-precision motivation to interpret
	 * tiny vectors. The viewer's scale multiplier is the user's
	 * opt-in to magnify shifts.
	 *
	 * Green = Stage B converged, red = fell back to zero shift. */
	if (mpp_shift_viewer_is_open() && run->shifts && run->shifts->shifts) {
		const int fv = mpp_shift_viewer_get_frame();
		const double scale = mpp_shift_viewer_get_scale();
		const int M = run->shifts->num_aps;
		if (fv >= 0 && fv < run->shifts->num_frames && M == run->aps->count) {
			cairo_set_line_width(dd->cr, 1.2 / dd->zoom);
			const double dot_r = 1.8 / dd->zoom;
			for (int a = 0; a < M; ++a) {
				const mpp_ap_record_t *ap = &run->aps->records[a];
				const size_t k = (size_t) fv * M + a;
				const double sdy = run->shifts->shifts[2 * k + 0];
				const double sdx = run->shifts->shifts[2 * k + 1];
				const uint8_t ok = run->shifts->success[k];
				if (ok)
					cairo_set_source_rgba(dd->cr, 0.3, 1.0, 0.3, 0.9);
				else
					cairo_set_source_rgba(dd->cr, 1.0, 0.2, 0.2, 0.9);
				const double x0 = ap->x + dx;
				/* Same pdata-row→display-y flip as the box draw above;
				 * sdy is a pdata-row delta so its display contribution
				 * is also negated. */
				const double y0 = (double)(H - 1) - (double)(ap->y + dy);
				const double x1 = x0 + sdx * scale;
				const double y1 = y0 - sdy * scale;
				cairo_move_to(dd->cr, x0, y0);
				cairo_line_to(dd->cr, x1, y1);
				cairo_stroke(dd->cr);
				/* Arrowhead: filled triangle at (x1, y1) pointing along
				 * the arrow direction. Skip when shift is sub-pixel
				 * tiny (would be a degenerate triangle anyway). */
				const double mag = hypot(x1 - x0, y1 - y0);
				if (mag > 1.5 / dd->zoom) {
					const double ux = (x1 - x0) / mag;
					const double uy = (y1 - y0) / mag;
					const double head = 7.0 / dd->zoom;
					cairo_move_to(dd->cr, x1, y1);
					cairo_line_to(dd->cr, x1 - head * ux - 0.5 * head * uy,
					                       y1 - head * uy + 0.5 * head * ux);
					cairo_line_to(dd->cr, x1 - head * ux + 0.5 * head * uy,
					                       y1 - head * uy - 0.5 * head * ux);
					cairo_close_path(dd->cr);
					cairo_fill(dd->cr);
				}
				/* Always paint a small filled dot at the AP centre so
				 * zero/sub-pixel-shift APs are still visible (and the
				 * success/failure colour distinguishes converged-with-
				 * zero-shift from failed-fallback-to-zero). */
				cairo_arc(dd->cr, x0, y0, dot_r, 0, 2 * M_PI);
				cairo_fill(dd->cr);
			}
		}
	}
}

static void draw_brg_boxes(const draw_data_t* dd) {
	GSList *list;
	GdkRGBA gdk_color;
	for (list = com.grad_samples; list; list = list->next) {
		background_sample *sample = (background_sample *)list->data;
		if (sample && background_sample_is_valid(sample)) {
			int radius = (int) (background_sample_get_size(sample) / 2);
			point position = background_sample_get_position(sample);
			cairo_set_line_width(dd->cr, 1.5 / dd->zoom);
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_bkg_samples);
			cairo_set_source_rgba(dd->cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			cairo_rectangle(dd->cr, position.x - radius - 1, position.y - radius,
					radius * 2, radius * 2);
			cairo_stroke(dd->cr);
		}
	}
}

static void draw_in_progress_poly(const draw_data_t* dd) {
	if (!gui.drawing_polypoints) return;

	cairo_t *cr = dd->cr;
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	gdk_cairo_set_source_rgba(cr, &gui.poly_ink);

	GSList *list = gui.drawing_polypoints;
	const point* start = (point*) list->data;
	cairo_move_to(cr, start->x + 0.5, start->y + 0.5);

	// Build the complete path first
	for (list = list->next; list; list = list->next) {
		const point* position = (point*) list->data;
		cairo_line_to(cr, position->x + 0.5, position->y + 0.5);
	}

	// Now stroke the entire path at once
	cairo_stroke(cr);
}

static void draw_compass(const draw_data_t* dd) {
	int pos = com.pref.gui.position_compass;
	if (!pos) return; // User chose None
	fits *fit = gfit;
	cairo_t *cr = dd->cr;
	cairo_set_line_width(cr, 3.0 / dd->zoom);
	double ra0, dec0;
	double xN, yN, xE, yE;

	double xpos = fit->rx * 0.5;
	double ypos = fit->ry * 0.5;
	pix2wcs(fit, xpos, ypos, &ra0, &dec0);
	if (ra0 == -1) return; // checks implicitly that wcslib member exists
	if (90. - fabs(dec0) < 2.78e-3) return;// center is less than 10"off from a pole
	double len = (double)fit->ry / 20.;
	wcs2pix(fit, ra0, dec0 + 0.1, &xN, &yN);
	wcs2pix(fit, ra0 + 0.1, dec0, &xE, &yE);
	double angleN = -atan2(yN - ypos, xN - xpos);
	double angleE = -atan2(yE - ypos, xE - xpos);

	cairo_set_font_size(cr, len / 3);

	double xdraw, ydraw;
	double pos_values[5][2] = { { 0.5, 0.5 }, { 0.1, 0.1 }, { 0.9, 0.1 }, { 0.1, 0.9 }, { 0.9, 0.9 } };
	xdraw = pos_values[pos - 1][0] * fit->rx;
	ydraw = pos_values[pos - 1][1] * fit->ry;

	/* draw north line and filled-arrow*/
	cairo_set_source_rgba(cr, 1., 0., 0., 1.0);
	cairo_save(cr); // save the original transform
	cairo_translate(cr, xdraw, ydraw);
	cairo_rotate(cr, angleN);
	cairo_move_to(cr, 0., 0.);
	cairo_line_to(cr, len, 0.);
	cairo_stroke(cr);
	cairo_line_to(cr, 0.75 * len, -0.15 * len);
	cairo_line_to(cr, 0.75 * len, +0.15 * len);
	cairo_line_to(cr, len, 0.);
	cairo_fill(cr);
	cairo_move_to(cr, len * 1.3, 0.1 * len);
	cairo_rotate(cr, -angleN);
	cairo_show_text(cr, "N");
	cairo_restore(cr); // restore the original transform

	/* draw east line */
	if (dd->neg_view)
		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
	else
		cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);

	cairo_save(cr); // save the original transform
	cairo_translate(cr, xdraw, ydraw);
	cairo_rotate(cr, angleE);
	cairo_move_to(cr, 0., 0.);
	cairo_line_to(cr, len / 2.0, 0.);
	cairo_stroke(cr);
	cairo_move_to(cr, (len / 2) * 2.0, -0.1 * len);
	cairo_rotate(cr, -angleE);
	cairo_show_text(cr, "E");
	cairo_stroke(cr);
	cairo_restore(cr); // restore the original transform
}

static label_point *new_label_point(double height, const double *pix1, const double *pix2, const double *world, gboolean isRA, int border) {
	label_point *pt = g_new(label_point, 1);

	pt->x = pix1[0];
	pt->y = height - pix1[1];
	pt->ra = world[0];
	pt->dec = world[1];
	pt->angle = -atan2(pix2[1] - pix1[1], pix2[0] - pix1[0]);
	pt->isRA = isRA;
	pt->border = border;

	return pt;
}

static int has_pole(fits *fit) {
	if (!wcs2pix(fit, 0., 90., NULL, NULL))
		return 1;
	if (!wcs2pix(fit, 0., -90., NULL, NULL))
		return -1;
	return 0;
}

static gboolean get_line_intersection(double p0_x, double p0_y, double p1_x,
		double p1_y, double p2_x, double p2_y, double p3_x, double p3_y,
		double *i_x, double *i_y) {
	double s1_x, s1_y, s2_x, s2_y;
	s1_x = p1_x - p0_x;	s1_y = p1_y - p0_y;
	s2_x = p3_x - p2_x;	s2_y = p3_y - p2_y;
	double det = -s2_x * s1_y + s1_x * s2_y;
	if (fabs(det) < DBL_EPSILON)
		return FALSE;
	double s, t;
	s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / det;
	t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / det;
	if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
		// Collision detected
		if (i_x != NULL)
			*i_x = p0_x + (t * s1_x);
		if (i_y != NULL)
			*i_y = p0_y + (t * s1_y);
		return TRUE;
	}
	return FALSE; // No collision
}

static gint border_compare(const label_point *a, const label_point *b) {
	if (a->border > b->border) return 1;
	if (a->border < b->border) return -1;
	return 0;
}

static double ra_values[] = { 45, 30, 15, 10, 7.5, 5, 3.75, 2.5, 1.5, 1.25, 1, 3. / 4., 1.
		/ 2., 1. / 4., 1. / 6., 1. / 8., 1. / 12., 1. / 16., 1. / 24., 1. / 40., 1. / 48. };


static void draw_wcs_grid(const draw_data_t* dd) {
	if (!gui.show_wcs_grid) return;
	fits *fit = gfit;
	if (!has_wcs(fit)) return;
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1. / dd->zoom);
	cairo_set_font_size(cr, 12.0 / dd->zoom);
	double ra0, dec0;
	GList *ptlist = NULL;
	double world[2], pix[2], pix2[2], img[2];
	double phi, theta;
	int status;

	double width = (double) fit->rx;
	double height = (double) fit->ry;
	cairo_rectangle(cr, 0., 0., width, height); // to clip the grid
	cairo_clip(cr);
	/* get ra and dec of center of the image */
	center2wcs(fit, &ra0, &dec0);
	if (ra0 == -1.) return;
	dec0 *= (M_PI / 180.0);
	ra0  *= (M_PI / 180.0);
	double range = get_wcs_image_resolution(fit) * sqrt(pow((width / 2.0), 2) + pow((height / 2.0), 2)); // range in degrees, FROM CENTER
	double step;

	/* Compute borders in pixel for tags*/
	const double pixbox[5][2] = { { 0., 0. }, { width, 0. }, { width, height }, { 0., height }, { 0., 0. } };
	const double pixval[4] = { 0., width, height, 0. }; // bottom, right, top, left with ref bottom left
	int pixtype[4] = { 1, 0, 1, 0 }; // y, x, y, x
	int polesign = has_pole(fit);

	/* calculate DEC step size */
	if (range > 16.0) {
		step = 8.; //step DEC 08:00
	} else if (range > 8.0) {
		step = 4.; // step DEC 04:00
	} else if (range > 4.0) { // image FOV about >2*4/sqrt(2) so >5 degrees
		step = 2.; // step DEC 02:00
	} else if (range > 2.0) {
		step = 1.; // step DEC 01:00
	} else if (range > 1.0) {
		step = 0.5; // step DEC 00:30
	} else if (range > 0.5) {
		step = 0.25; // step DEC 00:15
	} else if (range > 0.3) {
		step = 1. / 6.; // 0.166666, step DEC 00:10
	} else {
		step = 1. / 12.; // step DEC 00:05
	}

	// calculate RA step size
	double step2 = min(45, step / (cos(dec0) + 0.000001)); // exact value for stepRA, but not well rounded
	int iter = 0;
	double stepRA;
	do { // select nice rounded values for ra_step
		stepRA = ra_values[iter];
		iter++;
	} while ((stepRA >= step2) && (iter < G_N_ELEMENTS(ra_values))); // repeat until compatible value is found in ra_values
	if (polesign) stepRA = 45.;

	// round image centers
	double centra = stepRA * round(ra0 * 180 / (M_PI * stepRA)); // rounded image centers
	double centdec = step * round(dec0 * 180 / (M_PI * step));

	// plot DEC grid
	cairo_set_source_rgb(cr, 0.8, 0.0, 0.0);
	double di = (polesign) ? 0. : centra - 6 * stepRA;
	do { // dec lines
		double dj = max(centdec - 6 * step, -90);
		do {
			double xa, ya, xb, yb, x1, x2, y1, y2;

			wcs2pix(fit, di, dj, &xa, &ya);
			x1 = round(xa - 1);
			y1 = round(height - ya);

			wcs2pix(fit, di, (dj + step), &xb, &yb);
			x2 = round(xb - 1);
			y2 = round(height - yb);

			if (((x1 >= 0) && (y1 >= 0) && (x1 < width) && (y1 < height))
					|| ((x2 >= 0) && (y2 >= 0) && (x2 < width) && (y2 < height))) {
				cairo_move_to(cr, x1, y1);
				cairo_line_to(cr, x2, y2);
				cairo_stroke(cr);
			}
			// check crossing
			if (!(((xa >= 0) && (ya >= 0) && (xa < width) && (ya < height))
						&& ((xb >= 0) && (yb >= 0) && (xb < width) && (yb < height)))) {
				for (int k = 0; k < 4; k ++) {
					if (get_line_intersection(xa, ya, xb, yb, pixbox[k][0], pixbox[k][1], pixbox[k+1][0], pixbox[k+1][1], NULL, NULL)) {
						world[0] = di;
						pix[pixtype[k]] = pixval[k];
						double latspan[2] = {dj, dj+step};
						status = wcsmix(fit->keywords.wcslib, pixtype[k], 1, latspan, 1.0, 0, world, &phi, &theta, img, pix);
						if(!status) {
							wcs2pix(fit, world[0], world[1] + 0.1, &pix2[0], &pix2[1]);
							ptlist = g_list_append(ptlist, new_label_point(height, pix, pix2, world, TRUE, k));
						}
						break;
					}
				}
			}
			dj = dj + step;
		} while (dj <= min(centdec + 6 * step, 90.));
		di = di + stepRA;
	} while (di <= ((polesign) ? 360. : centra + 6 * stepRA));

	// plot RA grid
	cairo_set_source_rgb(cr, 0.0, 0.5, 1.0);
	double dj = max(centdec - step * 6, -90);
	do { // ra lines
		di = (polesign) ? 0. : centra - 6 * stepRA;
		do {
			double xa, ya, xb, yb, x1, x2, y1, y2;

			wcs2pix(fit, di, dj, &xa, &ya);
			x1 = round(xa - 1);
			y1 = round(height - ya);

			wcs2pix(fit, (di + step), dj, &xb, &yb);
			x2 = round(xb - 1);
			y2 = round(height - yb);

			if (((x1 >= 0) && (y1 >= 0) && (x1 < width) && (y1 < height))
					|| ((x2 >= 0) && (y2 >= 0) && (x2 < width) && (y2 < height))) {
				cairo_move_to(cr, x1, y1);
				cairo_line_to(cr, x2, y2);
				cairo_stroke(cr);
			}
				// check crossing
			if (!(((xa >= 0) && (ya >= 0) && (xa < width) && (ya < height))
				&& ((xb >= 0) && (yb >= 0) && (xb < width) && (yb < height)))) {
				for (int k = 0; k < 4; k ++) {
					if (get_line_intersection(xa, ya, xb, yb, pixbox[k][0], pixbox[k][1], pixbox[k+1][0], pixbox[k+1][1], NULL, NULL)) {
						world[1] = dj;
						pix[pixtype[k]] = pixval[k];
						double lngspan[2] = {di, di+step};
						status = wcsmix(fit->keywords.wcslib, pixtype[k], 2, lngspan, 1.0, 0, world, &phi, &theta, img, pix);
						if(!status) {
							wcs2pix(fit, world[0] + 0.1, world[1], &pix2[0], &pix2[1]);
							ptlist = g_list_append(ptlist, new_label_point(height, pix, pix2, world, FALSE, k));
						}
						break;
					}
				}
			}
			di = di + step;
		} while (di <= ((polesign) ? 360. : centra + 6 * stepRA));
		dj = dj + step;
	} while (dj <= min(centdec + step * 6, 90));

	// Add crossings labels
	ptlist = g_list_sort(ptlist, (GCompareFunc) border_compare); // sort potential tags by increasing border number
	if (dd->neg_view) {
		cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
	} else {
		cairo_set_source_rgb(cr, 0.8, 0.8, 0.8);
	}
	GSList *existingtags = NULL;
	gchar *RAfmt = (stepRA < 1./4.) ? "%02dh%02dm%02ds" : "%02dh%02dm";
	for (GList *l = ptlist; l != NULL; l = l->next) {
		// getting the label
		SirilWorldCS *world_cs;
		label_point *pt = (label_point*) l->data;
		world_cs = siril_world_cs_new_from_a_d(pt->ra, pt->dec);
		if (world_cs) {
			gchar *tag = (pt->isRA) ? siril_world_cs_alpha_format(world_cs, RAfmt) : siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'");
			siril_world_cs_unref(world_cs);
			if (!g_slist_find_custom(existingtags, tag, (GCompareFunc) strcompare)) { // this tag has already been used - skipping
				existingtags = g_slist_append(existingtags, (gpointer) tag);
				cairo_text_extents_t te1, te2;
				cairo_text_extents(cr, tag, &te1); // getting the dimensions of the textbox
				cairo_save(cr); // save the orginal transform
				cairo_translate(cr, pt->x, pt->y);
				// add pi for angles larger than +/- pi/2
				if (pt->angle > M_PI_2)
					pt->angle -= M_PI;
				if (pt->angle < -M_PI_2)
					pt->angle += M_PI;
				double dx = 0., dy = 0.;
				switch (pt->border) { // shift to get back in the image
				case 0: // bottom
					if (pt->angle > 0.)
						dx -= te1.x_advance;
					break;
				case 1: // right
					dx -= te1.x_advance;
					if (pt->angle > 0.)
						dy += te1.height;
					break;
				case 2: // top
					dy += te1.height;
					if (pt->angle < 0.)
						dx -= te1.x_advance;
					break;
				case 3: // left
					if (pt->angle < 0.)
						dy += te1.height;
					break;
				default:
					break;
				}
				cairo_rotate(cr, pt->angle);
				cairo_move_to(cr, dx, dy);
				cairo_text_extents(cr, tag, &te2);
				cairo_show_text(cr, tag);
				cairo_stroke(cr);
				cairo_restore(cr); // restore the orginal transform
			} else {
				g_free(tag);
			}
		}
	}
	g_list_free_full(ptlist, (GDestroyNotify) g_free);
	g_slist_free_full(existingtags, (GDestroyNotify) g_free);

	draw_compass(dd);
}

static void draw_wcs_disto(const draw_data_t* dd) {
	if (!gui.show_wcs_disto) return;
	fits *fit = gfit;
	if (!has_wcs(fit) || !fit->keywords.wcslib->lin.dispre) return; // no platesolve or no distortions
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 3. / dd->zoom);
	cairo_set_source_rgb(cr, 0.8, 0.0, 0.0);

	int nbpoints = 20;
	double radius = (dd->zoom < 1. )?  3. : 3. / dd->zoom;

	double paceX = (double)fit->rx / (double)(nbpoints - 1);
	double paceY = (double)fit->ry / (double)(nbpoints - 1);
	double currX = 0.5; // coords in fits/wcs conventions start at (0.5, 0.5) for the bottom left corner
	double currY = 0.5;
	for (int i = 0; i < nbpoints; i++) {
		currY = 0.5;
		for (int j = 0; j < nbpoints; j++) {
			double rawcrd[2] = { currX, currY };
			double discrd[2] = {0., 0.};
			int status = disp2x(fit->keywords.wcslib->lin.dispre, rawcrd, discrd);
			if (!status) {
				double disX = (discrd[0] - rawcrd[0]) * 5. +  rawcrd[0];
				double disY = (discrd[1] - rawcrd[1]) * 5. +  rawcrd[1];
				double startX, startY, endX, endY;
				fits_to_display(currX, currY, &startX, &startY, fit->ry);
				fits_to_display(disX, disY, &endX, &endY, fit->ry);
				cairo_arc(cr, startX, startY, radius,  0., 2. * M_PI);
				cairo_fill(cr);
				cairo_move_to(cr, startX, startY);
				cairo_line_to(cr, endX, endY);
				cairo_stroke(cr);
			}
			currY += paceY;
		}
		currX += paceX;
	}
}

static gdouble x_circle(gdouble x, gdouble radius, gdouble angle) {
	return x + radius * cos(angle);
}

static gdouble y_circle(gdouble y, gdouble radius, gdouble angle) {
	return y + radius * sin(angle);
}

static void draw_annotates(const draw_data_t* dd) {
	if (!com.found_object) return;
	gdouble resolution = get_wcs_image_resolution(gfit);
	if (resolution <= 0) return;
	double width = (double) gfit->rx;
	double height = (double) gfit->ry;
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);

	cairo_set_line_width(cr, 1.0 / dd->zoom);
	cairo_rectangle(cr, 0., 0., width, height); // to clip the grid
	cairo_clip(cr);

	for (GSList *list = com.found_object; list; list = list->next) {
		CatalogObjects *object = (CatalogObjects *)list->data;
		gdouble radius = get_catalogue_object_radius(object);
		gdouble x = get_catalogue_object_x(object);
		gdouble y = get_catalogue_object_y(object);
		gdouble x1 = get_catalogue_object_x1(object);
		gdouble y1 = get_catalogue_object_y1(object);
		gchar *code = get_catalogue_object_code_pretty(object);
		guint catalog = get_catalogue_object_cat(object);
		gboolean revert = FALSE;
		double angle = ANGLE_TOP;
		double addoffset = 0.;
		GdkRGBA gdk_color;

		switch (catalog) {
		case CAT_AN_USER_DSO:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_dso_annotations);
			cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			break;
		case CAT_AN_USER_SSO:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_sso_annotations);
			cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			break;
		case CAT_AN_USER_TEMP:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_tmp_annotations);
			cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			revert = TRUE;
			angle = ANGLE_BOT;
			break;
		default:
		case 0:
			gdk_rgba_parse(&gdk_color, com.pref.gui.config_colors.color_std_annotations);
			if (dd->neg_view) {
				cairo_set_source_rgba(cr, 1.0 - gdk_color.red, 1.0 - gdk_color.green, 1.0 - gdk_color.blue, gdk_color.alpha);
			} else {
				cairo_set_source_rgba(cr, gdk_color.red, gdk_color.green, gdk_color.blue, gdk_color.alpha);
			}
			break;
		}

		radius = radius / resolution / 60.0;
		// radius now in pixels

		point offset = {5., revert ? 5. : -5.};
		if (catalog == CAT_AN_CONST) { // constellation line
			cairo_move_to(cr, x, y);
			cairo_line_to(cr, x1, y1);
			cairo_stroke(cr);
		} else if (radius < 0 || catalog == CAT_AN_CONST_NAME) {
			// objects we don't have an accurate location (LdN, Sh2)
		} else if (radius > 5) {
			cairo_arc(cr, x, y, radius, 0., 2. * M_PI);
			cairo_stroke(cr);
			if (code) {
				cairo_move_to(cr, x_circle(x, radius, angle), y_circle(y, radius, angle));
				offset.x = x_circle(x, radius * 1.3, angle) - x;
				offset.y = y_circle(y, radius * 1.3, angle) - y;
				cairo_line_to(cr, offset.x + x, offset.y + y);
				cairo_stroke(cr);
			}
		} else {
			/* it is punctual */
			cairo_move_to(cr, x, y - 15);
			cairo_line_to(cr, x, y - 5);
			cairo_stroke(cr);
			cairo_move_to(cr, x, y + 15);
			cairo_line_to(cr, x, y + 5);
			cairo_stroke(cr);
			cairo_move_to(cr, x - 15, y);
			cairo_line_to(cr, x - 5, y);
			cairo_stroke(cr);
			cairo_move_to(cr, x + 15, y);
			cairo_line_to(cr, x + 5, y);
			cairo_stroke(cr);
		}
		if (code) {
			gchar *name = code;
			gdouble size = 18 * (com.pref.gui.font_scale / 100.0);
			cairo_select_font_face(cr, "Liberation Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
			cairo_set_font_size(cr, size / dd->zoom);
			if (revert) {
				cairo_text_extents_t te;
				cairo_text_extents(cr, name, &te); // getting the dimensions of the textbox
				addoffset = te.height;
			}
			cairo_move_to(cr, x + offset.x, y + offset.y + addoffset);
			cairo_show_text(cr, name);
			cairo_stroke(cr);
		}

	}
}

static void draw_rgb_centers(const draw_data_t* dd) {
	if (!gui.comp_layer_centering) return;
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);

	double red = gui.comp_layer_centering->saturated_color.red;
	double green = gui.comp_layer_centering->saturated_color.green;
	double blue = gui.comp_layer_centering->saturated_color.blue;
	cairo_set_source_rgb(cr, red, green, blue);
	cairo_set_line_width(cr, 2.0 / dd->zoom);

	double size = 10. / dd->zoom;
	cairo_arc(cr, gui.comp_layer_centering->center.x, gui.comp_layer_centering->center.y, size * 0.5, 0., 2. * M_PI);
	cairo_stroke(cr);
}

static void draw_analysis(const draw_data_t* dd) {
	if (com.tilt) {
		cairo_t *cr = dd->cr;
		cairo_set_dash(cr, NULL, 0, 0);

		cairo_set_source_rgb(cr, 1.0, 0.8, 0.7);
		cairo_set_line_width(cr, 2.0 / dd->zoom);
		cairo_move_to(cr, com.tilt->pt[0].x, com.tilt->pt[0].y);
		cairo_line_to(cr, com.tilt->pt[1].x, com.tilt->pt[1].y);
		cairo_line_to(cr, com.tilt->pt[2].x, com.tilt->pt[2].y);
		cairo_line_to(cr, com.tilt->pt[3].x, com.tilt->pt[3].y);
		cairo_line_to(cr, com.tilt->pt[1].x, com.tilt->pt[1].y);
		cairo_move_to(cr, com.tilt->pt[3].x, com.tilt->pt[3].y);
		cairo_line_to(cr, com.tilt->pt[0].x, com.tilt->pt[0].y);
		cairo_line_to(cr, com.tilt->pt[2].x, com.tilt->pt[2].y);
		cairo_stroke(cr);

		/* draw text */
		cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
		int size = 20.0 / dd->zoom;
		cairo_set_font_size(cr, size);

		/* fwhm 1 */
		gchar *str = g_strdup_printf("%.2f", com.tilt->fwhm[0]);
		cairo_move_to(cr, com.tilt->pt[0].x, com.tilt->pt[0].y - size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm 2 */
		str = g_strdup_printf("%.2f", com.tilt->fwhm[1]);
		cairo_move_to(cr, com.tilt->pt[1].x, com.tilt->pt[1].y - size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm 3 */
		str = g_strdup_printf("%.2f", com.tilt->fwhm[2]);
		cairo_move_to(cr, com.tilt->pt[2].x, com.tilt->pt[2].y + size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm 4 */
		str = g_strdup_printf("%.2f", com.tilt->fwhm[3]);
		cairo_move_to(cr, com.tilt->pt[3].x, com.tilt->pt[3].y + size);
		cairo_show_text(cr, str);
		g_free(str);
		/* fwhm center */
		str = g_strdup_printf("%.2f", com.tilt->fwhm_centre);
		cairo_move_to(cr, gfit->rx / 2.0, (gfit->ry / 2.0) + size);
		cairo_show_text(cr, str);
		cairo_stroke(cr);
		g_free(str);
	}
}

static void draw_regframe(const draw_data_t* dd) {
	if (com.script || com.headless) return;
	if (!sequence_is_loaded()) return;
	if (com.seq.current == RESULT_IMAGE) return;
	if (!drawframe)
		image_display_init_statics();
	if (!gtk_toggle_button_get_active(drawframe)) return;
	int activelayer = gtk_combo_box_get_active(seqcombo);
	if (!layer_has_registration(&com.seq, activelayer)) return;
	if (com.seq.reg_invalidated) return;
	transformation_type min, max;
	guess_transform_from_seq(&com.seq, activelayer, &min, &max, FALSE);
	if (max <= IDENTITY_TRANSFORMATION) return;

	if (guess_transform_from_H(com.seq.regparam[activelayer][com.seq.reference_image].H) == NULL_TRANSFORMATION ||
			guess_transform_from_H(com.seq.regparam[activelayer][com.seq.current].H) == NULL_TRANSFORMATION)
		return; // reference or current image H matrix is null matrix

	regframe framing = { 0 };
	framing.pt[0].x = 0.;
	framing.pt[0].y = 0.;
	framing.pt[1].x = (double)com.seq.imgparam[com.seq.reference_image].rx;
	framing.pt[1].y = 0.;
	framing.pt[2].x = (double)com.seq.imgparam[com.seq.reference_image].rx;
	framing.pt[2].y = (double)com.seq.imgparam[com.seq.reference_image].ry;
	framing.pt[3].x = 0.;
	framing.pt[3].y = (double)com.seq.imgparam[com.seq.reference_image].ry;
	double cogx = 0., cogy = 0., cx, cy;
	for (int i = 0; i < 4; i++) {
		cvTransfPoint(&framing.pt[i].x, &framing.pt[i].y, com.seq.regparam[activelayer][com.seq.reference_image].H, com.seq.regparam[activelayer][com.seq.current].H, 1.);
		cogx += framing.pt[i].x;
		cogy += framing.pt[i].y;
	}
	cogx *= 0.25;
	cogy *= 0.25;
	cx = (com.seq.is_variable) ? (double)com.seq.imgparam[com.seq.current].rx * 0.5 : (double)com.seq.rx * 0.5;
	cy = (com.seq.is_variable) ? (double)com.seq.imgparam[com.seq.current].ry * 0.5 : (double)com.seq.ry * 0.5;

	cairo_t *cr = dd->cr;
	double size = 10. / dd->zoom;
	cairo_set_dash(cr, NULL, 0, 0);
	gboolean has_disto = seq_has_any_distortion(&com.seq);
	if (max <= SHIFT_TRANSFORMATION)
		cairo_set_source_rgb(cr, 0.0, 0.5, 1.0);
	else {
		if (has_disto)
			cairo_set_source_rgb(cr, 0., 1., 0.5);
		else
			cairo_set_source_rgb(cr, 1., 0., 0.);
	}

	cairo_set_line_width(cr, 2.0 / dd->zoom);
	// reference origin
	cairo_arc(cr, framing.pt[0].x, framing.pt[0].y, size * 0.5, 0., 2. * M_PI);
	cairo_stroke_preserve(cr);
	cairo_fill(cr);
	// reference frame
	cairo_move_to(cr, framing.pt[0].x, framing.pt[0].y);
	cairo_line_to(cr, framing.pt[1].x, framing.pt[1].y);
	cairo_line_to(cr, framing.pt[2].x, framing.pt[2].y);
	cairo_line_to(cr, framing.pt[3].x, framing.pt[3].y);
	cairo_line_to(cr, framing.pt[0].x, framing.pt[0].y);
	cairo_stroke(cr);

	// reference center
	cairo_arc(cr, cogx, cogy, size * 0.5, 0., 2. * M_PI);
	cairo_stroke(cr);
	// current center
	cairo_set_source_rgb(cr, 0., 1., 0.);
	cairo_move_to(cr, cx - size * 0.5, cy);
	cairo_rel_line_to(cr, size, 0.);
	cairo_stroke(cr);
	cairo_move_to(cr, cx, cy - size * 0.5);
	cairo_rel_line_to(cr, 0., size);
	cairo_stroke(cr);
}

void initialize_image_display() {
	int i;
	siril_log_debug("HD AutoStretch bitdepth: %d\n", com.pref.hd_bitdepth);
	gui.hd_remap_max = 1 << (guint) com.pref.hd_bitdepth;
	for (i = 0; i < MAXGRAYVPORT; i++) {
		memset(gui.remap_index[i], 0, sizeof(gui.remap_index[i]));
		last_pente = 0.f;
		last_mode = HISTEQ_DISPLAY;
		// only HISTEQ mode always computes the index, it's a good initializer here
	}
	cairo_matrix_init_identity(&gui.display_matrix);
	gui.surface_scale = 1.0;
}

/* this function calculates the "fit to window" zoom values, given the window
 * size in argument and the image size in gfit->
 * Should not be called before displaying the main gray window when using zoom to fit */
double get_zoom_val() {
	int window_width, window_height;
	if (gui.zoom_value > 0.)
		return gui.zoom_value;
	/* else if zoom is < 0, it means fit to window */
	window_width = gtk_widget_get_allocated_width(gui.view[RED_VPORT].drawarea);
	window_height = gtk_widget_get_allocated_height(gui.view[RED_VPORT].drawarea);
	if (gfit->rx == 0 || gfit->ry == 0 || window_height <= 1 || window_width <= 1)
		return 1.0;
	double wtmp = (double) window_width / (double) gfit->rx;
	double htmp = (double) window_height / (double) gfit->ry;
	return min(wtmp, htmp);
}

static void invalidate_image_render_cache(int vport) {
	/* the render cache is a surface containing the rendering of the image
	 * from the mapped buffers to the drawing area.
	 * If some of the parameters change, the cache must be invalidated to
	 * redraw the image: image content, image mapping, widget dimensions.
	 */
	g_mutex_lock(&gui.cairo_mutex);
	for (int i = 0; i < MAXVPORT; i++) {
		if (vport >= 0 && i != vport)
			continue;
		if (gui.view[i].disp_surface)
			cairo_surface_destroy(gui.view[i].disp_surface);
		gui.view[i].disp_surface = NULL;
		gui.view[i].view_height = -1;
		gui.view[i].view_width = -1;
	}
	g_mutex_unlock(&gui.cairo_mutex);
	//siril_log_debug("###\t\t\tcache surface invalidated\t\t\t###\n");
}

void adjust_vport_size_to_image() {
	if (com.script) return;
	double zoom = get_zoom_val();
	if (zoom <= 0.0) return;
	/* Init display matrix from current display state */
	cairo_matrix_t new_matrix;
	/*siril_log_debug("computing matrix for zoom %g and offset [%g, %g]\n",
			zoom, gui.display_offset.x, gui.display_offset.y);*/
	cairo_matrix_init(&new_matrix,
			zoom, 0, 0, zoom,
			gui.display_offset.x,
			gui.display_offset.y);
	if (memcmp(&new_matrix, &gui.display_matrix, sizeof(new_matrix))) {
		invalidate_image_render_cache(-1);
		gui.display_matrix = new_matrix;

		/* Compute the inverse display matrix used for coordinate transformation */
		gui.image_matrix = gui.display_matrix;
		cairo_matrix_invert(&gui.image_matrix);
		//siril_log_debug("  matrix changed\n");
	}
}

void copy_roi_into_gfit() {
	size_t npixels_roi = gui.roi.selection.w * gui.roi.selection.h;
	if (npixels_roi == 0 || com.script || com.python_command)
		return;
	g_rw_lock_writer_lock(&gfit->rwlock);
	size_t npixels_gfit = gfit->rx * gfit->ry;
	if (gui.roi.fit.type != gfit->type) {
		size_t roi_ndata = gui.roi.fit.rx * gui.roi.fit.ry * gui.roi.fit.naxes[2];
		if (gfit->type == DATA_FLOAT) {
			fit_replace_buffer(&gui.roi.fit, gui.roi.fit.bitpix == BYTE_IMG ? ushort8_buffer_to_float(gui.roi.fit.data, roi_ndata): ushort_buffer_to_float(gui.roi.fit.data, roi_ndata), DATA_FLOAT);
			if (is_preview_active()) {
				fits *roi_backup = get_roi_backup();
				fit_replace_buffer(roi_backup, roi_backup->bitpix == BYTE_IMG ? ushort8_buffer_to_float(roi_backup->data, roi_ndata) : ushort_buffer_to_float(roi_backup->data, roi_ndata), DATA_FLOAT);
			}
		} else {
			fit_replace_buffer(&gui.roi.fit, float_buffer_to_ushort(gui.roi.fit.fdata, roi_ndata), DATA_USHORT);
			if (gfit->bitpix == BYTE_IMG) {
				for (size_t i = 0 ; i < roi_ndata ; i++) {
					gui.roi.fit.data[i] >>= 8;
				}
				gui.roi.fit.bitpix = BYTE_IMG;
			}
			if (is_preview_active()) {
				fits *roi_backup = get_roi_backup();
				fit_replace_buffer(roi_backup, float_buffer_to_ushort(roi_backup->fdata, roi_ndata), DATA_USHORT);
				if (gfit->bitpix == BYTE_IMG) {
					for (size_t i = 0 ; i < roi_ndata ; i++) {
						roi_backup->data[i] >>= 8;
					}
					roi_backup->bitpix = BYTE_IMG;
				}
			}
		}
	}

	if (gui.roi.fit.type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif
		for (uint32_t c = 0 ; c < gui.roi.fit.naxes[2] ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				const float *rowindex = gui.roi.fit.fdata + (y * gui.roi.fit.rx) + (c * npixels_roi);
				float *destindex = gfit->fdata + (c * npixels_gfit) + ((gfit->ry - gui.roi.selection.y - y - 1) * gfit->rx) + gui.roi.selection.x;
				memcpy(destindex, rowindex, gui.roi.selection.w * sizeof(float));
			}
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif
		for (uint32_t c = 0 ; c < gui.roi.fit.naxes[2] ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				const WORD *rowindex = gui.roi.fit.data + (y * gui.roi.fit.rx) + (c * npixels_roi);
				WORD *destindex = gfit->data + (npixels_gfit * c) + ((gfit->ry - gui.roi.selection.y - y - 1) * gfit->rx) + gui.roi.selection.x;
				memcpy(destindex, rowindex, gui.roi.selection.w * sizeof(WORD));
			}
		}
	}
	g_rw_lock_writer_unlock(&gfit->rwlock);
}

void remap_all() {
	stf_computed = FALSE;
	if (gui.rendering_mode == HISTEQ_DISPLAY || gui.rendering_mode == STF_DISPLAY) {
		for (int i = 0; i < gfit->naxes[2]; i++) {
			remap(i);
		}
	} else {
		remap_all_vports();
	}
	if (gfit->naxis == 3)
		remaprgb();
}

void redraw(remap_type doremap) {
	if (com.script && !com.python_script) return;
	switch (doremap) {
		case REDRAW_OVERLAY:
			break;
		case REDRAW_IMAGE:
			invalidate_image_render_cache(-1);
			break;
		case REMAP_ALL:
			/* redraw the 9-panel mosaic dialog if needed */
			redraw_aberration_inspector();
			break;
		default:
			siril_log_debug("UNKNOWN REMAP\n\n");
	}
	request_gtk_redraw_of_cvport();
}

static gboolean redraw_idle(gpointer p) {
	redraw((remap_type)GPOINTER_TO_INT(p)); // draw stars
	return FALSE;
}

static gpointer redraw_idle_thread_func(gpointer data) {
	sample_mutex_lock();
	execute_idle_and_wait_for_it(redraw_idle, data);
	sample_mutex_unlock();
	return FALSE;
}

void queue_redraw_and_wait_for_it(remap_type doremap) {
	if (!com.script && !com.python_command && !com.headless) {
		GThread *thread = g_thread_new("redraw", redraw_idle_thread_func, GINT_TO_POINTER((int)doremap));
		g_thread_join(thread);
	}
}

gboolean redraw_mask_idle(gpointer p) {
	g_rw_lock_reader_lock(&gfit->rwlock);
	if (gfit->mask && gfit->mask->data)
		remap_mask(gfit->mask);
	g_rw_lock_reader_unlock(&gfit->rwlock);
	if (com.pref.gui.mask_tints_vports) {
		/* Reader lock released before this call: notify_gfit_data_modified()
		 * may call copy_roi_into_gfit() which acquires the writer lock. */
		notify_gfit_data_modified();
		redraw(REMAP_ALL); // need to remap all to tint the image vports correctly
	}
	return FALSE;
}

void queue_redraw_mask() {
	siril_add_idle(redraw_mask_idle, NULL);
}

void queue_redraw(remap_type doremap) {
	// request a redraw from another thread
	siril_add_idle(redraw_idle, GINT_TO_POINTER((int)doremap));
}


/* callback for GtkDrawingArea, draw event */
gboolean redraw_drawingarea(GtkWidget *widget, cairo_t *cr, gpointer data) {
	draw_data_t dd;
	static GAction *action_neg = NULL;
	if (action_neg == NULL) {
		image_display_init_statics();
		action_neg = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
	}
	// we need to identify which vport is being redrawn
	dd.vport = match_drawing_area_widget(widget, TRUE);
	if (dd.vport == -1) {
		fprintf(stderr, "Could not find the vport for the draw callback\n");
		return TRUE;
	}

	/* While generic_image_worker is running the remap buffers (gfit pixel
	 * data) are stale.  Repaint from the cached display surface so the
	 * previous correct frame stays visible, avoiding a grey flash.
	 *
	 * If disp_surface was invalidated - e.g. the user changed zoom level or
	 * resized the window during a long operation - fall through to a normal
	 * draw.  Zoom and resize do not require a remap: gui.view[].buf already
	 * holds the full-resolution remapped image and is safe to read from the
	 * GTK thread.  draw_main_image() will re-render buf into a fresh
	 * disp_surface at the new viewport geometry without touching gfit. */
	if (g_atomic_int_get(&gui.suppress_drawarea_redraw)) {
		g_mutex_lock(&gui.cairo_mutex);
		cairo_surface_t *cached = gui.view[dd.vport].disp_surface;
		if (cached) {
			cairo_set_source_surface(cr, cached, 0, 0);
			cairo_paint(cr);
			g_mutex_unlock(&gui.cairo_mutex);
			return FALSE;
		}
		g_mutex_unlock(&gui.cairo_mutex);
		/* disp_surface invalidated by viewport change; buf is still valid -
		 * fall through to rebuild disp_surface from buf */
	}

	/* catch and compute rendering data */
	dd.cr = cr;
	dd.window_width = gtk_widget_get_allocated_width(widget);
	dd.window_height = gtk_widget_get_allocated_height(widget);

	if (dd.window_width != gui.view[dd.vport].view_width ||
			dd.window_height != gui.view[dd.vport].view_height) {
		//siril_log_debug("draw area and disp surface size mismatch: %d,%d vs %d,%d\n",
		//		dd.window_width, dd.window_height,
		//		gui.view[dd.vport].view_width, gui.view[dd.vport].view_height);
		invalidate_image_render_cache(dd.vport);
	}

	dd.zoom = get_zoom_val();
	dd.image_width = gfit->rx;
	dd.image_height = gfit->ry;
	dd.filter = (dd.zoom < 1.0) ? CAIRO_FILTER_GOOD : CAIRO_FILTER_FAST;

	GVariant *state = g_action_get_state(action_neg);
	dd.neg_view = g_variant_get_boolean(state);
	g_variant_unref(state);

#if 0
	static struct timeval prevtime = { 0 };
	struct timeval now;
	gettimeofday(&now, NULL);
	int us = (now.tv_sec - prevtime.tv_sec) * 1000000 + now.tv_usec - prevtime.tv_usec;
	prevtime = now;

	GdkRectangle rect;
	if (!gdk_cairo_get_clip_rectangle(cr, &rect) || (rect.width == 0 && rect.height == 0)) {
		printf("nothing to redraw\n");
		return FALSE;
	}
	if (us < 320000)
		printf("redraw %d ms\t(at %d, %d of size %d x %d)\n", us/1000,
				rect.x, rect.y, rect.width, rect.height);
#endif


	adjust_vport_size_to_image();

	/* RGB or gray images */
	draw_main_image(&dd);

	/* selection rectangle */
	draw_selection(&dd);

	/* ROI */
	draw_roi(&dd);

	/* cut line */
	draw_cut_line(&dd);

	/* draw measurement line */
	draw_measurement_line(&dd);

	/* draw user polygons */
	draw_in_progress_poly(&dd);
	draw_user_polygons(&dd);

	/* detected stars and highlight the selected star */
	g_mutex_lock(&com.mutex);
	draw_stars(&dd);
	g_mutex_unlock(&com.mutex);

	/* celestial grid */
	draw_wcs_grid(&dd);

	/* distortions */
	draw_wcs_disto(&dd);

	/* detected objects */
	draw_annotates(&dd);

	/* analysis tool */
	draw_analysis(&dd);

	/* background removal gradient selection boxes */
	draw_brg_boxes(&dd);

	/* multipoint planetary alignment-point grid (active after Analyze) */
	draw_mpp_aps(&dd);

	/* registration framing*/
	draw_regframe(&dd);

	/* RGB composition center points */
	draw_rgb_centers(&dd);

	/* allow custom rendering */
	if (gui.draw_extra)
		gui.draw_extra(&dd);

	return FALSE;
}

point get_center_of_vport() {
	image_display_init_statics();
	GtkWidget *widget = imgdisp_drawing_r;

	guint window_width = gtk_widget_get_allocated_width(widget);
	guint window_height = gtk_widget_get_allocated_height(widget);

	point center = { window_width / 2., window_height / 2. };

	return center;
}

void add_image_and_label_to_cairo(cairo_t *cr, int vport) {
	draw_data_t dd;
	image_display_init_statics();
	GtkWidget *widget = imgdisp_drawing_r;
	GAction *action_neg = g_action_map_lookup_action(G_ACTION_MAP(imgdisp_app_win), "negative-view");
	GVariant *state = g_action_get_state(action_neg);

	dd.vport = vport;
	dd.cr = cr;
	dd.window_width = gtk_widget_get_allocated_width(widget);
	dd.window_height = gtk_widget_get_allocated_height(widget);
	dd.zoom = get_zoom_val();
	dd.image_width = gfit->rx;
	dd.image_height = gfit->ry;
	dd.filter = (dd.zoom < 1.0) ? CAIRO_FILTER_GOOD : CAIRO_FILTER_FAST;
	dd.neg_view = g_variant_get_boolean(state);
	g_variant_unref(state);

	/* RGB or gray images */
	draw_main_image(&dd);
	/* detected stars and highlight the selected star */
	g_mutex_lock(&com.mutex);
	draw_stars(&dd);
	g_mutex_unlock(&com.mutex);
	/* wcs_grid */
	draw_wcs_grid(&dd);
	/* detected objects */
	draw_annotates(&dd);
	/* analysis */
	draw_analysis(&dd);
	/* distortions */
	draw_wcs_disto(&dd);
}
