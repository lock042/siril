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

/*
 * GTK-free declarations of GUI functions called directly from core/io code.
 *
 * These functions are implemented in either:
 *   src/gui/       — real GTK implementations (GUI build, src_files_gui)
 *   headless_stubs.c — no-op / minimal stubs  (headless build, siril_headless_lib)
 *
 * Non-GUI source files (src_files) that need to call these functions should
 * include this header instead of the GTK-heavy gui/ headers where these
 * declarations used to live.  The header itself has no GTK dependency.
 *
 * Functions here should be progressively replaced by gui_iface slots as the
 * separation work continues.  When a function migrates to gui_iface, remove
 * it from this header and from headless_stubs.c.
 */

#ifndef SRC_CORE_GUI_CALLS_H_
#define SRC_CORE_GUI_CALLS_H_

#include <glib.h>
#include "core/siril.h"   /* fits, sequence, sliders_mode, main_tabs, … */

/* ── Slider / display-range ─────────────────────────────────────────────── */
void sliders_mode_set_state(sliders_mode sliders);
void set_cutoff_sliders_max_values(void);
void set_cutoff_sliders_values(void);
void set_display_mode(void);

/* ── Full-image redraw ───────────────────────────────────────────────────── */
void remap_all(void);

/* ── Histogram ───────────────────────────────────────────────────────────── */
void invalidate_gfit_histogram(void);
void update_gfit_histogram_if_needed(void);
void compute_histo_for_fit(fits *fit);
void refresh_histogram_if_visible(void);

/* ── Selection / ROI ─────────────────────────────────────────────────────── */
gpointer on_clear_roi(void);
void delete_selected_area(void);

/* ── Sequence list ───────────────────────────────────────────────────────── */
void fill_sequence_list(sequence *seq, int layer, gboolean as_idle);
void sequence_list_change_current(void);

/* ── Registration reference view ─────────────────────────────────────────── */
void enable_view_reference_checkbox(gboolean status);

/* ── Tab / panel ─────────────────────────────────────────────────────────── */
/* Prefer gui_iface.switch_to_tab() for new code. */
void control_window_switch_to_tab(main_tabs tab);

/* ── Misc. GUI state ─────────────────────────────────────────────────────── */
/* Schedules func(arg) on the GTK main thread and blocks; in headless mode
 * calls func(arg) directly.  Prefer gui_iface.execute_idle_sync() for new
 * code. */
void execute_idle_and_wait_for_it(gboolean (*idle_function)(gpointer), gpointer arg);
gboolean close_tab(gpointer user_data);
gboolean init_right_tab(gpointer user_data);
void initialize_display_mode(void);
void display_filename(void);
void update_display_fwhm(void);
void update_zoom_label(void);
void update_prepro_interface(gboolean allow_debayer);
void adjust_sellabel(void);
void adjust_reginfo(void);
void adjust_refimage(int n);
int  set_layers_for_registration(void);
void queue_redraw_mask(void);
gboolean set_precision_switch(gpointer user_data);
void set_GUI_CAMERA(void);
gboolean update_MenuItem(gpointer user_data);
void update_seqlist(int layer);
gpointer update_seq_gui_idle_thread_func(gpointer data);
void ensure_seqlist_dialog_closed(void);

/* ── ROI / mask ─────────────────────────────────────────────────────────── */
int  populate_roi(void);
void copy_roi_into_gfit(void);
void lock_roi_mutex(void);
void unlock_roi_mutex(void);
void show_or_hide_mask_tab(void);
gboolean show_or_hide_mask_tab_idle(gpointer p);

/* ── Dialog management ──────────────────────────────────────────────────── */
int number_of_dialogs(void);

/* ── Preview / display refresh ─────────────────────────────────────────── */
void clear_previews(void);
gboolean redraw_previews(gpointer user_data);
void check_gfit_profile_identical_to_monitor(void);

/* reset_conv_kernel: resets the kernel display widget — GUI only */
void reset_conv_kernel(void);

/* ── StarNet / Remixer integration ──────────────────────────────────────── */
#define CALL_FROM_STARNET 1
int toggle_remixer_window_visibility(int invocation, fits *fit_left, fits *fit_right);

/* ── HEIF multi-image dialog ─────────────────────────────────────────────
 * Declared in gui/utils.h; stub in headless_stubs.c. */
#ifdef HAVE_LIBHEIF
struct heif_context;
gboolean heif_dialog(struct heif_context *heif, uint32_t *selected_image);
#endif

#endif /* SRC_CORE_GUI_CALLS_H_ */
