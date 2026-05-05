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
 * Headless-build stubs for GUI functions still called directly from core/
 * processing files.  This file is compiled ONLY when -Dgui=disabled; in the
 * normal GUI build the real implementations in src/gui/ and
 * src/livestacking/ are linked instead.
 *
 * Every stub is a no-op that returns a safe default value (0 / NULL / FALSE).
 * These stubs exist so the headless CLI (siril-cli) can link cleanly; actual
 * GUI functionality is simply absent in headless mode.
 *
 * When a function here is later properly routed through gui_iface, its entry
 * can be removed from this file.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <gsl/gsl_histogram.h>

#include "core/siril.h"
#include "io/siril_plot.h"

/* Forward-declare opaque GTK/GUI types used in stubs below. */
typedef struct _GtkWidget   GtkWidget;
typedef struct mtf_data     mtf_data;
typedef struct ght_data     ght_data;

/* int-typed enum forward compat (callers use int-compatible enum) */
typedef int main_tabs;

/* ── Global data variables ──────────────────────────────────────────────── */

/* gui/photometric_cc.c's global; used by command.c */
gchar **spcc_mirrors = NULL;

/* ── execute_idle_and_wait_for_it ───────────────────────────────────────── */
/* Headless: just call the function directly on the current thread. */
void execute_idle_and_wait_for_it(gboolean (*idle_function)(gpointer), gpointer arg) {
	if (idle_function)
		idle_function(arg);
}

/* ── Progress bar ───────────────────────────────────────────────────────── */
void set_progress_bar_data(const char *text, double percent) { (void)text; (void)percent; }

/* ── Display / histogram ────────────────────────────────────────────────── */
void display_filename(void) {}
void update_display_fwhm(void) {}
void update_zoom_label(void) {}
void initialize_display_mode(void) {}
void set_display_mode(void) {}
void remap_all(void) {}
void invalidate_gfit_histogram(void) {}
void refresh_histogram_if_visible(void) {}
void update_gfit_histogram_if_needed(void) {}
gsl_histogram *computeHisto(fits *fit, int layer) { (void)fit; (void)layer; return NULL; }
void compute_histo_for_fit(fits *thefit) { (void)thefit; }
gsl_histogram *computeHisto_Selection(fits *fit, int layer, rectangle *selection) {
	(void)fit; (void)layer; (void)selection; return NULL;
}

/* ── Stars / PSF ────────────────────────────────────────────────────────── */
psf_star *add_star(fits *fit, int layer, int *index) {
	(void)fit; (void)layer; (void)index; return NULL;
}
void clear_stars_list(gboolean refresh_GUI) { (void)refresh_GUI; }
gboolean clear_stars_list_as_idle(gpointer user_data) { (void)user_data; return FALSE; }
gboolean DrawPSF(gpointer user_data) { (void)user_data; return FALSE; }
gboolean end_findstar(gpointer p) { (void)p; return FALSE; }

/* ── Sequence / registration UI ─────────────────────────────────────────── */
int update_sequences_list(const char *sequence_name_to_select) {
	(void)sequence_name_to_select; return 0;
}
void fill_sequence_list(sequence *seq, int layer, gboolean as_idle) {
	(void)seq; (void)layer; (void)as_idle;
}
void sequence_list_change_current(void) {}
void adjust_reginfo(void) {}
void adjust_refimage(int n) { (void)n; }
void adjust_sellabel(void) {}
void update_reg_interface(gboolean dont_change_reg_radio) { (void)dont_change_reg_radio; }
void update_stack_interface(gboolean dont_change_stack_type) { (void)dont_change_stack_type; }
void update_prepro_interface(gboolean allow_debayer) { (void)allow_debayer; }
void update_seqlist(int layer) { (void)layer; }
void update_seq_gui_idle_thread_func(gpointer data) { (void)data; }
int set_layers_for_registration(void) { return 0; }
void enable_view_reference_checkbox(gboolean status) { (void)status; }
gboolean end_register_idle(gpointer p) { (void)p; return FALSE; }
void ensure_seqlist_dialog_closed(void) {}
gboolean update_MenuItem(gpointer user_data) { (void)user_data; return FALSE; }

/* ── Image loading / display ────────────────────────────────────────────── */
gboolean end_image_loading(gpointer arg) { (void)arg; return FALSE; }
gboolean init_right_tab(gpointer user_data) { (void)user_data; return FALSE; }
void drawPlot(void) {}
void clear_previews(void) {}
gboolean redraw_previews(gpointer user_data) { (void)user_data; return FALSE; }
void queue_redraw_mask(void) {}
void show_or_hide_mask_tab(void) {}
gboolean show_or_hide_mask_tab_idle(gpointer p) { (void)p; return FALSE; }

/* ── Window / tab / dialog control ─────────────────────────────────────── */
void control_window_switch_to_tab(main_tabs tab) { (void)tab; }
gboolean close_tab(gpointer user_data) { (void)user_data; return FALSE; }
int number_of_dialogs(void) { return 0; }

/* ── Sliders / precision / debayer controls ─────────────────────────────── */
void set_cutoff_sliders_values(void) {}
void set_cutoff_sliders_max_values(void) {}
void sliders_mode_set_state(sliders_mode s) { (void)s; }
gboolean set_precision_switch(gpointer user_data) { (void)user_data; return FALSE; }
void enable_debayer(gboolean arg) { (void)arg; }
void update_debayer_button_status(gboolean new_state) { (void)new_state; }
gboolean set_kernel_size_in_gui(gpointer user_data) { (void)user_data; return FALSE; }
void set_GUI_CAMERA(void) {}

/* ── ROI ────────────────────────────────────────────────────────────────── */
void on_clear_roi(void) {}
int populate_roi(void) { return 0; }
void copy_roi_into_gfit(void) {}
void delete_selected_area(void) {}
void lock_roi_mutex(void) {}
void unlock_roi_mutex(void) {}

/* ── Preview / backup ───────────────────────────────────────────────────── */
gboolean is_preview_active(void) { return FALSE; }
void copy_gfit_to_backup(void) {}
int copy_backup_to_gfit(void) { return 0; }
fits *get_preview_gfit_backup(void) { return NULL; }

/* ── Script widgets ─────────────────────────────────────────────────────── */
gboolean script_widgets_idle(gpointer user_data) { (void)user_data; return FALSE; }
gpointer refresh_scripts_in_thread(gpointer user_data) { (void)user_data; return NULL; }

/* ── Photometry / plot ──────────────────────────────────────────────────── */
void notify_new_photometry(void) {}
void clear_all_photometry_and_plot(void) {}
gboolean create_new_siril_plot_window(gpointer p) { (void)p; return FALSE; }
gboolean save_siril_plot_to_clipboard(siril_plot_data *spl_data, int width, int height) {
	(void)spl_data; (void)width; (void)height; return FALSE;
}
void init_plot_colors(void) {}
void drawPlot_idle(void) {}
gchar *get_log_as_string(void) { return NULL; }

/* ── Astrometry / catsearch end idles ──────────────────────────────────── */
gboolean end_plate_solver(gpointer p) { (void)p; return FALSE; }
gboolean end_platesolve_sequence(gpointer p) { (void)p; return FALSE; }
gboolean end_process_catsearch(gpointer p) { (void)p; return FALSE; }

/* ── Keywords dialog ────────────────────────────────────────────────────── */
gboolean refresh_keywords_dialog(gpointer user_data) { (void)user_data; return FALSE; }

/* ── ICC profile ────────────────────────────────────────────────────────── */
void check_gfit_profile_identical_to_monitor(void) {}

/* ── Aberration inspector ───────────────────────────────────────────────── */
void compute_aberration_inspector(void) {}

/* ── MTF / GHT hooks ────────────────────────────────────────────────────── */
struct mtf_data *create_mtf_data(void) { return NULL; }
void destroy_mtf_data(void *args) { (void)args; }
gchar *mtf_log_hook(gpointer p, log_hook_detail detail) { (void)p; (void)detail; return NULL; }
int mtf_single_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	(void)args; (void)fit; (void)threads; return 0;
}
gchar *invmtf_log_hook(gpointer p, log_hook_detail detail) { (void)p; (void)detail; return NULL; }
int invmtf_single_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	(void)args; (void)fit; (void)threads; return 0;
}
struct ght_data *create_ght_data(void) { return NULL; }
void destroy_ght_data(void *args) { (void)args; }
gchar *ght_log_hook(gpointer p, log_hook_detail detail) { (void)p; (void)detail; return NULL; }
int ght_single_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	(void)args; (void)fit; (void)threads; return 0;
}

/* ── Cut profile ────────────────────────────────────────────────────────── */
void initialize_cut_struct(cut_struct *arg) { if (arg) memset(arg, 0, sizeof(*arg)); }
gboolean cut_struct_is_valid(cut_struct *arg) { (void)arg; return FALSE; }
void free_cut_args(cut_struct *arg) { (void)arg; }
void apply_cut_to_sequence(cut_struct *cut_args) { (void)cut_args; }
gpointer cut_profile(gpointer p) { (void)p; return NULL; }
gpointer tri_cut(gpointer p) { (void)p; return NULL; }
gpointer cfa_cut(gpointer p) { (void)p; return NULL; }
void reset_cut_gui_filedependent(gpointer user_data) { (void)user_data; }

/* ── MTF/GHT sequence wrappers ──────────────────────────────────────────── */
void apply_mtf_to_sequence(struct mtf_data *mtf_args) { (void)mtf_args; }
void apply_ght_to_sequence(struct ght_data *ght_args) { (void)ght_args; }

/* ── Compositing / remixer ──────────────────────────────────────────────── */
void clear_stars_list_as_idle_h(gpointer user_data) { (void)user_data; }
int toggle_remixer_window_visibility(int _invocation, fits *_fit_left, fits *_fit_right) {
	(void)_invocation; (void)_fit_left; (void)_fit_right; return 0;
}

/* ── SPCC / photometric ─────────────────────────────────────────────────── */
void initialize_spcc_mirrors(void) {}
int get_favourite_spccobject(GList *list, const gchar *favourite) {
	(void)list; (void)favourite; return 0;
}
int get_favourite_oscsensor(GList *list, const gchar *favourite) {
	(void)list; (void)favourite; return 0;
}

/* ── Denoise ────────────────────────────────────────────────────────────── */
gboolean denoise_apply_idle(gpointer p) { (void)p; return FALSE; }

/* ── Livestacking ───────────────────────────────────────────────────────── */
void livestacking_display(gchar *str, gboolean free_after_display) {
	(void)str; (void)free_after_display;
}
gboolean livestacking_first_result_idle(gpointer p) { (void)p; return FALSE; }
void livestacking_update_number_of_images(int nb, double total_exposure, double noise, const char *process_time) {
	(void)nb; (void)total_exposure; (void)noise; (void)process_time;
}

/* ── Polygon helpers ────────────────────────────────────────────────────── */
UserPolygon *find_polygon_by_id(int id) { (void)id; return NULL; }
void free_user_polygon(gpointer data) { (void)data; }
int get_unused_polygon_id(void) { return 0; }
void clear_user_polygons(void) {}
gboolean delete_user_polygon(int id) { (void)id; return FALSE; }
UserPolygon *deserialize_polygon(const uint8_t *data, size_t size) {
	(void)data; (void)size; return NULL;
}
uint8_t *serialize_polygon(UserPolygon *polygon, size_t *size) {
	(void)polygon; if (size) *size = 0; return NULL;
}
uint8_t *serialize_polygon_list(GSList *polygons, size_t *out_size) {
	(void)polygons; if (out_size) *out_size = 0; return NULL;
}

/* ── Misc UI helpers ────────────────────────────────────────────────────── */
void set_suggested(GtkWidget *widget) { (void)widget; }
void unset_suggested(GtkWidget *widget) { (void)widget; }
int match_drawing_area_widget(const GtkWidget *drawing_area, gboolean allow_rgb) {
	(void)drawing_area; (void)allow_rgb; return 0;
}
int select_vport(int vport) { return vport < 3 ? vport : 0; }
void update_display_fwhm_idle(void) {}
void set_suggested_ht(GtkWidget *w) { (void)w; }
OverrangeResponse apply_limits(fits *fit, double minval, double maxval, OverrangeResponse method) {
	(void)fit; (void)minval; (void)maxval; return method;
}
gchar *build_save_filename(gchar *prepend, gchar *ext, gboolean forsequence, gboolean add_time_stamp) {
	(void)prepend; (void)ext; (void)forsequence; (void)add_time_stamp; return NULL;
}
