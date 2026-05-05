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
 * Functions fall into two categories:
 *   A) Pure GUI — no-op stubs are correct: the function only updates widgets.
 *   B) Mixed — the real function does processing (memory, state, threads) AND
 *      GUI updates.  These stubs perform the non-GUI work and skip the GTK
 *      calls, so headless operation stays correct.
 *
 * When a function is properly routed through gui_iface, its entry here can be
 * removed.  headless_stubs.c acts as the inventory of remaining direct calls.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <gsl/gsl_histogram.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "algos/astrometry_solver.h"
#include "algos/star_finder.h"
#include "registration/registration.h"
#include "filters/mtf.h"
#include "filters/ght.h"
#include "io/sequence.h"
#include "io/siril_plot.h"

/* ── Forward declarations for types not fully defined without GTK ───────── */
typedef struct _GtkWidget   GtkWidget;

/* ── Forward declarations for functions defined later in this file ──────── */
void destroy_mtf_data(void *args);
void destroy_ght_data(void *args);

/* ── Global data variables ──────────────────────────────────────────────── */

/* gui/photometric_cc.c's global — needed by command.c SPCC support */
gchar **spcc_mirrors = NULL;

/* ── execute_idle_and_wait_for_it ───────────────────────────────────────── */
/*
 * Category B: in GUI mode this schedules func on the main GTK thread and
 * blocks.  Headless has no event loop, so we call func directly.
 */
void execute_idle_and_wait_for_it(gboolean (*idle_function)(gpointer), gpointer arg) {
	if (idle_function)
		idle_function(arg);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Category B — stubs that perform non-GUI cleanup
 * ══════════════════════════════════════════════════════════════════════════ */

/* ── Stars list ─────────────────────────────────────────────────────────── */

/*
 * Real clear_stars_list already guards all GTK calls with
 * "if (refresh_GUI && !com.headless)", so we reproduce only the
 * lock + free + state-reset path here.
 */
void clear_stars_list(gboolean refresh_GUI) {
	(void)refresh_GUI;
	g_rw_lock_writer_lock(&com.stars_lock);
	psf_star **stars = com.stars;
	if (stars) {
		com.stars = NULL;
		g_rw_lock_writer_unlock(&com.stars_lock);
		if (stars[0]) {
			if (stars[1] || !com.star_is_seqdata) {
				int i = 0;
				while (i < MAX_STARS && stars[i])
					free_psf(stars[i++]);
			}
			free(stars);
		}
	} else {
		g_rw_lock_writer_unlock(&com.stars_lock);
	}
	com.star_is_seqdata = FALSE;
}

gboolean clear_stars_list_as_idle(gpointer user_data) {
	(void)user_data;
	clear_stars_list(FALSE);
	return FALSE;
}

/* ── Photometry cleanup ─────────────────────────────────────────────────── */

void clear_all_photometry_and_plot(void) {
	for (int i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++)
		free_photometry_set(&com.seq, i);
	clear_stars_list(FALSE);
	/* reset_plot() and drawPlot() are GUI-only — skip */
}

/* ── Processing-thread end idles ────────────────────────────────────────── */

/* end_findstar: stop thread + free args */
gboolean end_findstar(gpointer p) {
	struct starfinder_data *args = (struct starfinder_data *)p;
	stop_processing_thread();
	if (args) {
		g_free(args->starfile);
		free(args);
	}
	return FALSE;
}

/* end_register_idle: stop thread + update seq state + free args */
gboolean end_register_idle(gpointer p) {
	struct registration_args *args = (struct registration_args *)p;
	stop_processing_thread();
	if (args) {
		args->seq->reg_invalidated = FALSE;
		free(args->new_seq_name);
		free(args);
	}
	return FALSE;
}

/* end_plate_solver: stop thread + free astrometry data */
gboolean end_plate_solver(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *)p;
	stop_processing_thread();
	free_astrometry_data(args);
	return FALSE;
}

/* end_platesolve_sequence: sequence cleanup + free */
gboolean end_platesolve_sequence(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *)p;
	if (args) {
		if (!check_seq_is_comseq(args->seq))
			free_sequence(args->seq, TRUE);
		else if (!args->retval && !args->has_output) {
			gchar *seqname;
			if (g_str_has_suffix(com.seq.seqname, ".seq"))
				seqname = g_strdup(com.seq.seqname);
			else
				seqname = g_strdup_printf("%s.seq", com.seq.seqname);
			set_seq(seqname);
			g_free(seqname);
		}
		free(p);
	}
	return end_generic(NULL);
}

/* denoise_apply_idle: stop thread + free args */
gboolean denoise_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	free_generic_img_args(args);
	/* close_dialog calls omitted — no GTK in headless */
	return FALSE;
}

/* ── Cut struct lifecycle ───────────────────────────────────────────────── */

void initialize_cut_struct(cut_struct *arg) {
	if (!arg) return;
	arg->fit = gfit;
	arg->seq = NULL;
	arg->imgnumber = -1;
	arg->cut_start.x = -1;  arg->cut_start.y = -1;
	arg->cut_end.x   = -1;  arg->cut_end.y   = -1;
	arg->cut_wn1.x   = -1;  arg->cut_wn1.y   = -1;
	arg->cut_wn2.x   = -1;  arg->cut_wn2.y   = -1;
	arg->cut_measure     = FALSE;
	arg->plot_as_wavenumber = FALSE;
	arg->wavenumber1     = -1;
	arg->wavenumber2     = -1;
	arg->plot_spectro_bg = FALSE;
	arg->bg_poly_order   = 3;
	arg->tri             = FALSE;
	arg->mode            = CUT_MONO;
	arg->width           = 1;
	arg->step            = 1;
	arg->display_graph   = TRUE;
	g_free(arg->filename);   arg->filename    = NULL;
	g_free(arg->title);      arg->title       = NULL;
	g_free(arg->user_title); arg->user_title  = NULL;
	arg->title_has_sequence_numbers = FALSE;
	arg->save_dat        = FALSE;
	arg->save_png_too    = FALSE;
	arg->pref_as         = gfit->keywords.wcsdata.pltsolvd;
	arg->vport           = -1;
	/* gui_function(reset_cut_gui, NULL) is a no-op in headless (com.headless) */
}

void free_cut_args(cut_struct *arg) {
	if (!arg) return;
	g_free(arg->filename);    arg->filename   = NULL;
	g_free(arg->title);       arg->title      = NULL;
	g_free(arg->user_title);  arg->user_title = NULL;
	/* In headless mode there is no gui.cut, so always free */
	free(arg);
}

/* ── MTF / GHT data structures (structs moved to filters/mtf.h, ght.h) ─── */

struct mtf_data *create_mtf_data(void) {
	struct mtf_data *data = calloc(1, sizeof(struct mtf_data));
	if (!data) return NULL;
	data->linked       = TRUE;
	data->destroy_fn   = destroy_mtf_data;
	return data;
}

void destroy_mtf_data(void *args) {
	if (!args) return;
	struct mtf_data *data = (struct mtf_data *)args;
	free(data->seqEntry);
	free(data);
}

struct ght_data *create_ght_data(void) {
	struct ght_data *data = calloc(1, sizeof(struct ght_data));
	if (!data) return NULL;
	data->destroy_fn = destroy_ght_data;
	return data;
}

void destroy_ght_data(void *args) {
	if (!args) return;
	struct ght_data *data = (struct ght_data *)args;
	free(data->seqEntry);
	free(data->params_ght);
	free(data);
}

/* ── MTF / GHT processing hooks (Category B: do real work, skip GTK) ────── */

int mtf_single_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	(void)threads;
	struct mtf_data *data = (struct mtf_data *)args->user;
	if (!data) return 1;
	if (data->linked)
		apply_linked_mtf_to_fits(fit, fit, data->params, TRUE);
	else
		apply_unlinked_mtf_to_fits(fit, fit, data->uparams);
	/* auto_display_compensation path accesses gui.icc.monitor — skip headless */
	return 0;
}

int invmtf_single_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	(void)threads;
	struct mtf_data *data = (struct mtf_data *)args->user;
	if (!data) return 1;
	apply_linked_pseudoinverse_mtf_to_fits(fit, fit, data->params, TRUE);
	return 0;
}

int ght_single_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	(void)threads;
	struct ght_data *data = (struct ght_data *)args->user;
	if (!data || !data->params_ght) return 1;
	if (data->params_ght->payne_colourstretchmodel == COL_SAT) {
		/* SAT path needs prepare-hook buffers — call the dedicated function */
		apply_sat_ght_to_fits(fit, data->params_ght, TRUE);
	} else {
		apply_linked_ght_to_fits(fit, fit, data->params_ght, TRUE);
	}
	/* auto_display_compensation accesses gui.icc.monitor — skip headless */
	return 0;
}

/* log hooks — return NULL in headless (callers handle NULL gracefully) */
gchar *mtf_log_hook(gpointer p, log_hook_detail detail) { (void)p; (void)detail; return NULL; }
gchar *invmtf_log_hook(gpointer p, log_hook_detail detail) { (void)p; (void)detail; return NULL; }
gchar *ght_log_hook(gpointer p, log_hook_detail detail) { (void)p; (void)detail; return NULL; }

/* ── SPCC mirror initialisation ─────────────────────────────────────────── */

void initialize_spcc_mirrors(void) {
	g_strfreev(spcc_mirrors);
	spcc_mirrors = g_new(gchar *, 2);
	spcc_mirrors[0] = g_strdup("https://zenodo.org/records/17988559/files");
	spcc_mirrors[1] = NULL;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Category A — pure-GUI no-op stubs
 * ══════════════════════════════════════════════════════════════════════════ */

/* Progress bar */
void set_progress_bar_data(const char *text, double percent) { (void)text; (void)percent; }

/* Display / histogram */
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

/* Stars / PSF */
psf_star *add_star(fits *fit, int layer, int *index) {
	(void)fit; (void)layer; (void)index; return NULL;
}
gboolean DrawPSF(gpointer user_data) { (void)user_data; return FALSE; }

/* Sequence / registration UI */
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
gboolean end_register_idle_noop(gpointer p) { (void)p; return FALSE; }
void ensure_seqlist_dialog_closed(void) {}
gboolean update_MenuItem(gpointer user_data) { (void)user_data; return FALSE; }

/* Image loading / display */
gboolean end_image_loading(gpointer arg) { (void)arg; return FALSE; }
gboolean init_right_tab(gpointer user_data) { (void)user_data; return FALSE; }
void drawPlot(void) {}
void clear_previews(void) {}
gboolean redraw_previews(gpointer user_data) { (void)user_data; return FALSE; }
void queue_redraw_mask(void) {}
void show_or_hide_mask_tab(void) {}
gboolean show_or_hide_mask_tab_idle(gpointer p) { (void)p; return FALSE; }

/* Window / tab / dialog control */
void control_window_switch_to_tab(int tab) { (void)tab; }
gboolean close_tab(gpointer user_data) { (void)user_data; return FALSE; }
int number_of_dialogs(void) { return 0; }

/* Sliders / precision / debayer controls */
void set_cutoff_sliders_values(void) {}
void set_cutoff_sliders_max_values(void) {}
void sliders_mode_set_state(sliders_mode s) { (void)s; }
gboolean set_precision_switch(gpointer user_data) { (void)user_data; return FALSE; }
void enable_debayer(gboolean arg) { (void)arg; }
void update_debayer_button_status(gboolean new_state) { (void)new_state; }
gboolean set_kernel_size_in_gui(gpointer user_data) { (void)user_data; return FALSE; }
void set_GUI_CAMERA(void) {}

/* ROI */
void on_clear_roi(void) {}
int populate_roi(void) { return 0; }
void copy_roi_into_gfit(void) {}
void delete_selected_area(void) {}
void lock_roi_mutex(void) {}
void unlock_roi_mutex(void) {}

/* Preview / backup — GUI-only concept; no-ops in headless */
gboolean is_preview_active(void) { return FALSE; }
void copy_gfit_to_backup(void) {}
int copy_backup_to_gfit(void) { return 0; }
fits *get_preview_gfit_backup(void) { return NULL; }

/* Script widgets */
gboolean script_widgets_idle(gpointer user_data) { (void)user_data; return FALSE; }
gpointer refresh_scripts_in_thread(gpointer user_data) { (void)user_data; return NULL; }

/* Photometry / plot */
void notify_new_photometry(void) {}
gboolean create_new_siril_plot_window(gpointer p) { (void)p; return FALSE; }
gboolean save_siril_plot_to_clipboard(siril_plot_data *spl_data, int width, int height) {
	(void)spl_data; (void)width; (void)height; return FALSE;
}
void init_plot_colors(void) {}
gchar *get_log_as_string(void) { return NULL; }

/* Astrometry end idles */
gboolean end_process_catsearch(gpointer p) { (void)p; return FALSE; }

/* Keywords dialog */
gboolean refresh_keywords_dialog(gpointer user_data) { (void)user_data; return FALSE; }

/* ICC profile */
void check_gfit_profile_identical_to_monitor(void) {}

/* Aberration inspector */
void compute_aberration_inspector(void) {}

/* sequence-level MTF/GHT — complex GUI+sequence setup; no-op in headless */
void apply_mtf_to_sequence(struct mtf_data *mtf_args) { (void)mtf_args; }
void apply_ght_to_sequence(struct ght_data *ght_args) { (void)ght_args; }

/* Cut profile sequence/tri/CFA wrappers */
gboolean cut_struct_is_valid(cut_struct *arg) { (void)arg; return FALSE; }
void apply_cut_to_sequence(cut_struct *cut_args) { (void)cut_args; }
gpointer cut_profile(gpointer p) { (void)p; return NULL; }
gpointer tri_cut(gpointer p) { (void)p; return NULL; }
gpointer cfa_cut(gpointer p) { (void)p; return NULL; }
void reset_cut_gui_filedependent(gpointer user_data) { (void)user_data; }

/* Compositing / remixer */
int toggle_remixer_window_visibility(int _invocation, fits *_fit_left, fits *_fit_right) {
	(void)_invocation; (void)_fit_left; (void)_fit_right; return 0;
}

/* SPCC / photometric */
int get_favourite_spccobject(GList *list, const gchar *favourite) {
	(void)list; (void)favourite; return 0;
}
int get_favourite_oscsensor(GList *list, const gchar *favourite) {
	(void)list; (void)favourite; return 0;
}

/* Livestacking display */
void livestacking_display(gchar *str, gboolean free_after_display) {
	(void)str; (void)free_after_display;
}
gboolean livestacking_first_result_idle(gpointer p) { (void)p; return FALSE; }
void livestacking_update_number_of_images(int nb, double total_exposure, double noise, const char *process_time) {
	(void)nb; (void)total_exposure; (void)noise; (void)process_time;
}

/* Polygon helpers */
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

/* Misc UI helpers */
void set_suggested(GtkWidget *widget) { (void)widget; }
void unset_suggested(GtkWidget *widget) { (void)widget; }
int match_drawing_area_widget(const GtkWidget *drawing_area, gboolean allow_rgb) {
	(void)drawing_area; (void)allow_rgb; return 0;
}
int select_vport(int vport) { return vport < 3 ? vport : 0; }
OverrangeResponse apply_limits(fits *fit, double minval, double maxval, OverrangeResponse method) {
	(void)fit; (void)minval; (void)maxval; return method;
}
gchar *build_save_filename(gchar *prepend, gchar *ext, gboolean forsequence, gboolean add_time_stamp) {
	(void)prepend; (void)ext; (void)forsequence; (void)add_time_stamp; return NULL;
}
