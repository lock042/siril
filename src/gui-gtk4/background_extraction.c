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

#include "core/siril.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "algos/background_extraction.h"
#include "algos/statistics.h"
#include "algos/demosaicing.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/siril_actions.h"

static GtkDropDown *bkg_poly_order_combo = NULL;
static GtkDropDown *bkg_correction_combo = NULL;
static GtkSpinButton *bkg_nb_samples_spin = NULL;
static GtkRange *bkg_tolerance_scale = NULL;
static GtkDropDown *bkg_interp_combo = NULL;
static GtkSpinButton *bkg_smoothing_spin = NULL;
static GtkCheckButton *bkg_dither_btn = NULL;
static GtkCheckButton *bkg_randomize_btn = NULL;
static GtkCheckButton *bkg_grad_descent_btn = NULL;
static GtkWidget *bkg_ok_button = NULL;
static GtkWidget *bkg_compute_bkg_button = NULL;
static GtkDropDown *bkg_view_combo = NULL;
static GtkWidget *bkg_view_box = NULL;
static GtkLabel *bkg_label_samples = NULL;
static GtkCheckButton *bkg_keep_samples_btn = NULL;
static GtkCheckButton *bkg_seq_btn = NULL;
static GtkEntry *bkg_seq_entry = NULL;
static GtkNotebook *bkg_notebook = NULL;

/* method selection + automatic (sample-free) model widgets */
static GtkDropDown *bkg_method_combo = NULL;
static GtkWidget *bkg_samples_expander = NULL;
static GtkWidget *bkg_auto_expander = NULL;
static GtkSpinButton *ag_scale_spin = NULL;
static GtkSpinButton *ag_smoothness_spin = NULL;
static GtkCheckButton *ag_protect_check = NULL;
static GtkSpinButton *ag_protect_threshold_spin = NULL;
static GtkSpinButton *ag_protect_amount_spin = NULL;
static GtkWidget *ag_protect_threshold_label = NULL;
static GtkWidget *ag_protect_amount_label = NULL;
static GtkCheckButton *ag_simplified_check = NULL;
static GtkSpinButton *ag_degree_spin = NULL;
static GtkWidget *ag_degree_label = NULL;
static GtkDropDown *ag_downsample_combo = NULL;

static void bkg_sync_method(void);

static void background_extraction_init_statics(void) {
	if (bkg_poly_order_combo) return;
	bkg_poly_order_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "box_background_order"));
	bkg_correction_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "box_background_correction"));
	bkg_nb_samples_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_background_nb_samples"));
	bkg_tolerance_scale = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_background_tolerance"));
	bkg_interp_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "background_extraction_combo"));
	bkg_smoothing_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_background_smoothing"));
	bkg_dither_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "bkg_dither_button"));
	bkg_randomize_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "bkg_randomize_button"));
	bkg_grad_descent_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "bkg_grad_descent_button"));
	bkg_ok_button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "background_ok_button"));
	bkg_compute_bkg_button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "bkg_compute_bkg"));
	bkg_view_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "bkg_view_combo"));
	bkg_view_box = GTK_WIDGET(gtk_builder_get_object(gui.builder, "bkg_view_box"));
	bkg_label_samples = GTK_LABEL(gtk_builder_get_object(gui.builder, "background_label_samples"));
	bkg_keep_samples_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "subsky_keep_samples"));
	bkg_seq_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkBkgSeq"));
	bkg_seq_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryBkgSeq"));
	bkg_notebook = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "bkg_notebook_inter"));
	bkg_method_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "background_method_combo"));
	bkg_samples_expander = GTK_WIDGET(gtk_builder_get_object(gui.builder, "bkg_samples_expander"));
	bkg_auto_expander = GTK_WIDGET(gtk_builder_get_object(gui.builder, "bkg_auto_expander"));
	ag_scale_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ag_scale"));
	ag_smoothness_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ag_smoothness"));
	ag_protect_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "bkg_ag_protect_check"));
	ag_protect_threshold_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ag_protect_threshold"));
	ag_protect_amount_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ag_protect_amount"));
	ag_protect_threshold_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_ag_protect_threshold"));
	ag_protect_amount_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_ag_protect_amount"));
	ag_simplified_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "bkg_ag_simplified_check"));
	ag_degree_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ag_degree"));
	ag_degree_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_ag_degree"));
	ag_downsample_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "bkg_ag_downsample_combo"));
}

static poly_order get_poly_order() {
	return gtk_drop_down_get_selected(bkg_poly_order_combo);
}

static background_correction get_correction_type() {
	return gtk_drop_down_get_selected(bkg_correction_combo);
}

static int get_nb_samples_per_line() {
	return gtk_spin_button_get_value_as_int(bkg_nb_samples_spin);
}

static double get_tolerance_value() {
	return gtk_range_get_value(bkg_tolerance_scale);
}

static background_interpolation get_interpolation_method() {
	return gtk_drop_down_get_selected(bkg_interp_combo);
}

static double get_smoothing_parameter() {
	return gtk_spin_button_get_value(bkg_smoothing_spin);
}

static gboolean is_dither_checked() {
	return siril_toggle_get_active(GTK_WIDGET(bkg_dither_btn));
}

static gboolean is_randomize_checked() {
	return siril_toggle_get_active(GTK_WIDGET(bkg_randomize_btn));
}

static gboolean is_grad_descent_checked() {
	return siril_toggle_get_active(GTK_WIDGET(bkg_grad_descent_btn));
}

static background_method get_background_method() {
	return gtk_drop_down_get_selected(bkg_method_combo);
}

void update_bkg_compute_button_sensitivity(void) {
	background_extraction_init_statics();
	if (!bkg_compute_bkg_button)
		return;
	gboolean can_compute = (get_background_method() == BACKGROUND_METHOD_AUTO)
			|| (com.grad_samples != NULL);
	gtk_widget_set_sensitive(bkg_compute_bkg_button, can_compute);
}

static int get_ag_downsample() {
	static const int vals[] = {8, 4, 2, 1};
	guint sel = gtk_drop_down_get_selected(ag_downsample_combo);
	if (sel > 3) sel = 1;
	return vals[sel];
}

static void fill_autograd_from_ui(struct autograd_data *ag) {
	ag->scale = gtk_spin_button_get_value(ag_scale_spin);
	ag->smoothness = gtk_spin_button_get_value(ag_smoothness_spin);
	ag->protect = siril_toggle_get_active(GTK_WIDGET(ag_protect_check));
	ag->protect_threshold = gtk_spin_button_get_value(ag_protect_threshold_spin);
	ag->protect_amount = gtk_spin_button_get_value(ag_protect_amount_spin);
	ag->simplified = siril_toggle_get_active(GTK_WIDGET(ag_simplified_check));
	ag->degree = gtk_spin_button_get_value_as_int(ag_degree_spin);
	ag->downsample = get_ag_downsample();
}

/* management of the graphical state and backup image.
 * The original image is held by the shared preview backup (copy_gfit_to_backup);
 * bkg_processed keeps a copy of the last corrected result so the "View" selector
 * can switch between processed / background model / original at will. */
static fits bkg_processed = { 0 };
static background_correction bkg_view_correction = BACKGROUND_CORRECTION_SUBTRACT;
static gboolean background_computed = FALSE;

static void background_startup() {
	copy_gfit_to_backup();
}

/* Called on the GTK thread once a background has been computed and applied to
 * gfit. Stashes the result and arms the view selector (reset to "Processed"). */
static void bkg_on_computed(void) {
	background_computed = TRUE;
	clearfits(&bkg_processed);
	copyfits(gfit, &bkg_processed, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	bkg_view_correction = get_correction_type();
	g_rw_lock_reader_lock(&gfit->rwlock);
	gfit_modified_update_gui();
	g_rw_lock_reader_unlock(&gfit->rwlock);
	gtk_widget_set_sensitive(bkg_ok_button, TRUE);
	if (bkg_view_box)
		gtk_widget_set_sensitive(bkg_view_box, TRUE);
	if (bkg_view_combo)
		gtk_drop_down_set_selected(bkg_view_combo, 0);   /* show the processed result */
}

gboolean end_background(gpointer p) {
	struct background_data *args = (struct background_data *)p;
	stop_processing_thread();
	if (args) {
		bkg_on_computed();
		free(args);
	}
	set_cursor_waiting(FALSE);
	return FALSE;
}

/* Idle function for generic_image_worker path of on_bkg_compute_bkg_clicked */
static gboolean background_idle(gpointer p) {
	stop_processing_thread();
	bkg_on_computed();
	free_generic_img_args((struct generic_img_args *)p);
	set_cursor_waiting(FALSE);
	return FALSE;
}

/************* CALLBACKS *************/

/* Data for the sample generation worker thread */
struct bkg_generate_args {
	int nb_of_samples;
	double tolerance;
	gboolean randomize;
	gboolean grad_descent;
	rectangle sel; /* copy of com.selection captured on the GTK thread at click time */
};

/* Idle: runs on the GTK thread after sample generation completes */
static gboolean bkg_generate_idle(gpointer p) {
	struct bkg_generate_args *args = (struct bkg_generate_args *)p;
	stop_processing_thread();
	if (!args) {
		/* generation failed — log tab was already switched in the worker */
		update_bkg_compute_button_sensitivity();
		redraw(REDRAW_OVERLAY);
		set_cursor_waiting(FALSE);
		return FALSE;
	}
	free(args);
	update_bkg_compute_button_sensitivity();
	redraw(REDRAW_OVERLAY);
	set_cursor_waiting(FALSE);
	return FALSE;
}

/* Worker: runs in the processing thread, reads gfit to compute samples */
static gpointer bkg_generate_worker(gpointer p) {
	struct bkg_generate_args *args = (struct bkg_generate_args *)p;
	const rectangle *sel = (args->sel.w > 0 && args->sel.h > 0) ? &args->sel : NULL;
	int retval = generate_background_samples(args->nb_of_samples, args->tolerance, args->randomize, args->grad_descent, sel);
	if (retval) {
		/* Pass NULL to signal failure; args is freed here */
		free(args);
		siril_add_idle(bkg_generate_idle, NULL);
	} else {
		siril_add_idle(bkg_generate_idle, args);
	}
	return GINT_TO_POINTER(retval);
}

void on_bkg_randomize_button_toggled(GtkCheckButton *button, gpointer user_data) {
	background_extraction_init_statics();
	GtkLabel *label = bkg_label_samples;
	if (siril_toggle_get_active(GTK_WIDGET(button)))
		gtk_label_set_text(label, _("Number of samples"));
	else
		gtk_label_set_text(label, _("Samples per line"));
}

void on_background_generate_clicked(GtkButton *button, gpointer user_data) {
	background_extraction_init_statics();
	GtkCheckButton *keep_all_button = bkg_keep_samples_btn;
	gboolean keep_all = siril_toggle_get_active(GTK_WIDGET(keep_all_button));

	struct bkg_generate_args *args = calloc(1, sizeof(struct bkg_generate_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	args->nb_of_samples = get_nb_samples_per_line();
	args->tolerance = keep_all ? -1.0 : get_tolerance_value();
	args->randomize = is_randomize_checked();
	args->grad_descent = is_grad_descent_checked();
	args->sel = com.selection; /* capture on GTK thread before worker starts */

	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(bkg_generate_worker, args)) {
		free(args);
		set_cursor_waiting(FALSE);
	}
}

void on_background_clear_all_clicked(GtkButton *button, gpointer user_data) {
	sample_mutex_lock();
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	sample_mutex_unlock();

	update_bkg_compute_button_sensitivity();
	redraw(REDRAW_OVERLAY);
	set_cursor_waiting(FALSE);
}

void on_bkg_compute_bkg_clicked(GtkButton *button, gpointer user_data) {
	background_method method = get_background_method();

	if (method == BACKGROUND_METHOD_AUTO) {
		set_cursor_waiting(TRUE);
		copy_backup_to_gfit();

		struct background_data *bkg_args = calloc(1, sizeof(struct background_data));
		bkg_args->destroy_fn = free_background_data;
		bkg_args->method = BACKGROUND_METHOD_AUTO;
		bkg_args->threads = com.max_thread;
		bkg_args->from_ui = TRUE;
		bkg_args->correction = get_correction_type();
		bkg_args->dither = is_dither_checked();
		bkg_args->fit = gfit;
		fill_autograd_from_ui(&bkg_args->autograd);

		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		args->fit = gfit;
		args->mem_ratio = 6.0f;
		args->image_hook = remove_gradient_image_hook;
		args->log_hook = remove_gradient_log_hook;
		args->idle_function = background_idle;
		args->description = _("Automatic gradient removal");
		args->verbose = TRUE;
		/* Compute only previews the result; the undo state is saved on Apply.
		 * Without this the generic worker would push an undo entry right away,
		 * so a Close after Compute would leave a spurious undo step. */
		args->skip_generic_undo = TRUE;
		args->user = bkg_args;
		if (!start_in_new_thread(generic_image_worker, args)) {
			free_generic_img_args(args);
		}
		return;
	}

	if (com.grad_samples == NULL) {
		return;
	}
	set_cursor_waiting(TRUE);
	copy_backup_to_gfit();

	background_correction correction = get_correction_type();
	poly_order degree = get_poly_order();
	gboolean use_dither = is_dither_checked();
	double smoothing = get_smoothing_parameter();
	background_interpolation interpolation_method = get_interpolation_method();

	// Check if the image has a Bayer CFA pattern
	sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit->keywords.bayer_pattern);
	gboolean is_cfa = gfit->naxes[2] == 1 && pattern >= BAYER_FILTER_MIN && pattern <= BAYER_FILTER_MAX;

	struct background_data *bkg_args = calloc(1, sizeof(struct background_data));
	bkg_args->destroy_fn = free_background_data;
	bkg_args->threads = com.max_thread;
	bkg_args->from_ui = TRUE;
	bkg_args->correction = correction;
	bkg_args->interpolation_method = interpolation_method;
	bkg_args->degree = (poly_order)degree;
	bkg_args->smoothing = smoothing;
	bkg_args->dither = use_dither;
	bkg_args->fit = gfit;
	bkg_args->is_cfa = is_cfa;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = gfit;
	args->mem_ratio = 2.0f;
	args->image_hook = remove_gradient_image_hook;
	args->log_hook = remove_gradient_log_hook;
	args->idle_function = background_idle;
	args->description = _("Background extraction");
	args->verbose = TRUE;
	/* Compute only previews the result; the undo state is saved on Apply. */
	args->skip_generic_undo = TRUE;
	args->user = bkg_args;
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
	}
}

void on_background_ok_button_clicked(GtkButton *button, gpointer user_data) {
	background_extraction_init_statics();
	GtkCheckButton *seq_button = bkg_seq_btn;
	if (siril_toggle_get_active(GTK_WIDGET(seq_button)) && sequence_is_loaded()
			&& get_background_method() == BACKGROUND_METHOD_AUTO) {
		struct background_data *args = calloc(1, sizeof(struct background_data));
		args->method = BACKGROUND_METHOD_AUTO;
		args->correction = get_correction_type();
		args->dither = is_dither_checked();
		args->threads = com.max_thread;
		fill_autograd_from_ui(&args->autograd);
		set_cursor_waiting(TRUE);
		args->seqEntry = strdup(gtk_editable_get_text(GTK_EDITABLE(bkg_seq_entry)));
		if (args->seqEntry && args->seqEntry[0] == '\0') {
			free(args->seqEntry);
			args->seqEntry = strdup("bkg_");
		}
		args->seq = &com.seq;
		siril_toggle_set_active(GTK_WIDGET(seq_button), FALSE);
		apply_background_extraction_to_sequence(args);
	} else if (siril_toggle_get_active(GTK_WIDGET(seq_button)) && sequence_is_loaded()) {
		struct background_data *args = calloc(1, sizeof(struct background_data));
		args->nb_of_samples = get_nb_samples_per_line();
		args->tolerance = get_tolerance_value();
		args->correction = get_correction_type();
		args->degree = get_poly_order();
		args->smoothing = get_smoothing_parameter();
		args->dither = is_dither_checked();
		args->interpolation_method = get_interpolation_method();

		if (args->interpolation_method == BACKGROUND_INTER_POLY && args->degree > BACKGROUND_POLY_1) {
			int confirm = siril_confirm_dialog(_("Polynomial order seems too high."),
					_("You are about to process a sequence of preprocessed files with "
						"a polynomial degree greater than 1. This is unlikely because such "
						"gradients are often linear and a correction with a polynomial "
						"function of degree 1 is probably enough."), _("Extract Background"));
			if (!confirm) {
				free(args);
				set_cursor_waiting(FALSE);
				return;
			}
		} else if (args->interpolation_method != BACKGROUND_INTER_POLY) {
			int confirm = siril_confirm_dialog(_("Using wrong interpolation method"),
					_("You are about to process a sequence of preprocessed files with an RBF algorithm. "
						"This algorithm may not be very well suited for automated processing "
						"and we advise you to use the polynomial algorithm with a "
						"degree order of 1."), _("Extract Background"));
			if (!confirm) {
				free(args);
				set_cursor_waiting(FALSE);
				return;
			}
		}

		set_cursor_waiting(TRUE);

		args->seqEntry = strdup(gtk_editable_get_text(GTK_EDITABLE(bkg_seq_entry)));
		if (args->seqEntry && args->seqEntry[0] == '\0') {
			free(args->seqEntry);
			args->seqEntry = strdup("bkg_");
		}
		args->seq = &com.seq;
		/* now we uncheck the button */
		siril_toggle_set_active(GTK_WIDGET(seq_button), FALSE);
		apply_background_extraction_to_sequence(args);
	} else {
		if (background_computed) {
			/* the view selector may be showing the background model or the
			 * original, so restore the corrected result into gfit before saving */
			if (bkg_view_combo && gtk_drop_down_get_selected(bkg_view_combo) != 0
					&& bkg_processed.naxes[0] == gfit->naxes[0]) {
				copyfits(&bkg_processed, gfit, CP_COPYA, -1);
				notify_gfit_data_modified();
			}
			undo_save_state(get_preview_gfit_backup(), _("Background extraction (Correction: %s)"),
					bkg_view_correction == BACKGROUND_CORRECTION_DIVIDE ? "Division" : "Subtraction");
			background_computed = FALSE;
			clear_backup();
			siril_close_dialog("background_extraction_dialog");
		} else {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("No Background model computed"),
					_("You must first compute the background model."));
		}
	}
}

void apply_background_cancel() {
	siril_close_dialog("background_extraction_dialog");
}

void on_background_close_button_clicked(GtkButton *button, gpointer user_data) {
	apply_background_cancel();
}

gboolean bge_hide_on_delete(GtkWidget *widget) {
	apply_background_cancel();
	return TRUE;
}

void on_background_extraction_dialog_hide(GtkWidget *widget, gpointer user_data) {
	sample_mutex_lock();
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	sample_mutex_unlock();
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	redraw(REDRAW_OVERLAY);

	if (background_computed) {
		siril_preview_hide();
		background_computed = FALSE;
	} else {
		clear_backup();
	}
	clearfits(&bkg_processed);
	gtk_widget_set_sensitive(bkg_ok_button, FALSE);
	if (bkg_view_box)
		gtk_widget_set_sensitive(bkg_view_box, FALSE);
	if (bkg_view_combo)
		gtk_drop_down_set_selected(bkg_view_combo, 0);
}

void on_background_extraction_dialog_show(GtkWidget *widget, gpointer user_data) {
	background_extraction_init_statics();
	mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	background_startup();
	bkg_sync_method();   /* show the panel/overlay matching the selected method */
}

void on_background_extraction_combo_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	background_extraction_init_statics();
	gtk_notebook_set_current_page(bkg_notebook, gtk_drop_down_get_selected(combo));
}

/* grey out the automatic-model settings that don't apply with the current checkboxes */
static void ag_sync_sensitivity(void) {
	gboolean protect = siril_toggle_get_active(GTK_WIDGET(ag_protect_check));
	gtk_widget_set_sensitive(GTK_WIDGET(ag_protect_threshold_spin), protect);
	gtk_widget_set_sensitive(GTK_WIDGET(ag_protect_amount_spin), protect);
	gtk_widget_set_sensitive(ag_protect_threshold_label, protect);
	gtk_widget_set_sensitive(ag_protect_amount_label, protect);
	gboolean simpl = siril_toggle_get_active(GTK_WIDGET(ag_simplified_check));
	gtk_widget_set_sensitive(GTK_WIDGET(ag_degree_spin), simpl);
	gtk_widget_set_sensitive(ag_degree_label, simpl);
}

/* show the parameter panel matching the selected method; the automatic model
 * places no samples, so the interactive sample overlay is turned off for it */
static void bkg_sync_method(void) {
	gboolean is_auto = get_background_method() == BACKGROUND_METHOD_AUTO;
	gtk_widget_set_visible(bkg_samples_expander, !is_auto);
	gtk_widget_set_visible(bkg_auto_expander, is_auto);
	if (is_auto)
		ag_sync_sensitivity();
	update_bkg_compute_button_sensitivity();
	mouse_status = is_auto ? MOUSE_ACTION_SELECT_REG_AREA : MOUSE_ACTION_DRAW_SAMPLES;
	redraw(REDRAW_OVERLAY);
}

void on_background_method_combo_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj; (void)pspec; (void)user_data;
	background_extraction_init_statics();
	bkg_sync_method();
}

void on_bkg_ag_protect_toggled(GtkCheckButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	background_extraction_init_statics();
	ag_sync_sensitivity();
}

void on_bkg_ag_simplified_toggled(GtkCheckButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	background_extraction_init_statics();
	ag_sync_sensitivity();
}

/* View selector: after a background has been computed, the user can display the
 * corrected result, the background model that was removed, or the original.
 * The background model is reconstructed on the fly from the original (kept in
 * the preview backup) and the processed result, so it needs no extra storage
 * from the compute path and works identically for every removal method.
 * Reconstruction is exact up to a cosmetic brightness offset:
 *   subtract: bg = original - processed + level
 *   divide:   bg = original / processed * level
 * with `level` the per-channel mean of the original, so the shown background
 * sits at the image's own brightness. */
static void bkg_build_background_view(void) {
	fits *orig = get_preview_gfit_backup();
	if (!orig || (!orig->fdata && !orig->data))
		return;
	if (gfit->naxes[0] != orig->naxes[0] || gfit->naxes[1] != orig->naxes[1]
			|| gfit->naxes[2] != orig->naxes[2] || gfit->type != orig->type)
		return;
	if (bkg_processed.naxes[0] != orig->naxes[0] || bkg_processed.type != orig->type)
		return;

	const size_t npix = gfit->naxes[0] * gfit->naxes[1];
	const int nchan = gfit->naxes[2];
	const gboolean divide = (bkg_view_correction == BACKGROUND_CORRECTION_DIVIDE);

	for (int c = 0; c < nchan; c++) {
		if (gfit->type == DATA_FLOAT) {
			const float *o = orig->fpdata[c];
			const float *p = bkg_processed.fpdata[c];
			float *out = gfit->fpdata[c];
			double sum = 0.0;
			for (size_t i = 0; i < npix; i++) sum += o[i];
			const float level = (float)(sum / npix);
			for (size_t i = 0; i < npix; i++) {
				float v = divide ? o[i] / (p[i] > 1e-6f ? p[i] : 1e-6f) * level
								 : o[i] - p[i] + level;
				out[i] = v < 0.f ? 0.f : (v > 1.f ? 1.f : v);
			}
		} else {
			const WORD *o = orig->pdata[c];
			const WORD *p = bkg_processed.pdata[c];
			WORD *out = gfit->pdata[c];
			double sum = 0.0;
			for (size_t i = 0; i < npix; i++) sum += o[i];
			const double level = sum / npix;
			for (size_t i = 0; i < npix; i++) {
				double v = divide ? (double)o[i] / (p[i] > 0 ? p[i] : 1) * level
								  : (double)o[i] - p[i] + level;
				out[i] = v < 0.0 ? 0 : (v > USHRT_MAX ? USHRT_MAX : (WORD)(v + 0.5));
			}
		}
	}
	invalidate_stats_from_fit(gfit);
}

void on_bkg_view_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj; (void)pspec; (void)user_data;
	background_extraction_init_statics();
	if (!background_computed)
		return;
	switch (gtk_drop_down_get_selected(bkg_view_combo)) {
		case 0:  /* processed result */
			copyfits(&bkg_processed, gfit, CP_COPYA, -1);
			break;
		case 1:  /* background model (reconstructed) */
			bkg_build_background_view();
			break;
		default: /* original image */
			copy_backup_to_gfit();
			break;
	}
	notify_gfit_data_modified();
	gfit_modified_update_gui();
}

void on_checkBkgSeq_toggled(GtkCheckButton *button, gpointer user_data) {
	background_extraction_init_statics();
	GtkWidget *ok = bkg_ok_button;
	if (siril_toggle_get_active(GTK_WIDGET(button))) {
		gtk_widget_set_sensitive(ok, TRUE);
	} else {
		gboolean have_model = get_background_method() == BACKGROUND_METHOD_AUTO
				? background_computed : ((com.grad_samples != NULL) && background_computed);
		gtk_widget_set_sensitive(ok, have_model);
	}
}
