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
#include <string.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "registration/registration.h"
#include "registration/mpp_config.h"
#include "stacking/sum.h"
#include "stacking/stacking.h"

#include "gui/stacking.h"

static struct stacking_args stackparam = { 0 };

static struct filtering_tuple stackfilters[MAX_FILTERS];
/* Values for seq filtering for % or k value*/
static float filter_initvals[] = {90., 3.}; // max %, max k
static float filter_maxvals[] = {100., 5.}; // max %, max k
static float filter_increments[] = {1., 0.1}; // spin button steps for % and k
// update_adjustment passed here as a static (instead of function parameter like in registration)
/* Forward decls for the drizzle-combo CFA-aware filters; the definitions
 * live further down with on_combo_mpp_drizzle_changed but start_stacking
 * installs the cell-data-func at widget init. */
static void mpp_drizzle_combo_cell_func(GtkCellLayout *, GtkCellRenderer *,
                                         GtkTreeModel *, GtkTreeIter *, gpointer);
static void mpp_drizzle_combo_validate(GtkComboBox *);

// in order not to mess up all the calls to update_stack_interface
static int update_adjustment = -1;
static char *filter_tooltip_text[] = {
	N_("Percents of the images of the sequence."),
	N_("Number of standard deviations for rejection of worst images by k-sigma clipping algorithm.")
};

stack_method stacking_methods[] = {
	stack_summing_generic, stack_mean_with_rejection, stack_median, stack_addmax, stack_addmin,
	stack_mpp_handler
};

void initialize_stacking_methods() {
	init_stacking_args(&stackparam);
	GtkComboBoxText *stackcombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboboxstack_methods"));
	GtkComboBoxText *rejectioncombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comborejection"));
	GtkComboBoxText *weightingcombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboweighing"));
	GtkSpinButton *low = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
	GtkSpinButton *high = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
	gtk_combo_box_set_active(GTK_COMBO_BOX(stackcombo), com.pref.stack.method);
	gtk_combo_box_set_active(GTK_COMBO_BOX(rejectioncombo), com.pref.stack.rej_method);
	gtk_combo_box_set_active(GTK_COMBO_BOX(weightingcombo), com.pref.stack.weighting_method);
	switch (gtk_combo_box_get_active(GTK_COMBO_BOX(rejectioncombo))) {
	case PERCENTILE:
		gtk_spin_button_set_value(low, com.pref.stack.percentile_low);
		gtk_spin_button_set_value(high, com.pref.stack.percentile_high);
		break;
	case LINEARFIT:
		gtk_spin_button_set_value(low, com.pref.stack.linear_low);
		gtk_spin_button_set_value(high, com.pref.stack.linear_high);
		break;
	case SIGMA:
	case MAD:
	case SIGMEDIAN:
	case WINSORIZED:
		gtk_spin_button_set_value(low, com.pref.stack.sigma_low);
		gtk_spin_button_set_value(high, com.pref.stack.sigma_high);
		break;
	case GESDT:
		gtk_spin_button_set_value(low, 0.3);
		gtk_spin_button_set_value(high, 0.05);
		break;
	default:
		return;
	}
}

/* starts a summing operation using data stored in the stackparam structure
 * function is not reentrant but can be called again after it has returned and the thread is running */
static void start_stacking() {
	gchar *error = NULL;
	static GtkComboBox *method_combo = NULL, *rejec_combo = NULL, *norm_combo = NULL, *weighing_combo;
	static GtkEntry *output_file = NULL;
	static GtkToggleButton *overwrite = NULL, *force_norm = NULL, *max_framing = NULL,
					*fast_norm = NULL, *rejmaps = NULL, *merge_rejmaps = NULL, 
					*upscale_at_stacking = NULL, *overlap_norm = NULL, *force32b = NULL;
	static GtkSpinButton *sigSpin[2] = {NULL, NULL}, *feather_dist = NULL;
	static GtkWidget *norm_to_max = NULL, *RGB_equal = NULL, *blend_frame = NULL;
	static GtkComboBox *mpp_drizzle_combo = NULL, *mpp_driz_kernel_combo = NULL;
	static GtkSpinButton *mpp_stack_percent = NULL, *mpp_stack_frames = NULL,
	                     *mpp_bg_fraction = NULL, *mpp_bg_blend = NULL,
	                     *mpp_pixfrac = NULL;

	if (method_combo == NULL) {
		method_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "comboboxstack_methods"));
		output_file = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryresultfile"));
		overwrite = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "checkbutoverwrite"));
		sigSpin[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
		sigSpin[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
		rejec_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "comborejection"));
		norm_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combonormalize"));
		weighing_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "comboweighing"));
		force_norm = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "checkforcenorm"));
		fast_norm = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "checkfastnorm"));
		max_framing = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "check_maximize_framing"));
		norm_to_max = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max"));
		RGB_equal = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal"));
		rejmaps = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "rejmaps_checkbutton"));
		merge_rejmaps = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
		upscale_at_stacking = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "check_upscale_at_stacking"));
		feather_dist = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_stack_feather_dist"));
		blend_frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_blend_frame"));
		overlap_norm = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "check_norm_overlap"));
		force32b = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "check_force32b"));
		/* STACK_MPP widgets */
		mpp_drizzle_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_mpp_drizzle"));
		mpp_driz_kernel_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_mpp_driz_kernel"));
		/* Install the per-row sensitivity callback so stsci-* greys out
		 * for CFA input and bayer-* greys out for mono input. */
		GList *cells = gtk_cell_layout_get_cells(GTK_CELL_LAYOUT(mpp_drizzle_combo));
		if (cells && cells->data) {
			gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(mpp_drizzle_combo),
			                                   GTK_CELL_RENDERER(cells->data),
			                                   mpp_drizzle_combo_cell_func, NULL, NULL);
		}
		g_list_free(cells);
		mpp_stack_percent = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_stack_percent"));
		mpp_stack_frames  = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_stack_frames"));
		mpp_bg_fraction   = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_bg_fraction"));
		mpp_bg_blend      = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_bg_blend"));
		mpp_pixfrac       = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_pixfrac"));
	}

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	stackparam.sig[0] = (float) gtk_spin_button_get_value(sigSpin[0]);
	stackparam.sig[1] = (float) gtk_spin_button_get_value(sigSpin[1]);
	stackparam.type_of_rejection = gtk_combo_box_get_active(rejec_combo);
	stackparam.create_rejmaps = gtk_toggle_button_get_active(rejmaps);
	stackparam.merge_lowhigh_rejmaps = gtk_toggle_button_get_active(merge_rejmaps);
	stackparam.normalize = gtk_combo_box_get_active(norm_combo);
	stackparam.force_norm = gtk_toggle_button_get_active(force_norm);
	stackparam.output_norm = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(norm_to_max)) && gtk_widget_is_visible(norm_to_max);
	stackparam.maximize_framing = gtk_toggle_button_get_active(max_framing) && gtk_widget_is_visible(GTK_WIDGET(max_framing));
	stackparam.upscale_at_stacking = gtk_toggle_button_get_active(upscale_at_stacking) && gtk_widget_is_visible(GTK_WIDGET(upscale_at_stacking));
	stackparam.coeff.offset = NULL;
	stackparam.coeff.mul = NULL;
	stackparam.coeff.scale = NULL;
	stackparam.method =	stacking_methods[gtk_combo_box_get_active(method_combo)];
	gboolean weighing_is_enabled = gtk_widget_get_visible(GTK_WIDGET(weighing_combo));
	if (weighing_is_enabled) {
		int weight_type = gtk_combo_box_get_active(weighing_combo);
		if (weight_type == NOISE_WEIGHT) {
			int norm_type = gtk_combo_box_get_active(norm_combo);
			if (norm_type == NO_NORM) {
				siril_log_color_message(_("Weighting by noise is allowed only if normalization has been activated, ignoring weights.\n"), "red");
				weight_type = NO_WEIGHT;
			} else
				stackparam.weighting_type = norm_type;
		} else {
			stackparam.weighting_type = weight_type;
		}
	}
	stackparam.equalizeRGB = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(RGB_equal)) && gtk_widget_is_visible(RGB_equal)  && (gtk_combo_box_get_active(norm_combo) != NO_NORM);
	stackparam.lite_norm = gtk_toggle_button_get_active(fast_norm);
	stackparam.feather_dist = (int)gtk_spin_button_get_value(feather_dist) * gtk_widget_get_visible(blend_frame);
	stackparam.overlap_norm = gtk_toggle_button_get_active(overlap_norm) * gtk_widget_get_visible(blend_frame);
	if (stackparam.overlap_norm && stackparam.weighting_type == NOISE_WEIGHT) {
		siril_log_color_message(_("Weighting by noise cannot be used with overlap normalization, ignoring weights.\n"), "red");
		stackparam.weighting_type = NO_WEIGHT;
	}

	stackparam.use_32bit_output = gtk_toggle_button_get_active(force32b) || evaluate_stacking_should_output_32bits(stackparam.method,
			&com.seq, stackparam.nb_images_to_stack, &error);
	if (error) {
		siril_log_color_message(error, "red");
		return;
	}

	// ensure we have no normalization if not supported by the stacking method
	if (stackparam.method != stack_median && stackparam.method != stack_mean_with_rejection)
		stackparam.normalize = NO_NORM;
	stackparam.seq = &com.seq;
	stackparam.reglayer = get_registration_layer(stackparam.seq);
	/* STACK_MPP doesn't go through the homography-aware regdata path — it
	 * reads its own per-AP per-frame shifts from the .mpp sidecar. The
	 * "registration data with more than simple shifts" prompt would fire
	 * for a sidecar's quality-only regdata so we skip it for STACK_MPP. */
	if (stackparam.method != stack_mpp_handler) {
		// checking regdata is absent, or if present, is only shift
		if (!test_regdata_is_valid_and_shift(stackparam.seq, stackparam.reglayer)) {
			int confirm = siril_confirm_dialog(_("Registration data found"),
				_("Stacking has detected registration data with more than simple shifts.\n"
				"Normally, you should apply existing registration before stacking."),
				_("Stack anyway"));
			if (!confirm)
				return;
		}
	}

	if (stackparam.overlap_norm && stackparam.nb_images_to_stack > 20) {
		gchar *onorm_msg = g_strdup_printf(_("You have chosen to compute normalization on overlaps with more than %d images.\n"
			"This option should normally be used to stitch stacked mosaic tiles, not subs.\n"
			"If you proceed, execution may be slow"), MAX_IMAGES_FOR_OVERLAP);
		int confirm = siril_confirm_dialog(_("Large number of images"),
			onorm_msg,
			_("Stack anyway"));
		g_free(onorm_msg);
		if (!confirm)
			return;
	}

	/* Do not display that cause it uses the generic function that already
	 * displays this text
	 */
	if (stackparam.method != &stack_summing_generic)
		siril_log_color_message(_("Stacking: processing...\n"), "green");
	gettimeofday(&stackparam.t_start, NULL);
	set_cursor_waiting(TRUE);

	stackparam.output_overwrite = gtk_toggle_button_get_active(overwrite);
	stackparam.output_filename = gtk_entry_get_text(output_file);

	/* Stacking. Result is in gfit if success */
	struct stacking_args *params = calloc(1, sizeof(struct stacking_args));
	stacking_args_deep_copy(&stackparam, params);

	/* For STACK_MPP, capture the stack-side widget values into a fresh
	 * mpp_config_t and hang it off params. stack_mpp_handler reads from
	 * params->mpp_cfg and stacking_args_deep_free will release it.
	 * Allocated AFTER deep_copy so stackparam never owns the pointer
	 * (avoids a double-free if the user clicks Stack twice in quick
	 * succession). */
	if (stackparam.method == stack_mpp_handler) {
		mpp_config_t *cfg = calloc(1, sizeof(*cfg));
		mpp_config_defaults(cfg);
		const int driz_idx = gtk_combo_box_get_active(mpp_drizzle_combo);
		switch (driz_idx) {
			/* Index → (mode, factor) — must match the GtkComboBoxText item
			 * order in siril.ui. 0..3 are the Phase 5a bicubic entries
			 * (legacy); 4..7 are STScI / Bayer variants from Phase 5b. */
			case 0: cfg->drizzle_mode = MPP_DRIZZLE_OFF;     cfg->drizzle_factor = 1; break;
			case 1: cfg->drizzle_mode = MPP_DRIZZLE_BICUBIC; cfg->drizzle_factor = 3; break;  /* 1.5x */
			case 2: cfg->drizzle_mode = MPP_DRIZZLE_BICUBIC; cfg->drizzle_factor = 2; break;
			case 3: cfg->drizzle_mode = MPP_DRIZZLE_BICUBIC; cfg->drizzle_factor = 3; break;
			case 4: cfg->drizzle_mode = MPP_DRIZZLE_STSCI;   cfg->drizzle_factor = 2; break;
			case 5: cfg->drizzle_mode = MPP_DRIZZLE_STSCI;   cfg->drizzle_factor = 3; break;
			case 6: cfg->drizzle_mode = MPP_DRIZZLE_BAYER;   cfg->drizzle_factor = 2; break;
			case 7: cfg->drizzle_mode = MPP_DRIZZLE_BAYER;   cfg->drizzle_factor = 3; break;
			default: cfg->drizzle_mode = MPP_DRIZZLE_OFF;    cfg->drizzle_factor = 1; break;
		}
		cfg->alignment_points_frame_percent          = gtk_spin_button_get_value_as_int(mpp_stack_percent);
		cfg->alignment_points_frame_number           = gtk_spin_button_get_value_as_int(mpp_stack_frames);
		cfg->stack_frames_background_fraction        = gtk_spin_button_get_value(mpp_bg_fraction);
		cfg->stack_frames_background_blend_threshold = gtk_spin_button_get_value(mpp_bg_blend);
		cfg->drizzle_pixfrac = gtk_spin_button_get_value(mpp_pixfrac);
		const int k = gtk_combo_box_get_active(mpp_driz_kernel_combo);
		cfg->drizzle_kernel = (k >= 0 && k <= MPP_KERNEL_LANCZOS3) ? k : MPP_KERNEL_TURBO;
		params->mpp_cfg = cfg;
	}

	if (!start_in_new_thread(stack_function_handler, params)) {
		stacking_args_deep_free(params);
	}
}

void on_seqstack_button_clicked (GtkButton *button, gpointer user_data){
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_stacking();
}

/* Show / hide the STScI- and Bayer-only widgets (pixfrac spinner + label,
 * driz kernel combo + label) based on the drizzle-mode combo selection.
 * Indices 0..3 are the bicubic family (Off / 1.5x / 2x / 3x); 4..7 are
 * the dobox-driven modes that consume pixfrac and kernel. */
void on_combo_mpp_drizzle_changed(GtkComboBox *combo, gpointer user_data) {
	(void) user_data;
	static GtkWidget *lbl_pf = NULL, *spin_pf = NULL;
	static GtkWidget *lbl_kn = NULL, *combo_kn = NULL;
	if (!lbl_pf) {
		lbl_pf   = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_mpp_pixfrac"));
		spin_pf  = GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_mpp_pixfrac"));
		lbl_kn   = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_mpp_driz_kernel"));
		combo_kn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combo_mpp_driz_kernel"));
	}
	const int idx = gtk_combo_box_get_active(combo);
	const gboolean dobox_mode = (idx >= 4);   /* STScI 2x/3x or Bayer 2x/3x */
	gtk_widget_set_visible(lbl_pf,   dobox_mode);
	gtk_widget_set_visible(spin_pf,  dobox_mode);
	gtk_widget_set_visible(lbl_kn,   dobox_mode);
	gtk_widget_set_visible(combo_kn, dobox_mode);
}

/* Cell-renderer callback: gray out drizzle-mode combo items that are
 * incompatible with the loaded sequence's channel layout. STScI modes
 * (indices 4–5) require mono input — debayer-then-STScI is rejected at
 * the command-handler level because classical debayer artefacts would
 * be amplified by drizzle. Bayer modes (indices 6–7) require CFA input.
 * The bicubic family (0–3) is always available.
 *
 * Sensitive=FALSE makes the row unclickable in the dropdown but keeps it
 * visible so the user sees the available palette of options. The
 * mpp_drizzle_combo_validate function (called from update_stack_interface)
 * resets the active selection if the previously-selected mode became
 * invalid on a sequence-state change. */
static void mpp_drizzle_combo_cell_func(GtkCellLayout *cell_layout,
                                         GtkCellRenderer *cell,
                                         GtkTreeModel *tree_model,
                                         GtkTreeIter *iter,
                                         gpointer data) {
	(void)cell_layout; (void)data;
	GtkTreePath *path = gtk_tree_model_get_path(tree_model, iter);
	if (!path) return;
	const int idx = gtk_tree_path_get_indices(path)[0];
	gtk_tree_path_free(path);

	const gboolean seqloaded = sequence_is_loaded();
	const gboolean is_cfa  = seqloaded && com.seq.nb_layers > 1;
	const gboolean is_mono = seqloaded && com.seq.nb_layers == 1;

	gboolean sensitive = TRUE;
	if (idx == 4 || idx == 5) sensitive = !is_cfa;   /* stsci-* needs mono */
	if (idx == 6 || idx == 7) sensitive = !is_mono;  /* bayer-* needs CFA  */
	g_object_set(cell, "sensitive", sensitive, NULL);
}

/* Validate the currently-selected drizzle mode against the loaded
 * sequence. If the selection became invalid (e.g., user had stsci-2x
 * selected when no sequence was loaded, then opened a CFA SER), fall
 * back to "Off" so they don't hit a rejection at stack time. Called
 * from update_stack_interface on every sequence-state change. */
static void mpp_drizzle_combo_validate(GtkComboBox *combo) {
	if (!combo) return;
	const gboolean seqloaded = sequence_is_loaded();
	const gboolean is_cfa  = seqloaded && com.seq.nb_layers > 1;
	const gboolean is_mono = seqloaded && com.seq.nb_layers == 1;
	const int idx = gtk_combo_box_get_active(combo);
	gboolean invalid = FALSE;
	if ((idx == 4 || idx == 5) && is_cfa)  invalid = TRUE;
	if ((idx == 6 || idx == 7) && is_mono) invalid = TRUE;
	if (invalid) gtk_combo_box_set_active(combo, 0);
	/* Force the dropdown to re-evaluate cell sensitivities for the new
	 * sequence state. cell_data_func reads sequence state on each render. */
	gtk_widget_queue_draw(GTK_WIDGET(combo));
}

void on_comboboxstack_methods_changed (GtkComboBox *box, gpointer user_data) {
	static GtkNotebook* notebook = NULL;
	if (!notebook)
		notebook = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook4"));
	com.pref.stack.method = gtk_combo_box_get_active(box);

	gtk_notebook_set_current_page(notebook, com.pref.stack.method);
	update_stack_interface(TRUE);
}

void on_combonormalize_changed (GtkComboBox *box, gpointer user_data) {
	static GtkWidget *widgetnormalize = NULL, *force_norm = NULL, *fast_norm = NULL;
	if (!widgetnormalize) {
		widgetnormalize = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combonormalize"));
		force_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkforcenorm"));
		fast_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm"));
	}
	gtk_widget_set_sensitive(force_norm, gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
	gtk_widget_set_sensitive(fast_norm, gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
}

void on_comboweighing_changed (GtkComboBox *box, gpointer user_data) {
	com.pref.stack.weighting_method = gtk_combo_box_get_active(box);
}

void on_stack_siglow_button_value_changed(GtkSpinButton *button, gpointer user_data) {
	rejection type_of_rejection = gtk_combo_box_get_active((GtkComboBox *)user_data);

	switch (type_of_rejection) {
	case PERCENTILE:
		com.pref.stack.percentile_low = gtk_spin_button_get_value(button);
		break;
	case LINEARFIT:
		com.pref.stack.linear_low = gtk_spin_button_get_value(button);
		break;
	case SIGMA:
	case MAD:
	case SIGMEDIAN:
	case WINSORIZED:
		com.pref.stack.sigma_low = gtk_spin_button_get_value(button);
		break;
	case GESDT:
		break;
	default:
		return;
	}
}

void on_stack_sighigh_button_value_changed(GtkSpinButton *button, gpointer user_data) {
	rejection type_of_rejection = gtk_combo_box_get_active((GtkComboBox *)user_data);

	switch (type_of_rejection) {
	case PERCENTILE:
		com.pref.stack.percentile_high = gtk_spin_button_get_value(button);
		break;
	case LINEARFIT:
		com.pref.stack.linear_high = gtk_spin_button_get_value(button);
		break;
	case SIGMA:
	case MAD:
	case SIGMEDIAN:
	case WINSORIZED:
		com.pref.stack.sigma_high = gtk_spin_button_get_value(button);
		break;
	case GESDT:
		break;
	default:
		return;
	}
}

void on_comborejection_changed(GtkComboBox *box, gpointer user_data) {
	rejection type_of_rejection = gtk_combo_box_get_active(box);
	static GtkWidget *labellow = NULL, *labelhigh = NULL;
	static GtkWidget *siglow = NULL, *sighigh = NULL;
	static GtkWidget *rejmaps = NULL, *merge_rejmaps = NULL;

	if (!labellow) {
		labellow = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_low"));
		labelhigh = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_high"));
		siglow = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
		sighigh = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
		rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rejmaps_checkbutton"));
		merge_rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
	}

	g_signal_handlers_block_by_func(GTK_SPIN_BUTTON(siglow), on_stack_siglow_button_value_changed, NULL);
	g_signal_handlers_block_by_func(GTK_SPIN_BUTTON(sighigh), on_stack_sighigh_button_value_changed, NULL);
	/* set default values */
	switch (type_of_rejection) {
		case NO_REJEC:
			gtk_widget_set_visible(siglow, FALSE);
			gtk_widget_set_visible(rejmaps, FALSE);
			gtk_widget_set_visible(sighigh, FALSE);
			gtk_widget_set_visible(labellow, FALSE);
			gtk_widget_set_visible(labelhigh, FALSE);
			gtk_widget_set_visible(merge_rejmaps, FALSE);
			break;
		case PERCENTILE:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Low clipping factor for the percentile clipping rejection algorithm."));
			gtk_widget_set_tooltip_text(sighigh, _("High clipping factor for the percentile clipping rejection algorithm."));
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), com.pref.stack.percentile_low);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), com.pref.stack.percentile_high);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 1.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 1.0);
			gtk_label_set_text(GTK_LABEL(labellow), _("Percentile low: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("Percentile high: "));
			break;
		case LINEARFIT:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Tolerance for low pixel values of the linear fit clipping algorithm, in sigma units."));
			gtk_widget_set_tooltip_text(sighigh, _("Tolerance for high pixel values of the linear fit clipping algorithm, in sigma units."));
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 10.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), com.pref.stack.linear_low);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), com.pref.stack.linear_high);
			gtk_label_set_text(GTK_LABEL(labellow), _("Linear low: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("Linear high: "));
			break;
		default:
		case SIGMA:
		case MAD:
		case SIGMEDIAN:
		case WINSORIZED:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Low clipping factor for the sigma clipping rejection algorithm."));
			gtk_widget_set_tooltip_text(sighigh, _("High clipping factor for the sigma clipping rejection algorithm."));
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 10.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), com.pref.stack.sigma_low);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), com.pref.stack.sigma_high);
			gtk_label_set_text(GTK_LABEL(labellow), _("Sigma low: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("Sigma high: "));
			break;
		case GESDT:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Expected fraction of maximum outliers for the Generalized Extreme Studentized Deviate Test algorithm."));
			gtk_widget_set_tooltip_text(sighigh, _("Probability of making a false positive for the Generalized Extreme Studentized Deviate Test algorithm. "
					"Increasing this value will reject more pixels."));
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 1.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 1.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), 0.3);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), 0.05);
			gtk_label_set_text(GTK_LABEL(labellow), _("ESD Outliers: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("ESD Significance: "));
	}

	g_signal_handlers_unblock_by_func(GTK_SPIN_BUTTON(siglow), on_stack_siglow_button_value_changed, NULL);
	g_signal_handlers_unblock_by_func(GTK_SPIN_BUTTON(sighigh), on_stack_sighigh_button_value_changed, NULL);

	com.pref.stack.rej_method = gtk_combo_box_get_active(box);
}

void on_rejmaps_toggled(GtkToggleButton *button, gpointer user_data) {
	static GtkWidget *merge = NULL;
	if (!merge) merge = GTK_WIDGET(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
	gtk_widget_set_sensitive(merge, gtk_toggle_button_get_active(button));
}

void on_check_maximize_framing_toggled(GtkToggleButton *button, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_stacksel_changed(GtkComboBox *widget, gpointer user_data) {
	const gchar *caller = gtk_buildable_get_name(GTK_BUILDABLE (widget));
	if (g_str_has_prefix(caller, "filter_type")) {
		update_adjustment = (int)g_ascii_strtod(caller + 11, NULL) - 1; // filter_type1, 2 or 3 to be parsed as 0, 1 or 2
	}
	update_stack_interface(TRUE);
}

void on_spinbut_percent_change(GtkSpinButton *spinbutton, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_filter_add1_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter2 = NULL, *stackspin2 = NULL, *filter_add2 = NULL;
	static GtkWidget *filter_rem2 = NULL, *labelfilter2 = NULL, *filter_type2 = NULL;
	if (!combofilter2) {
		combofilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter2"));
		stackspin2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin2"));
		filter_add2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_add2"));
		filter_rem2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem2"));
		labelfilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter2"));
		filter_type2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type2"));
	}
	gtk_widget_set_visible(combofilter2, TRUE);
	gtk_widget_set_visible(stackspin2, TRUE);
	gtk_widget_set_visible(filter_add2, TRUE);
	gtk_widget_set_visible(filter_rem2, TRUE);
	gtk_widget_set_visible(labelfilter2, TRUE);
	gtk_widget_set_visible(filter_type2, TRUE);
	update_stack_interface(TRUE);
}

void on_filter_add2_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter3 = NULL, *stackspin3 = NULL, *filter_rem3 = NULL;
	static GtkWidget *labelfilter3 = NULL, *filter_type3 = NULL;
	if (!combofilter3) {
		combofilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter3"));
		stackspin3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin3"));
		filter_rem3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem3"));
		labelfilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter3"));
		filter_type3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type3"));
	}
	gtk_widget_set_visible(combofilter3, TRUE);
	gtk_widget_set_visible(stackspin3, TRUE);
	gtk_widget_set_visible(filter_rem3, TRUE);
	gtk_widget_set_visible(labelfilter3, TRUE);
	gtk_widget_set_visible(filter_type3, TRUE);
	update_stack_interface(TRUE);
}

void on_filter_rem2_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter2 = NULL, *stackspin2 = NULL, *filter_add2 = NULL;
	static GtkWidget *filter_rem2 = NULL, *labelfilter2 = NULL, *filter_type2 = NULL;
	if (!combofilter2) {
		combofilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter2"));
		stackspin2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin2"));
		filter_add2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_add2"));
		filter_rem2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem2"));
		labelfilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter2"));
		filter_type2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type2"));
	}
	gtk_widget_set_visible(combofilter2, FALSE);
	gtk_widget_set_visible(stackspin2, FALSE);
	gtk_widget_set_visible(filter_add2, FALSE);
	gtk_widget_set_visible(filter_rem2, FALSE);
	gtk_widget_set_visible(labelfilter2, FALSE);
	gtk_widget_set_visible(filter_type2, FALSE);
	update_stack_interface(TRUE);
}

void on_filter_rem3_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter3 = NULL, *stackspin3 = NULL, *filter_rem3 = NULL;
	static GtkWidget *labelfilter3 = NULL, *filter_type3 = NULL;
	if (!combofilter3) {
		combofilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter3"));
		stackspin3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin3"));
		filter_rem3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem3"));
		labelfilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter3"));
		filter_type3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type3"));
	}
	gtk_widget_set_visible(combofilter3, FALSE);
	gtk_widget_set_visible(stackspin3, FALSE);
	gtk_widget_set_visible(filter_rem3, FALSE);
	gtk_widget_set_visible(labelfilter3, FALSE);
	gtk_widget_set_visible(filter_type3, FALSE);
	update_stack_interface(TRUE);
}

void get_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter) {
	int filter, guifilter, channel = 0, type;
	gboolean is_ksig = FALSE;
	double percent = 0.0;
	static GtkComboBox *filter_combo[] = {NULL, NULL, NULL};
	static GtkAdjustment *stackadj[] = {NULL, NULL, NULL};
	static GtkWidget *spin[] = {NULL, NULL, NULL};
	static GtkWidget *ksig[] = {NULL, NULL, NULL};
	if (!spin[0]) {
		spin[0] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin1"));
		spin[1] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin2"));
		spin[2] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin3"));
		stackadj[0] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[0]));
		stackadj[1] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[1]));
		stackadj[2] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[2]));
		filter_combo[0] = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combofilter1"));
		filter_combo[1] = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combofilter2"));
		filter_combo[2] = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combofilter3"));
		ksig[0] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type1"));
		ksig[1] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type2"));
		ksig[2] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type3"));
	}
	for (filter = 0, guifilter = 0; guifilter < 3; guifilter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[guifilter]))) {
			continue;
		}

		type = gtk_combo_box_get_active(filter_combo[guifilter]);
		if (type != ALL_IMAGES && type != SELECTED_IMAGES) {
			channel = get_registration_layer(&com.seq);
			percent = gtk_adjustment_get_value(stackadj[guifilter]);
			is_ksig = gtk_combo_box_get_active(GTK_COMBO_BOX(ksig[guifilter]));
		}
		if (update_adjustment == filter) {
			g_signal_handlers_block_by_func(stackadj[filter], on_stacksel_changed, NULL);
			gtk_adjustment_set_upper(stackadj[filter], filter_maxvals[is_ksig]);
			gtk_adjustment_set_value(stackadj[filter], filter_initvals[is_ksig]);
			gtk_adjustment_set_step_increment(stackadj[filter], filter_increments[is_ksig]);
			gtk_widget_set_tooltip_text(spin[filter], filter_tooltip_text[is_ksig]);
			g_signal_handlers_unblock_by_func(stackadj[filter], on_stacksel_changed, NULL);
			update_adjustment = -1;
		}
		switch (type) {
			default:
			case ALL_IMAGES:
				stackfilters[filter].filter = seq_filter_all;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				gtk_widget_set_visible(ksig[guifilter], FALSE);
				break;
			case SELECTED_IMAGES:
				stackfilters[filter].filter = seq_filter_included;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				gtk_widget_set_visible(ksig[guifilter], FALSE);
				break;
			case BEST_PSF_IMAGES:
				stackfilters[filter].filter = seq_filter_fwhm;
				stackfilters[filter].param = compute_highest_accepted_fwhm(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_WPSF_IMAGES:
				stackfilters[filter].filter = seq_filter_weighted_fwhm;
				stackfilters[filter].param = compute_highest_accepted_weighted_fwhm(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_ROUND_IMAGES:
				stackfilters[filter].filter = seq_filter_roundness;
				stackfilters[filter].param = compute_lowest_accepted_roundness(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_BKG_IMAGES:
				stackfilters[filter].filter = seq_filter_background;
				stackfilters[filter].param = compute_highest_accepted_background(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_NBSTARS_IMAGES:
				stackfilters[filter].filter = seq_filter_nbstars;
				stackfilters[filter].param = compute_lowest_accepted_nbstars(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_QUALITY_IMAGES:
				stackfilters[filter].filter = seq_filter_quality;
				stackfilters[filter].param = compute_lowest_accepted_quality(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
		}
		filter++;
	}
	stackfilters[filter].filter = NULL;

	if (filter == 1) {
		*filtering_criterion = stackfilters[0].filter;
		*filtering_parameter = stackfilters[0].param;
	} else {
		*filtering_criterion = create_multiple_filter_from_list(stackfilters);
		*filtering_parameter = 0.0;
	}
}

static void update_filter_label() {
	static GtkComboBox *filter_combo[3] = { NULL };
	static GtkLabel *filter_label[3] = { NULL };
	if (!filter_combo[0]) {
		filter_combo[0] = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combofilter1"));
		filter_combo[1] = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combofilter2"));
		filter_combo[2] = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combofilter3"));
		filter_label[0] = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter1"));
		filter_label[1] = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter2"));
		filter_label[2] = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter3"));
	}

	for (int filter = 0; filter < 3; filter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[filter]))) {
			break;
		}

		int type = gtk_combo_box_get_active(filter_combo[filter]);
		double param = stackfilters[filter].param;
		gchar *filter_str;
		if (param == DBL_MIN || param == DBL_MAX || param == 0.0) {
			if (type == ALL_IMAGES || type == SELECTED_IMAGES)
				filter_str = g_strdup("");
			else filter_str = g_strdup("N/A");
		} else {
			switch (type) {
			default:
			case ALL_IMAGES:
			case SELECTED_IMAGES:
				filter_str = g_strdup("");
				break;
			case BEST_PSF_IMAGES:
			case BEST_WPSF_IMAGES:
				filter_str = g_strdup_printf("< %.2lf", param);
				break;
			case BEST_BKG_IMAGES :
				filter_str = (param < 1.) ? g_strdup_printf("< %.5lf", param) : g_strdup_printf("< %d", (int)param);
				break;
			case BEST_ROUND_IMAGES:
			case BEST_QUALITY_IMAGES:
				filter_str = g_strdup_printf("> %.3lf", param);
				break;
			case BEST_NBSTARS_IMAGES:
				filter_str = g_strdup_printf("> %d", (int)param);
				break;
			}
		}
		gtk_label_set_text(filter_label[filter], filter_str);
		g_free(filter_str);
	}
}

/* Activates or not the stack button if there are 2 or more selected images,
 * all data related to stacking is set in stackparam, except the method itself,
 * determined at stacking start.
 */
void update_stack_interface(gboolean dont_change_stack_type) {
	static GtkWidget *go_stack = NULL, *widgetnormalize = NULL, *force_norm =
			NULL, *output_norm = NULL, *RGB_equal = NULL, *fast_norm = NULL, *max_framing = NULL,
			*upscale_at_stacking = NULL, *blend_frame = NULL, *overlap_norm = NULL;
	static GtkComboBox *method_combo = NULL, *filter_combo = NULL;
	static GtkLabel *result_label = NULL;
	static GtkExpander *stack_expander_method = NULL, *stack_expander_output = NULL;
	gchar *labelbuffer;

	if(!go_stack) {
		go_stack = GTK_WIDGET(gtk_builder_get_object(gui.builder, "gostack_button"));
		filter_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combofilter1"));
		method_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "comboboxstack_methods"));
		widgetnormalize = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combonormalize"));
		force_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkforcenorm"));
		fast_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm"));
		result_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "stackfilter_label"));
		output_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max"));
		RGB_equal = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal"));
		max_framing = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_maximize_framing"));
		upscale_at_stacking = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_upscale_at_stacking"));
		blend_frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_blend_frame"));
		overlap_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_norm_overlap"));
		stack_expander_method = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "stack_expander_method"));
		stack_expander_output = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "stack_expander_output"));
	}
	gboolean seqloaded = sequence_is_loaded();
	gtk_widget_set_sensitive(GTK_WIDGET(stack_expander_method), seqloaded);
	gtk_expander_set_expanded(stack_expander_method, seqloaded);
	gtk_widget_set_sensitive(GTK_WIDGET(stack_expander_output), seqloaded);
	gtk_expander_set_expanded(stack_expander_output, seqloaded);
	/* Refresh the mpp drizzle-mode combo against the new sequence's
	 * channel layout: reset the active selection if it's no longer
	 * valid, and let the cell renderer grey out inapplicable options. */
	mpp_drizzle_combo_validate(GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_mpp_drizzle")));
	if (!seqloaded) {
		return;
	}
	stackparam.seq = &com.seq;

	if (!dont_change_stack_type && stackparam.seq->selnum < stackparam.seq->number) {
		g_signal_handlers_block_by_func(filter_combo, on_stacksel_changed, NULL);
		gtk_combo_box_set_active(filter_combo, SELECTED_IMAGES);
		g_signal_handlers_unblock_by_func(filter_combo, on_stacksel_changed, NULL);
	}

	/* If a .mpp sidecar already exists for this sequence (left by a previous
	 * "Multipoint Registration" run, or by the CLI register_mpp command),
	 * default the stack-method combo to STACK_MPP — the only method that can
	 * consume the sidecar. dont_change_stack_type guards the case where the
	 * user has already manually picked something. */
	if (!dont_change_stack_type && stackparam.seq->seqname) {
		gchar *mpp_path = g_strdup_printf("%s.mpp", stackparam.seq->seqname);
		if (g_file_test(mpp_path, G_FILE_TEST_EXISTS)) {
			gtk_combo_box_set_active(method_combo, STACK_MPP);
		}
		g_free(mpp_path);
	}
	gboolean can_reframe = layer_has_usable_registration(&com.seq, get_registration_layer(&com.seq));
	gboolean can_upscale = can_reframe && !com.seq.is_variable && !com.seq.is_drizzle;
	gboolean must_reframe = can_reframe && com.seq.is_variable;

	int stack_method = gtk_combo_box_get_active(method_combo);
	switch (stack_method) {
	default:
	case STACK_SUM:
	case STACK_MAX:
	case STACK_MIN:
		gtk_widget_set_sensitive(widgetnormalize, FALSE);
		gtk_widget_set_sensitive(force_norm, FALSE);
		gtk_widget_set_sensitive(fast_norm, FALSE);
		gtk_widget_set_visible(output_norm, FALSE);
		gtk_widget_set_visible(RGB_equal, FALSE);
		gtk_widget_set_visible(max_framing, can_reframe); // only shown if applicable
		if (can_reframe) {
			gtk_widget_set_sensitive(max_framing, !must_reframe);
			if (must_reframe)
				gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(max_framing), TRUE);
		}
		gtk_widget_set_visible(upscale_at_stacking, can_upscale); // only shown if applicable
		break;
	case STACK_MEAN:
	case STACK_MEDIAN:
		gtk_widget_set_sensitive(widgetnormalize, TRUE);
		gtk_widget_set_sensitive(force_norm,
				gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
		gtk_widget_set_sensitive(fast_norm,
				gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
		gtk_widget_set_visible(output_norm, TRUE);
		gtk_widget_set_visible(RGB_equal, TRUE);
		gtk_widget_set_visible(max_framing, can_reframe && stack_method != STACK_MEDIAN); // only shown if applicable and not for median
		gtk_widget_set_visible(upscale_at_stacking, can_upscale && stack_method != STACK_MEDIAN); // only shown if applicable and not for median
		if (can_reframe && stack_method != STACK_MEDIAN) {
			gtk_widget_set_sensitive(max_framing, !must_reframe);
			if (must_reframe)
				gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(max_framing), TRUE);
		}
	}

	gboolean can_feather = can_reframe && stack_method == STACK_MEAN;
	gtk_widget_set_visible(blend_frame, can_feather);
	if (can_feather) {
		gtk_widget_set_sensitive(overlap_norm, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(max_framing)));
		if (!gtk_widget_get_sensitive(overlap_norm))
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(overlap_norm), FALSE);
	}

	if (com.seq.reference_image == -1)
		com.seq.reference_image = sequence_find_refimage(&com.seq);
	stackparam.ref_image = com.seq.reference_image;

	get_sequence_filtering_from_gui(
			&stackparam.filtering_criterion, &stackparam.filtering_parameter);

	if (stackparam.description)
		g_free(stackparam.description);
	stackparam.description = describe_filter(stackparam.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);

	update_filter_label();

	stackparam.nb_images_to_stack = compute_nb_filtered_images(&com.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);
	labelbuffer = g_strdup_printf(_("Stacking %d images of the %d of the sequence"),
			stackparam.nb_images_to_stack, com.seq.number);
	gtk_label_set_text(result_label, labelbuffer);
	g_free(labelbuffer);

	if (stackparam.nb_images_to_stack >= 2) {
		stack_fill_list_of_unfiltered_images(&stackparam);
		gtk_widget_set_sensitive(go_stack, TRUE);
	} else {
		gtk_widget_set_sensitive(go_stack, FALSE);
	}
}
