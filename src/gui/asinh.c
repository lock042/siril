/*
 * GTK callbacks, preview logic and idle functions for the Asinh Stretch dialog.
 * Pure processing entry points live in src/filters/asinh.c.
 */

#include <math.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/undo.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "filters/asinh.h"
#include "gui/asinh.h"

static GtkSpinButton *asinh_spin_stretch = NULL, *asinh_spin_black = NULL;
static GtkToggleButton *asinh_toggle_rgb = NULL, *asinh_preview_btn = NULL;
static GtkComboBox *asinh_clipmode_combo = NULL;
static GtkWidget *asinh_clip_settings_widget = NULL;

static void asinh_dialog_init_statics(void) {
	if (asinh_spin_stretch) return;
	asinh_spin_stretch = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_asinh"));
	asinh_spin_black = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "black_point_spin_asinh"));
	asinh_toggle_rgb = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_RGBspace"));
	asinh_clipmode_combo = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "asinh_clipmode"));
	asinh_preview_btn = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "asinh_preview"));
	asinh_clip_settings_widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, "asinh_clip_settings"));
}

// Original ICC profile, in case we don't apply a stretch and need to revert
static cmsHPROFILE original_icc = NULL;
static gboolean single_image_stretch_applied = FALSE;

/* Idle function for final application */
static gboolean asinh_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		single_image_stretch_applied = TRUE;
		gfit_modified_update_gui();
	}
	free_generic_img_args(args);
	return FALSE;
}

/* Helper function to get current widget values */
static void get_asinh_values(float *stretch, float *black, gboolean *rgb_space, clip_mode_t *clipmode) {
	if (stretch)
		*stretch = gtk_spin_button_get_value(asinh_spin_stretch);
	if (black)
		*black = gtk_spin_button_get_value(asinh_spin_black);
	if (rgb_space)
		*rgb_space = gtk_toggle_button_get_active(asinh_toggle_rgb);
	if (clipmode)
		*clipmode = (clip_mode_t)gtk_combo_box_get_active(asinh_clipmode_combo);
}

/* Create and launch asinh processing */
static int asinh_process_with_worker(gboolean for_preview) {
	asinh_params *params = calloc(1, sizeof(asinh_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	get_asinh_values(&params->beta, &params->offset, &params->human_luminance, &params->clip_mode);

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return 1;
	}

	args->fit = gui.roi.active ? &gui.roi.fit : gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = asinh_image_hook;
	args->idle_function = for_preview ? NULL : asinh_apply_idle;
	args->description = _("Asinh stretch");
	args->verbose = !for_preview;
	args->user = params;
	args->mask_aware = TRUE;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = gui.roi.active;
	if (!for_preview)
		args->populate_roi_on_complete = TRUE;

	if (for_preview)
		generic_image_worker(args);
	else
		start_in_new_thread(generic_image_worker, args);
	return 0;
}

/* Update preview using the worker */
static int asinh_update_preview() {
	if (gtk_toggle_button_get_active(asinh_preview_btn)) {
		copy_backup_to_gfit();
		return asinh_process_with_worker(TRUE);
	}
	return 0;
}

void asinh_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(asinh_preview_btn);
	notify_update((gpointer) param);
}

static void asinh_startup() {
	add_roi_callback(asinh_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();
}

static void asinh_close(gboolean revert, gboolean revert_icc_profile) {
	set_cursor_waiting(TRUE);

	if (revert) {
		float stretch_value, black_value;
		get_asinh_values(&stretch_value, &black_value, NULL, NULL);

		if (stretch_value != 0.0f || black_value != 0.0f) {
			copy_backup_to_gfit();
			notify_gfit_data_modified();
			gfit_modified_update_gui();
		}
	} else {
		invalidate_stats_from_fit(gfit);

		float stretch_value, black_value;
		get_asinh_values(&stretch_value, &black_value, NULL, NULL);

		fits undo_fit = {0};
		memcpy(&undo_fit, get_preview_gfit_backup(), sizeof(fits));
		undo_fit.icc_profile = original_icc;
		undo_fit.color_managed = original_icc != NULL;

		undo_save_state(&undo_fit,
				_("Asinh Transformation: (stretch=%6.1lf, bp=%7.5lf)"),
				stretch_value, black_value);
	}

	roi_supported(FALSE);
	remove_roi_callback(asinh_change_between_roi_and_image);
	if (revert_icc_profile && !single_image_stretch_applied) {
		if (gfit->icc_profile)
			cmsCloseProfile(gfit->icc_profile);
		gfit->icc_profile = copyICCProfile(original_icc);
		color_manage(gfit, gfit->icc_profile != NULL);
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

/* Keep the existing asinhlut functions with clip_mode as parameter */
static void apply_asinh_changes() {
	float stretch_value, black_value;
	gboolean rgb_space;
	get_asinh_values(&stretch_value, &black_value, &rgb_space, NULL);

	gboolean changed = (stretch_value != 1.0f) || (black_value != 0.0f) || !rgb_space;
	asinh_close(!changed, !changed);
}

void apply_asinh_cancel() {
	asinh_close(TRUE, TRUE);
	siril_close_dialog("asinh_dialog");
}

void on_asinh_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_asinh_cancel();
}

gboolean asinh_hide_on_delete(GtkWidget *widget) {
	apply_asinh_cancel();
	return TRUE;
}

void on_asinh_dialog_show(GtkWidget *widget, gpointer user_data) {
	asinh_dialog_init_statics();
	gtk_widget_set_visible(asinh_clip_settings_widget, (gfit->naxes[2] == 3));

	if (gui.rendering_mode == LINEAR_DISPLAY)
		setup_stretch_sliders();

	if (original_icc)
		cmsCloseProfile(original_icc);
	original_icc = copyICCProfile(gfit->icc_profile);
	icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_STRETCH);
	single_image_stretch_applied = FALSE;

	if (single_image_is_loaded()) {
		if (original_icc) {
			cmsCloseProfile(original_icc);
			original_icc = copyICCProfile(gfit->icc_profile);
		}
		icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_STRETCH);
	} else {
		if (original_icc) {
			cmsCloseProfile(original_icc);
			original_icc = NULL;
		}
	}

	asinh_startup();

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(asinh_toggle_rgb, TRUE);
	gtk_spin_button_set_value(asinh_spin_stretch, 0.0f);
	gtk_spin_button_set_value(asinh_spin_black, 0.0f);
	gtk_spin_button_set_increments(asinh_spin_black, 0.001, 0.01);
	set_notify_block(FALSE);
}

void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	if (gtk_toggle_button_get_active(asinh_preview_btn))
		copy_backup_to_gfit();

	remove_roi_callback(asinh_change_between_roi_and_image);

	asinh_params *params = calloc(1, sizeof(asinh_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		add_roi_callback(asinh_change_between_roi_and_image);
		return;
	}

	get_asinh_values(&params->beta, &params->offset, &params->human_luminance, &params->clip_mode);

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		add_roi_callback(asinh_change_between_roi_and_image);
		return;
	}

	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = asinh_image_hook;
	args->idle_function = asinh_apply_idle;
	args->description = _("Asinh stretch");
	args->verbose = TRUE;
	args->user = params;
	args->log_hook = asinh_log_hook;
	args->max_threads = com.max_thread;
	args->mask_aware = TRUE;
	args->for_preview = FALSE;
	args->for_roi = FALSE;
	args->custom_undo = TRUE;
	args->populate_roi_on_complete = TRUE;

	start_in_new_thread(generic_image_worker, args);

	apply_asinh_changes();
	siril_close_dialog("asinh_dialog");
}

gboolean on_asinh_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_asinh_changes();
	return FALSE;
}

void on_asinh_undo_clicked(GtkButton *button, gpointer user_data) {
	set_notify_block(TRUE);
	gtk_toggle_button_set_active(asinh_toggle_rgb, TRUE);
	gtk_spin_button_set_value(asinh_spin_black, 0);
	gtk_spin_button_set_value(asinh_spin_stretch, 0);
	gtk_combo_box_set_active(asinh_clipmode_combo, 2);
	gtk_toggle_button_set_active(asinh_preview_btn, TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(asinh_preview_btn);
	notify_update((gpointer) param);
}

void on_asinh_parameter_changed(GtkWidget *widget, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(asinh_preview_btn);
	notify_update((gpointer) param);
}

void on_asinh_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(button)) {
		cancel_and_wait_for_preview();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = asinh_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}

void on_asinh_clipmode_changed(GtkComboBox *combo, gpointer user_data) {
	on_asinh_parameter_changed(GTK_WIDGET(combo), user_data);
}
