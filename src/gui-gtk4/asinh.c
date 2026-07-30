/*
 * GTK callbacks, preview logic and idle functions for the Asinh Stretch dialog.
 * Pure processing entry points live in src/filters/asinh.c.
 */

#include <math.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/op_descriptors.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/undo.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_replay.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/siril_preview.h"
#include "filters/asinh.h"
#include "gui-gtk4/asinh.h"

static GtkSpinButton *asinh_spin_stretch = NULL, *asinh_spin_black = NULL;
static GtkCheckButton *asinh_toggle_rgb = NULL, *asinh_preview_btn = NULL;
static GtkDropDown *asinh_clipmode_combo = NULL;
static GtkWidget *asinh_clip_settings_widget = NULL;
static GtkWidget *asinh_amend_note = NULL;

/* Amend mode (convergence C4): the dialog edits an existing history
 * record.  gfit holds the pre-record state installed by
 * nde_amend_preview_start; previews run against it through the ordinary
 * backup machinery, and Apply/Cancel route through
 * nde_amend_preview_end instead of a worker run + undo save. */
static gboolean asinh_amend_mode = FALSE;

static void asinh_dialog_init_statics(void) {
	if (asinh_spin_stretch) return;
	asinh_spin_stretch = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_asinh"));
	asinh_spin_black = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "black_point_spin_asinh"));
	asinh_toggle_rgb = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_RGBspace"));
	asinh_clipmode_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "asinh_clipmode"));
	asinh_preview_btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "asinh_preview"));
	asinh_clip_settings_widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, "asinh_clip_settings"));
	asinh_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "asinh_amend_note"));
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
		*rgb_space = siril_toggle_get_active(GTK_WIDGET(asinh_toggle_rgb));
	if (clipmode)
		*clipmode = (clip_mode_t)gtk_drop_down_get_selected(asinh_clipmode_combo);
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
	args->op = &op_desc_asinh;
	args->idle_function = for_preview ? NULL : asinh_apply_idle;
	args->verbose = !for_preview;
	args->user = params;
	args->mask_aware = TRUE;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = gui.roi.active;

	if (for_preview)
		generic_image_worker(args);
	else
		start_in_new_thread(generic_image_worker, args);
	return 0;
}

/* Update preview using the worker */
static int asinh_update_preview() {
	if (siril_toggle_get_active(GTK_WIDGET(asinh_preview_btn))) {
		copy_backup_to_gfit();
		return asinh_process_with_worker(TRUE);
	}
	return 0;
}

void asinh_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(asinh_preview_btn));
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
			/* Remap the Cairo display buffers from the restored gfit.
			 * Without this the display continues to show the preview
			 * because gfit_modified_update_gui() only queues a redraw. */
			notify_gfit_data_modified();
			gfit_modified_update_gui();
		}
	} else {
		invalidate_stats_from_fit(gfit);

		/* Read the full applied parameter set (the same values the OK
		 * handler's worker run used) so the provenance record round-trips
		 * exactly. */
		asinh_params applied = { 0 };
		get_asinh_values(&applied.beta, &applied.offset,
				&applied.human_luminance, &applied.clip_mode);

		gchar *summary = g_strdup_printf(
				_("Asinh Transformation: (stretch=%6.1f, bp=%7.5f)"),
				applied.beta, applied.offset);

		/* on_asinh_ok_clicked ran the worker with skip_generic_undo, so
		 * the worker recorded nothing: this is the sole commit point.
		 * Capture BEFORE the save (worker order) and tag the fresh entry. */
		/* NDE baseline (phase 2): the preview backup holds the pre-stretch pixels. */
		fits *nde_base = get_preview_gfit_backup();
		if (nde_base)
			nde_checkpoint_baseline_ensure(nde_base, nde_checkpoint_active_item_id());
		gint64 rid = nde_capture_from_descriptor(&op_desc_asinh, &applied, summary, gfit, TRUE);

		/* One entry captures the pre-stretch pixels (from the preview
		 * backup) and the pre-stretch ICC profile together, so a single
		 * Ctrl-Z reverts the whole operation. */
		if (!undo_save_state_with_icc(get_preview_gfit_backup(), original_icc,
				"%s", summary))
			undo_tag_top_nde_record(rid);
		g_free(summary);
	}

	roi_supported(FALSE);
	remove_roi_callback(asinh_change_between_roi_and_image);
	if (revert_icc_profile && !single_image_stretch_applied) {
		current_image_set_icc_profile(original_icc
			? copyICCProfile(original_icc) : NULL);
		current_image_color_manage(original_icc != NULL);
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
	if (asinh_amend_mode) {
		/* Leave amend mode: drop the pre-record backup and restore the
		 * true pixels (nothing was committed). */
		cancel_pending_update();
		cancel_and_wait_for_preview();
		clear_backup();
		nde_amend_preview_end(FALSE, NULL);
		asinh_amend_mode = FALSE;
		siril_close_dialog("asinh_dialog");
		return;
	}
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

/* Fill the widgets from the amended record's current parameters through
 * the op's own deserializer — the same struct the normal apply builds. */
static void asinh_prefill_from_amend(void) {
	asinh_params *p = op_desc_asinh.deserialize(nde_amend_preview_params(),
	                                            nde_amend_preview_op_version());
	if (!p)
		return;
	siril_toggle_set_active(GTK_WIDGET(asinh_toggle_rgb), p->human_luminance);
	gtk_spin_button_set_value(asinh_spin_stretch, p->beta);
	gtk_spin_button_set_value(asinh_spin_black, p->offset);
	gtk_drop_down_set_selected(asinh_clipmode_combo, p->clip_mode);
	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
}

void on_asinh_dialog_show(GtkWidget *widget, gpointer user_data) {
	asinh_dialog_init_statics();
	gtk_widget_set_visible(asinh_clip_settings_widget, (gfit->naxes[2] == 3));
	gtk_widget_set_visible(asinh_amend_note, asinh_amend_mode);

	if (gui.rendering_mode == LINEAR_DISPLAY)
		setup_stretch_sliders();

	if (asinh_amend_mode) {
		/* gfit already shows the pre-record state.  No ICC juggling (an
		 * amend only changes pixels) and no ROI (previews are full-image
		 * against the pre-record backup). */
		if (original_icc) {
			cmsCloseProfile(original_icc);
			original_icc = NULL;
		}
		single_image_stretch_applied = FALSE;
		if (gui.roi.active)
			on_clear_roi();
		copy_gfit_to_backup();   /* backup := pre-record state */

		set_notify_block(TRUE);
		siril_toggle_set_active(GTK_WIDGET(asinh_preview_btn), TRUE);
		asinh_prefill_from_amend();
		set_notify_block(FALSE);

		/* One preview tick so the dialog opens showing the record applied
		 * with its current parameters. */
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = asinh_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
		return;
	}

	if (original_icc)
		cmsCloseProfile(original_icc);
	original_icc = copyICCProfile(current_icc_profile());
	icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_STRETCH);
	single_image_stretch_applied = FALSE;

	if (single_image_is_loaded()) {
		if (original_icc) {
			cmsCloseProfile(original_icc);
			original_icc = copyICCProfile(current_icc_profile());
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
	siril_toggle_set_active(GTK_WIDGET(asinh_toggle_rgb), TRUE);
	gtk_spin_button_set_value(asinh_spin_stretch, 0.0f);
	gtk_spin_button_set_value(asinh_spin_black, 0.0f);
	gtk_spin_button_set_increments(asinh_spin_black, 0.001, 0.01);
	set_notify_block(FALSE);
}

void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	if (asinh_amend_mode) {
		/* Serialize the widget state through the SAME params struct and
		 * serializer the normal apply uses, then route it through the
		 * amend path.  preview_end restores the true pixels first (the
		 * on-screen preview is discarded), then replays the edited
		 * history — its tail replay starts from the pre-record state
		 * this preview session deposited in the cache. */
		cancel_pending_update();
		cancel_and_wait_for_preview();
		asinh_params applied = { 0 };
		get_asinh_values(&applied.beta, &applied.offset,
				&applied.human_luminance, &applied.clip_mode);
		gchar *blob = op_desc_asinh.serialize(&applied);
		clear_backup();
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		asinh_amend_mode = FALSE;
		siril_close_dialog("asinh_dialog");
		return;
	}

	if (!check_ok_if_cfa())
		return;

	if (siril_toggle_get_active(GTK_WIDGET(asinh_preview_btn)))
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
	args->op = &op_desc_asinh;
	args->idle_function = asinh_apply_idle;
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->mask_aware = TRUE;
	args->for_preview = FALSE;
	args->for_roi = FALSE;
	args->skip_generic_undo = TRUE;

	start_in_new_thread(generic_image_worker, args);

	apply_asinh_changes();
	siril_close_dialog("asinh_dialog");
}

gboolean on_asinh_dialog_close(GtkWindow *dialog, gpointer user_data) {
	(void)dialog; (void)user_data;
	apply_asinh_changes();
	return FALSE;
}

void on_asinh_undo_clicked(GtkButton *button, gpointer user_data) {
	set_notify_block(TRUE);
	if (asinh_amend_mode) {
		/* Reset means "back to the record's recorded parameters" here. */
		asinh_prefill_from_amend();
	} else {
		siril_toggle_set_active(GTK_WIDGET(asinh_toggle_rgb), TRUE);
		gtk_spin_button_set_value(asinh_spin_black, 0);
		gtk_spin_button_set_value(asinh_spin_stretch, 0);
		gtk_drop_down_set_selected(asinh_clipmode_combo, 2);
	}
	siril_toggle_set_active(GTK_WIDGET(asinh_preview_btn), TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(asinh_preview_btn));
	notify_update((gpointer) param);
}

void on_asinh_parameter_changed(GtkWidget *widget, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(asinh_preview_btn));
	notify_update((gpointer) param);
}

void on_asinh_preview_toggled(GtkCheckButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!siril_toggle_get_active(GTK_WIDGET(button))) {
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

void on_asinh_clipmode_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	on_asinh_parameter_changed(GTK_WIDGET(combo), user_data);
}

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void asinh_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	asinh_amend_mode = TRUE;
	siril_open_dialog("asinh_dialog");
}

gboolean asinh_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed.  A synchronous refusal is
	 * logged by the core and simply leaves everything unchanged. */
	nde_amend_preview_start(record_id, asinh_amend_ready, NULL);
	return TRUE;
}
