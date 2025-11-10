/*
 * Refactored asinh stretch using generic_image_worker
 */

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "algos/statistics.h"
#include "core/arithm.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"

#include "asinh.h"

// Original ICC profile, in case we don't apply a stretch and need to revert
static cmsHPROFILE original_icc = NULL;
static gboolean single_image_stretch_applied = FALSE;

/* The actual asinh processing hook */
int asinh_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	asinh_params *params = (asinh_params *)args->user;
	if (!params)
		return 1;

	return asinhlut(fit, params->beta, params->offset, params->human_luminance, params->clip_mode);
}

/* Idle function for preview updates */
static gboolean asinh_preview_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	if (args->retval == 0) {
		notify_gfit_modified();
	}

	// Free the params
	if (args->user)
		free(args->user);

	stop_processing_thread();
	free(args);
	return FALSE;
}

/* Idle function for final application */
static gboolean asinh_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	asinh_params *params = (asinh_params *)args->user;

	if (args->retval == 0) {
		single_image_stretch_applied = TRUE;
		populate_roi();
		notify_gfit_modified();
	}

	// Free the params
	if (params)
		free(params);

	stop_processing_thread();
	free(args);
	return FALSE;
}

/* Helper function to get current widget values */
static void get_asinh_values(float *stretch, float *black, gboolean *rgb_space, clip_mode_t *clipmode) {
	if (stretch)
		*stretch = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_asinh")));
	if (black)
		*black = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("black_point_spin_asinh")));
	if (rgb_space)
		*rgb_space = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace")));
	if (clipmode)
		*clipmode = (clip_mode_t)gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("asinh_clipmode")));
}

/* Create and launch asinh processing */
static int asinh_process_with_worker(gboolean for_preview) {
	// Allocate parameters
	asinh_params *params = calloc(1, sizeof(asinh_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Get current values from widgets
	get_asinh_values(&params->beta, &params->offset, &params->human_luminance, &params->clip_mode);

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		return 1;
	}

	// Set the fit based on whether ROI is active
	args->fit = gui.roi.active ? &gui.roi.fit : gfit;
	args->mem_ratio = 1.0f; // asinh needs roughly 1x image size for working memory
	args->image_hook = asinh_image_hook;
	args->idle_function = for_preview ? asinh_preview_idle : asinh_apply_idle;
	args->description = _("Asinh stretch");
	args->verbose = !for_preview; // Only verbose for final application
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = gui.roi.active;

	start_in_new_thread(generic_image_worker, args);
	return 0;
}

/* Update preview using the worker */
static int asinh_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")))) {
		copy_backup_to_gfit();
		return asinh_process_with_worker(TRUE);
	}
	return 0;
}

void asinh_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
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
			notify_gfit_modified();
		}
	} else {
		// Save undo state when applying (not reverting)
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
int asinhlut_ushort(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clip_mode) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	const gboolean do_channel[3] = { TRUE, TRUE, TRUE };
	float m_CB = 1.f;
	float norm = get_normalized_value(fit);
	float invnorm = 1.0f / norm;
	float asinh_beta = asinh(beta);
	float factor_red = human_luminance ? 0.2126f : 0.3333f;
	float factor_green = human_luminance ? 0.7152f : 0.3333f;
	float factor_blue = human_luminance ? 0.0722f : 0.3333f;

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	float globalmax = -FLT_MAX;
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->rx * 16)
#endif
		for (i = 0; i < n; i++) {
			blend_data data = { .sf = { 0.f }, .tf = { 0.f }, .do_channel = do_channel };
			float val[3] = { buf[RLAYER][i] * invnorm, buf[GLAYER][i] * invnorm, buf[BLAYER][i] * invnorm };
			float prime[3];
			for (int chan = 0 ; chan < 3 ; chan++)
				prime[chan] = max(0, (val[chan] - offset) / (1.0f - offset));

			float x = factor_red * prime[RLAYER] + factor_green * prime[GLAYER] + factor_blue * prime[BLAYER];
			float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * x) / (x * asinh_beta);
			for (int chan = 0 ; chan < 3 ; chan++)
				data.sf[chan] = min(1.0f, max(0.0f, (prime[chan] * k)));
			if (clip_mode == RGBBLEND) {
				for (int chan = 0 ; chan < 3 ; chan++)
					data.tf[chan] = (prime[chan] == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * prime[chan]) / (prime[chan] * asinh_beta);
			}
			float maxval = -FLT_MAX;
			switch (clip_mode) {
				case CLIP:
					for (int chan = 0 ; chan < 3 ; chan++)
						buf[chan][i] = do_channel[chan] ? roundf_to_WORD(norm * fmaxf(0.f, fminf(data.sf[chan], 1.f))) : roundf_to_WORD(fmaxf(0.f, val[chan]));
					break;
				case RESCALE:
					for (int chan = 0 ; chan < 3 ; chan++) {
						if (do_channel[chan] && data.sf[chan] > maxval)
							maxval = data.sf[chan];
					}
					if (maxval > 1.f) {
						float invmaxval = 1.f / maxval;
						for (int chan = 0 ; chan < 3 ; chan++) {
							data.sf[chan] *= invmaxval;
							buf[chan][i] = do_channel[chan] ? roundf_to_WORD(norm * fmaxf(0.f, data.sf[chan])) : roundf_to_WORD(norm * fmaxf(0.f, val[chan]));
						}
					} else {
						for (int chan = 0 ; chan < 3 ; chan++) {
							buf[chan][i] = do_channel[chan] ? roundf_to_WORD(norm * fmaxf(0.f, data.sf[chan])) : roundf_to_WORD(norm * fmaxf(0.f, val[chan]));
						}
					}
					break;
				case RGBBLEND:;
					float out[3] = {val[0], val[1], val[2]};
					rgbblend(&data, &out[RLAYER], &out[GLAYER], &out[BLAYER], m_CB);
					for (int chan = 0 ; chan < 3 ; chan++) {
						buf[chan][i] = roundf_to_WORD(norm * out[chan]);
					}
					break;
				case RESCALEGLOBAL:
					for (int chan = 0 ; chan < 3 ; chan++) {
						if (do_channel[chan] && data.sf[chan] > maxval)
							maxval = data.sf[chan];
					}
					if (maxval > globalmax)
						globalmax = maxval;
					break;
			}
		}
		if (clip_mode == RESCALEGLOBAL) {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0 ; i < n ; i++) {
				blend_data data = { .sf = { 0.f }, .tf = { 0.f }, .do_channel = do_channel };
				float val[3] = { buf[RLAYER][i] * invnorm, buf[GLAYER][i] * invnorm, buf[BLAYER][i] * invnorm };
				float prime[3];
				for (int chan = 0 ; chan < 3 ; chan++)
					prime[chan] = max(0, (val[chan] - offset) / (1.0f - offset));

				float x = factor_red * prime[RLAYER] + factor_green * prime[GLAYER] + factor_blue * prime[BLAYER];
				float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * x) / (x * asinh_beta);
				for (int chan = 0 ; chan < 3 ; chan++) {
					data.sf[chan] = min(1.0f, max(0.0f, (prime[chan] * k))) / globalmax;
					buf[chan][i] = (do_channel[chan]) ? roundf_to_WORD(norm * fmaxf(0.f, data.sf[chan])) : roundf_to_WORD(fmaxf(0.f, val[chan]));
				}
			}
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->rx * 16)
#endif
		for (i = 0; i < n; i++) {
			float x = buf[RLAYER][i] * invnorm;
			float xprime = max(0, (x - offset) / (1.0f - offset));
			float k = (xprime == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f,(xprime) * k)));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int asinhlut_float(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clip_mode) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	const gboolean do_channel[3] = { TRUE, TRUE, TRUE };

	float m_CB = 1.f;
	float asinh_beta = asinhf(beta);
	float factor_red = human_luminance ? 0.2126f : 0.3333f;
	float factor_green = human_luminance ? 0.7152f : 0.3333f;
	float factor_blue = human_luminance ? 0.0722f : 0.3333f;

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	float globalmax = -FLT_MAX;
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (i = 0; i < n; i++) {
			blend_data data = { .sf = { 0.f }, .tf = { 0.f }, .do_channel = do_channel };
			float val[3] = { buf[RLAYER][i], buf[GLAYER][i], buf[BLAYER][i] };
			float prime[3];
			for (int chan = 0 ; chan < 3 ; chan++)
				prime[chan] = max(0, (val[chan] - offset) / (1.0f - offset));

			float x = factor_red * prime[RLAYER] + factor_green * prime[GLAYER] + factor_blue * prime[BLAYER];
			float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * x) / (x * asinh_beta);
			for (int chan = 0 ; chan < 3 ; chan++)
				data.sf[chan] = min(1.0f, max(0.0f, (prime[chan] * k)));
			if (clip_mode == RGBBLEND) {
				for (int chan = 0 ; chan < 3 ; chan++)
					data.tf[chan] = (prime[chan] == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * prime[chan]) / (prime[chan] * asinh_beta);
			}
			float maxval = -FLT_MAX;
			switch (clip_mode) {
				case CLIP:
					for (int chan = 0 ; chan < 3 ; chan++)
						buf[chan][i] = do_channel[chan] ? fmaxf(0.f, fminf(data.sf[chan], 1.f)) : fmaxf(0.f, val[chan]);
					break;
				case RESCALE:
					for (int chan = 0 ; chan < 3 ; chan++) {
						if (do_channel[chan] && data.sf[chan] > maxval)
							maxval = data.sf[chan];
					}
					if (maxval > 1.f) {
						for (int chan = 0 ; chan < 3 ; chan++) {
							data.sf[chan] /= maxval;
							buf[chan][i] = do_channel[chan] ? fmaxf(0.f, data.sf[chan]) : fmaxf(0.f, val[chan]);
						}
					}
					break;
				case RGBBLEND:
					rgbblend(&data, &buf[RLAYER][i], &buf[GLAYER][i], &buf[BLAYER][i], m_CB);
					break;
				case RESCALEGLOBAL:
					for (int chan = 0 ; chan < 3 ; chan++) {
						if (do_channel[chan] && data.sf[chan] > maxval)
							maxval = data.sf[chan];
					}
					if (maxval > globalmax)
						globalmax = maxval;
					break;
			}
		}
		if (clip_mode == RESCALEGLOBAL) {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0 ; i < n ; i++) {
				blend_data data = { .sf = { 0.f }, .tf = { 0.f }, .do_channel = do_channel };
				float val[3] = { buf[RLAYER][i], buf[GLAYER][i], buf[BLAYER][i] };
				float prime[3];
				for (int chan = 0 ; chan < 3 ; chan++)
					prime[chan] = max(0, (val[chan] - offset) / (1.0f - offset));

				float x = factor_red * prime[RLAYER] + factor_green * prime[GLAYER] + factor_blue * prime[BLAYER];
				float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * x) / (x * asinh_beta);
				for (int chan = 0 ; chan < 3 ; chan++) {
					data.sf[chan] = min(1.0f, max(0.0f, (prime[chan] * k))) / globalmax;
					buf[chan][i] = (do_channel[chan]) ? fmaxf(0.f, data.sf[chan]) : fmaxf(0.f, val[chan]);
				}
			}
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (i = 0; i < n; i++) {
			float x, k;
			x = buf[RLAYER][i];
			float xprime = max(0.0f, (x - offset) / (1.0f - offset));
			k = (xprime == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = min(1.0f, max(0.0f,(xprime * k)));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int asinhlut(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clip_mode) {
	if (fit->type == DATA_USHORT)
		return asinhlut_ushort(fit, beta, offset, human_luminance, clip_mode);
	if (fit->type == DATA_FLOAT)
		return asinhlut_float(fit, beta, offset, human_luminance, clip_mode);
	return 1;
}

int command_asinh(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clipmode) {
	return asinhlut(fit, beta, offset, human_luminance, clipmode);
}

/* Keep all the GTK callbacks - now they don't set static variables */
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
	GtkSpinButton *spin_stretch = GTK_SPIN_BUTTON(lookup_widget("spin_asinh"));
	GtkSpinButton *spin_black_p = GTK_SPIN_BUTTON(lookup_widget("black_point_spin_asinh"));
	GtkToggleButton *toggle_rgb = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace"));
	GtkWidget *clipmode = lookup_widget("asinh_clip_settings");
	gtk_widget_set_visible(clipmode, (gfit->naxes[2] == 3));

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
	gtk_toggle_button_set_active(toggle_rgb, TRUE);
	gtk_spin_button_set_value(spin_stretch, 0.0f);
	gtk_spin_button_set_value(spin_black_p, 0.0f);
	gtk_spin_button_set_increments(spin_black_p, 0.001, 0.01);
	set_notify_block(FALSE);
}

void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	// Always process full image when Apply is clicked
	// even if ROI is active or preview is on
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")))) {
		// Preview is on, need to copy backup to gfit first
		copy_backup_to_gfit();
	}

	// Temporarily disable ROI callback to avoid interference
	remove_roi_callback(asinh_change_between_roi_and_image);

	// Allocate parameters for full image processing
	asinh_params *params = calloc(1, sizeof(asinh_params));
	if (!params) {
		PRINT_ALLOC_ERR;
		add_roi_callback(asinh_change_between_roi_and_image);
		return;
	}

	// Get current values from widgets
	get_asinh_values(&params->beta, &params->offset, &params->human_luminance, &params->clip_mode);

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free(params);
		add_roi_callback(asinh_change_between_roi_and_image);
		return;
	}

	// Always process the full gfit, not ROI
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = asinh_image_hook;
	args->idle_function = asinh_apply_idle;
	args->description = _("Asinh stretch");
	args->verbose = TRUE;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	start_in_new_thread(generic_image_worker, args);

	// Dialog will be closed and cleanup will happen in apply_asinh_changes
	apply_asinh_changes();
	siril_close_dialog("asinh_dialog");
}

gboolean on_asinh_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_asinh_changes();
	return FALSE;
}

void on_asinh_undo_clicked(GtkButton *button, gpointer user_data) {
	set_notify_block(TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace")), TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("black_point_spin_asinh")), 0);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_asinh")), 0);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("asinh_clipmode")), 2);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")), TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
	notify_update((gpointer) param);
}

/* Generic callback for all parameter changes that trigger preview update */
void on_asinh_parameter_changed(GtkWidget *widget, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
	notify_update((gpointer) param);
}

void on_asinh_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(button)) {
		waiting_for_thread();
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
