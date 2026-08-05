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

#include <assert.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/colors.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/op_descriptors.h"
#include "core/OS_utils.h"
#include "io/image_format_flis.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "algos/statistics.h"

#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/histogram.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/colors.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/flis_gui.h"
#include "core/nde_replay.h"

/* Amend mode (convergence C5b): the CCM dialog edits a color.ccm history
 * record.  gfit holds the pre-record state installed by
 * nde_amend_preview_start.  CCM has no preview machinery (dialogs.c
 * has_preview=FALSE), so there is no backup to arm or clear: the display shows
 * the installed pre-record state while the form is open, and Apply/Cancel route
 * through nde_amend_preview_end.  ICC is NOT touched in amend mode (skip-ICC
 * rule): converting the installed pre-record pixels would make the preview
 * diverge from the replayed commit. */
static gboolean ccm_amend_mode = FALSE;

void on_button_bkg_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_h"));
	}

	gtk_spin_button_set_value(selection_black_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_black_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_black_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_black_value[3], com.selection.h);
}

void initialize_calibration_interface() {
	static GtkAdjustment *selection_black_adjustment[4] = { NULL, NULL, NULL,
		NULL };
	static GtkAdjustment *selection_white_adjustment[4] = { NULL, NULL, NULL,
		NULL };

	if (!selection_black_adjustment[0]) {
		selection_black_adjustment[0] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_x"));
		selection_black_adjustment[1] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_y"));
		selection_black_adjustment[2] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_w"));
		selection_black_adjustment[3] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_h"));
	}
	if (!selection_white_adjustment[0]) {
		selection_white_adjustment[0] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_x"));
		selection_white_adjustment[1] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_y"));
		selection_white_adjustment[2] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_w"));
		selection_white_adjustment[3] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_h"));
	}
	gtk_adjustment_set_upper(selection_black_adjustment[0], gfit->rx);
	gtk_adjustment_set_upper(selection_black_adjustment[1], gfit->ry);
	gtk_adjustment_set_upper(selection_black_adjustment[2], gfit->rx);
	gtk_adjustment_set_upper(selection_black_adjustment[3], gfit->ry);
	gtk_adjustment_set_value(selection_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_black_adjustment[3], 0);

	gtk_adjustment_set_upper(selection_white_adjustment[0], gfit->rx);
	gtk_adjustment_set_upper(selection_white_adjustment[1], gfit->ry);
	gtk_adjustment_set_upper(selection_white_adjustment[2], gfit->rx);
	gtk_adjustment_set_upper(selection_white_adjustment[3], gfit->ry);
	gtk_adjustment_set_value(selection_white_adjustment[0], 0);
	gtk_adjustment_set_value(selection_white_adjustment[1], 0);
	gtk_adjustment_set_value(selection_white_adjustment[2], 0);
	gtk_adjustment_set_value(selection_white_adjustment[3], 0);
}


/* What the direct-apply path used to do inline after mutating gfit.  The
 * worker already handles undo, the NDE record, populate_roi() and
 * notify_gfit_data_modified(); what is left is the dialog's own tidying. */
static gboolean end_colour_tool(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *) p;
	stop_processing_thread();
	gfit_modified_update_gui();
	if (!args->retval)
		delete_selected_area();
	update_gfit_histogram_if_needed();
	redraw(REDRAW_ALL);
	gui_function(redraw_previews, NULL);
	free_generic_img_args(args);
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_button_bkg_neutralization_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	rectangle black_selection;
	int width, height;

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_h"));
	}
	width = (int) gtk_spin_button_get_value(selection_black_value[2]);
	height = (int) gtk_spin_button_get_value(selection_black_value[3]);

	if ((!width) || (!height)) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}
	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);

	/* Through the generic worker, like every other pixel op.  This one used to
	 * apply in place here and hand-roll its own undo, baseline and capture —
	 * an oversight, and the reason its provenance was opaque: without a
	 * descriptor there was nowhere to hang a serializer. */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	struct bkg_neutral_data *p = new_bkg_neutral_data();
	if (!p) {
		free(args);
		PRINT_ALLOC_ERR;
		return;
	}
	p->black_selection = black_selection;
	args->fit = gfit;
	args->op = &op_desc_bkg_neutral;
	args->user = p;
	args->verbose = TRUE;
	args->idle_function = end_colour_tool;
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

void on_button_white_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_white_value[4] = { NULL, NULL, NULL, NULL };

	if (!selection_white_value[0]) {
		selection_white_value[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_x"));
		selection_white_value[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_y"));
		selection_white_value[2] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_w"));
		selection_white_value[3] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_h"));
	}

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the white reference area"));
		return;
	}

	gtk_spin_button_set_value(selection_white_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_white_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_white_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_white_value[3], com.selection.h);
}

/* Read the dialog into the op's params.  The pixel work moved to
 * color_calib_image_hook (algos/colors.c) — it could not be serialized while it
 * lived here reading GtkRanges mid-calculation. */
static void fill_calib_params_from_widgets(struct color_calib_data *p) {
	static GtkRange *scale_white_balance[3] = { NULL, NULL, NULL };
	static GtkRange *scaleLimit[2] = { NULL, NULL };

	if (scale_white_balance[RLAYER] == NULL) {
		scale_white_balance[RLAYER] = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_r"));
		scale_white_balance[GLAYER] = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_g"));
		scale_white_balance[BLAYER] = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_b"));

		scaleLimit[0] = GTK_RANGE(gtk_builder_get_object(gui.builder, "lowWhiteColorCalibScale"));
		scaleLimit[1] = GTK_RANGE(gtk_builder_get_object(gui.builder, "upWhiteColorCalibScale"));
	}

	if (p->is_manual) {
		p->kw[RLAYER] = gtk_range_get_value(scale_white_balance[RLAYER]);
		p->kw[GLAYER] = gtk_range_get_value(scale_white_balance[GLAYER]);
		p->kw[BLAYER] = gtk_range_get_value(scale_white_balance[BLAYER]);
		/* Manual coefficients are the effective ones already; the hook has
		 * nothing left to work out. */
		p->have_effective = TRUE;
	} else {
		p->low = gtk_range_get_value(scaleLimit[0]);
		p->high = gtk_range_get_value(scaleLimit[1]);
	}
}

void on_calibration_apply_button_clicked(GtkButton *button, gpointer user_data) {
	rectangle black_selection, white_selection;
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	static GtkSpinButton *selection_white_value[4] = { NULL, NULL, NULL, NULL };

	/* The worker times and logs the operation itself (args->verbose). */
	static GtkCheckButton *manual = NULL;
	if (!manual) manual = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_manual_calibration"));
	gboolean is_manual = siril_toggle_get_active(GTK_WIDGET(manual));

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_bkg_h"));
	}

	if (!selection_white_value[0]) {
		selection_white_value[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_x"));
		selection_white_value[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_y"));
		selection_white_value[2] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_w"));
		selection_white_value[3] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_white_h"));
	}

	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);

	if ((!black_selection.w || !black_selection.h) && !is_manual) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	white_selection.x = gtk_spin_button_get_value(selection_white_value[0]);
	white_selection.y = gtk_spin_button_get_value(selection_white_value[1]);
	white_selection.w = gtk_spin_button_get_value(selection_white_value[2]);
	white_selection.h = gtk_spin_button_get_value(selection_white_value[3]);

	if ((!white_selection.w || !white_selection.h) && !is_manual) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the white reference area"));
		return;
	}

	/* Group path: a selected layer group that composites to colour gets the
	 * calibration applied to its composite and distributed into per-layer
	 * factors (flis_group_calibration_hook), instead of the generic worker
	 * refusing the group selection. */
	gint gid = is_current_image_flis() ? flis_panel_selected_group_id() : 0;
	if (gid) {
		struct color_calib_data widget_vals = { 0 };
		widget_vals.is_manual = is_manual;
		fill_calib_params_from_widgets(&widget_vals);

		flis_group_t *grp = flis_group_get_by_id(gid);
		GSList *members = grp ? flis_group_get_layers(grp) : NULL;
		if (!members) {
			siril_log_error(_("Colour calibration: the selected group has no layers\n"));
			return;
		}
		if (undo_save_flis_multi_layer(members, _("Colour calibration (group)")))
			siril_log_warning(_("Colour calibration: could not save undo state\n"));
		g_slist_free(members);

		struct flis_group_calib_args *p = new_flis_group_calib_args();
		if (!p) {
			PRINT_ALLOC_ERR;
			return;
		}
		p->kind = FLIS_GROUP_CALIB_CC;
		p->group_id = gid;
		p->manual = is_manual;
		p->manual_kw[0] = widget_vals.kw[0];
		p->manual_kw[1] = widget_vals.kw[1];
		p->manual_kw[2] = widget_vals.kw[2];
		p->white_sel = white_selection;
		p->black_sel = black_selection;
		p->low = widget_vals.low;
		p->high = widget_vals.high;

		struct generic_layer_args *gargs = calloc(1, sizeof(*gargs));
		if (!gargs) {
			free_flis_group_calib_args(p);
			PRINT_ALLOC_ERR;
			return;
		}
		gargs->layer = NULL;	/* multi-layer op */
		gargs->layer_hook = flis_group_calibration_hook;
		gargs->user = p;
		gargs->description = g_strdup(_("Colour calibration (group)"));
		gargs->invalidate_flags = FLIS_INV_ALL;
		set_cursor_waiting(TRUE);
		if (!start_in_new_thread(generic_layer_worker, gargs))
			free_generic_layer_args(gargs);
		set_cursor_waiting(FALSE);
		return;
	}

	/* Through the generic worker, like every other pixel op — see the
	 * background-neutralization button above. */
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	struct color_calib_data *p = new_color_calib_data();
	if (!p) {
		free(args);
		PRINT_ALLOC_ERR;
		return;
	}
	p->is_manual = is_manual;
	p->white_selection = white_selection;
	p->black_selection = black_selection;
	fill_calib_params_from_widgets(p);
	args->fit = gfit;
	args->op = &op_desc_color_calib;
	args->user = p;
	args->verbose = TRUE;
	args->idle_function = end_colour_tool;
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}

void on_calibration_close_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("color_calibration");
}

gboolean calibration_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("color_calibration");
	return TRUE;
}

void on_checkbutton_manual_calibration_toggled(GtkCheckButton *togglebutton,
		gpointer user_data) {
	static GtkWidget *cc_box_red = NULL, *scale_r = NULL, *cc_box_green = NULL;
	static GtkWidget *scale_g = NULL, *cc_box_blue = NULL, *scale_b = NULL;
	if (!cc_box_red) {
		cc_box_red = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cc_box_red"));
		scale_r = GTK_WIDGET(gtk_builder_get_object(gui.builder, "scale_r"));
		cc_box_green = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cc_box_green"));
		scale_g = GTK_WIDGET(gtk_builder_get_object(gui.builder, "scale_g"));
		cc_box_blue = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cc_box_blue"));
		scale_b = GTK_WIDGET(gtk_builder_get_object(gui.builder, "scale_b"));
	}
	gtk_widget_set_sensitive(cc_box_red, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
	gtk_widget_set_sensitive(scale_r, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
	gtk_widget_set_sensitive(cc_box_green, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
	gtk_widget_set_sensitive(scale_g, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
	gtk_widget_set_sensitive(cc_box_blue, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
	gtk_widget_set_sensitive(scale_b, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
}

void negative_processing() {
	set_cursor_waiting(TRUE);
	/* Direct-apply of the same op as the neg command's hook (arith.neg): no
	 * worker runs here, so this is the sole commit point.  pos_to_neg is void
	 * — capture after it and tag the entry saved just above.  The op is
	 * serializer-less/paramless, so params NULL yields a Tier B record. */
	int undo_err = undo_save_state(gfit, _("Negative Transformation"));
	siril_log_info(_("Negative Transformation\n"));
	/* NDE baseline (phase 2): direct-apply — snapshot pre-op pixels first. */
	nde_checkpoint_baseline_ensure(gfit, nde_capture_target_item());
	pos_to_neg(gfit);
	gint64 rid = nde_capture_from_descriptor(&op_desc_neg, NULL,
			_("Negative Transformation"), gfit, FALSE);
	if (!undo_err)
		undo_tag_top_nde_record(rid);
	invalidate_stats_from_fit(gfit);
	invalidate_gfit_histogram();
	update_gfit_histogram_if_needed();
	notify_gfit_data_modified();
	redraw(REDRAW_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}
/**********************************************************************/

void on_extract_channel_button_close_clicked(GtkButton *button,
		gpointer user_data) {
	siril_close_dialog("extract_channel_dialog");
}

void on_combo_extract_colors_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	static GtkLabel *label_c1 = NULL, *label_c2 = NULL, *label_c3 = NULL;
	if (!label_c1) {
		label_c1 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_extract_c1"));
		label_c2 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_extract_c2"));
		label_c3 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_extract_c3"));
	}
	switch(gtk_drop_down_get_selected(box)) {
		default:
		case 0: // RGB
			gtk_label_set_text(label_c1, _("Red: "));
			gtk_label_set_text(label_c2, _("Green: "));
			gtk_label_set_text(label_c3, _("Blue: "));
			break;
		case 1: // HSL
			gtk_label_set_text(label_c1, _("Hue: "));
			gtk_label_set_text(label_c2, _("Saturation: "));
			gtk_label_set_text(label_c3, _("Lightness: "));
			break;
		case 2: // HSV
			gtk_label_set_text(label_c1, _("Hue: "));
			gtk_label_set_text(label_c2, _("Saturation: "));
			gtk_label_set_text(label_c3, _("Value: "));
			break;
		case 3: // CIE L*a*b*
			gtk_label_set_text(label_c1, "L*: ");
			gtk_label_set_text(label_c2, "a*: ");
			gtk_label_set_text(label_c3, "b*: ");
	}
}

void on_extract_channel_button_ok_clicked(GtkButton *button, gpointer user_data) {
	static GtkEntry *channel_extract_entry[3] = { NULL, NULL, NULL };
	static GtkDropDown *combo_extract_channel = NULL;

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	struct extract_channels_data *args = calloc(1, sizeof(struct extract_channels_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}

	if (combo_extract_channel == NULL) {
		combo_extract_channel = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_extract_colors"));
		channel_extract_entry[0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "Ch1_extract_channel_entry"));
		channel_extract_entry[1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "Ch2_extract_channel_entry"));
		channel_extract_entry[2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "Ch3_extract_channel_entry"));
	}

	args->type = gtk_drop_down_get_selected(combo_extract_channel);
	/* IDs aligned with extract_channel_dialog.ui combo_extract_colors items */
	static const char *const extract_color_ids[] = {"RGB", "HSL", "HSV", "CIE L*a*b*"};
	args->str_type = (args->type < G_N_ELEMENTS(extract_color_ids))
			? extract_color_ids[args->type] : NULL;

	if (args->type != EXTRACT_RGB) {
		// Not RGB, so we need to value_check the image to avoid out-of-range pixels
		if (!value_check(gfit)) {
			siril_log_error(_("Error in value_check(). This should not happen...\n"));
			return;
		}
	}

	args->channel[0] = args->channel[1] = args->channel[2] = NULL;

	for (int i = 0; i < 3; i++) {
	    const gchar *text = gtk_editable_get_text(GTK_EDITABLE(channel_extract_entry[i]));
	    if (text && *text) {
	        args->channel[i] = g_strdup_printf("%s%s", text, com.pref.ext);
	    }
	}

	args->fit = calloc(1, sizeof(fits));
	set_cursor_waiting(TRUE);
	if (copyfits(gfit, args->fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_log_message(_("Could not copy the input image, aborting.\n"));
		clearfits(args->fit);
		free(args->fit);
		free(args->channel[0]);
		free(args->channel[1]);
		free(args->channel[2]);
		free(args);
	} else {
		copy_fits_metadata(gfit, args->fit);
		if (!start_in_new_thread(extract_channels, args)) {
			clearfits(args->fit);
			free(args->fit);
			free(args->channel[0]);
			free(args->channel[1]);
			free(args->channel[2]);
			free(args);
		}
	}
}

void update_button_sensitivity(GtkWidget *entry, gpointer user_data) {
    GtkWidget *button = GTK_WIDGET(user_data);
    static GtkEntry *channel_extract_entry[3] = { NULL };
    if (!channel_extract_entry[0]) {
        channel_extract_entry[0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "Ch1_extract_channel_entry"));
        channel_extract_entry[1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "Ch2_extract_channel_entry"));
        channel_extract_entry[2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "Ch3_extract_channel_entry"));
    }

    gboolean has_text = FALSE;

    for (int i = 0; i < 3; i++) {
        const gchar *text = gtk_editable_get_text(GTK_EDITABLE(channel_extract_entry[i]));
        if (text && *text) {
            has_text = TRUE;
            break;
        }
    }

    gtk_widget_set_sensitive(button, has_text);
}


/* Helper function to read matrix from GUI */
static void get_ccm_values(ccm matrix, float *power) {
	static GtkEntry *entry_m[3][3] = { { NULL } };
	static GtkSpinButton *spin_power = NULL;
	if (!entry_m[0][0]) {
		entry_m[0][0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m00"));
		entry_m[0][1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m01"));
		entry_m[0][2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m02"));
		entry_m[1][0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m10"));
		entry_m[1][1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m11"));
		entry_m[1][2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m12"));
		entry_m[2][0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m20"));
		entry_m[2][1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m21"));
		entry_m[2][2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m22"));
		spin_power = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ccm_power"));
	}
	matrix[0][0] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[0][0])), NULL);
	matrix[0][1] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[0][1])), NULL);
	matrix[0][2] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[0][2])), NULL);
	matrix[1][0] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[1][0])), NULL);
	matrix[1][1] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[1][1])), NULL);
	matrix[1][2] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[1][2])), NULL);
	matrix[2][0] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[2][0])), NULL);
	matrix[2][1] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[2][1])), NULL);
	matrix[2][2] = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry_m[2][2])), NULL);

	if (power) {
		*power = gtk_spin_button_get_value(spin_power);
	}
}

void on_ccm_apply_clicked(GtkButton* button, gpointer user_data) {
	static GtkCheckButton *btn = NULL;
	static GtkEntry *ccmSeqEntry = NULL;
	static GtkWidget *ccm_restore_icc = NULL;
	if (!btn) {
		btn = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "check_apply_seq_ccm"));
		ccmSeqEntry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryCCMSeq"));
		ccm_restore_icc = GTK_WIDGET(gtk_builder_get_object(gui.builder, "ccm_restore_icc"));
	}
	if (ccm_amend_mode) {
		/* Serialize the widget state through the SAME struct and serializer
		 * the normal apply uses, then route it through the amend path.  No
		 * ICC juggling and no sequence path in amend mode. */
		struct ccm_data applied = { 0 };
		get_ccm_values(applied.matrix, &applied.power);
		gchar *blob = op_desc_ccm.serialize(&applied);
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		ccm_amend_mode = FALSE;
		siril_close_dialog("ccm_dialog");
		return;
	}

	struct ccm_data *args = new_ccm_data();
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}

	get_ccm_values(args->matrix, &args->power);

	gboolean seq_toggle = siril_toggle_get_active(GTK_WIDGET(btn));

	if (seq_toggle && sequence_is_loaded()) {
		args->seqEntry = strdup(gtk_editable_get_text(GTK_EDITABLE(ccmSeqEntry)));
		args->seq = &com.seq;
		apply_ccm_to_sequence(args);
	} else {
		if (seq_toggle) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No sequence is loaded"));
			free_ccm_data(args);
			return;
		}

		// Check for ICC profile warning
		if (current_icc_profile() && current_image_color_managed()) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("ICC Profile"),
				_("This image has an attached ICC profile. Applying the CCM will invalidate the "
				"ICC profile therefore color management will be disabled. When you have completed low-level color manipulation and returned the image "
				"to the color space described by its ICC profile you can re-enable it using the button at the bottom of this dialog."));
			color_manage(gfit, FALSE);
			gtk_widget_set_sensitive(ccm_restore_icc, TRUE);
		}

		// Free the args structure as we're using the worker
		ccm temp_matrix;
		float temp_power = args->power;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				temp_matrix[i][j] = args->matrix[i][j];
			}
		}
		free_ccm_data(args);

		// Process using worker
		set_cursor_waiting(TRUE);
		if (ccm_process_with_worker(temp_matrix, temp_power) == 0) {
			invalidate_stats_from_fit(gfit);
		}
		set_cursor_waiting(FALSE);
	}
}

void on_ccm_restore_icc_clicked(GtkButton *button, gpointer user_data) {
	if (current_icc_profile()) {
		color_manage(gfit, TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(button), FALSE);
	}
}

static void update_ccm_matrix(ccm matrix) {
	static GtkEntry *entry_m[3][3] = { { NULL } };
	if (!entry_m[0][0]) {
		entry_m[0][0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m00"));
		entry_m[0][1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m01"));
		entry_m[0][2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m02"));
		entry_m[1][0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m10"));
		entry_m[1][1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m11"));
		entry_m[1][2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m12"));
		entry_m[2][0] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m20"));
		entry_m[2][1] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m21"));
		entry_m[2][2] = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_m22"));
	}
	gchar buf[G_ASCII_DTOSTR_BUF_SIZE+1];
	gtk_editable_set_text(GTK_EDITABLE(entry_m[0][0]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[0][0]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[0][1]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[0][1]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[0][2]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[0][2]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[1][0]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[1][0]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[1][1]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[1][1]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[1][2]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[1][2]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[2][0]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[2][0]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[2][1]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[2][1]));
	gtk_editable_set_text(GTK_EDITABLE(entry_m[2][2]), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[2][2]));
}

void on_ccm_reset_clicked(GtkButton* button, gpointer user_data) {
	ccm matrix = { { 0.f } };
	matrix[0][0] = matrix[1][1] = matrix[2][2] = 1.0f;
	update_ccm_matrix(matrix);
}

void on_combo_ccm_preset_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	static GtkSpinButton *spin_ccm_power = NULL;
	if (!spin_ccm_power) spin_ccm_power = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ccm_power"));
	ccm matrix;
	float power;
	int index = gtk_drop_down_get_selected(combo);
	switch (index) {
		case 0:
			// Custom - does nothing
			return;
		case 1: // Linear Rec.709 to XYZ
			matrix[0][0] = 0.4124564f;
			matrix[0][1] = 0.3575761f;
			matrix[0][2] = 0.1804375f;
			matrix[1][0] = 0.2126729f;
			matrix[1][1] = 0.7151522f;
			matrix[1][2] = 0.0721750f;
			matrix[2][0] = 0.0193339f;
			matrix[2][1] = 0.1191920f;
			matrix[2][2] = 0.9503041f;
			power = 1.f;
			break;
		case 3:
			matrix[0][0] = 0.4124564f;
			matrix[0][1] = 0.3575761f;
			matrix[0][2] = 0.1804375f;
			matrix[1][0] = 0.2126729f;
			matrix[1][1] = 0.7151522f;
			matrix[1][2] = 0.0721750f;
			matrix[2][0] = 0.0193339f;
			matrix[2][1] = 0.1191920f;
			matrix[2][2] = 0.9503041f;
			power = 2.2f;
			break;
		case 2: // XYZ to Linear Rec.709
			matrix[0][0] = 3.2404542f;
			matrix[0][1] = -1.5371385f;
			matrix[0][2] = -0.4985314f;
			matrix[1][0] = -0.9692660f;
			matrix[1][1] = 1.8760108f;
			matrix[1][2] = 0.0415560f;
			matrix[2][0] = 0.0556434f;
			matrix[2][1] = -0.2040259f;
			matrix[2][2] = 1.0572253f;
			power = 1.f;
			break;
		case 4:
			matrix[0][0] = 3.2404542f;
			matrix[0][1] = -1.5371385f;
			matrix[0][2] = -0.4985314f;
			matrix[1][0] = -0.9692660f;
			matrix[1][1] = 1.8760108f;
			matrix[1][2] = 0.0415560f;
			matrix[2][0] = 0.0556434f;
			matrix[2][1] = -0.2040259f;
			matrix[2][2] = 1.0572253f;
			power = 1.f / 2.2f;
			break;
		default:
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), _("This case is not handled yet"));
			return;
	}
	update_ccm_matrix(matrix);
	gtk_spin_button_set_value(spin_ccm_power, power);
}

void on_ccm_close_clicked(GtkButton* button, gpointer user_data) {
	if (ccm_amend_mode) {
		/* Close does not apply — leave amend mode without committing (restore
		 * the true pixels). */
		nde_amend_preview_end(FALSE, NULL);
		ccm_amend_mode = FALSE;
	}
	siril_close_dialog("ccm_dialog");
}

gboolean ccm_hide_on_delete(GtkWidget *widget) {
	if (ccm_amend_mode) {
		nde_amend_preview_end(FALSE, NULL);
		ccm_amend_mode = FALSE;
	}
	siril_close_dialog("ccm_dialog");
	return TRUE;
}

/* Fill the widgets from the amended record's current parameters through the
 * op's own deserializer — the same struct the normal apply reads. */
static void ccm_prefill_from_amend(void) {
	struct ccm_data *p = op_desc_ccm.deserialize(nde_amend_preview_params(),
	                                            nde_amend_preview_op_version());
	if (!p)
		return;
	static GtkSpinButton *spin_ccm_power = NULL;
	if (!spin_ccm_power)
		spin_ccm_power = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_ccm_power"));
	update_ccm_matrix(p->matrix);
	gtk_spin_button_set_value(spin_ccm_power, p->power);
	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
}

void on_ccm_dialog_show(GtkWidget *widget, gpointer user_data) {
	static GtkWidget *ccm_seq_check = NULL, *ccm_seq_entry = NULL, *ccm_restore_icc = NULL;
	static GtkWidget *ccm_amend_note = NULL;
	if (!ccm_amend_note) {
		ccm_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "ccm_amend_note"));
		ccm_seq_check = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_apply_seq_ccm"));
		ccm_seq_entry = GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryCCMSeq"));
		ccm_restore_icc = GTK_WIDGET(gtk_builder_get_object(gui.builder, "ccm_restore_icc"));
	}
	gtk_widget_set_visible(ccm_amend_note, ccm_amend_mode);
	if (ccm_amend_mode) {
		/* gfit already shows the pre-record state.  Sequence controls and the
		 * ICC-restore affordance are meaningless in amend mode. */
		if (gui.roi.active)
			on_clear_roi();
		gtk_widget_set_visible(ccm_seq_check, FALSE);
		gtk_widget_set_visible(ccm_seq_entry, FALSE);
		gtk_widget_set_sensitive(ccm_restore_icc, FALSE);
		set_notify_block(TRUE);
		ccm_prefill_from_amend();
		set_notify_block(FALSE);
	} else {
		gtk_widget_set_visible(ccm_seq_check, TRUE);
		gtk_widget_set_visible(ccm_seq_entry, TRUE);
	}
}

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void ccm_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	ccm_amend_mode = TRUE;
	siril_open_dialog("ccm_dialog");
}

gboolean ccm_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, ccm_amend_ready, NULL);
	return TRUE;
}
