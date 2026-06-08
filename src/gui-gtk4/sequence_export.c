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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/gui_iface.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/stacking.h"
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/save_dialog.h"
#include "io/sequence.h"
#include "io/sequence_export.h"
#include "registration/registration.h"
#ifdef HAVE_LIBTIFF
#include "io/Astro-TIFF.h"
#endif

void on_buttonExportSeq_clicked(GtkButton *button, gpointer user_data) {
	int selected = gtk_drop_down_get_selected(GTK_DROP_DOWN(GTK_WIDGET(gtk_builder_get_object(gui.builder, "comboExport"))));
	GtkEntry *entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryExportSeq")));
	const char *bname = gtk_editable_get_text(GTK_EDITABLE(entry));
	gboolean normalize = siril_toggle_get_active(GTK_WIDGET(GTK_WIDGET(gtk_builder_get_object(gui.builder, "exportNormalize"))));

	if (bname[0] == '\0') {
		widget_set_class(GTK_WIDGET(entry), "warning", "");
		return;
	}
	if (selected == -1) return;

	// checking regdata is absent, or if present, is only shift
	if (!test_regdata_is_valid_and_shift(&com.seq, get_registration_layer(&com.seq))) {
		int confirm = gui_iface.confirm_dialog(_("Registration data found"),
			_("Export has detected registration data with more than simple shifts.\n"
			"Normally, you should apply existing registration before exporting."),
			_("Export anyway"));
		if (!confirm)
			return;
	}

	struct exportseq_args *args = calloc(1, sizeof(struct exportseq_args));
	args->seq = &com.seq;
	get_sequence_filtering_from_gui(&args->filtering_criterion, &args->filtering_parameter);
	args->basename = g_str_to_ascii(bname, NULL);
	args->output = (export_format)selected;
	args->normalize = normalize;
	args->resample = FALSE;
	args->crop = com.selection.w && com.selection.h;
	if (args->crop)
		memcpy(&args->crop_area, &com.selection, sizeof(rectangle));

	if (args->output == EXPORT_AVI || args->output == EXPORT_MP4 || args->output == EXPORT_MP4_H265 || args->output == EXPORT_WEBM_VP9) {
		GtkEntry *fpsEntry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviFps")));
		args->film_fps = round_to_int(g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(fpsEntry)), NULL));
		if (args->film_fps <= 0) args->film_fps = 1;
	}
	if (args->output == EXPORT_MP4 || args->output == EXPORT_MP4_H265 || args->output == EXPORT_WEBM_VP9) {
		GtkAdjustment *adjQual = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder,"adjustmentqualscale"));
		GtkCheckButton *checkResize = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkAviResize")));
		args->film_quality = (int)gtk_adjustment_get_value(adjQual);
		args->resample = siril_toggle_get_active(GTK_WIDGET(checkResize));
		if (args->resample) {
			GtkEntry *widthEntry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviWidth")));
			GtkEntry *heightEntry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviHeight")));
			args->dest_width = g_ascii_strtoll(gtk_editable_get_text(GTK_EDITABLE(widthEntry)), NULL, 10);
			args->dest_height = g_ascii_strtoll(gtk_editable_get_text(GTK_EDITABLE(heightEntry)), NULL, 10);
			if (args->dest_height == 0 || args->dest_width == 0) {
				siril_log_message(_("Width or height cannot be null. Not resizing.\n"));
				siril_toggle_set_active(GTK_WIDGET(checkResize), FALSE);
				args->resample = FALSE;
			} else if (args->dest_height == args->seq->ry && args->dest_width == args->seq->rx) {
				siril_toggle_set_active(GTK_WIDGET(checkResize), FALSE);
				args->resample = FALSE;
			}
		}
	}
	else if (args->output == EXPORT_FITS || args->output == EXPORT_TIFF) {
		// add a trailing '_' for multiple-files sequences
		args->basename = format_basename(args->basename, TRUE);
		if (args->output == EXPORT_TIFF) {
#ifdef HAVE_LIBTIFF
			args->tiff_compression = get_tiff_compression();
#endif
		}
	}
	// Display a useful warning because I always forget to remove selection
	if (args->crop) {
		gboolean confirm = gui_iface.confirm_dialog(_("Export cropped sequence?"),
				_("An active selection was detected. The exported sequence will only contain data within the drawn selection. You can confirm the crop or cancel it. "
						"If you choose to click on cancel, the exported sequence will contain all data."), _("Confirm Crop"));
		args->crop = confirm;
	}

	gui_iface.set_busy(TRUE);
	if (!sequence_export_start(args)) {
		/* args already freed by sequence_export_start on failure */
	}
}

void on_comboExport_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	GtkWidget *avi_options = GTK_WIDGET(gtk_builder_get_object(gui.builder, "boxAviOptions"));
	GtkWidget *checkAviResize = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkAviResize"));
	GtkWidget *quality = GTK_WIDGET(gtk_builder_get_object(gui.builder, "exportQualScale"));
	int output_type = gtk_drop_down_get_selected(box);
	gtk_widget_set_visible(avi_options, output_type >= EXPORT_AVI);
	gtk_widget_set_visible(quality, output_type >= EXPORT_MP4);
	gtk_widget_set_sensitive(checkAviResize, output_type >= EXPORT_MP4);
}

void on_checkAviResize_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	GtkWidget *heightEntry = GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviHeight"));
	GtkWidget *widthEntry = GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviWidth"));
	GtkWidget *combo_export_preset = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combo_export_preset"));
	gtk_widget_set_sensitive(heightEntry, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
	gtk_widget_set_sensitive(widthEntry, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
	gtk_widget_set_sensitive(combo_export_preset, siril_toggle_get_active(GTK_WIDGET(togglebutton)));
}

void update_export_crop_label() {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "exportLabel")));
	if (com.selection.w && com.selection.h)
		gtk_label_set_text(label, _("Cropping to selection"));
	else gtk_label_set_text(label, _("Select area to crop"));
}

void on_entryExportSeq_changed(GtkEditable *editable, gpointer user_data){
	gchar *name = (gchar *)gtk_editable_get_text(GTK_EDITABLE(editable));
	if (*name != 0) {
		if (check_if_seq_exist(name, !g_str_has_suffix(name, ".ser"))) {
			set_icon_entry(GTK_ENTRY(editable), "dialog-warning");
		} else {
			set_icon_entry(GTK_ENTRY(editable), NULL);
		}
	} else {
		set_icon_entry(GTK_ENTRY(editable), NULL);
	}
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data);

void on_entryAviWidth_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_height;
	GtkEntry *heightEntry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviHeight")));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.ry / (double) com.seq.rx;
	width = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(editable)),NULL);
	height = ratio * width;
	c_height = g_strdup_printf("%d", (int)(height));

	g_signal_handlers_block_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	gtk_editable_set_text(GTK_EDITABLE(heightEntry), c_height);
	g_signal_handlers_unblock_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	g_free(c_height);
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_width;
	GtkEntry *widthEntry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviWidth")));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.rx / (double) com.seq.ry;
	height = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(editable)), NULL);
	width = ratio * height;
	c_width = g_strdup_printf("%d", (int)(width));

	g_signal_handlers_block_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	gtk_editable_set_text(GTK_EDITABLE(widthEntry), c_width);
	g_signal_handlers_unblock_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	g_free(c_width);
}

static gboolean g_signal_handlers_is_blocked_by_func(gpointer instance, GFunc func, gpointer data) {
	return g_signal_handler_find(instance,
			G_SIGNAL_MATCH_FUNC | G_SIGNAL_MATCH_DATA
					| G_SIGNAL_MATCH_UNBLOCKED, 0, 0, NULL, func, data) == 0;
}

void on_combo_export_preset_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	int preset = gtk_drop_down_get_selected(box);
	GtkWidget *widthEntry = GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviWidth"));
	GtkWidget *heightEntry = GTK_WIDGET(gtk_builder_get_object(gui.builder, "entryAviHeight"));

	switch(preset) {
	default:
	case PRESET_RATIO:
		if (g_signal_handlers_is_blocked_by_func(GTK_ENTRY(widthEntry), (GFunc) on_entryAviWidth_changed, NULL)) {
			g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
			g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		}
		gtk_widget_set_sensitive(widthEntry, TRUE);
		gtk_widget_set_sensitive(heightEntry, TRUE);
		break;
	case PRESET_FREE:
		if (!g_signal_handlers_is_blocked_by_func(GTK_ENTRY(widthEntry), (GFunc) on_entryAviWidth_changed, NULL)) {
			g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
			g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		}
		gtk_widget_set_sensitive(widthEntry, TRUE);
		gtk_widget_set_sensitive(heightEntry, TRUE);
		break;
	case PRESET_FULLHD:
		g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		gtk_editable_set_text(GTK_EDITABLE(widthEntry), "1920");
		gtk_editable_set_text(GTK_EDITABLE(heightEntry), "1080");
		gtk_widget_set_sensitive(widthEntry, FALSE);
		gtk_widget_set_sensitive(heightEntry, FALSE);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		break;
	case PRESET_4K:
		g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		gtk_editable_set_text(GTK_EDITABLE(widthEntry), "3840");
		gtk_editable_set_text(GTK_EDITABLE(heightEntry), "2160");
		gtk_widget_set_sensitive(widthEntry, FALSE);
		gtk_widget_set_sensitive(heightEntry, FALSE);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		break;
	case PRESET_8K:
		g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		gtk_editable_set_text(GTK_EDITABLE(widthEntry), "7680");
		gtk_editable_set_text(GTK_EDITABLE(heightEntry), "4320");
		gtk_widget_set_sensitive(widthEntry, FALSE);
		gtk_widget_set_sensitive(heightEntry, FALSE);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
	}
}
