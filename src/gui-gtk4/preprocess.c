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

/* GTK callbacks and widget-reading helpers for the preprocessing panel.
 * Moved here from core/preprocess.c to keep core/ free of GTK dependencies. */

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/proto.h"
#include "gui-gtk4/utils.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "core/preprocess.h"
#include "gui-gtk4/gui_state.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "filters/cosmetic_correction.h"
#include "algos/statistics.h"
#include "core/arithm.h"

/* Reads widget values to populate a preprocessing_data struct.
 * Returns TRUE on error. */
static gboolean test_for_master_files(struct preprocessing_data *args) {
	GtkCheckButton *tbutton;
	GtkEntry *entry;
	gboolean has_error = FALSE;
	args->bias_level = FLT_MAX;
	args->bad_pixel_map_file = NULL;
	fits reffit = { 0 };
	gboolean isseq = (args->seq != NULL);
	if (isseq) {
		// loading the sequence reference image's metadata in case it's needed
		int image_to_load = sequence_find_refimage(args->seq);
		if (seq_read_frame_metadata(args->seq, image_to_load, &reffit)) {
			siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
			has_error = TRUE;
		}
	} else {
		reffit = *gfit;
	}

	tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "useoffset_button")));
	if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
		const char *filename;
		entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "offsetname_entry")));
		filename = gtk_editable_get_text(GTK_EDITABLE(entry));
		if (filename[0] == '\0') {
			siril_toggle_set_active(GTK_WIDGET(tbutton), FALSE);
		} else {
			const char *error = NULL;
			if (filename[0] == '=') { // offset is specified as a level not a file
				gui_iface.set_progress(PROGRESS_NONE, _("Checking offset level..."));
				int offsetlevel = evaluateoffsetlevel(filename + 1, gfit);
				if (!offsetlevel) {
					error = _("NOT USING OFFSET: the offset value could not be parsed");
					args->use_bias = FALSE;
				} else {
					siril_log_message(_("Synthetic offset: Level = %d\n"),offsetlevel);
					int maxlevel = (gfit->orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
					if ((offsetlevel > maxlevel) || (offsetlevel < -maxlevel) ) {   // not excluding all neg values here to allow defining a pedestal
						error = _("NOT USING OFFSET: the offset value is not consistent with image bitdepth");
						args->use_bias = FALSE;
					} else {
						args->bias_level = (float)offsetlevel;
						args->bias_level *= (gfit->orig_bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE; //converting to [0 1] to use with soper
						args->use_bias = TRUE;
					}
				}
			} else {
				int status;
				gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
				if (status) {
					error = _("NOT USING OFFSET: could not parse the expression\n");
					has_error = TRUE;
				} else {
					args->bias = calloc(1, sizeof(fits));
					gui_iface.set_progress(PROGRESS_NONE, _("Opening offset image..."));
					if (!readfits(expression, args->bias, NULL, !com.pref.force_16bit)) {
						if (args->bias->naxes[2] != gfit->naxes[2]) {
							error = _("NOT USING OFFSET: number of channels is different");
						} else if (args->bias->naxes[0] != gfit->naxes[0] ||
								args->bias->naxes[1] != gfit->naxes[1]) {
							error = _("NOT USING OFFSET: image dimensions are different");
						} else {
							args->use_bias = TRUE;
							// if input is 8b, we assume 32b master needs to be rescaled
							if ((args->bias->type == DATA_FLOAT) && (gfit->orig_bitpix == BYTE_IMG)) {
								soper(args->bias, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
							}
						}
					} else
						error = _("NOT USING OFFSET: cannot open the file");
				}
				g_free(expression);
			}
			if (error) {
				siril_log_error("%s\n", error);
				gui_iface.set_progress(PROGRESS_DONE, error);
				if (args->bias)
					free(args->bias);
				args->use_bias = FALSE;
				has_error = TRUE;
			}
		}
	}

	tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "usedark_button")));
	if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
		const char *filename;
		entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "darkname_entry")));
		filename = gtk_editable_get_text(GTK_EDITABLE(entry));
		if (filename[0] == '\0') {
			siril_toggle_set_active(GTK_WIDGET(tbutton), FALSE);
		} else {
			const char *error = NULL;
			int status;
			gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
			if (status) {
				error = _("NOT USING DARK: could not parse the expression");
				has_error = TRUE;
			} else {
				gui_iface.set_progress(PROGRESS_NONE, _("Opening dark image..."));
				args->dark = calloc(1, sizeof(fits));
				if (!readfits(expression, args->dark, NULL, !com.pref.force_16bit)) {
					if (args->dark->naxes[2] != gfit->naxes[2]) {
						error = _("NOT USING DARK: number of channels is different");
					} else if (args->dark->naxes[0] != gfit->naxes[0] ||
							args->dark->naxes[1] != gfit->naxes[1]) {
						error = _("NOT USING DARK: image dimensions are different");
					} else {
						args->use_dark = TRUE;
						// if input is 8b, we assume 32b master needs to be rescaled
						if ((args->dark->type == DATA_FLOAT) && (gfit->orig_bitpix == BYTE_IMG)) {
							soper(args->dark, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
						}
					}

				} else error = _("NOT USING DARK: cannot open the file");
			}
			if (error) {
				siril_log_error("%s\n", error);
				gui_iface.set_progress(PROGRESS_DONE, error);
				if (args->dark)
					clearfits(args->dark);
				args->dark = NULL; // in order to be sure it is freed
				args->use_dark = FALSE;
				has_error = TRUE;
			}
			g_free(expression);

		}

		if (args->use_dark) {
			// dark optimization
			int optim = gtk_drop_down_get_selected(GTK_DROP_DOWN(GTK_WIDGET(gtk_builder_get_object(gui.builder, "comboDarkOptimize"))));
			args->use_dark_optim = optim != 0;
			args->use_exposure = optim == 2;
			const char *error = NULL;
			if ((gfit->orig_bitpix == BYTE_IMG) && args->use_dark_optim) {
				error = _("Dark optimization: This process cannot be applied to 8b images");
			}
			if (error) {
				siril_log_error("%s\n", error);
				gui_iface.set_progress(PROGRESS_DONE, error);
				clearfits(args->dark);
				gtk_editable_set_text(GTK_EDITABLE(entry), "");
				args->use_dark_optim = FALSE;
				has_error = TRUE;
			}

			// cosmetic correction
			tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "cosmEnabledCheck")));
			args->use_cosmetic_correction = siril_toggle_get_active(GTK_WIDGET(tbutton));

			if (args->use_cosmetic_correction) {
				tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkSigCold")));
				if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
					GtkSpinButton *sigCold = GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spinSigCosmeColdBox")));
					args->sigma[0] = gtk_spin_button_get_value(sigCold);
				} else args->sigma[0] = -1.0;

				tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkSigHot")));
				if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
					GtkSpinButton *sigHot = GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spinSigCosmeHotBox")));
					args->sigma[1] = gtk_spin_button_get_value(sigHot);
				} else args->sigma[1] = -1.0;

				/* Using Bad Pixel Map ? */
				const gchar *bad_pixel_f;
				entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixelmap_entry")));
				bad_pixel_f = gtk_editable_get_text(GTK_EDITABLE(entry));
				/* test for file */
				if (bad_pixel_f[0] != '\0') {
					args->bad_pixel_map_file = g_file_new_for_path(bad_pixel_f);
					if (!check_for_cosme_file_sanity(args->bad_pixel_map_file)) {
						siril_log_debug("cosme file sanity check failed...\n");
						has_error = TRUE;
					}
				}
			}
		}
	}

	/* now we want to know which cosmetic correction we choose */
	GtkStack *stack = GTK_STACK(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_cc")));
	GtkWidget *w = gtk_stack_get_visible_child(stack);
	args->cc_from_dark = TRUE; // default value
	if (w) {
		/* GTK4: gtk_stack_get_pages returns GtkSelectionModel* (which
		 * implements GListModel). */
		GtkSelectionModel *pages_sel = gtk_stack_get_pages(stack);
		GListModel *pages = G_LIST_MODEL(pages_sel);
		guint n = g_list_model_get_n_items(pages);
		for (guint i = 0; i < n; i++) {
			GtkStackPage *page = g_list_model_get_item(pages, i);
			if (gtk_stack_page_get_child(page) == w) {
				args->cc_from_dark = (i == 0);
				g_object_unref(page);
				break;
			}
			g_object_unref(page);
		}
		g_object_unref(pages_sel);
	}

	tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "useflat_button")));
	if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
		const char *filename;
		entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "flatname_entry")));
		filename = gtk_editable_get_text(GTK_EDITABLE(entry));
		if (filename[0] == '\0') {
			siril_toggle_set_active(GTK_WIDGET(tbutton), FALSE);
		} else {
			const char *error = NULL;
			int status;
			gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
			if (status) {
				error = _("NOT USING FLAT: could not parse the expression");
				has_error = TRUE;
			} else {
				gui_iface.set_progress(PROGRESS_NONE, _("Opening flat image..."));
				args->flat = calloc(1, sizeof(fits));
				if (!readfits(expression, args->flat, NULL, !com.pref.force_16bit)) {
					if (args->flat->naxes[2] != gfit->naxes[2]) {
						error = _("NOT USING FLAT: number of channels is different");
					} else if (args->flat->naxes[0] != gfit->naxes[0] ||
							args->flat->naxes[1] != gfit->naxes[1]) {
						error = _("NOT USING FLAT: image dimensions are different");
					} else {
						args->use_flat = TRUE;
						// no need to deal with bitdepth conversion as flat is just a division (unlike darks which need to be on same scale)
					}

				} else error = _("NOT USING FLAT: cannot open the file");
			}
			g_free(expression); // expression not used again after here, free before it falls out of scope
			if (error) {
				siril_log_error("%s\n", error);
				gui_iface.set_progress(PROGRESS_DONE, error);
				if (args->flat)
					free(args->flat);
				args->use_flat = FALSE;
				has_error = TRUE;
			}

			if (args->use_flat) {
				GtkCheckButton *autobutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkbutton_auto_evaluate")));
				args->autolevel = siril_toggle_get_active(GTK_WIDGET(autobutton));
				if (!args->autolevel) {
					GtkEntry *norm_entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entry_flat_norm")));
					args->normalisation = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(norm_entry)), NULL);
				}
			}
		}
	}
	if (isseq)
		clearfits(&reffit);
	return has_error;
}

void on_prepro_button_clicked(GtkButton *button, gpointer user_data) {
	if (!single_image_is_loaded() && !sequence_is_loaded())
		return;
	if (!single_image_is_loaded() && processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	GtkEntry *entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "preproseqname_entry")));
	GtkCheckButton *fix_xtrans = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "fix_xtrans_af")));
	GtkCheckButton *CFA = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "cosmCFACheck")));
	GtkCheckButton *debayer = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkButton_pp_dem")));
	GtkCheckButton *equalize_cfa = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkbutton_equalize_cfa")));
	GtkCheckButton *preprocess_excluded = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "toggle_preprocess_excluded")));
	GtkDropDown *output_type = GTK_DROP_DOWN(GTK_WIDGET(gtk_builder_get_object(gui.builder, "prepro_output_type_combo")));

	struct preprocessing_data *args = calloc(1, sizeof(struct preprocessing_data));
	if (test_for_master_files(args)) {
		siril_log_error(_("Some errors have been detected, Please check the logs.\n"));
		free(args);
		return;
	}

	siril_log_info(_("Preprocessing...\n"));

	// set output filename (preprocessed file name prefix)
	args->ppprefix = strdup(gtk_editable_get_text(GTK_EDITABLE(entry)));

	args->ignore_exclusion = siril_toggle_get_active(GTK_WIDGET(preprocess_excluded));
	args->is_cfa = siril_toggle_get_active(GTK_WIDGET(CFA));
	args->debayer = siril_toggle_get_active(GTK_WIDGET(debayer));
	args->equalize_cfa = siril_toggle_get_active(GTK_WIDGET(equalize_cfa));
	args->fix_xtrans = siril_toggle_get_active(GTK_WIDGET(fix_xtrans));

	/****/

	if (sequence_is_loaded()) {
		args->is_sequence = TRUE;
		args->seq = &com.seq;
		args->output_seqtype = gtk_drop_down_get_selected(output_type);
		if (args->output_seqtype < 0 || args->output_seqtype > SEQ_FITSEQ)
			args->output_seqtype = SEQ_REGULAR;
		args->allow_32bit_output = !com.pref.force_16bit && args->output_seqtype != SEQ_SER;
		gui_iface.set_busy(TRUE);
		gui_iface.show_panel("output_logs", TRUE);
		start_sequence_preprocessing(args);
	} else {
		int retval;
		args->is_sequence = FALSE;
		args->output_seqtype = SEQ_REGULAR;
		args->allow_32bit_output = !com.pref.force_16bit;
		gui_iface.set_busy(TRUE);
		gui_iface.show_panel("output_logs", TRUE);

		retval = calibrate_single_image(args);

		free(args);

		if (retval)
			gui_iface.set_progress(PROGRESS_NONE, _("Error in preprocessing."));
		else {
			gui_iface.set_progress(PROGRESS_RESET, PROGRESS_TEXT_RESET);
			gui_iface.invalidate_histogram();
			gui_iface.open_single_image_from_gfit();
		}
		gui_iface.set_busy(FALSE);
	}
}

void on_GtkButtonEvaluateCC_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;
	GtkLabel *label[2];
	GtkWidget *widget[2];
	const char *filename;
	gchar *str[2];
	double sig[2];
	long icold = 0L, ihot = 0L;
	double rate, total;
	fits fit = { 0 };
	fits reffit = { 0 };
	int status;
	gboolean isseq = FALSE;
	if (sequence_is_loaded()) {
		// loading the sequence reference image's metadata
		int image_to_load = sequence_find_refimage(&com.seq);
		if (seq_read_frame_metadata(&com.seq, image_to_load, &reffit)) {
			siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
			return;
		}
		isseq = TRUE;
	} else {
		reffit = *gfit;
	}

	gui_iface.set_busy(TRUE);
	sig[0] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spinSigCosmeCold"))));
	sig[1] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spinSigCosmeHot"))));
	widget[0] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "GtkLabelColdCC"));
	widget[1] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "GtkLabelHotCC"));
	label[0] = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "GtkLabelColdCC")));
	label[1] = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "GtkLabelHotCC")));
	entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "darkname_entry")));
	filename = gtk_editable_get_text(GTK_EDITABLE(entry));
	gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
	if (status || readfits(expression, &fit, NULL, !com.pref.force_16bit)) {
		str[0] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		str[1] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		gtk_label_set_markup(label[0], str[0]);
		gtk_label_set_markup(label[1], str[1]);
		gui_iface.set_busy(FALSE);
		g_free(expression);
		return;
	}
	find_deviant_pixels(&fit, sig, &icold, &ihot, TRUE);
	printf("sig[0]=%lf\n", sig[0]);
	printf("sig[1]=%lf\n", sig[1]);
	total = fit.rx * fit.ry;
	clearfits(&fit);
	rate = (double)icold / total;
	/* 1% of cold pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[0] = g_markup_printf_escaped("<span foreground=\"red\">%ld px</span>", icold);
		gtk_widget_set_tooltip_text(widget[0], _("This value may be too high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[0] = g_markup_printf_escaped("%ld px", icold);
		gtk_widget_set_tooltip_text(widget[0], "");
	}
	gtk_label_set_markup(label[0], str[0]);

	rate = (double)ihot / total;
	/* 1% of hot pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[1] = g_markup_printf_escaped("<span foreground=\"red\">%ld px</span>", ihot);
		gtk_widget_set_tooltip_text(widget[1], _("This value may be too high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[1] = g_markup_printf_escaped("%ld px", ihot);
		gtk_widget_set_tooltip_text(widget[1], "");
	}
	gtk_label_set_markup(label[1], str[1]);
	g_free(str[0]);
	g_free(str[1]);
	gui_iface.set_busy(FALSE);
	g_free(expression);
	if (isseq)
		clearfits(&reffit);
}
