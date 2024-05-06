/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/demosaicing.h"
#include "algos/siril_wcs.h"
#include "drizzle/cdrizzleutil.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "gui/registration.h"
#include "gui/sequence_list.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "stacking/stacking.h"

static gboolean keep_noout_state = FALSE;

#undef DEBUG
static char *tooltip_text[] = {
	N_("<b>1-2-3 Stars Registration</b>: This is the simplest method to register deep-sky images. "
		"Images are aligned using shifting, if you pick one star, "
		"or shifting + rotation if 2 or 3 stars are selected.\n"
		"If only shifts are computed, they are saved in the seq file and the aligned sequence does not need to be exported. "
		"Images will be shifted pixel-wise during stacking step.\n"
		"If rotation is also computed, then the aligned images need to be exported using Apply Existing Registration."),
	N_("<b>Global Star Alignment</b>: This is a more powerful and accurate algorithm (but also "
		"slower) to perform deep-sky images. The global matching is based on triangle "
		"similarity method for automatically identify common stars in each image. A new "
		"sequence is created with the prefix of your choice (r_ by default)."),
	N_("<b>Two-Pass Global Star Alignment</b>: The global star alignment is done in two passes, "
		"allowing the reference frame to be chosen from detected star information instead of "
		"automatically choosing the first frame of the sequence."),
	N_("<b>Image Pattern Alignment</b>: This is a simple registration by translation method "
		"using cross correlation in the spatial domain. This method is fast and is used to "
		"register planetary movies. It can also be used for some deep-sky images registration. "
		"Shifts at pixel precision are saved in seq file."),
	N_("<b>KOMBAT</b>: This simple algorithm tries to locate a single pattern on images and to "
		"align them accordingly. Only translation is taken into account yet."),
	N_("<b>Comet/Asteroid Registration</b>: This algorithm is dedicated to the comet and asteroid "
		"registration. It is necessary to have timestamps stored in FITS header and to load a "
		"sequence of star aligned images. This methods makes a translation of a certain number "
		"of pixels depending on the timestamp of each images and the global shift of the "
		"object between the first and the last image."),
	N_("<b>Apply Astrometric Registration</b>: This algorithm computes the transforms between plate-solved images "
	    " of a sequence and applies them"),
	N_("<b>Apply existing registration</b>: This is not an algorithm but rather a commodity to "
		"apply previously computed registration data stored in the sequence file. The "
		"interpolation method and simplified drizzle can be selected in the Output "
		"Registration section and it can be applied on selected images only, to avoid saving "
		"unnecessary images.")
};

static char *reg_frame_registration[] = {
	"framing-default.svg",
	"framing-max.svg",
	"framing-min.svg",
	"framing-cog.svg"
};

/*Possible values for max stars combo box
Needs to be consistent with list in comboreg_maxstars*/
static int maxstars_values[] = { 100, 200, 500, 1000, 2000 };

/* Values for seq filtering for % or k value*/
static float filter_initvals[] = {90., 3.}; // max %, max k
static float filter_maxvals[] = {100., 5.}; // max %, max k
static float filter_increments[] = {1., 0.1}; // spin button steps for % and k
static char *filter_tooltip_text[] = {
	N_("Percents of the images of the sequence."),
	N_("Number of standard deviations for rejection of worst images by k-sigma clipping algorithm.")
};

/* callback for the selected area event */
void _reg_selected_area_callback() {
	if (!com.headless)
		update_reg_interface(TRUE);
}

int populate_drizzle_data(struct driz_args_t *driz) {
	driz->use_flats = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("driz_use_flats")));
	driz->scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_driz_scale")));
	driz->weight_scale = 1.f; // Not used for now
	driz->kernel = (enum e_kernel_t) gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_driz_kernel")));
	driz->pixel_fraction = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_driz_dropsize")));
	if (driz->use_flats) {
		fits reffit = { 0 };
		GtkEntry *entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		const gchar *flat_filename = gtk_entry_get_text(entry);
		gchar *error = NULL;
		int status;
		gchar *expression = path_parse(&reffit, flat_filename, PATHPARSE_MODE_READ, &status);
		if (status) {
			error = _("NOT USING FLAT: could not parse the expression");
			driz->use_flats = FALSE;
		} else {
			free(expression);
			if (flat_filename[0] == '\0') {
				siril_log_message(_("Error: no master flat specified in the preprocessing tab.\n"));
				free(driz);
				return 1;
			} else {
				set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
				driz->flat = calloc(1, sizeof(fits));
				if (!readfits(flat_filename, driz->flat, NULL, TRUE)) {
					if (driz->flat->naxes[2] != com.seq.nb_layers) {
						error = _("NOT USING FLAT: number of channels is different");
					} else if (driz->flat->naxes[0] != com.seq.rx ||
							driz->flat->naxes[1] != com.seq.ry) {
						error = _("NOT USING FLAT: image dimensions are different");
					} else {
						// no need to deal with bitdepth conversion as readfits has already forced conversion to float
						siril_log_message(_("Master flat read for use as initial pixel weight\n"));
					}
				} else error = _("NOT USING FLAT: cannot open the file");
				if (error) {
					siril_log_color_message("%s\n", "red", error);
					set_progress_bar_data(error, PROGRESS_DONE);
					if (driz->flat) {
						clearfits(driz->flat);
						free(driz->flat);
					}
					free(driz);
					return 1;
				}
			}
		}
	}
	return 0;
}

void on_drizzleCheckButton_toggled(GtkToggleButton* button, gpointer user_data) {
	gboolean state = gtk_toggle_button_get_active(button);
	gtk_widget_set_visible(lookup_widget("box_drizzle_controls"), state);
	if (state) {
		gtk_notebook_set_current_page(GTK_NOTEBOOK(lookup_widget("notebook_registration")), REG_PAGE_APPLYREG);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("upscaleCheckButton")), FALSE);
	}
	gtk_widget_set_visible(lookup_widget("interp_box"), !state);
	gtk_widget_set_visible(lookup_widget("toggle_reg_clamp"), !state);
	gtk_widget_set_sensitive(lookup_widget("upscaleCheckButton"), !state);
	gtk_widget_set_visible(lookup_widget("regNoOutput"), FALSE);

}

void on_upscaleCheckButton_toggled(GtkToggleButton* button, gpointer user_data) {
	gboolean state = gtk_toggle_button_get_active(button);
	if (state) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("drizzleCheckButton")), FALSE);
	}
	GtkWidget *regNoOut = lookup_widget("regNoOutput");
	if (gtk_widget_get_visible(regNoOut) && state) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(regNoOut), FALSE);
		gtk_widget_set_sensitive(regNoOut, FALSE);
	} else {
		gtk_widget_set_sensitive(regNoOut, TRUE);
	}
}

static struct registration_method *reg_methods[NUMBER_OF_METHODS + 1];

void initialize_registration_methods() {
	GtkComboBoxText *regcombo;
	int i = 0, j = 0;
	GString *tip;
	gchar *ctip;

	reg_methods[i++] = new_reg_method(_("1-2-3 Stars Registration (deep-sky)"),
			&register_3stars, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Global Star Alignment (deep-sky)"),
			&register_star_alignment, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Two-Pass Global Star Alignment (deep-sky)"),
			&register_multi_step_global, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Image Pattern Alignment (planetary - full disk)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
	reg_methods[i++] = new_reg_method(_("KOMBAT (planetary surfaces or full disk)"),
			&register_kombat, REQUIRES_ANY_SELECTION, REGTYPE_PLANETARY);
	reg_methods[i++] = new_reg_method(_("Comet/Asteroid Registration"),
			&register_comet, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Apply Astrometric Registration"),
			&register_astrometric, REQUIRES_NO_SELECTION, REGTYPE_APPLY);
	reg_methods[i++] = new_reg_method(_("Apply Existing Registration"),
			&register_apply_reg, REQUIRES_NO_SELECTION, REGTYPE_APPLY);
	reg_methods[i] = NULL;

	tip = g_string_new ("");
	for (j = 0; j < i; j ++) {
		tip = g_string_append(tip, _(tooltip_text[j]));
		if (j < i - 1)
			tip = g_string_append(tip, "\n\n");
	}
	ctip = g_string_free (tip, FALSE);
	gtk_widget_set_tooltip_markup(lookup_widget("comboboxregmethod"), ctip);
	g_free(ctip);

	/* fill comboboxregmethod */
	regcombo = GTK_COMBO_BOX_TEXT(lookup_widget("comboboxregmethod"));
	gtk_combo_box_text_remove_all(regcombo);
	i = 0;
	while (reg_methods[i] != NULL) {
		gtk_combo_box_text_append_text(regcombo, reg_methods[i]->name);
		siril_log_message(_("Loading registration method: %s\n"),
				reg_methods[i]->name);
		i++;
	}
	if (i > 0) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(regcombo), com.pref.gui.reg_settings);
	}

	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("ComboBoxRegInter")), com.pref.gui.reg_interpolation);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_reg_clamp")), com.pref.gui.reg_clamping);
	gtk_widget_set_sensitive(lookup_widget("toggle_reg_clamp"),
			com.pref.gui.reg_interpolation == OPENCV_LANCZOS4 || com.pref.gui.reg_interpolation == OPENCV_CUBIC);

	/* register to the new area selected event */
	register_selection_update_callback(_reg_selected_area_callback);
}

struct registration_method *get_selected_registration_method() {
	GtkComboBoxText *regcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(gui.builder, "comboboxregmethod"));
	int index = 0;

	gchar *text = gtk_combo_box_text_get_active_text (regcombo);
	while (reg_methods[index] && text != NULL) {
		if (!strcmp(reg_methods[index]->name, text))
			break;
		index++;
	}
	g_free(text);
	return reg_methods[index];
}

void on_comboboxregmethod_changed(GtkComboBox *box, gpointer user_data) {
	int index = 0;
	gchar *text = gtk_combo_box_text_get_active_text (GTK_COMBO_BOX_TEXT(box));

	while (reg_methods[index] && text != NULL) {
		if (!strcmp(reg_methods[index]->name, text))
			break;
		index++;
	}
	g_free(text);

	com.pref.gui.reg_settings = index;
	reset_3stars();
	update_reg_interface(TRUE);
}

void on_toggle_reg_clamp_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean active = gtk_toggle_button_get_active(button);

	com.pref.gui.reg_clamping = active;
}

void on_ComboBoxRegInter_changed(GtkComboBox *box, gpointer user_data) {
	com.pref.gui.reg_interpolation = gtk_combo_box_get_active(box);
	gtk_widget_set_sensitive(lookup_widget("toggle_reg_clamp"),
			com.pref.gui.reg_interpolation == OPENCV_LANCZOS4
					|| com.pref.gui.reg_interpolation == OPENCV_CUBIC);
}

void on_comboreg_transfo_changed(GtkComboBox *box, gpointer user_data) {
	GtkAdjustment *register_minpairs = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "register_minpairs"));
	double val = gtk_adjustment_get_value(register_minpairs);

	switch(gtk_combo_box_get_active(box)) {
	case SHIFT_TRANSFORMATION:
	case SIMILARITY_TRANSFORMATION:
	case AFFINE_TRANSFORMATION:
		gtk_adjustment_set_lower (register_minpairs, 3);
		break;
	case HOMOGRAPHY_TRANSFORMATION:
		gtk_adjustment_set_lower (register_minpairs, 4);
		if (val < 4)
			gtk_adjustment_set_value(register_minpairs, 4);
		break;
	default:
		printf("on_comboreg_transfo_changed: Value not handled.\n");
	}
}

void on_comboreg_sel_all_combobox_changed(GtkComboBox *box, gpointer user_data) {
	update_reg_interface(TRUE);
}

/****************************************************************/

static struct filtering_tuple regfilters[MAX_FILTERS] = { 0 };

static void update_filters_registration(int update_adjustment);

void on_regsel_changed(GtkComboBox *widget, gpointer user_data) {
	int filter = -1;
	const gchar *caller = gtk_buildable_get_name(GTK_BUILDABLE (widget));
	if (g_str_has_prefix(caller, "filter_type"))
		filter = (int)g_ascii_strtod(caller + 11, NULL) - 4; // filter_type4, 5 or 6 to be parsed as 0, 1 or 2
	update_filters_registration(filter);
}

void on_regspinbut_percent_change(GtkSpinButton *spinbutton, gpointer user_data) {
	update_filters_registration(-1);
}

void on_filter_add4_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter5"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin5"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_add5"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem5"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter5"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_type5"), TRUE);
	update_filters_registration(-1);
}

void on_filter_add5_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter6"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin6"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem6"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter6"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_type6"), TRUE);
	update_filters_registration(-1);
}

void on_filter_rem5_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter5"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin5"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_add5"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem5"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter5"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_type5"), FALSE);
	update_filters_registration(-1);
}

void on_filter_rem6_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter6"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin6"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem6"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter6"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_type6"), FALSE);
	update_filters_registration(-1);
}

void get_reg_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter, int update_adjustment) {
	int filter, guifilter, channel = 0, type;
	gboolean is_ksig = FALSE;
	double percent = 0.0;
	static GtkComboBox *filter_combo[3] = { NULL };
	static GtkAdjustment *stackadj[3] = { NULL };
	static GtkWidget *spin[3] = { NULL };
	static GtkWidget *ksig[3] = { NULL };
	if (!spin[0]) {
		spin[0] = lookup_widget("stackspin4");
		spin[1] = lookup_widget("stackspin5");
		spin[2] = lookup_widget("stackspin6");
		stackadj[0] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[0]));
		stackadj[1] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[1]));
		stackadj[2] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[2]));
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter4"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter5"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter6"));
		ksig[0] = lookup_widget("filter_type4");
		ksig[1] = lookup_widget("filter_type5");
		ksig[2] = lookup_widget("filter_type6");
	}
	for (filter = 0, guifilter = 0; guifilter < 3; guifilter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[guifilter]))) {
			continue;
		}

		type = gtk_combo_box_get_active(filter_combo[guifilter]);
		if (type != ALL_IMAGES && type != SELECTED_IMAGES) {
			channel = get_registration_layer(&com.seq);
			percent = gtk_adjustment_get_value(stackadj[guifilter]);
			is_ksig =  gtk_combo_box_get_active(GTK_COMBO_BOX(ksig[guifilter]));
		}
		if (update_adjustment == filter) {
			g_signal_handlers_block_by_func(stackadj[filter], on_regsel_changed, NULL);
			gtk_adjustment_set_upper(stackadj[filter], filter_maxvals[is_ksig]);
			gtk_adjustment_set_value(stackadj[filter], filter_initvals[is_ksig]);
			gtk_adjustment_set_step_increment(stackadj[filter], filter_increments[is_ksig]);
			gtk_widget_set_tooltip_text(spin[filter], filter_tooltip_text[is_ksig]);
			g_signal_handlers_unblock_by_func(stackadj[filter], on_regsel_changed, NULL);
		}
		switch (type) {
			default:
			case ALL_IMAGES:
				regfilters[filter].filter = seq_filter_all;
				regfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				gtk_widget_set_visible(ksig[guifilter], FALSE);
				break;
			case SELECTED_IMAGES:
				regfilters[filter].filter = seq_filter_included;
				regfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				gtk_widget_set_visible(ksig[guifilter], FALSE);
				break;
			case BEST_PSF_IMAGES:
				regfilters[filter].filter = seq_filter_fwhm;
				regfilters[filter].param = compute_highest_accepted_fwhm(
						&com.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_WPSF_IMAGES:
				regfilters[filter].filter = seq_filter_weighted_fwhm;
				regfilters[filter].param = compute_highest_accepted_weighted_fwhm(
						&com.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_ROUND_IMAGES:
				regfilters[filter].filter = seq_filter_roundness;
				regfilters[filter].param = compute_lowest_accepted_roundness(
						&com.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_BKG_IMAGES:
				regfilters[filter].filter = seq_filter_background;
				regfilters[filter].param = compute_highest_accepted_background(
						&com.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_NBSTARS_IMAGES:
				regfilters[filter].filter = seq_filter_nbstars;
				regfilters[filter].param = compute_lowest_accepted_nbstars(
						&com.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_QUALITY_IMAGES:
				regfilters[filter].filter = seq_filter_quality;
				regfilters[filter].param = compute_lowest_accepted_quality(
						&com.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
		}
		filter++;
	}
	regfilters[filter].filter = NULL;

	if (filter == 1) {
		*filtering_criterion = regfilters[0].filter;
		*filtering_parameter = regfilters[0].param;
	} else {
		*filtering_criterion = create_multiple_filter_from_list(regfilters);
		*filtering_parameter = 0.0;
	}
}

static void update_filter_label(seq_image_filter filtering_criterion, double filtering_parameter) {
	static GtkComboBox *filter_combo[3] = { NULL };
	static GtkLabel *filter_label[3] = { NULL };
	static GtkLabel *result_label = NULL;
	if (!filter_combo[0]) {
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter4"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter5"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter6"));
		filter_label[0] = GTK_LABEL(lookup_widget("labelfilter4"));
		filter_label[1] = GTK_LABEL(lookup_widget("labelfilter5"));
		filter_label[2] = GTK_LABEL(lookup_widget("labelfilter6"));
		result_label = GTK_LABEL(lookup_widget("regfilter_label"));
	}

	for (int filter = 0; filter < 3; filter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[filter]))) {
			break;
		}

		int type = gtk_combo_box_get_active(filter_combo[filter]);
		double param = regfilters[filter].param;
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

	int nb_images_to_stack = compute_nb_filtered_images(&com.seq,
			filtering_criterion, filtering_parameter);
	gchar *labelbuffer = g_strdup_printf(_("Registering %d images of the %d of the sequence"),
			nb_images_to_stack, com.seq.number);
	gtk_label_set_text(result_label, labelbuffer);
	g_free(labelbuffer);
}

static void update_filters_registration(int update_adjustment) {
	if (!sequence_is_loaded())
		return;
	siril_debug_print("updating registration filters GUI\n");
	seq_image_filter criterion;
	double param;
	get_reg_sequence_filtering_from_gui(&criterion, &param, update_adjustment);
	update_filter_label(criterion, param);
}

static gboolean check_framing() {
	// TODO: need to cache
	framing_type framingmethod = (framing_type)gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboreg_framing")));
	GtkLabel *labelreginfo = GTK_LABEL(lookup_widget("labelregisterinfo"));
	if (framingmethod == FRAMING_MAX && com.seq.type == SEQ_FITSEQ) {
		gtk_label_set_text(labelreginfo, _("Max framing not allowed with fitseq, change to regular FITS images"));
		return FALSE;
	}
	// should not happen that often as the process checks the ref image is platesolved, which cannot happen for SER
	// can still be a case if ref image alone is solved through the GUI
	if (com.seq.type == SEQ_SER) {
		gtk_label_set_text(labelreginfo, _("Astrometric registration not allowed for SER sequences, change method"));
		return FALSE;
	}
	gtk_label_set_text(labelreginfo, "");
	gtk_widget_set_tooltip_text(GTK_WIDGET(labelreginfo), "");
	return TRUE;
}
/* Selects the "register all" or "register selected" according to the number of
 * selected images, if argument is false.
 * Verifies that enough images are selected and an area is selected.
 */
void update_reg_interface(gboolean dont_change_reg_radio) {
	static GtkWidget *go_register = NULL, *follow = NULL, *cumul_data = NULL,
	*noout = NULL, *toggle_reg_clamp = NULL, *onlyshift = NULL, *filter_box = NULL, *manualreg = NULL,
	*interpolation_algo = NULL, *undistort_check = NULL, *scale_box = NULL,
	*x2upscale = NULL, *go_estimate = NULL, *drizzle_checkbox = NULL;
	static GtkLabel *labelreginfo = NULL;
	static GtkComboBox *reg_all_sel_box = NULL, *reglayer = NULL, *filter_combo_init = NULL;
	static GtkNotebook *notebook_reg = NULL;
	int nb_images_reg; /* the number of images to register */
	struct registration_method *method = NULL;
	gboolean selection_is_done;
	gboolean has_reg, ready;
	int nbselstars = 0;

	if (!go_register) {
		go_register = lookup_widget("goregister_button");
		go_estimate = lookup_widget("proj_estimate");
		follow = lookup_widget("followStarCheckButton");
		onlyshift = lookup_widget("onlyshift_checkbutton");
		reg_all_sel_box = GTK_COMBO_BOX(lookup_widget("reg_sel_all_combobox"));
		labelreginfo = GTK_LABEL(lookup_widget("labelregisterinfo"));
		notebook_reg = GTK_NOTEBOOK(lookup_widget("notebook_registration"));
		cumul_data = lookup_widget("check_button_comet");
		noout = lookup_widget("regNoOutput");
		reglayer = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
		filter_combo_init = GTK_COMBO_BOX(lookup_widget("combofilter4"));
		toggle_reg_clamp = lookup_widget("toggle_reg_clamp");
		filter_box = lookup_widget("seq_filters_box_reg");
		manualreg = lookup_widget("manualreg_expander");
		interpolation_algo = lookup_widget("ComboBoxRegInter");
		scale_box = lookup_widget("reg_scaling_box");
		undistort_check = lookup_widget("reg_undistort");
		x2upscale = lookup_widget("upscaleCheckButton");
		drizzle_checkbox = lookup_widget("drizzleCheckButton");
	}

	if (!dont_change_reg_radio) {
		if (com.seq.selnum < com.seq.number) {
			gtk_combo_box_set_active(reg_all_sel_box, 1);
			gtk_combo_box_set_active(filter_combo_init, 1);
		} else
			gtk_combo_box_set_active(reg_all_sel_box, 0);
	}

	selection_is_done = (com.selection.w > 0 && com.selection.h > 0);

	/* initialize default */
	gtk_notebook_set_current_page(notebook_reg, REG_PAGE_MISC);
	gtk_widget_set_visible(cumul_data, FALSE);

	/* lock manual registration if sequence is of variable image size*/
	gtk_widget_set_sensitive(manualreg, !com.seq.is_variable);
	gtk_expander_set_expanded(GTK_EXPANDER(manualreg), !com.seq.is_variable);
	gtk_widget_set_tooltip_text(manualreg, (!com.seq.is_variable) ? "" : _("not available for sequences with variable image sizes"));

	/* getting the selected registration method */
	method = get_selected_registration_method();
	if (!method) {
		siril_log_color_message(_("Failed to determine registration method...\n"), "red");
		return;
	}

	/* show the appropriate frame selection widgets */
	gboolean isapplyreg = method->type == REGTYPE_APPLY;
	if (!isapplyreg)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(drizzle_checkbox), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(drizzle_checkbox), isapplyreg);
	gtk_widget_set_sensitive(GTK_WIDGET(drizzle_checkbox), method->method_ptr == &register_apply_reg); // TODO: remove when we allow drizzle with astrometric
	gtk_widget_set_visible(GTK_WIDGET(reg_all_sel_box), !isapplyreg);
	gtk_widget_set_visible(filter_box, isapplyreg);
	gtk_widget_set_visible(GTK_WIDGET(filter_combo_init), isapplyreg);
	if (isapplyreg) {
		if (!dont_change_reg_radio && com.seq.selnum < com.seq.number) {
			gtk_combo_box_set_active(filter_combo_init, SELECTED_IMAGES);
		}
		update_filters_registration(-1);
	}

	/* number of registered image */
	nb_images_reg = gtk_combo_box_get_active(reg_all_sel_box) == 0 ? com.seq.number : com.seq.selnum;

	/* registration data exists for the selected layer */
	has_reg = layer_has_registration(&com.seq, gtk_combo_box_get_active(reglayer));

	if (method && nb_images_reg > 1 && (selection_is_done || method->sel == REQUIRES_NO_SELECTION) && (has_reg || method->type != REGTYPE_APPLY) ) {
		if (method->method_ptr == &register_star_alignment || method->method_ptr == &register_multi_step_global) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_GLOBAL);
		} else if (method->method_ptr == &register_comet) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_COMET);
		} else if (method->method_ptr == &register_3stars) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_3_STARS);
		} else if (method->method_ptr == &register_kombat) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_KOMBAT);
		} else if (method->method_ptr == &register_apply_reg || method->method_ptr == &register_astrometric) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_APPLYREG);
			gtk_widget_set_visible(go_estimate, method->method_ptr == &register_astrometric);
		}
		ready = TRUE;
		if (method->method_ptr == &register_3stars) {
			ready = _3stars_check_selection(); // checks that the right image is loaded based on doall and dofollow
		} else if (gfit.naxes[2] == 1 && gfit.keywords.bayer_pattern[0] != '\0') {
			sensor_pattern pattern = get_bayer_pattern(&gfit);
			if (pattern <= BAYER_FILTER_MAX) {
				gtk_label_set_text(labelreginfo, _("Supported Bayer pattern detected"));
				gtk_widget_set_tooltip_text(GTK_WIDGET(labelreginfo), _("This sequence can be registered with the Bayer pattern intact based on registering the green layer with red and blue pixels interpolated"));
			} else {
				gtk_label_set_text(labelreginfo, _("Unsupported CFA pattern detected"));
				gtk_widget_set_tooltip_text(GTK_WIDGET(labelreginfo), _("This sequence cannot be registered with the CFA pattern intact. You must debayer it prior to registration"));
				ready = FALSE;
			}		
		} else if (method->type == REGTYPE_APPLY && sequence_is_loaded()) {
			ready = check_framing();
		} else {
			gtk_label_set_text(labelreginfo, "");
			gtk_widget_set_tooltip_text(GTK_WIDGET(labelreginfo), "");
		}		// the 3 stars method has special GUI requirements
		if (method->method_ptr == &register_3stars) {
			if (!ready)
				gtk_widget_set_sensitive(go_register, FALSE);
			else
				nbselstars = _3stars_check_registration_ready();
		} else {
			gtk_widget_set_sensitive(go_register, ready);
			gtk_widget_set_sensitive(go_estimate, ready);
		}
		if (method->method_ptr == &register_astrometric && sequence_is_loaded() && !has_wcs(&gfit)) {
			gtk_label_set_text(labelreginfo, _("Platesolve the sequence first"));
			gtk_widget_set_sensitive(go_register, FALSE);
			gtk_widget_set_sensitive(go_estimate, FALSE);
		}

		gtk_widget_set_visible(follow, method->method_ptr == &register_3stars);
		gtk_widget_set_visible(onlyshift, method->method_ptr == &register_3stars);
		gint interpolation_item = gtk_combo_box_get_active(GTK_COMBO_BOX(interpolation_algo));
		gtk_widget_set_sensitive(toggle_reg_clamp, (method->method_ptr == &register_apply_reg ||
		method->method_ptr == &register_star_alignment || method->method_ptr == &register_astrometric ||
		(method->method_ptr == &register_3stars && nbselstars > 1))
		&& (interpolation_item == OPENCV_CUBIC || interpolation_item == OPENCV_LANCZOS4));
		gtk_widget_set_visible(cumul_data, method->method_ptr == &register_comet);
	} else {
		gtk_widget_set_sensitive(go_register, FALSE);
		if (nb_images_reg <= 1 && !selection_is_done) {
			if (sequence_is_loaded()) {
				if (method->sel == REQUIRES_NO_SELECTION) {
					gtk_label_set_text(labelreginfo, _("Select images in the sequence"));
				} else {
					gtk_label_set_text(labelreginfo, _("Select an area in image first, and select images in the sequence"));
				}
			}
			else {
				gtk_label_set_text(labelreginfo, _("Load a sequence first"));
			}
		} else if (nb_images_reg <= 1) {
			gtk_label_set_text(labelreginfo, _("Select images in the sequence"));
		} else if (method->type == REGTYPE_APPLY && !has_reg){
			gtk_label_set_text(labelreginfo, _("Select a layer with existing registration"));
		} else {
			gtk_label_set_text(labelreginfo, _("Select an area in image first"));
		}
	}
	/* we temporary save value as keep_noout_state will be changed in the callback */
	gboolean save_state = keep_noout_state;
	// for now, methods which do not save images but only shift in seq files are constrained to this option (no_output is true and unsensitive)

	gboolean is_astrometric = method->method_ptr == &register_astrometric;
	gtk_widget_set_visible(undistort_check, is_astrometric);
	gtk_widget_set_visible(scale_box, is_astrometric);

	if (((method->method_ptr == &register_comet) ||
			(method->method_ptr == &register_kombat) ||
			(method->method_ptr == &register_shift_dft) ||
			(method->method_ptr == &register_multi_step_global) ||
			(method->method_ptr == &register_3stars && nbselstars <= 1))) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(noout), TRUE);
		gtk_widget_set_sensitive(noout, FALSE);
		gtk_widget_set_visible(noout, TRUE);
		gtk_widget_set_visible(x2upscale, TRUE);
		gtk_widget_set_sensitive(x2upscale, FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(x2upscale), FALSE);
	} else if (method->method_ptr == &register_apply_reg ||
				method->method_ptr == &register_astrometric ) { // cannot have no output with apply registration/astrometric method
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(noout), FALSE);
		gtk_widget_set_sensitive(noout, FALSE);
		gtk_widget_set_visible(noout, FALSE);
		gtk_widget_set_visible(x2upscale, !is_astrometric);
		gtk_widget_set_sensitive(x2upscale, !is_astrometric);
	} else {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(noout), save_state);
		gtk_widget_set_sensitive(noout, TRUE);
		gtk_widget_set_visible(noout, TRUE);
		gtk_widget_set_visible(x2upscale, TRUE);
		gtk_widget_set_sensitive(x2upscale, TRUE);
	}
	keep_noout_state  = save_state;

}


/* callback for no output button */
void on_regNoOutput_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *Algo = lookup_widget("ComboBoxRegInter");
	GtkWidget *clamping = lookup_widget("toggle_reg_clamp");
	GtkWidget *Prefix = lookup_widget("regseqname_entry");
	GtkWidget *x2upscale = lookup_widget("upscaleCheckButton");

	gboolean toggled = gtk_toggle_button_get_active(togglebutton);

	gtk_widget_set_sensitive(Algo, !toggled);
	gtk_widget_set_sensitive(clamping, !toggled);
	gtk_widget_set_sensitive(Prefix, !toggled);
	gtk_widget_set_sensitive(x2upscale, !toggled);
	if (toggled)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(x2upscale), FALSE);

	keep_noout_state = toggled;
}

void on_regfollowStar_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	update_reg_interface(TRUE);
}

void on_shiftonly_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	gboolean toggled = gtk_toggle_button_get_active(togglebutton);
	GtkWidget *noout = lookup_widget("regNoOutput");
	gtk_widget_set_sensitive(noout, !toggled);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(noout), toggled);
}

static int fill_registration_structure_from_GUI(struct registration_args *reg_args) {
	char *msg;
	struct registration_method *method;
	GtkToggleButton *follow, *matchSel, *x2upscale, *cumul, *onlyshift, *undistort;
	GtkComboBox *cbbt_layers, *reg_all_sel_box;
	GtkComboBoxText *ComboBoxRegInter, *ComboBoxTransfo, *ComboBoxMaxStars, *ComboBoxFraming;
	GtkSpinButton *minpairs, *percent_moved, *scaling_spin;

	if (!reserve_thread()) {	// reentrant from here
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	if (!com.seq.regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		// means that a call to seq_check_basic_data() or
		// check_or_allocate_regparam() is missing somewhere else
		return 1;
	}

	method = get_selected_registration_method();

	if (com.selection.w <= 0 && com.selection.h <= 0
			&& method->sel != REQUIRES_NO_SELECTION) {
		msg = siril_log_message(
				_("All prerequisites are not filled for registration. Select a rectangle first.\n"));
		siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), msg);
		return 1;
	}

	control_window_switch_to_tab(OUTPUT_LOGS);

	/* filling the arguments for registration */
	follow = GTK_TOGGLE_BUTTON(lookup_widget("followStarCheckButton"));
	onlyshift = GTK_TOGGLE_BUTTON(lookup_widget("onlyshift_checkbutton"));
	matchSel = GTK_TOGGLE_BUTTON(lookup_widget("checkStarSelect"));
	x2upscale = GTK_TOGGLE_BUTTON(lookup_widget("upscaleCheckButton"));
	cbbt_layers = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
	ComboBoxRegInter = GTK_COMBO_BOX_TEXT(lookup_widget("ComboBoxRegInter"));
	cumul = GTK_TOGGLE_BUTTON(lookup_widget("check_button_comet"));
	minpairs = GTK_SPIN_BUTTON(lookup_widget("spinbut_minpairs"));
	percent_moved = GTK_SPIN_BUTTON(lookup_widget("spin_kombat_percent"));
	ComboBoxMaxStars = GTK_COMBO_BOX_TEXT(lookup_widget("comboreg_maxstars"));
	ComboBoxTransfo = GTK_COMBO_BOX_TEXT(lookup_widget("comboreg_transfo"));
	ComboBoxFraming = GTK_COMBO_BOX_TEXT(lookup_widget("comboreg_framing"));
	reg_all_sel_box = GTK_COMBO_BOX(GTK_COMBO_BOX_TEXT(lookup_widget("reg_sel_all_combobox")));
	scaling_spin =GTK_SPIN_BUTTON(lookup_widget("reg_scaling_spin"));
	undistort =  GTK_TOGGLE_BUTTON(lookup_widget("reg_undistort"));

	reg_args->clamp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_reg_clamp")));

	reg_args->func = method->method_ptr;
	reg_args->seq = &com.seq;
	reg_args->reference_image = sequence_find_refimage(&com.seq);
	reg_args->follow_star = gtk_toggle_button_get_active(follow);
	reg_args->matchSelection = gtk_toggle_button_get_active(matchSel);
	reg_args->no_output = keep_noout_state;
	reg_args->x2upscale = gtk_toggle_button_get_active(x2upscale);
	reg_args->cumul = gtk_toggle_button_get_active(cumul);
	reg_args->prefix = strdup( gtk_entry_get_text(GTK_ENTRY(lookup_widget("regseqname_entry"))));
	reg_args->min_pairs = gtk_spin_button_get_value_as_int(minpairs);
	reg_args->percent_moved = (float) gtk_spin_button_get_value(percent_moved) / 100.f;
	int starmaxactive = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxMaxStars));
	reg_args->max_stars_candidates = (starmaxactive == -1) ? MAX_STARS_FITTED : maxstars_values[starmaxactive];
	if (method->method_ptr != register_3stars)
		reg_args->type = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxTransfo));
	else {
		reg_args->type = (gtk_toggle_button_get_active(onlyshift)) ? SHIFT_TRANSFORMATION : SIMILARITY_TRANSFORMATION;
		reg_args->no_output = (gtk_toggle_button_get_active(onlyshift)) ? TRUE : keep_noout_state;
	}
	reg_args->framing = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxFraming));
	reg_args->undistort = gtk_toggle_button_get_active(undistort);
	reg_args->astrometric_scale = (float)gtk_spin_button_get_value(scaling_spin);
#ifndef HAVE_CV44
	if (reg_args->type == SHIFT_TRANSFORMATION && method->method_ptr != register_3stars) {
		siril_log_color_message(_("Shift-only registration is only possible with OpenCV 4.4\n"), "red");
		free(reg_args->prefix);
		return 1;
	}
#endif
	if (method->method_ptr == register_apply_reg || method->method_ptr == register_astrometric) {
		get_reg_sequence_filtering_from_gui(
				&reg_args->filtering_criterion, &reg_args->filtering_parameter, -1);
	} else {
		reg_args->filters.filter_included = gtk_combo_box_get_active(reg_all_sel_box);
	}
	if ((method->method_ptr == register_star_alignment || method->method_ptr == register_multi_step_global) &&
		reg_args->matchSelection && reg_args->seq->is_variable) {
		siril_log_color_message(_("Cannot use area selection on a sequence with variable image sizes\n"), "red");
		return 1;
	}

	if ((method->method_ptr == register_star_alignment || method->method_ptr == register_multi_step_global) &&
		!reg_args->matchSelection) {
		delete_selected_area(); // otherwise it is enforced
	}

	/* We check that available disk space is enough when
	the registration method produces a new sequence
	*/
	if (!reg_args->no_output && method->method_ptr == register_star_alignment) {

		int nb_frames = reg_args->filters.filter_included ? reg_args->seq->selnum : reg_args->seq->number;
		gint64 size = seq_compute_size(reg_args->seq, nb_frames, get_data_type(reg_args->seq->bitpix));
		if (reg_args->x2upscale)
			size *= 4;
		if (test_available_space(size)) {
			siril_log_color_message(_("Not enough space to save the output images, aborting\n"), "red");
			return 1;
		}
	} else if (method->method_ptr == register_comet) {
		pointf velocity = get_velocity();
		if ((velocity.x == 0.0 && velocity.y == 0.0)
				|| isinf(velocity.x) || isinf(velocity.y)) {
			msg = siril_log_color_message(_("The object is not moving, please check your registration data.\n"), "red");
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), msg);
			return 1;
		}
	}
	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	reg_args->interpolation = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxRegInter));
	get_the_registration_area(reg_args, method);	// sets selection
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE; // only TRUE for some methods. Will be updated in these cases

	if (method->method_ptr == register_star_alignment) { // seqpplyreg case is dealt with in the sanity checks of the method
		if (reg_args->interpolation == OPENCV_NONE && (reg_args->x2upscale || com.seq.is_variable)) {
			siril_log_color_message(_("When interpolation is set to None, the images must be of same size and no upscaling can be applied. Aborting\n"), "red");
			return 1;
		}
		if (reg_args->interpolation == OPENCV_NONE && (reg_args->type > SHIFT_TRANSFORMATION)) {
			siril_log_color_message(_("When interpolation is set to None, the transformation can only be set to Shift. Aborting\n"), "red");
			return 1;
		}
	}
	if (((method->method_ptr == register_star_alignment || method->method_ptr == register_3stars || method->method_ptr == register_apply_reg || method->method_ptr == register_astrometric) &&
		(reg_args->interpolation == OPENCV_AREA || reg_args->interpolation == OPENCV_LINEAR || reg_args->interpolation == OPENCV_NEAREST || reg_args->interpolation == OPENCV_NONE)) ||
		reg_args->no_output)
		reg_args->clamp = FALSE;

	if (method->method_ptr != register_3stars)
		clear_stars_list(TRUE); //to avoid problems with com.stars later on in the process

	return 0;
}

/* callback for 'Go register' button, GTK thread */
void on_seqregister_button_clicked(GtkButton *button, gpointer user_data) {

	struct registration_args *reg_args = calloc(1, sizeof(struct registration_args));

	char *msg;
	if (fill_registration_structure_from_GUI(reg_args)) {
		free(reg_args);
		unreserve_thread();
		return;
	}
	fits fit_ref = { 0 };
	int ret = seq_read_frame_metadata(reg_args->seq, reg_args->reference_image, &fit_ref);
	if (ret) {
		siril_log_message(_("Error: unable to read reference frame metadata\n"));
		free(reg_args);
		unreserve_thread();
		return;
	}
	reg_args->bayer = (fit_ref.keywords.bayer_pattern[0] != '\0');
	clearfits(&fit_ref);

	struct registration_method *method = get_selected_registration_method();

	const gchar *caller = gtk_buildable_get_name(GTK_BUILDABLE(button));
	if (!g_strcmp0(caller, "proj_estimate"))
		reg_args->no_output = TRUE;
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("drizzleCheckButton")))) {
		reg_args->driz = calloc(1, sizeof(struct driz_args_t));
		if (populate_drizzle_data(reg_args->driz)) {
			free(reg_args);
			return;
		}
	}

	msg = siril_log_color_message(_("Registration: processing using method: %s\n"),
			"green", method->name);
	msg[strlen(msg) - 1] = '\0';

	if (!reg_args->driz && reg_args->clamp && !reg_args->no_output)
		siril_log_message(_("Interpolation clamping active\n"));
	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_reserved_thread(register_thread_func, reg_args);
}

// end of registration, GTK thread. Executed when started from the GUI and in
// the graphical command line but not from a script (headless mode)
gboolean end_register_idle(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	stop_processing_thread();

	if (!args->retval) {
		if (!args->load_new_sequence && sequence_is_loaded()) {
			int chan = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboboxreglayer")));
			update_seqlist(chan);
			fill_sequence_list(args->seq, chan, FALSE);
			set_layers_for_registration();	// update display of available reg data
			seq_load_image(args->seq, args->seq->reference_image, TRUE);
			redraw(REDRAW_OVERLAY); // plot registration frame
		}
		else {
			check_seq();
			update_sequences_list(args->new_seq_name);
		}
	}
	set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);
	args->seq->reg_invalidated = FALSE;
	drawPlot();
	update_stack_interface(TRUE);
	adjust_sellabel();

	set_cursor_waiting(FALSE);
	if (args->func == &register_3stars) reset_3stars();

	free(args->new_seq_name);
	if (!check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);
	free(args->driz);
	free(args);
	return FALSE;
}

void on_comboreg_framing_changed(GtkComboBox *box, gpointer user_data) {
	gchar *name;
	GtkImage *image = GTK_IMAGE(lookup_widget("framing-image"));
	int i = gtk_combo_box_get_active(box);

	if (i >= 0 && i < G_N_ELEMENTS(reg_frame_registration)) {
		name = g_strdup_printf("/org/siril/ui/pixmaps/%s", reg_frame_registration[i]);
		gtk_image_set_from_resource(image, name);

		g_free(name);
	}
	GtkWidget *go_register = lookup_widget("goregister_button");
	GtkWidget *go_estimate = lookup_widget("proj_estimate");
	gboolean ready = check_framing();
	gtk_widget_set_sensitive(go_register, ready);
	gtk_widget_set_sensitive(go_estimate, ready);
}
