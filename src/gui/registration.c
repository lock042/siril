/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "gui/registration.h"
#include "gui/sequence_list.h"
#include "gui/stacking.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "stacking/stacking.h"
#include "opencv/opencv.h"

#undef DEBUG
static char *tooltip_text[] = {
	N_("<b>Global Star Alignment</b>: This is the go-to algorithm for deep-sky images"
		"The global matching is based on triangle similarity method to automatically identify "
		"common stars between images. A new sequence is created with the prefix of your choice "
		"(r_ by default), unless you tick the 2pass check box. In that case, registration "
		"data will be saved to the sequence and reference will be picked based on star detection "
		" information. When choosing 2pass, you will then need to export the sequence using "
		" the Apply Existing Registration method"),
	N_("<b>1-2-3 Stars Registration</b>: This is the simplest method to register deep-sky images. "
		"Images are aligned using shifting, if you pick one star, "
		"or shifting + rotation if 2 or 3 stars are selected.\n"
		"Transformations are saved to the input sequence file."),
	N_("<b>Image Pattern Alignment</b>: This is a simple registration by translation method "
		"using cross correlation in the spatial domain. This method is fast and is used to "
		"register planetary movies. It can also be used for some deep-sky images registration. "
		"Shifts at pixel precision are saved to the input sequence file."),
	N_("<b>KOMBAT</b>: This simple algorithm tries to locate a single pattern on images and to "
		"align them accordingly. Only translation is taken into account. "
		"Shifts at pixel precision are saved to the input sequence file."),
	N_("<b>Comet/Asteroid Registration</b>: This algorithm is dedicated to the comet and asteroid "
		"registration. It is necessary to have timestamps stored in FITS header and to load a "
		"sequence of star aligned images. This methods makes a translation of a certain number "
		"of pixels depending on the timestamp of each image and the global shift of the "
		"object between the first and the last image. If some registration data already exists,"
		" the shifts are composed with the existing transformations. A new sequence is created, "
		"with the prefix of your choice (comet_ by default), but the images are not duplicated, they "
		"are symlinked to the input images."),
	N_("<b>Apply existing registration</b>: This is not an algorithm but rather a commodity to "
		"apply previously computed registration data stored in the sequence file. The "
		"interpolation method or drizzling can be selected in the Output "
		"Registration section and it can be applied on selected/filtered images only, to avoid saving "
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
static struct filtering_tuple regfilters[MAX_FILTERS] = { 0 };

// comet specific statics
static pointf velocity = { 0.f, 0.f };
static GDateTime *t_of_image_1 = NULL;
static GDateTime *t_of_image_2 = NULL;
static point pos_of_image1 = { 0 };
static point pos_of_image2 = { 0 };

static struct registration_method *reg_methods[NUMBER_OF_METHODS + 1];

// Statics declarations
static GtkAdjustment *register_minpairs = NULL;
static GtkBox *seq_filters_box_reg = NULL, *reg_wcsfilechooser_box = NULL;
static GtkButton *filter_add4 = NULL, *filter_add5 = NULL, *filter_rem5 = NULL, *filter_rem6 = NULL, *proj_estimate = NULL, *goregister_button = NULL, *reg_wcsfile_button = NULL;
static GtkComboBoxText *comboboxregmethod = NULL, *comboboxreglayer = NULL, *comboreg_maxstars = NULL, *comboreg_transfo = NULL, *reg_sel_all_combobox = NULL, *combofilter4 = NULL, *filter_type4 = NULL, *combofilter5 = NULL, *filter_type5 = NULL, *combofilter6 = NULL, *filter_type6 = NULL, *comboreg_framing = NULL, *ComboBoxRegInter = NULL, *combo_driz_kernel = NULL, *comboreg_undistort = NULL;
static GtkEntry *entry1_x_comet = NULL, *entry2_x_comet = NULL, *entry1_y_comet = NULL, *entry2_y_comet = NULL, *regseqname_entry = NULL, *flatname_entry = NULL, *reg_wcsfile_entry = NULL, *cometseqname_entry = NULL;
static GtkExpander *autoreg_expander = NULL, *manualreg_expander = NULL;
static GtkFrame *output_reg_frame = NULL;
static GtkGrid *grid_reg_framing = NULL, *grid_interp_controls = NULL, *grid_drizzle_controls = NULL, *grid_reg_wcs = NULL;
static GtkImage *framing_image = NULL;
static GtkLabel *label1_comet = NULL, *regfilter_label = NULL, *labelfilter4 = NULL, *labelfilter5 = NULL, *labelfilter6 = NULL, *labelregisterinfo = NULL, *labelRegRef = NULL, *estimate_label = NULL;
static GtkNotebook *notebook_registration = NULL;
static GtkSpinButton *spinbut_minpairs = NULL, *spin_kombat_percent = NULL, *stackspin4 = NULL, *stackspin5 = NULL, *stackspin6 = NULL, *reg_scaling_spin = NULL, *spin_driz_dropsize = NULL, *spinbut_shiftx = NULL, *spinbut_shifty = NULL;
static GtkStack *interp_drizzle_stack = NULL;
static GtkStackSwitcher *interp_drizzle_stack_switcher = NULL;
static GtkToggleButton *checkStarSelect = NULL, *reg_2pass = NULL, *followStarCheckButton = NULL, *onlyshift_checkbutton = NULL, *toggle_reg_clamp = NULL, *driz_use_flats = NULL, *checkbutton_displayref = NULL, *toggle_reg_manual1 = NULL, *toggle_reg_manual2 = NULL;
GtkWindow *control_window = NULL;

// additional statics
static GtkComboBox *filter_combo[3] = { NULL };
static GtkAdjustment *stackadj[3] = { NULL };
static GtkWidget *spin[3] = { NULL };
static GtkWidget *ksig[3] = { NULL };
static GtkLabel *filter_label[3] = { NULL };

static GList *switcher_buttons = NULL;

/****************************************************************/
/* Initialization                                               */
/****************************************************************/

static void registration_init_statics() {
	if (register_minpairs == NULL) {
		// GtkAdjustment
		register_minpairs = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "register_minpairs"));
		// GtkBox
		seq_filters_box_reg = GTK_BOX(gtk_builder_get_object(gui.builder, "seq_filters_box_reg"));
		reg_wcsfilechooser_box = GTK_BOX(gtk_builder_get_object(gui.builder, "reg_wcsfilechooser_box"));
		// GtkButton
		filter_add4 = GTK_BUTTON(gtk_builder_get_object(gui.builder, "filter_add4"));
		filter_add5 = GTK_BUTTON(gtk_builder_get_object(gui.builder, "filter_add5"));
		filter_rem5 = GTK_BUTTON(gtk_builder_get_object(gui.builder, "filter_rem5"));
		filter_rem6 = GTK_BUTTON(gtk_builder_get_object(gui.builder, "filter_rem6"));
		proj_estimate = GTK_BUTTON(gtk_builder_get_object(gui.builder, "proj_estimate"));
		goregister_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "goregister_button"));
		reg_wcsfile_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "reg_wcsfile_button"));
		// GtkComboBoxText
		comboboxregmethod = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboboxregmethod"));
		comboboxreglayer = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboboxreglayer"));
		comboreg_maxstars = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboreg_maxstars"));
		comboreg_transfo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboreg_transfo"));
		reg_sel_all_combobox = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "reg_sel_all_combobox"));
		combofilter4 = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "combofilter4"));
		filter_type4 = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "filter_type4"));
		combofilter5 = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "combofilter5"));
		filter_type5 = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "filter_type5"));
		combofilter6 = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "combofilter6"));
		filter_type6 = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "filter_type6"));
		comboreg_framing = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboreg_framing"));
		ComboBoxRegInter = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "ComboBoxRegInter"));
		combo_driz_kernel = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "combo_driz_kernel"));
		comboreg_undistort = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "comboreg_undistort"));
		// GtkEntry
		entry1_x_comet = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry1_x_comet"));
		entry2_x_comet = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry2_x_comet"));
		entry1_y_comet = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry1_y_comet"));
		entry2_y_comet = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry2_y_comet"));
		cometseqname_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "cometseqname_entry"));
		regseqname_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "regseqname_entry"));
		flatname_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "flatname_entry"));
		reg_wcsfile_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "reg_wcsfile_entry"));
		// GtkExpander
		autoreg_expander = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "autoreg_expander"));
		manualreg_expander = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "manualreg_expander"));
		// GtkFrame
		output_reg_frame = GTK_FRAME(gtk_builder_get_object(gui.builder, "output_reg_frame"));
		// GtkGrid
		grid_reg_framing = GTK_GRID(gtk_builder_get_object(gui.builder, "grid_reg_framing"));
		grid_interp_controls = GTK_GRID(gtk_builder_get_object(gui.builder, "grid_interp_controls"));
		grid_drizzle_controls = GTK_GRID(gtk_builder_get_object(gui.builder, "grid_drizzle_controls"));
		grid_reg_wcs = GTK_GRID(gtk_builder_get_object(gui.builder, "grid_reg_wcs"));
		// GtkImage
		framing_image = GTK_IMAGE(gtk_builder_get_object(gui.builder, "framing-image"));
		// GtkLabel
		label1_comet = GTK_LABEL(gtk_builder_get_object(gui.builder, "label1_comet"));
		regfilter_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "regfilter_label"));
		labelfilter4 = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter4"));
		labelfilter5 = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter5"));
		labelfilter6 = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter6"));
		labelregisterinfo = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelregisterinfo"));
		labelRegRef = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelRegRef"));
		estimate_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "estimate_label"));
		// GtkNotebook
		notebook_registration = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook_registration"));
		// GtkSpinButton
		spinbut_minpairs = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbut_minpairs"));
		spin_kombat_percent = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_kombat_percent"));
		stackspin4 = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stackspin4"));
		stackspin5 = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stackspin5"));
		stackspin6 = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stackspin6"));
		reg_scaling_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "reg_scaling_spin"));
		spin_driz_dropsize = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_driz_dropsize"));
		spinbut_shiftx = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbut_shiftx"));
		spinbut_shifty = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbut_shifty"));
		// GtkStack
		interp_drizzle_stack = GTK_STACK(gtk_builder_get_object(gui.builder, "interp_drizzle_stack"));
		// GtkStackSwitcher
		interp_drizzle_stack_switcher = GTK_STACK_SWITCHER(gtk_builder_get_object(gui.builder, "interp_drizzle_stack_switcher"));
		// GtkToggleButton
		checkStarSelect = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "checkStarSelect"));
		reg_2pass = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "reg_2pass"));
		followStarCheckButton = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "followStarCheckButton"));
		onlyshift_checkbutton = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "onlyshift_checkbutton"));
		toggle_reg_clamp = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "toggle_reg_clamp"));
		driz_use_flats = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "driz_use_flats"));
		checkbutton_displayref = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_displayref"));
		toggle_reg_manual1 = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "toggle_reg_manual1"));
		toggle_reg_manual2 = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "toggle_reg_manual2"));
		control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));

		// additional statics
		spin[0] = GTK_WIDGET(stackspin4);
		spin[1] = GTK_WIDGET(stackspin5);
		spin[2] = GTK_WIDGET(stackspin6);
		stackadj[0] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[0]));
		stackadj[1] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[1]));
		stackadj[2] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[2]));
		filter_combo[0] = GTK_COMBO_BOX(combofilter4);
		filter_combo[1] = GTK_COMBO_BOX(combofilter5);
		filter_combo[2] = GTK_COMBO_BOX(combofilter6);
		ksig[0] = GTK_WIDGET(filter_type4);
		ksig[1] = GTK_WIDGET(filter_type5);
		ksig[2] = GTK_WIDGET(filter_type6);
		filter_label[0] = labelfilter4;
		filter_label[1] = labelfilter5;
		filter_label[2] = labelfilter6;

		// switcher buttons list
		switcher_buttons = gtk_container_get_children(GTK_CONTAINER(interp_drizzle_stack_switcher));
	}
}

static void _reg_selected_area_callback() {
	if (!com.headless)
		update_reg_interface(TRUE);
}

void initialize_registration_methods() {
	int j = 0;
	GString *tip;
	gchar *ctip;
	registration_init_statics();

	reg_methods[REG_GLOBAL] = new_reg_method(_("Global Star Alignment (deep-sky)"),
			&register_star_alignment, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[REG_3STARS] = new_reg_method(_("1-2-3 Stars Registration (deep-sky)"),
			&register_3stars, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[REG_DFT] = new_reg_method(_("Image Pattern Alignment (planetary - full disk)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
	reg_methods[REG_KOMBAT] = new_reg_method(_("KOMBAT (planetary surfaces or full disk)"),
			&register_kombat, REQUIRES_ANY_SELECTION, REGTYPE_PLANETARY);
	reg_methods[REG_COMET] = new_reg_method(_("Comet/Asteroid Registration"),
			&register_comet, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[REG_APPLY] = new_reg_method(_("Apply Existing Registration"),
			&register_apply_reg, REQUIRES_NO_SELECTION, REGTYPE_APPLY);
	// we register 2-pass but we won't add it to the combo/tooltip
	reg_methods[REG_2PASS] = new_reg_method(_("Two-Pass Global Star Alignment (deep-sky)"),
			&register_multi_step_global, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[NUMBER_OF_METHODS] = NULL;

	tip = g_string_new ("");
	for (j = 0; j < NUMBER_OF_METHODS - 1; j ++) {
		tip = g_string_append(tip, _(tooltip_text[j]));
		if (j < NUMBER_OF_METHODS - 2)
			tip = g_string_append(tip, "\n\n");
	}
	ctip = g_string_free (tip, FALSE);
	gtk_widget_set_tooltip_markup(GTK_WIDGET(comboboxregmethod), ctip);
	g_free(ctip);

	/* fill comboboxregmethod */
	gtk_combo_box_text_remove_all(comboboxregmethod);
	for (j = 0; j < NUMBER_OF_METHODS - 1; j ++) {
		gtk_combo_box_text_append_text(comboboxregmethod, reg_methods[j]->name);
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(comboboxregmethod), com.pref.gui.reg_settings);

	gtk_combo_box_set_active(GTK_COMBO_BOX(ComboBoxRegInter), com.pref.gui.reg_interpolation);
	gtk_toggle_button_set_active(toggle_reg_clamp, com.pref.gui.reg_clamping);
	gtk_widget_set_sensitive(GTK_WIDGET(toggle_reg_clamp),
			com.pref.gui.reg_interpolation == OPENCV_LANCZOS4 || com.pref.gui.reg_interpolation == OPENCV_CUBIC);

	/* register to the new area selected event */
	register_selection_update_callback(_reg_selected_area_callback);
}

/****************************************************************/
/* GUI elements callbacks                                       */
/****************************************************************/

void on_reg_2pass_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	gboolean toggled = gtk_toggle_button_get_active(togglebutton);
	GtkWidget *outputframe = GTK_WIDGET(output_reg_frame);
	gtk_widget_set_sensitive(outputframe, !toggled);
	update_reg_interface(TRUE);
}

void on_regfollowStar_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	update_reg_interface(TRUE);
}

void on_comboboxregmethod_changed(GtkComboBox *box, gpointer user_data) {
	int index = 0;
	gchar *text = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(box));

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
	gtk_widget_set_sensitive(GTK_WIDGET(toggle_reg_clamp),
			com.pref.gui.reg_interpolation == OPENCV_LANCZOS4
					|| com.pref.gui.reg_interpolation == OPENCV_CUBIC);
}

void on_comboreg_transfo_changed(GtkComboBox *box, gpointer user_data) {
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

void on_comboreg_framing_changed(GtkComboBox *box, gpointer user_data) {
	gchar *name;
	int i = gtk_combo_box_get_active(box);

	if (i >= 0 && i < G_N_ELEMENTS(reg_frame_registration)) {
		name = g_strdup_printf("/org/siril/ui/pixmaps/%s", reg_frame_registration[i]);
		gtk_image_set_from_resource(framing_image, name);

		g_free(name);
	}
	// Re-check framing is OK for the seq type and reg method
	update_reg_interface(TRUE);
}

void on_comboreg_undistort_changed(GtkComboBox *box, gpointer user_data) {
	update_reg_interface(TRUE);
	gtk_widget_grab_focus(GTK_WIDGET(box));
}

void on_reg_wcsfile_button_clicked(GtkButton *button, gpointer user_data) {
		SirilWidget *widgetdialog;
		GtkFileChooser *dialog = NULL;
		GtkWindow *control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_local_only(dialog, FALSE);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		gtk_filter_add(dialog, _("WCS Files (*.wcs)"),
				"*.wcs", FALSE);
		gtk_filter_add(dialog, _("FITS Files (*.fit, *.fits, *.fts, *.fit.fz, *.fits.fz, *.fts.fz)"),
				"*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS;*.fit.fz;*.FIT.fz;*.fits.fz;*.FITS.fz;*.fts.fz;*.FTS.fz", gui.file_ext_filter == TYPEFITS);
		gint res = siril_dialog_run(widgetdialog);
		if (res == GTK_RESPONSE_ACCEPT) {
			gchar *file = siril_file_chooser_get_filename(dialog);
			gtk_entry_set_text(reg_wcsfile_entry, file);
			gtk_editable_set_position(GTK_EDITABLE(reg_wcsfile_entry), -1);
			g_free(file);
		}
		siril_widget_destroy(widgetdialog);
		update_reg_interface(TRUE);
}

gboolean on_switcher_stack_clicked(GtkWidget *widget,
	GdkEventButton *event, gpointer user_data) {
	update_reg_interface(TRUE);
	return TRUE;
}

/****************************************************************/
/* comet specific callbacks                                       */
/****************************************************************/

static void update_velocity() {
	velocity = compute_velocity(t_of_image_1, t_of_image_2, pos_of_image1, pos_of_image2);
	gchar *v_txt = g_strdup_printf("Δx: %.2lf, Δy: %.2lf", velocity.x, -velocity.y);
	gtk_label_set_text(label1_comet, v_txt);
	g_free(v_txt);
}

static void update_comet_entry(point pt, GtkEntry *entry_x, GtkEntry *entry_y) {
	gchar *txt_x, *txt_y;
	txt_x = g_strdup_printf("%7.2f", pt.x);
	txt_y = g_strdup_printf("%7.2f", pt.y);
	gtk_entry_set_text(entry_x, txt_x);
	gtk_entry_set_text(entry_y, txt_y);
	g_free(txt_x);
	g_free(txt_y);
}

void on_button_comet_clicked(GtkButton *button, gpointer p) {
	psf_star *result = NULL;
	int layer = gtk_combo_box_get_active(GTK_COMBO_BOX(comboboxreglayer));
	const gchar *caller = gtk_buildable_get_name(GTK_BUILDABLE(button));
	gboolean first = !g_strcmp0(caller, "button1_comet");

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		psf_error error = PSF_NO_ERR;
		result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE, FALSE, NULL, FALSE, com.pref.starfinder_conf.profile, &error);
		if (result && (result->x0 <= 0. || result->x0 >= com.selection.w || result->y0 <= 0. || result->x0 >= com.selection.h) && error != PSF_NO_ERR) { // we check result is inside the selection box
			siril_log_color_message(_("Comet PSF center is out of the box, will use selection center instead\n"), "salmon");
			free_psf(result);
			result = NULL;
		}
		point pos;
		if (result) {
			pos.x = result->x0 + com.selection.x;
			pos.y = com.selection.y + com.selection.h - result->y0;
		} else { // if psf fails we use the center of selection
			pos.x = 0.5 * (double)com.selection.w + com.selection.x;
			pos.y = com.selection.y + com.selection.h - 0.5 * (double)com.selection.h;
		}
		free_psf(result);
		if (layer_has_registration(&com.seq, layer) &&
				guess_transform_from_H(com.seq.regparam[layer][com.seq.reference_image].H) > NULL_TRANSFORMATION &&
				guess_transform_from_H(com.seq.regparam[layer][com.seq.current].H) > NULL_TRANSFORMATION) {
			cvTransfPoint(&pos.x, &pos.y, com.seq.regparam[layer][com.seq.current].H, com.seq.regparam[layer][com.seq.reference_image].H, 1.);
		}
		if (first) // comet1_clicked
			pos_of_image1 = pos;
		else
			pos_of_image2 = pos;
		if (!gfit.keywords.date_obs) {
			siril_message_dialog(GTK_MESSAGE_ERROR,
					_("There is no timestamp stored in the file"),
					_("Siril cannot perform the registration without date information in the file."));
		} else {
			GDateTime *time = g_date_time_ref(gfit.keywords.date_obs);
			if (!time) {
				siril_message_dialog(GTK_MESSAGE_ERROR,
						_("Unable to convert DATE-OBS to a valid date"),
						_("Siril cannot convert the DATE-OBS keyword into a valid date needed in the alignment."));
			}
			if (first) {
				if (t_of_image_1) {
					g_date_time_unref(t_of_image_1);
				}
				t_of_image_1 = time;
				update_comet_entry(pos_of_image1, entry1_x_comet, entry1_y_comet);
			} else {
				if (t_of_image_2) {
					g_date_time_unref(t_of_image_2);
				}
				t_of_image_2 = time;
				update_comet_entry(pos_of_image2, entry2_x_comet, entry2_y_comet);
			}
		}
		set_cursor_waiting(FALSE);
	}
	update_reg_interface(TRUE);
}

void on_entry_comet_changed(GtkEditable *editable, gpointer user_data) {
	pos_of_image1.x = g_ascii_strtod(gtk_entry_get_text(entry1_x_comet), NULL);
	pos_of_image1.y = g_ascii_strtod(gtk_entry_get_text(entry1_y_comet), NULL);
	pos_of_image2.x = g_ascii_strtod(gtk_entry_get_text(entry2_x_comet), NULL);
	pos_of_image2.y = g_ascii_strtod(gtk_entry_get_text(entry2_y_comet), NULL);
	update_velocity();
}

/****************************************************************/
/* Filters                                                      */
/****************************************************************/

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
	gtk_widget_set_visible(GTK_WIDGET(combofilter5), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(stackspin5), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(filter_add5), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(filter_rem5), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(labelfilter5), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(filter_type5), TRUE);
	update_filters_registration(-1);
}

void on_filter_add5_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(GTK_WIDGET(combofilter6), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(stackspin6), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(filter_rem6), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(labelfilter6), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(filter_type6), TRUE);
	update_filters_registration(-1);
}

void on_filter_rem5_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(GTK_WIDGET(combofilter5), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(stackspin5), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(filter_add5), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(filter_rem5), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(labelfilter5), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(filter_type5), FALSE);
	update_filters_registration(-1);
}

void on_filter_rem6_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(GTK_WIDGET(combofilter6), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(stackspin6), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(filter_rem6), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(labelfilter6), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(filter_type6), FALSE);
	update_filters_registration(-1);
}

static void update_filter_label(seq_image_filter filtering_criterion, double filtering_parameter) {
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
	gtk_label_set_text(regfilter_label, labelbuffer);
	g_free(labelbuffer);
}

static void get_reg_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter, int update_adjustment) {
	int filter, guifilter, channel = 0, type;
	gboolean is_ksig = FALSE;
	double percent = 0.0;

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

static void update_filters_registration(int update_adjustment) {
	if (!sequence_is_loaded())
		return;
	siril_debug_print("updating registration filters GUI\n");
	seq_image_filter criterion;
	double param;
	get_reg_sequence_filtering_from_gui(&criterion, &param, update_adjustment);
	update_filter_label(criterion, param);
}

/****************************************************************/
/* UI setup and collection                                      */
/****************************************************************/

// Checkers
static gboolean check_applyreg(regmethod_index index) {
	if (index != REG_APPLY)
		return TRUE;
	framing_type framingmethod = (framing_type)gtk_combo_box_get_active(GTK_COMBO_BOX(comboreg_framing));
	if (framingmethod == FRAMING_MAX && (com.seq.type == SEQ_FITSEQ || com.seq.type == SEQ_SER)) {
		gtk_label_set_text(labelregisterinfo, _("Max framing not allowed with FITSEQ or SER, change to regular FITS images"));
		return FALSE;
	}
	return TRUE;
}

static gboolean check_comet(regmethod_index index) {
	if (index != REG_COMET)
		return TRUE;
	if ((velocity.x == 0.0 && velocity.y == 0.0)
				|| isinf(velocity.x) || isinf(velocity.y)) {
		gtk_label_set_text(labelregisterinfo, _("Pick object in two images"));
		return FALSE;
	}
	return TRUE;
}

static gboolean check_3stars(regmethod_index index) {
	if (index != REG_3STARS)
		return TRUE;
	if (_3stars_check_selection()) {// checks that the right image is loaded based on doall and dofollow
		int nbselstars = _3stars_get_number_selected_stars();
		return nbselstars > 0;
	}
	return FALSE;
}

static gboolean check_disto(disto_source index) {
	gchar *label = NULL;
	gchar *tooltip = NULL;
	const gchar *text = gtk_entry_get_text(reg_wcsfile_entry);
	gboolean status =  validate_disto_params(&gfit, text, index, &tooltip, &label);
	if (!status) {
		gtk_label_set_text(labelregisterinfo, label);
		gtk_widget_set_tooltip_text(GTK_WIDGET(labelregisterinfo), tooltip);
		g_free(label);
		g_free(tooltip);
	}
	return status;
}
// Helpers
struct registration_method *get_selected_registration_method(int *index) {
	*index = REG_UNDEF;
	int ind = 0;
	gchar *text = gtk_combo_box_text_get_active_text(comboboxregmethod);
	while (reg_methods[ind] && text != NULL) {
		if (!strcmp(reg_methods[ind]->name, text))
			break;
		ind++;
	}
	g_free(text);
	if (ind == NUMBER_OF_METHODS)
		return NULL;

	if (ind == REG_GLOBAL && gtk_toggle_button_get_active(reg_2pass)) // 2 pass method is global + 2pass button checked
		ind = REG_2PASS;
	*index = ind;
	return reg_methods[ind];
}

/* Selects the "register all" or "register selected" according to the number of
 * selected images, if argument is false.
 * Verifies that enough images are selected and an area is selected.
 */
void update_reg_interface(gboolean dont_change_reg_radio) {
	int nb_images_reg; /* the number of images to register */
	struct registration_method *method = NULL;
	gboolean selection_is_done;
	gboolean has_reg, has_output, has_drizzle, must_have_drizzle, must_have_interp, samesizeseq_required, check_bayer_ok, ready, seqloaded;
	gboolean isapplyreg, is_global, is_star_align;
	sensor_pattern pattern = BAYER_FILTER_NONE;
	disto_source disto_source_index = DISTO_UNDEF;

	seqloaded = sequence_is_loaded();
	gtk_widget_set_sensitive(GTK_WIDGET(manualreg_expander), seqloaded);
	if (!seqloaded)
		gtk_expander_set_expanded(manualreg_expander, FALSE); // no need to clutter the interface, but we don't want to reset if user has expanded it
	gtk_widget_set_sensitive(GTK_WIDGET(autoreg_expander), seqloaded);
	gtk_expander_set_expanded(autoreg_expander, seqloaded);
	if (!seqloaded) // no need to go further, hide all and return
		return;

	if (!dont_change_reg_radio) {
		if (com.seq.selnum < com.seq.number) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(reg_sel_all_combobox), 1);
			gtk_combo_box_set_active(GTK_COMBO_BOX(combofilter4), 1);
		} else
			gtk_combo_box_set_active(GTK_COMBO_BOX(reg_sel_all_combobox), 0);
	}

	selection_is_done = (com.selection.w > 0 && com.selection.h > 0);

	/* number of registered image */
	nb_images_reg = gtk_combo_box_get_active(GTK_COMBO_BOX(reg_sel_all_combobox)) == 0 ? com.seq.number : com.seq.selnum;

	/* getting the selected registration method */
	regmethod_index regindex = REG_UNDEF;
	method = get_selected_registration_method(&regindex);
	if (!method) {
		siril_log_color_message(_("Failed to determine registration method...\n"), "red");
		return;
	}

	/* specific methods */
	isapplyreg = regindex == REG_APPLY;
	is_global = regindex == REG_GLOBAL;
	is_star_align = is_global || regindex == REG_2PASS;

	/* registration data exists for the selected layer */
	has_reg = layer_has_registration(&com.seq, gtk_combo_box_get_active(GTK_COMBO_BOX(comboboxreglayer)));
	/* registration will produce output */
	has_output = isapplyreg || is_global;
	/* drizzle is selected */
	has_drizzle = has_output && gtk_stack_get_visible_child(interp_drizzle_stack) == GTK_WIDGET(grid_drizzle_controls);
	/* must enforce drizzle/interp */
	must_have_drizzle = has_output && gfit.naxes[2] == 1 && gfit.keywords.bayer_pattern[0] != '\0';
	must_have_interp  = has_output && gfit.naxes[2] == 3;


	/* initialize default */
	gtk_notebook_set_current_page(notebook_registration, REG_PAGE_MISC);
	gtk_label_set_text(labelregisterinfo, "");
	gtk_widget_set_tooltip_text(GTK_WIDGET(labelregisterinfo), "");

	/* variable image sizes => disable manual reg */
	gtk_widget_set_sensitive(GTK_WIDGET(manualreg_expander), !com.seq.is_variable);
	if (com.seq.is_variable)
		gtk_expander_set_expanded(manualreg_expander, FALSE); // no need to clutter the interface
	gtk_widget_set_tooltip_text(GTK_WIDGET(manualreg_expander), (!com.seq.is_variable) ? "" : _("not available for sequences with variable image sizes"));

	/* show the appropriate frame selection widgets */
	gtk_widget_set_visible(GTK_WIDGET(reg_sel_all_combobox), !isapplyreg);
	gtk_widget_set_visible(GTK_WIDGET(seq_filters_box_reg), isapplyreg);
	gtk_widget_set_visible(GTK_WIDGET(combofilter4), isapplyreg);
	if (isapplyreg) {
		if (!dont_change_reg_radio && com.seq.selnum < com.seq.number) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(combofilter4), SELECTED_IMAGES);
		}
		update_filters_registration(-1);
	}
	/* update star_align widgets*/
	if (is_star_align) {
		gtk_widget_set_visible(GTK_WIDGET(grid_reg_wcs), !com.seq.is_variable);
		if (!com.seq.is_variable) {
			disto_source_index = gtk_combo_box_get_active(GTK_COMBO_BOX(comboreg_undistort));
			gtk_widget_set_visible(GTK_WIDGET(reg_wcsfilechooser_box), disto_source_index == DISTO_FILE);
		}
	}

	/* show the appropriate outputregframe widgets */
	gtk_widget_set_visible(GTK_WIDGET(output_reg_frame), isapplyreg || is_global);
	gtk_widget_set_sensitive(GTK_WIDGET(output_reg_frame), isapplyreg || is_global);
	gtk_widget_set_visible(GTK_WIDGET(proj_estimate), isapplyreg);
	gtk_widget_set_visible(GTK_WIDGET(notebook_registration), !isapplyreg);
	gtk_widget_set_visible(GTK_WIDGET(grid_reg_framing), isapplyreg);
	if (must_have_drizzle) {
		gtk_stack_set_visible_child(interp_drizzle_stack, GTK_WIDGET(grid_drizzle_controls));
		has_drizzle = TRUE;
	}
	if (must_have_interp) {
		gtk_stack_set_visible_child(interp_drizzle_stack, GTK_WIDGET(grid_interp_controls));
		has_drizzle = FALSE;
	}
	gtk_widget_set_sensitive(GTK_WIDGET(interp_drizzle_stack_switcher), !(must_have_drizzle || must_have_interp));
	set_switcher_buttons_colors(switcher_buttons, (has_drizzle) ? 1 : 0);
	if (has_output && !has_drizzle) {
		gint interpolation_item = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxRegInter));
		gtk_widget_set_sensitive(GTK_WIDGET(toggle_reg_clamp), interpolation_item == OPENCV_CUBIC || interpolation_item == OPENCV_LANCZOS4);
	}

	// show the relevant notebook page
	if (is_star_align) {
		gtk_notebook_set_current_page(notebook_registration, REG_PAGE_GLOBAL);
	} else if (regindex == REG_COMET) {
		gtk_notebook_set_current_page(notebook_registration, REG_PAGE_COMET);
	} else if (regindex == REG_3STARS) {
		gtk_notebook_set_current_page(notebook_registration, REG_PAGE_3_STARS);
	} else if (regindex == REG_KOMBAT) {
		gtk_notebook_set_current_page(notebook_registration, REG_PAGE_KOMBAT);
	}

	// if not debayered, check that the bayer pattern is known
	// at this stage it does not need to be correct, we just need to know there's one to identify if it's a valid bayer pattern
	// TODO: we can't force drizzle as we do for debayer because here, we rely on the header having a bayer_pattern keyword
	// This would need either to add a force bayer button in the registration tab (ugly)
	// or direct the user to correcting the bayer pattern in the images headers (prefered)
	if (gfit.naxes[2] == 1 && gfit.keywords.bayer_pattern[0] != '\0') {
		pattern = get_cfa_pattern_index_from_string(gfit.keywords.bayer_pattern);
	}

	// checking bayer status is ok
	check_bayer_ok = !has_output || //if no output, we don't need to check
					 (gfit.naxes[2] == 1 && gfit.keywords.bayer_pattern[0] == '\0') || // mono sequence
					 (gfit.naxes[2] == 1 && has_drizzle && pattern > BAYER_FILTER_MIN) || // bayer-drizzle sequence
					 gfit.naxes[2] == 3 || // debayered sequence
					 com.seq.type == SEQ_SER || // SER can be debayered on-the-fly
					 (has_drizzle && gfit.naxes[2] == 1); // drizzle or bayer-drizzle will be applied to produce the output

	// checking if it requires same size sequence
	samesizeseq_required = (regindex >= REG_3STARS && regindex <= REG_KOMBAT) || ((is_star_align || isapplyreg) && has_drizzle);

	// performing all checks
	ready = nb_images_reg > 1 &&
			check_bayer_ok &&
			(!samesizeseq_required || (samesizeseq_required && !com.seq.is_variable)) &&
			(selection_is_done || method->sel == REQUIRES_NO_SELECTION) &&
			(!isapplyreg || has_reg) && // must have reg data if applyreg
			check_applyreg(regindex) &&
			check_comet(regindex) &&
			check_3stars(regindex) &&
			check_disto(disto_source_index);

	if (!ready) { // all the other cases not set by the checkers
		if (method->sel > REQUIRES_NO_SELECTION && !selection_is_done ) {
			gtk_label_set_text(labelregisterinfo, _("Select an area in image first"));
		} else if (nb_images_reg <= 1) {
			gtk_label_set_text(labelregisterinfo, _("Select images in the sequence"));
		} else if (regindex == REG_APPLY && !has_reg){
			gtk_label_set_text(labelregisterinfo, _("Select a layer with existing registration"));
		} else if (samesizeseq_required && com.seq.is_variable) {
			gtk_label_set_text(labelregisterinfo, _("not available for sequences with variable image sizes"));
		} else if (pattern > BAYER_FILTER_MAX || pattern < BAYER_FILTER_MIN) {
			gtk_label_set_text(labelregisterinfo, _("Unsupported CFA pattern detected"));
			gtk_widget_set_tooltip_text(GTK_WIDGET(labelregisterinfo), _("This sequence cannot be registered with the CFA pattern intact. You must debayer it prior to registration"));
		} else if (!check_bayer_ok) {
			gtk_label_set_text(labelregisterinfo, _("CFA pattern detected"));
			gtk_widget_set_tooltip_text(GTK_WIDGET(labelregisterinfo), _("This sequence cannot be registered with the CFA pattern intact. Select drizzle or debayer it prior to registration"));
		}
	}
	gtk_label_set_text(estimate_label, "");
	gtk_widget_set_sensitive(GTK_WIDGET(goregister_button), ready);
	gtk_widget_set_sensitive(GTK_WIDGET(proj_estimate), ready);
}

// Functions to collect GUI data
static int populate_drizzle_data(struct driz_args_t *driz, sequence *seq) {
	driz->use_flats = gtk_toggle_button_get_active(driz_use_flats);
	driz->weight_scale = 1.f; // Not used for now
	driz->kernel = (enum e_kernel_t) gtk_combo_box_get_active(GTK_COMBO_BOX(combo_driz_kernel));
	driz->pixel_fraction = gtk_spin_button_get_value(spin_driz_dropsize);
	sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit.keywords.bayer_pattern); // no need to get the validated version
	driz->is_bayer = gfit.naxes[2] == 1 && pattern > BAYER_FILTER_MIN && pattern < BAYER_FILTER_MAX;
	int status;
	gchar *error = NULL;
	gchar *expression = NULL;
	if (driz->use_flats) {
		fits reffit = { 0 };
		if (seq_read_frame_metadata(seq, seq->reference_image, &reffit)) {
			siril_log_color_message(_("NOT USING FLAT: Could not load reference image\n"), "red");
			free(driz);
			clearfits(&reffit);
			return 1;
		}
		const gchar *flat_filename = gtk_entry_get_text(flatname_entry);
		expression = path_parse(&reffit, flat_filename, PATHPARSE_MODE_READ, &status);
		clearfits(&reffit);
		if (status) {
			error = _("NOT USING FLAT: could not parse the expression");
			free(driz);
			return 1;
		} else {
			if (expression[0] == '\0') {
				siril_log_message(_("Error: no master flat specified in the preprocessing tab.\n"));
				free(driz);
				return 1;
			} else {
				set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
				driz->flat = calloc(1, sizeof(fits));
				if (!readfits(expression, driz->flat, NULL, TRUE)) {
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

static int fill_registration_structure_from_GUI(struct registration_args *regargs) {
	struct registration_method *method;

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
	control_window_switch_to_tab(OUTPUT_LOGS);

	regmethod_index regindex = REG_UNDEF;
	method = get_selected_registration_method(&regindex);
	if (!method) {
		return 1;
	}

	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	regargs->layer = gtk_combo_box_get_active(GTK_COMBO_BOX(comboboxreglayer));
	get_the_registration_area(regargs, method);	// sets selection
	regargs->run_in_thread = TRUE;
	regargs->load_new_sequence = FALSE; // only TRUE for some methods. Will be updated in these cases

	/* specific methods */
	gboolean isapplyreg = regindex == REG_APPLY;
	gboolean is_global = regindex == REG_GLOBAL;
	gboolean is_star_align = is_global || regindex == REG_2PASS;
	/* registration will produce output */
	gboolean has_output_images = isapplyreg || is_global;
	gboolean has_drizzle = has_output_images && gtk_stack_get_visible_child(interp_drizzle_stack) == GTK_WIDGET(grid_drizzle_controls);

	regargs->func = method->method_ptr;
	regargs->seq = &com.seq;
	regargs->reference_image = sequence_find_refimage(&com.seq);
	regargs->no_output = !has_output_images && regindex != REG_COMET; // comet produces a new sequence with symlinks to previous images
	if (regindex == REG_3STARS) {
		regargs->follow_star = gtk_toggle_button_get_active(followStarCheckButton);
		regargs->type = (gtk_toggle_button_get_active(onlyshift_checkbutton)) ? SHIFT_TRANSFORMATION : SIMILARITY_TRANSFORMATION;
	}
	if (is_star_align) {
		regargs->min_pairs = gtk_spin_button_get_value_as_int(spinbut_minpairs);
		int starmaxactive = gtk_combo_box_get_active(GTK_COMBO_BOX(comboreg_maxstars));
		regargs->max_stars_candidates = (starmaxactive == -1) ? MAX_STARS_FITTED : maxstars_values[starmaxactive];
		regargs->type = gtk_combo_box_get_active(GTK_COMBO_BOX(comboreg_transfo));
		regargs->matchSelection = gtk_toggle_button_get_active(checkStarSelect);
		if (regargs->matchSelection && regargs->seq->is_variable) {
			siril_log_color_message(_("Cannot use area selection on a sequence with variable image sizes\n"), "red");
			return 1;
		}
		if (!regargs->matchSelection) {
			delete_selected_area(); // otherwise it is enforced
		}
	}
	if (regindex == REG_KOMBAT) {
		regargs->percent_moved = (float)gtk_spin_button_get_value(spin_kombat_percent) / 100.f;
	}
	if (regindex == REG_COMET) {
		regargs->velocity = velocity;
		regargs->prefix = g_strdup(gtk_entry_get_text(cometseqname_entry)); //to create the .seq file
	}
	if (regindex == REG_APPLY) {
		regargs->framing = gtk_combo_box_get_active(GTK_COMBO_BOX(comboreg_framing));
		get_reg_sequence_filtering_from_gui(
				&regargs->filtering_criterion, &regargs->filtering_parameter, -1);
	} else {
		regargs->filters.filter_included = gtk_combo_box_get_active(GTK_COMBO_BOX(reg_sel_all_combobox));
	}

	if (has_output_images) {
		regargs->prefix = strdup(gtk_entry_get_text(regseqname_entry));
		regargs->output_scale = (float)gtk_spin_button_get_value(reg_scaling_spin);
		if (has_drizzle) {
			regargs->driz = calloc(1, sizeof(struct driz_args_t));
			if (populate_drizzle_data(regargs->driz, regargs->seq)) {
				return 1;
			}
		} else { // interpolation
			regargs->interpolation = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxRegInter));
			regargs->clamp = gtk_toggle_button_get_active(toggle_reg_clamp) && (regargs->interpolation == OPENCV_CUBIC || regargs->interpolation == OPENCV_LANCZOS4);
		}
	}

	regargs->undistort = DISTO_UNDEF;
	if (is_star_align) {
		regargs->undistort = gtk_combo_box_get_active(GTK_COMBO_BOX(comboreg_undistort));
		if (regargs->undistort > DISTO_UNDEF) {
			regargs->distoparam.index = regargs->undistort;
			if (regargs->undistort == DISTO_FILE) {
				regargs->distoparam.filename = g_strdup(gtk_entry_get_text(reg_wcsfile_entry));
			}
		}
	}
	if (isapplyreg && seq_has_any_distortion(regargs->seq)) {
		regargs->undistort = regargs->seq->distoparam[regargs->layer].index;
		regargs->distoparam = regargs->seq->distoparam[regargs->layer];
	}

	//Checks

#ifndef HAVE_CV44
	if (regargs->type == SHIFT_TRANSFORMATION && is_star_align) {
		siril_log_color_message(_("Shift-only registration is only possible with OpenCV 4.4\n"), "red");
		free(regargs->prefix);
		return 1;
	}
#endif

	/* We check that available disk space is enough when
	the registration method produces a new sequence with images
	*/
	if (has_output_images) {
		int nb_frames = regargs->filters.filter_included ? regargs->seq->selnum : regargs->seq->number;
		gint64 size = seq_compute_size(regargs->seq, nb_frames, get_data_type(regargs->seq->bitpix));
		if (regargs->output_scale != 1.f)
			size = (int64_t)(regargs->output_scale * regargs->output_scale * (float)size);
		if (test_available_space(size)) {
			siril_log_color_message(_("Not enough space to save the output images, aborting\n"), "red");
			return 1;
		}
	}

	if (regindex == REG_GLOBAL && regargs->interpolation == OPENCV_NONE) { // seqpplyreg case is dealt with in the sanity checks of the method
		if (regargs->output_scale != 1.f || com.seq.is_variable) {
			siril_log_color_message(_("When interpolation is set to None, the images must be of same size and no scaling can be applied. Aborting\n"), "red");
			return 1;
		}
		if (regargs->type > SHIFT_TRANSFORMATION) {
			siril_log_color_message(_("When interpolation is set to None, the transformation can only be set to Shift. Aborting\n"), "red");
			return 1;
		}
		if (regargs->undistort) {
			siril_log_color_message(_("When interpolation is set to None, distortions must be set to None as well. Aborting\n"), "red");
			return 1;
		}
	}

	if (regindex != REG_3STARS)
		clear_stars_list(TRUE); //to avoid problems with com.stars later on in the process

	return 0;
}

int get_registration_layer_from_GUI(const sequence *seq) {
	int reglayer = gtk_combo_box_get_active(GTK_COMBO_BOX(comboboxreglayer));
	if (!seq || !seq->regparam || !seq->regparam[reglayer] || seq->nb_layers < 0 || seq->nb_layers <= reglayer)
		return -1;
	return reglayer;
}

/* callback for 'Go register' button, GTK thread */
void on_seqregister_button_clicked(GtkButton *button, gpointer user_data) {

	struct registration_args *regargs = calloc(1, sizeof(struct registration_args));

	char *msg;
	if (fill_registration_structure_from_GUI(regargs)) {
		free(regargs);
		unreserve_thread();
		return;
	}
	fits fit_ref = { 0 };
	int ret = seq_read_frame_metadata(regargs->seq, regargs->reference_image, &fit_ref);
	if (ret) {
		siril_log_message(_("Error: unable to read reference frame metadata\n"));
		free(regargs);
		unreserve_thread();
		return;
	}
	regargs->bayer = (fit_ref.keywords.bayer_pattern[0] != '\0');
	clearfits(&fit_ref);

	regmethod_index regindex = REG_UNDEF;
	struct registration_method *method = get_selected_registration_method(&regindex);

	const gchar *caller = gtk_buildable_get_name(GTK_BUILDABLE(button));
	if (!g_strcmp0(caller, "proj_estimate"))
		regargs->no_output = TRUE;

	msg = siril_log_color_message(_("Registration: processing using method: %s\n"),
			"green", method->name);
	msg[strlen(msg) - 1] = '\0';

	if (regargs->clamp)
		siril_log_message(_("Interpolation clamping active\n"));
	set_progress_bar_data(msg, PROGRESS_RESET);

	if (!start_in_reserved_thread(register_thread_func, regargs)) {
		free(regargs);
		unreserve_thread();
	}
}

// end of registration, GTK thread. Executed when started from the GUI and in
// the graphical command line but not from a script (headless mode)
gboolean end_register_idle(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	stop_processing_thread();

	if (!args->retval) {
		if (!args->load_new_sequence && sequence_is_loaded()) {
			int chan = gtk_combo_box_get_active(GTK_COMBO_BOX(comboboxreglayer));
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

	if (args->func == &register_apply_reg && args->no_output) { // Estimate was pressed, we need to update the label
		gchar *downscale = (args->output_scale != 1.f) ? g_strdup_printf(_(" (assuming a scaling factor of %.2f)"), args->output_scale) : g_strdup("");
		gchar *scalemsg = g_strdup_printf(_("Output image: %d x %d pixels%s\n"), args->framingd.roi_out.w, args->framingd.roi_out.h, downscale);
		gtk_label_set_text(estimate_label, scalemsg);
		g_free(downscale);
		g_free(scalemsg);
		if (!args->retval)
			control_window_switch_to_tab(REGISTRATION); // if there are some warnings we stay on the Console tab
	}

	free(args->new_seq_name);
	free(args);
	return FALSE;
}
