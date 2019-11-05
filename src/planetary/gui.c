/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include <float.h>
#include <assert.h>
#include "core/siril.h"
#include "core/proto.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "stacking/stacking.h"
#include "planetary.h"
#include "gui.h"
#include "algos/statistics.h"
#include "registration/registration.h"
#include "gui/gui.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/zones.h"
#include "gui/plot.h"

static gboolean seqimage_range_mouse_pressed = FALSE;
static double lowest_accepted_quality = 0.0;
static double percent_images_to_keep = 100.0;

static double get_overlapamout() {
	double val;
	GtkAdjustment *overlapadj = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment_overlap"));
	val = gtk_adjustment_get_value(overlapadj);

	return val;
}

/* the planetary mode uses a slider to change displayed image whereas the
 * deep-sky mode uses a spin button (on_imagenumberspin_output) */
void on_seqimage_changed(GtkRange *range, gpointer user_data) {
	if (!sequence_is_loaded() || seqimage_range_mouse_pressed) return;
	int index = (int)gtk_range_get_value(range);
	seq_load_image(&com.seq, index, TRUE);
}

gboolean on_seqimage_button_press(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	seqimage_range_mouse_pressed = TRUE;
	return FALSE;
}

gboolean on_seqimage_button_release(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	seqimage_range_mouse_pressed = FALSE;
	on_seqimage_changed(GTK_RANGE(widget), NULL);
	return FALSE;
}

void on_bestimage_changed(GtkRange *range, gpointer user_data) {
	if (!sequence_is_loaded()) return;
	int layer = get_registration_layer(&com.seq);
	percent_images_to_keep = gtk_range_get_value(range);
	lowest_accepted_quality = compute_lowest_accepted_quality(&com.seq, layer,
			percent_images_to_keep);
	plot_set_filtering_threshold(lowest_accepted_quality);
}

void activate_mpp_processing_button() {
	static GtkWidget *stack_button = NULL;
	int status = sequence_is_loaded() && refimage_is_set() &&
		com.stacking_zones && com.stacking_zones[0].centre.x >= 0;
	// we should also check that there are at least two images to stack
	if (!stack_button)
	       stack_button = lookup_widget("gostack_button");
	gtk_widget_set_sensitive(stack_button, status);
}

// add the new zone index in the graphical list
void update_zones_list() {
	static GtkComboBox *combo = NULL;
	static GtkComboBoxText *combotext = NULL;
	if (!combo) {
		combo = GTK_COMBO_BOX(lookup_widget("selectedzonecombo"));
		combotext = GTK_COMBO_BOX_TEXT(combo);
	}
	if (!com.stacking_zones) return;
	int previous_selection = gtk_combo_box_get_active(combo);

	gtk_combo_box_text_remove_all(combotext);
	gtk_combo_box_text_append_text(combotext, "global");

	int i = 0;
	char buf[30];
	while (get_zone_available(i)) {
		i++;
		snprintf(buf, 30, "%s %d", _("zone"), i);
		gtk_combo_box_text_append_text(combotext, buf);
	}
	if (i > previous_selection)
		gtk_combo_box_set_active(combo, previous_selection);
	else 	gtk_combo_box_set_active(combo, 0);	// global
}

gboolean end_mpp_analysis(gpointer arg) {
	update_zones_list();
	return end_generic(arg);
}

/* First step: running the global registration and building the reference image */
void on_planetary_analysis_clicked(GtkButton *button, gpointer user_data) {
	if (!sequence_is_loaded()) return;
	GtkComboBox *cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));
	GtkSpinButton *kernelsz = GTK_SPIN_BUTTON(
			gtk_builder_get_object(builder, "kernelspin"));
	GtkToggleButton *save_precomp = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "store_precomp_checkbutton"));

	struct registration_args *reg_args;
	reg_args = calloc(1, sizeof(struct registration_args));
	reg_args->func = get_selected_registration_method()->method_ptr;
	reg_args->seq = &com.seq;
	reg_args->reference_image = sequence_find_refimage(&com.seq);
	reg_args->process_all_frames = FALSE;
	reg_args->follow_star = FALSE;
	reg_args->match_selection = FALSE;
	reg_args->translation_only = TRUE;
	reg_args->x2upscale = FALSE;
	reg_args->kernel_size = gtk_spin_button_get_value_as_int(kernelsz);
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE;
	reg_args->use_caching = gtk_toggle_button_get_active(save_precomp);

	char *msg = siril_log_color_message(_("Starting sequence analysis\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	set_cursor_waiting(TRUE);

	start_in_new_thread(sequence_analysis_thread_func, reg_args);
}

/* clicking on registration, only start registration */
gboolean on_mpreg_button_clicked(GtkButton *button, gpointer user_data) {
	GtkComboBox *cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));
	GtkSpinButton *kernelsz = GTK_SPIN_BUTTON(
			gtk_builder_get_object(builder, "kernelspin"));
	GtkToggleButton *save_precomp = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "store_precomp_checkbutton"));
	struct mpr_args *args = malloc(sizeof(struct mpr_args));
	args->seq = &com.seq;
	args->layer = gtk_combo_box_get_active(cbbt_layers);
	args->filtering_criterion = seq_filter_included;
	args->use_caching = gtk_toggle_button_get_active(save_precomp);
	args->kernel_size = gtk_spin_button_get_value_as_int(kernelsz);
	start_in_new_thread(the_multipoint_analysis, args);
	return FALSE;
}

/* clicking on stacking, start everything */
gboolean on_planetary_processing_button_clicked(GtkButton *button, gpointer user_data) {
	GtkComboBox *cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));
	GtkEntry *output_file = GTK_ENTRY(gtk_builder_get_object(builder, "entryresultfile"));
	GtkToggleButton *overwrite = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "checkbutoverwrite"));
	GtkSpinButton *kernelsz = GTK_SPIN_BUTTON(
			gtk_builder_get_object(builder, "kernelspin"));
	GtkToggleButton *save_precomp = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "store_precomp_checkbutton"));
	GtkComboBox *ap_method = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "ap_method_combo"));

	struct mpr_args *args = malloc(sizeof(struct mpr_args));
	args->seq = &com.seq;
	args->layer = gtk_combo_box_get_active(cbbt_layers);
	args->filtering_criterion = seq_filter_quality;
	args->filtering_parameter = lowest_accepted_quality;
	args->filtering_percent = percent_images_to_keep;
	args->ap_method = gtk_combo_box_get_active(ap_method);
	//args->using_homography = FALSE;	// TRUE is not working yet, translation only
	args->nb_closest_AP = 5;	// min(this, nb_AP) will be used
	args->max_distance = 250.0;	// AP farther than this will be ignored
	//args->own_distance_f = 0.5;	// unused
	args->output_filename = strdup(gtk_entry_get_text(output_file));
	args->output_overwrite = gtk_toggle_button_get_active(overwrite);
	args->global_image = NULL;
	args->use_caching = gtk_toggle_button_get_active(save_precomp);
	args->kernel_size = gtk_spin_button_get_value_as_int(kernelsz);
	start_in_new_thread(the_multipoint_processing, args);
	return FALSE;
}

/* when the reglayer is changed we update the status of the reference image */
void on_comboboxreglayer_planetary_changed(GtkComboBox *widget, gpointer user_data) {
	update_refimage_on_layer_change(&com.seq, gtk_combo_box_get_active(widget));
	activate_mpp_processing_button();
	display_refimage_if_needed();
}

void on_showref_check_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (refimage_is_set() && gtk_toggle_button_get_active(togglebutton)) {
		if (single_image_is_loaded()) {
		       if (!strcmp(com.uniq->filename, get_refimage_filename()))
			       return;
		       close_single_image();
		}
		clearfits(&gfit);
		const fits *refimage = get_refimage();
		copyfits(refimage, &gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		// code from end_stacking
		com.uniq = calloc(1, sizeof(single));
		com.uniq->comment = strdup(_("Stacking result image"));
		com.uniq->nb_layers = gfit.naxes[2];
		com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
		com.uniq->fit = &gfit;
		com.uniq->filename = strdup(get_refimage_filename());

		com.seq.current = RESULT_IMAGE;
		adjust_cutoff_from_updated_gfit();
		set_sliders_value_to_gfit();
		initialize_display_mode();

		sliders_mode_set_state(com.sliders);
		set_cutoff_sliders_max_values();

		set_display_mode();
		display_filename();
	} else {
		if (!single_image_is_loaded() || !sequence_is_loaded()) return;
		if (com.seq.current != RESULT_IMAGE) return;
		// close_single_image does not display a sequence image, it had
		// to be modified as below
		if (com.uniq->filename)
			free(com.uniq->filename);
		if (com.uniq->comment)
			free(com.uniq->comment);
		free(com.uniq);
		com.uniq = NULL;

		int image_to_load = sequence_find_refimage(&com.seq);
		seq_load_image(&com.seq, image_to_load, TRUE);
	}
	redraw(com.cvport, REMAP_ALL);
}

void display_refimage_if_needed() {
	on_showref_check_toggled(GTK_TOGGLE_BUTTON(lookup_widget("showref_check")), NULL);
}

static WORD Compute_threshold(fits *fit, int layer, WORD *norm) {
	WORD threshold;
	imstats *stat;

	stat = statistics(NULL, -1, fit, layer, NULL, STATS_BASIC);
	if (!stat) {
		//siril_log_message(_("Error: statistics computation failed.\n"));
		return 0;
	}
	threshold = (WORD) stat->median + 1 * (WORD) stat->sigma;
	*norm = (WORD) stat->normValue;
	free_stats(stat);

	return threshold;
}

static void auto_add_stacking_zones(fits *fit, int layer, int size) {
	WORD **image;
	WORD threshold, norm;
	int nx = fit->rx;
	int ny = fit->ry;
	int k, y, nbzone = 0;
	int areaX0 = 0;
	int areaY0 = 0;
	int areaX1 = nx;
	int areaY1 = ny;
	// TODO wich value should we take in radius ??
	int radius = 50;
	double overlap = get_overlapamout();

	threshold = Compute_threshold(fit, layer, &norm);

	/* FILL real image upside-down */
	image = malloc(ny * sizeof(WORD *));
	if (image == NULL) {
		fprintf(stderr, "Memory allocation failed\n");
		return;
	}
	for (k = 0; k < ny; k++)
		image[ny - k - 1] = fit->pdata[layer] + k * nx;

	for (y = radius + areaY0; y < areaY1 - radius; y++) {
		int x;
		for (x = radius + areaX0; x < areaX1 - radius; x++) {
			WORD pixel = image[y][x];
			if (pixel > threshold/* && pixel < norm*/) {
				int yy, xx;
				gboolean bingo = TRUE;
				WORD neighbor;
				for (yy = y - 1; yy <= y + 1; yy++) {
					for (xx = x - 1; xx <= x + 1; xx++) {
						if (xx == x && yy == y)
							continue;
						neighbor = image[yy][xx];
						if (neighbor > pixel) {
							bingo = FALSE;
							break;
						} else if (neighbor == pixel) {
							if ((xx <= x && yy <= y) || (xx > x && yy < y)) {
								bingo = FALSE;
								break;
							}
						}
					}
				}
				if (bingo) {
					nbzone++;
					if (!zone_is_too_close(x, y, overlap))
						add_stacking_zone(x, y, size);
				}
			}
		}
	}
	if (nbzone) {
		update_zones_list();
		redraw(com.cvport, REMAP_ALL);
	}
	free(image);
}

void on_autoposition_button_clicked(GtkButton *button, gpointer user_data) {
	GtkAdjustment *sizeadj = GTK_ADJUSTMENT(
			gtk_builder_get_object(builder, "adjustment_zonesize"));
	int half_size, size = round_to_int(gtk_adjustment_get_value(sizeadj));
	assert(size > 0);
	half_size = (size - 1) / 2;

	remove_stacking_zones();

	// TODO: set layer
	auto_add_stacking_zones(&gfit, 0, half_size);
}

void on_check_autopos_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {

	GtkWidget *overlapamount = lookup_widget("overlapamount");
	GtkWidget *autoposition_button = lookup_widget("autoposition_button");
	gboolean is_active = gtk_toggle_button_get_active(togglebutton);

	gtk_widget_set_sensitive(overlapamount, is_active);
	gtk_widget_set_sensitive(autoposition_button, is_active);
}

void on_selectedzonecombo_changed(GtkComboBox *widget, gpointer user_data) {
	int zone_id = gtk_combo_box_get_active(widget);
	com.stacking_zone_focus = -1;

	if (zone_id < 0) return;
	if (zone_id == 0) {
		drawPlot();
		redraw(com.cvport, REMAP_NONE);		// un-highlight prev selected zone
		return;
	}

	// show graph for a zone and highlight it
	int i = 0;
	while (com.stacking_zones[i].centre.x >= 0 && i < com.stacking_zones_size - 1) {
		i++;
		if (i == zone_id) {
			com.stacking_zone_focus = i-1;
			redraw(com.cvport, REMAP_NONE);		// show the new selected zone
			//fprintf(stdout, "setting focus to zone %d\n", i);
			drawPlot();
		}
	}
}

