/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/remixer.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/siril_preview.h"
#include "gui/histogram.h"
#include "core/undo.h"
#include "core/arithm.h"
#include "core/siril_app_dirs.h"
#include "filters/ght.h"

#include "gui/histogram.h"
#include <gsl/gsl_histogram.h>

// Invocation: 1 if called directly from starnet GUI,
// 2 if called from Siril menu
static int invocation = 0;

static gdouble leftD = 0.0, rightD = 0.0;
static gdouble leftB = 0.0, rightB = 0.0;
static gdouble leftLP = 0.0, rightLP = 0.0;
static gdouble leftSP = 0.0, rightSP = 0.0;
static gdouble leftHP = 1.0, rightHP = 1.0;
static gdouble leftBP = 0.0, rightBP = 0.0;
static gdouble mastermixer = 0.5;
static double finalstretch = 0.0;
static int type_left = STRETCH_PAYNE_NORMAL, type_right = STRETCH_PAYNE_NORMAL;
static int colour_left = COL_INDEP, colour_right = COL_INDEP;
static fits fit_left;
static fits fit_right;
static fits fit_left_calc;
static fits fit_right_calc;

static ght_params params_left, params_right, params_histo_left, params_histo_right;
static ght_compute_params cp_histo_left = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static ght_compute_params cp_histo_right = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

static gboolean left_changed = FALSE, right_changed = FALSE;
static gboolean leftBP_changed = FALSE, rightBP_changed = FALSE;
static gchar* filename_left;
static gchar* filename_right;

static gboolean advanced_interface = FALSE;
static gboolean permit_calculation = FALSE;
static gboolean remixer_show_preview;
static gboolean left_loaded = FALSE;
static gboolean right_loaded =FALSE;

////////////////////////////////////////////
// Remixer histogram functionality        //
// Some callouts to histogram.h           //
// Many functions adaptedfrom histogram.c //
// to support 2 sets of histograms        //
////////////////////////////////////////////

static gsl_histogram *remix_histlayers_left[MAXVPORT];
static gsl_histogram *remix_histlayers_right[MAXVPORT];
static gsl_histogram *remix_histlayers_backup_left[MAXVPORT];
static gsl_histogram *remix_histlayers_backup_right[MAXVPORT];

static double histo_color_r[] = { 1.0, 0.0, 0.0, 0.0 };
static double histo_color_g[] = { 0.0, 1.0, 0.0, 0.0 };
static double histo_color_b[] = { 0.0, 0.0, 1.0, 0.0 };

// Applies stretch to histogram
static void apply_remix_histo(gsl_histogram *histo, float norm,
		ght_params *params, ght_compute_params *cp) {

	size_t int_norm = (size_t)norm;
	gsl_histogram *mtf_histo = gsl_histogram_alloc(int_norm + 1);

	gsl_histogram_set_ranges_uniform(mtf_histo, 0, norm);

#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
#endif
	for (size_t i = 0; i < int_norm + 1; i++) {
		float ght;
		float binval = gsl_histogram_get(histo, i);
		float pxl = ((float)i / norm);

		ght = GHT(pxl, params->B, params->D, params->LP, params->SP, params->HP, params->BP, params->stretchtype, cp) * norm;
		gsl_histogram_accumulate(mtf_histo, (double) ght, (double) binval);
	}
#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp critical
#endif
#endif
	gsl_histogram_memcpy(histo, mtf_histo);
	gsl_histogram_free(mtf_histo);
}

static void draw_remix_curve(cairo_t *cr, int width, int height, ght_params *params, ght_compute_params *cp) {
	// draw curve
	int k;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, 0.98, 0.5, 0.45);
	GHTsetup(cp, params->B, params->D, params->LP, params->SP, params->HP, params->stretchtype);
	for (k = 0; k < width + 1; k++) {
		float x = k / (float) width;
		float y = (float) GHT((double) x, params->B, params->D, params->LP, params->SP, params->HP, 0.0, params->stretchtype, cp);
		cairo_line_to(cr, k, height * (1 - y));
	}
	cairo_stroke(cr);
}

void display_remix_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height, double zoomH, double zoomV, gboolean isOrig) {
	if (width <= 0) return;
	int current_bin;
	size_t norm = gsl_histogram_bins(histo) - 1;
	gboolean remix_log_scale = TRUE;

	float vals_per_px = (float)norm / (float)width;	// size of a bin
	size_t i, nb_orig_bins = gsl_histogram_bins(histo);

	// We need to store the binned histogram in order to find the binned maximum
	static gfloat *displayed_values = NULL;
	static int nb_bins_allocated = 0;
	/* we create a bin for each pixel in the displayed width.
	 * nb_bins_allocated is thus equal to the width of the image */
	if (nb_bins_allocated != width) {
		gfloat *tmp;
		nb_bins_allocated = width;
		tmp = realloc(displayed_values, nb_bins_allocated * sizeof(gfloat));
		if (!tmp) {
			if (displayed_values != NULL) {
				g_free(displayed_values);
				displayed_values = NULL;
			}
			PRINT_ALLOC_ERR;
//			histo_close(TRUE);
			return;
		}
		displayed_values = tmp;
		memset(displayed_values, 0, nb_bins_allocated);
	}
	if (gfit.naxis == 2)
		cairo_set_source_rgb(cr, 255.0, 255.0, 255.0);
	else
		cairo_set_source_rgb(cr, histo_color_r[layer], histo_color_g[layer],
				histo_color_b[layer]);
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.5);
	if (isOrig)
		cairo_set_line_width(cr, 0.5);

	// first loop builds the bins and finds the maximum
	i = 0;
	float graph_height = 0.f;
	current_bin = 0;
	do {
		float bin_val = 0.f;
		while (i < nb_orig_bins
				&& (float)i / vals_per_px <= (float)current_bin + 0.5f) {
			bin_val += (float)gsl_histogram_get(histo, i);
			i++;
		}
		if (remix_log_scale && bin_val != 0.f) {
			bin_val = logf(bin_val);
		}
		displayed_values[current_bin] = bin_val;
		if (bin_val > graph_height)	// check for maximum
			graph_height = bin_val;
		current_bin++;
	} while (i < nb_orig_bins && current_bin < nb_bins_allocated);
	for (i = 0; i < nb_bins_allocated; i++) {
		float bin_height = height - height * displayed_values[i] / graph_height;
		cairo_line_to(cr, i, bin_height);
	}
	cairo_stroke(cr);
}

gboolean redraw_remix_histo_left(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width, height;

	width = gtk_widget_get_allocated_width(lookup_widget("remix_histo_left"));
	height = gtk_widget_get_allocated_height(lookup_widget("remix_histo_left"));

	if (height == 1)
		return FALSE;
	erase_histo_display(cr, width, height);
	draw_grid(cr, width, height);
		for (i = 0; i < MAXVPORT; i++) {
		if (remix_histlayers_left[i]) {
			display_remix_histo(remix_histlayers_backup_left[i], cr, i, width, height, 1.0, 1.0, TRUE);
			display_remix_histo(remix_histlayers_left[i], cr, i, width, height, 1.0, 1.0, FALSE);
		}
	draw_remix_curve(cr, width, height, &params_histo_left, &cp_histo_left);
	}
	return TRUE;
}

gboolean redraw_remix_histo_right(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width, height;

	width = gtk_widget_get_allocated_width(lookup_widget("remix_histo_right"));
	height = gtk_widget_get_allocated_height(lookup_widget("remix_histo_right"));

	if (height == 1)
		return FALSE;
	erase_histo_display(cr, width, height);
	draw_grid(cr, width, height);

	for (i = 0; i < MAXVPORT; i++) {
		if (remix_histlayers_right[i]) {
			display_remix_histo(remix_histlayers_backup_right[i], cr, i, width, height, 1.0, 1.0, TRUE);
			display_remix_histo(remix_histlayers_right[i], cr, i, width, height, 1.0, 1.0, FALSE);
		}
		draw_remix_curve(cr, width, height, &params_histo_right, &cp_histo_right);
	}
	return TRUE;
}

static void update_remix_histo_left() {
	if (!remix_histlayers_left[0]) return;
	float norm = (float)gsl_histogram_bins(remix_histlayers_left[0]) - 1;
	params_histo_left = (ght_params) { leftB, leftD, leftLP, leftSP, leftHP, leftBP, type_left, colour_left, TRUE, TRUE, TRUE };

	GHTsetup(&cp_histo_left, params_histo_left.B, params_histo_left.D, params_histo_left.LP, params_histo_left.SP, params_histo_left.HP, params_histo_left.stretchtype);
	for (size_t i = 0; i < fit_left.naxes[2]; i++) {
		gsl_histogram_memcpy(remix_histlayers_left[i], remix_histlayers_backup_left[i]);
		apply_remix_histo(remix_histlayers_left[i], norm, &params_histo_left, &cp_histo_left);
	}
	if (leftBP_changed) {
		params_histo_left.stretchtype = STRETCH_LINEAR;
		for (size_t i = 0; i < fit_left.naxes[2]; i++) {
			apply_remix_histo(remix_histlayers_left[i], norm, &params_histo_left, &cp_histo_left);
			}
	}
	gtk_widget_queue_draw(lookup_widget("remix_histo_left"));
}

static void update_remix_histo_right() {
	if (!remix_histlayers_right[0]) return;
	float norm = (float)gsl_histogram_bins(remix_histlayers_right[0]) - 1;
	params_histo_right = (ght_params) { rightB, rightD, rightLP, rightSP, rightHP, rightBP, type_right, colour_right, TRUE, TRUE, TRUE };

	int nlayers = fit_right.naxes[2];
	GHTsetup(&cp_histo_right, params_histo_right.B, params_histo_right.D, params_histo_right.LP, params_histo_right.SP, params_histo_right.HP, params_histo_right.stretchtype);
	for (size_t i = 0; i < nlayers; i++) {
		gsl_histogram_memcpy(remix_histlayers_right[i], remix_histlayers_backup_right[i]);
		apply_remix_histo(remix_histlayers_right[i], norm, &params_histo_right, &cp_histo_right);
	}
	if (rightBP_changed) {
		params_histo_right.stretchtype = STRETCH_LINEAR;
		for (size_t i = 0; i < fit_right.naxes[2]; i++) {
			apply_remix_histo(remix_histlayers_right[i], norm, &params_histo_right, &cp_histo_right);
		}
	}
	gtk_widget_queue_draw(lookup_widget("remix_histo_right"));
}

static void set_remix_histogram(gsl_histogram *histo, int layer, int side) {
	g_assert(layer >= 0 && layer < MAXVPORT);
	if (side == 0) {
		if (remix_histlayers_left[layer])
			gsl_histogram_free(remix_histlayers_left[layer]);
		remix_histlayers_left[layer] = histo;
	} else if (side == 1) {
		if (remix_histlayers_right[layer])
			gsl_histogram_free(remix_histlayers_right[layer]);
		remix_histlayers_right[layer] = histo;
	}
}

void compute_remix_histo_for_fit(fits *fit, int side) {
	int nb_layers = 3;
	if (fit->naxis == 2)
		nb_layers = 1;
	for (int i = 0; i < nb_layers; i++) {
		if (side == 0) {
			if (!remix_histlayers_left[i])
				set_remix_histogram(computeHisto(fit, i), i, side);
		} else if (side == 1) {
			if (!remix_histlayers_right[i])
				set_remix_histogram(computeHisto(fit, i), i, side);
		}
		else return;
	}
}

static void remix_histo_startup_left() {
	compute_remix_histo_for_fit(&fit_left, 0);
	for (int i = 0; i < fit_left.naxes[2]; i++) {
		remix_histlayers_backup_left[i] = gsl_histogram_clone(remix_histlayers_left[i]);
	}
}

static void remix_histo_startup_right() {
	compute_remix_histo_for_fit(&fit_right, 1);
	for (int i = 0; i < fit_right.naxes[2]; i++) {
		remix_histlayers_backup_right[i] = gsl_histogram_clone(remix_histlayers_right[i]);
	}
}

void close_histograms() {
	params_histo_left = (ght_params) { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0, 0, TRUE, TRUE, TRUE };
	params_histo_right = (ght_params) { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0, 0, TRUE, TRUE, TRUE };

	if (remix_histlayers_left[0]) {
		for (int i = 0; i < gfit.naxes[2]; i++) {
			gsl_histogram_free(remix_histlayers_left[i]);
			remix_histlayers_left[i] = NULL;
		}
	}
	if (remix_histlayers_right[0]) {
		for (int i = 0; i < gfit.naxes[2]; i++) {
			gsl_histogram_free(remix_histlayers_right[i]);
			remix_histlayers_right[i] = NULL;
		}
	}
	if (remix_histlayers_backup_left[0]) {
		for (int i = 0; i < gfit.naxes[2]; i++) {
			gsl_histogram_free(remix_histlayers_backup_left[i]);
			remix_histlayers_backup_left[i] = NULL;
		}
	}
	if (remix_histlayers_backup_right[0]) {
		for (int i = 0; i < gfit.naxes[2]; i++) {
			gsl_histogram_free(remix_histlayers_backup_right[i]);
			remix_histlayers_backup_right[i] = NULL;
		}
	}
}
/////////////////////////
/////////////////////////
/////////////////////////

static void remixer_startup() {
//	copy_gfit_to_backup();
	return;
}

void initialise_image() {
	clear_stars_list(TRUE);
	com.seq.current = UNRELATED_IMAGE;
	com.uniq = calloc(1, sizeof(single));
	com.uniq->comment = strdup(_("Star recomposition"));
	com.uniq->filename = strdup(_("Unsaved star recomposition result"));
	com.uniq->fileexist = FALSE;
	com.uniq->nb_layers = gfit.naxes[2];
	com.uniq->fit = &gfit;

	initialize_display_mode();
	update_zoom_label();
	sliders_mode_set_state(gui.sliders);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();
	set_display_mode();
	update_prepro_interface(TRUE);
	adjust_sellabel();
	display_filename();	// display filename in gray window
	set_precision_switch(); // set precision on screen
	/* update menus */
	update_MenuItem();
	close_tab();
}

gboolean check_images_match(fits *fit1, fits *fit2) {
	if (fit1->rx != fit2->rx) return FALSE;
	if (fit1->ry != fit2->ry) return FALSE;
	if (fit1->type != fit2->type) return FALSE;
	return TRUE;
}

// imoper_scaled(a, b, oper, factor): a = a oper (factor * b)
int imoper_scaled(fits *a, fits *b, image_operator oper, float factor) {
	size_t n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	float *result;

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		return 1;
	}

	if (a->type == DATA_FLOAT) {
		result = a->fdata;
	}
	else if (a->type == DATA_USHORT) {
		result = malloc(n * sizeof(float));
		if (!result) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}
	else return 1;

	for (size_t i = 0; i < n; ++i) {
		float aval = a->type == DATA_USHORT ? ushort_to_float_bitpix(a, a->data[i]) : a->fdata[i];
		float bval = b->type == DATA_USHORT ? ushort_to_float_bitpix(b, b->data[i]) : b->fdata[i];
		switch (oper) {
			case OPER_ADD:
				result[i] = aval + (bval * factor);
				break;
			case OPER_SUB:
				result[i] = aval - (bval * factor);
				break;
			case OPER_MUL:
				result[i] = aval * bval * factor;
				break;
			case OPER_DIV:
				if (bval == 0.0f)
					result[i] = 0.0f;
				else result[i] = aval / bval * factor;
		}
	}
	if (a->type == DATA_USHORT) {
		for (size_t i = 0; i < n ; i++)
			a->data[i] = float_to_ushort_range(result[i]);
		free(result);
	} else invalidate_stats_from_fit(a);
	return 0;
}

int remixer() {
	// Processing chain
	// (applies to each side)
	// * Apply GHT stretch according to chosen parameters
	// * Apply linear stretch according to chosen blackpoint
	// (combining)
	// Add together (mastermixer * right image) + ((1-mastermixer) * left image)
	// Renormalize to [0,1]

	// Are we allowed to proceed?
	if(!permit_calculation)
		return 1;

	params_left = (ght_params) { leftB, leftD, leftLP, leftSP, leftHP, leftBP, type_left, colour_left, TRUE, TRUE, TRUE };
	params_right = (ght_params) { rightB, rightD, rightLP, rightSP, rightHP, rightBP, type_right, colour_right, TRUE, TRUE, TRUE };
	ght_compute_params cp_left = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	ght_compute_params cp_right = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	// Process left image
	if (left_loaded && (left_changed || leftBP_changed)) {
		apply_linked_ght_to_fits(&fit_left, &fit_left_calc, params_left, cp_left, TRUE);
	}
		// Now do the linear BP shift, if needed. The only parameter that matters is BP so
		// we just need to change the stretch type, no need to recompute params.
	if (left_loaded && (leftBP_changed || (left_changed && leftBP != 0.0))) {
		params_left.stretchtype = STRETCH_LINEAR;
		apply_linked_ght_to_fits(&fit_left_calc, &fit_left_calc, params_left, cp_left, TRUE);
		leftBP_changed = FALSE;
	}
	left_changed = FALSE;

	// Process right image
	if (right_loaded && (right_changed || rightBP_changed)) {
		apply_linked_ght_to_fits(&fit_right, &fit_right_calc, params_right, cp_right, TRUE);
	}
		// As above, BP shift if required.
	if (right_loaded && (rightBP_changed || (right_changed && rightBP != 0.0))) {
		params_right.stretchtype = STRETCH_LINEAR;
		apply_linked_ght_to_fits(&fit_right_calc, &fit_right_calc, params_right, cp_right, TRUE);
		rightBP_changed = FALSE;
	}
	right_changed = FALSE;

	// Combine images together
//	clearfits(&gfit);
//	copyfits(&fit_left_calc, &gfit, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
	const size_t ndata = gfit.naxes[0] * gfit.naxes[1] * gfit.naxes[2];
	if (gfit.data)
		memset(gfit.data, 0, ndata * sizeof(WORD));
	if (gfit.fdata)
		memset(gfit.fdata, 0, ndata * sizeof(float));

	// Scale left and right images according to master mixer slider
	// imoper_scaled ( a, b, oper, factor ) = a oper factor * b
	if (left_loaded && (mastermixer != 1.0)) {
		imoper_scaled(&gfit, &fit_left_calc, OPER_ADD, (1.0 - mastermixer));
	}
	if (right_loaded && (mastermixer != 0.0)) {
		imoper_scaled(&gfit, &fit_right_calc, OPER_ADD, mastermixer);
	}

	if (finalstretch != 0.0) {
		ght_compute_params cp_final = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		ght_params params_final = { (finalstretch - 5.0), 4.0, 0.0, 0.0, 1.0, 0.0, STRETCH_PAYNE_NORMAL, COL_HUMANLUM, TRUE, TRUE, TRUE };
		apply_linked_ght_to_fits(&gfit, &gfit, params_final, cp_final, TRUE);
	}
	// If 16bit preference is set, check the images are 16bit
	if (com.pref.force_16bit && gfit.type == DATA_FLOAT)
		fit_replace_buffer(&gfit, float_buffer_to_ushort(gfit.fdata, ndata), DATA_USHORT);

	notify_gfit_modified();

	return 0;
}

static int remixer_update_preview() {
	remixer();
	return 0;
}

void reset_values() {
	leftD = leftB = leftLP = leftSP = leftBP = 0.0;
	leftHP = 1.0;
	rightD = rightB = rightLP = rightSP = rightBP = 0.0;
	rightHP = 1.0;
	type_left = type_right = STRETCH_PAYNE_NORMAL;
	colour_left = colour_right = COL_INDEP;
	mastermixer = 0.5;
	finalstretch = 0.0;
}

void reset_filechoosers() {
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_left")));
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_right")));
}

void reset_controls_and_values() {
	GtkSpinButton *spin_remix_D_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_D_left"));
	GtkSpinButton *spin_remix_D_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_D_right"));
	GtkSpinButton *spin_remix_B_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_B_left"));
	GtkSpinButton *spin_remix_B_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_B_right"));
	GtkSpinButton *spin_remix_LP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_left"));
	GtkSpinButton *spin_remix_LP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_right"));
	GtkSpinButton *spin_remix_SP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_left"));
	GtkSpinButton *spin_remix_SP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_right"));
	GtkSpinButton *spin_remix_HP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_left"));
	GtkSpinButton *spin_remix_HP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_right"));
	GtkSpinButton *spin_remix_BP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_BP_left"));
	GtkSpinButton *spin_remix_BP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_BP_right"));
	GtkSpinButton *spin_remix_finalstretch = GTK_SPIN_BUTTON(lookup_widget("spin_finalstretch"));
	GtkSpinButton *spin_remix_mastermixer = GTK_SPIN_BUTTON(lookup_widget("spin_mastermixer"));

	reset_values();

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_remix_D_left, leftD);
	gtk_spin_button_set_value(spin_remix_D_right, rightD);
	gtk_spin_button_set_value(spin_remix_B_left, leftB);
	gtk_spin_button_set_value(spin_remix_B_right, rightB);
	gtk_spin_button_set_value(spin_remix_LP_left, leftLP);
	gtk_spin_button_set_value(spin_remix_LP_right, rightLP);
	gtk_spin_button_set_value(spin_remix_SP_left, leftSP);
	gtk_spin_button_set_value(spin_remix_SP_right, rightSP);
	gtk_spin_button_set_value(spin_remix_HP_left, leftHP);
	gtk_spin_button_set_value(spin_remix_HP_right, rightHP);
	gtk_spin_button_set_value(spin_remix_BP_left, leftBP);
	gtk_spin_button_set_value(spin_remix_BP_right, rightBP);
	gtk_spin_button_set_value(spin_remix_finalstretch, finalstretch);
	gtk_spin_button_set_value(spin_remix_mastermixer, mastermixer);
	// Set left_type, right_type, left_colour, right_colour;
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_type_left")), 0);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_type_right")), 0);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_colour_left")), COL_INDEP);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_colour_right")), COL_INDEP);
	set_notify_block(FALSE);
}

static void remixer_close() {
	close_histograms();
	invalidate_stats_from_fit(&gfit);
	clearfits(&fit_left);
	clearfits(&fit_right);
	clearfits(&fit_left_calc);
	clearfits(&fit_right_calc);
	reset_controls_and_values();
	reset_filechoosers();
	permit_calculation = FALSE;
	gtk_button_set_label(GTK_BUTTON(lookup_widget("remix_advanced")), _("Advanced"));
	gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("remix_advanced")), _("Show advanced stretch options."));
	advanced_interface = FALSE;
	siril_close_dialog("dialog_star_remix");
}

void apply_remix_cancel() {
	set_cursor_waiting(TRUE);
	remixer_close();
	if (left_loaded || right_loaded) {
		close_single_image();
	}
	set_cursor_waiting(FALSE);
}

/*** callbacks **/

void on_dialog_star_remix_show(GtkWidget *widget, gpointer user_data) {
	remixer_startup();
	reset_controls_and_values();
	remixer_show_preview = TRUE;

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

int toggle_remixer_window_visibility(int _invocation, const fits* _fit_left, const fits* _fit_right) {
	invocation = _invocation;
	if (gtk_widget_get_visible(lookup_widget("dialog_star_remix"))) {
		set_cursor_waiting(TRUE);
		reset_controls_and_values();
		remixer_close();
		set_cursor_waiting(FALSE);
		siril_close_dialog("dialog_star_remix");
	} else {
		left_loaded = FALSE;
		right_loaded = FALSE;
		if (invocation == CALL_FROM_STARNET) {
			fit_left = *_fit_left;
			remix_histo_startup_left();
			copyfits(&fit_left, &fit_left_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
			close_single_image();
			copyfits(&fit_left, &gfit, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
			left_loaded = TRUE; // Mark LHS image as loaded
			left_changed = TRUE; // Force update on initial draw
			permit_calculation = TRUE;
			fit_right = *_fit_right;
			remix_histo_startup_right();
			copyfits(&fit_right, &fit_right_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
			right_loaded = TRUE; // Mark RHS image as loaded
			right_changed = TRUE; // Force update on initial draw
			initialise_image();

			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("remix_filechooser_left")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("remix_filechooser_right")), FALSE);
			update_image *param = malloc(sizeof(update_image));
			param->update_preview_fn = remixer_update_preview;
			param->show_preview = TRUE;
			notify_update((gpointer) param);

		} else {
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("remix_filechooser_left")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("remix_filechooser_right")), TRUE);
			gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_left")), com.wd);
			gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_right")), com.wd);
		}
		// Set eyedropper icons to light or dark according to theme
		gchar *image;
		GtkWidget *v, *w;
		if (com.pref.gui.combo_theme == 0) {
			image = g_build_filename(siril_get_system_data_dir(), "pixmaps", "eyedropper_dark.svg", NULL);
			v = gtk_image_new_from_file(image);
			w = gtk_image_new_from_file(image);
		} else {
			image = g_build_filename(siril_get_system_data_dir(), "pixmaps", "eyedropper.svg", NULL);
			v = gtk_image_new_from_file(image);
			w = gtk_image_new_from_file(image);
		}
		gtk_button_set_image(GTK_BUTTON(lookup_widget("eyedropper_SP_left")), v);
		gtk_button_set_image(GTK_BUTTON(lookup_widget("eyedropper_SP_right")), w);
		gtk_widget_show(v);
		gtk_widget_show(w);
		g_free(image);

		// Hide the advanced widgets, these can be show using the Advanced button for full control
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtStretchTypecontrols3")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols2")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols2")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtHPcontrols2")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols2")), FALSE);

		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtStretchTypecontrols1")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols1")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols1")), FALSE);
//		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_SP_right")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols1")), FALSE);

		siril_open_dialog("dialog_star_remix");
	}
	return 0;
}

void on_remix_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_remix_cancel();
}

void on_remix_reset_clicked(GtkButton *button, gpointer user_data) {
	reset_controls_and_values();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_apply_clicked(GtkButton *button, gpointer user_data) {
	close_histograms();
	remixer_close();
	set_cursor_waiting(FALSE);
}

/*** adjusters **/
void on_spin_remix_D_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftD = expm1(gtk_spin_button_get_value(button));
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_D_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightD = expm1(gtk_spin_button_get_value(button));
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_B_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftB = gtk_spin_button_get_value(button);
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}
void on_spin_remix_B_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightB = gtk_spin_button_get_value(button);
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_LP_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftLP = gtk_spin_button_get_value(button);
	if (leftLP > leftSP) {
		leftLP = leftSP;
		gtk_spin_button_set_value(button, leftLP);
	}
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_LP_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightLP = gtk_spin_button_get_value(button);
	if (rightLP > rightSP) {
		rightLP = rightSP;
		gtk_spin_button_set_value(button, rightLP);
	}
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_SP_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	GtkSpinButton *spin_remix_LP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_left"));
	GtkSpinButton *spin_remix_HP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_left"));
	left_changed = TRUE;
	leftSP = gtk_spin_button_get_value(button);
	if (leftSP < leftLP) {
		gtk_spin_button_set_value(spin_remix_LP_left, leftSP);
		leftLP = leftSP;
	}
	if (leftSP > leftHP) {
		leftHP = leftSP;
		gtk_spin_button_set_value(spin_remix_HP_left, leftSP);
	}
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}
void on_spin_remix_SP_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	GtkSpinButton *spin_remix_LP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_right"));
	GtkSpinButton *spin_remix_HP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_right"));
	right_changed = TRUE;
	rightSP = gtk_spin_button_get_value(button);
	if (rightSP < rightLP) {
		rightLP = rightSP;
		gtk_spin_button_set_value(spin_remix_LP_right, rightSP);
	}
	if (rightSP > rightHP) {
		rightHP = rightSP;
		gtk_spin_button_set_value(spin_remix_HP_right, rightSP);
	}
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_HP_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftHP = gtk_spin_button_get_value(button);
	if (leftHP < leftSP) {
		leftHP = leftSP;
		gtk_spin_button_set_value(button, leftHP);
	}
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}
void on_spin_remix_HP_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightHP = gtk_spin_button_get_value(button);
	if (rightHP < rightSP) {
		rightHP = rightSP;
		gtk_spin_button_set_value(button, rightHP);
	}
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_BP_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	leftBP_changed = TRUE;
	leftBP = gtk_spin_button_get_value(button);
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_BP_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	rightBP_changed = TRUE;
	rightBP = gtk_spin_button_get_value(button);
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_mastermixer_value_changed(GtkSpinButton *button, gpointer user_data) {
	mastermixer = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_finalstretch_value_changed(GtkSpinButton *button, gpointer user_data) {
	finalstretch = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_colour_left_changed(GtkComboBox *combo, gpointer user_data) {
	left_changed = TRUE;
	colour_left = gtk_combo_box_get_active(combo);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_colour_right_changed(GtkComboBox *combo, gpointer user_data) {
	right_changed = TRUE;
	colour_right = gtk_combo_box_get_active(combo);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_type_left_changed(GtkComboBox *combo, gpointer user_data) {
	left_changed = TRUE;
	type_left = gtk_combo_box_get_active(combo);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_type_right_changed(GtkComboBox *combo, gpointer user_data) {
	right_changed = TRUE;
	type_right = gtk_combo_box_get_active(combo);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_filechooser_left_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	filename_left = g_strdup(gtk_file_chooser_get_filename(filechooser));
	if (readfits(filename_left, &fit_left, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		right_loaded = FALSE;
		return;
	}
	remix_histo_startup_left();
	if(fit_left.type == DATA_FLOAT && com.pref.force_16bit) {
		const size_t ndata = fit_left.naxes[0] * fit_left.naxes[1] * fit_left.naxes[2];
		fit_replace_buffer(&fit_left, float_buffer_to_ushort(fit_left.fdata, ndata), DATA_USHORT);
	}
	if (right_loaded) {
		if(!check_images_match(&fit_left, &fit_right)) {
			siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: images do not match"),
				_("Image width, height and bit depth must match"));
			clearfits(&fit_left);
			gtk_file_chooser_unselect_all(filechooser);
			left_loaded = FALSE;
			return;
		}
		else {
			left_loaded = TRUE;
		}
	} else {
		close_single_image();
		clearfits(&gfit);
		copyfits(&fit_left, &gfit, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
		initialise_image();
		left_loaded = TRUE;
		clearfits(&fit_right_calc);
		copyfits(&fit_left, &fit_right_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
	}
	if (left_loaded || right_loaded)
		permit_calculation = TRUE;
	else
		permit_calculation = FALSE;
	clearfits(&fit_left_calc);
	copyfits(&fit_left, &fit_left_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
	left_changed = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_filechooser_right_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	filename_right = g_strdup(gtk_file_chooser_get_filename(filechooser));
	if (readfits(filename_right, &fit_right, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
				_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		right_loaded = FALSE;
		return;
	}
	remix_histo_startup_right();
	if(fit_right.type == DATA_FLOAT && com.pref.force_16bit) {
		const size_t ndata = fit_right.naxes[0] * fit_right.naxes[1] * fit_right.naxes[2];
		fit_replace_buffer(&fit_right, float_buffer_to_ushort(fit_right.fdata, ndata), DATA_USHORT);
	}
	if (left_loaded) {
		if(!check_images_match(&fit_left, &fit_right)) {
			siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: images do not match"),
				_("Image width, height and bit depth must match"));
			clearfits(&fit_right);
			gtk_file_chooser_unselect_all(filechooser);
			right_loaded = FALSE;
			return;
		}
		else {
			right_loaded = TRUE;
		}
	} else {
		close_single_image();
		clearfits(&gfit);
		copyfits(&fit_right, &gfit, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
		initialise_image();
		right_loaded = TRUE;
		clearfits(&fit_left_calc);
		copyfits(&fit_right, &fit_left_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
	}
	if (left_loaded || right_loaded)
		permit_calculation = TRUE;
	else
		permit_calculation = FALSE;
	clearfits(&fit_right_calc);
	copyfits(&fit_right, &fit_right_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
	right_changed = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_eyedropper_SP_left_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_remix_LP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_left"));
	GtkSpinButton *spin_remix_HP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_left"));
	int chan, channels = fit_left.naxes[2];
	imstats* stats[3];
	double ref = 0.0;
	double norm = 1.0;
	if (get_preview_gfit_backup()->type == DATA_USHORT)
		norm = USHRT_MAX_DOUBLE;

	if (!com.selection.w || !com.selection.h) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the point or area to set SP"));
		return;
	}
	for (chan = 0; chan < fit_left.naxes[2]; chan++) {
		stats[chan] = statistics(NULL, -1, &fit_left, chan, &com.selection, STATS_BASIC, MULTI_THREADED);
		if (!stats[chan]) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return;
		}
		ref += stats[chan]->mean;
		free_stats(stats[chan]);
	}
	ref /= channels;
	ref /= norm;
	leftSP = ref;
	if (leftSP < leftLP) {
		gtk_spin_button_set_value(spin_remix_LP_left, leftSP);
		rightLP = rightSP;
	}
	if (leftSP > leftHP) {
		gtk_spin_button_set_value(spin_remix_HP_left, leftSP);
		rightHP = rightSP;
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_left")),leftSP);
	left_changed = TRUE;
//	if (remix_histos_visible())
		update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE; // no need of preview button. This is always in preview
	notify_update((gpointer) param);
}

void on_eyedropper_SP_right_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_remix_LP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_right"));
	GtkSpinButton *spin_remix_HP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_right"));
	int chan, channels = fit_right.naxes[2];
	imstats* stats[3];
	double ref = 0.0;
	double norm = 1.0;
	if (get_preview_gfit_backup()->type == DATA_USHORT)
		norm = USHRT_MAX_DOUBLE;
	if (!com.selection.w || !com.selection.h) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the point or area to set SP"));
		return;
	}
	for (chan = 0; chan < fit_right.naxes[2]; chan++) {
		stats[chan] = statistics(NULL, -1, &fit_right, chan, &com.selection, STATS_BASIC, MULTI_THREADED);
		if (!stats[chan]) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return;
		}
		ref += stats[chan]->mean;
		free_stats(stats[chan]);
	}
	ref /= channels;
	ref /= norm;
	rightSP = ref;
	if (rightSP < rightLP) {
		gtk_spin_button_set_value(spin_remix_LP_right, rightSP);
		rightLP = rightSP;
	}
	if (rightSP > rightHP) {
		gtk_spin_button_set_value(spin_remix_HP_right, rightSP);
		rightHP = rightSP;
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_right")),leftSP);
	right_changed = TRUE;
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE; // no need of preview button. This is always in preview
	notify_update((gpointer) param);
}

void on_remix_advanced_clicked(GtkButton *button, gpointer user_data) {
	// Show all the widgets
	if (!advanced_interface) {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtStretchTypecontrols3")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols2")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols2")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtHPcontrols2")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols2")), TRUE);

		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtStretchTypecontrols1")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols1")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols1")), TRUE);
//		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_SP_right")), TRUE);

		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols1")), TRUE);
		gtk_button_set_label(GTK_BUTTON(lookup_widget("remix_advanced")), _("Simple"));
		gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("remix_advanced")), _("Show only simple stretch options. Note that any parameters that have been set while in advanced mode will remain set."));
		advanced_interface = TRUE;
	}
	else {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtStretchTypecontrols3")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols2")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols2")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtHPcontrols2")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols2")), FALSE);

		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtStretchTypecontrols1")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols1")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols1")), FALSE);
//		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_SP_right")), FALSE);

		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols1")), FALSE);
		gtk_button_set_label(GTK_BUTTON(lookup_widget("remix_advanced")), _("Advanced"));
		gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("remix_advanced")), _("Show advanced stretch options."));
		advanced_interface = FALSE;
	}
}

void on_dialog_star_remix_close() {
	apply_remix_cancel();
}

void on_remixer_help_button_clicked(GtkButton *button, gpointer user_data) {
	siril_open_dialog("remixer_help");
}

void on_remixer_help_okay_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("remixer_help");
}

void on_remix_histos_button_clicked(GtkButton *button, gpointer user_data) {
	siril_open_dialog("remixer_histos");
}
