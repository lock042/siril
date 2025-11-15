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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "algos/statistics.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
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
#include "core/arithm.h"
#include "core/siril_log.h"
#include "filters/ght.h"

#include <gsl/gsl_histogram.h>

// Invocation: 1 if called directly from starnet GUI,
// 2 if called from Siril menu
static int invocation = 0;

static clip_mode_t clip_mode = RGBBLEND;
static float leftD = 0.0f, rightD = 0.0f;
static float leftB = 0.0f, rightB = 0.0f;
static float leftLP = 0.0f, rightLP = 0.0f;
static float leftSP = 0.0f, rightSP = 0.0f;
static float leftHP = 1.0f, rightHP = 1.0f;
static float leftBP = 0.0f, rightBP = 0.0f;
static int type_left = STRETCH_PAYNE_NORMAL, type_right = STRETCH_ASINH;
static int colour_left = COL_HUMANLUM, colour_right = COL_HUMANLUM;
static fits fit_left;
static fits fit_right;
static fits fit_left_calc;
static fits fit_right_calc;

static ght_params params_left, params_right, params_histo_left, params_histo_right;
static ght_compute_params cp_histo_left = { 0.0f };
static ght_compute_params cp_histo_right = { 0.0f };

static gboolean left_changed = FALSE, right_changed = FALSE;
static gboolean leftBP_changed = FALSE, rightBP_changed = FALSE;
static gchar* filename_left;
static gchar* filename_right;

static gboolean advanced_interface = FALSE;
static gboolean permit_calculation = FALSE;
static gboolean remixer_show_preview;
static gboolean left_loaded = FALSE;
static gboolean right_loaded = FALSE;
static gboolean remix_log_scale = FALSE;

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
		float y = GHT(x, params->B, params->D, params->LP, params->SP, params->HP, 0.0, params->stretchtype, cp);
		cairo_line_to(cr, k, height * (1 - y));
	}
	cairo_stroke(cr);
}

void display_remix_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height, double zoomH, double zoomV, gboolean isOrig) {
	if (width <= 0) return;
	int current_bin;
	size_t norm = gsl_histogram_bins(histo) - 1;

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
	if (!graph_height)
		return;
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
	params_histo_left = (ght_params) { leftB, leftD, leftLP, leftSP, leftHP, leftBP, type_left, colour_left, TRUE, TRUE, TRUE, clip_mode };

	GHTsetup(&cp_histo_left, params_histo_left.B, params_histo_left.D, params_histo_left.LP, params_histo_left.SP, params_histo_left.HP, params_histo_left.stretchtype);
	for (size_t i = 0; i < fit_left.naxes[2]; i++) {
		gsl_histogram_memcpy(remix_histlayers_left[i], remix_histlayers_backup_left[i]);
		apply_remix_histo(remix_histlayers_left[i], norm, &params_histo_left, &cp_histo_left);
	}
	if (leftBP_changed || leftBP != 0.f) {
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
	params_histo_right = (ght_params) { rightB, rightD, rightLP, rightSP, rightHP, rightBP, type_right, colour_right, TRUE, TRUE, TRUE, clip_mode };

	int nlayers = fit_right.naxes[2];
	GHTsetup(&cp_histo_right, params_histo_right.B, params_histo_right.D, params_histo_right.LP, params_histo_right.SP, params_histo_right.HP, params_histo_right.stretchtype);
	for (size_t i = 0; i < nlayers; i++) {
		gsl_histogram_memcpy(remix_histlayers_right[i], remix_histlayers_backup_right[i]);
		apply_remix_histo(remix_histlayers_right[i], norm, &params_histo_right, &cp_histo_right);
	}
	if (rightBP_changed || rightBP != 0.f) {
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

void close_histograms(gboolean clear_left, gboolean clear_right) {
	params_histo_left = (ght_params) { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, STRETCH_PAYNE_NORMAL, COL_INDEP, TRUE, TRUE, TRUE };
	params_histo_right = (ght_params) { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, STRETCH_PAYNE_NORMAL, COL_INDEP, TRUE, TRUE, TRUE };

	if (clear_left) {
		if (remix_histlayers_left[0]) {
			for (int i = 0; i < gfit.naxes[2]; i++) {
				gsl_histogram_free(remix_histlayers_left[i]);
				remix_histlayers_left[i] = NULL;
			}
		}
		if (remix_histlayers_backup_left[0]) {
			for (int i = 0; i < gfit.naxes[2]; i++) {
				gsl_histogram_free(remix_histlayers_backup_left[i]);
				remix_histlayers_backup_left[i] = NULL;
			}
		}
	}
	if (clear_right) {
		if (remix_histlayers_right[0]) {
			for (int i = 0; i < gfit.naxes[2]; i++) {
				gsl_histogram_free(remix_histlayers_right[i]);
				remix_histlayers_right[i] = NULL;
			}
		}
		if (remix_histlayers_backup_right[0]) {
			for (int i = 0; i < gfit.naxes[2]; i++) {
				gsl_histogram_free(remix_histlayers_backup_right[i]);
				remix_histlayers_backup_right[i] = NULL;
			}
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
	if (!create_uniq_from_gfit(strdup(_("Unsaved star recomposition result")), FALSE))
		com.uniq->comment = strdup(_("Star recomposition"));

	initialize_display_mode();
	update_zoom_label();
	sliders_mode_set_state(gui.sliders);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();
	set_display_mode();
	update_prepro_interface(TRUE);
	adjust_sellabel();
	display_filename();	// display filename in gray window
	gui_function(set_precision_switch, NULL); // set precision on screen
	/* update menus */
	gui_function(update_MenuItem, NULL);
	gui_function(close_tab, NULL);
	gui_function(init_right_tab, NULL);
}

static gboolean check_images_match(const fits *fit1, const fits *fit2) {
	if (fit1->rx != fit2->rx) return FALSE;
	if (fit1->ry != fit2->ry) return FALSE;
	if (fit1->type != fit2->type) return FALSE;
	if (fit1->naxes[2] != fit2->naxes[2]) return FALSE;
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
	// Combine in CIE LCh colorspace
	// Renormalize to [0,1]

	// Are we allowed to proceed?
	if(!permit_calculation)
		return 1;

	params_left = (ght_params) { leftB, leftD, leftLP, leftSP, leftHP, leftBP, type_left, colour_left, TRUE, TRUE, TRUE, clip_mode };
	params_right = (ght_params) { rightB, rightD, rightLP, rightSP, rightHP, rightBP, type_right, colour_right, TRUE, TRUE, TRUE, clip_mode };

	// Process left image
	if (left_loaded && (left_changed || leftBP_changed)) {
		apply_linked_ght_to_fits(&fit_left, &fit_left_calc, &params_left, TRUE);
	}
		// Now do the linear BP shift, if needed. The only parameter that matters is BP so
		// we just need to change the stretch type, no need to recompute params.
	if (left_loaded && (leftBP_changed || (left_changed && leftBP != 0.0f))) {
		params_left.stretchtype = STRETCH_LINEAR;
		apply_linked_ght_to_fits(&fit_left_calc, &fit_left_calc, &params_left, TRUE);
		leftBP_changed = FALSE;
	}
	left_changed = FALSE;

	// Process right image
	if (right_loaded && (right_changed || rightBP_changed)) {
		apply_linked_ght_to_fits(&fit_right, &fit_right_calc, &params_right, TRUE);
	}
		// As above, BP shift if required.
	if (right_loaded && (rightBP_changed || (right_changed && rightBP != 0.0f))) {
		params_right.stretchtype = STRETCH_LINEAR;
		apply_linked_ght_to_fits(&fit_right_calc, &fit_right_calc, &params_right, TRUE);
		rightBP_changed = FALSE;
	}
	right_changed = FALSE;

	// Combine images together

	const size_t ndata = gfit.naxes[0] * gfit.naxes[1] * gfit.naxes[2];
	if (gfit.data)
		memset(gfit.data, 0, ndata * sizeof(WORD));
	if (gfit.fdata)
		memset(gfit.fdata, 0, ndata * sizeof(float));

	size_t npixels = gfit.naxes[0] * gfit.naxes[1];
	const float norm = USHRT_MAX_SINGLE;
	const float invnorm = 1.f / norm;
	if (gfit.naxes[2] == 1) {
		switch (gfit.type) {
			case DATA_FLOAT:
				if (gfit.fdata) {
					for (size_t i = 0 ; i < npixels ; i++) {
						gfit.fdata[i] = fit_left_calc.fdata[i] + fit_right_calc.fdata[i] - fit_left_calc.fdata[i] * fit_right_calc.fdata[i];
					}
				}
				break;
			case DATA_USHORT:
				if (gfit.data) {
					for (size_t i = 0 ; i < npixels ; i++) {
						gfit.data[i] = roundf_to_WORD(norm * (1.f - (1.f - (fit_left_calc.data[i]*invnorm))*(1.f - (fit_right_calc.data[i]*invnorm))));
					}
				}
				break;
			default:
				break;
		}
	} else {
		switch (gfit.type) {
			case DATA_FLOAT:
				if (left_loaded && !right_loaded) {
					memcpy(gfit.fdata, fit_left_calc.fdata, npixels * 3 * sizeof(float));
				} else if (right_loaded && !left_loaded) {
					memcpy(gfit.fdata, fit_right_calc.fdata, npixels * 3 * sizeof(float));
				} else {
					screen(&fit_left_calc, &fit_right_calc, &gfit, !com.pref.force_16bit, com.max_thread);
				}
				break;
			case DATA_USHORT:
				if (left_loaded && !right_loaded) {
					memcpy(gfit.data, fit_left_calc.data, npixels * 3 * sizeof(WORD));
				} else if (right_loaded && !left_loaded) {
					memcpy(gfit.data, fit_right_calc.data, npixels * 3 * sizeof(WORD));
				} else {
					screen(&fit_left_calc, &fit_right_calc, &gfit, !com.pref.force_16bit, com.max_thread);
				}
				break;
			default:
				break;
		}
	}

	notify_gfit_modified();

	return 0;
}

static int remixer_update_preview() {
	remixer();
	return 0;
}

void reset_values_left() {
	leftD = leftB = leftLP = leftSP = leftBP = 0.0f;
	leftHP = 1.0f;
	type_left = STRETCH_PAYNE_NORMAL;
	colour_left = COL_HUMANLUM;
}

void reset_values_right() {
	rightD = rightB = rightLP = rightSP = rightBP = 0.0f;
	rightHP = 1.0f;
	type_right = STRETCH_ASINH;
	colour_right = COL_HUMANLUM;
}

void reset_values() {
	reset_values_left();
	reset_values_right();
}

void reset_filechoosers() {
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_left")));
	gtk_file_chooser_unselect_all(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_right")));
}

void reset_left() {
	GtkSpinButton *spin_remix_D_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_D_left"));
	GtkSpinButton *spin_remix_B_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_B_left"));
	GtkSpinButton *spin_remix_LP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_left"));
	GtkSpinButton *spin_remix_SP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_left"));
	GtkSpinButton *spin_remix_HP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_left"));
	GtkSpinButton *spin_remix_BP_left = GTK_SPIN_BUTTON(lookup_widget("spin_remix_BP_left"));
	reset_values_left();
	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_remix_D_left, leftD);
	gtk_spin_button_set_value(spin_remix_B_left, leftB);
	gtk_spin_button_set_value(spin_remix_LP_left, leftLP);
	gtk_spin_button_set_value(spin_remix_SP_left, leftSP);
	gtk_spin_button_set_value(spin_remix_HP_left, leftHP);
	gtk_spin_button_set_value(spin_remix_BP_left, leftBP);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_type_left")), STRETCH_PAYNE_NORMAL);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_colour_left")), COL_HUMANLUM);
	set_notify_block(FALSE);
}

void reset_right() {
	GtkSpinButton *spin_remix_D_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_D_right"));
	GtkSpinButton *spin_remix_B_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_B_right"));
	GtkSpinButton *spin_remix_LP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_LP_right"));
	GtkSpinButton *spin_remix_SP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_right"));
	GtkSpinButton *spin_remix_HP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_HP_right"));
	GtkSpinButton *spin_remix_BP_right = GTK_SPIN_BUTTON(lookup_widget("spin_remix_BP_right"));
	reset_values_right();
	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_remix_D_right, rightD);
	gtk_spin_button_set_value(spin_remix_B_right, rightB);
	gtk_spin_button_set_value(spin_remix_LP_right, rightLP);
	gtk_spin_button_set_value(spin_remix_SP_right, rightSP);
	gtk_spin_button_set_value(spin_remix_HP_right, rightHP);
	gtk_spin_button_set_value(spin_remix_BP_right, rightBP);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_type_right")), STRETCH_ASINH);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("remix_colour_right")), COL_HUMANLUM);
	set_notify_block(FALSE);
}

void reset_controls_and_values() {
	reset_left();
	reset_right();
}

static void remixer_close() {
	close_histograms(TRUE, TRUE);
	invalidate_stats_from_fit(&gfit);
	clearfits(&fit_left);
	clearfits(&fit_right);
	clearfits(&fit_left_calc);
	clearfits(&fit_right_calc);
	reset_controls_and_values();
	reset_filechoosers();
	permit_calculation = FALSE;
	gtk_button_set_label(GTK_BUTTON(lookup_widget("remix_advanced")), _("Advanced"));
	gtk_widget_set_tooltip_text(lookup_widget("remix_advanced"), _("Show advanced stretch options."));
	advanced_interface = FALSE;

	siril_close_dialog("dialog_star_remix");
}

void apply_remix_cancel() {
	remixer_close();
	if (left_loaded || right_loaded || (!single_image_is_loaded() && !sequence_is_loaded())) {
		close_single_image();
	}
	set_cursor_waiting(FALSE);
}

/*** callbacks **/

void on_dialog_star_remix_show(GtkWidget *widget, gpointer user_data) {
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_left")), "filefilter_fits", "FITS files");
	siril_set_file_filter(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_right")), "filefilter_fits", "FITS files");
	remixer_startup();
	reset_controls_and_values();
	remixer_show_preview = TRUE;

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

int toggle_remixer_window_visibility(int _invocation, fits* _fit_left, fits* _fit_right) {
	invocation = _invocation;
	if (gtk_widget_get_visible(lookup_widget("dialog_star_remix"))) {
		set_cursor_waiting(TRUE);
		reset_controls_and_values();

		remixer_close();
		set_cursor_waiting(FALSE);
		siril_close_dialog("dialog_star_remix");
	} else {
		if (gui.rendering_mode == LINEAR_DISPLAY)
			setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535
		left_loaded = FALSE;
		right_loaded = FALSE;
		remix_log_scale = (com.pref.gui.display_histogram_mode == LOG_DISPLAY ? TRUE : FALSE);
		GtkToggleButton *toggle_log = GTK_TOGGLE_BUTTON(lookup_widget("toggle_remixer_log_histograms"));
		gtk_toggle_button_set_active(toggle_log, remix_log_scale);
		if (invocation == CALL_FROM_STARNET) {
			copyfits(_fit_left, &fit_left, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
			copy_fits_metadata(_fit_left, &fit_left);
			clearfits(_fit_left);
			free(_fit_left);
			close_histograms(TRUE, TRUE);
			remix_histo_startup_left();
			copyfits(&fit_left, &fit_left_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
			close_single_image();
			copyfits(&fit_left, &gfit, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
			left_loaded = TRUE; // Mark LHS image as loaded
			left_changed = TRUE; // Force update on initial draw
			permit_calculation = TRUE;
			copyfits(_fit_right, &fit_right, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
			copy_fits_metadata(_fit_right, &fit_right);
			remix_histo_startup_right();
			clearfits(_fit_right);
			free(_fit_right);
			copyfits(&fit_right, &fit_right_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
			right_loaded = TRUE; // Mark RHS image as loaded
			right_changed = TRUE; // Force update on initial draw
			merge_fits_headers_to_result(&gfit, FALSE, &fit_left, &fit_right, NULL);
			// Avoid doubling STACKCNT and LIVETIME as we are merging starless and star parts of a single image
			gfit.keywords.stackcnt = fit_left.keywords.stackcnt;
			gfit.keywords.livetime = fit_left.keywords.livetime;
			// Update the header
			update_fits_header(&gfit);
			initialise_image();
			GtkWidget *clip = lookup_widget("remixer_clip_mode_settings");
			gtk_widget_set_visible(clip, (fit_left.naxes[2] == 3));
			gtk_widget_set_visible(lookup_widget("remix_filechooser_left"), FALSE);
			gtk_widget_set_visible(lookup_widget("remix_filechooser_right"), FALSE);
			update_image *param = malloc(sizeof(update_image));
			param->update_preview_fn = remixer_update_preview;
			param->show_preview = TRUE;
			notify_update((gpointer) param);

		} else {
			gtk_widget_set_visible(lookup_widget("remix_filechooser_left"), TRUE);
			gtk_widget_set_visible(lookup_widget("remix_filechooser_right"), TRUE);
			gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_left")), com.wd);
			gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(lookup_widget("remix_filechooser_right")), com.wd);
		}
		// Set eyedropper icons to light or dark according to theme
		GtkWidget *v = NULL, *w = NULL;
		if (com.pref.gui.combo_theme == 0) {
			v = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper_dark.svg");
			w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper_dark.svg");
		} else {
			v = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper.svg");
			w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper.svg");
		}
		gtk_button_set_image(GTK_BUTTON(lookup_widget("eyedropper_SP_left")), v);
		gtk_button_set_image(GTK_BUTTON(lookup_widget("eyedropper_SP_right")), w);
		gtk_widget_show(v);
		gtk_widget_show(w);

		// Hide the advanced widgets, these can be show using the Advanced button for full control
		gtk_widget_set_visible(lookup_widget("ghtStretchTypecontrols3"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtColourmodelcontrols2"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtLPcontrols2"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtHPcontrols2"), FALSE);

		gtk_widget_set_visible(lookup_widget("ghtStretchTypecontrols1"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtColourmodelcontrols1"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtLPcontrols1"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtBcontrols1"), FALSE);

		siril_open_dialog("dialog_star_remix");
	}
	return 0;
}

gboolean on_remix_close_clicked(GtkButton *button, gpointer user_data) {
	close_histograms(TRUE, TRUE);
	remixer_close();
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_remix_reset_left_clicked(GtkButton *button, gpointer user_data) {
	reset_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_reset_right_clicked(GtkButton *button, gpointer user_data) {
	reset_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_apply_left_clicked(GtkButton *button, gpointer user_data) {
	switch (fit_left.type) {
		case DATA_FLOAT:
			memcpy(fit_left.fdata, fit_left_calc.fdata, fit_left.rx * fit_left.ry * fit_left.naxes[2] * sizeof(float));
			break;
		case DATA_USHORT:
			memcpy(fit_left.data, fit_left_calc.data, fit_left.rx * fit_left.ry * fit_left.naxes[2] * sizeof(WORD));
			break;
		default:
			break;
	}
	reset_left();
	close_histograms(TRUE, FALSE);
	remix_histo_startup_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);

}

void on_remix_apply_right_clicked(GtkButton *button, gpointer user_data) {
	switch (fit_right.type) {
		case DATA_FLOAT:
			memcpy(fit_right.fdata, fit_right_calc.fdata, fit_right.rx * fit_right.ry * fit_right.naxes[2] * sizeof(float));
			break;
		case DATA_USHORT:
			memcpy(fit_right.data, fit_right_calc.data, fit_right.rx * fit_right.ry * fit_right.naxes[2] * sizeof(WORD));
			break;
		default:
			break;
	}
	reset_right();
	close_histograms(FALSE, TRUE);
	remix_histo_startup_right();
}

/*** adjusters **/
void on_spin_remix_D_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftD = expm1f((float)gtk_spin_button_get_value(button));
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_D_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightD = expm1f((float)gtk_spin_button_get_value(button));
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_B_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftB = (float) gtk_spin_button_get_value(button);
	if (fabsf(leftB) < 1.e-3f) {
		leftB = 0.f;
		gtk_spin_button_set_value(button, 0.f);
	}
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}
void on_spin_remix_B_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightB = (float) gtk_spin_button_get_value(button);
	if (fabsf(rightB) < 1.e-3f) {
		rightB = 0.f;
		gtk_spin_button_set_value(button, 0.f);
	}	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_LP_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftLP = (float) gtk_spin_button_get_value(button);
	if (leftLP > leftSP) {
		leftLP = leftSP;
		gtk_spin_button_set_value(button, (double) leftLP);
	}
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_LP_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightLP = (float) gtk_spin_button_get_value(button);
	if (rightLP > rightSP) {
		rightLP = rightSP;
		gtk_spin_button_set_value(button, (float) rightLP);
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
	leftSP = (float) gtk_spin_button_get_value(button);
	if (leftSP < leftLP) {
		gtk_spin_button_set_value(spin_remix_LP_left, (double) leftSP);
		leftLP = leftSP;
	}
	if (leftSP > leftHP) {
		leftHP = leftSP;
		gtk_spin_button_set_value(spin_remix_HP_left, (double) leftSP);
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
	rightSP = (float) gtk_spin_button_get_value(button);
	if (rightSP < rightLP) {
		rightLP = rightSP;
		gtk_spin_button_set_value(spin_remix_LP_right, (double) rightSP);
	}
	if (rightSP > rightHP) {
		rightHP = rightSP;
		gtk_spin_button_set_value(spin_remix_HP_right, (double) rightSP);
	}
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_HP_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	left_changed = TRUE;
	leftHP = (float) gtk_spin_button_get_value(button);
	if (leftHP < leftSP) {
		leftHP = leftSP;
		gtk_spin_button_set_value(button, (double) leftHP);
	}
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}
void on_spin_remix_HP_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	right_changed = TRUE;
	rightHP = (float) gtk_spin_button_get_value(button);
	if (rightHP < rightSP) {
		rightHP = rightSP;
		gtk_spin_button_set_value(button, (double) rightHP);
	}
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_spin_remix_BP_left_value_changed(GtkSpinButton *button, gpointer user_data) {
	leftBP_changed = TRUE;
	leftBP = (float) gtk_spin_button_get_value(button);
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
	if (leftBP == 0) {
		gtk_widget_set_sensitive(lookup_widget("eyedropper_SP_left"), TRUE);
		gtk_widget_set_tooltip_text(lookup_widget("eyedropper_SP_left"), _("This dropper allows you to set the background stretch symmetry point  based on the average value of the current selection. This helps in picking out faint detail above the dark background."));
	} else {
		gtk_widget_set_sensitive(lookup_widget("eyedropper_SP_left"), FALSE);
		gtk_widget_set_tooltip_text(lookup_widget("eyedropper_SP_left"), _("Eyedropper disabled with unapplied BP shift"));
	}
}

void on_spin_remix_BP_right_value_changed(GtkSpinButton *button, gpointer user_data) {
	rightBP_changed = TRUE;
	rightBP = (float) gtk_spin_button_get_value(button);
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
	if (rightBP == 0) {
		gtk_widget_set_sensitive(lookup_widget("eyedropper_SP_right"), TRUE);
		gtk_widget_set_tooltip_text(lookup_widget("eyedropper_SP_right"), _("This dropper allows you to set the background stretch symmetry point  based on the average value of the current selection. This helps in picking out faint detail above the dark background."));
	} else {
		gtk_widget_set_sensitive(lookup_widget("eyedropper_SP_right"), FALSE);
		gtk_widget_set_tooltip_text(lookup_widget("eyedropper_SP_right"), _("Eyedropper disabled with unapplied BP shift"));
	}
}

void on_remix_colour_left_changed(GtkComboBox *combo, gpointer user_data) {
	left_changed = TRUE;
	colour_left = gtk_combo_box_get_active(combo);
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_colour_right_changed(GtkComboBox *combo, gpointer user_data) {
	right_changed = TRUE;
	colour_right = gtk_combo_box_get_active(combo);
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_type_left_changed(GtkComboBox *combo, gpointer user_data) {
	left_changed = TRUE;
	type_left = gtk_combo_box_get_active(combo);
	if (type_left == STRETCH_ASINH || type_left == STRETCH_INVASINH)
		gtk_widget_set_visible(lookup_widget("ghtBcontrols2"), FALSE);
	else
		gtk_widget_set_visible(lookup_widget("ghtBcontrols2"), TRUE);
	update_remix_histo_left();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_type_right_changed(GtkComboBox *combo, gpointer user_data) {
	right_changed = TRUE;
	type_right = gtk_combo_box_get_active(combo);
	if (type_right == STRETCH_ASINH || type_right == STRETCH_INVASINH)
		gtk_widget_set_visible(lookup_widget("ghtBcontrols1"), FALSE);
	else
		gtk_widget_set_visible(lookup_widget("ghtBcontrols1"), TRUE);
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}

void on_remix_filechooser_left_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	close_histograms(TRUE, FALSE);
	if (left_loaded) {
		clearfits(&fit_left);
		clearfits(&fit_left_calc);
		left_loaded = FALSE;
	}
	filename_left = siril_file_chooser_get_filename(filechooser);
	if (readfits(filename_left, &fit_left, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		right_loaded = FALSE;
		return;
	}
	if (right_loaded) {
		if(!check_images_match(&fit_left, &fit_right)) {
			siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: images do not match"),
				_("Image width, height, bit depth and number of channels must match"));
			clearfits(&fit_left);
			gtk_file_chooser_unselect_all(filechooser);
			left_loaded = FALSE;
			return;
		}
		else {
			left_loaded = TRUE;
			/* We are the second image to be loaded so we convert to the first
			* image's color profile */
			if (!profiles_identical(fit_left.icc_profile, fit_right.icc_profile)) {
				if (fit_right.icc_profile) {
					siril_colorspace_transform(&fit_left, fit_right.icc_profile);
					siril_log_message(_("Color profiles did not match: left-hand image has been converted to right-hand image color profile.\n"));
				} else {
					fit_right.icc_profile = copyICCProfile(fit_left.icc_profile);
					siril_log_message(_("Color profiles did not match: right-hand image has been assigned the left-hand image color profile.\n"));
				}
			}
			merge_fits_headers_to_result(&gfit, FALSE, &fit_left, &fit_right, NULL);
			if (fit_left.keywords.filter[0] != '\0' && fit_right.keywords.filter[0] != '\0' && strlen(fit_left.keywords.filter) >= 8 && strlen(fit_right.keywords.filter) >= 8) {
				gchar *temp_l = g_malloc(strlen(fit_left.keywords.filter) - 7);
				g_strlcpy(temp_l, fit_left.keywords.filter, strlen(fit_left.keywords.filter) - 8);
				gchar *temp_r = g_malloc(strlen(fit_right.keywords.filter) - 7);
				g_strlcpy(temp_r, fit_right.keywords.filter, strlen(fit_right.keywords.filter) - 8);
				if (strlen(fit_left.keywords.filter) == 8)
					temp_l[0] = '\0';
				if (strlen(fit_right.keywords.filter) == 8)
					temp_r[0] = '\0';
				if (!strcmp(temp_l, temp_r)) {
					if (fit_left.keywords.livetime >= fit_right.keywords.livetime) {
						gfit.keywords.livetime = fit_left.keywords.livetime;
						gfit.keywords.stackcnt = fit_left.keywords.stackcnt;
					} else {
						gfit.keywords.livetime = fit_right.keywords.livetime;
						gfit.keywords.stackcnt = fit_right.keywords.stackcnt;
					}
				}
				g_free(temp_l);
				g_free(temp_r);
			} // Otherwise the addition from merge_fits_headers_to_result is probably the best
			  // value to use as the image is probably being combined from 2 different stacks
		}
	} else {
		check_profile_correct(&fit_left);
		close_single_image();
		close_sequence(FALSE);
		clearfits(&gfit);
		copyfits(&fit_left, &gfit, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
		siril_log_message(_("Setting the output image ICC profile to the working color space.\n"));
		initialise_image();
		left_loaded = TRUE;
		clearfits(&fit_right_calc);
		copyfits(&fit_left, &fit_right_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
	}
	remix_histo_startup_left();
	GtkWidget *clip = lookup_widget("remixer_clip_mode_settings");
	gtk_widget_set_visible(clip, (fit_left.naxes[2] == 3));
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
	close_histograms(FALSE, TRUE);
	if (right_loaded) {
		clearfits(&fit_right);
		clearfits(&fit_right_calc);
		right_loaded = FALSE;
	}
	filename_right = siril_file_chooser_get_filename(filechooser);
	if (readfits(filename_right, &fit_right, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
				_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser);
		right_loaded = FALSE;
		return;
	}
	if (left_loaded) {
		if(!check_images_match(&fit_left, &fit_right)) {
			siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: images do not match"),
				_("Image width, height, bit depth and number of channels must match"));
			clearfits(&fit_right);
			gtk_file_chooser_unselect_all(filechooser);
			right_loaded = FALSE;
			return;
		}
		else {
			right_loaded = TRUE;
			/* We are the second image to be loaded so we convert to the first
			* image's color profile */
			if (!profiles_identical(fit_left.icc_profile, fit_right.icc_profile)) {
				if (fit_left.icc_profile) {
					siril_colorspace_transform(&fit_right, fit_left.icc_profile);
					siril_log_message(_("Color profiles did not match: right-hand image has been converted to left-hand image color profile.\n"));
				} else {
					fit_left.icc_profile = copyICCProfile(fit_right.icc_profile);
					siril_log_message(_("Color profiles did not match: left-hand image has been assigned the right-hand image color profile.\n"));
				}
			}
			if (!profiles_identical(fit_left.icc_profile, fit_right.icc_profile)) {
				siril_colorspace_transform(&fit_right, fit_left.icc_profile);
				siril_log_message(_("Color profiles did not match: right-hand image has been converted to left-hand image color profile.\n"));
			}
			merge_fits_headers_to_result(&gfit, FALSE, &fit_left, &fit_right, NULL);
			if (fit_left.keywords.filter[0] != '\0' && fit_right.keywords.filter[0] != '\0' && strlen(fit_left.keywords.filter) >= 8 && strlen(fit_right.keywords.filter) >= 8) {
				gchar *temp_l = g_malloc(strlen(fit_left.keywords.filter) - 7);
				g_strlcpy(temp_l, fit_left.keywords.filter, strlen(fit_left.keywords.filter) - 8);
				gchar *temp_r = g_malloc(strlen(fit_right.keywords.filter) - 7);
				g_strlcpy(temp_r, fit_right.keywords.filter, strlen(fit_right.keywords.filter) - 8);
				if (strlen(fit_left.keywords.filter) == 8)
					temp_l[0] = '\0';
				if (strlen(fit_right.keywords.filter) == 8)
					temp_r[0] = '\0';
				if (!strcmp(temp_l, temp_r)) {
					if (fit_left.keywords.livetime >= fit_right.keywords.livetime) {
						gfit.keywords.livetime = fit_left.keywords.livetime;
						gfit.keywords.stackcnt = fit_left.keywords.stackcnt;
					} else {
						gfit.keywords.livetime = fit_right.keywords.livetime;
						gfit.keywords.stackcnt = fit_right.keywords.stackcnt;
					}
				}
				g_free(temp_l);
				g_free(temp_r);
			} // Otherwise the addition from merge_fits_headers_to_result is probably the best
			  // value to use as the image is probably being combined from 2 different stacks
		}
	} else {
		check_profile_correct(&fit_right);
		close_single_image();
		close_sequence(FALSE);
		clearfits(&gfit);
		copyfits(&fit_right, &gfit, (CP_ALLOC | CP_COPYA | CP_FORMAT), 0);
		initialise_image();
		right_loaded = TRUE;
		clearfits(&fit_left_calc);
		copyfits(&fit_right, &fit_left_calc, (CP_ALLOC | CP_INIT | CP_FORMAT), 0);
	}
	remix_histo_startup_right();
	GtkWidget *clip = lookup_widget("remixer_clip_mode_settings");
	gtk_widget_set_visible(clip, (fit_right.naxes[2] == 3));
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
	if (!left_loaded) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no starless image loaded"),
				_("Load an image using the left-hand file chooser before attempting to set SP using the eyedropper"));
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
	leftSP = (float) ref;
	if (leftSP < leftLP) {
		gtk_spin_button_set_value(spin_remix_LP_left, (double) leftSP);
		rightLP = rightSP;
	}
	if (leftSP > leftHP) {
		gtk_spin_button_set_value(spin_remix_HP_left, (double) leftSP);
		rightHP = rightSP;
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_left")),leftSP);
	left_changed = TRUE;
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
	if (!right_loaded) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no starless image loaded"),
				_("Load an image using the right-hand file chooser before attempting to set SP using the eyedropper"));
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
	rightSP = (float) ref;
	if (rightSP < rightLP) {
		gtk_spin_button_set_value(spin_remix_LP_right, (double) rightSP);
		rightLP = rightSP;
	}
	if (rightSP > rightHP) {
		gtk_spin_button_set_value(spin_remix_HP_right, (double) rightSP);
		rightHP = rightSP;
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_remix_SP_right")),rightSP);
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
		gtk_widget_set_visible(lookup_widget("ghtStretchTypecontrols3"), TRUE);
		gtk_widget_set_visible(lookup_widget("ghtColourmodelcontrols2"), TRUE);
		gtk_widget_set_visible(lookup_widget("ghtLPcontrols2"), TRUE);
		gtk_widget_set_visible(lookup_widget("ghtHPcontrols2"), TRUE);
//		gtk_widget_set_visible(lookup_widget("ghtBPcontrols2"), TRUE);

		gtk_widget_set_visible(lookup_widget("ghtStretchTypecontrols1"), TRUE);
		gtk_widget_set_visible(lookup_widget("ghtColourmodelcontrols1"), TRUE);
		gtk_widget_set_visible(lookup_widget("ghtLPcontrols1"), TRUE);

//		gtk_widget_set_visible(lookup_widget("ghtBPcontrols1"), TRUE);
		gtk_button_set_label(GTK_BUTTON(lookup_widget("remix_advanced")), _("Simple"));
		gtk_widget_set_tooltip_text(lookup_widget("remix_advanced"), _("Show only simple stretch options. Note that any parameters that have been set while in advanced mode will remain set."));
		advanced_interface = TRUE;
	}
	else {
		gtk_widget_set_visible(lookup_widget("ghtStretchTypecontrols3"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtColourmodelcontrols2"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtLPcontrols2"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtHPcontrols2"), FALSE);
//		gtk_widget_set_visible(lookup_widget("ghtBPcontrols2"), FALSE);

		gtk_widget_set_visible(lookup_widget("ghtStretchTypecontrols1"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtColourmodelcontrols1"), FALSE);
		gtk_widget_set_visible(lookup_widget("ghtLPcontrols1"), FALSE);

//		gtk_widget_set_visible(lookup_widget("ghtBPcontrols1"), FALSE);
		gtk_button_set_label(GTK_BUTTON(lookup_widget("remix_advanced")), _("Advanced"));
		gtk_widget_set_tooltip_text(lookup_widget("remix_advanced"), _("Show advanced stretch options."));
		advanced_interface = FALSE;
	}
}

void on_toggle_remixer_log_histograms_toggled(GtkToggleButton *button, gpointer user_data) {
	remix_log_scale = gtk_toggle_button_get_active(button);
	update_remix_histo_left();
	update_remix_histo_right();
}

void on_dialog_star_remix_close() {
	apply_remix_cancel();
}

void on_remixer_help_okay_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("remixer_help");
}

void on_remix_histos_button_clicked(GtkButton *button, gpointer user_data) {
	siril_open_dialog("remixer_histos");
}

void on_remixer_clipmode_changed(GtkComboBox *combo, gpointer user_data) {
	left_changed = TRUE;
	right_changed = TRUE;
	clip_mode = (clip_mode_t) gtk_combo_box_get_active(combo);
	update_remix_histo_left();
	update_remix_histo_right();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = remixer_update_preview;
	param->show_preview = TRUE;
	notify_update((gpointer) param);
}
