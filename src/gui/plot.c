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

#include "plot.h"

#include <cairo.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "core/sleef.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/sequence_list.h"
#include "registration/registration.h"
#include "kplot.h"
#include "algos/PSF.h"
#include "io/ser.h"
#include "io/sequence.h"
#include "io/gnuplot_i.h"
#include "gui/PSF_list.h"
#include "opencv/opencv.h"

// TODO: Would probably be more efficient to cache plot surface and add selection as an overlay

#define XLABELSIZE 15
#define PLOT_SLIDER_THICKNESS 10. // thickness of the sliders in pixels
#define SIDE_MARGIN 12. // the margin in pixels top and left to define the start of selectable zone (allows to write info atop of selection)
#define SEL_TOLERANCE 3. // toerance in pixels for grabbing the selection borders

static GtkWidget *drawingPlot = NULL, *sourceCombo = NULL, *combo = NULL,
		*varCurve = NULL, *buttonClearAll = NULL,
		*buttonClearLatest = NULL, *arcsec = NULL, *julianw = NULL,
		*comboX = NULL, *layer_selector = NULL, *buttonSavePrt = NULL, *buttonSaveCSV = NULL;
static pldata *plot_data;
static struct kpair ref, curr;
static gboolean use_photometry = FALSE, requires_seqlist_update = FALSE;
static char *ylabel = NULL;
static gchar *xlabel = NULL;
static enum photometry_source photometry_selected_source = FWHM;
static enum registration_source registration_selected_source = r_FWHM;
static enum registration_source X_selected_source = r_FRAME;
static int julian0 = 0;
static int reglayer = 0;
static gboolean is_fwhm = TRUE;
static gnuplot_ctrl *gplot = NULL;
static gboolean is_arcsec = FALSE;
static gboolean force_Julian = FALSE;
static plot_draw_data_t pdd = { 0 };
static char *regfmt32[] = { "%0.2f", "%0.2f", "%0.2f", "%0.4f", "%0.0f", "%0.1f", "%0.1f", "%0.3f", "%0.0f" };
static char *regfmt16[] = { "%0.2f", "%0.2f", "%0.2f", "%0.0f", "%0.0f", "%0.1f", "%0.1f", "%0.3f", "%0.0f" };
static char *phtfmt32[] = { "%0.2f", "%0.2f", "%0.2f", "%0.2f", "%0.4f", "%0.1f", "%0.1f", "%0.2f"};
static char *phtfmt16[] = { "%0.2f", "%0.2f", "%0.2f", "%0.2f", "%0.0f", "%0.1f", "%0.1f", "%0.2f"};
static GtkMenu *menu = NULL;
static GtkMenuItem *menu_item1 = NULL, *menu_item2 = NULL, *menu_item3 = NULL;

static void formatX(double v, char *buf, size_t bufsz) {
	char *fmt;
	if (use_photometry && julian0 && force_Julian) {
		fmt = "%0.5f";
	} else {
		fmt = (gfit.type == DATA_FLOAT) ? regfmt32[X_selected_source] : regfmt16[X_selected_source];
	}
	// size of buf is 128
	// https://kristaps.bsd.lv/kplot/kplot_alloc.3.html
	snprintf(buf, 128, fmt, v);
}
static void formatY(double v, char *buf, size_t bufsz) {
	char *fmt;
	if (use_photometry) {
		fmt = (gfit.type == DATA_FLOAT) ? phtfmt32[photometry_selected_source] : phtfmt16[photometry_selected_source];
	} else {
		fmt = (gfit.type == DATA_FLOAT) ? regfmt32[registration_selected_source] : regfmt16[registration_selected_source];
	}
	snprintf(buf, 128, fmt, v);
}

static void update_ylabel();
static void set_colors(struct kplotcfg *cfg);
static void free_colors(struct kplotcfg *cfg);
void on_JulianPhotometry_toggled(GtkToggleButton *button, gpointer user_data);
void on_plotCombo_changed(GtkComboBox *box, gpointer user_data);
void on_plotComboX_changed(GtkComboBox *box, gpointer user_data);


static const gchar *photometry_labels[] = {
	N_("FWHM"),
	N_("Roundness"),
	N_("Amplitude"),
	N_("Magnitude"),
	N_("Background"),
	N_("X Position"),
	N_("Y Position"),
	N_("SNR")
};

static const gchar *registration_labels[] = {
	N_("FWHM"),
	N_("Roundness"),
	N_("wFWHM"),
	N_("Background"),
	N_("# Stars"),
	N_("X Position"),
	N_("Y Position"),
	N_("Quality"),
	N_("Frame"),
};

pldata *alloc_plot_data(int size) {
	pldata *plot = malloc(sizeof(pldata));
	plot->frame = calloc(size, sizeof(double));
	plot->julian = calloc(size, sizeof(double));
	plot->data = calloc(size, sizeof(struct kpair));
	plot->err = calloc(size, sizeof(struct kpair));
	if (!plot->frame || !plot->julian || !plot->data || !plot->err) {
		PRINT_ALLOC_ERR;
		free(plot->frame);
		free(plot->julian);
		free(plot->data);
		free(plot->err);
		free(plot);
		return NULL;
	}
	plot->nb = size;
	plot->next = NULL;
	return plot;
}


static void convert_surface_to_plot_coords(gdouble x, gdouble y, double *xpos, double *ypos) {
	*xpos = max(min(pdd.scale.x * (x - pdd.offset.x) + pdd.pdatamin.x, pdd.pdatamax.x), pdd.pdatamin.x);
	*ypos = max(min(pdd.scale.y * (pdd.range.y - y + pdd.offset.y) + pdd.pdatamin.y, pdd.pdatamax.y), pdd.pdatamin.y);
}

// static void convert_plot_to_surface_coords(double x, double y, double *xpos, double *ypos) {
// 	*xpos = 1./ pdd.scale.x * (x - pdd.pdatamin.x) + pdd.offset.x;
// 	*ypos = 1./ pdd.scale.y * (pdd.pdatamin.y - y) + pdd.range.y + pdd.offset.y;
// }

static void reset_plot_zoom() {
	pdd.xrange[0] = 0.;
	pdd.xrange[1] = 1.;
	pdd.yrange[0] = 0.;
	pdd.yrange[1] = 1.;
	pdd.marker_grabbed = MARKER_NONE;
	pdd.action = SELACTION_NONE;
	pdd.border_grabbed = SELBORDER_NONE;
	pdd.slider_grabbed = SLIDER_NONE;
	pdd.selection = (rectangled){0., 0., 0., 0.};
}

static gboolean selection_is_active() {
	return pdd.selection.w > 0. && pdd.selection.h > 0.;
}

static gboolean is_inside_selection(double x, double y) {
	if (!selection_is_active()) return FALSE;
	if (x >= pdd.selection.x + SEL_TOLERANCE &&
		x <= pdd.selection.x + pdd.selection.w - SEL_TOLERANCE &&
		y >= pdd.selection.y + SEL_TOLERANCE && 
		y <= pdd.selection.y + pdd.selection.h - SEL_TOLERANCE) 
			return TRUE;
	return FALSE;
}

static gboolean is_inside_grid(double x, double y) {
	if (x <= pdd.offset.x + pdd.range.x + SEL_TOLERANCE &&
		x >= pdd.offset.x - SEL_TOLERANCE &&
		y <= pdd.offset.y + pdd.range.y + SEL_TOLERANCE &&
		y >= pdd.offset.y - SEL_TOLERANCE)
			return TRUE;
	return FALSE;
}

static enum border_type is_over_selection_border(double x, double y) {
	if (!selection_is_active()) return SELBORDER_NONE;
	if (x >= pdd.selection.x && x <= pdd.selection.x + pdd.selection.w && y >= pdd.selection.y - SEL_TOLERANCE && y <= pdd.selection.y + SEL_TOLERANCE) return SELBORDER_TOP;
	if (x >= pdd.selection.x && x <= pdd.selection.x + pdd.selection.w && y >= pdd.selection.y + pdd.selection.h - SEL_TOLERANCE && y <= pdd.selection.y + pdd.selection.h + SEL_TOLERANCE) return SELBORDER_BOTTOM;
	if (y >= pdd.selection.y && y <= pdd.selection.y + pdd.selection.h && x >= pdd.selection.x - SEL_TOLERANCE && x <= pdd.selection.x + SEL_TOLERANCE) return SELBORDER_LEFT;
	if (y >= pdd.selection.y && y <= pdd.selection.y + pdd.selection.h && x >= pdd.selection.x + pdd.selection.w - SEL_TOLERANCE && x <= pdd.selection.x + pdd.selection.w + SEL_TOLERANCE) return SELBORDER_RIGHT;
	return SELBORDER_NONE;
}

// Returns TRUE if within all the plotting surface except
// the full-length bottom and right rectangles that enclose sliders
// This should avoid conflicts with sliders interactions
static gboolean is_inside_selectable_zone(double x, double y) {
	if (x <= pdd.surf_w - PLOT_SLIDER_THICKNESS && x >= SIDE_MARGIN && y <= pdd.surf_h - PLOT_SLIDER_THICKNESS && y >= SIDE_MARGIN) return TRUE;
	return FALSE;
}

static gboolean is_inside_slider(double x, double y, enum slider_type slider_t) {
	switch (slider_t) {
		// some margins are included to make sure we can grab
		case SLIDER_X:
			return (x <= pdd.offset.x + pdd.range.x + SEL_TOLERANCE && x >= pdd.offset.x - SEL_TOLERANCE && y >= pdd.surf_h - PLOT_SLIDER_THICKNESS);
		case SLIDER_Y:
			return (y <= pdd.offset.y + pdd.range.y + SEL_TOLERANCE && y >= pdd.offset.y - SEL_TOLERANCE && x >= pdd.surf_w - PLOT_SLIDER_THICKNESS);
		default:
			return FALSE;
	}
}

static gboolean is_over_marker(double x, double y, enum marker_type marker_t) {
	switch (marker_t) {
		case MARKER_X_MIN:
			return (fabs(x - (pdd.offset.x + pdd.range.x * pdd.xrange[0])) <= PLOT_SLIDER_THICKNESS * 0.5 && fabs(y - (pdd.surf_h - PLOT_SLIDER_THICKNESS * 0.5)) <= PLOT_SLIDER_THICKNESS * 0.5);
		case MARKER_X_MAX:
			return (fabs(x - (pdd.offset.x + pdd.range.x * pdd.xrange[1])) <= PLOT_SLIDER_THICKNESS * 0.5 && fabs(y - (pdd.surf_h - PLOT_SLIDER_THICKNESS * 0.5)) <= PLOT_SLIDER_THICKNESS * 0.5);
		case MARKER_Y_MIN:
			return (fabs(x - (pdd.surf_w - PLOT_SLIDER_THICKNESS * 0.5)) <= PLOT_SLIDER_THICKNESS * 0.5 && fabs(y - (pdd.offset.y + pdd.range.y * (1. - pdd.yrange[0]))) <= PLOT_SLIDER_THICKNESS * 0.5);
		case MARKER_Y_MAX:
			return (fabs(x - (pdd.surf_w - PLOT_SLIDER_THICKNESS * 0.5)) <= PLOT_SLIDER_THICKNESS * 0.5 && fabs(y - (pdd.offset.y + pdd.range.y * (1. - pdd.yrange[1]))) <= PLOT_SLIDER_THICKNESS * 0.5);
		default:
			return FALSE;
	}
}

static int get_closest_marker(double x, double y, enum slider_type slider_t, double *valrange) {
	double val = (slider_t == SLIDER_X) ? (x - pdd.offset.x) / pdd.range.x : (pdd.offset.y + pdd.range.y - y) / pdd.range.y;
	// should not be outside of [0, 1] but just in case
	val = min(1., val);
	val = max(0., val);
	return (int)(fabs(val - valrange[0]) > fabs(val - valrange[1]));
}

static void find_range_from_pos(double x, double y, enum slider_type slider_t, int index, double *valrange) {
	double val = (slider_t == SLIDER_X) ? (x - pdd.offset.x) / pdd.range.x : (pdd.offset.y + pdd.range.y - y) / pdd.range.y;
	val = min(1., val);
	val = max(0., val);
	double otherval = (index == 0) ? valrange[1] : valrange[0];
	// making sure the values are correctly ordered
	valrange[0] = min(val, otherval);
	valrange[1] = max(val, otherval);
}

static double find_rangeval_from_pos(double v, enum slider_type slider_t) {
	return (slider_t == SLIDER_X) ? (v - pdd.offset.x) / pdd.range.x : (pdd.offset.y + pdd.range.y - v) / pdd.range.y;
}

static void update_slider(enum slider_type slider_t, double valmin, double valmax) {
	switch (slider_t) {
		case SLIDER_X:
			pdd.xrange[0] = valmin;
			pdd.xrange[1] = valmax;
			break;
		case SLIDER_Y:
			pdd.yrange[0] = valmin;
			pdd.yrange[1] = valmax;
			break;
		default:
			break;
	}
	drawPlot();
}

static gboolean get_index_of_frame(double x, double y, gboolean check_index_incl, double *index, double *xpos, double *ypos) {
	int closestframe = -1;
	double mindist = DBL_MAX;
	pldata *plot = plot_data;
	convert_surface_to_plot_coords(x, y, xpos, ypos);
	// double testx, testy;
	// convert_plot_to_surface_coords(pdd.datamin.x, pdd.datamin.y, &testx, &testy);

	double invrangex = 1./(pdd.pdatamax.x - pdd.pdatamin.x);
	double invrangey = 1./(pdd.pdatamax.y - pdd.pdatamin.y);
	*index = *xpos;

	while (plot) {
		for (int j = 0; j < plot->nb; j++) {
			double dist = xpow((*index - plot->data[j].x) * invrangex, 2) + xpow((*ypos - plot->data[j].y) * invrangey, 2);
			if (dist < mindist) {
				mindist = dist;
				closestframe = plot->frame[j];
			}
		}
		plot = plot->next;
	}
	*index = (mindist < 0.0004) ? closestframe : -1; // only set index if distance between cursor and a point is small enough (2% of scales)

	if (check_index_incl && (*index >= 0 && *index <= pdd.pdatamax.x)) return com.seq.imgparam[(int)*index - 1].incl;
	return TRUE;

}

static void plot_draw_all_sliders(cairo_t *cr) {
	double color = (com.pref.gui.combo_theme == 0) ? 1.0 : 0.0;
	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, color, color, color);
	// x-slider
	cairo_rectangle(cr, pdd.offset.x, pdd.surf_h - PLOT_SLIDER_THICKNESS, pdd.range.x, PLOT_SLIDER_THICKNESS);
	cairo_stroke(cr);
	// y-slider
	cairo_rectangle(cr, pdd.surf_w - PLOT_SLIDER_THICKNESS, pdd.offset.y, PLOT_SLIDER_THICKNESS, pdd.range.y);
	cairo_stroke(cr);
}

static void plot_draw_slider_fill(cairo_t *cr, enum slider_type slider_t) {
	double color = 0.5;
	cairo_set_source_rgb(cr, color, color, color);
	switch (slider_t) {
		default:
		case SLIDER_X:
			cairo_rectangle(cr, pdd.offset.x + pdd.range.x * pdd.xrange[0], pdd.surf_h - PLOT_SLIDER_THICKNESS, pdd.range.x * (pdd.xrange[1] - pdd.xrange[0]), PLOT_SLIDER_THICKNESS);
			break;
		case SLIDER_Y:
			cairo_rectangle(cr, pdd.surf_w - PLOT_SLIDER_THICKNESS, pdd.offset.y + pdd.range.y * (1. - pdd.yrange[1]), PLOT_SLIDER_THICKNESS, pdd.range.y * (pdd.yrange[1] - pdd.yrange[0]));
			break;
	}
	cairo_fill(cr);
}

static void plot_draw_all_sliders_fill(cairo_t *cr) {
	plot_draw_slider_fill(cr, SLIDER_X);
	plot_draw_slider_fill(cr, SLIDER_Y);
}


static void plot_draw_marker(cairo_t *cr, enum marker_type marker_t) {
	cairo_set_source_rgb(cr, 0.8, 0., 0.);
	switch (marker_t) {
		default:
		case MARKER_X_MIN:
			cairo_arc(cr, pdd.offset.x + pdd.range.x * pdd.xrange[0], pdd.surf_h - PLOT_SLIDER_THICKNESS * 0.5, PLOT_SLIDER_THICKNESS * 0.5, 0., 2. * M_PI);
			break;
		case MARKER_X_MAX:
			cairo_arc(cr, pdd.offset.x + pdd.range.x * pdd.xrange[1], pdd.surf_h - PLOT_SLIDER_THICKNESS * 0.5, PLOT_SLIDER_THICKNESS * 0.5, 0., 2. * M_PI);
			break;
		case MARKER_Y_MIN:
			cairo_arc(cr, pdd.surf_w - PLOT_SLIDER_THICKNESS * 0.5, pdd.offset.y + pdd.range.y * (1. - pdd.yrange[0]), PLOT_SLIDER_THICKNESS * 0.5, 0., 2. * M_PI);
			break;
		case MARKER_Y_MAX:
			cairo_arc(cr, pdd.surf_w - PLOT_SLIDER_THICKNESS * 0.5, pdd.offset.y + pdd.range.y * (1. - pdd.yrange[1]), PLOT_SLIDER_THICKNESS * 0.5, 0., 2. * M_PI);
			break;
	}
	cairo_fill(cr);
}

static void plot_draw_all_markers(cairo_t *cr) {
	for (int i = MARKER_X_MIN; i <= MARKER_Y_MAX; i++)
		plot_draw_marker(cr, i);
}

static void plot_draw_selection(cairo_t *cr){
	if (pdd.selection.h == 0. || pdd.selection.w == 0.) return;
	double dash_format[] = { 4.0, 2.0 };
	double color = (com.pref.gui.combo_theme == 0) ? 0.8 : 0.2;
	cairo_set_source_rgb(cr, color, color, color);
	cairo_set_dash(cr, dash_format, 2, 0);
	cairo_set_line_width(cr, 1.);
	cairo_rectangle(cr, pdd.selection.x,  pdd.selection.y,
						 pdd.selection.w,  pdd.selection.h);
	cairo_stroke(cr);
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_move_to(cr, pdd.selection.x,  pdd.selection.y - 2);
	double xmin, ymin, xmax, ymax;
	convert_surface_to_plot_coords(pdd.selection.x, pdd.selection.y, &xmin, &ymax);
	convert_surface_to_plot_coords(pdd.selection.x + pdd.selection.w, pdd.selection.y + pdd.selection.h, &xmax, &ymin);
	gchar fmt[256] = { 0 };
	gchar buffer[256] = { 0 };
	if (use_photometry) {
		if (julian0 && force_Julian) {
			g_sprintf(fmt, "Nb: %s (for V star) - %s: [ %s , %s ] - %s: [ %s , %s ]", "\%d", xlabel,
			"%0.5f",
			"%0.5f",
			ylabel,
			(gfit.type == DATA_FLOAT) ? phtfmt32[photometry_selected_source] : phtfmt16[photometry_selected_source],
			(gfit.type == DATA_FLOAT) ? phtfmt32[photometry_selected_source] : phtfmt16[photometry_selected_source]);
		} else {
			g_sprintf(fmt, "Nb: %s (for V star) - %s: [ %s , %s ] - %s: [ %s , %s ]", "\%d", xlabel,
			(gfit.type == DATA_FLOAT) ? regfmt32[X_selected_source] : regfmt16[X_selected_source],
			(gfit.type == DATA_FLOAT) ? regfmt32[X_selected_source] : regfmt16[X_selected_source],
			ylabel,
			(gfit.type == DATA_FLOAT) ? phtfmt32[photometry_selected_source] : phtfmt16[photometry_selected_source],
			(gfit.type == DATA_FLOAT) ? phtfmt32[photometry_selected_source] : phtfmt16[photometry_selected_source]);
		}
	} else {
		g_sprintf(fmt, "Nb: %s - %s: [ %s , %s ] - %s: [ %s , %s ]", "\%d", xlabel,
		(gfit.type == DATA_FLOAT) ? regfmt32[X_selected_source] : regfmt16[X_selected_source],
		(gfit.type == DATA_FLOAT) ? regfmt32[X_selected_source] : regfmt16[X_selected_source],
		ylabel,
		(gfit.type == DATA_FLOAT) ? regfmt32[registration_selected_source] : regfmt16[registration_selected_source],
		(gfit.type == DATA_FLOAT) ? regfmt32[registration_selected_source] : regfmt16[registration_selected_source]);
	}
	g_sprintf(buffer, fmt, pdd.nbselected, xmin, xmax, ymin, ymax);
	cairo_show_text(cr, buffer);
	cairo_stroke(cr);
}

static void build_registration_dataset(sequence *seq, int layer, int ref_image,
		pldata *plot) {
	int i, j;
	double fwhm;
	double dx, dy;
	double cx,cy;
	gboolean Href_is_invalid;
	cx = (seq->is_variable) ? (double)seq->imgparam[ref_image].rx * 0.5 : (double)seq->rx * 0.5;
	cy = (seq->is_variable) ? (double)seq->imgparam[ref_image].ry * 0.5 : (double)seq->ry * 0.5;
	Homography Href = seq->regparam[layer][ref_image].H;
	Href_is_invalid = (guess_transform_from_H(Href) == -2) ? TRUE : FALSE;
	pdd.datamin = (point){ DBL_MAX, DBL_MAX};
	pdd.datamax = (point){ -DBL_MAX, -DBL_MAX};

	for (i = 0, j = 0; i < plot->nb; i++) {
		if (!seq->imgparam[i].incl)
			continue;
		switch (X_selected_source) {
			case r_ROUNDNESS:
				plot->data[j].x = seq->regparam[layer][i].roundness;
				break;
			case r_FWHM:
				if (is_arcsec) {
					double bin = gfit.unbinned ? (double) gfit.binning_x : 1.0;
					convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][i].fwhm, bin, (double) gfit.pixel_size_x, gfit.focal_length, &fwhm);
				} else {
					fwhm = seq->regparam[layer][i].fwhm;
				}
				plot->data[j].x = fwhm;
				break;
			case r_X_POSITION:
			case r_Y_POSITION:
				// compute the center of image i in the axes of the reference frame
				dx = (seq->is_variable) ? (double)seq->imgparam[i].rx * 0.5 : (double)seq->rx * 0.5;
				dy = (seq->is_variable) ? (double)seq->imgparam[i].ry * 0.5 : (double)seq->ry * 0.5;
				if (Href_is_invalid || guess_transform_from_H(seq->regparam[layer][i].H) == -2) {
					plot->data[j].x = 0;
					break;
				}
				cvTransfPoint(&dx, &dy, seq->regparam[layer][i].H, Href);
				plot->data[j].x = (X_selected_source == r_X_POSITION) ? dx - cx : dy - cy;
				break;
			case r_WFWHM:
				if (is_arcsec) {
					double bin = gfit.unbinned ? (double) gfit.binning_x : 1.0;
					convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][i].weighted_fwhm, bin, (double) gfit.pixel_size_x, gfit.focal_length, &fwhm);
				} else {
					fwhm = seq->regparam[layer][i].weighted_fwhm;
				}
				plot->data[j].x = fwhm;
				break;
			case r_QUALITY:
				plot->data[j].x = seq->regparam[layer][i].quality;
				break;
			case r_BACKGROUND:
				plot->data[j].x = seq->regparam[layer][i].background_lvl;
				break;
			case r_NBSTARS:
				plot->data[j].x = seq->regparam[layer][i].number_of_stars;
				break;
			case r_FRAME:
				plot->data[j].x = (double) i + 1;
				break;
			default:
				break;
		}
		plot->data[j].x = (isnan(plot->data[j].x)) ? 0.0 : plot->data[j].x;
		switch (registration_selected_source) {
			case r_ROUNDNESS:
				plot->data[j].y = seq->regparam[layer][i].roundness;
				break;
			case r_FWHM:
				if (is_arcsec) {
					double bin = gfit.unbinned ? (double) gfit.binning_x : 1.0;
					convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][i].fwhm, bin, (double) gfit.pixel_size_x, gfit.focal_length, &fwhm);
				} else {
					fwhm = seq->regparam[layer][i].fwhm;
				}
				plot->data[j].y = fwhm;
				break;
			case r_X_POSITION:
			case r_Y_POSITION:
				// compute the center of image i in the axes of the reference frame
				dx = (seq->is_variable) ? (double)seq->imgparam[i].rx * 0.5 : (double)seq->rx * 0.5;
				dy = (seq->is_variable) ? (double)seq->imgparam[i].ry * 0.5 : (double)seq->ry * 0.5;
				if (Href_is_invalid || guess_transform_from_H(seq->regparam[layer][i].H) == -2) {
					plot->data[j].y = 0;
					break;
				}
				cvTransfPoint(&dx, &dy, seq->regparam[layer][i].H, Href);
				plot->data[j].y = (registration_selected_source == r_X_POSITION) ? dx - cx : dy - cy;
				break;
			case r_WFWHM:
				if (is_arcsec) {
					double bin = gfit.unbinned ? (double) gfit.binning_x : 1.0;
					convert_single_fwhm_to_arcsec_if_possible(seq->regparam[layer][i].weighted_fwhm, bin, (double) gfit.pixel_size_x, gfit.focal_length, &fwhm);
				} else {
					fwhm = seq->regparam[layer][i].weighted_fwhm;
				}
				plot->data[j].y = fwhm;
				break;
			case r_QUALITY:
				plot->data[j].y = seq->regparam[layer][i].quality;
				break;
			case r_BACKGROUND:
				plot->data[j].y = seq->regparam[layer][i].background_lvl;
				break;
			case r_NBSTARS:
				plot->data[j].y = seq->regparam[layer][i].number_of_stars;
				break;
			default:
				break;
		}
		plot->data[j].y = (isnan(plot->data[j].y)) ? 0.0 : plot->data[j].y;
		plot->frame[j] =  (double) i + 1;
		if (i == ref_image) {
			ref.x = plot->data[j].x;
			ref.y = plot->data[j].y;
		}
		if ((i == seq->current) & com.seq.imgparam[i].incl) {
			curr.x = plot->data[j].x;
			curr.y = plot->data[j].y;
		}
		// caching the data range
		if (pdd.datamin.x > plot->data[j].x) pdd.datamin.x = plot->data[j].x;
		if (pdd.datamax.x < plot->data[j].x) pdd.datamax.x = plot->data[j].x;
		if (pdd.datamin.y > plot->data[j].y) pdd.datamin.y = plot->data[j].y;
		if (pdd.datamax.y < plot->data[j].y) pdd.datamax.y = plot->data[j].y;
		j++;
	}
	plot->nb = j;
	pdd.pdatamin.x = pdd.datamin.x + (pdd.datamax.x - pdd.datamin.x) * pdd.xrange[0];
	pdd.pdatamax.x = pdd.datamin.x + (pdd.datamax.x - pdd.datamin.x) * pdd.xrange[1];
	pdd.pdatamin.y = pdd.datamin.y + (pdd.datamax.y - pdd.datamin.y) * pdd.yrange[0];
	pdd.pdatamax.y = pdd.datamin.y + (pdd.datamax.y - pdd.datamin.y) * pdd.yrange[1];
}

static void set_x_photometry_values(sequence *seq, pldata *plot, int image_index, int point_index) {
	if (seq->imgparam[image_index].date_obs) {
		double julian;
		GDateTime *tsi = g_date_time_ref(seq->imgparam[image_index].date_obs);
		if (seq->exposure > 0.0) {
			GDateTime *new_dt = g_date_time_add_seconds(tsi, seq->exposure * 0.5);
			julian = date_time_to_Julian(new_dt);
			g_date_time_unref(new_dt);
		} else {
			julian = date_time_to_Julian(tsi);
		}

		plot->julian[point_index] = julian - (double)julian0;

		g_date_time_unref(tsi);
	} else {
		plot->julian[point_index] = (double) image_index + 1; // should not happen
		siril_debug_print("no DATE-OBS information for frame %d\n", image_index);
	}
	plot->frame[point_index] = (double) image_index + 1;

	if (julian0 && force_Julian) {
		plot->data[point_index].x = plot->julian[point_index];
	} else {
		plot->data[point_index].x = plot->frame[point_index];
	}

	plot->err[point_index].x = plot->data[point_index].x;
}

static void build_photometry_dataset(sequence *seq, int dataset, int size,
		int ref_image, pldata *plot) {
	int i, j;
	double offset = -1001.0;
	double fwhm;
	psf_star **psfs = seq->photometry[dataset], *ref_psf;
	if (seq->reference_star >= 0 && !seq->photometry[seq->reference_star])
		seq->reference_star = -1;

	for (i = 0, j = 0; i < size; i++) {
		if (!seq->imgparam[i].incl || !psfs[i])
			continue;
		if (!julian0 && !xlabel) {
			if (seq->imgparam[i].date_obs) {
				GDateTime *ts0 = g_date_time_ref(seq->imgparam[i].date_obs);
				if (seq->exposure > 0.0) {
					GDateTime *new_dt = g_date_time_add_seconds(ts0, seq->exposure * 0.5);
					julian0 = (int) date_time_to_Julian(new_dt);
					g_date_time_unref(new_dt);
				} else {
					julian0 = (int) date_time_to_Julian(ts0);
				}
				g_date_time_unref(ts0);
				//siril_debug_print("julian0 set to %d\n", julian0);
			}
			if (julian0 && force_Julian) {
				xlabel = malloc(XLABELSIZE * sizeof(char));
				g_snprintf(xlabel, XLABELSIZE, "(JD) %d +", julian0);
			} else {
				xlabel = _("Frames");
			}
		}
		set_x_photometry_values(seq, plot, i, j);

		switch (photometry_selected_source) {
			case ROUNDNESS:
				plot->data[j].y = psfs[i]->fwhmy / psfs[i]->fwhmx;
				break;
			case FWHM:
				if (is_arcsec) {
					fwhm_to_arcsec_if_needed(&gfit, psfs[i]);
					fwhm = psfs[i]->fwhmx_arcsec < 0 ? psfs[i]->fwhmx : psfs[i]->fwhmx_arcsec;
				} else {
					fwhm_to_pixels(psfs[i]);
					fwhm = psfs[i]->fwhmx;
				}
				plot->data[j].y = fwhm;
				break;
			case AMPLITUDE:
				plot->data[j].y = psfs[i]->A;
				break;
			case MAGNITUDE:
				if (!psfs[i]->phot_is_valid)
					continue;
				plot->data[j].y = psfs[i]->mag;
				plot->err[j].y = psfs[i]->s_mag;

				if (seq->reference_star >= 0) {
					/* we have a reference star for the sequence,
					 * with photometry data */
					ref_psf = seq->photometry[seq->reference_star][i];
					if (ref_psf)
						offset = seq->reference_mag - ref_psf->mag;
				} else if (com.magOffset > 0.0)
					offset = com.magOffset;

				/* apply the absolute apparent magnitude offset */
				if (offset > -1000.0)
					plot->data[j].y += offset;
				break;
			case BACKGROUND:
				plot->data[j].y = psfs[i]->B;
				break;
			case X_POSITION:
				plot->data[j].y = psfs[i]->xpos;
				break;
			case Y_POSITION:
				plot->data[j].y = psfs[i]->ypos;
				break;
			case SNR:
				plot->data[j].y = psfs[i]->SNR;
				break;
			default:
				break;
		}

		/* we'll just take the reference image point from the last data set rendered */
		if (i == ref_image) {
			ref.x = plot->data[j].x;
			ref.y = plot->data[j].y;
		}
		if (i == seq->current && com.seq.imgparam[i].incl) {
			curr.x = plot->data[j].x;
			curr.y = plot->data[j].y;
		}
		// caching the data range
		if (pdd.datamin.x > plot->data[j].x) pdd.datamin.x = plot->data[j].x;
		if (pdd.datamax.x < plot->data[j].x) pdd.datamax.x = plot->data[j].x;
		if (pdd.datamin.y > plot->data[j].y) pdd.datamin.y = plot->data[j].y;
		if (pdd.datamax.y < plot->data[j].y) pdd.datamax.y = plot->data[j].y;
		j++;
	}
	plot->nb = j;
	pdd.pdatamin.x = pdd.datamin.x + (pdd.datamax.x - pdd.datamin.x) * pdd.xrange[0];
	pdd.pdatamax.x = pdd.datamin.x + (pdd.datamax.x - pdd.datamin.x) * pdd.xrange[1];
	pdd.pdatamin.y = pdd.datamin.y + (pdd.datamax.y - pdd.datamin.y) * pdd.yrange[0];
	pdd.pdatamax.y = pdd.datamin.y + (pdd.datamax.y - pdd.datamin.y) * pdd.yrange[1];
}

static double get_error_for_time(pldata *plot, double time) {
	/* when data are sorted we need to check order by matching
	 * timestamps in order to sort uncertainties as well
	 */
	for (int k = 0; k < plot->nb; k++) {
		if (plot->err[k].x == time)
			return plot->err[k].y;
	}
	return 0.0;
}

// call after having filled the plot data of the multiple stars with either
// generate_magnitude_data() from the command or build_photometry_dataset() from the GUI
// the first will be the target
int light_curve(pldata *plot, sequence *seq, gchar *filename) {
	int i, j, nbImages = 0;
	double *vmag, *err, *x, *real_x;
	gboolean use_gnuplot = gnuplot_is_available();
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the light curve data will be produced in %s but no image will be created.\n"), filename);
	}
	if (!seq->photometry[0]) {
		siril_log_color_message(_("No photometry data found, error\n"), "red");
		return -1;
	}

	/* get number of valid frames for each star */
	int ref_valid_count[MAX_SEQPSF] = { 0 };
	gboolean ref_valid[MAX_SEQPSF] = { FALSE };
	for (i = 0; i < plot->nb; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)
			continue;
		++nbImages;
		for (int r = 1; r < MAX_SEQPSF && seq->photometry[r]; r++) {
			if (seq->photometry[r][i] && seq->photometry[r][i]->phot_is_valid)
				ref_valid_count[r]++;
		}
	}
	siril_debug_print("we have %d images with a valid photometry for the variable star\n", nbImages);
	if (nbImages < 1)
		return -1;

	int nb_ref_stars = 0;
	// select reference stars that are only available at least 3/4 of the time
	for (int r = 1; r < MAX_SEQPSF && seq->photometry[r]; r++) {
		ref_valid[r] = ref_valid_count[r] >= nbImages * 3 / 4;
		siril_debug_print("reference star %d has %d/%d valid measures, %s\n", r, ref_valid_count[r], nbImages, ref_valid[r] ? "including" : "discarding");
		if (ref_valid[r])
			nb_ref_stars++;
	}

	if (nb_ref_stars == 0) {
		siril_log_color_message(_("The reference stars are not good enough, probably out of the configured valid pixel range, cannot calibrate the light curve\n"), "red");
		return -1;
	}
	if (nb_ref_stars < 1)
		siril_log_color_message(_("Only one reference star was validated, this will not result in an accurate light curve. Try to add more reference stars or check the configured valid pixel range\n"), "salmon");
	else siril_log_message(_("Using %d stars to calibrate the light curve\n"), nb_ref_stars);

	vmag = calloc(nbImages, sizeof(double));
	err = calloc(nbImages, sizeof(double));
	x = calloc(nbImages, sizeof(double));
	real_x = calloc(nbImages, sizeof(double));
	if (!vmag || !err || !x || !real_x) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	// i is index in dataset, j is index in output
	for (i = 0, j = 0; i < plot->nb; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)
			continue;
		double cmag = 0.0, cerr = 0.0;

		x[j] = plot->data[j].x;			// relative date
		real_x[j] = x[j] + (double)julian0;	// absolute date
		vmag[j] = plot->data[j].y;		// magnitude
		err[j] = get_error_for_time(plot, x[j]);// error of the magnitude

		int nb_ref = 0;
		pldata *cur_plot = plot->next;
		/* First data plotted are variable data, others are references
		 * Variable is done above, now we compute references */
		for (int r = 1; r < MAX_SEQPSF && seq->photometry[r]; r++) {
			if (ref_valid[r] && seq->photometry[r][i] && seq->photometry[r][i]->phot_is_valid) {
				/* variable data, inversion of Pogson's law
				 * Flux = 10^(-0.4 * mag)
				 */
				cmag += pow(10, -0.4 * cur_plot->data[j].y);
				cerr += get_error_for_time(cur_plot, x[j]);
				++nb_ref;
			}
			cur_plot = cur_plot->next;
		}
		/* Converting back to magnitude */
		if (nb_ref == nb_ref_stars) {
			/* we consider an image to be invalid if all references are not valid,
			 * because it changes the mean otherwise and makes nonsense data */
			cmag = -2.5 * log10(cmag / nb_ref);
			cerr = (cerr / nb_ref) / sqrt((double) nb_ref);

			vmag[j] = vmag[j] - cmag;
			err[j] = fmin(9.999, sqrt(err[j] * err[j] + cerr * cerr));
			j++;
		}
	}
	int nb_valid_images = j;

	siril_log_message(_("Calibrated data for %d points of the light curve, %d excluded because of invalid calibration\n"), nb_valid_images, plot->nb - nb_valid_images);

	/*  data are computed, now plot the graph. */

	if (use_gnuplot) {
		/* First, close the graph if already exists */
		if (gplot != NULL) {
			gnuplot_close(gplot);
		}

		if ((gplot = gnuplot_init())) {
			/* Plotting light curve */
			gnuplot_set_title(gplot, _("Light Curve"));
			gnuplot_set_xlabel(gplot, xlabel);
			gnuplot_reverse_yaxis(gplot);
			gnuplot_setstyle(gplot, "errorbars");
			gnuplot_plot_xyyerr(gplot, x, vmag, err, nb_valid_images, "", 0);
		}
		else siril_log_message(_("Communicating with gnuplot failed, still creating the data file\n"));
	}

	/* Exporting data in a dat file */
	if ((gnuplot_write_xyyerr_dat(filename, real_x, vmag, err, nb_valid_images, "JD_UT V-C err"))) {
		if (com.script)
			siril_log_color_message(_("Failed to create the light curve data file %s\n"), "red", filename);
		else siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Something went wrong while saving plot"));
	} else {
		siril_log_message(_("%s has been saved.\n"), filename);
	}

	free(vmag);
	free(err);
	free(x);
	free(real_x);
	return 0;
}

static int exportCSV(pldata *plot, sequence *seq, gchar *filename) {
	GError *error = NULL;

	GFile *file = g_file_new_for_path(filename);
	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			fprintf(stderr, "exportCSV: Cannot export\n");
		}
		g_object_unref(file);
		return 1;
	}

	if (use_photometry) {
		pldata *tmp_plot = plot;
		for (int i = 0, j = 0; i < plot->nb; i++) {
			if (!seq->imgparam[i].incl)
				continue;
			int x = 0;
			double date = tmp_plot->data[j].x;
			if (julian0 && force_Julian) {
				date += julian0;
			}
			gchar *buffer = g_strdup_printf("%.10lf", date);
			if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
				g_warning("%s\n", error->message);
				g_free(buffer);
				g_clear_error(&error);
				g_object_unref(output_stream);
				g_object_unref(file);
				return 1;
			}
			g_free(buffer);
			while (x < MAX_SEQPSF && seq->photometry[x]) {
				buffer = g_strdup_printf(", %g", tmp_plot->data[j].y);
				if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
					g_warning("%s\n", error->message);
					g_free(buffer);
					g_clear_error(&error);
					g_object_unref(output_stream);
					g_object_unref(file);
					return 1;
				}
				tmp_plot = tmp_plot->next;
				++x;
				g_free(buffer);
			}
			if (!g_output_stream_write_all(output_stream, "\n", 1, NULL, NULL, &error)) {
				g_warning("%s\n", error->message);
				g_free(buffer);
				g_clear_error(&error);
				g_object_unref(output_stream);
				g_object_unref(file);
				return 1;
			}
			tmp_plot = plot;
			j++;
		}
	} else {
		for (int i = 0, j = 0; i < plot->nb; i++) {
			if (!seq->imgparam[i].incl)
				continue;
			double date = plot->data[j].x;
			if (julian0) {
				date += julian0;
			}
			gchar *buffer = g_strdup_printf("%.10lf, %g\n", date, plot->data[j].y);
			if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
				g_warning("%s\n", error->message);
				g_free(buffer);
				g_clear_error(&error);
				g_object_unref(output_stream);
				g_object_unref(file);
				return 1;
			}
			j++;
			g_free(buffer);
		}
	}
	siril_log_message(_("%s has been saved.\n"), filename);
	g_object_unref(output_stream);
	g_object_unref(file);

	return 0;
}

void free_plot_data() {
	pldata *plot = plot_data;
	while (plot) {
		pldata *next = plot->next;
		free(plot->julian);
		free(plot->frame);
		free(plot->data);
		free(plot->err);
		free(plot);
		plot = next;
	}
	julian0 = 0;
	xlabel = NULL;
	plot_data = NULL;
	free(pdd.selected);
	pdd.selected = NULL;
}

static void set_sensitive(GtkCellLayout *cell_layout,
		GtkCellRenderer *cell,
		GtkTreeModel *tree_model,
		GtkTreeIter *iter,
		gpointer data) {
	gboolean sensitive = TRUE;

	if (!use_photometry) {
		GtkTreePath* path = gtk_tree_model_get_path(tree_model, iter);
		if (!path) return;
		gint *index = gtk_tree_path_get_indices(path); // search by index to avoid translation problems
		if (index) {
			if (!is_fwhm) {
				sensitive = (index[0] == r_FRAME || index[0] == r_QUALITY ||
						index[0] == r_X_POSITION || index[0] == r_Y_POSITION);
			} else {
				sensitive = (index[0] == r_FRAME || index[0] == r_FWHM ||
						index[0] == r_WFWHM || index[0] == r_ROUNDNESS ||
						index[0] == r_BACKGROUND || index[0] == r_NBSTARS ||
						index[0] == r_X_POSITION || index[0] == r_Y_POSITION);
			}
		}
		gtk_tree_path_free(path);
	}
	g_object_set(cell, "sensitive", sensitive, NULL);
}

static void fill_plot_statics() {
	if (drawingPlot == NULL) {
		drawingPlot = lookup_widget("DrawingPlot");
		combo = lookup_widget("plotCombo");
		comboX = lookup_widget("plotComboX");
		varCurve = lookup_widget("varCurvePhotometry");
		buttonSaveCSV = lookup_widget("ButtonSaveCSV");
		buttonSavePrt = lookup_widget("ButtonSavePlot");
		arcsec = lookup_widget("arcsecPhotometry");
		julianw = lookup_widget("JulianPhotometry");
		sourceCombo = lookup_widget("plotSourceCombo");
		buttonClearAll = lookup_widget("clearAllPhotometry");
		buttonClearLatest = lookup_widget("clearLastPhotometry");
		layer_selector = lookup_widget("seqlist_dialog_combo");
	}
}

static void validate_combos() {
	fill_plot_statics();
	use_photometry = gtk_combo_box_get_active(GTK_COMBO_BOX(sourceCombo));
	if (!use_photometry) {
		reglayer = gtk_combo_box_get_active(GTK_COMBO_BOX(layer_selector));
		if (!(com.seq.regparam) || !(com.seq.regparam[reglayer]))
			reglayer = get_registration_layer(&com.seq);
	}
	gtk_widget_set_visible(varCurve, use_photometry);
	gtk_widget_set_visible(buttonSaveCSV, TRUE);
	gtk_widget_set_visible(buttonSavePrt, TRUE);
	g_signal_handlers_block_by_func(julianw, on_JulianPhotometry_toggled, NULL);
	gtk_widget_set_visible(julianw, use_photometry);
	g_signal_handlers_unblock_by_func(julianw, on_JulianPhotometry_toggled, NULL);

	g_signal_handlers_block_by_func(combo, on_plotCombo_changed, NULL);
	g_signal_handlers_block_by_func(comboX, on_plotComboX_changed, NULL);
	gtk_combo_box_text_remove_all(GTK_COMBO_BOX_TEXT(combo));
	int i = 0;
	if (use_photometry) {
		while (i < (G_N_ELEMENTS(photometry_labels))) {
			gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo), _(photometry_labels[i]));
			i++;
		}
		gtk_combo_box_set_active(GTK_COMBO_BOX(comboX), r_FRAME);
		X_selected_source = r_FRAME;
		gtk_widget_set_sensitive(comboX, FALSE);
	} else {
		while (i < (G_N_ELEMENTS(registration_labels)) - 1) { // do not write 'Frame' as possible Y value
			gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo), _(registration_labels[i]));
			i++;
		}
		gtk_combo_box_set_active(GTK_COMBO_BOX(sourceCombo), 0);
		gtk_combo_box_set_active(GTK_COMBO_BOX(comboX), X_selected_source);
		gtk_widget_set_sensitive(comboX, TRUE);
		if ((!is_fwhm) && registration_selected_source < r_X_POSITION) {
			registration_selected_source = r_QUALITY;
		}
		if ((is_fwhm) && registration_selected_source > r_Y_POSITION) {
			registration_selected_source = r_FWHM;
		}
	}
	gtk_cell_layout_clear(GTK_CELL_LAYOUT(combo));
	GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
	gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, TRUE);
	gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", 0, NULL);
	gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(combo), renderer, set_sensitive, NULL, NULL);

	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), (use_photometry) ? photometry_selected_source : registration_selected_source);

	gtk_cell_layout_clear(GTK_CELL_LAYOUT(comboX));
	GtkCellRenderer *rendererX = gtk_cell_renderer_text_new();
	gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(comboX), rendererX, TRUE);
	gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(comboX), rendererX, "text", 0, NULL);
	gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(comboX), rendererX, set_sensitive, NULL, NULL);

	g_signal_handlers_unblock_by_func(combo, on_plotCombo_changed, NULL);
	g_signal_handlers_unblock_by_func(comboX, on_plotComboX_changed, NULL);
}

void on_plotSourceCombo_changed(GtkComboBox *box, gpointer user_data) {
	validate_combos();
	requires_seqlist_update = TRUE;
	reset_plot_zoom();
	drawPlot();
}

void reset_plot() {
	free_plot_data();
	reset_plot_zoom();
	int layer;
	if (sourceCombo) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(sourceCombo), 0);
		gtk_combo_box_set_active(GTK_COMBO_BOX(combo), registration_selected_source); //remove?
		gtk_combo_box_set_active(GTK_COMBO_BOX(comboX), r_FRAME);
		gtk_widget_set_sensitive(comboX, TRUE);
		gtk_widget_set_sensitive(sourceCombo, FALSE);
		gtk_widget_set_visible(varCurve, FALSE);
		gtk_widget_set_sensitive(buttonSaveCSV, FALSE);
		gtk_widget_set_sensitive(buttonSavePrt, FALSE);
		gtk_widget_set_visible(julianw, FALSE);
		gtk_widget_set_sensitive(buttonClearLatest, FALSE);
		gtk_widget_set_sensitive(buttonClearAll, FALSE);
		layer = get_registration_layer(&com.seq);
		update_seqlist(layer);
	}
}

static int comparex(const void *a, const void *b) {
	struct kpair datax_a = * ((struct kpair *) a);
	struct kpair datax_b = * ((struct kpair *) b);

	if (datax_a.x > datax_b.x) return 1;
	if (datax_a.x < datax_b.x) return -1;
	return 0;
}

static int comparey(const void *a, const void *b) {
	struct kpair datax_a = * ((struct kpair *) a);
	struct kpair datax_b = * ((struct kpair *) b);

	if (datax_a.y == 0.) return 1;	// push zeros at the back
	if (datax_b.y == 0.) return -1;	// push zeros at the back
	if (datax_a.y > datax_b.y) return 1;
	if (datax_a.y < datax_b.y) return -1;
	return 0;
}

static int comparey_desc(const void *a, const void *b) {
	struct kpair datax_a = * ((struct kpair *) b); // inverted here
	struct kpair datax_b = * ((struct kpair *) a);

	if (datax_a.y == 0.) return -1;	// push zeros at the back
	if (datax_b.y == 0.) return 1;	// push zeros at the back
	if (datax_a.y > datax_b.y) return 1;
	if (datax_a.y < datax_b.y) return -1;
	return 0;
}

void drawPlot() {
	int ref_image;
	sequence *seq;

	validate_combos();

	seq = &com.seq;
	free_plot_data();

	if (seq->reference_image == -1)
		ref_image = 0;
	else ref_image = seq->reference_image;

	gboolean arcsec_is_ok = (gfit.focal_length > 0.0 && gfit.pixel_size_x > 0.f
		&& gfit.pixel_size_y > 0.f && gfit.binning_x > 0
		&& gfit.binning_y > 0);
	int current_selected_source = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	if (!arcsec_is_ok) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(arcsec), FALSE);
		is_arcsec = FALSE;
	}

	ref.x = -DBL_MAX;
	ref.y = -DBL_MAX;
	curr.x = -DBL_MAX;
	curr.y = -DBL_MAX;

	if (use_photometry) {
		// photometry data display
		pldata *plot;
		gtk_widget_set_visible(arcsec, current_selected_source == FWHM && arcsec_is_ok);
		update_ylabel();

		plot = alloc_plot_data(seq->number);
		plot_data = plot;
		// datamin/max must be init here to be common to all the stars
		pdd.datamin = (point){ DBL_MAX, DBL_MAX};
		pdd.datamax = (point){ -DBL_MAX, -DBL_MAX};
		for (int i = 0; i < MAX_SEQPSF && seq->photometry[i]; i++) {
			if (i > 0) {
				plot->next = alloc_plot_data(seq->number);
				plot = plot->next;
			}

			build_photometry_dataset(seq, i, seq->number, ref_image, plot);
			qsort(plot->data, plot->nb, sizeof(struct kpair), comparex);
		}
		if (requires_seqlist_update) { // update seq list if combo or arcsec changed
			update_seqlist(reglayer);
			requires_seqlist_update = FALSE;
		}
	} else {
		// registration data display
		if (!seq->regparam || reglayer < 0 || !seq->regparam[reglayer])
			return;

		is_fwhm = (seq->regparam[reglayer][ref_image].fwhm > 0.0f);
		gtk_widget_set_visible(arcsec, (current_selected_source == r_FWHM || current_selected_source == r_WFWHM || X_selected_source == r_FWHM || X_selected_source == r_WFWHM) && arcsec_is_ok);

		update_ylabel();
		/* building data array */
		plot_data = alloc_plot_data(seq->number);
		build_registration_dataset(seq, reglayer, ref_image, plot_data);
		if (requires_seqlist_update) { // update seq list if combo or arcsec changed
			update_seqlist(reglayer);
			requires_seqlist_update = FALSE;
		}
	}
	gtk_widget_set_sensitive(julianw, julian0);
	gtk_widget_queue_draw(drawingPlot);
}

static void set_filter(GtkFileChooser *dialog, const gchar *format) {
	GtkFileFilter *f = gtk_file_filter_new();
	gchar *name = g_strdup_printf(_("Output files (*%s)"), format);
	gchar *pattern = g_strdup_printf("*%s", format);
	gtk_file_filter_set_name(f, name);
	gtk_file_filter_add_pattern(f, pattern);
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);

	g_free(name);
	g_free(pattern);
}

static void save_dialog(const gchar *format, int (export_function)(pldata *, sequence *, gchar *)) {
	GtkWindow *control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
	SirilWidget *widgetdialog = siril_file_chooser_save(control_window, GTK_FILE_CHOOSER_ACTION_SAVE);
	GtkFileChooser *dialog = GTK_FILE_CHOOSER(widgetdialog);

	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	gtk_file_chooser_set_current_name(dialog, format);
	set_filter(dialog, format);

	gint res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = gtk_file_chooser_get_filename(dialog);
		export_function(plot_data, &com.seq, file);

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
}

void on_ButtonSaveCSV_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	save_dialog(".csv", exportCSV);
	set_cursor_waiting(FALSE);
}

void on_varCurvePhotometry_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	save_dialog(".dat", light_curve);
	set_cursor_waiting(FALSE);
}

void on_clearLatestPhotometry_clicked(GtkButton *button, gpointer user_data) {
	int i;
	for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++)
		;
	if (i != 0) {
		i--;
		free_photometry_set(&com.seq, i);
	}
	if (i == 0) {
		reset_plot();
		clear_stars_list(TRUE);
	}
	drawPlot();
}

void clear_all_photometry_and_plot() {
	for (int i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++) {
		free_photometry_set(&com.seq, i);
	}
	reset_plot();
	clear_stars_list(TRUE);
	drawPlot();
}

void on_clearAllPhotometry_clicked(GtkButton *button, gpointer user_data) {
	clear_all_photometry_and_plot();
}

void drawing_the_graph(GtkWidget *widget, cairo_t *cr, gboolean for_saving) {
	guint width, height;
	struct kplotcfg cfgplot;
	struct kdatacfg cfgdata;
	struct kdata *d1, *ref_d, *mean_d, *curr_d;

	if (!plot_data || !widget)
		return;
	pldata *plot = plot_data;
	d1 = ref_d = mean_d = NULL;

	double color = (com.pref.gui.combo_theme == 0) ? 0.0 : 1.0;

	kplotcfg_defaults(&cfgplot);
	kdatacfg_defaults(&cfgdata);
	set_colors(&cfgplot);
	cfgplot.ticlabel = TICLABEL_LEFT | TICLABEL_BOTTOM;
	cfgplot.border = BORDER_ALL;
	cfgplot.borderline.clr.type = KPLOTCTYPE_RGBA;
	cfgplot.borderline.clr.rgba[0] = 0.5;
	cfgplot.borderline.clr.rgba[1] = 0.5;
	cfgplot.borderline.clr.rgba[2] = 0.5;
	cfgplot.borderline.clr.rgba[3] = 1.0;
	cfgplot.xaxislabel = xlabel == NULL ? _("Frames") : xlabel;
	cfgplot.xtics = 5;
	cfgplot.yaxislabel = ylabel;
	cfgplot.yaxislabelrot = M_PI_2 * 3.0;
	cfgplot.xticlabelpad = cfgplot.yticlabelpad = 10.0;
	cfgplot.xticlabelfmt = formatX;
	cfgplot.yticlabelfmt = formatY;
	cfgdata.point.radius = 10;
	// binding the extrema to the sliders
	cfgplot.extrema = 0x0F;
	cfgplot.extrema_xmin = pdd.pdatamin.x;
	cfgplot.extrema_xmax = pdd.pdatamax.x;
	cfgplot.extrema_ymin = pdd.pdatamin.y;
	cfgplot.extrema_ymax = pdd.pdatamax.y;

	struct kplot *p = kplot_alloc(&cfgplot);

	// data plots
	int nb_graphs = 0;

	while (plot) {
		d1 = kdata_array_alloc(plot->data, plot->nb);
		if (X_selected_source == r_FRAME) {
			kplot_attach_data(p, d1,
					plot_data->nb <= 100 ? KPLOT_LINESPOINTS : KPLOT_LINES,
					NULL);
		} else {
			kplot_attach_data(p, d1, KPLOT_POINTS, NULL);
		}
		plot = plot->next;
		nb_graphs++;
	}

	/* mean and min/max */
	double mean = kdata_ymean(d1);
	int min_data = kdata_xmin(d1, NULL);
	int max_data = kdata_xmax(d1, NULL);

	if (nb_graphs == 1 && plot_data->nb > 0) {
		if (!use_photometry && (registration_selected_source == r_FWHM || registration_selected_source == r_WFWHM ||
					registration_selected_source == r_ROUNDNESS || registration_selected_source == r_QUALITY ||
					registration_selected_source == r_BACKGROUND || registration_selected_source == r_NBSTARS)) {
			if (X_selected_source == r_FRAME) {
				struct kpair *sorted_data;
				sorted_data = calloc(plot_data->nb, sizeof(struct kpair));
				for (int i = 0; i < plot_data->nb; i++) {
					sorted_data[i].x = plot_data->data[i].x;
					sorted_data[i].y = plot_data->data[i].y;
				}
				qsort(sorted_data, plot_data->nb, sizeof(struct kpair),
						(registration_selected_source == r_ROUNDNESS || registration_selected_source == r_QUALITY ||
						 registration_selected_source == r_NBSTARS) ? comparey_desc : comparey);
				double imin = pdd.pdatamin.x;
				double imax = pdd.pdatamax.x;
				double pace = (imax - imin) / ((double)plot_data->nb - 1.);
				for (int i = 0; i < plot_data->nb; i++) {
					sorted_data[i].x = imin + (double)i * pace;
				}
				d1 = kdata_array_alloc(sorted_data, plot_data->nb);
				kplot_attach_data(p, d1, KPLOT_LINES, NULL);
				free(sorted_data);
			}
		} else if (use_photometry){
			int nb_data = max_data - min_data + 1;
			struct kpair *avg = calloc(nb_data, sizeof(struct kpair));
			for (int i = 0, j = min_data; i < nb_data; i++, j++) {
				avg[i].x = plot_data->data[j].x;
				avg[i].y = mean;
			}
			mean_d = kdata_array_alloc(avg, nb_data);
			kplot_attach_data(p, mean_d, KPLOT_LINES, NULL);	// mean plot
			free(avg);
		}

		if (ref.x > -DBL_MAX && ref.y > -DBL_MAX) {
			ref_d = kdata_array_alloc(&ref, 1);
			kplot_attach_data(p, ref_d, KPLOT_POINTS, &cfgdata);	// ref image dot
		}
		if (curr.x > -DBL_MAX && curr.y > -DBL_MAX) {
			curr_d = kdata_array_alloc(&curr, 1);
			kplot_attach_data(p, curr_d, KPLOT_MARKS, &cfgdata);	// ref image dot
		}
	}

	width =  gtk_widget_get_allocated_width(widget);
	height = gtk_widget_get_allocated_height(widget);
	pdd.surf_w = (double)width;
	pdd.surf_h = (double)height;
	cairo_surface_t *surface = NULL;

	if (for_saving) {
		surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
		cr = cairo_create(surface);
	}

	cairo_set_source_rgb(cr, color, color, color);
	cairo_rectangle(cr, 0.0, 0.0, width, height);
	cairo_fill(cr);
	kplot_draw(p, width, height, cr);

	if (for_saving) {
		gchar *timestamp, *filename;
		timestamp = build_timestamp_filename();
		filename = g_strdup_printf("%s.png", timestamp);
		g_free(timestamp);

		cairo_surface_write_to_png(surface, filename);
		cairo_surface_destroy(surface);
		siril_log_message(_("Saved plot to image %s\n"), filename);
		g_free(filename);
	} else {
		// caching more data
		pdd.range = (point){ get_dimx(),  get_dimy()};
		pdd.scale = (point){ (pdd.pdatamax.x - pdd.pdatamin.x) / pdd.range.x, (pdd.pdatamax.y - pdd.pdatamin.y) / pdd.range.y};
		pdd.offset = (point){ get_offsx(),  get_offsy()};
		// dealing with selection here after plot specifics have been updated. Otherwise change of scale is flawed (when arsec/julian state are changed)
		pdd.selected = calloc(com.seq.number, sizeof(gboolean));
		if (selection_is_active()) {
			double xmin, ymin, xmax, ymax;
			int n = 0;
			convert_surface_to_plot_coords(pdd.selection.x, pdd.selection.y, &xmin, &ymax);
			convert_surface_to_plot_coords(pdd.selection.x + pdd.selection.w, pdd.selection.y + pdd.selection.h, &xmax, &ymin);
			int i, j;
			for (i = 0, j = 0; i < com.seq.number; i++) {
				if (!com.seq.imgparam[i].incl) continue;
				if (plot_data->data[j].x >= xmin && plot_data->data[j].x <= xmax && plot_data->data[j].y >= ymin && plot_data->data[j].y <= ymax) {
					n++;
					pdd.selected[i] = TRUE;
				}
				j++;
			}
			pdd.nbselected = n;
		} else {
			pdd.nbselected = 0;
		}
		//drawing the sliders and markers
		plot_draw_all_sliders(cr);
		plot_draw_all_sliders_fill(cr);
		plot_draw_all_markers(cr);
		plot_draw_selection(cr);
	}
	free_colors(&cfgplot);
	kplot_free(p);
	kdata_destroy(d1);
	kdata_destroy(ref_d);
	if (mean_d)
		kdata_destroy(mean_d);
}

gboolean on_DrawingPlot_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	drawing_the_graph(widget, cr, FALSE);
	return FALSE;
}

void on_ButtonSavePlot_clicked(GtkButton *button, gpointer user_data) {
	drawing_the_graph((GtkWidget *)user_data, NULL, TRUE);
}

void on_plotCombo_changed(GtkComboBox *box, gpointer user_data) {
	if (use_photometry) {
		photometry_selected_source = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	} else {
		registration_selected_source = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	}
	requires_seqlist_update = TRUE;
	pdd.selection = (rectangled){0., 0., 0., 0.};
	update_slider(SLIDER_Y, 0., 1.);
}

void on_plotComboX_changed(GtkComboBox *box, gpointer user_data) {
	X_selected_source = gtk_combo_box_get_active(GTK_COMBO_BOX(comboX));
	pdd.selection = (rectangled){0., 0., 0., 0.};
	update_slider(SLIDER_X, 0., 1.);
}

void on_arcsecPhotometry_toggled(GtkToggleButton *button, gpointer user_data) {
	requires_seqlist_update = TRUE;
	is_arcsec = gtk_toggle_button_get_active(button);
	drawPlot();
}

void on_JulianPhotometry_toggled(GtkToggleButton *button, gpointer user_data) {
	force_Julian = gtk_toggle_button_get_active(button);
	drawPlot();
}

static void update_ylabel() {
	int current_selected_source = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	gtk_widget_set_sensitive(buttonSavePrt, TRUE);
	if (use_photometry) {
		gtk_widget_set_sensitive(buttonSaveCSV, TRUE);
		gtk_widget_set_sensitive(varCurve, current_selected_source == MAGNITUDE);
		switch (current_selected_source) {
			case ROUNDNESS:
				ylabel = _("Star roundness (1 is round)");
				break;
			case FWHM:
				ylabel = is_arcsec ? _("FWHM ('')") : _("FWHM (px)");
				break;
			case AMPLITUDE:
				ylabel = _("Amplitude");
				break;
			case MAGNITUDE:
				if (com.magOffset > 0.0 || com.seq.reference_star >= 0)
					ylabel = _("Star magnitude (absolute)");
				else
					ylabel = _("Star magnitude (relative, use setmag)");
				break;
			case BACKGROUND:
				ylabel = _("Background value");
				break;
			case X_POSITION:
				ylabel = _("Star position on X axis");
				break;
			case Y_POSITION:
				ylabel = _("Star position on Y axis");
				break;
			case SNR:
				ylabel = _("SNR (dB)");
				break;
			default:
				break;
		}
	} else {
		gtk_widget_set_sensitive(buttonSaveCSV, TRUE);
		switch (current_selected_source) {
			case r_ROUNDNESS:
				ylabel = _("Star roundness (1 is round)");
				break;
			case r_FWHM:
				ylabel = (is_arcsec) ? _("FWHM ('')") : _("FWHM (px)");
				break;
			case r_WFWHM:
				ylabel = (is_arcsec) ? _("wFWHM ('')") : _("wFWHM (px)");
				break;
			case r_X_POSITION:
				ylabel = _("Shift on X axis");
				break;
			case r_Y_POSITION:
				ylabel = _("Shift on Y axis");
				break;
			case r_QUALITY:
				ylabel = _("Quality");
				break;
			case r_BACKGROUND:
				ylabel = _("Background value");
				break;
			case r_NBSTARS:
				ylabel = _("Number of stars");
				break;
			default:
				break;
		}
		switch (X_selected_source) {
			case r_ROUNDNESS:
				xlabel = _("Star roundness (1 is round)");
				break;
			case r_FWHM:
				xlabel = (is_arcsec) ? _("FWHM ('')") : _("FWHM (px)");
				break;
			case r_WFWHM:
				xlabel = (is_arcsec) ? _("wFWHM ('')") : _("wFWHM (px)");
				break;
			case r_X_POSITION:
				xlabel = _("Shift on X axis");
				break;
			case r_Y_POSITION:
				xlabel = _("Shift on Y axis");
				break;
			case r_QUALITY:
				xlabel = _("Quality");
				break;
			case r_BACKGROUND:
				xlabel = _("Background value");
				break;
			case r_NBSTARS:
				xlabel = _("Number of stars");
				break;
			case r_FRAME:
				xlabel = _("Frames");
				break;
			default:
				break;
		}
	}
}

/* initialize the colors of each star for photometry, outside
 * on_DrawingPlot_draw because it may not be called if the panel is hidden, but
 * the circles in the main image display still use these colors.
 */
void init_plot_colors() {
	struct kplotcfg cfgplot;
	set_colors(&cfgplot);
	/* copy graph colours for star highlight */
	for (int i = 0; i < cfgplot.clrsz; i++) {
		com.seq.photometry_colors[i][0] = cfgplot.clrs[i].rgba[0];
		com.seq.photometry_colors[i][1] = cfgplot.clrs[i].rgba[1];
		com.seq.photometry_colors[i][2] = cfgplot.clrs[i].rgba[2];
	}
}

void notify_new_photometry() {
	control_window_switch_to_tab(PLOT);
	init_plot_colors();
	gtk_widget_set_sensitive(sourceCombo, TRUE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(sourceCombo), 1);
	gtk_widget_set_sensitive(buttonClearLatest, TRUE);
	gtk_widget_set_sensitive(buttonClearAll, TRUE);
	gtk_widget_set_sensitive(comboX, FALSE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(comboX), r_FRAME);
}

/* copied from subprojects/kplot/draw.c, in which the max is 7, not MAX_SEQPSF.
 * We assign the same color to all others */
static void set_colors(struct kplotcfg *cfg) {
	int i;
	cfg->clrsz = MAX_SEQPSF;
	cfg->clrs = calloc(cfg->clrsz, sizeof(struct kplotccfg));
	for (i = 0; i < cfg->clrsz; i++) {
		cfg->clrs[i].type = KPLOTCTYPE_RGBA;
		cfg->clrs[i].rgba[3] = 1.0;
		if (i > 6) {
			cfg->clrs[i].rgba[0] = 0x00;
			cfg->clrs[i].rgba[1] = 0xaa;
			cfg->clrs[i].rgba[2] = 0xbb;
		}
	}
	cfg->clrs[0].rgba[0] = 0x94 / 255.0;
	cfg->clrs[0].rgba[1] = 0x04 / 255.0;
	cfg->clrs[0].rgba[2] = 0xd3 / 255.0;
	cfg->clrs[1].rgba[0] = 0x00 / 255.0;
	cfg->clrs[1].rgba[1] = 0x9e / 255.0;
	cfg->clrs[1].rgba[2] = 0x73 / 255.0;
	cfg->clrs[2].rgba[0] = 0x56 / 255.0;
	cfg->clrs[2].rgba[1] = 0xb4 / 255.0;
	cfg->clrs[2].rgba[2] = 0xe9 / 255.0;
	cfg->clrs[3].rgba[0] = 0xe6 / 255.0;
	cfg->clrs[3].rgba[1] = 0x9f / 255.0;
	cfg->clrs[3].rgba[2] = 0x00 / 255.0;
	cfg->clrs[4].rgba[0] = 0xf0 / 255.0;
	cfg->clrs[4].rgba[1] = 0xe4 / 255.0;
	cfg->clrs[4].rgba[2] = 0x42 / 255.0;
	cfg->clrs[5].rgba[0] = 0x00 / 255.0;
	cfg->clrs[5].rgba[1] = 0x72 / 255.0;
	cfg->clrs[5].rgba[2] = 0xb2 / 255.0;
	cfg->clrs[6].rgba[0] = 0xe5 / 255.0;
	cfg->clrs[6].rgba[1] = 0x1e / 255.0;
	cfg->clrs[6].rgba[2] = 0x10 / 255.0;
}

static void free_colors(struct kplotcfg *cfg) {
	free(cfg->clrs);
}

gboolean on_DrawingPlot_motion_notify_event(GtkWidget *widget,
		GdkEventMotion *event, gpointer user_data) {

	if (!plot_data) return FALSE;
	if (!com.seq.imgparam) return FALSE;

	gtk_widget_set_has_tooltip(widget, FALSE);

	double x = (double)event->x;
	double y = (double)event->y;
	if (pdd.action == SELACTION_SELECTING) {
		double x1, x2, y1, y2;
		if (x <= pdd.selection.x) {
			x1 = x;
			x2 = pdd.selection.x + pdd.selection.w;
		} else {
			x1 = pdd.selection.x;
			x2 = x;
		}
		if (y <= pdd.selection.y) {
			y1 = y;
			y2 = pdd.selection.y + pdd.selection.h;
		} else {
			y1 = pdd.selection.y;
			y2 = y;
		}
		x1 = max(SIDE_MARGIN, x1);
		x2 = min(pdd.surf_w - PLOT_SLIDER_THICKNESS, x2);
		y1 = max(SIDE_MARGIN, y1);
		y2 = min(pdd.surf_h - PLOT_SLIDER_THICKNESS, y2);
		//siril_debug_print("%.0f %.0f\n", x1, x2);
		pdd.selection = (rectangled){x1, y1, x2 - x1, y2 - y1};
		drawPlot();
		return TRUE;
	}
	if (pdd.action == SELACTION_RESIZING) {
		double x1, x2, y1, y2;
		x1 = pdd.selection.x;
		y1 = pdd.selection.y;
		x2 = pdd.selection.x + pdd.selection.w;
		y2 = pdd.selection.y + pdd.selection.h;
		switch (pdd.border_grabbed) {
			case SELBORDER_LEFT:
				x1 = min(x , pdd.selection.x + pdd.selection.w);
				x2 = max(x , pdd.selection.x + pdd.selection.w);
				break;
			case SELBORDER_RIGHT:
				x1 = min(x , pdd.selection.x);
				x2 = max(x , pdd.selection.x);
				break;
			case SELBORDER_TOP:
				y1 = min(y , pdd.selection.y + pdd.selection.h);
				y2 = max(y , pdd.selection.y + pdd.selection.h);
				break;
			case SELBORDER_BOTTOM:
				y1 = min(y , pdd.selection.y);
				y2 = max(y , pdd.selection.y);
				break;
			default:
				break;
		}
		x1 = max(SIDE_MARGIN, x1);
		x2 = min(pdd.surf_w - PLOT_SLIDER_THICKNESS, x2);
		y1 = max(SIDE_MARGIN, y1);
		y2 = min(pdd.surf_h - PLOT_SLIDER_THICKNESS, y2);
		pdd.selection = (rectangled){x1, y1, x2 - x1, y2 - y1};
		drawPlot();
		return TRUE;
	}
	if (pdd.action == SELACTION_MOVING) {
		double dx, dy;
		dx = x - pdd.start.x;
		dy = y - pdd.start.y;
		if (pdd.selection.x + dx < SIDE_MARGIN) {
			dx = -pdd.selection.x + SIDE_MARGIN;
			x = dx + pdd.start.x;
		}
		if (pdd.selection.y + dy < SIDE_MARGIN) {
			dy = -pdd.selection.y + SIDE_MARGIN;
			y = dy + pdd.start.y;
		}
		if (pdd.selection.x + pdd.selection.w + dx > pdd.surf_w - PLOT_SLIDER_THICKNESS) {
			dx = - (pdd.selection.x + pdd.selection.w) - PLOT_SLIDER_THICKNESS + pdd.surf_w;
			x = dx + pdd.start.x;
		}
		if (pdd.selection.y + pdd.selection.h + dy > pdd.surf_h - PLOT_SLIDER_THICKNESS) {
			dy = - (pdd.selection.y + pdd.selection.h) - PLOT_SLIDER_THICKNESS + pdd.surf_h;
			y = dy + pdd.start.y;
		}
		pdd.start.x = x;
		pdd.start.y = y;
		pdd.selection = (rectangled){pdd.selection.x + dx, pdd.selection.y + dy, pdd.selection.w, pdd.selection.h};
		drawPlot();
		return TRUE;
	}
	for (int i = SLIDER_X; i <= SLIDER_Y; i++) {
		for (int j = 0; j < 2; j++) {
			if (pdd.marker_grabbed == 2 * i + j) {
				double *valrange = (i == 0) ? pdd.xrange : pdd.yrange;
				find_range_from_pos(x, y, i, j, valrange);
				update_slider(i, valrange[0], valrange[1]);
				return TRUE;
			} else if (pdd.slider_grabbed == i) {
				double dv, v1, v2;
				if (i == SLIDER_X) {
					dv = x - pdd.start.x;
					v1 = pdd.xrange[0] * pdd.range.x + pdd.offset.x;
					v2 = pdd.xrange[1] * pdd.range.x + pdd.offset.x;
					if (v1 + dv <= pdd.offset.x) {
						dv = pdd.offset.x - v1;
						x = dv + pdd.start.x;
					}
					if (v2 + dv >= pdd.offset.x + pdd.range.x) {
						dv = pdd.offset.x + pdd.range.x - v2;
						x = dv + pdd.start.x;
					}
				} else {
					dv = pdd.start.y - y;
					v1 = (1. - pdd.yrange[0]) * pdd.range.y + pdd.offset.y; // y axis is reversed
					v2 = (1. - pdd.yrange[1]) * pdd.range.y + pdd.offset.y; 
					if (v2 - dv <= pdd.offset.y) {
						dv = v2 - pdd.offset.y;
						y = pdd.start.y - dv;
					}
					if (v1 - dv >= pdd.offset.y + pdd.range.y) {
						dv = -pdd.offset.y - pdd.range.y + v1;
						y = pdd.start.y - dv;
					}
				}
				pdd.start.x = x;
				pdd.start.y = y;
				if (i == SLIDER_X) {
					v1 += dv;
					v2 += dv;
					pdd.xrange[0] = find_rangeval_from_pos(v1, i);
					pdd.xrange[1] = find_rangeval_from_pos(v2, i);
					update_slider(i, pdd.xrange[0], pdd.xrange[1]);
				} else {
					v1 -= dv;
					v2 -= dv;
					pdd.yrange[0] = find_rangeval_from_pos(v1, i);
					pdd.yrange[1] = find_rangeval_from_pos(v2, i);
					update_slider(i, pdd.yrange[0], pdd.yrange[1]);
				}
				return TRUE;
			} else if (is_over_marker(x, y, 2 * i + j)) {
				set_cursor("grab");
				return TRUE;
			}
		}
	}
	enum border_type border = is_over_selection_border(x, y);
	if (border > SELBORDER_NONE) {
		if (border <= SELBORDER_BOTTOM)
			set_cursor("n-resize");
		else
			set_cursor("w-resize");
		return TRUE;
	} else if (is_inside_selection(x, y)) {
		set_cursor("all-scroll");
	} else {
		set_cursor("tcross");
	}
	if (is_inside_grid(x, y)) {
		double index, xpos, ypos;
		gboolean getvals = get_index_of_frame(x, y, FALSE, &index, &xpos, &ypos);
		gchar *tooltip_text;
		if (getvals) {
			if (index > 0) {
				tooltip_text = g_strdup_printf("X pos: %0.3f\nY pos: %0.3f\nFrame#%d", xpos, ypos,(int)index);
			} else {
				tooltip_text = g_strdup_printf("X pos: %0.3f\nY pos: %0.3f", xpos, ypos);
			}
			gtk_widget_set_tooltip_text(widget, tooltip_text);
			g_free(tooltip_text);
			return TRUE;
		}
	}
	set_cursor("tcross");
	return TRUE;
}

gboolean on_DrawingPlot_enter_notify_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (plot_data && pdd.marker_grabbed == MARKER_NONE && pdd.border_grabbed == SELBORDER_NONE && pdd.slider_grabbed == SLIDER_NONE && pdd.action == SELACTION_NONE)
		set_cursor("tcross");
	return TRUE;
}

gboolean on_DrawingPlot_leave_notify_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (get_thread_run()) {
		set_cursor_waiting(TRUE);
	} else {
		/* trick to get default cursor */
		if (pdd.marker_grabbed == MARKER_NONE && pdd.border_grabbed == SELBORDER_NONE && pdd.slider_grabbed == SLIDER_NONE && pdd.action == SELACTION_NONE)
			set_cursor_waiting(FALSE);
	}
	return TRUE;
}

static void do_popup_singleframemenu(GtkWidget *my_widget, GdkEventButton *event) {
	if (!event) {
		return;
	}

	double index, xpos, ypos;
	gboolean getvals = get_index_of_frame((double)event->x,(double)event->y, TRUE, &index, &xpos, &ypos);
	if (!getvals) return;
	if (index < 0) return;

	if (!menu) {
		menu = GTK_MENU(lookup_widget("menu_plot"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
		menu_item1 = GTK_MENU_ITEM(lookup_widget("menu_plot_item1"));
		menu_item2 = GTK_MENU_ITEM(lookup_widget("menu_plot_item2"));
		menu_item3 = GTK_MENU_ITEM(lookup_widget("menu_plot_item3"));
	}
	gchar *str = g_strdup_printf(_("Exclude Frame %d"), (int)index);
	gtk_menu_item_set_label(menu_item1, str);
	gchar *str2 = g_strdup_printf(_("Show Frame %d"), (int)index);
	gtk_menu_item_set_label(menu_item2, str2);
	gtk_menu_item_set_label(menu_item3, "");
	gtk_widget_set_sensitive(GTK_WIDGET(menu_item3), FALSE);
	g_free(str);
	g_free(str2);

#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
	int button = event->button;
	int event_time = event->time;

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button,
			event_time);
#endif
}

static void do_popup_selectionmenu(GtkWidget *my_widget, GdkEventButton *event) {
	if (!event) {
		return;
	}

	if (!menu) {
		menu = GTK_MENU(lookup_widget("menu_plot"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
		menu_item1 = GTK_MENU_ITEM(lookup_widget("menu_plot_item1"));
		menu_item2 = GTK_MENU_ITEM(lookup_widget("menu_plot_item2"));
		menu_item3 = GTK_MENU_ITEM(lookup_widget("menu_plot_item3"));
	}
	gtk_menu_item_set_label(menu_item1, _("Zoom to selection"));
	gtk_menu_item_set_label(menu_item2, _("Only keep points within selection"));
	gtk_menu_item_set_label(menu_item3, _("Exclude selected points"));
	gtk_widget_set_sensitive(GTK_WIDGET(menu_item3), TRUE);

#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
	int button = event->button;
	int event_time = event->time;

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button,
			event_time);
#endif
}

gboolean on_DrawingPlot_button_press_event(GtkWidget *widget,
	GdkEventButton *event, gpointer user_data) {
	if (!plot_data) return FALSE;
	if (!com.seq.imgparam) return FALSE;
	double x = (double)event->x;
	double y = (double)event->y;
	// if a zone is selected, right-click pops out menu to zoom/keep/exclude in bulk
	if (event->button == GDK_BUTTON_SECONDARY) {
		if (selection_is_active()) {
			do_popup_selectionmenu(widget, event);
			return TRUE;
		}
		// open or exclude image (if close enough to a data point)
		if (is_inside_grid(x, y) && event->button == GDK_BUTTON_SECONDARY) {
			do_popup_singleframemenu(widget, event);
			return TRUE;
		}
	}
	
	if (is_inside_selectable_zone(x, y) && event->button == GDK_BUTTON_PRIMARY) {
		enum border_type border = is_over_selection_border(x, y);
		if (event->type == GDK_DOUBLE_BUTTON_PRESS) {  // double-click resets zoom
			reset_plot_zoom();
			return TRUE;
		} else if (is_inside_selection(x, y)) { // start moving selection
			pdd.action = SELACTION_MOVING;
			pdd.start = (point){x, y};
			return TRUE;
		} else if (border > SELBORDER_NONE) { // start resizing selection
			pdd.action = SELACTION_RESIZING;
			pdd.border_grabbed = border;
			return TRUE;
		} else { // start drawing selection
			pdd.action = SELACTION_SELECTING;
			pdd.selection = (rectangled){x, y, 0., 0.};
			return TRUE;
		}
	}

	for (int i = SLIDER_X; i <= SLIDER_Y; i++) {
		if (is_inside_slider(x, y, i)) {
			if (event->button == GDK_BUTTON_PRIMARY) {
				// double - click on slider resets both markers
				if (event->type == GDK_DOUBLE_BUTTON_PRESS) {
					update_slider(i, 0., 1.);
					return TRUE;
				}
				for (int j = 0; j < 2; j++) {
					if (is_over_marker(x, y, 2 * i + j) && pdd.marker_grabbed == MARKER_NONE) {
						pdd.marker_grabbed = 2 * i + j;
						set_cursor("grabbing");
						return TRUE;
					}
				}
				//otherwise, it's just clicked once - we find the closest marker
				// In case of double-click, that's called once first... we decided to live with that
				double *valrange = (i == 0) ? &pdd.xrange[0] : &pdd.yrange[0];
				int j = get_closest_marker(x, y, i, valrange);
				find_range_from_pos(x, y, i, j, valrange);
				update_slider(i, valrange[0], valrange[1]);
				return TRUE;
			} else if (event->button == GDK_BUTTON_SECONDARY) {
				(i == SLIDER_X) ? set_cursor("w-resize") : set_cursor("n-resize");
				pdd.slider_grabbed = i;
				pdd.start = (point){x, y};
				return TRUE;
			}
		}
	}
	return TRUE;
}

gboolean on_DrawingPlot_button_release_event(GtkWidget *widget,
	GdkEventButton *event, gpointer user_data) {
	if (plot_data) {
		set_cursor("tcross");
	} else {
		set_cursor_waiting(FALSE);
	}
	pdd.marker_grabbed = MARKER_NONE;
	pdd.action = SELACTION_NONE;
	pdd.border_grabbed = SELBORDER_NONE;
	pdd.slider_grabbed = SLIDER_NONE;
	if (pdd.selection.w < 1. || pdd.selection.h < 1. )
		pdd.selection = (rectangled){0., 0., 0., 0.};
	drawPlot();
	return TRUE;
}

static signed long extract_int_from_label(const gchar *str) {
	gchar *p = (gchar *)str;
	while (*p) {
		if (g_ascii_isdigit(*p) || ((*p == '-' || *p == '+') && g_ascii_isdigit(*(p + 1)))) {
			// Found a number
			return g_ascii_strtoll(p, &p, 10); // return number
		} else {
			// Otherwise, move on to the next character.
			p++;
		}
	}
	return -1;
}

void on_menu_plot_item1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (!selection_is_active()) { // Exclude single frame
		const gchar *label = gtk_menu_item_get_label(menuitem);
		gint index;

		index = extract_int_from_label(label);
		if (index > 0) {
			index--;

			exclude_single_frame(index);
			update_seqlist(use_photometry ? 0 : reglayer);
		}
	} else { // Zoom to selection
		double xmin, xmax, ymin, ymax;
		convert_surface_to_plot_coords(pdd.selection.x, pdd.selection.y, &xmin, &ymax);
		convert_surface_to_plot_coords(pdd.selection.x + pdd.selection.w, pdd.selection.y + pdd.selection.h, &xmax, &ymin);
		pdd.xrange[0] = (xmin - pdd.datamin.x) / (pdd.datamax.x - pdd.datamin.x);
		pdd.xrange[1] = (xmax - pdd.datamin.x) / (pdd.datamax.x - pdd.datamin.x);
		pdd.yrange[0] = (ymin - pdd.datamin.y) / (pdd.datamax.y - pdd.datamin.y);
		pdd.yrange[1] = (ymax - pdd.datamin.y) / (pdd.datamax.y - pdd.datamin.y);
		pdd.selection = (rectangled){0., 0., 0., 0.};
		update_slider(SLIDER_X, pdd.xrange[0], pdd.xrange[1]);
		update_slider(SLIDER_Y, pdd.yrange[0], pdd.yrange[1]);
	}
}

void on_menu_plot_item2_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (!selection_is_active()) { // Open single frame
		const gchar *label = gtk_menu_item_get_label(menuitem);
		gint index;

		index = extract_int_from_label(label);
		if (index > 0) {
			siril_open_dialog("seqlist_dialog");
			update_seqlist(use_photometry ? 0 : reglayer);
			sequence_list_select_row_from_index(index - 1, TRUE);
		}
	} else {
		select_unselect_frames_from_list(pdd.selected, TRUE);
		reset_plot_zoom();
		drawPlot();
	}
}
void on_menu_plot_item3_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (!selection_is_active()) { // Open single frame
		return;
	} else {
		select_unselect_frames_from_list(pdd.selected, FALSE);
		reset_plot_zoom();
		drawPlot();
	}
}

