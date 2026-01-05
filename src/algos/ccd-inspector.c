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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "gui/image_display.h"
#include "gui/dialogs.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"

#include "ccd-inspector.h"

static void draw_polygon(float rx, float ry, float m1, float m2, float m3, float m4, float mcentre) {
	float r1, r2, r3, r4;
	pointf c = { rx / 2.f, ry / 2.f };
	float m = (m1 + m2 + m3 + m4) / 4.f;
	float diag = sqrtf(rx * rx + ry * ry) / 4.f;

	/* now we compute the four radius. */
	r1 = diag * (((m1 - m) / m) + 1);
	r2 = diag * (((m2 - m) / m) + 1);
	r3 = diag * (((m3 - m) / m) + 1);
	r4 = diag * (((m4 - m) / m) + 1);

	com.tilt = malloc(sizeof(sensor_tilt));

	com.tilt->pt[0].x = c.x + (r1 * sin(7.0 * M_PI / 4.0));
	com.tilt->pt[0].y = ry - (c.y + (r1 * cos(7.0 * M_PI / 4.0)));
	com.tilt->fwhm[0] = m1;

	com.tilt->pt[1].x = c.x + (r2 * sin(M_PI / 4));
	com.tilt->pt[1].y = ry - (c.y + (r2 * cos(M_PI / 4.0)));
	com.tilt->fwhm[1] = m2;

	com.tilt->pt[2].x = c.x + (r3 * sin(5.0 * M_PI / 4));
	com.tilt->pt[2].y = ry - (c.y + (r3 * cos(5.0 * M_PI / 4.0)));
	com.tilt->fwhm[2] = m3;

	com.tilt->pt[3].x = c.x + (r4 * sin(3.0 * M_PI / 4));
	com.tilt->pt[3].y = ry - (c.y + (r4 * cos(3.0 * M_PI / 4.0)));
	com.tilt->fwhm[3] = m4;

	com.tilt->fwhm_centre = mcentre;

	redraw(REDRAW_OVERLAY);
}

void clear_sensor_tilt() {
	free(com.tilt);
	com.tilt = NULL;
}

static int compute_tilt_values(fits *fit, int nbstars, psf_star **stars, float *m, float *m1, float *m2, float *m3, float *m4, float *mr1, float *mr2) {
	int ret = 1;
	int i1 = 0, i2 = 0, i3 = 0, i4 = 0, ir1 = 0, ir2 = 0;
	pointf center = {fit->rx / 2.f, fit->ry / 2.f };

	float r = sqrtf(center.x * center.x + center.y * center.y);
	float r1 = 0.25f * r;
	float r2 = 0.75f * r;

	float *f = malloc(nbstars * sizeof(float));

	float *f1 = calloc(nbstars, sizeof(float));
	float *f2 = calloc(nbstars, sizeof(float));
	float *f3 = calloc(nbstars, sizeof(float));
	float *f4 = calloc(nbstars, sizeof(float));

	float *fr1 = calloc(nbstars, sizeof(float));
	float *fr2 = calloc(nbstars, sizeof(float));

	for (int i = 0; i < nbstars; i++) {
		float x = (float) stars[i]->xpos;
		float y = (float) stars[i]->ypos;

		/* global */
		f[i] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		/* check for 4 areas */
		if ((x < center.x) && (y < center.y)) {
			f1[i1++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if ((x > center.x) && (y < center.y)) {
			f2[i2++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if ((x < center.x) && (y > center.y)) {
			f3[i3++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if ((x > center.x) && (y > center.y)) {
			f4[i4++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		}
		/* check for off-axis aberration */
		if (((x - center.x) * (x - center.x) + (y - center.y) * (y - center.y)) < (r1 * r1)) {
			fr1[ir1++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if (((x - center.x) * (x - center.x) + (y - center.y) * (y - center.y)) > (r2 * r2)) {
			fr2[ir2++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		}
	}

	if ((i1 != 0) && (i2 != 0) && (i3 != 0) && (i4 != 0) && (ir1 != 0) && (ir2 != 0)) {
		quicksort_f(f, nbstars);
		quicksort_f(f1, i1);
		quicksort_f(f2, i2);
		quicksort_f(f3, i3);
		quicksort_f(f4, i4);
		quicksort_f(fr1, ir1);
		quicksort_f(fr2, ir2);

		*m = siril_stats_trmean_from_sorted_data(0.25, f, 1, nbstars);
		*m1 = siril_stats_trmean_from_sorted_data(0.25, f1, 1, i1);
		*m2 = siril_stats_trmean_from_sorted_data(0.25, f2, 1, i2);
		*m3 = siril_stats_trmean_from_sorted_data(0.25, f3, 1, i3);
		*m4 = siril_stats_trmean_from_sorted_data(0.25, f4, 1, i4);

		*mr1 = siril_stats_trmean_from_sorted_data(0.25, fr1, 1, ir1);
		*mr2 = siril_stats_trmean_from_sorted_data(0.25, fr2, 1, ir2);

		ret = 0;
	}

	free(f);

	free(f1);
	free(f2);
	free(f3);
	free(f4);

	free(fr1);
	free(fr2);

	return ret;
}

int draw_sensor_tilt(fits *fit) {
	int nbstars = 0;

	float m = 0;
	float m1 = 0;
	float m2 = 0;
	float m3 = 0;
	float m4 = 0;
	float mr1 = 0;
	float mr2 = 0;

	delete_selected_area();

	image im = { .fit = fit, .from_seq = NULL, .index_in_seq = -1 };
	psf_star **stars = peaker(&im, select_vport(gui.cvport), &com.pref.starfinder_conf, &nbstars, NULL, FALSE, FALSE, MAX_STARS_FITTED, com.pref.starfinder_conf.profile, com.max_thread);

	if (!compute_tilt_values(fit, nbstars, stars, &m, &m1, &m2, &m3, &m4, &mr1, &mr2)) {
		float best = min(min(m1, m2), min(m3, m4));
		float worst = max(max(m1, m2), max(m3, m4));

		float ref = (m1 + m2 + m3 + m4) / 4.f;

		draw_polygon((float) fit->rx, (float) fit->ry, m1, m2, m3, m4, mr1);

		siril_log_message(_("Stars: %d, Truncated mean[FWHM]: %.2f, Sensor tilt[FWHM]: %.2f (%.0f%%), Off-axis aberration[FWHM]: %.2f\n"),
				nbstars, m, worst - best, roundf(((worst - best) / ref) * 100.f), mr2 - mr1);
	}

	free_fitted_stars(stars);
	return 0;
}

static int compute_tilt_to_image(image *im, struct tilt_data *t_args) {
	int nbstars = 0;
	int layer = im->fit->naxes[2] > 1 ? GLAYER : RLAYER;

	psf_star **stars = peaker(im, layer, &com.pref.starfinder_conf, &nbstars, NULL, FALSE, FALSE, MAX_STARS_FITTED, com.pref.starfinder_conf.profile, com.max_thread);

	float m = 0;
	float m1 = 0;
	float m2 = 0;
	float m3 = 0;
	float m4 = 0;
	float mr1 = 0;
	float mr2 = 0;

	if (!compute_tilt_values(im->fit, nbstars, stars, &m, &m1, &m2, &m3, &m4, &mr1, &mr2)) {
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			t_args->m += m;
			t_args->m1 += m1;
			t_args->m2 += m2;
			t_args->m3 += m3;
			t_args->m4 += m4;
			t_args->mr1 += mr1;
			t_args->mr2 += mr2;

			t_args->nbstars += nbstars;
		}

		free_fitted_stars(stars);
	}

	return 0;
}


/** Tilt on sequence **/

static int tilt_finalize_hook(struct generic_seq_args *args) {
	struct tilt_data *t_args = (struct tilt_data*) args->user;

	t_args->m /= args->seq->selnum;
	t_args->m1 /= args->seq->selnum;
	t_args->m2 /= args->seq->selnum;
	t_args->m3 /= args->seq->selnum;
	t_args->m4 /= args->seq->selnum;
	t_args->mr1 /= args->seq->selnum;
	t_args->mr2 /= args->seq->selnum;

	t_args->nbstars /= args->seq->selnum;

	float best = min(min(t_args->m1, t_args->m2), min(t_args->m3, t_args->m4));
	float worst = max(max(t_args->m1, t_args->m2), max(t_args->m3, t_args->m4));

	float ref = (t_args->m1 + t_args->m2 + t_args->m3 + t_args->m4) / 4.f;

	if (t_args->draw_polygon) {
		draw_polygon((float) gfit.rx, (float) gfit.ry, t_args->m1, t_args->m2, t_args->m3, t_args->m4, t_args->mr1);
	}

	siril_log_message(_("Stars: %d, Truncated mean[FWHM]: %.2f, Sensor tilt[FWHM]: %.2f (%.0f%%), Off-axis aberration[FWHM]: %.2f\n"),
			t_args->nbstars, t_args->m, worst - best, roundf(((worst - best) / ref) * 100.f), t_args->mr2 - t_args->mr1);

	free(t_args);

	return 0;
}

static int tilt_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct tilt_data *t_args = (struct tilt_data*) args->user;

	image im = { .fit = fit, .from_seq = args->seq, .index_in_seq = i };
	return compute_tilt_to_image(&im, t_args);
}

void apply_tilt_to_sequence(struct tilt_data *tilt_args) {
	struct generic_seq_args *args = create_default_seqargs(tilt_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = tilt_args->seq->selnum;
	args->prepare_hook = NULL;
	args->finalize_hook = tilt_finalize_hook;
	args->image_hook = tilt_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Tilt evaluation");
	args->has_output = FALSE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = NULL;
	args->load_new_sequence = FALSE;
	args->user = tilt_args;

	tilt_args->fit = NULL;	// not used here

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(tilt_args);
		free_generic_seq_args(args, TRUE);
	}
}


/**** show edges features **********/
static cairo_surface_t *edge_surface = NULL;
int image_width = -1;
int image_height = -1;
static char *edge_w[] = {
		"left_top",
		"center_top",
		"right_top",
		"left_center",
		"center_center",
		"right_center",
		"left_bottom",
		"center_bottom",
		"right_bottom"
};

static void set_edge_square(gchar **panel) {
	int cvport = gfit.naxes[2] > 1 ? RGB_VPORT : RED_VPORT;

	struct image_view *view = &gui.view[cvport];

	siril_open_dialog("edge_dialog");
	if (edge_surface)
		cairo_surface_destroy(edge_surface);

	image_width = gfit.rx;
	image_height = gfit.ry;
	/* New surface as we modify it */
	edge_surface = cairo_image_surface_create_for_data(view->buf, CAIRO_FORMAT_RGB24, image_width, image_height, view->full_surface_stride);

	if (cairo_surface_status(edge_surface) != CAIRO_STATUS_SUCCESS) {
		cairo_surface_destroy(edge_surface);
		edge_surface = NULL;
		return;
	}
	int widget_size = com.pref.analysis.mosaic_window / 3;
	double scale = (double) com.pref.analysis.mosaic_panel / widget_size;
	if (scale < 1.0) scale = 1.0;
	cairo_surface_set_device_scale(edge_surface, scale, scale);
	image_width = (int) ((double)image_width / scale);
	image_height = (int) ((double) image_height / scale);


	for (int i = 0; i < G_N_ELEMENTS(edge_w); i++)
		gtk_widget_queue_draw(lookup_widget(panel[i]));
}

gboolean on_left_top_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, 0, 0);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_center_top_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, (area_width - image_width) * 0.5, 0);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_right_top_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, area_width - image_width, 0);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_left_center_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, 0, (area_height - image_height) * 0.5);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_center_center_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, (area_width - image_width) * 0.5, (area_height - image_height) * 0.5);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_right_center_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, area_width - image_width, (area_height - image_height) * 0.5);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_left_bottom_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, 0, area_height - image_height);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_center_bottom_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, (area_width - image_width) * 0.5, area_height - image_height);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

gboolean on_right_bottom_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int area_width = gtk_widget_get_allocated_width (widget);
	int area_height = gtk_widget_get_allocated_height (widget);

	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	cairo_set_source_surface(cr, edge_surface, area_width - image_width, area_height - image_height);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	return FALSE;
}

void compute_aberration_inspector() {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		int widget_size = com.pref.analysis.mosaic_window / 3;

		gtk_widget_set_size_request(lookup_widget("left_top"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("center_top"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("right_top"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("left_center"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("center_center"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("right_center"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("left_bottom"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("center_bottom"), widget_size, widget_size);
		gtk_widget_set_size_request(lookup_widget("right_bottom"), widget_size, widget_size);
		set_edge_square(edge_w);
	}
}

void redraw_aberration_inspector() {
	if (!gtk_widget_is_visible(lookup_widget("edge_dialog"))) return;
	compute_aberration_inspector();
}
