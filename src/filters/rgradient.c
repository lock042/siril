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
#include "core/arithm.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/PSF.h"
#include "gui/image_display.h"
#include "gui/histogram.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/registration_preview.h"
#include "opencv/opencv.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

#include "rgradient.h"

static void to_polar(int x, int y, point center, double *r, double *theta) {
	double dx = x - center.x;
	double dy = y - center.y;
	*r = sqrt(dx * dx + dy * dy);
	*theta = atan2(dy, dx);
}

static void to_cartesian(double r, double theta, point center, point *p) {
	p->x = center.x + r * cos(theta);
	p->y = center.y + r * sin(theta);
}

static gboolean end_rgradient_filter(gpointer p) {
	struct rgradient_filter_data *args = (struct rgradient_filter_data *) p;
	stop_processing_thread();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);

	free(args);
	return FALSE;
}

gpointer rgradient_filter(gpointer p) {
	struct timeval t_start, t_end;

	siril_log_color_message(_("Rotational gradient: processing...\n"), "green");
	gettimeofday(&t_start, NULL);

	struct rgradient_filter_data *args = (struct rgradient_filter_data *) p;

	gboolean was_ushort;
	fits imA = { 0 }, imB = { 0 };
	int retval = 0;
	const point center = {args->xc, args->yc};
	const double dAlpha = M_PI / 180.0 * args->da;

	int cur_nb = 0; // only used for progress bar
	const double total = args->fit->ry * args->fit->naxes[2];	// only used for progress bar
	set_progress_bar_data(_("Rotational gradient in progress..."), PROGRESS_RESET);

	was_ushort = args->fit->type == DATA_USHORT;

	/* convenient transformation to not inverse y sign */
	fits_flip_top_to_bottom(args->fit);

	if (was_ushort) {
		const long n = args->fit->naxes[0] * args->fit->naxes[1] * args->fit->naxes[2];
		float *newbuf = ushort_buffer_to_float(args->fit->data, n);
		if (!newbuf) { retval = 1; goto end_rgradient; }
		fit_replace_buffer(args->fit, newbuf, DATA_FLOAT);
	}

	retval = copyfits(args->fit, &imA, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (retval) { retval = 1; goto end_rgradient; }

	retval = copyfits(&imA, &imB, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (retval) { retval = 1; goto end_rgradient; }

	int layer, y;
	float global_min = FLT_MAX;

	// Processing the image layers
	for (layer = 0; layer < args->fit->naxes[2]; layer++) {
		float *gbuf = args->fit->fpdata[layer];
		float *Abuf = imA.fpdata[layer];
		float *Bbuf = imB.fpdata[layer];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(min:global_min)
#endif
		for (y = 0; y < args->fit->ry; y++) {
			size_t i = y * args->fit->rx;
#ifdef _OPENMP
			#pragma omp critical
			#endif
			{
				set_progress_bar_data(NULL, cur_nb / total);
				cur_nb++;
			}

			for (int x = 0; x < args->fit->rx; x++) {
				float buf = gbuf[i] + gbuf[i];

				double r, theta;
				point delta1, delta2;

				to_polar(x, y, center, &r, &theta);

				// Positive differential
				to_cartesian(r - args->dR, theta + dAlpha, center, &delta1);
				buf -= bilinear(Abuf, args->fit->rx, args->fit->ry, delta1.x, delta1.y);

				// Negative differential
				to_cartesian(r - args->dR, theta - dAlpha, center, &delta2);
				buf -= bilinear(Bbuf, args->fit->rx, args->fit->ry, delta2.x, delta2.y);

				gbuf[i] = buf > 1.f ? 1.f : buf;

				// Update global minimum value
				if (gbuf[i] < global_min) {
					global_min = gbuf[i];
				}

				i++;
			}
		}
	}

	retval = soper(args->fit, global_min, OPER_SUB, TRUE);
	if (retval) { retval = 1; goto end_rgradient; }

end_rgradient:
	fits_flip_top_to_bottom(args->fit);
	set_progress_bar_data(_("Rotational gradient complete."), PROGRESS_DONE);

	clearfits(&imA);
	clearfits(&imB);
	if (!retval) {
		if (was_ushort) {
			const long n = args->fit->naxes[0] * args->fit->naxes[1] * args->fit->naxes[2];
			WORD *newbuf = float_buffer_to_ushort(args->fit->fdata, n);
			if (!newbuf)
				retval = 1;
			else fit_replace_buffer(args->fit, newbuf, DATA_USHORT);
		}
		set_progress_bar_data(_("Rotational gradient complete."), PROGRESS_DONE);
	}

	invalidate_stats_from_fit(args->fit);
	update_gfit_histogram_if_needed();
	siril_add_idle(end_rgradient_filter, args);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	return GINT_TO_POINTER(retval);
}

/// GUI

static double get_xc() {
	GtkEntry *entry;

	entry = GTK_ENTRY(lookup_widget("entry_rgradient_xc"));

	return g_ascii_strtod(gtk_entry_get_text(entry), NULL);
}

static double get_yc() {
	GtkEntry *entry;

	entry = GTK_ENTRY(lookup_widget("entry_rgradient_yc"));

	return g_ascii_strtod(gtk_entry_get_text(entry), NULL);
}

static double get_dR() {
	GtkRange *range;

	range = GTK_RANGE(lookup_widget("scale_radial_rgradient"));

	return gtk_range_get_value(range);
}

static double get_da() {
	GtkRange *range;

	range = GTK_RANGE(lookup_widget("scale_rot_rgradient"));

	return gtk_range_get_value(range);
}

////// CALLBACKS

void on_rgradient_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("rgradient_dialog");
}

void on_rgradient_Apply_clicked(GtkButton *button, gpointer user_data) {
	if(!check_ok_if_cfa())
		return;
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (!single_image_is_loaded()) return;

	if (gfit.orig_bitpix == BYTE_IMG) {
		siril_log_color_message(_("This process cannot be applied to 8b images\n"), "red");
		return;
	}

	struct rgradient_filter_data *args = malloc(sizeof(struct rgradient_filter_data));
	args->xc = get_xc();
	args->yc = get_yc();
	args->dR = get_dR();
	args->da = get_da();
	args->fit = &gfit;

	if ((args->xc >= args->fit->rx) || (args->yc >= args->fit->ry)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong center coordinates"),
				_("The coordinates cannot be greater than the size of the image. Please change their values and retry."));
	} else {

	set_cursor_waiting(TRUE);

	undo_save_state(&gfit, _("RGradient: (dR=%5.2lf, dA=%4.2lf, xc=%7.1lf, yc=%7.1lf)"),
			args->dR, args->da, args->xc, args->yc);

	start_in_new_thread(rgradient_filter, args);
	}
}

void on_button_rgradient_selection_clicked(GtkButton *button, gpointer user_data) {
	if (com.selection.h && com.selection.w) {
		psf_star *result = psf_get_minimisation(&gfit, 0, &com.selection, FALSE, NULL, TRUE, PSF_GAUSSIAN, NULL);
		if (result) {
			gchar *x0 = g_strdup_printf("%.3lf", result->x0 + com.selection.x);
			gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_rgradient_xc")), x0);
			gchar *y0 = g_strdup_printf("%.3lf", com.selection.y + com.selection.h - result->y0);
			gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_rgradient_yc")), y0);
			g_free(x0);
			g_free(y0);
			free_psf(result);
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Center coordinate selection error"),
				_("No valid PSF found within selection."));
		}
	}
}

