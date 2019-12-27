/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "io/single_image.h"

#include "median.h"


void on_menuitem_medianfilter_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	if (single_image_is_loaded())
		siril_open_dialog("Median_dialog");
}


void on_Median_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("Median_dialog");
}

void on_Median_Apply_clicked(GtkButton *button, gpointer user_data) {
	int combo_size = gtk_combo_box_get_active(
			GTK_COMBO_BOX(
					gtk_builder_get_object(builder, "combo_ksize_median")));
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_median")));
	int iterations = round_to_int(gtk_spin_button_get_value(GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "median_button_iterations"))));

	if (get_thread_run()) {
		siril_log_message(
				_(	"Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct median_filter_data *args = malloc(sizeof(struct median_filter_data));

	switch (combo_size) {
	default:
	case 0:
		args->ksize = 3;
		break;
	case 1:
		args->ksize = 5;
		break;
	case 2:
		args->ksize = 7;
		break;
	case 3:
		args->ksize = 9;
		break;
	case 4:
		args->ksize = 11;
		break;
	case 5:
		args->ksize = 13;
		break;
	case 6:
		args->ksize = 15;
		break;
	}
	undo_save_state(&gfit, "Processing: Median Filter (filter=%dx%d px)",
			args->ksize, args->ksize);

	args->fit = &gfit;
	args->amount = amount;
	args->iterations = iterations;
	set_cursor_waiting(TRUE);
	start_in_new_thread(median_filter, args);

}

/*****************************************************************************
 *                M E D I A N     I M A G E     F I L T E R S                *
 ****************************************************************************/

/* get the median of the neighbors of pixel (xx, yy), including itself if
 * include_self is TRUE. radius is 1 for a 3x3, 2 for a 5x5 and so on.
 * w and h are the size of the image passed in buf.
 */
double get_median_ushort(WORD *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	WORD *values;
	double median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(WORD));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				// ^ limit to image bounds ^
				// v exclude centre pixel v
				if (include_self || x != xx || y != yy) {
					values[n++] = buf[x + y * w];
				}
			}
		}
	}
	median = quickmedian(values, n);
	free(values);
	return median;
}

double get_median_float(float *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	float *values;
	double median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(float));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				// ^ limit to image bounds ^
				// v exclude centre pixel v
				if (include_self || x != xx || y != yy) {
					values[n++] = buf[x + y * w];
				}
			}
		}
	}
	median = quickmedian_float(values, n);
	free(values);
	return median;
}

double get_median_gsl(gsl_matrix *mat, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	double *values, median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(double));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				// ^ limit to image bounds ^
				// v exclude centre pixel v
				if (include_self || x != xx || y != yy) {
					values[n++] = gsl_matrix_get(mat, y, x);
				}
			}
		}
	}
	median = quickmedian_double(values, n);
	free(values);
	return median;
}


/*****************************************************************************
 *                      M E D I A N     F I L T E R                          *
 ****************************************************************************/

static gboolean end_median_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

static gpointer median_filter_ushort(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	int progress = 0, x, y, layer, iter = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	double total, norm = (double)get_normalized_value(args->fit);
	struct timeval t_start, t_end;
	int radius = (args->ksize - 1) / 2;

	g_assert(args->ksize % 2 == 1 && args->ksize > 1);
	g_assert(nx > 0 && ny > 0);
	total = ny * args->fit->naxes[2] * args->iterations;

	char *msg = siril_log_color_message(_("Median Filter: processing...\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	do {
		for (layer = 0; layer < args->fit->naxes[2]; layer++) {
			WORD *data = args->fit->pdata[layer];
			for (y = 0; y < ny; y++) {
				int pix_idx = y * nx;
				if (!get_thread_run()) break;
				if (!(y % 16))	// every 16 iterations
					set_progress_bar_data(NULL, (double)progress / total);
				progress++;
				for (x = 0; x < nx; x++) {
					double median = get_median_ushort(data, x, y, nx, ny, radius, FALSE, TRUE);
					if (args->amount != 1.0) {
						double pixel = args->amount * (median / norm);
						pixel += (1.0 - args->amount)
							* ((double)data[pix_idx] / norm);
						data[pix_idx] = round_to_WORD(pixel * norm);
					} else {
						data[pix_idx] = round_to_WORD(median);
					}
					pix_idx++;
				}
			}
		}
		iter++;
	} while (iter < args->iterations && get_thread_run());
	invalidate_stats_from_fit(args->fit);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	set_progress_bar_data(_("Median filter applied"), PROGRESS_DONE);
	siril_add_idle(end_median_filter, args);

	return GINT_TO_POINTER(0);
}

static gpointer median_filter_float(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	int progress = 0, x, y, layer, iter = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	double total;
	struct timeval t_start, t_end;
	int radius = (args->ksize - 1) / 2;

	g_assert(args->ksize % 2 == 1 && args->ksize > 1);
	g_assert(nx > 0 && ny > 0);
	total = ny * args->fit->naxes[2] * args->iterations;

	char *msg = siril_log_color_message(_("Median Filter: processing...\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	do {
		for (layer = 0; layer < args->fit->naxes[2]; layer++) {
			float *data = args->fit->fpdata[layer];
			for (y = 0; y < ny; y++) {
				int pix_idx = y * nx;
				if (!get_thread_run()) break;
				if (!(y % 16))	// every 16 iterations
					set_progress_bar_data(NULL, (double)progress / total);
				progress++;
				for (x = 0; x < nx; x++) {
					double median = get_median_float(data, x, y, nx, ny, radius, FALSE, TRUE);
					if (args->amount != 1.0) {
						double pixel = args->amount * median;
						pixel += (1.0 - args->amount) * (double)data[pix_idx];
						data[pix_idx] = (float)pixel;
					} else {
						data[pix_idx] = (float)median;
					}
					pix_idx++;
				}
			}
		}
		iter++;
	} while (iter < args->iterations && get_thread_run());
	invalidate_stats_from_fit(args->fit);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	set_progress_bar_data(_("Median filter applied"), PROGRESS_DONE);
	siril_add_idle(end_median_filter, args);

	return GINT_TO_POINTER(0);
}

/* The function smoothes an image using the median filter with the
 * ksize x ksize aperture. Each channel of a multi-channel image is
 * processed independently. In-place operation is supported. */
gpointer median_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	if (args->fit->type == DATA_USHORT)
		return median_filter_ushort(p);
	if (args->fit->type == DATA_FLOAT)
		return median_filter_float(p);
	siril_add_idle(end_median_filter, args);
	return GINT_TO_POINTER(1);
}
