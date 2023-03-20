#include <math.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui/image_interactions.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/gnuplot_i.h"

int sign(double x) {
	return x < 0. ? -1 : x > 0. ? 1 : 0;
}

float interpf(fits* fit, float x, float y, int chan) {
	if (chan >= fit->naxes[2])
		return -9999.f;
	int w = fit->rx;
	int h = fit->ry;
	int npixels = w * h;
	int x0 = (int)(x);
	int x1 = (int)(x) + 1;
	int y0 = (int)(y);
	int y1 = (int)(y) + 1;
	float val00 = fit->fdata[x0 + y0 * w + npixels * chan];
	float val01 = fit->fdata[x1 + y0 * w + npixels * chan];
	float val10 = fit->fdata[x0 + y1 * w + npixels * chan];
	float val11 = fit->fdata[x1 + y1 * w + npixels * chan];
	float interp1 = (x - x0) * val00 + (x1 - x) * val01;
	float interp2 = (x - x0) * val10 + (x1 - x) * val11;
	float interp = (y - y0) * interp1 + (y1 - y) * interp2;
	return interp;
}

float interpw(fits* fit, float x, float y, int chan) {
	if (chan >= fit->naxes[2])
		return -9999.f;
	int w = fit->rx;
	int h = fit->ry;
	int npixels = w * h;
	int x0 = (int)(x);
	int x1 = (int)(x) + 1;
	int y0 = (int)(y);
	int y1 = (int)(y) + 1;
	float val00 = (float) fit->data[x0 + y0 * w + npixels * chan];
	float val01 = (float) fit->data[x1 + y0 * w + npixels * chan];
	float val10 = (float) fit->data[x0 + y1 * w + npixels * chan];
	float val11 = (float) fit->data[x1 + y1 * w + npixels * chan];
	float interp1 = (x - x0) * val00 + (x1 - x) * val01;
	float interp2 = (x - x0) * val10 + (x1 - x) * val11;
	float interp = (y - y0) * interp1 + (y1 - y) * interp2;
	return interp;
}

float interp(fits *fit, float x, float y, int chan) {
	float val;
	switch (fit->type) {
		case DATA_FLOAT:
			val = interpf(fit, x, y, chan);
			return val;
			break;
		case DATA_USHORT:
			return interpw(fit, x, y, chan);
			break;
		default:
			return -9999.f;
			break;
	}
}

gpointer cut_profile(gpointer p) {
	cut_args *args = (cut_args *) p;
	int retval = 0;
	char *filename = "cut.dat";
	gboolean use_gnuplot = gnuplot_is_available();
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), filename);
	}
	pointi delta;
	delta.x = args->finish.x - args->start.x;
	delta.y = args->finish.y - args->start.y;
	printf("sx: %d, sy: %d, fx: %d, fy:%d\n", args->start.x, args->start.y, args->finish.x, args->finish.y);
	printf("dx: %d dy: %d\n", delta.x, delta.y);
	double *x = NULL, *r = NULL, *g = NULL, *b = NULL;
	double length = sqrtf(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.f) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;

	r = malloc(nbr_points * sizeof(double));
	if (gfit.naxes[2] > 1) {
		g = malloc(nbr_points * sizeof(double));
		b = malloc(nbr_points * sizeof(double));
	}
	x = malloc(nbr_points * sizeof(double));
	int w = gfit.rx;
	int h = gfit.ry;
	for (int i = 0 ; i < nbr_points ; i++) {
		x[i] = i * point_spacing;
		if (abs(point_spacing_x == 1.f)) { // Horizontal, no interpolation
			if (gfit.type == DATA_FLOAT)
				r[i] = gfit.fdata[args->start.x + (i * sign(point_spacing_x)) + args->start.y * w];
			else
				r[i] = gfit.data[args->start.x + (i * sign(point_spacing_x)) + args->start.y * w];
		} else if (abs(point_spacing_y == 1.f)) { // Vertical, no interpolation
			if (gfit.type == DATA_FLOAT)
				r[i] = gfit.fdata[args->start.x + (args->start.y + i * sign(point_spacing_y)) * w];
			else
				r[i] = gfit.data[args->start.x + i + (args->start.y + i * sign(point_spacing_y)) * w];
		} else { // Neither horizontal nor vertical: interpolate
			r[i] = interp(&gfit, (double) (args->start.x + point_spacing_x * i), (double) (args->start.y + point_spacing_y * i), 0);
		}
		if (gfit.naxes[2] > 1) {
			if (abs(point_spacing_x == 1.f)) { // Horizontal, no interpolation
				if (gfit.type == DATA_FLOAT) {
					g[i] = gfit.fdata[args->start.x + (i * sign(point_spacing_x)) + args->start.y * w + w * h];
					b[i] = gfit.fdata[args->start.x + (i * sign(point_spacing_x)) + args->start.y * w + w * h * 2];
				} else {
					g[i] = gfit.data[args->start.x + (i * sign(point_spacing_x)) + args->start.y * w + w * h];
					b[i] = gfit.data[args->start.x + (i * sign(point_spacing_x)) + args->start.y * w + w * h * 2];
				}
			} else if (abs(point_spacing_y == 1.f)) { // Vertical, no interpolation
				if (gfit.type == DATA_FLOAT) {
					g[i] = gfit.fdata[args->start.x + (args->start.y + i * sign(point_spacing_y)) * w + w * h];
					b[i] = gfit.fdata[args->start.x + (args->start.y + i * sign(point_spacing_y)) * w + w * h * 2];
				} else {
					g[i] = gfit.data[args->start.x + i + (args->start.y + i * sign(point_spacing_y)) * w + w * h];
					b[i] = gfit.data[args->start.x + i + (args->start.y + i * sign(point_spacing_y)) * w + w * h * 2];
				}
			} else { // Neither horizontal nor vertical: interpolate (simple bilinear interpolation)
				r[i] = interp(&gfit, (double) (args->start.x + point_spacing_x * i), (double) (args->start.y + point_spacing_y * i), 0);
				g[i] = interp(&gfit, (double) (args->start.x + point_spacing_x * i), (double) (args->start.y + point_spacing_y * i), 1);
				b[i] = interp(&gfit, (double) (args->start.x + point_spacing_x * i), (double) (args->start.y + point_spacing_y * i), 2);
			}
		}
	}
	if (gfit.naxes[2] == 3 && args->mode == MONO) {
		for (int i = 0 ; i < nbr_points ; i++) {
			r[i] = (r[i] + g[i] + b[i]) / 3.0;
		}
	}
	if (gfit.naxes[2] == 1 || args->mode == MONO)
		retval = gnuplot_write_xy_dat(filename, x, r, nbr_points, "x L");
	else
		retval = gnuplot_write_xrgb_dat(filename, x, r, g, b, nbr_points, "x R G B");
	if (retval) {
		if (com.script)
			siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", filename);
	} else {
		siril_log_message(_("%s has been saved.\n"), filename);
	}
	if (use_gnuplot) {
		gnuplot_ctrl *gplot = gnuplot_init();
		if (gplot) {
			/* Plotting cut profile */
			gchar *title = g_strdup_printf(_("Data Cut Profile"));
			gnuplot_set_title(gplot, title);
			gchar *xlabel = g_strdup_printf(_("Distance along cut"));
			gnuplot_set_xlabel(gplot, xlabel);
			gnuplot_setstyle(gplot, "lines");
			if (args->display_graph) {
				if (gfit.naxes[2] == 1)
					gnuplot_plot_xy_from_datfile(gplot, filename, "L");
				else
					gnuplot_plot_xrgb_from_datfile(gplot, filename);
			} else {
				gchar *imagename = replace_ext(filename, ".png");
				if (gfit.naxes[2] == 1)
					gnuplot_plot_xy_datfile_to_png(gplot, filename, "test", imagename);
				else
					gnuplot_plot_xrgb_datfile_to_png(gplot, filename, "test", imagename);
				g_free(imagename);
			}
			g_free(title);
			g_free(xlabel);
//			gnuplot_close(gplot);
		}
		else siril_log_message(_("Communicating with gnuplot failed\n"));
	}

END:
	// Clean up
	free(x);
	x = NULL;
	free(r);
	r = NULL;
	free(g);
	g = NULL;
	free(b);
	b = NULL;
	free(args);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

//// GUI callbacks ////

void on_cut_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton* cut_mono = (GtkToggleButton*)lookup_widget("cut_radio_mono");
	GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
	GtkToggleButton* cut_spectroscopy = (GtkToggleButton*)lookup_widget("cut_radio_spectroscopy");
	cut_args *cut_data = malloc(sizeof(cut_args));
	cut_data->start.x = com.cut_start.x;
	cut_data->finish.x = com.cut_point.x;
	// Reverse the y coordinates to look up the correct point in gfit.
	// This can't be done in the image_interactions callbacks otherwise
	// the cut preview line is drawn with the wrong orientation.
	cut_data->start.y = gfit.ry - 1 - com.cut_start.y;
	cut_data->finish.y = gfit.ry - 1 - com.cut_point.y;
	if (gtk_toggle_button_get_active(cut_mono))
		cut_data->mode = MONO;
	else if (gtk_toggle_button_get_active(cut_color))
		cut_data->mode = COLOR;
	else if (gtk_toggle_button_get_active(cut_spectroscopy))
		cut_data->mode = SPECTROSCOPY;
	cut_data->display_graph = TRUE;
	start_in_new_thread(cut_profile, cut_data);
	siril_close_dialog("cut_dialog");
}

void on_cut_cancel_button_clicked(GtkButton *button, gpointer user_data) {
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		siril_close_dialog("cut_dialog");
}

void on_cut_dialog_show(GtkButton *button, gpointer user_data) {
}

void on_cut_button_toggled(GtkToggleButton *button, gpointer user_data) {
	if (gtk_toggle_tool_button_get_active(button)) {
		mouse_status = MOUSE_ACTION_CUT_SELECT;
	} else {
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		siril_close_dialog("cut_dialog");
	}
}
