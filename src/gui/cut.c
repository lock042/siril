#include <math.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui/image_interactions.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/gnuplot_i.h"
#include "gui/image_display.h"

double wavenumber1 = 0.0, wavenumber2 = 0.0;
int width = 1;

cut_args cut_data = { 0 };

int sign(double x) {
	return x < 0. ? -1 : x > 0. ? 1 : 0;
}

gboolean spectroscopy_selections_are_valid() {
	gboolean a = (wavenumber1 != wavenumber2) && (wavenumber1 > 0.0) && (wavenumber2 > 0.0);
	gboolean b = (com.cut_wn1.x >= 0) && (com.cut_wn1.y >= 0) && (com.cut_wn2.x >= 0) && (com.cut_wn2.y >= 0) && (com.cut_wn1.x < gfit.rx) && (com.cut_wn1.y < gfit.ry) && (com.cut_wn2.x < gfit.rx) && (com.cut_wn2.y < gfit.ry);
	gboolean c = (!((com.cut_wn1.x == com.cut_wn2.x) && (com.cut_wn1.y == com.cut_wn2.y)));
	return a && b && c;
}

double interpf(fits* fit, double x, double y, int chan) {
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
	return (double) interp;
}

double interpw(fits* fit, double x, double y, int chan) {
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
	return (double) interp;
}

double interp(fits *fit, double x, double y, int chan, int num, double dx, double dy) {
	double val = 0.0;
	int hw = (num - 1) / 2;
	for (int i = -hw ; i < hw + 1 ; i++) {
		switch (fit->type) {
			case DATA_FLOAT:
				val += interpf(fit, x + (i * dy), y + (i * dx), chan);
				break;
			case DATA_USHORT:
				val += interpw(fit, x + (i * dy), y + (i * dx), chan);
				break;
			default:
				return -9999.0;
				break;
		}
	}
	val /= num;
	return val;
}

double nointerpf(fits *fit, int x, int y, int chan) {
	int w = gfit.rx;
	int h = gfit.ry;
	double val = (double) gfit.fdata[x + y * w + w * h * chan];
	return val;
}

double nointerpw(fits *fit, int x, int y, int chan) {
	int w = gfit.rx;
	int h = gfit.ry;
	double val = (double) gfit.data[x + y * w + w * h * chan];
	return val;
}

double nointerp(fits *fit, int x, int y, int chan, int num, int dx, int dy) {
	double val = 0.0;
	int hw = (num - 1) / 2;
	for (int i = -hw ; i < hw + 1 ; i++) {
		switch (fit->type) {
			case DATA_FLOAT:
				// Note the mismatch of x and dy etc is intentional, we are integrating
				// perpendicular to the cut line
				val += (double) nointerpf(fit, x + (i * dy), y + (i * dx), chan);
				break;
			case DATA_USHORT:
				val += (double) nointerpw(fit, x + (i * dy), y + (i * dx), chan);
				break;
			default:
				return -9999.0;
				break;
		}
	}
	val /= num;
	return val;
}

void calc_zero_and_offset(double *zero, double *spectro_spacing) {
	point wndelta = { (double) com.cut_wn2.x - com.cut_wn1.x , (double) com.cut_wn2.y - com.cut_wn1.y };
	double wndiff_dist = sqrt(wndelta.x * wndelta.x + wndelta.y * wndelta.y);
	double wndiff = wavenumber2 - wavenumber1;
	*spectro_spacing = wndiff / wndiff_dist;
	wndelta.x = com.cut_wn1.x - com.cut_start.x;
	wndelta.y = com.cut_wn1.y - com.cut_start.x;
	wndiff_dist = sqrt(wndelta.x * wndelta.x + wndelta.y * wndelta.y);
	*zero = wavenumber1 - wndiff_dist * *spectro_spacing;
	return;
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
	gboolean hv = ((point_spacing_x == 1.f) || (point_spacing_y == 1.f) || (point_spacing_x == -1.f) || (point_spacing_y == -1.f));

	r = malloc(nbr_points * sizeof(double));
	if (gfit.naxes[2] > 1) {
		g = malloc(nbr_points * sizeof(double));
		b = malloc(nbr_points * sizeof(double));
	}
	x = malloc(nbr_points * sizeof(double));
	gboolean xscale = spectroscopy_selections_are_valid();
	double zero = 0.0, spectro_spacing = 1.0;
	if (xscale)
		calc_zero_and_offset(&zero, &spectro_spacing);
	for (int i = 0 ; i < nbr_points ; i++) {
		if (xscale) {
			x[i] = zero + i * spectro_spacing;
		} else {
			x[i] = i * point_spacing;
		}
		if (hv) {
			// Horizontal / vertical, no interpolation
			r[i] = nointerp(&gfit, args->start.x + point_spacing_x * i, args->start.y + point_spacing_y * i, 0, width, (int) point_spacing_x, (int) point_spacing_y);
		} else {
			// Neither horizontal nor vertical: interpolate
			r[i] = interp(&gfit, (double) (args->start.x + point_spacing_x * i),
						  (double) (args->start.y + point_spacing_y * i), 0,
						  width, point_spacing_x, point_spacing_y);
		}
		if (gfit.naxes[2] > 1) {
			if (abs(point_spacing_x == 1.f) || abs(point_spacing_y == 1.f)) { // Horizontal, no interpolation
				g[i] = nointerp(&gfit, args->start.x + point_spacing_x * i, args->start.y + point_spacing_y * i, 1,
								width, (int) point_spacing_x, (int) point_spacing_y);
				b[i] = nointerp(&gfit, args->start.x + point_spacing_x * i, args->start.y + point_spacing_y * i, 2,
								width,  (int) point_spacing_x, (int) point_spacing_y);
			} else { // Neither horizontal nor vertical: interpolate (simple bilinear interpolation)
				g[i] = interp(&gfit, (double) (args->start.x + point_spacing_x * i), (double) (args->start.y + point_spacing_y * i), 1,
							  width, point_spacing_x, point_spacing_y);
				b[i] = interp(&gfit, (double) (args->start.x + point_spacing_x * i), (double) (args->start.y + point_spacing_y * i), 2,
							  width, point_spacing_x, point_spacing_y);
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
		gnuplot_ctrl *gplot = gnuplot_init(TRUE);
		if (gplot) {
			/* Plotting cut profile */
			gchar *title = g_strdup_printf(_("Data Cut Profile"));
			gnuplot_set_title(gplot, title);
			gchar *xlabel = NULL;
			if (xscale)
				xlabel = g_strdup_printf(_("Wavenumber"));
			else
				xlabel = g_strdup_printf(_("Distance along cut"));
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
			gnuplot_close(gplot);
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
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

//// GUI callbacks ////

void on_cut_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton* cut_mono = (GtkToggleButton*)lookup_widget("cut_radio_mono");
	GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
	GtkToggleButton* cut_spectroscopy = (GtkToggleButton*)lookup_widget("cut_radio_spectroscopy");
	cut_data.start.x = com.cut_start.x;
	cut_data.finish.x = com.cut_point.x;
	cut_data.start.y = gfit.ry - 1 - com.cut_start.y;
	cut_data.finish.y = gfit.ry - 1 - com.cut_point.y;
	if (gtk_toggle_button_get_active(cut_mono))
		cut_data.mode = MONO;
	else if (gtk_toggle_button_get_active(cut_color))
		cut_data.mode = COLOR;
	else if (gtk_toggle_button_get_active(cut_spectroscopy))
		cut_data.mode = SPECTROSCOPY;
	cut_data.display_graph = TRUE;
	start_in_new_thread(cut_profile, &cut_data);
}

void on_cut_close_button_clicked(GtkButton *button, gpointer user_data) {
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	GtkToggleToolButton *toolbutton = (GtkToggleToolButton*) lookup_widget("cut_button");
	gtk_toggle_tool_button_set_active(toolbutton, FALSE);
	siril_close_dialog("cut_dialog");
}

void match_adjustments_to_gfit() {
	GtkAdjustment *sxa = (GtkAdjustment*) lookup_adjustment("adj_cut_xstart");
	GtkAdjustment *fxa = (GtkAdjustment*) lookup_adjustment("adj_cut_xfinish");
	GtkAdjustment *sya = (GtkAdjustment*) lookup_adjustment("adj_cut_ystart");
	GtkAdjustment *fya = (GtkAdjustment*) lookup_adjustment("adj_cut_yfinish");
	gtk_adjustment_set_upper(sxa, gfit.rx);
	gtk_adjustment_set_upper(fxa, gfit.rx);
	gtk_adjustment_set_upper(sya, gfit.ry);
	gtk_adjustment_set_upper(fya, gfit.ry);
}

void on_cut_manual_coords_button_clicked(GtkButton* button, gpointer user_data) {
	match_adjustments_to_gfit();
	g_signal_handlers_block_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
	GtkWidget *cut_coords_dialog = lookup_widget("cut_coords_dialog");
	if (!gtk_widget_is_visible(cut_coords_dialog))
		siril_open_dialog("cut_coords_dialog");
}

void on_cut_spectroscopic_button_clicked(GtkButton* button, gpointer user_data) {
	g_signal_handlers_block_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
	GtkWidget *cut_spectroscopy_dialog = lookup_widget("cut_spectroscopy_dialog");
	if (!gtk_widget_is_visible(cut_spectroscopy_dialog))
		siril_open_dialog("cut_spectroscopy_dialog");
	mouse_status = MOUSE_ACTION_NONE;

}

void on_cut_spectro_cancel_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("cut_spectroscopy_dialog");
}

void on_cut_coords_cancel_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("cut_coords_dialog");
}

void on_cut_coords_dialog_hide(GtkWindow *window, gpointer user_data) {
	siril_open_dialog("cut_dialog");
	g_signal_handlers_unblock_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
}

void on_cut_spectroscopy_dialog_hide(GtkWindow *window, gpointer user_data) {
	siril_open_dialog("cut_dialog");
	mouse_status = MOUSE_ACTION_CUT_SELECT;
	g_signal_handlers_unblock_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
}

void on_cut_coords_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	int sx = (int) gtk_spin_button_get_value(startx);
	int sy = (int) gtk_spin_button_get_value(starty);
	int fx = (int) gtk_spin_button_get_value(finishx);
	int fy = (int) gtk_spin_button_get_value(finishy);
	siril_debug_print("start (%d, %d) finish (%d, %d)\n", sx, sy, fx, fy);
	com.cut_start.x = sx;
	com.cut_start.y = sy;
	com.cut_point.x = fx;
	com.cut_point.y = fy;
	redraw(REDRAW_OVERLAY);
	siril_close_dialog("cut_coords_dialog");
}

void on_cut_wavenumber1_clicked(GtkButton *button, gpointer user_data) {
	mouse_status = MOUSE_ACTION_CUT_WN1;
}
void on_cut_wavenumber2_clicked(GtkButton *button, gpointer user_data) {
	mouse_status = MOUSE_ACTION_CUT_WN2;
}

void on_cut_spectro_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton* cut_spin_wavenumber1 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber1");
	GtkSpinButton* cut_spin_wavenumber2 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber2");
	GtkSpinButton* cut_spin_width = (GtkSpinButton*) lookup_widget("cut_spin_width");
	wavenumber1 = gtk_spin_button_get_value(cut_spin_wavenumber1);
	wavenumber2 = gtk_spin_button_get_value(cut_spin_wavenumber2);
	width = (int) gtk_spin_button_get_value(cut_spin_width);
	siril_close_dialog("cut_spectroscopy_dialog");
}
