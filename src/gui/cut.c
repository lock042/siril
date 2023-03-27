/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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
#include "algos/extraction.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "gui/message_dialog.h"
#include "gui/image_interactions.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "algos/PSF.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/gnuplot_i.h"
#include "gui/image_display.h"

static gboolean no_gnuplot_warning_given = FALSE;

char* profile_tmpfile()
{
    static char const * tmp_filename_template = "profile_tmpdatafile_XXXXXX";
    char *              tmp_filename = NULL;
    char const        * tmp_dir = g_get_tmp_dir();
#ifndef _WIN32
    int                 unx_fd;
#endif // #ifndef _WIN32

/* Due to a Windows behavior and Mingw temp file name,
 * we escapes the special characters by inserting a '\' before them */
#ifdef _WIN32
    gchar *tmp = g_build_filename(tmp_dir, tmp_filename_template, NULL);
    tmp_filename = g_strescape(tmp, NULL);
    g_free(tmp);
#else
    tmp_filename = g_build_filename(tmp_dir, tmp_filename_template, NULL);
#endif

#ifdef _WIN32
    if (_mktemp(tmp_filename) == NULL)
    {
        return NULL;
    }
#else // #ifdef _WIN32
    unx_fd = mkstemp(tmp_filename);
    if (unx_fd == -1)
    {
        free(tmp_filename);
        return NULL;
    }
    close(unx_fd);

#endif // #ifdef _WIN32

    return tmp_filename;
}

void initialize_com_cut() {
	com.cut.fit = &gfit;
	com.cut.cut_start.x = -1;
	com.cut.cut_start.y = -1;
	com.cut.cut_end.x = -1;
	com.cut.cut_end.y = -1;
	com.cut.cut_wn1.x = -1;
	com.cut.cut_wn1.y = -1;
	com.cut.cut_wn2.x = -1;
	com.cut.cut_wn2.y = -1;
	com.cut.cut_measure = FALSE;
	com.cut.wavenumber1 = -1;
	com.cut.wavenumber2 = -1;
	com.cut.tri = FALSE;
	com.cut.mode = MONO;
	com.cut.width = 1;
	com.cut.step = 1;
	com.cut.display_graph = TRUE;
	com.cut.filename = NULL;
	com.cut.save_dat = FALSE;
}

int sign(double x) {
	return x < 0. ? -1 : x > 0. ? 1 : 0;
}

void measure_line(point start, point finish) {
	int deg = -1;
	point delta = { finish.x - start.x, finish.y - start.y };
	double pixdist = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (pixdist == 0.) return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	gboolean unit_is_as = (com.cut.fit->focal_length > 0.0) && (com.cut.fit->pixel_size_x > 0.0) && (com.cut.fit->pixel_size_y == com.cut.fit->pixel_size_x);
	if (unit_is_as) {
		double bin_X = com.pref.binning_update ? (double) com.cut.fit->binning_x : 1.0;
		double conversionfactor = (((3600.0 * 180.0) / M_PI) / 1.0E3 * (double) com.cut.fit->pixel_size_x / com.cut.fit->focal_length) * bin_X;
		double asdist = pixdist * conversionfactor;
		if (asdist < 60.0) {
			siril_log_message(_("Measurement: %.1f\"\n"), asdist);
		} else {
			int min = (int) asdist / 60;
			double sec = asdist - (min * 60);
			if (asdist < 3600) {
				siril_log_message(_("Measurement: %d\' %.1f\"\n"), min, sec);
			} else {
				deg = (int) asdist / 3600;
				min -= (deg * 60);
				siril_log_message(_("Measurement: %dº %d\' %.0f\"\n"), deg, min, sec);
			}
		}
	} else {
		siril_log_message(_("Measurement: %.1f px\n"), pixdist);
	}
	if (deg > 10.0)
		siril_log_color_message(_("Warning: angular measurement > 10º. Error is > 1%\n"), "salmon");
}


gboolean spectroscopy_selections_are_valid() {
	gboolean a = (com.cut.wavenumber1 != com.cut.wavenumber2) && (com.cut.wavenumber1 > 0.0) && (com.cut.wavenumber2 > 0.0);
	gboolean b = (com.cut.cut_wn1.x >= 0) && (com.cut.cut_wn1.y >= 0) && (com.cut.cut_wn2.x >= 0) && (com.cut.cut_wn2.y >= 0) && (com.cut.cut_wn1.x < com.cut.fit->rx) && (com.cut.cut_wn1.y < com.cut.fit->ry) && (com.cut.cut_wn2.x < com.cut.fit->rx) && (com.cut.cut_wn2.y < com.cut.fit->ry);
	gboolean c = (!((com.cut.cut_wn1.x == com.cut.cut_wn2.x) && (com.cut.cut_wn1.y == com.cut.cut_wn2.y)));
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
	int w = fit->rx;
	int h = fit->ry;
	double val = (double) fit->fdata[x + y * w + w * h * chan];
	return val;
}

double nointerpw(fits *fit, int x, int y, int chan) {
	int w = fit->rx;
	int h = fit->ry;
	double val = (double) fit->data[x + y * w + w * h * chan];
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

void calc_zero_and_spacing(double *zero, double *spectro_spacing) {
	point wndelta = { (double) com.cut.cut_wn2.x - com.cut.cut_wn1.x , (double) com.cut.cut_wn2.y - com.cut.cut_wn1.y };
	double wndiff_dist = sqrt(wndelta.x * wndelta.x + wndelta.y * wndelta.y);
	double wndiff = com.cut.wavenumber2 - com.cut.wavenumber1;
	*spectro_spacing = wndiff / wndiff_dist;

	// To calculate the zero we will work from whichever of x or y has the biggest difference
	double z2_z1 = wndelta.y > wndelta.x ? wndelta.y : wndelta.x;
	double z1_z0 = wndelta.y > wndelta.x ? com.cut.cut_wn1.y - com.cut.cut_start.y : com.cut.cut_wn1.x - com.cut.cut_start.x;
	*zero = com.cut.wavenumber1 - ( ( z1_z0 * wndiff ) / z2_z1 );
	// printf("zero %.3f spacing %.3f\n", *zero, *spectro_spacing);
	return;
}

gpointer cut_profile(gpointer p) {
	int retval = 0;
	gnuplot_ctrl *gplot = NULL;
	gboolean tmpfile = FALSE;
	double starty = com.cut.fit->ry - 1 - com.cut.cut_start.y;
	double endy = com.cut.fit->ry - 1 - com.cut.cut_end.y;
	gboolean use_gnuplot = gnuplot_is_available();
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), com.cut.filename);
	} else {
		gplot = gnuplot_init(TRUE);
	}
	if (!com.cut.filename) {
		siril_debug_print("Generating filename ");
		if (com.cut.save_dat) {
		siril_debug_print("for storage: ");
			if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
				siril_debug_print("from current image: ");
				gchar* temp = g_path_get_basename(com.uniq->filename);
				gchar* temp2 = g_strdup_printf("%s%s", build_timestamp_filename(), temp);
				com.cut.filename = replace_ext(temp2, ".dat");
				g_free(temp);
				g_free(temp2);
				siril_debug_print("%s\n", com.cut.filename);
			} /*else if (sequence_is_loaded()) {
				com.cut.filename = g_strdup_printf("%s%s%.5d", build_timestamp_filename(), args->seq->seqname, args->imgnumber + 1);
			}*/ else {
				com.cut.filename = g_strdup_printf("%s_image", build_timestamp_filename());
				siril_debug_print("%s\n", com.cut.filename);
			}
		} else {
			siril_debug_print("for temporary use: ");
			gchar* temp = profile_tmpfile();
			com.cut.filename = g_strdup_printf("%s.dat", temp);
			tmpfile = TRUE;
			g_free(temp);
			siril_debug_print("%s\n", com.cut.filename);
		}
	}
	point delta;
	delta.x = com.cut.cut_end.x - com.cut.cut_start.x;
	delta.y = endy - starty;
	double *x = NULL, *r = NULL, *g = NULL, *b = NULL;
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.f) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	printf("delta %f %f length %f nbr_points %d \n", delta.x, delta.y, length, nbr_points);
	gboolean hv = ((point_spacing_x == 1.f) || (point_spacing_y == 1.f) || (point_spacing_x == -1.f) || (point_spacing_y == -1.f));

	r = malloc(nbr_points * sizeof(double));
	if (com.cut.fit->naxes[2] > 1) {
		g = malloc(nbr_points * sizeof(double));
		b = malloc(nbr_points * sizeof(double));
	}
	x = malloc(nbr_points * sizeof(double));
	gboolean xscale = spectroscopy_selections_are_valid();
	double zero = 0.0, spectro_spacing = 1.0;
	if (xscale)
		calc_zero_and_spacing(&zero, &spectro_spacing);
	for (int i = 0 ; i < nbr_points ; i++) {
		if (xscale) {
			x[i] = zero + i * spectro_spacing;
		} else {
			x[i] = i * point_spacing;
		}
		if (hv) {
			// Horizontal / vertical, no interpolation
			r[i] = nointerp(com.cut.fit, com.cut.cut_start.x + point_spacing_x * i, starty + point_spacing_y * i, 0, com.cut.width, (int) point_spacing_x, (int) point_spacing_y);
		} else {
			// Neither horizontal nor vertical: interpolate
			r[i] = interp(com.cut.fit, (double) (com.cut.cut_start.x + point_spacing_x * i),
						  (double) (starty + point_spacing_y * i), 0,
						  com.cut.width, point_spacing_x, point_spacing_y);
		}
		if (com.cut.fit->naxes[2] > 1) {
			if (abs(point_spacing_x == 1.f) || abs(point_spacing_y == 1.f)) { // Horizontal, no interpolation
				g[i] = nointerp(com.cut.fit, com.cut.cut_start.x + point_spacing_x * i, starty + point_spacing_y * i, 1,
								com.cut.width, (int) point_spacing_x, (int) point_spacing_y);
				b[i] = nointerp(com.cut.fit, com.cut.cut_start.x + point_spacing_x * i, starty + point_spacing_y * i, 2,
								com.cut.width,  (int) point_spacing_x, (int) point_spacing_y);
			} else { // Neither horizontal nor vertical: interpolate (simple bilinear interpolation)
				g[i] = interp(com.cut.fit, (double) (com.cut.cut_start.x + point_spacing_x * i), (double) (starty + point_spacing_y * i), 1,
							  com.cut.width, point_spacing_x, point_spacing_y);
				b[i] = interp(com.cut.fit, (double) (com.cut.cut_start.x + point_spacing_x * i), (double) (starty + point_spacing_y * i), 2,
							  com.cut.width, point_spacing_x, point_spacing_y);
			}
		}
	}
	if (com.cut.fit->naxes[2] == 3 && com.cut.mode == MONO) {
		for (int i = 0 ; i < nbr_points ; i++) {
			r[i] = (r[i] + g[i] + b[i]) / 3.0;
		}
	}
	if (com.cut.fit->naxes[2] == 1 || com.cut.mode == MONO)
		retval = gnuplot_write_xy_dat(com.cut.filename, x, r, nbr_points, "x L");
	else
		retval = gnuplot_write_xrgb_dat(com.cut.filename, x, r, g, b, nbr_points, "x R G B");
	if (retval) {
		if (com.script)
			siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", com.cut.filename);
	} else {
		siril_log_message(_("%s has been saved.\n"), com.cut.filename);
	}
	if (use_gnuplot) {
		if (gplot) {
			// Add tmpfile to handle so it is deleted automatically
			/* Plotting cut profile */
			gchar *xlabel = NULL, *title = NULL;
			if (xscale) {
				title = g_strdup_printf(_("Spectrogram"));
				xlabel = g_strdup_printf(_("Wavenumber"));
			} else {
				title = g_strdup_printf(_("Data Cut Profile"));
				xlabel = g_strdup_printf(_("Distance along cut"));
			}
			gnuplot_set_title(gplot, title);
			gnuplot_set_xlabel(gplot, xlabel);
			gnuplot_setstyle(gplot, "lines");
			if (com.cut.display_graph) {
				if (com.cut.fit->naxes[2] == 1)
					gnuplot_plot_xy_from_datfile(gplot, com.cut.filename);
					else
					gnuplot_plot_xrgb_from_datfile(gplot, com.cut.filename);
			} else {
				gchar *imagename = replace_ext(g_path_get_basename(com.cut.filename), ".png");
				if (com.cut.fit->naxes[2] == 1)
					gnuplot_plot_xy_datfile_to_png(gplot, com.cut.filename, "test", imagename);
				else
					gnuplot_plot_xrgb_datfile_to_png(gplot, com.cut.filename, "test", imagename);
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
	if (tmpfile) {
		g_unlink(com.cut.filename);
	}
	g_free(com.cut.filename);
	com.cut.filename = NULL;
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

gpointer tri_cut(gpointer p) {
	int retval = 0;
	gboolean tmpfile = FALSE;
	double starty = com.cut.fit->ry - 1 - com.cut.cut_start.y;
	double endy = com.cut.fit->ry - 1 - com.cut.cut_end.y;
	GtkSpinButton *spin_step = (GtkSpinButton*) lookup_widget("cut_tricut_step");
	double step = gtk_spin_button_get_value(spin_step);
	gboolean use_gnuplot = gnuplot_is_available();
	gnuplot_ctrl *gplot = NULL;
	if (use_gnuplot) {
		gplot = gnuplot_init(TRUE);
	} else {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), com.cut.filename);
	}
	if (!com.cut.filename) {
		siril_debug_print("Generating filename ");
		if (com.cut.save_dat) {
		siril_debug_print("for storage: ");
			if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
				siril_debug_print("from current image: ");
				gchar* temp = g_path_get_basename(com.uniq->filename);
				gchar* temp2 = g_strdup_printf("%s%s", build_timestamp_filename(), temp);
				com.cut.filename = replace_ext(temp2, ".dat");
				g_free(temp);
				g_free(temp2);
				siril_debug_print("%s\n", com.cut.filename);
			} /*else if (sequence_is_loaded()) {
				com.cut.filename = g_strdup_printf("%s%s%.5d", build_timestamp_filename(), args->seq->seqname, args->imgnumber + 1);
			}*/ else {
				com.cut.filename = g_strdup_printf("%s_image", build_timestamp_filename());
				siril_debug_print("%s\n", com.cut.filename);
			}
		} else {
			siril_debug_print("for temporary use: ");
			gchar* temp = profile_tmpfile();
			com.cut.filename = g_strdup_printf("%s.dat", temp);
			g_free(temp);
			tmpfile = TRUE;
			siril_debug_print("%s\n", com.cut.filename);
		}
	}
	point delta;
	delta.x = com.cut.cut_end.x - com.cut.cut_start.x;
	delta.y = endy - starty;
	double *x = NULL, *r[3] = { 0 };
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.f) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	gboolean hv = ((point_spacing_x == 1.f) || (point_spacing_y == 1.f) || (point_spacing_x == -1.f) || (point_spacing_y == -1.f));
	for (int i = 0 ; i < 3 ; i++)
		r[i] = malloc(nbr_points * sizeof(double));
	x = malloc(nbr_points * sizeof(double));
	for (int offset = -1 ; offset < 2 ; offset++) {
		double offstartx = com.cut.cut_start.x - (offset * point_spacing_y * step);
		double offstarty = starty + (offset * point_spacing_x * step);
		for (int i = 0 ; i < nbr_points ; i++) {
			x[i] = i * point_spacing;
			if (hv) {
				// Horizontal / vertical, no interpolation
				r[offset+1][i] = nointerp(com.cut.fit, offstartx + point_spacing_x * i, offstarty + point_spacing_y * i, 0, com.cut.width, (int) point_spacing_x, (int) point_spacing_y);
			} else {
				// Neither horizontal nor vertical: interpolate
				r[offset+1][i] = interp(com.cut.fit, (double) (offstartx + point_spacing_x * i),
							(double) (offstarty + point_spacing_y * i), 0,
							com.cut.width, point_spacing_x, point_spacing_y);
			}
			if (com.cut.fit->naxes[2] > 1) {
				if (abs(point_spacing_x == 1.) || abs(point_spacing_y == 1.)) { // Horizontal, no interpolation
					r[offset+1][i] += nointerp(com.cut.fit, offstartx + point_spacing_x * i, offstarty + point_spacing_y * i, 1,
									com.cut.width, (int) point_spacing_x, (int) point_spacing_y);
					r[offset+1][i] += nointerp(com.cut.fit, offstartx + point_spacing_x * i, offstarty + point_spacing_y * i, 2,
									com.cut.width,  (int) point_spacing_x, (int) point_spacing_y);
				} else { // Neither horizontal nor vertical: interpolate (simple bilinear interpolation)
					r[offset+1][i] += interp(com.cut.fit, (double) (offstartx + point_spacing_x * i), (double) (offstarty + point_spacing_y * i), 1,
								com.cut.width, point_spacing_x, point_spacing_y);
					r[offset+1][i] += interp(com.cut.fit, (double) (offstartx + point_spacing_x * i), (double) (offstarty + point_spacing_y * i), 2,
								com.cut.width, point_spacing_x, point_spacing_y);
				}
			}
		}
		if (com.cut.fit->naxes[2] == 3) {
			for (int i = 0 ; i < nbr_points ; i++) {
				r[offset+1][i] /= 3.0;
			}
		}
	}
	gchar *titletext = g_strdup_printf("x L(-%dpx) L L(+%dpx)", (int) step, (int) step);
	retval = gnuplot_write_xrgb_dat(com.cut.filename, x, r[0], r[1], r[2], nbr_points, titletext);
	g_free(titletext);
	if (retval) {
		if (com.script)
			siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", com.cut.filename);
	} else {
		siril_log_message(_("%s has been saved.\n"), com.cut.filename);
	}
	if (use_gnuplot) {
		if (gplot) {
			/* Plotting cut profile */
			gchar *xlabel = NULL, *title = NULL;
			title = g_strdup_printf(_("Data Cut Profile"));
			xlabel = g_strdup_printf(_("Distance along cut / px"));
			gnuplot_set_title(gplot, title);
			gnuplot_set_xlabel(gplot, xlabel);
			gnuplot_setstyle(gplot, "lines");
			if (com.cut.display_graph) {
				gnuplot_plot_xrgb_from_datfile(gplot, com.cut.filename);
			} else {
				gchar *imagename = replace_ext(g_path_get_basename(com.cut.filename), ".png");
				gnuplot_plot_xrgb_datfile_to_png(gplot, com.cut.filename, title, imagename);
				g_free(imagename);
			}
			com.cut.filename = NULL;
			g_free(title);
			g_free(xlabel);
		}
	}
	else siril_log_message(_("Communicating with gnuplot failed\n"));

END:
	gnuplot_close(gplot);
	// Clean up
	if (tmpfile)
		g_unlink(com.cut.filename);
	g_free(com.cut.filename);
	com.cut.filename = NULL;
	free(x);
	x = NULL;
	for (int i = 0 ; i < 3 ; i++) {
		free(r[i]);
		r[i] = NULL;
	}
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

gpointer cfa_cut(gpointer p) {
	int retval = 0, ret = 0;
	gboolean tmpfile = FALSE;
	gboolean use_gnuplot = gnuplot_is_available();
	gnuplot_ctrl *gplot = NULL;

	// Split com.cut.fit into 4 x Bayer sub-patterns cfa[0123]
	fits cfa[4];
	memset(cfa, 0, 4 * sizeof(fits));
	if (com.cut.fit->type == DATA_USHORT) {
		if ((ret = split_cfa_ushort(com.cut.fit, &cfa[0], &cfa[1], &cfa[2], &cfa[3]))) {
			siril_log_color_message(_("Error: failed to split FITS into CFA sub-patterns.\n"), "red");
			return GINT_TO_POINTER(1);
		}
	} else {
		if ((ret = split_cfa_float(com.cut.fit, &cfa[0], &cfa[1], &cfa[2], &cfa[3]))) {
			siril_log_color_message(_("Error: failed to split FITS into CFA sub-patterns.\n"), "red");
			return GINT_TO_POINTER(1);
		}
	}

	if (use_gnuplot) {
		gplot = gnuplot_init(TRUE);
	} else {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), com.cut.filename);
	}

	if (!com.cut.filename) {
		siril_debug_print("Generating filename ");
		if (com.cut.save_dat) {
		siril_debug_print("for storage: ");
			if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
				siril_debug_print("from current image: ");
				gchar* temp = g_path_get_basename(com.uniq->filename);
				gchar* temp2 = g_strdup_printf("%s%s", build_timestamp_filename(), temp);
				com.cut.filename = replace_ext(temp2, ".dat");
				g_free(temp);
				g_free(temp2);
				siril_debug_print("%s\n", com.cut.filename);
			} /*else if (sequence_is_loaded()) {
				com.cut.filename = g_strdup_printf("%s%s%.5d", build_timestamp_filename(), args->seq->seqname, args->imgnumber + 1);
			}*/ else {
				com.cut.filename = g_strdup_printf("%s_image", build_timestamp_filename());
				siril_debug_print("%s\n", com.cut.filename);
			}
		} else {
			siril_debug_print("for temporary use: ");
			gchar* temp = profile_tmpfile();
			com.cut.filename = g_strdup_printf("%s.dat", temp);
			g_free(temp);
			tmpfile = TRUE;
			siril_debug_print("%s\n", com.cut.filename);
		}
	}

	// Coordinates of profile start and endpoints in CFA space
	point cut_start_cfa = { (int) (com.cut.cut_start.x / 2) , (int) (com.cut.cut_start.y / 2) };
	point cut_end_cfa = { (int) (com.cut.cut_end.x / 2) , (int) (com.cut.cut_end.y / 2) };
	double starty = (com.cut.fit->ry / 2) - 1 - cut_start_cfa.y;
	double endy = (com.cut.fit->ry / 2) - 1 - cut_end_cfa.y;
	point delta = { cut_end_cfa.x - cut_start_cfa.x, endy - starty };
	double *x = NULL, *r[4] = { 0 };
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.f) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	gboolean hv = ((point_spacing_x == 1.f) || (point_spacing_y == 1.f) || (point_spacing_x == -1.f) || (point_spacing_y == -1.f));
	for (int i = 0 ; i < 4 ; i++)
		r[i] = malloc(nbr_points * sizeof(double));
	x = malloc(nbr_points * sizeof(double));
	for (int j = 0 ; j < 4 ; j++) {
		for (int i = 0 ; i < nbr_points ; i++) {
			x[i] = i * point_spacing;
			if (hv) {
				// Horizontal / vertical, no interpolation
				r[j][i] = nointerp(&cfa[j], cut_start_cfa.x + point_spacing_x * i, starty + point_spacing_y * i,
								   0, 1, (int) point_spacing_x, (int) point_spacing_y);
			} else {
				// Neither horizontal nor vertical: interpolate
				r[j][i] = interp(&cfa[j], (double) (cut_start_cfa.x + point_spacing_x * i), (double) (starty + point_spacing_y * i),
								 0, 1, point_spacing_x, point_spacing_y);
			}
		}
	}
	gchar *titletext = g_strdup("x CFA0 CFA1 CFA2 CFA3");
	retval = gnuplot_write_xcfa_dat(com.cut.filename, x, r[0], r[1], r[2], r[3], nbr_points, titletext);
	g_free(titletext);
	if (retval) {
		if (com.script)
			siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", com.cut.filename);
	} else {
		siril_log_message(_("%s has been saved.\n"), com.cut.filename);
	}
	if (use_gnuplot) {
		if (gplot) {
			/* Plotting cut profile */
			gchar *xlabel = NULL, *title = NULL;
			title = g_strdup_printf(_("Data Cut Profile"));
			xlabel = g_strdup_printf(_("Distance along cut / px"));
			gnuplot_set_title(gplot, title);
			gnuplot_set_xlabel(gplot, xlabel);
			gnuplot_setstyle(gplot, "lines");
			if (com.cut.display_graph) {
				gnuplot_plot_xcfa_from_datfile(gplot, com.cut.filename);
			} else {
				gchar *imagename = replace_ext(g_path_get_basename(com.cut.filename), ".png");
				gnuplot_plot_xcfa_datfile_to_png(gplot, com.cut.filename, title, imagename);
				g_free(imagename);
			}
			com.cut.filename = NULL;
			g_free(title);
			g_free(xlabel);
		}
	}
	else siril_log_message(_("Communicating with gnuplot failed\n"));

END:
	gnuplot_close(gplot);
	// Clean up
	if (tmpfile)
		g_unlink(com.cut.filename);
	g_free(com.cut.filename);
	com.cut.filename = NULL;
	free(x);
	x = NULL;
	for (int i = 0 ; i < 4 ; i++) {
		free(r[i]);
		r[i] = NULL;
		clearfits(&cfa[i]);
	}
	siril_debug_print("Finished CFA plot\n");
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

static void update_spectro_coords() {
	GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	gtk_spin_button_set_value(startx, com.cut.cut_start.x);
	gtk_spin_button_set_value(starty, com.cut.cut_start.y);
	gtk_spin_button_set_value(finishx, com.cut.cut_end.x);
	gtk_spin_button_set_value(finishy, com.cut.cut_end.y);
}

//// GUI callbacks ////

void on_cut_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
	com.cut.fit = &gfit;
	if (com.cut.tri) {
		siril_debug_print("Tri-profile\n");
		start_in_new_thread(tri_cut, NULL);
	} else if (com.cut.cfa) {
		if ((com.cut.fit->naxes[2] > 1) || (com.cut.fit->bayer_pattern == NULL))
			siril_message_dialog(GTK_MESSAGE_ERROR,
					_("Invalid file for this profiling mode"),
					_("CFA mode requires a mono image file with a Bayer pattern (before debayering is carried out)."));
		else {
			siril_debug_print("CFA profiling\n");
			start_in_new_thread(cfa_cut, NULL);
		}
	} else {
		if (gtk_toggle_button_get_active(cut_color))
			com.cut.mode = COLOR;
		else
			com.cut.mode = MONO;
		com.cut.display_graph = TRUE;
		start_in_new_thread(cut_profile, NULL);
	}
}

void on_cut_close_button_clicked(GtkButton *button, gpointer user_data) {
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	GtkToggleToolButton *toolbutton = (GtkToggleToolButton*) lookup_widget("cut_button");
	gtk_toggle_tool_button_set_active(toolbutton, FALSE);
	siril_close_dialog("cut_dialog");
}

void match_adjustments_to_fit() {
	GtkAdjustment *sxa = (GtkAdjustment*) lookup_adjustment("adj_cut_xstart");
	GtkAdjustment *fxa = (GtkAdjustment*) lookup_adjustment("adj_cut_xfinish");
	GtkAdjustment *sya = (GtkAdjustment*) lookup_adjustment("adj_cut_ystart");
	GtkAdjustment *fya = (GtkAdjustment*) lookup_adjustment("adj_cut_yfinish");
	gtk_adjustment_set_upper(sxa, com.cut.fit->rx);
	gtk_adjustment_set_upper(fxa, com.cut.fit->rx);
	gtk_adjustment_set_upper(sya, com.cut.fit->ry);
	gtk_adjustment_set_upper(fya, com.cut.fit->ry);
}

void on_cut_manual_coords_button_clicked(GtkButton* button, gpointer user_data) {
	match_adjustments_to_fit();
	g_signal_handlers_block_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
	GtkWidget *cut_coords_dialog = lookup_widget("cut_coords_dialog");
	if (com.cut.cut_start.x != -1) // If there is a cut line already made, show the
							   // endpoint coordinates in the dialog
		update_spectro_coords();
	if (!gtk_widget_is_visible(cut_coords_dialog))
		siril_open_dialog("cut_coords_dialog");
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
}

void on_cut_spectroscopic_button_clicked(GtkButton* button, gpointer user_data) {
	g_signal_handlers_block_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
	GtkWidget *cut_spectroscopy_dialog = lookup_widget("cut_spectroscopy_dialog");
	if (!gtk_widget_is_visible(cut_spectroscopy_dialog))
		siril_open_dialog("cut_spectroscopy_dialog");
	mouse_status = MOUSE_ACTION_NONE;

}

void on_cut_dialog_show(GtkWindow *dialog, gpointer user_data) {
	GtkWidget* colorbutton = lookup_widget("cut_radio_color");
	gtk_widget_set_sensitive(colorbutton, (com.cut.fit->naxes[2] == 3));
}

void on_cut_spectro_cancel_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("cut_spectroscopy_dialog");
}

void on_cut_coords_cancel_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("cut_coords_dialog");
}

void on_cut_coords_dialog_hide(GtkWindow *window, gpointer user_data) {
	siril_open_dialog("cut_dialog");
	mouse_status = MOUSE_ACTION_CUT_SELECT;
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
	com.cut.cut_start.x = sx;
	com.cut.cut_start.y = sy;
	com.cut.cut_end.x = fx;
	com.cut.cut_end.y = fy;
	redraw(REDRAW_OVERLAY);
}

void on_cut_wavenumber1_clicked(GtkButton *button, gpointer user_data) {
	set_cursor("crosshair");
	mouse_status = MOUSE_ACTION_CUT_WN1;
}
void on_cut_wavenumber2_clicked(GtkButton *button, gpointer user_data) {
	set_cursor("crosshair");
	mouse_status = MOUSE_ACTION_CUT_WN2;
}

void on_cut_spin_width_value_changed(GtkSpinButton *button, gpointer user_data) {
	int n = (int) gtk_spin_button_get_value(button);
	if (!(n % 2)) {
		n++;
		gtk_spin_button_set_value(button, n);
	}
}

void on_cut_spectro_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton* cut_spin_wavenumber1 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber1");
	GtkSpinButton* cut_spin_wavenumber2 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber2");
	GtkSpinButton* cut_spin_width = (GtkSpinButton*) lookup_widget("cut_spin_width");
	com.cut.wavenumber1 = gtk_spin_button_get_value(cut_spin_wavenumber1);
	com.cut.wavenumber2 = gtk_spin_button_get_value(cut_spin_wavenumber2);
	com.cut.width = (int) gtk_spin_button_get_value(cut_spin_width);
	siril_close_dialog("cut_spectroscopy_dialog");
}

void on_start_select_from_star_clicked(GtkToolButton *button, gpointer user_data) {
	psf_star *result = NULL;
	int layer = 0; // Detect stars in layer 0 (mono or red) as this will always be present

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		result = psf_get_minimisation(com.cut.fit, layer, &com.selection, FALSE, NULL, FALSE, com.pref.starfinder_conf.profile, NULL);
		if (result) {
			com.cut.cut_start.x = result->x0 + com.selection.x;
			com.cut.cut_start.y = com.selection.y + com.selection.h - result->y0;
			GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
			GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
			gtk_spin_button_set_value(startx, com.cut.cut_start.x);
			gtk_spin_button_set_value(starty, com.cut.cut_start.y);
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR,
						_("No star detected"),
						_("Siril cannot set the start coordinate as no star has been detected in the selection"));
		}
		free_psf(result);
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR,
					_("No selection"),
					_("Siril cannot set the start coordinate as no selection is made"));
	}
}

void on_end_select_from_star_clicked(GtkToolButton *button, gpointer user_data) {
	psf_star *result = NULL;
	int layer = 0; // Detect stars in layer 0 (mono or red) as this will always be present

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		result = psf_get_minimisation(com.cut.fit, layer, &com.selection, FALSE, NULL, FALSE, com.pref.starfinder_conf.profile, NULL);
		if (result) {
			com.cut.cut_end.x = result->x0 + com.selection.x;
			com.cut.cut_end.y = com.selection.y + com.selection.h - result->y0;
			GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
			GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
			gtk_spin_button_set_value(finishx, com.cut.cut_end.x);
			gtk_spin_button_set_value(finishy, com.cut.cut_end.y);
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR,
						_("No star detected"),
						_("Siril cannot set the start coordinate as no star has been detected in the selection"));
		}
		free_psf(result);
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR,
					_("No selection"),
					_("Siril cannot set the start coordinate as no selection is made"));
	}
}

void on_cut_measure_profile_toggled(GtkToggleButton *button, gpointer user_data) {
	com.cut.cut_measure = gtk_toggle_button_get_active(button);
}

void on_cut_coords_measure_button_clicked(GtkButton *button, gpointer user_data) {
	measure_line(com.cut.cut_start, com.cut.cut_end);
}

void on_cut_tricut_step_value_changed(GtkSpinButton *button, gpointer user_data) {
	com.cut.step = gtk_spin_button_get_value(button);
	redraw(REDRAW_OVERLAY);
}

void on_cut_tri_cut_toggled(GtkToggleButton *button, gpointer user_data) {
	com.cut.tri = gtk_toggle_button_get_active(button);
	if (com.cut.tri)
		com.cut.cfa = FALSE;
	redraw(REDRAW_OVERLAY);
}

void on_cut_cfa_toggled(GtkToggleButton *button, gpointer user_data) {
	com.cut.cfa = gtk_toggle_button_get_active(button);
	if (com.cut.cfa) {
		siril_debug_print("CFA mode enabled\n");
		com.cut.tri = FALSE;
	}
	redraw(REDRAW_OVERLAY);
}

void on_cut_save_checkbutton_toggled(GtkToggleButton *button, gpointer user_data) {
	com.cut.save_dat = gtk_toggle_button_get_active(button);
}

void on_cut_button_toggled(GtkToggleToolButton *button, gpointer user_data) {
	if ((!gnuplot_is_available()) && (!(no_gnuplot_warning_given))) {
		siril_message_dialog(GTK_MESSAGE_WARNING,
					_("GNUplot not available"),
					_("This functionality relies on GNUplot to generate graphs. Without it, profile data files can still be produced but you will need to use external software to display them."));
		no_gnuplot_warning_given = TRUE;
	}
}
