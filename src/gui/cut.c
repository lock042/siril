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
#include "algos/demosaicing.h"
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

void initialize_cut_struct(cut_struct *arg) {
	arg->fit = &gfit;
	arg->seq = NULL;
	arg->imgnumber = -1;
	arg->cut_start.x = -1;
	arg->cut_start.y = -1;
	arg->cut_end.x = -1;
	arg->cut_end.y = -1;
	arg->cut_wn1.x = -1;
	arg->cut_wn1.y = -1;
	arg->cut_wn2.x = -1;
	arg->cut_wn2.y = -1;
	arg->cut_measure = FALSE;
	arg->wavenumber1 = -1;
	arg->wavenumber2 = -1;
	arg->tri = FALSE;
	arg->mode = MONO;
	arg->width = 1;
	arg->step = 1;
	arg->display_graph = TRUE;
	arg->filename = NULL;
	arg->save_dat = FALSE;
}

int sign(double x) {
	return x < 0. ? -1 : x > 0. ? 1 : 0;
}

void measure_line(cut_struct *arg, point start, point finish) {
	int deg = -1;
	point delta = { finish.x - start.x, finish.y - start.y };
	double pixdist = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (pixdist == 0.) return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	gboolean unit_is_as = (arg->fit->focal_length > 0.0) && (arg->fit->pixel_size_x > 0.0) && (arg->fit->pixel_size_y == arg->fit->pixel_size_x);
	if (unit_is_as) {
		double bin_X = com.pref.binning_update ? (double) arg->fit->binning_x : 1.0;
		double conversionfactor = (((3600.0 * 180.0) / M_PI) / 1.0E3 * (double) arg->fit->pixel_size_x / arg->fit->focal_length) * bin_X;
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


gboolean spectroscopy_selections_are_valid(cut_struct *arg) {
	gboolean a = (arg->wavenumber1 != arg->wavenumber2) && (arg->wavenumber1 > 0.0) && (arg->wavenumber2 > 0.0);
	gboolean b = (arg->cut_wn1.x >= 0) && (arg->cut_wn1.y >= 0) && (arg->cut_wn2.x >= 0) && (arg->cut_wn2.y >= 0) && (arg->cut_wn1.x < arg->fit->rx) && (arg->cut_wn1.y < arg->fit->ry) && (arg->cut_wn2.x < arg->fit->rx) && (arg->cut_wn2.y < arg->fit->ry);
	gboolean c = (!((arg->cut_wn1.x == arg->cut_wn2.x) && (arg->cut_wn1.y == arg->cut_wn2.y)));
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

void calc_zero_and_spacing(cut_struct *arg, double *zero, double *spectro_spacing) {
	point wndelta = { (double) arg->cut_wn2.x - arg->cut_wn1.x , (double) arg->cut_wn2.y - arg->cut_wn1.y };
	double wndiff_dist = sqrt(wndelta.x * wndelta.x + wndelta.y * wndelta.y);
	double wndiff = arg->wavenumber2 - arg->wavenumber1;
	*spectro_spacing = wndiff / wndiff_dist;

	// To calculate the zero we will work from whichever of x or y has the biggest difference
	double z2_z1 = wndelta.y > wndelta.x ? wndelta.y : wndelta.x;
	double z1_z0 = wndelta.y > wndelta.x ? arg->cut_wn1.y - arg->cut_start.y : arg->cut_wn1.x - arg->cut_start.x;
	*zero = arg->wavenumber1 - ( ( z1_z0 * wndiff ) / z2_z1 );
	return;
}

gpointer cut_profile(gpointer p) {
	cut_struct* arg = (cut_struct*) p;
	int retval = 0;
	gnuplot_ctrl *gplot = NULL;
	gboolean tmpfile = FALSE;
	char *filename = NULL, *tempfilename = NULL, *imagefilename = NULL;
	gchar *temp = NULL;
	double starty = arg->fit->ry - 1 - arg->cut_start.y;
	double endy = arg->fit->ry - 1 - arg->cut_end.y;
	gboolean use_gnuplot = gnuplot_is_available();
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), arg->filename);
	} else {
		gplot = gnuplot_init();
	}
	if (arg->filename) {
		filename = arg->filename;
	} else {
		if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
			temp = g_path_get_basename(com.uniq->filename);
			filename = g_strdup_printf("profile_%s_%s.dat", temp, build_timestamp_filename());
			g_free(temp);
		} else if (arg->seq) {
			filename = g_strdup_printf("profile_%s%.5d.dat", arg->seq->seqname, arg->imgnumber + 1);
		} else {
			filename = g_strdup_printf("profile_%s", build_timestamp_filename());
			siril_debug_print("%s\n", arg->filename);
		}
		tempfilename = profile_tmpfile();
		if (!(arg->save_dat || arg->seq))
			tmpfile = TRUE;
	}
	temp = g_path_get_basename(filename);
	imagefilename = replace_ext(temp, ".png");
	g_free(temp);
	temp = NULL;
	if (tmpfile) {
		g_free(filename);
		filename = g_strdup(tempfilename);
		free(tempfilename);
	}

	point delta;
	delta.x = arg->cut_end.x - arg->cut_start.x;
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
	if (arg->fit->naxes[2] > 1) {
		g = malloc(nbr_points * sizeof(double));
		b = malloc(nbr_points * sizeof(double));
	}
	x = malloc(nbr_points * sizeof(double));
	gboolean xscale = spectroscopy_selections_are_valid(arg);
	double zero = 0.0, spectro_spacing = 1.0;
	if (xscale)
		calc_zero_and_spacing(arg, &zero, &spectro_spacing);
	for (int i = 0 ; i < nbr_points ; i++) {
		if (xscale) {
			x[i] = zero + i * spectro_spacing;
		} else {
			x[i] = i * point_spacing;
		}
		if (hv) {
			// Horizontal / vertical, no interpolation
			r[i] = nointerp(arg->fit, arg->cut_start.x + point_spacing_x * i, starty + point_spacing_y * i, 0, arg->width, (int) point_spacing_x, (int) point_spacing_y);
		} else {
			// Neither horizontal nor vertical: interpolate
			r[i] = interp(arg->fit, (double) (arg->cut_start.x + point_spacing_x * i),
						  (double) (starty + point_spacing_y * i), 0,
						  arg->width, point_spacing_x, point_spacing_y);
		}
		if (arg->fit->naxes[2] > 1) {
			if (abs(point_spacing_x == 1.f) || abs(point_spacing_y == 1.f)) { // Horizontal, no interpolation
				g[i] = nointerp(arg->fit, arg->cut_start.x + point_spacing_x * i, starty + point_spacing_y * i, 1,
								arg->width, (int) point_spacing_x, (int) point_spacing_y);
				b[i] = nointerp(arg->fit, arg->cut_start.x + point_spacing_x * i, starty + point_spacing_y * i, 2,
								arg->width,  (int) point_spacing_x, (int) point_spacing_y);
			} else { // Neither horizontal nor vertical: interpolate (simple bilinear interpolation)
				g[i] = interp(arg->fit, (double) (arg->cut_start.x + point_spacing_x * i), (double) (starty + point_spacing_y * i), 1,
							  arg->width, point_spacing_x, point_spacing_y);
				b[i] = interp(arg->fit, (double) (arg->cut_start.x + point_spacing_x * i), (double) (starty + point_spacing_y * i), 2,
							  arg->width, point_spacing_x, point_spacing_y);
			}
		}
	}
	if (arg->fit->naxes[2] == 3 && arg->mode == MONO) {
		for (int i = 0 ; i < nbr_points ; i++) {
			r[i] = (r[i] + g[i] + b[i]) / 3.0;
		}
	}
	if (arg->fit->naxes[2] == 1 || arg->mode == MONO)
		retval = gnuplot_write_xy_dat(filename, x, r, nbr_points, "x L");
	else
		retval = gnuplot_write_xrgb_dat(filename, x, r, g, b, nbr_points, "x R G B");
	if (retval) {
		siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", filename);
		goto END;
	} else {
		if (!tmpfile)
			siril_log_message(_("%s has been saved.\n"), filename);
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
			if (arg->display_graph) {
				if (arg->fit->naxes[2] == 1)
					gnuplot_plot_xy_from_datfile(gplot, filename);
					else
					gnuplot_plot_xrgb_from_datfile(gplot, filename);
			} else {
				if (arg->fit->naxes[2] == 1) {
					gnuplot_plot_xy_datfile_to_png(gplot, filename, title, imagefilename);
				} else {
					gnuplot_plot_xrgb_datfile_to_png(gplot, filename, title, imagefilename);
				}
				siril_log_message(_("%s has been saved.\n"), imagefilename);
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
		siril_debug_print("Unlinking temporary file %s\n", filename);
		if (g_unlink(filename))
			siril_debug_print("Error in g_unlink()\n");
	}
	g_free(arg->filename);
	arg->filename = NULL;
	free(imagefilename);
	imagefilename = NULL;
	g_free(filename);
	filename = NULL;
	free(x);
	x = NULL;
	free(r);
	r = NULL;
	free(g);
	g = NULL;
	free(b);
	b = NULL;
	gboolean in_sequence = (arg->seq != NULL);
	if (arg != &gui.cut)
		free(arg);
	if (!in_sequence)
		siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

gpointer tri_cut(gpointer p) {
	cut_struct* arg = (cut_struct*) p;
	int retval = 0;
	gboolean tmpfile = FALSE;
	char *filename = NULL, *tempfilename = NULL, *imagefilename = NULL;
	gchar *temp = NULL;
	double starty = arg->fit->ry - 1 - arg->cut_start.y;
	double endy = arg->fit->ry - 1 - arg->cut_end.y;
	GtkSpinButton *spin_step = (GtkSpinButton*) lookup_widget("cut_tricut_step");
	double step = gtk_spin_button_get_value(spin_step);
	gboolean use_gnuplot = gnuplot_is_available();
	gnuplot_ctrl *gplot = NULL;
	if (use_gnuplot) {
		gplot = gnuplot_init();
	} else {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), arg->filename);
	}
	if (arg->filename) {
		filename = arg->filename;
	} else {
		if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
			temp = g_path_get_basename(com.uniq->filename);
			filename = g_strdup_printf("profile_%s_%s.dat", temp, build_timestamp_filename());
			g_free(temp);
			temp = NULL;
		} else if (arg->seq) {
			filename = g_strdup_printf("profile_%s%.5d.dat", arg->seq->seqname, arg->imgnumber + 1);
		} else {
			filename = g_strdup_printf("profile_%s", build_timestamp_filename());
			siril_debug_print("%s\n", arg->filename);
		}
		tempfilename = profile_tmpfile();
		if (!(arg->save_dat || arg->seq))
			tmpfile = TRUE;
	}
	temp = g_path_get_basename(filename);
	imagefilename = replace_ext(temp, ".png");
	g_free(temp);
	temp = NULL;
	if (tmpfile) {
		g_free(filename);
		filename = g_strdup(tempfilename);
		free(tempfilename);
	}
	point delta;
	delta.x = arg->cut_end.x - arg->cut_start.x;
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
		double offstartx = arg->cut_start.x - (offset * point_spacing_y * step);
		double offstarty = starty + (offset * point_spacing_x * step);
		for (int i = 0 ; i < nbr_points ; i++) {
			x[i] = i * point_spacing;
			if (hv) {
				// Horizontal / vertical, no interpolation
				r[offset+1][i] = nointerp(arg->fit, offstartx + point_spacing_x * i, offstarty + point_spacing_y * i, 0, arg->width, (int) point_spacing_x, (int) point_spacing_y);
			} else {
				// Neither horizontal nor vertical: interpolate
				r[offset+1][i] = interp(arg->fit, (double) (offstartx + point_spacing_x * i),
							(double) (offstarty + point_spacing_y * i), 0,
							arg->width, point_spacing_x, point_spacing_y);
			}
			if (arg->fit->naxes[2] > 1) {
				if (abs(point_spacing_x == 1.) || abs(point_spacing_y == 1.)) { // Horizontal, no interpolation
					r[offset+1][i] += nointerp(arg->fit, offstartx + point_spacing_x * i, offstarty + point_spacing_y * i, 1,
									arg->width, (int) point_spacing_x, (int) point_spacing_y);
					r[offset+1][i] += nointerp(arg->fit, offstartx + point_spacing_x * i, offstarty + point_spacing_y * i, 2,
									arg->width,  (int) point_spacing_x, (int) point_spacing_y);
				} else { // Neither horizontal nor vertical: interpolate (simple bilinear interpolation)
					r[offset+1][i] += interp(arg->fit, (double) (offstartx + point_spacing_x * i), (double) (offstarty + point_spacing_y * i), 1,
								arg->width, point_spacing_x, point_spacing_y);
					r[offset+1][i] += interp(arg->fit, (double) (offstartx + point_spacing_x * i), (double) (offstarty + point_spacing_y * i), 2,
								arg->width, point_spacing_x, point_spacing_y);
				}
			}
		}
		if (arg->fit->naxes[2] == 3) {
			for (int i = 0 ; i < nbr_points ; i++) {
				r[offset+1][i] /= 3.0;
			}
		}
	}
	gchar *titletext = g_strdup_printf("x L(-%dpx) L L(+%dpx)", (int) step, (int) step);
	retval = gnuplot_write_xrgb_dat(filename, x, r[0], r[1], r[2], nbr_points, titletext);
	g_free(titletext);
	if (retval) {
		siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", arg->filename);
		goto END;
	} else {
		if(!tmpfile)
			siril_log_message(_("%s has been saved.\n"), filename);
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
			if (arg->display_graph) {
				gnuplot_plot_xrgb_from_datfile(gplot, filename);
			} else {
				gnuplot_plot_xrgb_datfile_to_png(gplot, filename, title, imagefilename);
				siril_log_message(_("%s has been saved.\n"), imagefilename);
			}
			g_free(title);
			g_free(xlabel);
		}
	}
	else siril_log_message(_("Communicating with gnuplot failed\n"));

END:
	gnuplot_close(gplot);
	// Clean up
	if (tmpfile) {
		siril_debug_print("Unlinking temporary file %s\n", filename);
		if (g_unlink(filename))
			siril_debug_print("Error in g_unlink()\n");
	}
	g_free(filename);
	filename = NULL;
	free(imagefilename);
	imagefilename = NULL;
	free(x);
	x = NULL;
	for (int i = 0 ; i < 3 ; i++) {
		free(r[i]);
		r[i] = NULL;
	}
	gboolean in_sequence = (arg->seq != NULL);
	if (arg != &gui.cut)
		free(arg);

	if (!in_sequence)
		siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

gpointer cfa_cut(gpointer p) {
	cut_struct* arg = (cut_struct*) p;
	int retval = 0, ret = 0;
	gboolean tmpfile = FALSE;
	double *x = NULL, *r[4] = { 0 };
	char *filename = NULL, *tempfilename = NULL, *imagefilename = NULL;
	gchar *temp = NULL;
	fits cfa[4];
	memset(cfa, 0, 4 * sizeof(fits));
	gboolean use_gnuplot = gnuplot_is_available();
	gnuplot_ctrl *gplot = NULL;
	if (use_gnuplot) {
		gplot = gnuplot_init();
	} else {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), arg->filename);
	}

	sensor_pattern pattern = get_cfa_pattern_index_from_string(arg->fit->bayer_pattern);
	if ((arg->fit->naxes[2] > 1) || ((!(pattern == BAYER_FILTER_RGGB || pattern == BAYER_FILTER_GRBG || pattern == BAYER_FILTER_BGGR || pattern == BAYER_FILTER_GBRG)))) {
		siril_log_color_message(_("Error: CFA mode cannot be used with color images or mono images with no Bayer pattern.\n"), "red");
		retval = 1;
		goto END;
	}

	// Split arg->fit into 4 x Bayer sub-patterns cfa[0123]
	if (arg->fit->type == DATA_USHORT) {
		if ((ret = split_cfa_ushort(arg->fit, &cfa[0], &cfa[1], &cfa[2], &cfa[3]))) {
			siril_log_color_message(_("Error: failed to split FITS into CFA sub-patterns.\n"), "red");
			retval = 1;
			goto END;
		}
	} else {
		if ((ret = split_cfa_float(arg->fit, &cfa[0], &cfa[1], &cfa[2], &cfa[3]))) {
			siril_log_color_message(_("Error: failed to split FITS into CFA sub-patterns.\n"), "red");
			retval = 1;
			goto END;
		}
	}
	if (arg->filename) {
		filename = arg->filename;
	} else {
		if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
			temp = g_path_get_basename(com.uniq->filename);
			filename = g_strdup_printf("profile_%s_%s.dat", temp, build_timestamp_filename());
			g_free(temp);
		} else if (arg->seq) {
			filename = g_strdup_printf("profile_%s%.5d.dat", arg->seq->seqname, arg->imgnumber + 1);
		} else {
			filename = g_strdup_printf("profile_%s", build_timestamp_filename());
			siril_debug_print("%s\n", arg->filename);
		}
		tempfilename = profile_tmpfile();
		if (!(arg->save_dat || arg->seq))
			tmpfile = TRUE;
	}
	temp = g_path_get_basename(filename);
	imagefilename = replace_ext(temp, ".png");
	g_free(temp);
	temp = NULL;
	if (tmpfile) {
		g_free(filename);
		filename = g_strdup(tempfilename);
		free(tempfilename);
	}

	// Coordinates of profile start and endpoints in CFA space
	point cut_start_cfa = { (int) (arg->cut_start.x / 2) , (int) (arg->cut_start.y / 2) };
	point cut_end_cfa = { (int) (arg->cut_end.x / 2) , (int) (arg->cut_end.y / 2) };
	double starty = (arg->fit->ry / 2) - 1 - cut_start_cfa.y;
	double endy = (arg->fit->ry / 2) - 1 - cut_end_cfa.y;
	point delta = { cut_end_cfa.x - cut_start_cfa.x, endy - starty };
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
	retval = gnuplot_write_xcfa_dat(filename, x, r[0], r[1], r[2], r[3], nbr_points, titletext);
	g_free(titletext);
	if (retval) {
		siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", filename);
		retval = 1;
		goto END;
	} else {
		if (!tmpfile)
			siril_log_message(_("%s has been saved.\n"), filename);
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
			if (arg->display_graph) {
				gnuplot_plot_xcfa_from_datfile(gplot, filename);
			} else {
				gnuplot_plot_xcfa_datfile_to_png(gplot, filename, title, imagefilename);
				siril_log_message(_("%s has been saved.\n"), imagefilename);
			}
			g_free(title);
			g_free(xlabel);
		}
	}
	else siril_log_message(_("Communicating with gnuplot failed\n"));

END:
	gnuplot_close(gplot);
	// Clean up
	if (tmpfile) {
		siril_debug_print("Unlinking temporary file %s\n", filename);
		if (g_unlink(filename))
			siril_debug_print("Error in g_unlink()\n");
	}
	g_free(filename);
	arg->filename = NULL;
	free(imagefilename);
	imagefilename = NULL;
	free(x);
	x = NULL;
	for (int i = 0 ; i < 4 ; i++) {
		free(r[i]);
		r[i] = NULL;
		clearfits(&cfa[i]);
	}
	gboolean in_sequence = (arg->seq != NULL);
	if (arg != &gui.cut)
		free(arg);
	if (!in_sequence)
		siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

static void update_spectro_coords() {
	GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	gtk_spin_button_set_value(startx, gui.cut.cut_start.x);
	gtk_spin_button_set_value(starty, gui.cut.cut_start.y);
	gtk_spin_button_set_value(finishx, gui.cut.cut_end.x);
	gtk_spin_button_set_value(finishy, gui.cut.cut_end.y);
}

//// GUI callbacks ////

void on_cut_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
	gui.cut.fit = &gfit;
	sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit.bayer_pattern);
	gui.cut.seq = NULL;
	if (gui.cut.tri) {
		siril_debug_print("Tri-profile\n");
		start_in_new_thread(tri_cut, &gui.cut);
	} else if (gui.cut.cfa) {
		if ((gfit.naxes[2] > 1) || ((!(pattern == BAYER_FILTER_RGGB || pattern == BAYER_FILTER_GRBG || pattern == BAYER_FILTER_BGGR || pattern == BAYER_FILTER_GBRG))))
			siril_message_dialog(GTK_MESSAGE_ERROR,
					_("Invalid file for this profiling mode"),
					_("CFA mode requires a mono image file with a Bayer pattern (before debayering is carried out)."));
		else {
			siril_debug_print("CFA profiling\n");
			start_in_new_thread(cfa_cut, &gui.cut);
		}
	} else {
		if (gtk_toggle_button_get_active(cut_color))
			gui.cut.mode = COLOR;
		else
			gui.cut.mode = MONO;
		gui.cut.display_graph = TRUE;
		start_in_new_thread(cut_profile, &gui.cut);
	}
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
	gtk_adjustment_set_upper(sxa, gui.cut.fit->rx);
	gtk_adjustment_set_upper(fxa, gui.cut.fit->rx);
	gtk_adjustment_set_upper(sya, gui.cut.fit->ry);
	gtk_adjustment_set_upper(fya, gui.cut.fit->ry);
}

void on_cut_manual_coords_button_clicked(GtkButton* button, gpointer user_data) {
	match_adjustments_to_gfit();
	g_signal_handlers_block_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
	GtkWidget *cut_coords_dialog = lookup_widget("cut_coords_dialog");
	if (gui.cut.cut_start.x != -1) // If there is a cut line already made, show the
							   // endpoint coordinates in the dialog
		update_spectro_coords(&gui.cut);
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
	gtk_widget_set_sensitive(colorbutton, (gui.cut.fit->naxes[2] == 3));
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
	gui.cut.cut_start.x = sx;
	gui.cut.cut_start.y = sy;
	gui.cut.cut_end.x = fx;
	gui.cut.cut_end.y = fy;
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
	gui.cut.width = n;
}

void on_cut_spectro_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton* cut_spin_wavenumber1 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber1");
	GtkSpinButton* cut_spin_wavenumber2 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber2");
	GtkSpinButton* cut_spin_width = (GtkSpinButton*) lookup_widget("cut_spin_width");
	gui.cut.wavenumber1 = gtk_spin_button_get_value(cut_spin_wavenumber1);
	gui.cut.wavenumber2 = gtk_spin_button_get_value(cut_spin_wavenumber2);
	gui.cut.width = (int) gtk_spin_button_get_value(cut_spin_width);
	siril_close_dialog("cut_spectroscopy_dialog");
}

void on_start_select_from_star_clicked(GtkToolButton *button, gpointer user_data) {
	psf_star *result = NULL;
	int layer = 0; // Detect stars in layer 0 (mono or red) as this will always be present

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		result = psf_get_minimisation(gui.cut.fit, layer, &com.selection, FALSE, NULL, FALSE, com.pref.starfinder_conf.profile, NULL);
		set_cursor_waiting(FALSE);
		if (result) {
			gui.cut.cut_start.x = result->x0 + com.selection.x;
			gui.cut.cut_start.y = com.selection.y + com.selection.h - result->y0;
			GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
			GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
			gtk_spin_button_set_value(startx, gui.cut.cut_start.x);
			gtk_spin_button_set_value(starty, gui.cut.cut_start.y);
			redraw(REDRAW_OVERLAY);
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
		result = psf_get_minimisation(gui.cut.fit, layer, &com.selection, FALSE, NULL, FALSE, com.pref.starfinder_conf.profile, NULL);
		set_cursor_waiting(FALSE);
		if (result) {
			gui.cut.cut_end.x = result->x0 + com.selection.x;
			gui.cut.cut_end.y = com.selection.y + com.selection.h - result->y0;
			GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
			GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
			gtk_spin_button_set_value(finishx, gui.cut.cut_end.x);
			gtk_spin_button_set_value(finishy, gui.cut.cut_end.y);
			redraw(REDRAW_OVERLAY);
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
	gui.cut.cut_measure = gtk_toggle_button_get_active(button);
}

void on_cut_coords_measure_button_clicked(GtkButton *button, gpointer user_data) {
	measure_line(&gui.cut, gui.cut.cut_start, gui.cut.cut_end);
}

void on_cut_tricut_step_value_changed(GtkSpinButton *button, gpointer user_data) {
	int n = (int) gtk_spin_button_get_value(button);
	if (!(n % 2)) {
		n++;
		gtk_spin_button_set_value(button, n);
	}
	gui.cut.step = n;
	redraw(REDRAW_OVERLAY);
}

void on_cut_tri_cut_toggled(GtkToggleButton *button, gpointer user_data) {
	gui.cut.tri = gtk_toggle_button_get_active(button);
	if (gui.cut.tri) {
		gui.cut.cfa = FALSE;
	}
	printf("cfa: %d tri: %d\n", gui.cut.cfa, gui.cut.tri);
	redraw(REDRAW_OVERLAY);
}

void on_cut_cfa_toggled(GtkToggleButton *button, gpointer user_data) {
	gui.cut.cfa = gtk_toggle_button_get_active(button);
	if (gui.cut.cfa) {

		gui.cut.tri = FALSE;
	}
	printf("cfa: %d tri: %d\n", gui.cut.cfa, gui.cut.tri);
	redraw(REDRAW_OVERLAY);
}

void on_cut_save_checkbutton_toggled(GtkToggleButton *button, gpointer user_data) {
	gui.cut.save_dat = gtk_toggle_button_get_active(button);
}

void on_cut_button_toggled(GtkToggleToolButton *button, gpointer user_data) {
	if ((!gnuplot_is_available()) && (!(no_gnuplot_warning_given))) {
		siril_message_dialog(GTK_MESSAGE_WARNING,
					_("GNUplot not available"),
					_("This functionality relies on GNUplot to generate graphs. Without it, profile data files can still be produced but you will need to use external software to display them."));
		no_gnuplot_warning_given = TRUE;
	}
}

//// Sequence Processing ////

static int cut_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	cut_struct* arg = (cut_struct*) args->user;
	unsigned int MB_per_input_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_input_image, NULL, &MB_avail);
	if (arg->cfa)
		required = 2 * MB_per_input_image;
	else
		required = MB_per_input_image;
	// TODO: this doesn't account for the memory taken up by each copy of the GNUplot
	// process, which feels dangerous...
	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		limit = thread_limit;
	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);
		g_free(mem_per_thread);
		g_free(mem_available);
	}
	return limit;
}

gboolean cut_idle_function(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	cut_struct *data = (cut_struct *) args->user;
	free(data);
	gboolean retval = end_generic_sequence(args);
	return retval;
}

int cut_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	cut_struct *cut_args = (cut_struct*) args->user;
	int ret = 0;
	cut_struct *private_args = malloc(sizeof(cut_struct));
	memcpy(private_args, cut_args, sizeof(cut_struct));
	private_args->fit = fit;
	private_args->imgnumber = i;
	if (private_args->cfa)
		ret = GPOINTER_TO_INT(cfa_cut(private_args));
	else if (private_args->tri)
		ret = GPOINTER_TO_INT(tri_cut(private_args));
	else
		ret = GPOINTER_TO_INT(cut_profile(private_args));
	printf("image_worker ret: %d\n", ret);
	return ret;
}

void apply_cut_to_sequence(cut_struct* cut_args) {
	struct generic_seq_args *args = create_default_seqargs(cut_args->seq);
	args->seq = cut_args->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = args->seq->selnum;
	args->compute_mem_limits_hook = cut_compute_mem_limits;
	args->image_hook = cut_image_hook;
	args->idle_function = cut_idle_function;
	args->description = _("Intensity Profile");
	args->has_output = FALSE;
	args->load_new_sequence = FALSE;
	args->force_ser_output = FALSE;
	args->stop_on_error = FALSE;
	args->user = cut_args;

	start_in_new_thread(generic_sequence_worker, args);
}

