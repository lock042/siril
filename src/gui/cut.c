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

static char* profile_tmpfile()
{
    static char const * tmp_filename_template = "profile_tmpdatafile_XXXXXX";
    char *              tmp_filename = NULL;
    char const        * tmp_dir = g_get_tmp_dir();
#ifndef _WIN32
    int                 unx_fd;
#endif // #ifndef _WIN32

/* Due to a Windows behavior and Mingw temp file name,
 * we escape the special characters by inserting a '\' before them */
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

void reset_cut_gui_filedependent() { // Searated out to avoid having to repeat too much after opening a new file
	GtkWidget *colorbutton = (GtkWidget*) lookup_widget("cut_radio_color");
	GtkWidget *cfabutton = (GtkWidget*) lookup_widget("cut_cfa");
	gtk_widget_set_sensitive(colorbutton, (gfit.naxes[2] == 3));
	sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit.bayer_pattern);
	gboolean cfa_disabled = ((gfit.naxes[2] > 1) || ((!(pattern == BAYER_FILTER_RGGB || pattern == BAYER_FILTER_GRBG || pattern == BAYER_FILTER_BGGR || pattern == BAYER_FILTER_GBRG))));
		gtk_widget_set_sensitive(cfabutton, !cfa_disabled);
}

static void reset_cut_gui() {
	GtkToggleButton *radio_mono = (GtkToggleButton*) lookup_widget("cut_radio_mono");
	gtk_toggle_button_set_active(radio_mono, TRUE);
	GtkToggleButton *measure = (GtkToggleButton*) lookup_widget("cut_measure_profile");
	gtk_toggle_button_set_active(measure, FALSE);
	GtkToggleButton *save_dat = (GtkToggleButton*) lookup_widget("cut_save_checkbutton");
	gtk_toggle_button_set_active(save_dat, FALSE);
	GtkToggleButton *save_png = (GtkToggleButton*) lookup_widget("cut_save_png");
	gtk_toggle_button_set_active(save_png, FALSE);
	GtkSpinButton *cut_startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton *cut_starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton *cut_finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton *cut_finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	gtk_spin_button_set_value(cut_startx, -1);
	gtk_spin_button_set_value(cut_starty, -1);
	gtk_spin_button_set_value(cut_finishx, -1);
	gtk_spin_button_set_value(cut_finishy, -1);
	GtkSpinButton *cut_width = (GtkSpinButton*) lookup_widget("cut_spin_width");
	GtkSpinButton *cut_wn1 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber1");
	GtkSpinButton *cut_wn2 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber2");
	gtk_spin_button_set_value(cut_width, 1);
	gtk_spin_button_set_value(cut_wn1, -1);
	gtk_spin_button_set_value(cut_wn2, -1);
	GtkLabel *wn1x = (GtkLabel*) lookup_widget("label_wn1_x");
	GtkLabel *wn1y = (GtkLabel*) lookup_widget("label_wn1_y");
	GtkLabel *wn2x = (GtkLabel*) lookup_widget("label_wn2_x");
	GtkLabel *wn2y = (GtkLabel*) lookup_widget("label_wn2_y");
	gtk_label_set_text(wn1x, "");
	gtk_label_set_text(wn1y, "");
	gtk_label_set_text(wn2x, "");
	gtk_label_set_text(wn2y, "");
	reset_cut_gui_filedependent();
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
	arg->mode = CUT_MONO;
	arg->width = 1;
	arg->step = 1;
	arg->display_graph = TRUE;
	arg->filename = NULL;
	arg->save_dat = FALSE;
	arg->save_png_too = FALSE;
	arg->vport = -1;
	if (!com.script)
		reset_cut_gui();
}

void free_cut_args(cut_struct *arg) {
	if (arg->filename) {
		g_free(arg->filename);
	}
	if (arg != &gui.cut)
		free(arg);
	return;
}

gboolean cut_struct_is_valid(cut_struct *arg) {
	// define # layers
	int nb_layers = (arg->seq) ? arg->seq->nb_layers : arg->fit->naxes[2];
	// checking mutually exclusive choices
	if (arg->tri && arg->cfa) {
		siril_log_color_message(_("Error: CFA mode and tri-profile are mutually exclusive.\n"), "red");
		return FALSE;
	}
	if (arg->tri && arg->mode == CUT_COLOR) {
		siril_log_color_message(_("Error: color plot and tri-profile are mutually exclusive.\n"), "red");
		return FALSE;
	}
	if (nb_layers == 1 && arg->mode == CUT_COLOR) {
		siril_log_color_message(_("Error: mono image is loaded: color mode is unavailable.\n"), "red");
		return FALSE;
	}
	if (nb_layers == 1 && arg->vport > 0) {
		siril_log_color_message(_("Error: mono image is loaded: colored layers cannot be specified.\n"), "red");
		return FALSE;
	}
	if (arg->cfa && arg->mode == CUT_COLOR) {
		siril_log_color_message(_("Error: color plot and CFA mode are mutually exclusive.\n"), "red");
		return FALSE;
	}
		if (arg->cut_start.x == arg->cut_end.x && arg->cut_start.y == arg->cut_end.y) {
		siril_log_message(_("Error: start and finish points are the same.\n"));
		return FALSE;
	}
	if (arg->seq && arg->seq->is_variable) {
		siril_log_message(_("Error: variable image size in sequence.\n"));
		return FALSE;
	}

	// Check args are cromulent
	int rx = (arg->seq) ? arg->seq->rx : gfit.rx;
	int ry = (arg->seq) ? arg->seq->ry : gfit.ry;

	if (arg->cut_wn1.x > -1.0 && arg->cut_wn1.x == arg->cut_wn2.x && arg->cut_wn1.y == arg->cut_wn2.y) {
		siril_log_message(_("Error: wavenumber points are the same.\n"));
		return FALSE;
	}
	if (arg->cut_wn1.x > -1.0 && arg->wavenumber1 == -1.0) {
		siril_log_message(_("Error: wavenumber for -wn1= is not set.\n"));
		return FALSE;
	}
	if (arg->cut_wn2.x > -1.0 && arg->wavenumber2 == -1.0) {
		siril_log_message(_("Error: wavenumber for -wn2= is not set.\n"));
		return FALSE;
	}
	if (arg->wavenumber1 >=0.0 && (arg->cut_wn1.x < 0.0 || arg->cut_wn1.y < 0.0)) {
		siril_log_message(_("Error: wavenumber1 set but corresponding location is not set.\n"));
		return FALSE;
	}
	if (arg->wavenumber2 >=0.0 && (arg->cut_wn2.x < 0.0 || arg->cut_wn2.y < 0.0)) {
		siril_log_message(_("Error: wavenumber2 set but corresponding location is not set.\n"));
		return FALSE;
	}
	if ((arg->cut_wn1.x < 0.0 || arg->cut_wn2.y < 0.0 || arg->cut_wn2.x < 0.0 || arg->cut_wn2.y < 0.0) &&
		(!(arg->cut_wn1.x == -1.0 && arg->cut_wn1.y == -1.0 && arg->cut_wn2.x == -1.0 && arg->cut_wn2.y == -1.0))) {
		siril_log_message(_("Error: wavenumber point outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->cut_wn1.x > rx - 1 || arg->cut_wn2.x > rx - 1 || arg->cut_wn1.y > ry - 1 || arg->cut_wn2.y > ry - 1) {
		siril_log_message(_("Error: wavenumber point outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->cut_start.x > rx - 1 || arg->cut_start.y > ry - 1 || arg->cut_end.x > rx - 1 || arg->cut_end.y > ry - 1) {
		siril_log_message(_("Error: profile line endpoint outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->step != 1 && !arg->tri) {
		siril_log_message(_("Error: spacing is set but tri-profile mode is not selected.\n"));
		return FALSE;
	}
	if (arg->vport < 0 || arg->vport > nb_layers) {
		siril_log_message(_("Error: layer out of range.\n"));
		return FALSE;
	}
	return TRUE;

}
/*
static int sign(double x) {
	return x < 0. ? -1 : x > 0. ? 1 : 0;
}
*/
void measure_line(fits *fit, point start, point finish) {
	int deg = -1;
	point delta = { finish.x - start.x, finish.y - start.y };
	double pixdist = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (pixdist == 0.) return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	gboolean unit_is_as = (fit->focal_length > 0.0) && (fit->pixel_size_x > 0.0) && (fit->pixel_size_y == fit->pixel_size_x);
	if (unit_is_as) {
		double bin_X = com.pref.binning_update ? (double) fit->binning_x : 1.0;
		double conversionfactor = (((3600.0 * 180.0) / M_PI) / 1.0E3 * (double) fit->pixel_size_x / fit->focal_length) * bin_X;
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

static gboolean spectroscopy_selections_are_valid(cut_struct *arg) {
	gboolean a = (arg->wavenumber1 != arg->wavenumber2) && (arg->wavenumber1 > 0.0) && (arg->wavenumber2 > 0.0);
	gboolean b = (arg->cut_wn1.x >= 0) && (arg->cut_wn1.y >= 0) && (arg->cut_wn2.x >= 0) && (arg->cut_wn2.y >= 0) && (arg->cut_wn1.x < arg->fit->rx) && (arg->cut_wn1.y < arg->fit->ry) && (arg->cut_wn2.x < arg->fit->rx) && (arg->cut_wn2.y < arg->fit->ry);
	gboolean c = (!((arg->cut_wn1.x == arg->cut_wn2.x) && (arg->cut_wn1.y == arg->cut_wn2.y)));
	return a && b && c;
}

static double interpf(fits* fit, double x, double y, int chan) {
	if (chan >= fit->naxes[2])
		return -DBL_MAX;
	int w = fit->rx;
	int h = fit->ry;
	int npixels = w * h;
	int x0 = (int)(x);
	int x1 = (int)(x) + 1;
	int y0 = (int)(y);
	int y1 = (int)(y) + 1;
	if (x0 < 0 || x1 > w-1 || y0 < 0 || y1 > h-1)
		return NAN;
	float val00 = fit->fdata[x0 + y0 * w + npixels * chan];
	float val01 = fit->fdata[x1 + y0 * w + npixels * chan];
	float val10 = fit->fdata[x0 + y1 * w + npixels * chan];
	float val11 = fit->fdata[x1 + y1 * w + npixels * chan];
	float interp1 = (x - x0) * val00 + (x1 - x) * val01;
	float interp2 = (x - x0) * val10 + (x1 - x) * val11;
	float interp = (y - y0) * interp1 + (y1 - y) * interp2;
	return (double) interp;
}

static double interpw(fits* fit, double x, double y, int chan) {
	if (chan >= fit->naxes[2])
		return -DBL_MAX;
	int w = fit->rx;
	int h = fit->ry;
	int npixels = w * h;
	int x0 = (int)(x);
	int x1 = (int)(x) + 1;
	int y0 = (int)(y);
	int y1 = (int)(y) + 1;
	if (x0 < 0 || x1 > w-1 || y0 < 0 || y1 > h-1)
		return NAN;
	float val00 = (float) fit->data[x0 + y0 * w + npixels * chan];
	float val01 = (float) fit->data[x1 + y0 * w + npixels * chan];
	float val10 = (float) fit->data[x0 + y1 * w + npixels * chan];
	float val11 = (float) fit->data[x1 + y1 * w + npixels * chan];
	float interp1 = (x - x0) * val00 + (x1 - x) * val01;
	float interp2 = (x - x0) * val10 + (x1 - x) * val11;
	float interp = (y - y0) * interp1 + (y1 - y) * interp2;
	return (double) interp;
}

static double interp(fits *fit, double x, double y, int chan, int num, double dx, double dy) {
	double val = 0.0;
	int hw = (num - 1) / 2;
	for (int i = -hw ; i < hw + 1 ; i++) {
		switch (fit->type) {
			case DATA_FLOAT:
				val += interpf(fit, x + (i * dy), y + (i * dx), chan);
				break;
			case DATA_USHORT:

				val = interpw(fit, x + (i * dy), y + (i * dx), chan);
				break;
			default:
				return NAN;
				break;
		}
	}
	val /= num;
	return val;
}

static double nointerpf(fits *fit, int x, int y, int chan) {
	int w = fit->rx;
	int h = fit->ry;
	if (x < 0 || x > w-1 || y < 0 || y > h-1)
		return NAN;
	double val = (double) fit->fdata[x + y * w + w * h * chan];
	return val;
}

static double nointerpw(fits *fit, int x, int y, int chan) {
	int w = fit->rx;
	int h = fit->ry;
	if (x < 0 || x > w-1 || y < 0 || y > h-1)
		return NAN;
	double val = (double) fit->data[x + y * w + w * h * chan];
	return val;
}

static double nointerp(fits *fit, int x, int y, int chan, int num, int dx, int dy) {
	double val = 0.0;
	int hw = (num - 1) / 2;
	for (int i = -hw ; i < hw + 1 ; i++) {
		switch (fit->type) {
			case DATA_FLOAT:
				// Note the mismatch of x and dy etc is intentional, we are integrating
				// perpendicular to the cut line
				val += nointerpf(fit, x + (i * dy), y + (i * dx), chan);
				break;
			case DATA_USHORT:
				val += nointerpw(fit, x + (i * dy), y + (i * dx), chan);
				break;
			default:
				return NAN;
				break;
		}
	}
	val /= num;
	return val;
}

static void calc_zero_and_spacing(cut_struct *arg, double *zero, double *spectro_spacing) {
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

static void build_profile_filenames(cut_struct *arg, gchar **filename, gchar **imagefilename, gboolean *hastmpfile) {
	gchar *temp = NULL, *tempfilename = NULL;
	if (arg->filename) {
		*filename = g_strdup(arg->filename);
	} else {
		if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
			temp = g_path_get_basename(com.uniq->filename);
			*filename = g_strdup_printf("profile_%s_%s.dat", temp, build_timestamp_filename());
			g_free(temp);
			temp = NULL;
		} else if (arg->seq) {
			char* seq_image_canonical_name = calloc(256, 1);
			seq_get_image_filename(arg->seq, arg->imgnumber, seq_image_canonical_name);
			*filename = g_strdup_printf("profile_%s.dat", seq_image_canonical_name);
			free(seq_image_canonical_name);
		} else if (sequence_is_loaded()) { // when sequence is loaded in the GUI and we only apply cut to the currently loaded image
			char* seq_image_canonical_name = calloc(256, 1);
			seq_get_image_filename(&com.seq, com.seq.current, seq_image_canonical_name);
			*filename = g_strdup_printf("profile_%s.dat", seq_image_canonical_name);
			free(seq_image_canonical_name);
		} else {
			*filename = g_strdup_printf("profile_%s", build_timestamp_filename());
			siril_debug_print("%s\n", arg->filename);
		}
		tempfilename = profile_tmpfile();
		if (!(arg->save_dat || arg->seq))
			*hastmpfile = TRUE;
	}
	printf("Filename: %s\n", *filename);
	temp = g_path_get_basename(*filename);
	*imagefilename = replace_ext(temp, ".png");
	g_free(temp);
	temp = NULL;
	if (*hastmpfile) {
		g_free(*filename);
		*filename = g_strdup(tempfilename);
		free(tempfilename);
	}
}

gpointer cut_profile(gpointer p) {
	cut_struct* arg = (cut_struct*) p;
	gchar* legend = NULL;
	int retval = 0;
	gnuplot_ctrl *gplot = NULL;
	gboolean tmpfile = FALSE;
	gchar *filename = NULL, *imagefilename = NULL;
	double starty = arg->fit->ry - 1 - arg->cut_start.y;
	double endy = arg->fit->ry - 1 - arg->cut_end.y;

	build_profile_filenames(arg, &filename, &imagefilename, &tmpfile);
	gboolean use_gnuplot = gnuplot_is_available();
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), filename);
	} else {
		gplot = gnuplot_init();
	}

	point delta;
	delta.x = arg->cut_end.x - arg->cut_start.x;
	delta.y = endy - starty;
	double *x = NULL, *r = NULL, *g = NULL, *b = NULL;
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	gboolean hv = ((point_spacing_x == 1.) || (point_spacing_y == 1.) || (point_spacing_x == -1.) || (point_spacing_y == -1.));

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
			if (abs(point_spacing_x == 1.) || abs(point_spacing_y == 1.)) { // Horizontal, no interpolation
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
	if (arg->mode == CUT_MONO) {
		if (arg->fit->naxes[2] == 3) {
			// Now the vport is defined, set the data to a channel or even luminance
			if (arg->vport == 0) {
				printf("vport R\n");
				legend = g_strdup("x R");
			} else if (arg->vport == 1) {
				for (int i = 0 ; i < nbr_points ; i++)
					r[i] = g[i];
				printf("vport G\n");
				legend = g_strdup("x G");
			} else if (arg->vport == 2) {
				for (int i = 0 ; i < nbr_points ; i++)
					r[i] = b[i];
				printf("vport B\n");
				legend = g_strdup("x B");
			} else {
				for (int i = 0 ; i < nbr_points ; i++) {
					r[i] = (r[i] + g[i] + b[i]) / 3.0;
				}
				printf("vport L\n");
				legend = g_strdup("x L");
			}
		} else if (arg->fit->naxes[2] == 1) {
			legend = g_strdup("x Mono");
		}
		retval = gnuplot_write_xy_dat(filename, x, r, nbr_points, legend);
	} else {
		retval = gnuplot_write_xrgb_dat(filename, x, r, g, b, nbr_points, "x R G B");
	}
	printf("legend: %s\n", legend);
	if (retval) {
		siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", filename);
		goto END;
	} else {
		if (!(tmpfile || ((arg->seq) && !arg->save_dat)))
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
				if (arg->fit->naxes[2] == 1 || (arg->fit->naxes[2] == 3 && arg->mode == CUT_MONO)) {
					gnuplot_plot_xy_from_datfile(gplot, filename);
					if (arg->save_png_too) {
						gnuplot_plot_xy_datfile_colheader_to_png(gplot, filename, legend, imagefilename);
						siril_log_message(_("%s has been saved.\n"), imagefilename);
					}
				} else {
					gnuplot_plot_xrgb_from_datfile(gplot, filename);
					if (arg->save_png_too) {
						gnuplot_plot_xrgb_datfile_to_png(gplot, filename, imagefilename);
						siril_log_message(_("%s has been saved.\n"), imagefilename);
					}
				}
			} else {
				if (arg->fit->naxes[2] == 1) {
					gnuplot_plot_xy_datfile_colheader_to_png(gplot, filename, legend, imagefilename);
				} else {
					gnuplot_plot_xrgb_datfile_to_png(gplot, filename, imagefilename);
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
	if (tmpfile || (arg->seq != NULL && !arg->save_dat)) {
		siril_debug_print("Unlinking temporary file %s\n", filename);
		if (g_unlink(filename))
			siril_debug_print("Error in g_unlink()\n");
	}
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(filename);
	g_free(imagefilename);
	g_free(legend);
	free(x);
	free(r);
	free(g);
	free(b);
	arg->vport = -1;
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
	char *filename = NULL, *imagefilename = NULL;
	double starty = arg->fit->ry - 1 - arg->cut_start.y;
	double endy = arg->fit->ry - 1 - arg->cut_end.y;
	gboolean use_gnuplot = gnuplot_is_available();
	gnuplot_ctrl *gplot = NULL;

	build_profile_filenames(arg, &filename, &imagefilename, &tmpfile);
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), filename);
	} else {
		gplot = gnuplot_init();
	}

	point delta;
	delta.x = arg->cut_end.x - arg->cut_start.x;
	delta.y = endy - starty;
	double *x = NULL, *r[3] = { 0 };
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	gboolean hv = ((point_spacing_x == 1.) || (point_spacing_y == 1.) || (point_spacing_x == -1.) || (point_spacing_y == -1.));
	for (int i = 0 ; i < 3 ; i++)
		r[i] = malloc(nbr_points * sizeof(double));
	x = malloc(nbr_points * sizeof(double));
	for (int offset = -1 ; offset < 2 ; offset++) {
		double offstartx = arg->cut_start.x - (offset * point_spacing_y * arg->step);
		double offstarty = starty + (offset * point_spacing_x * arg->step);
		gboolean single_channel = (arg->vport == 0 || arg->vport == 1 || arg->vport == 2);
		gboolean redvport = arg->vport < 3 ? arg->vport : 0;
		for (int i = 0 ; i < nbr_points ; i++) {
			x[i] = i * point_spacing;
			if (hv) {
				// Horizontal / vertical, no interpolation
				r[offset+1][i] = nointerp(arg->fit, offstartx + point_spacing_x * i, offstarty + point_spacing_y * i, redvport, arg->width, (int) point_spacing_x, (int) point_spacing_y);
			} else {
				// Neither horizontal nor vertical: interpolate
				r[offset+1][i] = interp(arg->fit, (double) (offstartx + point_spacing_x * i),
							(double) (offstarty + point_spacing_y * i), redvport,
							arg->width, point_spacing_x, point_spacing_y);
			}
			if (arg->fit->naxes[2] > 1 && (!single_channel)) {
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
	gchar *titletext = NULL;
	if (arg->vport == 0) {
		if (arg->fit->naxes[2] == 3)
			titletext = g_strdup_printf("x R(-%dpx) R R(+%dpx)", (int) arg->step, (int) arg->step);
		else
			titletext = g_strdup_printf("x Mono(-%dpx) Mono Mono(+%dpx)", (int) arg->step, (int) arg->step);
	}
	else if (arg->vport == 1)
		titletext = g_strdup_printf("x G(-%dpx) G G(+%dpx)", (int) arg->step, (int) arg->step);
	else if (arg->vport == 2)
		titletext = g_strdup_printf("x B(-%dpx) B B(+%dpx)", (int) arg->step, (int) arg->step);
	else
		titletext = g_strdup_printf("x L(-%dpx) L L(+%dpx)", (int) arg->step, (int) arg->step);

	retval = gnuplot_write_xrgb_dat(filename, x, r[0], r[1], r[2], nbr_points, titletext);
	g_free(titletext);
	if (retval) {
		siril_log_color_message(_("Failed to create the cut data file %s\n"), "red", arg->filename);
		goto END;
	} else {
		if (!(tmpfile || ((arg->seq) && !arg->save_dat)))
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
				if (arg->save_png_too) {
					gnuplot_plot_xrgb_datfile_to_png(gplot, filename, imagefilename);
					siril_log_message(_("%s has been saved.\n"), imagefilename);
				}
			} else {
				gnuplot_plot_xrgb_datfile_to_png(gplot, filename, imagefilename);
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
	if (tmpfile || (arg->seq != NULL && !arg->save_dat)) {
		siril_debug_print("Unlinking temporary file %s\n", filename);
		if (g_unlink(filename))
			siril_debug_print("Error in g_unlink()\n");
	}
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(filename);
	g_free(imagefilename);
	free(x);
	for (int i = 0 ; i < 3 ; i++) {
		free(r[i]);
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
	char *filename = NULL,*imagefilename = NULL;
	fits cfa[4];
	memset(cfa, 0, 4 * sizeof(fits));
	gboolean use_gnuplot = gnuplot_is_available();
	gnuplot_ctrl *gplot = NULL;

	build_profile_filenames(arg, &filename, &imagefilename, &tmpfile);
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the brightness profile data will be produced in %s but no image will be created.\n"), filename);
	} else {
		gplot = gnuplot_init();
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

	// Coordinates of profile start and endpoints in CFA space
	point cut_start_cfa = { (int) (arg->cut_start.x / 2) , (int) (arg->cut_start.y / 2) };
	point cut_end_cfa = { (int) (arg->cut_end.x / 2) , (int) (arg->cut_end.y / 2) };
	double starty = (arg->fit->ry / 2) - 1 - cut_start_cfa.y;
	double endy = (arg->fit->ry / 2) - 1 - cut_end_cfa.y;
	point delta = { cut_end_cfa.x - cut_start_cfa.x, endy - starty };
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	gboolean hv = ((point_spacing_x == 1.) || (point_spacing_y == 1.) || (point_spacing_x == -1.) || (point_spacing_y == -1.));
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
		if (!(tmpfile || ((arg->seq) && !arg->save_dat)))
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
				if (arg->save_png_too) {
					gnuplot_plot_xcfa_datfile_to_png(gplot, filename, imagefilename);
					siril_log_message(_("%s has been saved.\n"), imagefilename);
				}
			} else {
				gnuplot_plot_xcfa_datfile_to_png(gplot, filename, imagefilename);
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
	if (tmpfile || (arg->seq != NULL && !arg->save_dat)) {
		siril_debug_print("Unlinking temporary file %s\n", filename);
		if (g_unlink(filename))
			siril_debug_print("Error in g_unlink()\n");
	}
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(filename);
	g_free(imagefilename);
	free(x);
	for (int i = 0 ; i < 4 ; i++) {
		free(r[i]);
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
void apply_cut_to_sequence(cut_struct* cut_args);

void on_cut_sequence_apply_from_gui() {
	GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
	cut_struct *arg = calloc(1, sizeof(cut_struct));
	memcpy(arg, &gui.cut, sizeof(cut_struct));
	arg->save_png_too = FALSE;
	arg->fit = NULL;
	arg->seq = &com.seq;
	if (gtk_toggle_button_get_active(cut_color))
		arg->mode = CUT_COLOR;
	else
		arg->mode = CUT_MONO;
	arg->display_graph = FALSE;
	arg->cut_measure = FALSE;
	control_window_switch_to_tab(OUTPUT_LOGS);
	// Check args are cromulent
	if (cut_struct_is_valid(arg))
		apply_cut_to_sequence(arg);
}

void on_cut_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton* apply_to_sequence = (GtkToggleButton*)lookup_widget("cut_apply_to_sequence");
	if (gtk_toggle_button_get_active(apply_to_sequence)) {
		if (sequence_is_loaded())
			on_cut_sequence_apply_from_gui();
		else
			siril_message_dialog(GTK_MESSAGE_ERROR,
					_("No sequence is loaded"),
					_("The Apply to sequence option is checked, but no sequence is loaded."));

	} else {
		GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
		gui.cut.fit = &gfit;
		gui.cut.seq = NULL;
		gui.cut.vport = gui.cvport;
		gui.cut.display_graph = TRUE;
		if (gui.cut.tri) {
			siril_debug_print("Tri-profile\n");
			start_in_new_thread(tri_cut, &gui.cut);
		} else if (gui.cut.cfa) {
			siril_debug_print("CFA profiling\n");
			start_in_new_thread(cfa_cut, &gui.cut);
		} else {
			if (gtk_toggle_button_get_active(cut_color)) {
				siril_debug_print("Color profiling\n");
				gui.cut.mode = CUT_COLOR;
			} else {
				siril_debug_print("Mono profiling\n");
				gui.cut.mode = CUT_MONO;
			}
			start_in_new_thread(cut_profile, &gui.cut);
		}
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
	GtkWidget* cfabutton = lookup_widget("cut_cfa");
	gtk_widget_set_sensitive(colorbutton, (gfit.naxes[2] == 3));
	sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit.bayer_pattern);
	gboolean cfa_disabled = ((gfit.naxes[2] > 1) || ((!(pattern == BAYER_FILTER_RGGB || pattern == BAYER_FILTER_GRBG || pattern == BAYER_FILTER_BGGR || pattern == BAYER_FILTER_GBRG))));
		gtk_widget_set_sensitive(cfabutton, !cfa_disabled);
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

void on_cut_save_png_toggled(GtkToggleButton *button, gpointer user_data) {
	gui.cut.save_png_too = gtk_toggle_button_get_active(button);
}

void on_cut_apply_to_sequence_toggled(GtkToggleButton *button, gpointer user_data) {
	// The sequence mode will always save PNGs, so we set the option at the same time.
	// It's not really necessary as the structure member isn't used, but it keeps
	// things consistent for the user.
	GtkToggleButton *pngbutton = (GtkToggleButton*) lookup_widget("cut_save_png");
	GtkWidget *measurebutton = (GtkWidget*) lookup_widget("cut_measure_profile");
	if (gtk_toggle_button_get_active(button)) {
		gtk_toggle_button_set_active(pngbutton, TRUE);
		gui.cut.save_png_too = TRUE;
		gtk_toggle_button_set_active((GtkToggleButton*) measurebutton, FALSE);
		gtk_widget_set_sensitive(measurebutton, FALSE);
		gui.cut.cut_measure = FALSE;
	} else {
		gtk_widget_set_sensitive(measurebutton, TRUE);
	}
}

void on_cut_coords_measure_button_clicked(GtkButton *button, gpointer user_data) {
	measure_line(&gfit, gui.cut.cut_start, gui.cut.cut_end);
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

static gboolean cut_idle_function(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	cut_struct *data = (cut_struct *) args->user;
	free(data);
	gboolean retval = end_generic_sequence(args);
	return retval;
}

static int cut_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	cut_struct *cut_args = (cut_struct*) args->user;
	int ret = 0;
	cut_struct *private_args = malloc(sizeof(cut_struct));
	memcpy(private_args, cut_args, sizeof(cut_struct));
	private_args->fit = fit;

	// Check points are not out of bounds. The >=0 criterion is checked earlier
	// but we have to check <fit.rx and <fit.ry for every image in case of
	// sequences with variable image sizes. (It probably doesn't make sense to
	// profile such a sequence but at least we should ensure we don't segfault.)

	if (private_args->cut_wn1.x > fit->rx - 1 || private_args->cut_wn2.x > fit->rx - 1 || private_args->cut_wn1.y > fit->ry - 1 || private_args->cut_wn2.y > fit->ry - 1) {
		siril_log_message(_("Error: wavenumber point outside image dimensions.\n"));
		return 1;
	}
	if (private_args->cut_start.x > fit->rx - 1 || private_args->cut_start.y > fit->ry - 1 || private_args->cut_end.x > fit->rx - 1 || private_args->cut_end.y > fit->ry - 1) {
		siril_log_message(_("Error: profile line endpoint outside image dimensions.\n"));
		return 1;
	}
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

