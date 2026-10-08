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

/*
 * Non-GUI lifecycle helpers for cut_struct (intensity profile / spectrogram).
 * The GTK dialog callbacks remain in gui/cut.c.
 */

#include <math.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/cut.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "core/siril_date.h"
#include "algos/extraction.h"
#include "algos/fitting.h"
#include "algos/demosaicing.h"
#include "algos/PSF.h"
#include "io/siril_plot.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"

void initialize_cut_struct(cut_struct *arg) {
	arg->fit = gfit;
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
	arg->plot_as_wavenumber = FALSE;
	arg->wavenumber1 = -1;
	arg->wavenumber2 = -1;
	arg->plot_spectro_bg = FALSE;
	arg->bg_poly_order = 3;
	arg->tri = FALSE;
	arg->mode = CUT_MONO;
	arg->width = 1;
	arg->step = 1;
	arg->display_graph = TRUE;
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(arg->title);
	arg->title = NULL;
	g_free(arg->user_title);
	arg->user_title = NULL;
	arg->title_has_sequence_numbers = FALSE;
	arg->save_dat = FALSE;
	arg->save_png_too = FALSE;
	arg->pref_as = gfit->keywords.wcsdata.pltsolvd;
	arg->vport = -1;
	/* gui_function(reset_cut_gui, NULL) is called by gui/cut.c when needed */
}

void free_cut_args(cut_struct *arg) {
	if (!arg) return;
	g_free(arg->filename);   arg->filename   = NULL;
	g_free(arg->title);      arg->title      = NULL;
	g_free(arg->user_title); arg->user_title = NULL;
	free(arg);
}

gboolean cut_struct_is_valid(cut_struct *arg) {
	int nb_layers = (arg->seq) ? arg->seq->nb_layers : arg->fit->naxes[2];
	if (arg->tri && arg->cfa) {
		siril_log_error(_("Error: CFA mode and tri-profile are mutually exclusive.\n"));
		return FALSE;
	}
	if (arg->tri && arg->mode == CUT_COLOR) {
		siril_log_error(_("Error: color plot and tri-profile are mutually exclusive.\n"));
		return FALSE;
	}
	if (nb_layers == 1 && arg->mode == CUT_COLOR) {
		siril_log_error(_("Error: mono image is loaded: color mode is unavailable.\n"));
		return FALSE;
	}
	if (nb_layers == 1 && arg->vport > 0) {
		siril_log_error(_("Error: mono image is loaded: colored layers cannot be specified.\n"));
		return FALSE;
	}
	if (arg->cfa && arg->mode == CUT_COLOR) {
		siril_log_error(_("Error: color plot and CFA mode are mutually exclusive.\n"));
		return FALSE;
	}
	if (arg->cut_start.x == arg->cut_end.x && arg->cut_start.y == arg->cut_end.y) {
		siril_log_error(_("Error: start and finish points are the same.\n"));
		return FALSE;
	}
	if (arg->seq && arg->seq->is_variable) {
		siril_log_error(_("Error: variable image size in sequence.\n"));
		return FALSE;
	}

	int rx = (arg->seq) ? arg->seq->rx : gfit->rx;
	int ry = (arg->seq) ? arg->seq->ry : gfit->ry;

	if (arg->cut_wn1.x > -1.0 && arg->cut_wn1.x == arg->cut_wn2.x && arg->cut_wn1.y == arg->cut_wn2.y) {
		siril_log_error(_("Error: wavenumber points are the same.\n"));
		return FALSE;
	}
	if (arg->cut_wn1.x > -1.0 && arg->wavenumber1 == -1.0) {
		siril_log_error(_("Error: wavenumber / wavelength for point 1 is not set.\n"));
		return FALSE;
	}
	if (arg->cut_wn2.x > -1.0 && arg->wavenumber2 == -1.0) {
		siril_log_error(_("Error: wavenumber / wavelength for point 2 is not set.\n"));
		return FALSE;
	}
	if (arg->wavenumber1 >= 0.0 && (arg->cut_wn1.x < 0.0 || arg->cut_wn1.y < 0.0)) {
		siril_log_error(_("Error: wavenumber / wavelength set for point 1 but corresponding location is not set.\n"));
		return FALSE;
	}
	if (arg->wavenumber2 >= 0.0 && (arg->cut_wn2.x < 0.0 || arg->cut_wn2.y < 0.0)) {
		siril_log_error(_("Error: wavenumber / wavelength set for point 2 but corresponding location is not set.\n"));
		return FALSE;
	}
	if ((arg->cut_wn1.x < 0.0 || arg->cut_wn2.y < 0.0 || arg->cut_wn2.x < 0.0 || arg->cut_wn2.y < 0.0) &&
		(!(arg->cut_wn1.x == -1.0 && arg->cut_wn1.y == -1.0 && arg->cut_wn2.x == -1.0 && arg->cut_wn2.y == -1.0))) {
		siril_log_error(_("Error: wavenumber / wavelength point outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->cut_wn1.x > rx - 1 || arg->cut_wn2.x > rx - 1 || arg->cut_wn1.y > ry - 1 || arg->cut_wn2.y > ry - 1) {
		siril_log_error(_("Error: wavenumber / wavelength point outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->cut_start.x > rx - 1 || arg->cut_start.y > ry - 1 || arg->cut_end.x > rx - 1 || arg->cut_end.y > ry - 1) {
		siril_log_error(_("Error: profile line endpoint outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->step != 1 && !arg->tri) {
		siril_log_error(_("Error: spacing is set but tri-profile mode is not selected.\n"));
		return FALSE;
	}
	if (arg->vport < 0 || arg->vport > nb_layers) {
		siril_log_error(_("Error: layer out of range.\n"));
		return FALSE;
	}
	if (arg->bg_poly_order < 1 || arg->bg_poly_order > 6) {
		siril_log_error(_("Error: background removal polynomial degree should be >= 1 and <= 6.\n"));
		return FALSE;
	}
	return TRUE;
}

gchar* cut_make_title(cut_struct *arg, gboolean spectro) {
	gchar* str = NULL;

	if (arg->user_title) {
		// Check for trailing brackets
		if (g_str_has_suffix(arg->user_title, "()")) {
			arg->title_has_sequence_numbers = TRUE;
			gchar* temp = g_malloc(strlen(arg->user_title) - 1);
			snprintf(temp, strlen(arg->user_title) - 1, "%s", arg->user_title);
			g_free(arg->title);
			arg->title = g_strdup(temp);
			g_free(temp);
		} else {
			arg->title_has_sequence_numbers = FALSE;
			arg->title = g_strdup(arg->user_title);
		}

		if (arg->title_has_sequence_numbers && arg->seq) {
			str = g_strdup_printf("%s(%d / %d)", arg->title, arg->imgnumber + 1, arg->seq->selnum);
		} else {
			str = g_strdup(arg->title);
		}
	} else {
		str = g_strdup(spectro ? _("Spectrogram") : _("Intensity Profile"));
	}
	return str;
}


double get_conversion_factor(fits *fit) {
	gboolean unit_is_as = (fit->keywords.focal_length > 0.0) && (fit->keywords.pixel_size_x > 0.0) && (fit->keywords.pixel_size_y == fit->keywords.pixel_size_x);
	double conversionfactor = -DBL_MAX;
	if (unit_is_as) {
		double bin_X = com.pref.binning_update ? (double) fit->keywords.binning_x : 1.0;
		conversionfactor = (((3600.0 * 180.0) / G_PI) / 1.0E3 * (double) fit->keywords.pixel_size_x / fit->keywords.focal_length) * bin_X;
	}
	return conversionfactor;
}

gboolean spectroscopy_selections_are_valid(cut_struct *arg) {
	gboolean a = (arg->wavenumber1 != arg->wavenumber2) && (arg->wavenumber1 > 0.0) && (arg->wavenumber2 > 0.0);
	gboolean b = (arg->cut_wn1.x >= 0) && (arg->cut_wn1.y >= 0) && (arg->cut_wn2.x >= 0) && (arg->cut_wn2.y >= 0) && (arg->cut_wn1.x < arg->fit->rx) && (arg->cut_wn1.y < arg->fit->ry) && (arg->cut_wn2.x < arg->fit->rx) && (arg->cut_wn2.y < arg->fit->ry);
	gboolean c = (!((arg->cut_wn1.x == arg->cut_wn2.x) && (arg->cut_wn1.y == arg->cut_wn2.y)));
	return a && b && c;
}

static double interp(fits *fit, double x, double y, int chan, int num, double dx, double dy) {
	double val = 0.0;
	int hw = (num - 1) / 2;
	for (int i = -hw ; i < hw + 1 ; i++) {
		switch (fit->type) {
			case DATA_FLOAT:
				val += bilinear(fit->fpdata[chan], fit->rx, fit->ry, x + ((double) i * dy), y + ((double) i * dx));
				break;
			case DATA_USHORT:
				val += bilinear_ushort(fit->pdata[chan], fit->rx, fit->ry, x + ((double) i * dy), y + ((double) i * dx));
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
	double wavelength1 = 10000000. / arg->wavenumber1;
	double wavelength2 = 10000000. / arg->wavenumber2;

	double wndiff = wavelength2 - wavelength1;
	*spectro_spacing = wndiff / wndiff_dist; // Spacing is in wavelength

	// To calculate the zero we will work from whichever of x or y has the biggest difference
	double z2_z1 = fabs(wndelta.y) > fabs(wndelta.x) ? wndelta.y : wndelta.x;
	double z1_z0 = fabs(wndelta.y) > fabs(wndelta.x) ? arg->cut_wn1.y - arg->cut_start.y : arg->cut_wn1.x - arg->cut_start.x;
	double z2_z0 = fabs(wndelta.y) > fabs(wndelta.x) ? arg->cut_wn2.y - arg->cut_start.y : arg->cut_wn2.x - arg->cut_start.x;
	*zero = wavelength1 - ( ( z1_z0 * wndiff ) / z2_z1 ); // Zero is in wavelength

	// Check whether the spacing needs to be positive or negative
	if ((fabs(z2_z0) < fabs(z1_z0) && wavelength1 < wavelength2))// || (fabs(z2_z0) > fabs(z1_z0) && wavelength1 > wavelength2))
		*spectro_spacing *= -1.0;
	return;
}

static void build_profile_filenames(cut_struct *arg, gchar **filename, gchar **imagefilename) {
	gchar *temp = NULL, *timestamp = NULL;
	if (arg->filename) {
		*filename = g_strdup(arg->filename);
	} else {
		if (arg->seq) {
			char seq_image_canonical_name[256] = "";
			seq_get_image_filename(arg->seq, arg->imgnumber, seq_image_canonical_name);
			*filename = g_strdup_printf("profile_%s.dat", seq_image_canonical_name);
		} else if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
			temp = g_path_get_basename(com.uniq->filename);
			timestamp = build_timestamp_filename();
			*filename = g_strdup_printf("profile_%s_%s.dat", temp, timestamp);
			g_free(temp);
		} else if (sequence_is_loaded()) { // when sequence is loaded in the GUI and we only apply cut to the currently loaded image
			char seq_image_canonical_name[256] = "";
			seq_get_image_filename(&com.seq, com.seq.current, seq_image_canonical_name);
			*filename = g_strdup_printf("profile_%s.dat", seq_image_canonical_name);
		} else {
			timestamp = build_timestamp_filename();
			*filename = g_strdup_printf("profile_%s", timestamp);
			siril_log_debug("%s\n", arg->filename);
		}
	}
	temp = g_path_get_basename(*filename);
	*imagefilename = replace_ext(temp, ".png");
	g_free(temp);
	g_free(timestamp);
	temp = NULL;
}

gpointer cut_profile(gpointer p) {
	cut_struct* arg = (cut_struct*) p;
	gchar *spl_legend = NULL;
	int retval = 0;
	gchar *filename = NULL, *imagefilename = NULL;
	double starty = arg->fit->ry - 1 - arg->cut_start.y;
	double endy = arg->fit->ry - 1 - arg->cut_end.y;
	double *x = NULL, *r = NULL, *g = NULL, *b = NULL;
	siril_plot_data *spl_data = init_siril_plot_data();
	if (!spl_data) {
		retval = 1;
		goto END;
	}

	build_profile_filenames(arg, &filename, &imagefilename);

	point delta;
	delta.x = arg->cut_end.x - arg->cut_start.x;
	delta.y = endy - starty;
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	gboolean xscale = spectroscopy_selections_are_valid(arg);
	// If pref_as is set andmetadata allows, represent spacing as arcsec
	double conversionfactor = get_conversion_factor(arg->fit);
	if (arg->pref_as && !xscale) {
		if (conversionfactor != -DBL_MAX) {
			point_spacing *= conversionfactor;
		}
	}
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	gboolean hv = ((point_spacing_x == 1.) || (point_spacing_y == 1.) || (point_spacing_x == -1.) || (point_spacing_y == -1.));

	r = malloc(nbr_points * sizeof(double));
	if (arg->fit->naxes[2] > 1) {
		g = malloc(nbr_points * sizeof(double));
		b = malloc(nbr_points * sizeof(double));
	}
	x = malloc(nbr_points * sizeof(double));
	double zero = 0.0, spectro_spacing = 1.0;
	if (xscale)
		calc_zero_and_spacing(arg, &zero, &spectro_spacing);
	for (int i = 0 ; i < nbr_points ; i++) {
		if (xscale) {
			x[i] = arg->plot_as_wavenumber ? 10000000. / (zero + i * spectro_spacing) : zero + i * spectro_spacing;
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
			if (fabs(point_spacing_x) == 1. || fabs(point_spacing_y) == 1.) { // Horizontal, no interpolation
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
				spl_legend = g_strdup("R");
			} else if (arg->vport == 1) {
				for (int i = 0 ; i < nbr_points ; i++)
					r[i] = g[i];
				spl_legend = g_strdup("G");
			} else if (arg->vport == 2) {
				for (int i = 0 ; i < nbr_points ; i++)
					r[i] = b[i];
				spl_legend = g_strdup("B");
			} else {
				for (int i = 0 ; i < nbr_points ; i++) {
					r[i] = (r[i] + g[i] + b[i]) / 3.0;
				}
				spl_legend = g_strdup("L");
			}
		} else if (arg->fit->naxes[2] == 1) {
			spl_legend = g_strdup("Mono");
		}
	} else {
		spl_legend = g_strdup("R");
	}
	/* Plotting cut profile */
	gchar *xlabel = NULL, *title = NULL;
	title = cut_make_title(arg, xscale); // must be freed with g_free()
	if (xscale) {
		xlabel = arg->plot_as_wavenumber ? g_strdup_printf(_("Wavenumber / cm^{-1}")) : g_strdup_printf(_("Wavelength / nm"));
	} else {
		if (arg->pref_as && conversionfactor != -DBL_MAX) {
			xlabel = g_strdup_printf(_("Distance along cut / arcsec"));
		} else {
			xlabel = g_strdup_printf(_("Distance along cut / px"));
		}
	}

	siril_plot_set_title(spl_data, title);
	siril_plot_set_xlabel(spl_data, xlabel);
	siril_plot_add_xydata(spl_data, spl_legend, nbr_points, x, r, NULL, NULL);
	siril_plot_set_nth_color(spl_data, 1, (double[3]){1., 0., 0.});
	siril_plot_set_savename(spl_data, "profile");
	if (arg->fit->naxes[2] == 3 && arg->mode != CUT_MONO) {
		siril_plot_add_xydata(spl_data, "G", nbr_points, x, g, NULL, NULL);
		siril_plot_add_xydata(spl_data, "B", nbr_points, x, b, NULL, NULL);
		siril_plot_set_nth_color(spl_data, 2, (double[3]){0., 1., 0.});
		siril_plot_set_nth_color(spl_data, 3, (double[3]){0., 0., 1.});
	}
	if (arg->save_dat)
		siril_plot_save_dat(spl_data, filename, FALSE);
	if (arg->save_png_too || !arg->display_graph)
		siril_plot_save_png(spl_data, imagefilename, 0, 0);
	if (!arg->display_graph) { // if not used for display we can free spl_data now
		free_siril_plot_data(spl_data);
		spl_data = NULL; // just in case we try to use it later on
	}

	g_free(title);
	g_free(xlabel);

END:
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(filename);
	g_free(imagefilename);
	g_free(spl_legend);
	g_free(arg->title);
	arg->title = NULL;
	free(x);
	free(r);
	free(g);
	free(b);
	arg->vport = -1;
	gboolean in_sequence = (arg->seq != NULL);
	if (in_sequence || !arg->display_graph) {
		free_siril_plot_data(spl_data);
	} else if (!in_sequence) {
		if (arg->display_graph)
			gui_iface.show_siril_plot(spl_data);
		siril_add_idle(end_generic, NULL);
	}
	free(arg);
	return GINT_TO_POINTER(retval);
}

gpointer tri_cut(gpointer p) {
	cut_struct* arg = (cut_struct*) p;
	int retval = 0;
	char *filename = NULL, *imagefilename = NULL;
	double starty = arg->fit->ry - 1 - arg->cut_start.y;
	double endy = arg->fit->ry - 1 - arg->cut_end.y;
	double *x = NULL, *r[3] = { 0 };
	gchar *spllabels[3] = { NULL };
	siril_plot_data *spl_data = init_siril_plot_data();
	if (!spl_data) {
		retval = 1;
		goto END;
	}

	build_profile_filenames(arg, &filename, &imagefilename);

	point delta;
	delta.x = arg->cut_end.x - arg->cut_start.x;
	delta.y = endy - starty;
	double length = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.) {
		retval = 1;
		free_siril_plot_data(spl_data);
		spl_data = NULL;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	gboolean do_spec = spectroscopy_selections_are_valid(arg);
	// If pref_as is set andmetadata allows, represent spacing as arcsec
	double conversionfactor = get_conversion_factor(arg->fit);
	if (arg->pref_as && !do_spec) {
		if (conversionfactor != -DBL_MAX) {
			point_spacing *= conversionfactor;
		}
	}
	double point_spacing_x = (double) delta.x / nbr_points;
	double point_spacing_y = (double) delta.y / nbr_points;
	gboolean hv = ((point_spacing_x == 1.) || (point_spacing_y == 1.) || (point_spacing_x == -1.) || (point_spacing_y == -1.));
	for (int i = 0 ; i < 3 ; i++)
		r[i] = malloc(nbr_points * sizeof(double));
	x = malloc(nbr_points * sizeof(double));
	double zero = 0.0, spectro_spacing = 1.0;
	if (do_spec)
		calc_zero_and_spacing(arg, &zero, &spectro_spacing);
	for (int i = 0 ; i < nbr_points ; i++) {
		if (do_spec) {
			x[i] = arg->plot_as_wavenumber ? 10000000. / (zero + i * spectro_spacing) : zero + i * spectro_spacing;
		} else {
			x[i] = i * point_spacing;
		}
	}
	for (int offset = -1 ; offset < 2 ; offset++) {
		double offstartx = arg->cut_start.x - (offset * point_spacing_y * arg->step);
		double offstarty = starty + (offset * point_spacing_x * arg->step);
		gboolean single_channel = (arg->vport == 0 || arg->vport == 1 || arg->vport == 2);
		gboolean redvport = arg->vport < 3 ? arg->vport : 0;
		for (int i = 0 ; i < nbr_points ; i++) {
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
				if (fabs(point_spacing_x) == 1. || fabs(point_spacing_y) == 1.) { // Horizontal, no interpolation
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

	// If spectroscopic options are set, use r[0] and r[2] as dark strips,
	// to calibrate the spectrum. (This intentionally clobbers r[0] and r[2])
	if (do_spec) {
		for (int i = 0 ; i < nbr_points ; i++) {
			// Average the dark strips on either side
			r[2][i] = (r[0][i] + r[2][i]) / 2.0;
			r[0][i] = (double) i;
		}
		int degree = arg->bg_poly_order;
		double *coeffs = calloc(degree + 1, sizeof(double));
		double *uncertainties = calloc(degree + 1, sizeof(double));
		double sigma;
		robust_polynomial_fit(r[0], r[2], nbr_points, degree, coeffs, uncertainties, NULL, &sigma);
		GString *text = g_string_new(_("Coefficients: y = "));
		for (int i = 0 ; i <= degree ; i++) {
			gchar *tmp;
			if (i == 0) {
				tmp = g_strdup_printf("%.2e (±%.2e) ", coeffs[0], uncertainties[0]);
			} else {
				if (coeffs[i] >= 0.0)
					tmp = g_strdup_printf("+ %.2e (±%.2e) * x^%d ", coeffs[i], uncertainties[i], i);
				else
					tmp = g_strdup_printf("- %.2e (±%.2e) * x^%d ", -coeffs[i], uncertainties[i], i);
			}
			text = g_string_append(text, tmp);
			g_free(tmp);
		}

		gchar *coeffs_text = g_string_free(text, FALSE);
		siril_log_message(_("Subtracting dark strips: polynomial fit of degree %d:\n"), degree);
		siril_log_status("%s\nFit σ: %.2e\n", coeffs_text, sigma);
		g_free(coeffs_text);

		for (int i = 0 ; i < nbr_points ; i++) {
			r[0][i] = evaluate_polynomial(coeffs, degree, (double) i);
			r[1][i] -= r[0][i];
		}
	}

	if (arg->vport == 0) {
		if (do_spec) {
			spllabels[1] = g_strdup(_("Intensity"));
		} else if (arg->fit->naxes[2] == 3) {
			spllabels[1] = g_strdup(_("R"));
		} else {
			spllabels[1] = g_strdup(_("Mono"));
		}
	}
	else if (arg->vport == 1) {
		spllabels[1] = g_strdup(_("G"));
	} else if (arg->vport == 2) {
		spllabels[1] = g_strdup(_("B"));
	} else {
		spllabels[1] = do_spec ? g_strdup(_("Intensity")) : g_strdup(_("L"));
	}
	if (do_spec) {
		spllabels[0] = g_strdup_printf(_("Fitted background"));
		spllabels[2] = g_strdup_printf(_("Measured background"));
	} else {
		spllabels[0] = g_strdup_printf("%s(-%dpx)", spllabels[1], (int) arg->step);
		spllabels[2] = g_strdup_printf("%s(+%dpx)", spllabels[1], (int) arg->step);
	}

	/* Plotting cut profile */
	gchar *xlabel = NULL, *title = NULL;
	title = cut_make_title(arg, FALSE); // must be freed with g_free()
	if (do_spec) {
		xlabel = arg->plot_as_wavenumber ? g_strdup_printf(_("Wavenumber / cm^{-1}")) : g_strdup_printf(_("Wavelength / nm"));
	} else if (arg->pref_as && conversionfactor != -DBL_MAX) {
		xlabel = g_strdup_printf(_("Distance along cut / arcsec"));
	} else {
		xlabel = g_strdup_printf(_("Distance along cut / px"));
	}

	siril_plot_set_title(spl_data, title);
	siril_plot_set_xlabel(spl_data, xlabel);
	siril_plot_set_savename(spl_data, "profile");
	if (!do_spec || arg->plot_spectro_bg)
		siril_plot_add_xydata(spl_data, spllabels[0], nbr_points, x, r[0], NULL, NULL);
	siril_plot_add_xydata(spl_data, spllabels[1], nbr_points, x, r[1], NULL, NULL);
	if (!do_spec || arg->plot_spectro_bg)
		siril_plot_add_xydata(spl_data, spllabels[2], nbr_points, x, r[2], NULL, NULL);
	if (arg->save_dat)
		siril_plot_save_dat(spl_data, filename, FALSE);
	if (arg->save_png_too || !arg->display_graph)
		siril_plot_save_png(spl_data, imagefilename, 0, 0);
	if (!arg->display_graph) { // if not used for display we can free spl_data now
		free_siril_plot_data(spl_data);
		spl_data = NULL; // just in case we try to use it later on
	}

	g_free(title);
	g_free(xlabel);

END:
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(filename);
	g_free(imagefilename);
	free(x);
	g_free(arg->title);
	arg->title = NULL;
	for (int i = 0 ; i < 3 ; i++) {
		free(r[i]);
		g_free(spllabels[i]);
	}
	gboolean in_sequence = (arg->seq != NULL);
	if (!in_sequence) {
		if (arg->display_graph)
			gui_iface.show_siril_plot(spl_data);
		siril_add_idle(end_generic, NULL);
	}
	free(arg);
	return GINT_TO_POINTER(retval);
}

gpointer cfa_cut(gpointer p) {
	cut_struct* arg = (cut_struct*) p;
	int retval = 0, ret = 0;
	double *x = NULL, *r[4] = { 0 };
	char *filename = NULL,*imagefilename = NULL;
	fits cfa[4];
	memset(cfa, 0, 4 * sizeof(fits));
	siril_plot_data *spl_data = init_siril_plot_data();
	if (!spl_data) {
		retval = 1;
		goto END;
	}
	build_profile_filenames(arg, &filename, &imagefilename);
	sensor_pattern pattern = get_validated_cfa_pattern(arg->fit, FALSE, FALSE);
	if (pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX) {
			siril_log_error(_("Error: failed to read CFA pattern or invalid found.\n"));
			retval = 1;
			free_siril_plot_data(spl_data);
			spl_data = NULL;
			goto END;	
	}
	const gchar* pattern_str = filter_pattern[pattern];

	// Split arg->fit into 4 x Bayer sub-patterns cfa[0123]
	if (arg->fit->type == DATA_USHORT) {
		if ((ret = split_cfa_ushort(arg->fit, &cfa[0], &cfa[1], &cfa[2], &cfa[3]))) {
			siril_log_error(_("Error: failed to split FITS into CFA sub-patterns.\n"));
			retval = 1;
			free_siril_plot_data(spl_data);
			spl_data = NULL;
			goto END;
		}
	} else {
		if ((ret = split_cfa_float(arg->fit, &cfa[0], &cfa[1], &cfa[2], &cfa[3]))) {
			siril_log_error(_("Error: failed to split FITS into CFA sub-patterns.\n"));
			retval = 1;
			free_siril_plot_data(spl_data);
			spl_data = NULL;
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
		free_siril_plot_data(spl_data);
		spl_data = NULL;
		goto END;
	}
	int nbr_points = (int) length;
	double point_spacing = length / nbr_points;
	// If pref_as is set andmetadata allows, represent spacing as arcsec
	double conversionfactor = get_conversion_factor(arg->fit);
	if (arg->pref_as) {
		if (conversionfactor != -DBL_MAX) {
			point_spacing *= conversionfactor;
		}
	}
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
	/* Plotting cut profile */
	gchar *xlabel = NULL, *title = NULL;
	title = cut_make_title(arg, FALSE); // must be freed with g_free()
	if (arg->pref_as && conversionfactor != -DBL_MAX) {
		xlabel = g_strdup_printf(_("Distance along cut / arcsec"));
	} else {
		xlabel = g_strdup_printf(_("Distance along cut / px"));
	}

	siril_plot_set_title(spl_data, title);
	siril_plot_set_xlabel(spl_data, xlabel);
	siril_plot_set_savename(spl_data, "profile");
	gboolean first_green = TRUE;
	for (int i = 0; i < 4; i++) {
		double color[3] = { 0.0, 0.0, 0.0 };
		const char *label;
		switch (pattern_str[i]) {
			case 'R':
				label = "Red";
				color[0] = 1.0;
				break;
			case 'G':
			//RGB(0.078, 0.392, 0.078)
			//RGB(0.196, 1.0, 0.196)
				label = first_green ? "Green1" : "Green2";
				color[0] = first_green ? 0.078 : 0.196;
				color[1] = first_green ? 0.392 : 1.000;
				color[2] = first_green ? 0.078 : 0.196;
				first_green = FALSE;
				break;
			case 'B':
				label = "Blue";
				color[2] = 1.0;
				break;
			default:
				label = "Error";
				break;
		}
		siril_plot_add_xydata(spl_data, label, nbr_points, x, r[i], NULL, NULL);
		siril_plot_set_nth_color(spl_data, i + 1, color);
	}
	if (arg->save_dat)
		siril_plot_save_dat(spl_data, filename, FALSE);
	if (arg->save_png_too || !arg->display_graph)
		siril_plot_save_png(spl_data, imagefilename, 0, 0);
	if (!arg->display_graph) { // if not used for display we can free spl_data now
		free_siril_plot_data(spl_data);
		spl_data = NULL; // just in case we try to use it later on
	}

	g_free(title);
	g_free(xlabel);
END:
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(filename);
	g_free(imagefilename);
	g_free(arg->title);
	arg->title = NULL;
	free(x);
	for (int i = 0 ; i < 4 ; i++) {
		free(r[i]);
		clearfits(&cfa[i]);
	}
	gboolean in_sequence = (arg->seq != NULL);
	if (!in_sequence) {
		if (arg->display_graph)
			gui_iface.show_siril_plot(spl_data);
		siril_add_idle(end_generic, NULL);
	}
	free(arg);
	return GINT_TO_POINTER(retval);
}

static int cut_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	cut_struct* arg = (cut_struct*) args->user;
	unsigned int MB_per_input_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_input_image, NULL, &MB_avail);
	if (arg->cfa)
		required = 2 * MB_per_input_image;
	else
		required = MB_per_input_image;
	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		limit = thread_limit;
	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_error(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"), args->description, mem_per_thread, mem_available);
		g_free(mem_per_thread);
		g_free(mem_available);
	}
	return limit;
}

static gboolean cut_idle_function(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	cut_struct *data = (cut_struct *) args->user;
	g_free(data->user_title);
	data->user_title = NULL;
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
	private_args->imgnumber = i;
	if (private_args->cfa)
		ret = GPOINTER_TO_INT(cfa_cut(private_args));
	else if (private_args->tri)
		ret = GPOINTER_TO_INT(tri_cut(private_args));
	else
		ret = GPOINTER_TO_INT(cut_profile(private_args));
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

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free_cut_args((cut_struct*) args->user);
		free_generic_seq_args(args, TRUE);
	}
}
