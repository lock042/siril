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
#include <gtk/gtk.h>
#include "gui/utils.h"

#include "comparison_stars.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "algos/search_objects.h"
#include "io/siril_catalogues.h"
#include "io/remote_catalogues.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"

#define BORDER_RATIO 0.10 // the amount of image that is considered at border

// area is not recovering another at least half their size (assumed square and identical)
static int area_is_unique(rectangle *area, rectangle *areas, int nb_areas) {
	int half_size = area->w / 2;
	for (int i = 0; i < nb_areas; i++) {
		if (abs(area->x - areas[i].x) < half_size && abs(area->y - areas[i].y) < half_size)
			return 0;
	}
	return 1;
}

static int get_photo_area_from_ra_dec(fits *fit, double ra, double dec, rectangle *ret_area) {
	double fx, fy;
	if (wcs2pix(fit, ra, dec, &fx, &fy)) {
		siril_debug_print("star is outside image\n");
		return 1;
	}
	double dx, dy;
	fits_to_display(fx, fy, &dx, &dy, fit->ry);
	// some margin needs to be allowed for star centroid movement
	double start = 1.5 * com.pref.phot_set.outer;
	double size = 3.0 * com.pref.phot_set.outer;
	rectangle area;
	area.x = round_to_int(dx - start);
	area.y = round_to_int(dy - start);
	area.w = size;
	area.h = size;
	if (area.x < 0 || area.y < 0 ||
			area.h <= 0 || area.w <= 0 ||
			area.x + area.w >= fit->rx ||
			area.y + area.h >= fit->ry) {
		siril_debug_print("star is outside image\n");
		return 1;
	}
	*ret_area = area;
	siril_debug_print("Pixel coordinates of a star: %.1f, %.1f\n", dx, dy);
	return 0;
}

/* com.pref.phot_set.outer needs to be set at this point */
int parse_nina_stars_file_using_WCS(struct light_curve_args *args, const char *file_path,
		gboolean use_comp1, gboolean use_comp2, fits *first) {
	/* The file is a CSV with at least these fields:
	 * Type,Name, Ra,Dec
	 *
	 * Type can be 'Target' for the variable star to analyse, 'Var' for variable stars to
	 * absolutely exclude as calibration reference, 'Comp1' are reference stars obtained
	 * from SIMBAD based on color, 'Comp2' are stars obtained from the AAVSO site.
	 * We just need this and Ra,Dec. xPos and yPos are in pixels, but we'll use wcs2pix
	 * to get them with our plate solve and our star fitting.
	 * There can be comments at the top for Siril-generated files, starting with a '#'.
	 * In the comments we put some metadata, like parameters that were used to select
	 * the stars of the list. We can parse them to keep them within the light curve.
	 */
	siril_catalogue *siril_cat = calloc(1, sizeof(siril_catalogue));
	siril_cat->cattype = CAT_COMPSTARS;
	siril_cat->columns = siril_catalog_columns(siril_cat->cattype);
	siril_catalog_load_from_file(siril_cat, file_path);
	
	// parsing the file header to get metadata
	if (siril_cat->header) { 
		gchar **tokens = g_strsplit(siril_cat->header, "\n", -1);
		int length = g_strv_length(tokens);
		for (int i = 0; i < length; i++) {
			if (g_str_has_prefix(tokens[i], "# Siril: ")) {
				args->metadata = calloc(1, sizeof(struct compstars_arg));
				/* parse the siril header */
				int nb = sscanf(tokens[i], "# Siril: %d stars, dVmag %lf, dBV %lf, max e_mag %lf",
						&args->metadata->nb_comp_stars,
						&args->metadata->delta_Vmag,
						&args->metadata->delta_BV,
						&args->metadata->max_emag);
				if (nb != 4) {
					free(args->metadata);
					args->metadata = NULL;
				}
			}
		}
		g_strfreev(tokens);
	}

	rectangle *areas = malloc(MAX_REF_STARS * sizeof(rectangle));
	areas[0].x = 0; areas[0].y = 0;
	int stars_count = 0;
	gboolean target_acquired = FALSE;

	for (int i = 0; i < siril_cat->nbitems; i++) {
		cat_item item = siril_cat->cat_items[i];
		if (!strcasecmp(item.type, "target")) {
			args->target_descr = g_strdup(item.name);
			if (!get_photo_area_from_ra_dec(first, item.ra, item.dec, &areas[0])) {
				target_acquired = TRUE;
				stars_count++;
				siril_log_message(_("Target star identified: %s\n"), item.name);
			} else {
				siril_log_color_message(_("There was a problem finding the target star in the image, cannot continue with the light curve\n"), "red");
				break;
			}
		} else if ((use_comp1 && !strcasecmp(item.type, "comp1")) || (use_comp2 && !strcasecmp(item.type, "comp2"))) {
			int index = target_acquired ? stars_count : stars_count + 1;
			if (!get_photo_area_from_ra_dec(first, item.ra, item.dec, &areas[index])) {
				if (area_is_unique(&areas[index], areas, index)) {
					stars_count++;
					siril_log_message(_("star %s [%s] added as a reference star\n"), item.name, item.type);
				}
				else siril_log_message(_("Star %s ignored because it was too close to another\n"), item.name);
			}
			else siril_log_message(_("Star %s could not be used because it's on the borders or outside\n"), item.name);
		}
		if (stars_count == MAX_REF_STARS)
			break;
	}

	if (target_acquired) {
		args->areas = areas;
		args->nb = stars_count;
	} else {
		free(areas);
	}
	siril_catalog_free(siril_cat);
	return !target_acquired;
}

// may not be called as an idle, do not call GTK code if !has_GUI
static gboolean end_compstars(gpointer p) {
	siril_debug_print("end_compstars\n");
	struct compstars_arg *args = (struct compstars_arg *) p;
	g_free(args->target_name);

	// display the comp star list
	clear_stars_list(args->has_GUI);
	if (args->has_GUI) {
		if (!args->retval)
			purge_temp_user_catalogue();
		if (!args->nina_file) {
			com.stars = args->comp_stars;
			siril_debug_print("displaying results as stars\n");
			for (int i = 0; i < args->nb_comp_stars; i++) {
				double x, y;
				if (!wcs2pix(&gfit, com.stars[i]->ra, com.stars[i]->dec, &x, &y)) {
					fits_to_display(x, y, &com.stars[i]->xpos, &com.stars[i]->ypos, gfit.ry);
				}
			}
		}
		else if (!args->retval) {
			load_csv_targets_to_temp(args->nina_file);
			/* display the new 'found_object' */
			GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
			refresh_found_objects();
			if (!gtk_toggle_tool_button_get_active(button)) {
				gtk_toggle_tool_button_set_active(button, TRUE);
			} else {
				refresh_found_objects();
				redraw(REDRAW_OVERLAY);
			}
		}
	}
	g_free(args->nina_file);

	siril_catalog_free(args->cat_stars);

	if (args->target_star)
		free_psf(args->target_star);
	g_free(args->AAVSO_chartid);
	g_free(args->AAVSO_uri);
	redraw(REDRAW_OVERLAY);
	free(args);
	return end_generic(NULL);
}


#define ONE_ARCSEC 0.000277778
/* determines if two stars are the same based on their coordinates */
static gboolean is_same_star(psf_star *s1, cat_item *s2) {
	return (fabs(s1->ra - s2->ra) < 2.0 * ONE_ARCSEC) &&
			(fabs(s1->dec - s2->dec) < 2.0 * ONE_ARCSEC);
}

int sort_compstars(struct compstars_arg *args) {
	if (!args->target_star || !args->cat_stars || args->cat_stars->nbitems <= 0) {
		siril_log_color_message(_("Not enough stars available\n"), "salmon");
		return 1;
	}
	siril_catalogue *siril_cat = args->cat_stars;
	siril_log_message(_("Sort parameters: delta_Vmag: %1.2lf, delta_BV: %1.2lf, max e_Vmag: %1.3lf\n"),
			args->delta_Vmag, args->delta_BV, args->max_emag);
	siril_log_message(_("Target star: Vmag = %2.2lf, B-V = %+2.2lf\n"),
			args->target_star->mag, args->target_star->Bmag - args->target_star->mag);

	int nb_phot_stars = 0;
	args->comp_stars = malloc((siril_cat->nbincluded + 1) * sizeof(psf_star *));
	double BV0 = fabs(args->target_star->mag - args->target_star->Bmag);	// B-V of the target star
	double xmin, xmax, ymin, ymax;
	xmin = (double)gfit.rx * BORDER_RATIO;
	xmax = (double)gfit.rx * (1. - BORDER_RATIO);
	ymin = (double)gfit.ry * BORDER_RATIO;
	ymax = (double)gfit.ry * (1. - BORDER_RATIO);


	for (int i = 0; i < args->nb_cat_stars; i++) {
		if (!siril_cat->cat_items[i].included) // included means inside the image
			continue;
		double d_mag = fabs(siril_cat->cat_items[i].mag - args->target_star->mag);
		// TODO: the fabs should not be there, it's B-V to get the color
		//double BVi = fabs(siril_cat->cat_items[i].mag - siril_cat->cat_items[i].bmag);
		double BVi = siril_cat->cat_items[i].bmag - siril_cat->cat_items[i].mag; // B-V index

		// Criteria #0: the star has to be within the image and far from the borders
		// (discards 10% of the width/height on both borders)
		if ((siril_cat->cat_items[i].x > xmin && siril_cat->cat_items[i].x < xmax && siril_cat->cat_items[i].y > ymin && siril_cat->cat_items[i].y < ymax) &&
				d_mag <= args->delta_Vmag &&		// Criteria #1: nearly same V magnitude
				// TODO: AAVSO guide says the compstars do not need to be the same color as target
				// but that their BV index should be between +0.3 and +1.0, preferably around +0.7, do we need to do smthg with that info?
				fabs(BVi - BV0) <= args->delta_BV &&	// Criteria #2: nearly same colors
				!is_same_star(args->target_star, &siril_cat->cat_items[i]) &&
				// TODO: AAVSO guide mentions compstars should have a low magnitude error, not sure it means in Bband
				// Here this criterion (which is not optional) discards most stars and forces to increase d_mag
				// Stars in AAVSO charts commonly have Vmag above 0.02
				// siril_cat->cat_items[i].e_mag < 0.015 && siril_cat->cat_items[i].e_bmag < 0.015) {	// Criteria #3: e_Vmag and e_Bmag < 0.015
				siril_cat->cat_items[i].e_mag < args->max_emag) {
			// TODO: we could use a siril_cat output, and provided csv_writer later on
			args->comp_stars[nb_phot_stars] = new_psf_star();
			args->comp_stars[nb_phot_stars]->xpos = siril_cat->cat_items[i].x;
			args->comp_stars[nb_phot_stars]->ypos = siril_cat->cat_items[i].y;
			args->comp_stars[nb_phot_stars]->mag = siril_cat->cat_items[i].mag;
			args->comp_stars[nb_phot_stars]->Bmag = siril_cat->cat_items[i].bmag;
			args->comp_stars[nb_phot_stars]->s_mag = siril_cat->cat_items[i].e_mag;
			args->comp_stars[nb_phot_stars]->s_Bmag = siril_cat->cat_items[i].e_bmag;
			args->comp_stars[nb_phot_stars]->ra = siril_cat->cat_items[i].ra;
			args->comp_stars[nb_phot_stars]->dec = siril_cat->cat_items[i].dec;
			if (siril_cat->cat_items[i].name)
				args->comp_stars[nb_phot_stars]->star_name = g_strdup(siril_cat->cat_items[i].name);
			else
				args->comp_stars[nb_phot_stars]->star_name = g_strdup_printf("%d", nb_phot_stars + 1);
			nb_phot_stars++;
		}
		i++;
	}
	args->comp_stars[nb_phot_stars] = NULL;
	// reallocating list to the required size
	psf_star **result = realloc(args->comp_stars, (nb_phot_stars + 1) * sizeof(psf_star *));
	if (!result) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	args->comp_stars = result;
	siril_log_message(_("%d comparison stars after sort.\n"), nb_phot_stars);
	args->nb_comp_stars = nb_phot_stars;
	return nb_phot_stars == 0;
}

static int get_catstars(struct compstars_arg *args) {
	if (!args->target_star) {
		siril_log_color_message(_("No variable star selected\n"), "salmon");
		return 1;
	}
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return 1;
	}
	g_assert(args->cat == CAT_APASS || args->cat == CAT_NOMAD);

	double ra, dec;
	center2wcs(&gfit, &ra, &dec);

	double resolution = get_wcs_image_resolution(&gfit);
	uint64_t sqr_radius = 0;
	double radius = 0.0;

	if (args->narrow_fov) {
		// Limited to the image smallest dimension, to avoid the corners with their potential vignettage
		radius = resolution * min(gfit.rx, gfit.ry) / 2.0;	// in degrees
	} else {
		// The whole Field of View
		sqr_radius = (gfit.rx * gfit.rx + gfit.ry * gfit.ry) / 4;
		radius = resolution * sqrt((double)sqr_radius);	// in degrees
	}

	// preparing the query
	siril_catalogue *siril_cat = siril_catalog_fill_from_fit(&gfit, args->cat, max(args->target_star->mag + 6.0, 17.0));
	siril_cat->radius = radius * 60.; // overwriting to account for narrow argument
	siril_cat->phot = TRUE;

	// and retrieving its results
	if (!siril_catalog_conesearch(siril_cat)) // returns the nb of stars
		return 1;
	sort_cat_items_by_mag(siril_cat);

	if (!siril_catalog_project_with_WCS(siril_cat, &gfit, FALSE)) {
		args->nb_cat_stars = siril_cat->nbitems;
		args->cat_stars = siril_cat;
	}
	return (int)(args->nb_cat_stars == 0);
}

static void write_nina_file(struct compstars_arg *args) {
	if (!args->nina_file)
		return;
	if (!g_strcmp0(args->nina_file, "auto")) {
		g_free(args->nina_file);
		args->nina_file = g_strdup_printf("%s_SirilstarList_%1.2lf_%1.2lf_%1.2lf_%s.csv",
				args->target_star->star_name, args->delta_Vmag, args->delta_BV, args->max_emag,
				catalog_to_str(args->cat));
	}
	FILE *fd = g_fopen(args->nina_file, "w+");
	if (!fd)
		return;
	siril_log_message(_("Creating csv output file %s\n"), args->nina_file);

	fprintf(fd, "# Sorted comparison stars for %s from %s according to the following criteria\n",
			args->target_star->star_name, catalog_to_str(args->cat));

	// if (args->cat == CAT_AAVSO && args->AAVSO_chartid && args->AAVSO_uri) {
	// 	fprintf(fd, "# AAVSO chartid: %s, image_uri: %s\n",
	// 			args->AAVSO_chartid, args->AAVSO_uri);
	// }

	if (args->delta_Vmag <= 0.0 && args->delta_BV <= 0.0)
		fprintf(fd, "# No criteria applied\n");
	else fprintf(fd, "# Siril: %d stars, dVmag %.2f, dBV %.2f, max e_mag %.2f\n",
			args->nb_comp_stars, args->delta_Vmag, args->delta_BV, args->max_emag);

	fprintf(fd, "Type,Name,Ra,Dec\n");

	fprintf(fd, "Target,%s,%lf,%lf\n", args->target_star->star_name,
			args->target_star->ra, args->target_star->dec);

	int i = 0;
	while (args->comp_stars[i]) {
		int retval = fprintf(fd, "Comp1,%s,%f,%f\n",
				args->comp_stars[i]->star_name /*? args->comp_stars[i]->star_name : "x"*/,
				args->comp_stars[i]->ra, args->comp_stars[i]->dec);
		if (retval < 0)
			break;
		i++;
	}

	fclose(fd);
}

void chk_compstars(struct compstars_arg *args) {
	if (!args->target_star) {
		siril_log_color_message(_("No variable star selected\n"), "salmon");
		return;
	}

	siril_log_message("d_mag and d_BV are discrepancies from Vmag and B-V of the target star\n");
	siril_log_message(_("e_Vmag and e_Bmag are photometric errors supplied from %s\n"), catalog_to_str(args->cat));
	siril_log_message("Index_nbr        V     B-V   d_mag    d_BV  e_Vmag e_Bmag\n");
	int i = 0;
	double BV0 = args->target_star->Bmag - args->target_star->mag;
	while (args->comp_stars[i]) {
		const char *title = "Comp star";
		double BVi = args->comp_stars[i]->Bmag - args->comp_stars[i]->mag;
		siril_log_message(_("%s %3d: %4.2lf, %+4.2lf, %+5.3lf, %+5.3lf, %4.3lf, %4.3lf\n"),
				title, i + 1, args->comp_stars[i]->mag, BVi,
				args->comp_stars[i]->mag - args->target_star->mag,
				BVi - BV0,
				args->comp_stars[i]->s_mag, args->comp_stars[i]->s_Bmag);
		i++;
	}

	if (args->nb_comp_stars) {
		write_nina_file(args);
	} else siril_log_message(_("No csv output file to create\n"));
}

// uses gfit
gpointer compstars_worker(gpointer arg) {
	int retval;
	struct compstars_arg *args = (struct compstars_arg *) arg;
	// 1. search for the variable star
	retval = cached_object_lookup(args->target_name, &args->target_star);
	if (retval)
		goto end;
	if (!args->target_star) {
		siril_log_color_message(_("Target star invalid.\n"), "red");
		retval = 1;
		goto end;
	}
	if (args->target_star->mag == 0.0 || args->target_star->Bmag == 0.0) {
		siril_log_color_message(_("Target star photometric information not available.\n"), "red");
		retval = 1;
		goto end;
	}

	// 2. get a catalogue of stars for the field
	retval = get_catstars(args);
	if (retval) goto end;

	// 3. extract those matching our photometry need
	retval = sort_compstars(args);
	if (retval) goto end;

	// generate NINA-style star list
	chk_compstars(args);

end:
	args->retval = retval;
	args->has_GUI = TRUE;
	if (!siril_add_idle(end_compstars, args)) {
		args->has_GUI = FALSE;
		end_compstars(args);	// we still need to free all
	}
	return GINT_TO_POINTER(retval);
}

/* generates subtitles for the dat file header and the plots based on what's available */
gchar *generate_lc_subtitle(struct compstars_arg *metadata, gboolean for_plot) {
	if (!metadata)
		return g_strdup("");
	GString *str = g_string_new("");
	gboolean first = TRUE;
	if (metadata->nb_comp_stars > 0 && metadata->delta_Vmag != 0.0 && metadata->delta_BV != 0.0 && metadata->max_emag != 0.0) {
		if (for_plot)
			g_string_append_printf(str,
				"\n<span size=\"small\">%d %s &#x03B4;<sub>Vmag</sub> = %.2f, &#x03B4;<sub>BV</sub> = %.2f, max e_mag = %.2f</span>",
				metadata->nb_comp_stars, _("stars within"), metadata->delta_Vmag, metadata->delta_BV, metadata->max_emag);
		else g_string_append_printf(str, "#%d %s delta Vmag = %.2f, delta BV = %.2f, max e_mag = %.2f",
					metadata->nb_comp_stars, _("stars within"), metadata->delta_Vmag, metadata->delta_BV, metadata->max_emag);
		first = FALSE;
	}
	if (metadata->AAVSO_chartid) {
		if (for_plot)
			g_string_append_printf(str, "<span size=\"small\">%sAAVSO chart %s</span>", first ? "\n" : " - ", metadata->AAVSO_chartid);
		else
			g_string_append_printf(str, "%sAAVSO chart %s", first ? "#" : " - ", metadata->AAVSO_chartid);
		first = FALSE;
	}
	if (!for_plot && !first)
		g_string_append(str, "\n");
	return g_string_free(str, FALSE);
}

