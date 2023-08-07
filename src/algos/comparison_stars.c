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
#include "io/remote_catalogues.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"

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

static gboolean is_inside_border(fits *fit, double ra, double dec, double offset_ratio) {
	return wcs2pix_off(fit, ra, dec, NULL, NULL, offset_ratio) == 0;
}

/* com.pref.phot_set.outer needs to be set at this point */
int parse_nina_stars_file_using_WCS(struct light_curve_args *args, const char *file_path,
		gboolean use_comp1, gboolean use_comp2, fits *first) {
	/* The file is a CSV with these fields:
	 * Type,Name,HFR,xPos,yPos,AvgBright,MaxBright,Background,Ra,Dec
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
	FILE *fd = g_fopen(file_path, "r");
	if (!fd) {
		siril_log_message(_("Could not open file %s: %s\n"), file_path, strerror(errno));
		return 1;
	}

	char buf[512];
	rectangle *areas = malloc(MAX_REF_STARS * sizeof(rectangle));
	areas[0].x = 0; areas[0].y = 0;
	int ra_index = -1, dec_index = -1, name_index = 2;
	int stars_count = 0;
	gboolean ready_to_parse = FALSE, target_acquired = FALSE;
	while (fgets(buf, 512, fd)) {
		if (buf[0] == '\0' || buf[0] == '\r' || buf[0] == '\n')
			continue;
		remove_trailing_eol(buf);
		if (buf[0] == '#') {
			if (g_str_has_prefix(buf, "# Siril: ")) {
				args->metadata = calloc(1, sizeof(struct compstars_arg));
				/* parse the siril header */
				int nb = sscanf(buf, "# Siril: %d stars, dVmag %lf, dBV %lf",
						&args->metadata->nb_comp_stars,
						&args->metadata->delta_Vmag,
						&args->metadata->delta_BV);
				if (nb == 0) {
					free(args->metadata);
					args->metadata = NULL;
				}
			}
			continue;
		}
		gchar **tokens = g_strsplit(buf, ",", -1);
		int length = g_strv_length(tokens);
		if (!tokens[0] || length <= ra_index || length <= dec_index) {
			siril_debug_print("malformed line: %s\n", buf);
			g_strfreev(tokens);
			continue;
		}
		gchar *type = tokens[0];
		if (!ready_to_parse) {
			if (!strcasecmp(type, "type")) {
				siril_debug_print("header from the NINA file: %s\n", buf);
				for (int i = 1; tokens[i]; i++) {
					if (!strcasecmp(tokens[i], "ra"))
						ra_index = i;
					else if (!strcasecmp(tokens[i], "dec"))
						dec_index = i;
					else if (!strcasecmp(tokens[i], "name"))
						name_index = i;
				}
				g_strfreev(tokens);

				if (ra_index < 1 || dec_index < 1) {
					siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
					target_acquired = FALSE;
					break;
				}
				siril_debug_print("Found RA and Dec indices in file: %d and %d\n", ra_index, dec_index);
				ready_to_parse = TRUE;
				continue;
			}
			siril_debug_print("malformed line: %s\n", buf);
			siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
			g_strfreev(tokens);
			target_acquired = FALSE;
			break;
		}

		if (!strcasecmp(type, "#")) continue;	// skip comment line

		if (!strcasecmp(type, "target")) {
			gchar *end1, *end2;
			double ra = g_ascii_strtod(tokens[ra_index], &end1);
			double dec = g_ascii_strtod(tokens[dec_index], &end2);
			if (end1 == tokens[ra_index] || end2 == tokens[dec_index]) {
				siril_debug_print("malformed line: %s\n", buf);
				siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
				g_strfreev(tokens);
				target_acquired = FALSE;
				break;
			}

			args->target_descr = g_strdup(tokens[name_index]);
			if (!get_photo_area_from_ra_dec(first, ra, dec, &areas[0])) {
				target_acquired = TRUE;
				stars_count++;
				siril_log_message(_("Target star identified: %s\n"), tokens[name_index]);
			} else {
				siril_log_message(_("There was a problem finding the target star in the image, cannot continue with the light curve\n"));
			}
		}
		else if (!strcasecmp(type, "var")) {
			// we don't use them for this, but we could add them in the
			// user catalogue for annotations, or a local database
		}
		else if (!strcasecmp(type, "comp1")) {
			if (!use_comp1) {
				siril_debug_print("ignoring comp1 star\n");
				g_strfreev(tokens);
				continue;
			}
			gchar *end1, *end2;
			double ra = g_ascii_strtod(tokens[ra_index], &end1);
			double dec = g_ascii_strtod(tokens[dec_index], &end2);
			if (end1 == tokens[ra_index] || end2 == tokens[dec_index]) {
				siril_debug_print("malformed line: %s\n", buf);
				siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
				g_strfreev(tokens);
				target_acquired = FALSE;
				break;
			}
			int index = target_acquired ? stars_count : stars_count + 1;
			if (!get_photo_area_from_ra_dec(first, ra, dec, &areas[index])) {
				if (area_is_unique(&areas[index], areas, index)) {
					stars_count++;
					siril_log_message(_("Star %s added as a reference star\n"), tokens[name_index]);
				}
				else siril_log_message(_("Star %s ignored because it was too close to another\n"), tokens[name_index]);
			}
			else siril_log_message(_("Star %s could not be used because it's on the borders or outside\n"), tokens[name_index]);
		}
		else if (!strcasecmp(type, "comp2")) {
			if (!use_comp2) {
				siril_debug_print("ignoring comp2 star\n");
				g_strfreev(tokens);
				continue;
			}
			gchar *end1, *end2;
			double ra = g_ascii_strtod(tokens[ra_index], &end1);
			double dec = g_ascii_strtod(tokens[dec_index], &end2);
			if (end1 == tokens[ra_index] || end2 == tokens[dec_index]) {
				siril_debug_print("malformed line: %s\n", buf);
				siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
				g_strfreev(tokens);
				target_acquired = FALSE;
				break;
			}
			int index = target_acquired ? stars_count : stars_count + 1;
			if (!get_photo_area_from_ra_dec(first, ra, dec, &areas[index])) {
				if (area_is_unique(&areas[index], areas, index)) {
					stars_count++;
					siril_log_message(_("Star %s added as a reference star\n"), tokens[name_index]);
				}
				else siril_log_message(_("Star %s ignored because it was too close to another\n"), tokens[name_index]);
			}
			else siril_log_message(_("Star %s could not be used because it's on the borders or outside\n"), tokens[name_index]);
		}
		else {
			siril_debug_print("malformed line: %s\n", buf);
		}
		g_strfreev(tokens);
		if (stars_count >= MAX_REF_STARS)
			break;
	}
	if (target_acquired) {
		args->areas = areas;
		args->nb = stars_count;
	} else {
		free(areas);
	}
	fclose(fd);
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
					com.stars[i]->xpos = x;
					com.stars[i]->ypos = gfit.ry - y - 1;
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

	for (int i = 0; i < args->nb_cat_stars; i++)
		if (args->cat_stars[i].star_name)
			g_free(args->cat_stars[i].star_name);
	free(args->cat_stars);

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
static gboolean is_same_star(psf_star *s1, psf_star *s2) {
	return (fabs(s1->ra - s2->ra) < 2.0 * ONE_ARCSEC) &&
			(fabs(s1->dec - s2->dec) < 2.0 * ONE_ARCSEC);
}

int sort_compstars(struct compstars_arg *args) {
	if (!args->target_star || !args->cat_stars || args->nb_cat_stars <= 0) {
		siril_log_color_message(_("Not enough stars available\n"), "salmon");
		return 1;
	}
	siril_log_message(_("Sort parameters: delta_Vmag: %1.2lf, delta_BV: %1.2lf\n"),
			args->delta_Vmag, args->delta_BV);
	siril_log_message(_("Target star: Vmag=%2.2lf, Bmag=%2.2lf\n"),
			args->target_star->mag, args->target_star->Bmag);
	int nb_phot_stars = 0;
	args->comp_stars = malloc((args->nb_cat_stars + 1) * sizeof(psf_star *));
	double BV0 = fabs(args->target_star->mag - args->target_star->Bmag);	// B-V of the target star

	for (int i = 0; i < args->nb_cat_stars; i++) {
		double d_mag = fabs(args->cat_stars[i].mag - args->target_star->mag);
		double BVi = fabs(args->cat_stars[i].mag - args->cat_stars[i].Bmag);

		// Criteria #0: the star has to be within the image and far from the borders
		// (discards 15% of the width/height on both borders)
		if (is_inside_border(&gfit, args->cat_stars[i].ra, args->cat_stars[i].dec, 0.15) &&
				d_mag <= args->delta_Vmag &&		// Criteria #1: nearly same V magnitude
				fabs(BVi - BV0) <= args->delta_BV &&	// Criteria #2: nearly same colors
				!is_same_star(args->target_star, &args->cat_stars[i]) &&
				args->cat_stars[i].s_mag < 0.015 && args->cat_stars[i].s_Bmag < 0.015) {	// Criteria #3: e_Vmag and e_Bmag < 0.015
			args->comp_stars[nb_phot_stars] = duplicate_psf(&args->cat_stars[i]);
			if (!args->comp_stars[nb_phot_stars]->star_name)
				args->comp_stars[nb_phot_stars]->star_name = g_strdup_printf("%d", nb_phot_stars+1);
			nb_phot_stars++;
		}
		i++;
	}
	args->comp_stars[nb_phot_stars] = NULL;

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
	SirilWorldCS *center = siril_world_cs_new_from_a_d(ra, dec);

	double resolution = get_wcs_image_resolution(&gfit);
	uint64_t sqr_radius = 0;
	double radius = 0.0;

	if (args->narrow_fov) {
		// Limited to ry, to avoid the corners with their potential vignettage
		radius = resolution * gfit.ry / 2.0;	// in degres
	} else {
		// The whole Field of View
		sqr_radius = (gfit.rx * gfit.rx + gfit.ry * gfit.ry) / 4;
		radius = resolution * sqrt((double)sqr_radius);	// in degrees
	}

	double limit_mag = max(args->target_star->mag + 6.0, 17.0);
	GFile *catalog_file = download_catalog(args->cat, center, radius * 60.0, limit_mag);
	siril_world_cs_unref(center);
	if (!catalog_file) {
		siril_log_message(_("Could not download the online star catalog.\n"));
		return 1;
	}
	siril_log_message(_("The %s catalog has been successfully downloaded.\n"),
			catalog_to_str(args->cat));

	load_catalog(catalog_file, TRUE, &args->cat_stars, &args->nb_cat_stars);

	g_object_unref(catalog_file);
	return 0;
}

static void write_nina_file(struct compstars_arg *args) {
	if (!args->nina_file)
		return;
	if (!g_strcmp0(args->nina_file, "auto")) {
		g_free(args->nina_file);
		args->nina_file = g_strdup_printf("%s_SirilstarList_%1.2lf_%1.2lf_%s.csv",
				args->target_star->star_name, args->delta_Vmag, args->delta_BV,
				catalog_to_str(args->cat));
	}
	FILE *fd = g_fopen(args->nina_file, "w+");
	if (!fd)
		return;
	siril_log_message(_("Creating csv output file %s\n"), args->nina_file);

	fprintf(fd, "# Sorted comparison stars for %s from %s according to the following criteria\n",
			args->target_star->star_name, catalog_to_str(args->cat));

	if (args->cat == CAT_AAVSO && args->AAVSO_chartid && args->AAVSO_uri) {
		fprintf(fd, "# AAVSO chartid: %s, image_uri: %s\n",
				args->AAVSO_chartid, args->AAVSO_uri);
	}

	if (args->delta_Vmag <= 0.0 && args->delta_BV <= 0.0)
		fprintf(fd, "# No criteria applied\n");
	else fprintf(fd, "# Siril: %d stars, dVmag %.2f, dBV %.2f\n",
			args->nb_comp_stars, args->delta_Vmag, args->delta_BV);

	fprintf(fd, "Type,Name,HFR,xPos,yPos,AvgBright,MaxBright,Background,Ra,Dec\n");

	fprintf(fd, "Target,%s,,,,,,,%lf,%lf\n", args->target_star->star_name,
			args->target_star->ra, args->target_star->dec);

	int i = 0;
	while (args->comp_stars[i]) {
		int retval = fprintf(fd, "Comp1,%s,,,,,,,%f,%f\n",
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
	siril_log_message(_("e_Vmag and e_Bmag are intrinsec errors supplied from %s\n"), catalog_to_str(args->cat));
	siril_log_message("Index_nbr / Vmag / Bmag / d_mag / d_BV / e_Vmag / e_Bmag\n");
	siril_log_message(_("Target star: %2.2lf, %2.2lf, 0, %1.3lf, %1.3lf, %1.3lf\n"),
			args->target_star->mag, args->target_star->Bmag,
			fabs(fabs(args->target_star->mag - args->target_star->Bmag) - 0.7),
			args->target_star->s_mag, args->target_star->s_Bmag);

	int i = 0;
	int Cmp_St = 0;
	while (args->comp_stars[i]) {
		const char *title = "Comp star";
		Cmp_St++;
		siril_log_message(_("%s  %d: %2.2lf, %2.2lf, %1.3lf, %1.3lf, %1.3lf, %1.3lf\n"),
				title, i+1, args->comp_stars[i]->mag, args->comp_stars[i]->Bmag,
				fabs(args->comp_stars[i]->mag - args->target_star->mag),
				fabs(fabs(args->comp_stars[i]->mag - args->comp_stars[i]->Bmag) - 0.7),
				args->comp_stars[i]->s_mag, args->comp_stars[i]->s_Bmag);
		i++;
	}

	if (Cmp_St) {
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
gchar *generate_lc_subtitle(struct compstars_arg *metadata, gboolean for_gnuplot) {
	if (!metadata)
		return g_strdup("");
	GString *str = for_gnuplot ? g_string_new("\\n{/*0.8 ") : g_string_new("");
	gboolean first = TRUE;
	if (metadata->nb_comp_stars > 0 && metadata->delta_Vmag != 0.0 && metadata->delta_BV != 0.0) {
		if (for_gnuplot) {
#ifdef _WIN32
			g_string_append_printf(str,
					"%d stars within delta_{Vmag} = %.2f, delta_{BV} = %.2f",
					metadata->nb_comp_stars, metadata->delta_Vmag, metadata->delta_BV);
#else
			g_string_append_printf(str,
					"%d stars within {/Symbol d}_{Vmag} = %.2f, {/Symbol d}_{BV} = %.2f",
					metadata->nb_comp_stars, metadata->delta_Vmag, metadata->delta_BV);
#endif
		}
		else g_string_append_printf(str, "%d stars within delta Vmag = %.2f, delta BV = %.2f",
					metadata->nb_comp_stars, metadata->delta_Vmag, metadata->delta_BV);
		first = FALSE;
	}
	if (metadata->AAVSO_chartid) {
		g_string_append_printf(str, "%sAAVSO chart %s", first ? "" : " - ", metadata->AAVSO_chartid);
	}
	if (for_gnuplot)
		g_string_append(str, "}");
	else g_string_append(str, "\n# "); // in that case it's the first line...

	return g_string_free(str, FALSE);
}

