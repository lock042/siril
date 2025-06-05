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
#include <gtk/gtk.h>
#include "gui/utils.h"

#include "comparison_stars.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/siril_wcs.h"
#include "algos/search_objects.h"
#include "io/siril_catalogues.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"


#define BORDER_RATIO 0.10 // the amount of image that is considered at border
#define MAX_VAR_CAT 3 // max number of variable stars catalogues
#define MAX_SEPARATION  2 * 0.000277778 // the max sepration to consider stars are the same
const siril_cat_index var_cat[MAX_VAR_CAT] = { CAT_VSX, CAT_GCVS, CAT_VARISUM }; // the variable stars catalogues, order to be consistent with discard_var_catalogues

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
	siril_to_display(fx, fy, &dx, &dy, fit->ry);
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
	siril_catalogue *siril_cat = siril_catalog_new(CAT_COMPSTARS);
	if (siril_catalog_load_from_file(siril_cat, file_path)) {
		siril_catalog_free(siril_cat);
		return 1;
	}

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
		cat_item *item = &siril_cat->cat_items[i];
		if (!g_ascii_strcasecmp(item->type, "target")) {
			args->target_descr = g_strdup(item->name);
			if (!get_photo_area_from_ra_dec(first, item->ra, item->dec, &areas[0])) {
				target_acquired = TRUE;
				stars_count++;
				siril_log_message(_("Target star identified: %s\n"), item->name);
			} else {
				siril_log_color_message(_("There was a problem finding the target star in the image, cannot continue with the light curve\n"), "red");
				break;
			}
		} else if ((use_comp1 && !g_ascii_strcasecmp(item->type, "comp1")) || (use_comp2 && !g_ascii_strcasecmp(item->type, "comp2"))) {
			int index = target_acquired ? stars_count : stars_count + 1;
			if (!get_photo_area_from_ra_dec(first, item->ra, item->dec, &areas[index])) {
				if (area_is_unique(&areas[index], areas, index)) {
					stars_count++;
					siril_log_message(_("star %s [%s] added as a reference star\n"), item->name, item->type);
				}
				else siril_log_message(_("Star %s ignored because it was too close to another\n"), item->name);
			}
			else siril_log_message(_("Star %s could not be used because it's on the borders or outside\n"), item->name);
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
	if (args->has_GUI && !args->retval) {
		purge_user_catalogue(CAT_AN_USER_TEMP); // we always clear user_temp in case there are two successive calls to findcompstars
		if (!load_siril_cat_to_temp(args->comp_stars)) {
			GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
			refresh_found_objects();
			if (!gtk_toggle_tool_button_get_active(button)) {
				gtk_toggle_tool_button_set_active(button, TRUE);
			} else {
				refresh_found_objects();
				redraw(REDRAW_OVERLAY);
			}
		}
	} else {
		siril_catalog_free(args->comp_stars); // we don't free args->compstars if it's being used for the annotations
	}
	g_free(args->nina_file);
	siril_catalog_free(args->cat_stars);
	if (args->target_star)
		siril_catalog_free_item(args->target_star);
	g_free(args->AAVSO_chartid);
	g_free(args->AAVSO_uri);
	redraw(REDRAW_OVERLAY);
	free(args);
	return end_generic(NULL);
}

void write_nina_file(struct compstars_arg *args) {
	if (!args->nina_file)
		return;
	if (!g_strcmp0(args->nina_file, "auto")) {
		g_free(args->nina_file);
		args->nina_file = g_strdup_printf("%s_SirilstarList_%1.2lf_%1.2lf_%1.2lf_%s.csv",
				args->target_star->name, args->delta_Vmag, args->delta_BV, args->max_emag,
				catalog_to_str(args->cat));
	}
	replace_spaces_from_str(args->nina_file, '_');
	replace_char_from_str(args->nina_file, '*', '-');

	siril_log_message(_("Creating csv output file %s\n"), args->nina_file);

	GString *header_lines = g_string_new("");
	g_string_append_printf(header_lines, "# Sorted comparison stars for %s from %s according to the following criteria\n",
			args->target_star->name, catalog_to_str(args->cat));

	// if (args->cat == CAT_AAVSO && args->AAVSO_chartid && args->AAVSO_uri) {
	// 	fprintf(fd, "# AAVSO chartid: %s, image_uri: %s\n",
	// 			args->AAVSO_chartid, args->AAVSO_uri);
	// }

	if (args->delta_Vmag <= 0.0 && args->delta_BV <= 0.0 && args->max_emag <= 0.0)
		g_string_append_printf(header_lines, "# No criteria applied");
	else
		g_string_append_printf(header_lines, "# Siril: %d stars, dVmag %.2f, dBV %.2f, max e_mag %.2f",
			args->nb_comp_stars, args->delta_Vmag, args->delta_BV, args->max_emag);
	args->comp_stars->header = g_string_free(header_lines, FALSE);
	if (!siril_catalog_write_to_file(args->comp_stars, args->nina_file))
		siril_log_color_message(_("Problem writing the comparison stars file\n"), "red");
}

/* determines if two stars are the same based on their coordinates */
static gboolean is_same_star(cat_item *s1, cat_item *s2){
	return (fabs(s1->ra - s2->ra) < MAX_SEPARATION) &&
			(fabs(s1->dec - s2->dec) < MAX_SEPARATION);
}

void fill_compstar_item(cat_item *item, double ra, double dec, float mag, gchar *name, const gchar *type) {
	item->ra = ra;
	item->dec = dec;
	item->mag = mag;
	if (name)
		item->name = g_strdup(name);
	item->type = g_strdup(type);
	item->included = TRUE;
}

// This function compares objects and return them sorted by distance
static int compare_items_by_dist(const void* item1, const void* item2) {
	compstar_dist *i1 = (compstar_dist*) item1;
	compstar_dist *i2 = (compstar_dist*) item2;
	if (i1->dist < i2->dist)
		return -1;
	if (i1->dist > i2->dist)
		return 1;
	return 0;
}

// Checks if a star belongs to a given catalog
static int is_var_star(cat_item *item, GList *siril_cats, int nb_disc[]) {
	if (!siril_cats) return 0;	// Does the catalog-to-be-discarded exist?
	while (siril_cats) {
		int n = 0;
		siril_catalogue *siril_cat = (siril_catalogue *)siril_cats->data;
		for (int i = 0; i < siril_cat->nbitems; i++) {
			if(siril_cat->cat_items[i].included && is_same_star(&siril_cat->cat_items[i], item)) {
				nb_disc[n] += 1;
				return TRUE;
			}
		}
		siril_cats = siril_cats->next;
		n++;
	}
	return FALSE;
}

int sort_compstars(struct compstars_arg *args) {
	if (!args->target_star || !args->cat_stars || args->cat_stars->nbitems <= 0) {
		siril_log_color_message(_("Not enough stars available\n"), "salmon");
		return 1;
	}
	siril_catalogue *siril_cat = args->cat_stars;
	siril_log_message(_("Filtering parameters: delta_Vmag: %1.2lf, delta_BV: %1.2lf, max e_Vmag: %1.3lf\n"),
			args->delta_Vmag, args->delta_BV, args->max_emag);
	siril_log_message(_("Target star: Vmag = %2.2lf, B-V = %+2.2lf\n"),
			args->target_star->mag, args->target_star->bmag - args->target_star->mag);

	// we will just measure distance to center, we will copy them later
	compstar_dist *sorter = calloc(args->cat_stars->nbincluded, sizeof(compstar_dist));
	int nb_phot_stars = 0;

	// prepare the reference values
	double BV0 = args->target_star->bmag - args->target_star->mag;	// B-V of the target star
	double xmin, xmax, ymin, ymax; // borders boundaries

	xmin = (double)args->fit->rx * BORDER_RATIO;
	xmax = (double)args->fit->rx * (1. - BORDER_RATIO);
	ymin = (double)args->fit->ry * BORDER_RATIO;
	ymax = (double)args->fit->ry * (1. - BORDER_RATIO);

	const gchar *startype = (args->cat == CAT_NOMAD || args->cat == CAT_APASS) ? "Comp1" : "Comp2";
	gboolean cat2discard = args->var_stars_cat != NULL;
	int nb_disc[MAX_VAR_CAT] = { 0 };
	for (int i = 0; i < args->cat_stars->nbitems; i++) {
		cat_item *item = &siril_cat->cat_items[i];
		if (!item->included || is_same_star(args->target_star, item)) // included means inside the image after wcs projection
			continue;
		double d_mag = fabs(item->mag - args->target_star->mag);
		double BVi = item->bmag - item->mag; // B-V index

		// Criteria #0: the star has to be within the image and far from the borders
		// (discards BORDER_RATIO of the width/height on both borders)
		if ((item->x > xmin && item->x < xmax && item->y > ymin && item->y < ymax) &&
				d_mag <= args->delta_Vmag &&		// Criteria #1: nearly same V magnitude
				fabs(BVi - BV0) <= args->delta_BV &&	// Criteria #2: nearly same colors
				((args->cat == CAT_APASS) ? (item->e_mag > 0. && item->e_mag <= args->max_emag) : TRUE) &&	// Criteria #3: e_mag smaller than threshold, for catalogues that have the info
				(!cat2discard || !is_var_star(item, args->var_stars_cat, nb_disc))) {// Criteria #4: not a variable star - we search here to avoid comparing long lists to long lists
			sorter[nb_phot_stars] = (compstar_dist){i, compute_coords_distance(siril_cat->center_ra, siril_cat->center_dec, item->ra, item->dec)};
			nb_phot_stars++;
		}
	}

	if (cat2discard) {
		GList *siril_cats = args->var_stars_cat;
		int i = 0;
		while (siril_cats) {
			siril_catalogue *siril_cat = (siril_catalogue *)siril_cats->data;
			if (nb_disc[i] > 0)
				siril_log_message(_("-> %d stars matching those criteria were found in %s and discarded\n"), nb_disc[i++], catalog_to_str(siril_cat->cat_index));
			siril_catalog_free(siril_cat);
			siril_cats = siril_cats->next;
		}
		g_list_free(args->var_stars_cat);
	}

	if (nb_phot_stars > 0) {
		// For now, we sort the stars by increasing radius from center
		// if we later on want to add capability to sort by mag instead
		// we will just add a flag in the args that does or not the sorting by rad
		// as args->cat_stars is sorted by mag
		qsort(sorter, nb_phot_stars, sizeof(compstar_dist), compare_items_by_dist);
		siril_log_message(_("Stars sorted by increasing distance (in arcmin) wrt. image center\n"));
		// preparing the output catalog
		args->comp_stars = siril_catalog_new(CAT_COMPSTARS);
		args->comp_stars->columns |= (1 << CAT_FIELD_MAG); // we add mag to write it in the output file (it is not a mandatory field at readout)
		// allocating final sorted list to the required size
		cat_item *result = calloc(nb_phot_stars + 1, sizeof(cat_item));
		// write the target star
		fill_compstar_item(&result[0], args->target_star->ra, args->target_star->dec, args->target_star->mag, args->target_star->name, "Target");
		// and write the stars sorted by radius
		siril_log_message("d_mag and d_BV are discrepancies from Vmag and B-V of the target star\n");
		siril_log_message(_("e_Vmag and e_Bmag are photometric errors supplied from %s\n"), catalog_to_str(args->cat));
		siril_log_message("Index_nbr        V     B-V   d_mag    d_BV  e_Vmag e_Bmag   Dist\n");
		for (int i = 0; i < nb_phot_stars; i++) {
			cat_item *item = &siril_cat->cat_items[sorter[i].index];
			gchar *name = (item->name) ? g_strdup((item->name)) : g_strdup_printf("%d", i + 1);
			fill_compstar_item(&result[i + 1], item->ra, item->dec, item->mag, name, startype);
			g_free(name);
			siril_log_message(_("Comp star %3d: %4.2lf, %+4.2lf, %+5.3lf, %+5.3lf, %4.3lf, %4.3lf, %5.2lf\n"),
					i + 1, item->mag, item->bmag - item->mag,
					item->mag - args->target_star->mag,
					item->bmag - item->mag - BV0,
					item->e_mag,
					item->e_bmag,
					sorter[i].dist * 60.);
		}
		args->comp_stars->cat_items = result;
		args->comp_stars->nbitems = nb_phot_stars + 1;
		args->comp_stars->nbincluded = nb_phot_stars + 1;
		args->nb_comp_stars = nb_phot_stars;
	}
	// free the intermediate sorting list
	free(sorter);
	// and log results
	siril_log_message(_("%d comparison stars after sort.\n"), nb_phot_stars);
	if (nb_phot_stars && args->nina_file)
		write_nina_file(args);
	else
		siril_log_message(_("No csv output file to create\n"));
	return nb_phot_stars == 0;
}

static gboolean is_variable_cat(siril_cat_index cat_index) {
	switch (cat_index) {
		case CAT_VARISUM:
		case CAT_VSX:
		case CAT_GCVS:
			return TRUE;
		default:
			return FALSE;
	}
}

static siril_catalogue *get_catstars(struct compstars_arg *args, siril_cat_index cat_index, double *fov_radius_deg) {
	gboolean is_variable = is_variable_cat(cat_index);

	double resolution = get_wcs_image_resolution(args->fit);
	uint64_t sqr_radius = 0;
	double radius = 0.0;

	if (args->narrow_fov) {
		// Limited to the image smallest dimension, to avoid the corners with their potential vignettage
		radius = resolution * min(args->fit->rx, args->fit->ry) / 2.0;	// in degrees
	} else {
		// The whole Field of View
		sqr_radius = (uint64_t)args->fit->rx * args->fit->rx + (uint64_t)args->fit->ry * args->fit->ry / 4;
		radius = resolution * sqrt((double)sqr_radius);	// in degrees
	}
	if (fov_radius_deg)
		*fov_radius_deg = radius;

	// preparing the query
	siril_catalogue *siril_cat = siril_catalog_fill_from_fit(args->fit, cat_index, max(args->target_star->mag + 6.0, 17.0));
	siril_cat->radius = radius * 60.; // overwriting to account for narrow argument
	siril_cat->phot = !is_variable;

	// and retrieving its results
	if (siril_catalog_conesearch(siril_cat) <= 0) {// returns the nb of stars
		siril_log_color_message(_("No stars retrieved from the catalog %s\n"), "red", catalog_to_str(cat_index));
		siril_catalog_free(siril_cat);
		return NULL;
	}

	// projecting to check if stars are within the image
	if (!siril_catalog_project_with_WCS(siril_cat, args->fit, FALSE, FALSE) && siril_cat->nbincluded > 0) {
		siril_log_message(_("-> %i %s stars found within the image from %s\n"),
		siril_cat->nbincluded, (is_variable) ? _("variable") : _("comparison"), catalog_to_str(cat_index));
		return siril_cat;
	} else {
		siril_catalog_free(siril_cat);
		return NULL;
	}
}

gpointer compstars_worker(gpointer p) {
	int retval;
	siril_log_color_message(_("Comparison stars: processing...\n"), "green");
	struct compstars_arg *args = (struct compstars_arg *) p;
	sky_object_query_args *query_args = NULL;
	//0. check pre-requisites
	if (!has_wcs(args->fit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		retval = 1;
		goto end;
	}
	g_assert(args->cat == CAT_APASS || args->cat == CAT_NOMAD);

	// 1. search for the variable star
	query_args = init_sky_object_query(); // for the reference star
	query_args->fit = args->fit;
	query_args->name = g_strdup(args->target_name);
	query_args->server = QUERY_SERVER_SIMBAD_PHOTO;
	retval = cached_object_lookup(query_args);
	if (retval)
		goto end;
	if (query_args->item->mag == 0.0 || query_args->item->bmag == 0.0) {
		siril_log_color_message(_("Target star photometric information not available.\n"), "red");
		retval = 1;
		goto end;
	}
	args->target_star = calloc(1, sizeof(cat_item));
	siril_catalogue_copy_item(query_args->item, args->target_star);
	if (!args->target_star) {
		siril_log_color_message(_("No variable star selected\n"), "salmon");
		retval = 1;
		goto end;
	}
	// 2. get a catalogue of stars for the field
	double radius;
	args->cat_stars = get_catstars(args, args->cat, &radius);
	if (!args->cat_stars) {
		retval = 1;
		siril_log_color_message(_("No comparison stars found in the image, aborting\n"), "red");
		goto end;
	}
	// and check the target star is relatively well centered - warn if not
	double dist = compute_coords_distance(args->cat_stars->center_ra, args->cat_stars->center_dec, args->target_star->ra, args->target_star->dec); // in degrees
	if (dist > 0.2 * radius) {
		siril_log_color_message("Target star is off the center of field of view by %.1f arcmin, photometry results may be impacted\n", "salmon", dist * 60.);
	}

	// 3. get variable stars catalogues if any
	int cat2discard = com.pref.phot_set.discard_var_catalogues;
	if (cat2discard > 0) {
		for (int i = 0; i < MAX_VAR_CAT; i++) {
			if (cat2discard & ( 1 << i )) {
				siril_catalogue *siril_cat_var = get_catstars(args, var_cat[i], NULL);
				if (siril_cat_var)
					args->var_stars_cat = g_list_prepend(args->var_stars_cat, siril_cat_var);
			}
		}
	}

	// 4. extract those matching our photometry need
	retval = sort_compstars(args);

end:
	args->retval = retval;
	args->has_GUI = TRUE;
	if (!siril_add_idle(end_compstars, args)) {
		args->has_GUI = FALSE;
		end_compstars(args);	// we still need to free all
	}
	free_sky_object_query(query_args);
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

