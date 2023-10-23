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
#include "algos/search_objects.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "io/annotation_catalogues.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "algos/astrometry_solver.h"
#include "io/remote_catalogues.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/image_display.h"


/* parse response from online catalogue lookups (search_in_online_catalogs()
 * for QUERY_SERVER_EPHEMCC and QUERY_SERVER_SIMBAD_PHOTO) and stores the
 * result in the local annotation catalogues and returns it in the argument if
 * non-NULL */
int parse_catalog_buffer(const gchar *buffer, sky_object_query_args *args) {
	gchar **token, *objname = NULL, *objtype = NULL;
	int nargs;
	gboolean check_for_duplicates = FALSE;
	annotations_cat target_cat = ANCAT_NONE;
	if (!buffer || buffer[0] == '\0')
		return 1;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	switch (args->server) {
	case (QUERY_SERVER_SIMBAD_PHOTO):
		// In case of a bad request...
		if (nargs <= 2 || g_strrstr(buffer, "No known catalog") || g_strrstr(buffer, "incorrect id")) {
			g_strfreev(token);
			args->retval = 1;
			break;
		}
		objname = g_shell_unquote(token[2], NULL);
		gint rank = 3;
		double ra = 0.0, dec = 0.0, pmra = 0.0, pmdec = 0.0;
		double Vmag = 0.0, Bmag = 0.0;
		while (token[rank]) {
			if (g_str_has_prefix(token[rank], "Coordinates(ICRS")) {
				gchar **fields = g_strsplit(token[rank], " ", -1);
				double hours = 0.0, min = 0.0, seconds = 0.0, degres = 0.0;
				// RA
				sscanf(fields[1], "%lf", &hours);
				sscanf(fields[2], "%lf", &min);
				sscanf(fields[3], "%lf", &seconds);
				ra = 360.0*hours/24.0 + 360.0*min/(24.0*60.0) + 360.0*seconds/(24.0*60.0*60.0);
				// DEC
				sscanf(fields[5], "%lf", &degres);
				sscanf(fields[6], "%lf", &min);
				sscanf(fields[7], "%lf", &seconds);
				if (degres < 0.0)
					dec = degres - min/60.0 - seconds/3600.0;
				else dec = degres + min/60.0 + seconds/3600.0;
				g_strfreev(fields);
			}
			// Finally, retrieve the B and V magnitudes
			else if (g_str_has_prefix(token[rank], "Flux B")){
				gchar **fields = g_strsplit(token[rank], " ", -1);
				sscanf(fields[3], "%lf", &Bmag);
				g_strfreev(fields);
			}
			else if (g_str_has_prefix(token[rank], "Flux V")){
				gchar **fields = g_strsplit(token[rank], " ", -1);
				sscanf(fields[3], "%lf", &Vmag);
				g_strfreev(fields);
			}
			else if (g_str_has_prefix(token[rank], "Proper motions")){
				gchar **fields = g_strsplit(token[rank], " ", -1);
				sscanf(fields[2], "%lf", &pmra);
				sscanf(fields[3], "%lf", &pmdec);
				g_strfreev(fields);
			}
			rank++;
		}
		cat_item *item = calloc(1, sizeof(cat_item));
		item->mag = Vmag;
		item->bmag = Bmag;
		item->ra = ra;
		item->dec = dec;
		item->pmra = pmra;
		item->pmdec = pmdec;
		item->name = g_strdup(objname);
		args->item = item;
		args->retval = 0;
		check_for_duplicates = TRUE;
		target_cat = USER_DSO_CAT_INDEX;
		break;
	// This is for a Miriade emphemcc search (SSO)
	case QUERY_SERVER_EPHEMCC:
		if (!g_str_has_prefix(buffer, "# Flag: 1")) {
			siril_log_color_message(_("IMCCE server returned:\n"), "red");
			int max_lines = min(3, nargs); // we display 3 lines max
			for (int i = 0; i < max_lines; i++)
				siril_log_color_message(_("%s\n"), "red", token[i]);
			g_strfreev(token);
			args->retval = 1;
			return args->retval;
		}
		// Extract the first word of the third line, surrounded by (":" + space) and ("space + "|")
		// "# Comet : PANSTARRS  (C/2017 T2) | ..."
		// "# Comet : P/Churyumov-Gerasimenko  (67P) | ..."
		// "# Asteroid: 103516 2000 BY4 | ..."
		// "# Planet : 4 Mars | ..."
		gchar *pattern = "#(.+?):(.+?)\\|";
		gchar **parts = g_regex_split_simple(pattern, token[2], G_REGEX_DEFAULT, 0);
		if (g_strv_length(parts) == 4) {
			objtype = g_strdup(g_strstrip(parts[1]));
			objname = g_strdup(g_strstrip(parts[2]));
		}
		if (!objname) {
			siril_log_message(_("Unsupported object type in response: %s\n"), token[2]);
			g_strfreev(token);
			g_strfreev(parts);
			g_free(objname);
			g_free(objtype);
			args->retval = 1;
			return 1;
		}
		siril_debug_print("Object name: %s\n", objname);
		// Then, retrieve the data
		gchar **fields = g_strsplit(token[4], ",", -1);
		guint n = g_strv_length(fields);
		if (n == 16) {
			double ra, dec, vra, vdec, mag;
			ra = parse_hms(fields[2]);
			dec = parse_dms(fields[3]);
			vra = g_strtod(fields[13], NULL) * 60.; // vra stored in arcsec/hr but given in arcsec/min
			vdec = g_strtod(fields[14], NULL) * 60.; // vdec stored in arcsec/hr but given in arcsec/min
			mag = g_strtod(fields[9], NULL);
			// and preparing to store
			cat_item *item = calloc(1, sizeof(cat_item));
			item->ra = ra;
			item->dec = dec;
			item->name = g_strdup(objname); // we use the offical name returned by the server
			item->alias = g_strdup(args->name); // we store the name that was queried
			item->sitelon = args->fit->sitelong;
			item->sitelat = args->fit->sitelat;
			item->siteelev = args->fit->siteelev;
			item->dateobs = date_time_to_Julian(args->fit->date_obs);
			item->type = g_strdup(objtype);
			item->vra = vra;
			item->vdec = vdec;
			item->mag = mag;
			args->item = item;
			check_for_duplicates = TRUE;
			target_cat = USER_SSO_CAT_INDEX;
			args->retval = 0;
		} else {
			siril_log_message(_("Could not parse the server response: %s\n"), token[3]);
			g_strfreev(token);
			g_strfreev(parts);
			g_free(objname);
			g_free(objtype);
			args->retval = 1;
			return 1;
		}
		g_strfreev(parts);
		g_free(objname);
		g_free(objtype);
		g_strfreev(fields);
	default:
		siril_debug_print("unknown query catalogue\n");
	}
	g_strfreev(token);

	if (!args->retval && args->item) {
		SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(args->item->ra, args->item->dec);
		gchar *alpha = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02.2lfs");
		gchar *delta = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02.2lf\"");
		siril_log_message(_("Found %s at coordinates: %s, %s\n"), args->item->name, alpha, delta);
		g_free(alpha);
		g_free(delta);
		siril_world_cs_unref(world_cs);
		if (target_cat != ANCAT_NONE) {
			add_item_in_catalogue(args->item, target_cat, check_for_duplicates);
			com.pref.gui.catalog[target_cat] = TRUE;	// and display it
		}
		return 0;
	}
	return 1;
}

/* this should be the new entry point for single object search from their name */
int cached_object_lookup(sky_object_query_args *args) {
	gboolean solarsystem = args->prefix != NULL;
	args->retval = 1;
	if (!solarsystem) {
		// local search first, if useful
		object_catalog cattype;
		args->item = search_in_annotations_by_name(args->name, &cattype);
		if (args->item) {
			SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(args->item->ra, args->item->dec);
			gchar *alpha = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02.2lfs");
			gchar *delta = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02.2lf\"");
			siril_log_message(_("Object %s (exact match: %s) was found in the %s catalogue at: %s, %s\n"),
					args->name, args->item->name, catalog_to_str(cattype),
					alpha, delta);
			g_free(alpha);
			g_free(delta);
			siril_world_cs_unref(world_cs);
		}
	} else {
		args->item = search_in_solar_annotations(args);
	}
	if (!args->item) {
		gchar *result = search_in_online_catalogs(args);
		if (result) {
			parse_catalog_buffer(result, args);
			free_fetch_result(result);
		} else {
			args->retval = 1;
		}
	} else {
		args->retval = 0;
	}
	if (!args->retval) {
		if (!is_inside(args->fit, args->item->ra, args->item->dec)) {
			siril_log_message(_("Object %s record was found but is not within the bounds of the image\n"), args->name);
			siril_catalog_free_item(args->item);
			args->item = NULL;
			args->retval = 1;
		}
	} else {
		siril_log_message(_("Object %s not found or encountered an error processing it\n"), args->name);
	}
	return args->retval;
}

void on_search_objects_entry_activate(GtkEntry *entry, gpointer user_data) {
	if (!has_wcs(&gfit))
		return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	sky_object_query_args *args = init_sky_object_query();
	args->name = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));
	args->fit = &gfit;
	start_in_new_thread(catsearch_worker, args);
}


/********** object lookup for astrometry, can we merge both? ************/

struct sky_object platedObject[RESOLVER_NUMBER];

/* instead of populating platedObject from a resolver buffer, use the passed object */
void add_plated_from_annotations(const cat_item *obj) {
	double ra = obj->ra;
	double dec = obj->dec;
	int resolver = RESOLVER_LOCAL;
	platedObject[resolver].world_cs = siril_world_cs_new_from_a_d(ra, dec);
	platedObject[resolver].imageCenter.x = ra;
	platedObject[resolver].imageCenter.y = dec;
	platedObject[resolver].south = (dec < 0.0);
	platedObject[resolver].name = g_strdup(obj->name);
}

/* parse the result from search_in_online_catalogs() for QUERY_SERVER_CDS,
 * QUERY_SERVER_VIZIER, QUERY_SERVER_SIMBAD, for object name to coordinates
 * conversion (used only by the astrometry GUI). See also parse_catalog_buffer().
 */
int parse_resolver_buffer(const char *buffer, struct sky_object *obj) {
	gchar **token, **fields;
	point center;
	int nargs, i = 0;
	resolver_t resolver = RESOLVER_UNSET;
	gboolean SIMBAD_alternative = FALSE;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	while (i < nargs) {
		if (g_strrstr (token[i], "=NED")) {
			resolver = RESOLVER_NED;
		} else if (g_strrstr (token[i], "=Simbad")) {
			resolver = RESOLVER_SIMBAD;
		} else if (g_str_has_prefix (token[i], "oid")) {
			resolver = RESOLVER_SIMBAD;
			SIMBAD_alternative = TRUE;
		} else if (g_strrstr(token[i], "=VizieR")) {
			resolver = RESOLVER_VIZIER;
		} else if (g_str_has_prefix (token[i], "%J ")) {
			fields = g_strsplit(token[i], " ", -1);
			sscanf(fields[1], "%lf", &center.x);
			sscanf(fields[2], "%lf", &center.y);
			if (resolver != RESOLVER_UNSET) {
				platedObject[resolver].world_cs = siril_world_cs_new_from_a_d(center.x, center.y);

				/* others */
				platedObject[resolver].imageCenter = center;
				platedObject[resolver].south = (center.y < 0.0);
			}
			g_strfreev(fields);
		} else if (g_str_has_prefix (token[i], "%I.0 ")) {
			if (resolver != RESOLVER_UNSET) {
				gchar *name = g_strstr_len(token[i], strlen(token[i]), "%I.0 ");
				gchar *realname;
				realname = g_strdup(name + 5);
				platedObject[resolver].name = realname;
			}
		} else if (g_str_has_prefix (token[i], "%I NAME ")) {
			if (resolver != RESOLVER_UNSET) {
				gchar *name = g_strstr_len(token[i], strlen(token[i]), "%I NAME ");
				gchar *realname;
				realname = g_strdup(name + 5 + 3);
				g_free(platedObject[resolver].name);
				platedObject[resolver].name = realname;
			}
		} else if (SIMBAD_alternative) {
			fields = g_strsplit(token[i], "\t", -1);
			guint n = g_strv_length(token);
			if (n > 2 && resolver != RESOLVER_UNSET) {
				sscanf(fields[1], "%lf", &center.x);
				sscanf(fields[2], "%lf", &center.y);
				gchar *realname = g_shell_unquote(fields[3], NULL);
				g_free(platedObject[resolver].name);
				platedObject[resolver].name = realname;
				platedObject[resolver].world_cs = siril_world_cs_new_from_a_d(center.x, center.y);
				platedObject[resolver].imageCenter = center;
				platedObject[resolver].south = (center.y < 0.0);
				// don't come back
				SIMBAD_alternative = FALSE;
			}
			g_strfreev(fields);
		}

		i++;
	}
	g_strfreev(token);
	return 0;
}

void free_Platedobject() {
	for (int i = 0; i < RESOLVER_NUMBER; i++) {
		if (platedObject[i].name) {
			siril_world_cs_unref(platedObject[i].world_cs);
			g_free(platedObject[i].name);
			platedObject[i].name = NULL;
		}
	}
}

gboolean has_nonzero_coords() {
	for (int i = 0; i < RESOLVER_NUMBER; i++){
		if (fabs(platedObject[i].imageCenter.x) > 0.000001) return TRUE;
		if (fabs(platedObject[i].imageCenter.y) > 0.000001) return TRUE;
	}
	return FALSE;
}

// returns a string describing the site coordinates on Earth in a format suited for queries
static gchar *retrieve_site_coord(fits *fit) {
	// if lat and lon are null, we still format them so as to have consistent outputs to parse
	GString *formatted_site_coord = g_string_new("");
	g_string_printf(formatted_site_coord, "+%f,+%f,%f", fit->sitelat, fit->sitelong, fit->siteelev);

	// if (fit->sitelat < 0.)
	// 	formatted_site_coord = g_string_append(formatted_site_coord, "%2D");
	// else formatted_site_coord = g_string_append(formatted_site_coord, "%2B");
	// g_string_append_printf(formatted_site_coord, "%f,", fabs(fit->sitelat));

	// if (fit->sitelong < 0.)
	// 	formatted_site_coord = g_string_append(formatted_site_coord, "%2D");
	// else formatted_site_coord = g_string_append(formatted_site_coord, "%2B");
	// g_string_append_printf(formatted_site_coord, "%f,", fabs(fit->sitelong));

	// g_string_append_printf(formatted_site_coord, "%f", fabs(fit->siteelev));

	return g_string_free(formatted_site_coord, FALSE);
}

// free the result with free_fetch_result
gchar *search_in_online_catalogs(sky_object_query_args *args) {
#ifndef HAVE_NETWORKING
	siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
	return NULL;
#else
	GString *string_url = NULL;
	gchar *name = args->name;
	switch(args->server) {
	case QUERY_SERVER_CDS:
		string_url = g_string_new(name);
		g_string_replace(string_url, "+", "%2B", 0);
		g_string_replace(string_url, "-", "%2D", 0);
		string_url = g_string_prepend(string_url, "/-oI/A?");
		string_url = g_string_prepend(string_url, CDSSESAME);
		siril_log_message(_("Searching %s in CDSESAME...\n"), name);
		break;
	case QUERY_SERVER_VIZIER:
		string_url = g_string_new(name);
		g_string_replace(string_url, "+", "%2B", 0);
		g_string_replace(string_url, "-", "%2D", 0);
		string_url = g_string_prepend(string_url, "/-oI/A?");
		string_url = g_string_prepend(string_url, VIZIERSESAME);
		siril_log_message(_("Searching %s in VIZIER...\n"), name);
		break;
	default:
	case QUERY_SERVER_SIMBAD:
		string_url = g_string_new(name);
		g_string_replace(string_url, "+", "%2B", 0);
		g_string_replace(string_url, "-", "%2D", 0);
		string_url = g_string_prepend(string_url, SIMBADSESAME);
		string_url = g_string_append(string_url, "';");
		siril_log_message(_("Searching %s in SIMBAD...\n"), name);
		break;
	case QUERY_SERVER_EPHEMCC:
		// see https://ssp.imcce.fr/webservices/miriade/api/ephemcc/
		if (!args->fit->date_obs) {
			siril_log_color_message(_("This command only works on images that have observation date information\n"), "red");
			g_string_free(string_url, TRUE);
			return NULL;
		}
		string_url = g_string_new(EPHEMCC);
		g_string_append_printf(string_url, "&-name=%s:%s", args->prefix, name);
		gchar *formatted_date = date_time_to_FITS_date(args->fit->date_obs);
		g_string_append_printf(string_url, "&-ep=%s", formatted_date);
		gchar *formatted_site = retrieve_site_coord(args->fit);
		g_string_append_printf(string_url, "&-observer=%s", formatted_site);
		siril_log_message(_("Searching for solar system object %s on observation date %s\n"),
				name, formatted_date);
		if (args->fit->sitelat == 0.0 && args->fit->sitelong == 0.0) {
			siril_log_color_message(_("No topocentric data available. Set to geocentric\n"), "salmon");
		} else {
			siril_log_message(_("at lat: %f, long: %f, alt: %f\n"), args->fit->sitelat,
				args->fit->sitelong, args->fit->siteelev);
		}
		g_free(formatted_site);
		g_free(formatted_date);
		break;
	case QUERY_SERVER_SIMBAD_PHOTO:
		// SIMBAD request to get the magnitudes (BVRIJ) in addition to coordinates for a particular star/object
		string_url = g_string_new(name);
		g_string_replace(string_url, "+", "%2B", 0);
		g_string_replace(string_url, "-", "%2D", 0);
		string_url = g_string_prepend(string_url, SIMBAD);
		siril_log_message(_("Searching %s in SIMBAD\n"), name);
		break;
	}

	gchar *url = g_string_free(string_url, FALSE);
	gchar *cleaned_url = url_cleanup(url);
	g_free(url);
	siril_debug_print("URL: %s\n", cleaned_url);

	char *result = fetch_url(cleaned_url);

	g_free(cleaned_url);
	return result;
#endif
}

gpointer catsearch_worker(gpointer p) {
	sky_object_query_args *args = (sky_object_query_args *)p;
	if(!args)
		return GINT_TO_POINTER(1);

	gchar **parts = g_strsplit(args->name, ":", -1);
	if (g_strv_length(parts) == 2) {
		g_free(args->name);
		args->name = g_strdup(parts[1]);
		args->prefix = g_strdup(parts[0]);
		args->server = QUERY_SERVER_EPHEMCC;
	} else {
		args->server = QUERY_SERVER_SIMBAD_PHOTO;
	}
	gboolean found_it = !cached_object_lookup(args);

	if (!com.script) { // we need to update GUI
		siril_add_idle(end_process_catsearch, args); // this will free the args
	} else {
		free_sky_object_query(args);
	}

	return GINT_TO_POINTER(!found_it);
}
