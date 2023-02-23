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


#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "algos/annotate.h"
#include "algos/siril_wcs.h"
#include "algos/search_objects.h"
#include "astrometry_solver.h"

/* parse results from online catalogue searches and store them in the local annotation catalogues */
int parse_buffer(const gchar *buffer, double lim_mag) {
	gchar **token, **part, *realname = NULL;
	int nargs;
	gboolean is_solar_system = FALSE;
	gboolean is_in_field = FALSE;
	SirilWorldCS *world_cs = NULL;
	gchar *objname = NULL;
	if (!buffer || buffer[0] == '\0')
		return 1;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	/* The goal here is to determine if a cone search has been performed or
	 * not, by searching the (valid) answer for paricular "word" */
	gboolean conesearch = FALSE;
	if (g_str_has_prefix(buffer, "# Flag:") && nargs >= 3) {
		// Header of the answer to be tested
		if (g_str_has_prefix(buffer, "# Flag: 0")) {
			siril_log_message("No Solar System object found in this FOV \n");
			g_strfreev(token);
			return 1;
		}

		// Is the answer valid for conesearch?
		// Extract the first word of the third line, surrounded by ("#" + space) and ("space + "|")
		// We are looking for "Num", which is typical for a consearch request
		char *start = strchr(token[2], '#');
		if (start && !strncmp(start, "# Num |", 7))
			conesearch = TRUE;
	}

	// This is for a request to SIMBAD but in photometry style
	if (g_str_has_prefix(buffer, "B	V	R	I	J")) {
		gchar **fields = g_strsplit(token[1], "\t", -1);
		guint n = g_strv_length(fields);
		if (n > 4) {
			double Bmag, Vmag; //	Rmag, Imag and Jmag are not used but are available, in case...
			if (sscanf(fields[0], "%lf", &Bmag) && sscanf(fields[1], "%lf", &Vmag)) {
				//For now, do nothing. Bmag and Vmag are aimed at beeing saved in a cstars, as target data. See the next MR...
				siril_log_message(_("Bmag= %lf  , Vmag= %lf\n"), Bmag, Vmag);
			}
		}
		g_strfreev(fields);
		g_strfreev(token);
		return 1;
	}

	if (g_str_has_prefix(buffer, "oid")) {
		gchar **fields = g_strsplit(token[1], "\t", -1);
		guint n = g_strv_length(fields);
		if (n > 3) {
			double ra, dec;
			if (sscanf(fields[1], "%lf", &ra) && sscanf(fields[2], "%lf", &dec)) {
				world_cs = siril_world_cs_new_from_a_d(ra, dec);
				realname = g_shell_unquote(fields[3], NULL);
			}
		}
		g_strfreev(fields);
	}

	/// Case of a Solar System Object in Search mode
	else if (g_str_has_prefix(buffer, "# Flag:") && !conesearch) {
		// First, retrieve the name of the object
		// example strings:
		// "# Comet : PANSTARRS  (C/2017 T2) | ..."
		// "# Comet : P/Churyumov-Gerasimenko  (67P) | ..."
		// "# Asteroid: 103516 2000 BY4 | ..."
		// "# Planet : 4 Mars | ..."

		// The result of the request is different between geocentric and topocentric
		// To retrieve the right RA and DEC, we have to introcuce a +1 offset in the position (see RA and DEC retrieve, below)
		int offset_ind = 0;
		if (gfit.sitelat != 0.) {
			offset_ind = 1;
		}

		// Is the found object valid?
		// Extract the first word of the third line, surrounded by (":" + space) and ("space + "|")
		char *start = strchr(token[2], ':');
		if (start && *(start+1) != '\0') {
			start += 2; // remove the space
			char *end = strchr(start, '|');
			if (end) {
				end -= 2; // remove the space
				objname = g_strndup(start, end-start+1);
			}
		}
		if (!objname) {
			siril_log_message(_("Unsupported object type in response: %s\n"), token[2]);
			g_strfreev(token);
			g_free(objname);
			return 1;
		}
		siril_debug_print("Object name: %s\n", objname);

		// Add date at the end
		gchar *formatted_date = date_time_to_date(gfit.date_obs);
		realname = g_strdup_printf("%s\\n(%s)", objname, formatted_date);
		g_free(formatted_date);

		// Then, retrieve the coordinates
		gchar **fields = g_strsplit(token[4], ",", -1);
		guint n = g_strv_length(fields);
		if (n > 2) {
			double hours = 0.0, min = 0.0, seconds = 0.0, degres = 0.0;
			double ra, dec;
			gboolean valid = FALSE;

			// For RA coordinates
			part = g_strsplit(fields[1 + offset_ind], " ", -1);
			if (g_strv_length(part) > 3) {
				sscanf(part[1], "%lf", &hours);
				sscanf(part[2], "%lf", &min);
				sscanf(part[3], "%lf", &seconds);
				ra = 360.0*hours/24.0 + 360.0*min/(24.0*60.0) + 360.0*seconds/(24.0*60.0*60.0);
				valid = TRUE;
			}
			g_strfreev(part);

			// For Dec coordinates
			if (valid) {
				part = g_strsplit(fields[2 + offset_ind], " ", -1);
				if (g_strv_length(part) > 3) {
					sscanf(part[1], "%lf", &degres);
					sscanf(part[2], "%lf", &min);
					sscanf(part[3], "%lf", &seconds);
					if (degres < 0.0)
						dec = degres - min/60.0 - seconds/3600.0;
					else dec = degres + min/60.0 + seconds/3600.0;
				}
				else valid = FALSE;
				g_strfreev(part);
			}

			if (valid) {
				world_cs = siril_world_cs_new_from_a_d(ra, dec);
				is_solar_system = TRUE;
				is_in_field = FALSE;
			}
		}
		g_free(objname);
		objname = NULL;
		g_strfreev(fields);
		fields = NULL;
	}

	/// Case of Solar System Objects in Conesearch mode
	else if (g_str_has_prefix(buffer, "# Flag:") && conesearch) {
		gint rank = 3;
		guint nbr_a = 0;
		while (token[rank]) {
			char *start = strchr(token[rank], '|');
			if (start && *(start+1) != '\0') {
				start += 2; // remove the space
				char *end = strchr(start, '|');
				if (end) {
					end -= 2; // remove the space
					g_free(objname); // in case it was set in a previous run through the loop
					objname = g_strndup(start, end-start+1);
				}
			}

			if (!objname) {
				siril_log_message("Unsupported object type in response: \n");
				g_strfreev(token);
				g_free(objname);
				return 1;
			}
			siril_debug_print("Object name: %s\n", objname);

			// Add date at the end
			gchar *formatted_date = date_time_to_date(gfit.date_obs);
			realname = g_strdup_printf("%s", objname);
			g_free(formatted_date);

			// Then, retrieve the coordinates
			gchar **fields = g_strsplit(token[rank], "|", -1);
			guint n = g_strv_length(fields);
			double mag = 0.0;
			double hours = 0.0, min = 0.0, seconds = 0.0, degres = 0.0;
			double ra = 0.0, dec =0;
			gboolean valid = FALSE;
			if (n > 2) {
				// For RA coordinates
				part = g_strsplit(fields[2], " ", -1);
				if (g_strv_length(part) > 3) {
					sscanf(part[1], "%lf", &hours);
					sscanf(part[2], "%lf", &min);
					sscanf(part[3], "%lf", &seconds);
					ra = 360.0*hours/24.0 + 360.0*min/(24.0*60.0) + 360.0*seconds/(24.0*60.0*60.0);
					valid = TRUE;
				}
				g_strfreev(part);

				// For Dec coordinates
				if (valid) {
					part = g_strsplit(fields[3], " ", -1);
					if (g_strv_length(part) > 3) {
						sscanf(part[1], "%lf", &degres);
						sscanf(part[2], "%lf", &min);
						sscanf(part[3], "%lf", &seconds);
						if (degres < 0.0)
							dec = degres - min/60.0 - seconds/3600.0;
						else dec = degres + min/60.0 + seconds/3600.0;
					}
					else valid = FALSE;
					g_strfreev(part);
				}

				// Then, retrieve the magnitude, used to sort the objects to be displayed
				if (valid) {
					if (fields[5]) sscanf(fields[5], "%lf", &mag);	// Get the magnitude
					else valid = FALSE;

					// Is mag enought to be shown? Intentionaly hard-limited to 20
					if (mag < lim_mag && is_inside(&gfit, ra, dec)) nbr_a+=1;
					else valid = FALSE;
				}
				if (valid) {
					world_cs = siril_world_cs_new_from_a_d(ra, dec);
					is_solar_system = TRUE;
					is_in_field = TRUE;
				}
			}
			g_strfreev(fields);

			if (world_cs && realname && valid) {
				gchar **display_name = g_strsplit(realname, "\\n", 2);
				gchar* alpha = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02ds");
				gchar* delta = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
				siril_log_message(_("Found %s at coordinates: %s, %s with magnitude: %0.1lf\n"),display_name[0],
						alpha, delta, mag);
				g_free(alpha);
				g_free(delta);

				if (!com.script && !com.headless) {	// Write in catalogue only if in GUI mode, not in script-mode or headless-mode
					com.pref.gui.catalog[6] = TRUE;	// enabling the user catalog in which it will be added
					add_object_in_catalogue(realname, world_cs, is_solar_system, is_in_field);
					siril_world_cs_unref(world_cs);
				}
				g_strfreev(display_name);
			}
			g_free(realname);
			rank+=1;
		}
		g_free(objname);
		g_strfreev(token);
		siril_log_message(_("Found %d Solar System Objects with mag < %0.1lf\n"), nbr_a, lim_mag);
		return 0;

	} else {
		// what case is that?
		int i = 0;
		while (i < nargs) {
			if (g_str_has_prefix(token[i], "%I.0 ") && !realname) {
				gchar *name = g_strstr_len(token[i], strlen(token[i]), "%I.0 ");
				realname = g_strdup(name + 5);

			} else if (g_str_has_prefix (token[i], "%J ") && !world_cs) {
				gchar **fields = g_strsplit(token[i], " ", -1);
				guint n = g_strv_length(fields);
				if (n > 2) {
					double ra, dec;
					if (sscanf(fields[1], "%lf", &ra) && sscanf(fields[2], "%lf", &dec))
						world_cs = siril_world_cs_new_from_a_d(ra, dec);
				}
				g_strfreev(fields);
			}
			i++;
		}
	}
	g_strfreev(token);

	if (world_cs && realname) {
		gchar **display_name = g_strsplit(realname, "\\n", 2);
		gchar *alpha = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02ds");
		gchar *delta = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
		siril_log_message(_("Found %s at coordinates: %s, %s\n"),display_name[0], alpha, delta);
		com.pref.gui.catalog[6] = TRUE;	// enabling the user catalog in which it will be added
		add_object_in_catalogue(realname, world_cs, is_solar_system, is_in_field);

		g_free(alpha);
		g_free(delta);
		g_free(realname);
		g_strfreev(display_name);
		siril_world_cs_unref(world_cs);
		return 0;
	}
	if (world_cs)
		siril_world_cs_unref(world_cs);
	return 1;
}

gchar *search_object(const gchar *name) {
	gboolean solarsystem = g_str_has_prefix(name, "a:") ||
		g_str_has_prefix(name, "c:") ||	g_str_has_prefix(name, "p:");
	if (!solarsystem) {
		gchar *result = search_in_online_catalogs(name, QUERY_SERVER_SIMBAD_PHOTO);
		parse_buffer(result, 20.0);
		// TODO: avoid doing the second request
		return search_in_online_catalogs(name, QUERY_SERVER_SIMBAD);
	}
	return search_in_online_catalogs(name, QUERY_SERVER_EPHEMCC);
}

void on_search_objects_entry_activate(GtkEntry *entry, gpointer user_data) {
	if (!has_wcs(&gfit)) return;
	// TODO: search in local annotation catalogues first
	const gchar *name = gtk_entry_get_text(GTK_ENTRY(entry));

	control_window_switch_to_tab(OUTPUT_LOGS);

	gchar *result = search_object(name);

	if (result && !parse_buffer(result, 20.0)) {
		GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
		if (!gtk_toggle_tool_button_get_active(button)) {
			gtk_toggle_tool_button_set_active(button, TRUE);
		} else {
			refresh_found_objects();
			redraw(REDRAW_OVERLAY);
		}
		gtk_entry_set_text(GTK_ENTRY(entry), "");
		siril_close_dialog("search_objects");
	}
	else siril_log_message(_("Object %s not found or encountered an error processing it\n"), name);
}
