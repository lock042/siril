/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/siril_world_cs.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "core/siril_date.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "algos/annotate.h"
#include "algos/siril_wcs.h"
#include "astrometry_solver.h"

static int parse_buffer(char *buffer) {
	gchar **token, **part, *realname = NULL;
	int nargs;
	gboolean is_solar_system = FALSE;
	SirilWorldCS *world_cs = NULL;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

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
	else if (g_str_has_prefix(buffer, "# Flag:")) { /// case of an asteroid
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

		gchar *objname = NULL;
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
			return FALSE;
		}
		siril_debug_print("Object name: %s\n", objname);
		// Add date at the end
		gchar *formatted_date = date_time_to_date(gfit.date_obs);
		realname = g_strdup_printf("%s\\n(%s)", objname, formatted_date);
		g_free(objname);
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
			}
		}
		g_strfreev(fields);

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
		add_object_in_catalogue(realname, world_cs, is_solar_system);

		g_free(alpha);
		g_free(delta);
		g_free(realname);
		g_strfreev(display_name);
		siril_world_cs_unref(world_cs);
		return 0;
	}
	return 1;
}

void on_search_objects_entry_activate(GtkEntry *entry, gpointer user_data) {
	if (!has_wcs(&gfit)) return;
	// TODO: search in local annotation catalogues first
	const gchar *name = gtk_entry_get_text(GTK_ENTRY(entry));
	gboolean solarsystem = g_str_has_prefix(name, "a:") || g_str_has_prefix(name, "c:") ||
		g_str_has_prefix(name, "p:");

	control_window_switch_to_tab(OUTPUT_LOGS);
	gchar *result = search_in_online_catalogs(name, solarsystem ? QUERY_SERVER_EPHEMCC : QUERY_SERVER_SIMBAD);

	if (result && !parse_buffer(result)) {
		GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
		if (!gtk_toggle_tool_button_get_active(button)) {
			gtk_toggle_tool_button_set_active(button, TRUE);
		} else {
			force_to_refresh_catalogue_list();
			redraw(REDRAW_OVERLAY);
		}
		gtk_entry_set_text(GTK_ENTRY(entry), "");
		siril_close_dialog("search_objects");
	}
	else siril_log_message(_("Object %s not found or encountered an error processing it\n"), name);
}
