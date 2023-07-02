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
#include "algos/annotate.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "io/remote_catalogues.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/image_display.h"

#define MAX_DATA_ITEMS 4

/* parse response from online catalogue lookups (search_in_online_catalogs()
 * for QUERY_SERVER_EPHEMCC and QUERY_SERVER_SIMBAD_PHOTO) and stores the
 * result in the local annotation catalogues and returns it in the argument if
 * non-NULL */
int parse_catalog_buffer(const gchar *buffer, psf_star **result) {
	gchar **token, **part, *realname = NULL;
	int nargs;
	gboolean check_for_duplicates = FALSE;
	annotations_cat target_cat = ANCAT_NONE;
	SirilWorldCS *world_cs = NULL;
	gchar *objname = NULL;
	if (!buffer || buffer[0] == '\0')
		return 1;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	// This is for a request to SIMBAD (DSO)
	if (g_str_has_prefix(buffer, "C.D.S.") && nargs > 2) {
		// In case of a sad request...
		if (g_strrstr(buffer, "No known catalog") || g_strrstr(buffer, "incorrect id")) return 1;
		// This case does not end up being stored in a annotation catalogue, it's
		// saved in an alternate structure that contains photometry and returns
		realname = g_shell_unquote(token[2], NULL);
		gint rank = 3;
		double ra = 0.0, dec = 0.0;
		double Vmag = 0.0, Bmag = 0.0;
		while (token[rank]) {
			if (g_str_has_prefix(token[rank], "Coordinates(ICRS")){
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

			rank++;
		}
		world_cs = siril_world_cs_new_from_a_d(ra, dec);

		if (result && (Vmag == 0.0 || Bmag == 0.0))
			siril_log_color_message(_("No photometric data available for this object.\n"), "red");
		psf_star *star = new_psf_star();
		star->mag = Vmag;
		star->Bmag = Bmag;
		star->BV = Bmag - Vmag;
		star->phot = NULL;
		star->ra = ra;
		star->dec = dec;
		star->star_name = g_strndup(realname, 20);

		if (!is_inside(&gfit, ra, dec))
			siril_log_color_message(_("Object not in the field of view...\n"), "salmon");
		else if (result)
			*result = star;

		check_for_duplicates = TRUE;
		target_cat = USER_DSO_CAT_INDEX;
	}

	// This is for a Miriade emphemcc search (SSO)
	else if (g_str_has_prefix(buffer, "# Flag:")) {
		// The result of the request is different between geocentric and topocentric
		// To retrieve the right RA and DEC, we have to introcuce a +1
		// offset in the position (see RA and DEC retrieve, below)
		int offset_ind = 0;
		// TODO: parse this information from the result maybe?
		if (gfit.sitelat != 0.0 && gfit.sitelong != 0.0)
			offset_ind = 1;

		// Extract the first word of the third line, surrounded by (":" + space) and ("space + "|")
		// "# Comet : PANSTARRS  (C/2017 T2) | ..."
		// "# Comet : P/Churyumov-Gerasimenko  (67P) | ..."
		// "# Asteroid: 103516 2000 BY4 | ..."
		// "# Planet : 4 Mars | ..."
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

		// Add date at the end of the name
		gchar *formatted_date = date_time_to_date(gfit.date_obs);
		realname = g_strdup_printf("%s\\n(%s)", objname, formatted_date);
		g_free(formatted_date);

		// Then, retrieve the coordinates
		gchar **fields = g_strsplit(token[4], ",", -1);
		guint n = g_strv_length(fields);
		if (n > 2 + offset_ind) {
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
				check_for_duplicates = FALSE;
				target_cat = USER_SSO_CAT_INDEX;
			}
		}
		g_free(objname);
		objname = NULL;
		g_strfreev(fields);
		fields = NULL;
	}
	else if (g_str_has_prefix(buffer, "!! ")) {
		// error response, for example:
		// "!! 'kelt16' this identifier has an incorrect format for catalog:\n\tkelt : Kilodegree Extremely Little Telescope\n\n query string: kelt16\n"
		siril_log_color_message(_("Failure in name resolution: %s\n"), "red", token[0]+3);
		if (nargs > 1)
			siril_log_color_message("%s\n", "red", token[1]);
	} else {
		// unsupported response format, should not happen
		siril_log_color_message(_("Error in name resolution: %s\n"), "red", token[0]);
	}
	g_strfreev(token);

	if (world_cs && realname) {
		gchar **display_name = g_strsplit(realname, "\\n", 2);
		gchar *alpha = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02.2lfs");
		gchar *delta = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02.2lf\"");
		siril_log_message(_("Found %s at coordinates: %s, %s\n"), display_name[0], alpha, delta);
		g_strfreev(display_name);
		g_free(alpha);
		g_free(delta);

		if (target_cat != ANCAT_NONE) {
			add_object_in_catalogue(realname, world_cs, check_for_duplicates, target_cat);
			com.pref.gui.catalog[target_cat] = TRUE;	// and display it
		}
		g_free(realname);
		siril_world_cs_unref(world_cs);
		return 0;
	}
	g_free(realname);
	if (world_cs)
		siril_world_cs_unref(world_cs);
	return 1;
}

/////////////////////////////////////////////////////////
// Defines the possible data contained in a Vizier catalogue
/////////////////////////////////////////////////////////

typedef struct {
	int CAT_TYP;		// Catalogue type, or number
	int ra_len;			// Lenght of the RA field in the returned stream
	int ra_ind;			// Index of the RA field in the returned stream
	int dec_len;		// Lenght of the DEC field in the returned stream
	int dec_ind;		// Index of the DEC field in the returned stream
	int name_len;		// Lenght of the NAME field in the returned stream
	int name_ind;		// Index of the NAME field in the returned stream
	int Vmag_len;		// Lenght of the Vmag field in the returned stream
	int Vmag_ind;		// Index of the Vmag field in the returned stream
	int Bmag_len;		// Lenght of the Bmag field in the returned stream
	int Bmag_ind;		// Index of the Bmag field in the returned stream
	int eVmag_len;		// Lenght of the eVmag field in the returned stream
	int eVmag_ind;		// Index of the eVmag field in the returned stream
	int eBmag_len;		// Lenght of the eBmag field in the returned stream
	int eBmag_ind;		// Index of the eBmag field in the returned stream
	int fmt;			// Data type for RA and DEC: 0 stands for Alphanumeric, 1 stands for Double
} catalog_plan;

static catalog_plan catalog_data = {0};


	//// APASS fields: RAJ2000 DEJ2000 Vmag Bmag e_Vmag e_Bmag recno
	//// NOMAD fields: RAJ2000 DEJ2000 Vmag Bmag
	//// GAIA3 fields: RAJ2000 DEJ2000 Gmag BPmag e_Gmag e_BPmag DR3Name
	// 1- For now, only the 4 first data are used, but could be extended to
	//    the next items defined in catalog_plan
	// 2- The keyword order (RA/DEC/Name/mg/...) must be kept as is.
	// 3- The catalogue order MUST match the one defined in the structure 'online_catalog' in "remote_catalogues.h"
	// 4- The keywords are specific to each catalogue and are explained in the Vizier pages
static char *catalog_data_label[][MAX_DATA_ITEMS] = {	
	{"RAJ2000", "DEJ2000", "recno", "VTmag"},				//CAT_TYCHO2
	{"RAJ2000", "DEJ2000", "NOMAD1", "Vmag"},				//CAT_NOMAD
	{"RAJ2000", "DEJ2000", "DR3Name", "Gmag"},				//CAT_GAIADR3
	{"RAJ2000", "DEJ2000", "PPMXL", "imag"},				//CAT_PPMXL
	{"RAJ2000", "DEJ2000", "HD", "Vmag"},					//CAT_BRIGHT_STARS
	{"RAJ2000", "DEJ2000", "Name", NULL},					//CAT_APASS No real interrest as this is the same as CAT_AAVSO with no name fields
	{"RAJ2000", "DEJ2000", "GCVS", "magMax"},				//CAT_GCVS
	{"RAJ2000", "DEJ2000", "Name", "max"},					//CAT_AAVSO_Var
	{"RAJ2000", "DEJ2000", "recno", "Vmag"},				//CAT_AAVSO
	{"RAJ2000", "DEJ2000", "PGC", ""},						//CAT_PGC
	{NULL, NULL, NULL, NULL}								// Dummy
};

// Designed to get the index/lenght (see catalog_plan) of the dataset determined by the keywords (see catalog_data_label)
int catalog_mapping(const gchar *buffer, int cata) {
	gchar **token = g_strsplit(buffer, "\n", -1);
	int nb_lines = g_strv_length(token);
	catalog_data.fmt = 0;
	size_t label_qty = sizeof(catalog_data_label[cata]) / (sizeof(catalog_data_label[cata][0]));

	for (int j =0; j < label_qty; j++) {		// Loop on items contained in the catalogue
		int cumul = 0;		// cumulative index
		int F77_len = 0;	// lenght of each item
		for (int i = 0; i < nb_lines; i++) {		// Loops on the lines
			gchar *line = token[i];

		// Prototype of a significant line:
		// #Column	Item	(F77 format statement)	[useless stuff]
		// Example:
		// #Column	RAJ2000	(F10.6)		[ucd=pos.eq.ra;meta.main]
		// #Column	DEJ2000	(F10.6)		[ucd=pos.eq.dec;meta.main]

			if (g_str_has_prefix(line, "#Column")) {		// A significant line begins with "#Column" in the header
				gchar **sub_line = g_strsplit(line, "\t", -1);

				char *start = strchr(sub_line[2], '(');		// extract lenght
				if (start && *(start+1) != '\0') {
					start += 1; // 
					gchar *tmp = NULL;
					char *end = strchr(start, '.');
					if (end) {
						end -= 1; // 
						tmp = g_strndup(start+1, end-start);
						F77_len = g_ascii_strtoll(tmp, NULL, 10);
					} else {
						end = strchr(start, ')');
						end -= 1; // 
						tmp = g_strndup(start+1, end-start);
						F77_len = g_ascii_strtoll(tmp, NULL, 10);
					}
					g_free(tmp);


					if (!g_strcmp0(sub_line[1], catalog_data_label[cata][j]) && catalog_data_label[cata][j]) {		// one item is found
						switch (j) {
							case 0:		//0 being the index of the 1st data, ie RA
								catalog_data.ra_ind = cumul;
								catalog_data.ra_len = F77_len;
								if (g_str_has_prefix(sub_line[2], "(F")) catalog_data.fmt = 1;	// RA/DEC format, 1: FLOAT, 0: CHAR. Depends on catalogue
								break;
							case 1:		//1 being the index of the 2nd data, ie DEC
								catalog_data.dec_ind = cumul;
								catalog_data.dec_len = F77_len;
								if (g_str_has_prefix(sub_line[2], "(F")) catalog_data.fmt = 1;	// RA/DEC format, 1: FLOAT, 0: CHAR. Depends on catalogue
								break;
							case 2:		//2 being the index of the 3rd data, ie Name/Label of the object
								catalog_data.name_ind = cumul;
								catalog_data.name_len = F77_len;
								break;
							case 3:		//3 being the index of the 4th data, ie Vmag
								catalog_data.Vmag_ind = cumul;
								catalog_data.Vmag_len = F77_len;
								break;
						}
					}
					cumul = F77_len + 1 + cumul;	
				}
				g_free(sub_line);
			}			
		}	
	}
	g_strfreev(token);
	return 0;
}

int parse_vizier_buffer(const gchar *buffer, double lim_mag, int cata) {
	gchar **token = g_strsplit(buffer, "\n", -1);
	int nb_lines = g_strv_length(token);
	SirilWorldCS *world_cs = NULL;
	purge_temp_user_catalogue();

	if (!buffer || nb_lines <= 25) {				// 25 because nearly the end of the INFO lines in the header...
		siril_log_color_message(_("Useless file. Check if server available. Aborted.\n"), "red");
		g_strfreev(token);
		return 0;
	}

	double ra = 0.0, dec = 0.0;
	double mag = 0.0;
	int cat = cata - 0;
	gboolean valid = FALSE;
	gboolean unlock = FALSE;
	gboolean first_blank = FALSE;
	guint nbr_a = 0;
	gchar *objname = NULL;
	gchar *magg = NULL;

// Determines the data mapping of he used catalogue
	catalog_mapping(buffer, cat);

// Discards heder lines
	for (int i = 0; i < nb_lines; i++) {		// Loops on the lines
		gchar *line = token[i];
		guint line_len = strlen (line);
		valid = FALSE;

		if (!unlock) {
			if (!line[0] && first_blank) {
				siril_log_color_message(_("No object from %s found in this FOV.\n"), "salmon", catalog_to_str(cat));
				return 0;
			}
			if (g_str_has_prefix(line, "#INFO	Error")) {		// case described here: https://gitlab.com/free-astro/siril/-/issues/1131#note_1446793482
				siril_log_color_message(_("The VizieR server returned a connetion error. Try later.\n"), "salmon");
				return 0;
			}		
			if (!line[0]) first_blank = TRUE;					// Loops until the begining is found
			if (g_str_has_prefix(line, "-----")) unlock = TRUE;	// Flag line is found, unlock and go ahead
			continue;
		}

		if (!line[0] && unlock) break;				// The first blank line after a good data block is the end of the file (avoids "is_blank()" and a new dependency)

// If we reach this point, the file is valid and the index/lenghts have been found
// Reading the desired data in a valid line
		for (int j= 0; j < line_len; j++) {			// Loop across one line
			if (j == catalog_data.ra_ind) {			// RA read
				gchar *temp = g_strndup(&line[j], catalog_data.ra_len);
				ra = catalog_data.fmt ? g_ascii_strtod(temp, NULL) : parse_hms(temp);
				g_free(temp);
			}
			if (j == catalog_data.dec_ind) {			// DEC read
				gchar *temp = g_strndup(&line[j], catalog_data.dec_len);
				dec = catalog_data.fmt ? g_ascii_strtod(temp, NULL) : parse_dms(temp);
				g_free(temp);
			}
			if (j == catalog_data.Vmag_ind) 			// Magitude read
				magg = g_strndup(&line[j], catalog_data.Vmag_len);

			if (j == catalog_data.name_ind) {		// Name read
				gchar *temp = g_strndup(&line[j], catalog_data.name_len);	
				if (cat == CAT_PGC)	objname = g_strconcat("PGC ", temp, NULL);	// Adding a prefix for PGC. Much nicer!!
				else if (cat == CAT_BRIGHT_STARS)	objname = g_strconcat("HD ", temp, NULL);	// Adding a prefix for BSC. Much nicer also!!
				else objname = g_strdup(temp);
				if (!temp) objname = "x";		// if noname, just "x", but should not happen
				g_free(temp);
			}
		}
		valid = !isnan(dec) && !isnan(ra);

		if (valid) {
			world_cs = siril_world_cs_new_from_a_d(ra, dec);
			gchar *end;

			if (cat == CAT_PGC)	{
				mag = -99.0;							// PGC doesn't contain magnitude data
			} else mag = g_ascii_strtod(magg, &end);

			if (end == magg)
				valid = FALSE;
			else if (mag < lim_mag && is_inside(&gfit, ra, dec))
				nbr_a++;
			else valid = FALSE;
		}

// And finally displays the result (or not)
		if (world_cs && objname && valid) {
			gchar* alpha = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02ds");
			gchar* delta = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
			siril_log_message(_("Found %s at coordinates: %s, %s with magnitude: %0.1lf\n"),
					objname, alpha, delta, mag);
			g_free(alpha);
			g_free(delta);

			if (!com.script && !com.headless) {	// Write in catalogue only if in GUI mode
				add_object_in_catalogue(objname, world_cs, FALSE, USER_TEMP_CAT_INDEX);
				com.pref.gui.catalog[USER_TEMP_CAT_INDEX] = TRUE;	// and display it
				siril_world_cs_unref(world_cs);
				world_cs = NULL;
			}
		}
	}
	g_free(magg);
	g_free(objname);
	g_strfreev(token);
	siril_log_message(_("Found %d Objects with mag < %0.1lf\n"), nbr_a, lim_mag);
	return nbr_a <= 0;
}


int parse_conesearch_buffer(const gchar *buffer, double lim_mag) {
	if (!buffer || buffer[0] == '\0')
		return 1;
	if (!g_str_has_prefix(buffer, "# Flag:"))
		return 1;
	if (g_str_has_prefix(buffer, "# Flag: 0")) {
		siril_log_message("No Solar System object found in this FOV\n");
		return 1;
	}

	gchar **token = g_strsplit(buffer, "\n", -1);
	int nb_lines = g_strv_length(token);

	if (nb_lines <= 3) {
		g_strfreev(token);
		return 1;
	}

	guint nbr_a = 0;
	for (int i = 3; i < nb_lines; i++) {
		gchar *line = token[i];
		// format is '# Num | Name | RA(h) | DE(deg) | Class | Mv | Err(arcsec) | d(arcsec)'
		gchar *objname = NULL;
		char *start = strchr(line, '|');
		if (start && *(start+1) != '\0') {
			start += 2; // remove the space
			char *end = strchr(start, '|');
			if (end) {
				end -= 2; // remove the space
				objname = g_strndup(start, end-start+1);
			}
		} else continue;

		if (!objname) {
			siril_debug_print("Unsupported object type in response: %s\n", line);
			continue;
		}
		siril_debug_print("Object name: %s\n", objname);

		// Then, retrieve the coordinates
		gchar **fields = g_strsplit(line, "|", -1);
		guint n_fields = g_strv_length(fields);
		double mag = 0.0;
		//double hours = 0.0, min = 0.0, seconds = 0.0, degres = 0.0;
		double ra = 0.0, dec = 0.0;
		gboolean valid = FALSE;
		SirilWorldCS *world_cs = NULL;
		if (n_fields < 6) {
			g_strfreev(fields);
			continue;
		}
		ra = parse_hms(fields[2]);
		valid = !isnan(ra);
		if (valid) {
			dec = parse_dms(fields[3]);
			valid = !isnan(dec);
		}

		// Then, retrieve the magnitude, used to sort the objects to be displayed
		if (valid) {
			world_cs = siril_world_cs_new_from_a_d(ra, dec);

			gchar *end;
			mag = g_ascii_strtod(fields[5], &end);
			if (end == fields[5])
				valid = FALSE;
			else if (mag < lim_mag && is_inside(&gfit, ra, dec))
				nbr_a++;
			else valid = FALSE;
		}
		g_strfreev(fields);

		if (world_cs && objname && valid) {
			gchar *alpha = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02.2lfs");
			gchar *delta = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02.2lf\"");
			siril_log_message(_("Found %s at coordinates: %s, %s with magnitude: %0.1lf\n"),
					objname, alpha, delta, mag);
			g_free(alpha);
			g_free(delta);

			if (!com.script && !com.headless) {	// Write in catalogue only if in GUI mode
				add_object_in_catalogue(objname, world_cs, FALSE, USER_TEMP_CAT_INDEX);
				com.pref.gui.catalog[USER_TEMP_CAT_INDEX] = TRUE;	// and display it
				siril_world_cs_unref(world_cs);
				world_cs = NULL;
			}
		}
	}
	g_strfreev(token);
	siril_log_message(_("Found %d Solar System Objects with mag < %0.1lf\n"), nbr_a, lim_mag);
	return nbr_a <= 0;
}

/* this should be the new entry point for single object search from their name */
// TODO: the time and telescope location should be an input too for SSO, or the
// fit that contains them instead of assumed gfit like now
int cached_object_lookup(const gchar *name, psf_star **opt_result) {
	gboolean solarsystem = g_str_has_prefix(name, "a:") ||
		g_str_has_prefix(name, "c:") ||	g_str_has_prefix(name, "p:");
	gboolean found_it = FALSE;
	/* if we need the result, it means we need more than what the
	 * annotation catalogues contain, like photometric magnitudes
	 */
	if (!solarsystem && !opt_result) {
		// local search first, if useful
		const CatalogObjects *obj = search_in_annotations_by_name(name);
		if (obj) {
			siril_log_message(_("Object %s was found in the %s catalogue, RA = %f, Dec = %f\n"),
					get_catalogue_object_code(obj), cat_index_to_name(get_catalogue_object_cat(obj)),
					get_catalogue_object_ra(obj), get_catalogue_object_dec(obj));
			found_it = TRUE;
		}
	}

	if (!found_it) {
		// remote search then
		gchar *result;
		if (solarsystem)
			result = search_in_online_catalogs(name, QUERY_SERVER_EPHEMCC);
		else result = search_in_online_catalogs(name, QUERY_SERVER_SIMBAD_PHOTO);
		if (result && !parse_catalog_buffer(result, opt_result))
			found_it = TRUE;
		free_fetch_result(result);
	}

	return !found_it;
}

void on_search_objects_entry_activate(GtkEntry *entry, gpointer user_data) {
	if (!has_wcs(&gfit)) return;
	control_window_switch_to_tab(OUTPUT_LOGS);

	const gchar *name = gtk_entry_get_text(GTK_ENTRY(entry));
	gboolean found_it = cached_object_lookup(name, NULL) == 0;

	if (found_it) {
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


/********** object lookup for astrometry, can we merge both? ************/

struct sky_object platedObject[RESOLVER_NUMBER];

/* instead of populating platedObject from a resolver buffer, use the passed object */
void add_plated_from_annotations(const CatalogObjects *obj) {
	double ra = get_catalogue_object_ra(obj);
	double dec = get_catalogue_object_dec(obj);
	int resolver = RESOLVER_LOCAL;
	platedObject[resolver].world_cs = siril_world_cs_new_from_a_d(ra, dec);
	platedObject[resolver].imageCenter.x = ra;
	platedObject[resolver].imageCenter.y = dec;
	platedObject[resolver].south = (dec < 0.0);
	platedObject[resolver].name = g_strdup(get_catalogue_object_code(obj));
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
	//  In case no localisation data is available, set the reference to geocentric, ie @500 in the request
	if (fit->sitelat == 0.0 && fit->sitelong == 0.0)
		return g_strdup("@500");

	GString *formatted_site_coord = g_string_new("");

	if (fit->sitelat < 0.)
		formatted_site_coord = g_string_append(formatted_site_coord, "%2D");
	else formatted_site_coord = g_string_append(formatted_site_coord, "%2B");
	g_string_append_printf(formatted_site_coord, "%f,", fabs(fit->sitelat));

	if (fit->sitelong < 0.)
		formatted_site_coord = g_string_append(formatted_site_coord, "%2D");
	else formatted_site_coord = g_string_append(formatted_site_coord, "%2B");
	g_string_append_printf(formatted_site_coord, "%f,", fabs(fit->sitelong));

	g_string_append_printf(formatted_site_coord, "%f", fabs(fit->siteelev));

	return g_string_free(formatted_site_coord, FALSE);
}

// warning: uses gfit for ephemcc
// free the result with free_fetch_result
gchar *search_in_online_catalogs(const gchar *object, query_server server) {
#ifndef HAVE_NETWORKING
	siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
	return NULL;
#else
	GString *string_url = NULL;
	gchar *name = g_utf8_strdown(object, -1);
	switch(server) {
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
		string_url = g_string_new(EPHEMCC);
		string_url = g_string_append(string_url, "-name=");
		string_url = g_string_append(string_url, name);
		string_url = g_string_append(string_url, "&-type=");
		string_url = g_string_append(string_url, "&-ep=");
		gchar *formatted_date = date_time_to_FITS_date(gfit.date_obs);
		string_url = g_string_append(string_url, formatted_date);
		string_url = g_string_append(string_url, "&-nbd=1");
		string_url = g_string_append(string_url, "&-tscale=UTC");
		string_url = g_string_append(string_url, "&-observer=");
		gchar *formatted_site = retrieve_site_coord(&gfit);
		string_url = g_string_append(string_url, formatted_site);
		string_url = g_string_append(string_url, "&-theory=INPOP");
		string_url = g_string_append(string_url, "&-teph=1");
		string_url = g_string_append(string_url, "&-tcoor=5");
		string_url = g_string_append(string_url, "&-oscelem=astorb");
		string_url = g_string_append(string_url, "&-mime=text/csv");
		string_url = g_string_append(string_url, "&-output=--jd");
		string_url = g_string_append(string_url, "&-from=Siril;");
		if (!gfit.date_obs) {
			siril_log_color_message(_("This command only works on images that have observation date information\n"), "red");
			g_string_free(string_url, TRUE);
			return NULL;
		}
		siril_log_message(_("Searching for solar system object %s on observation date %s\n"),
				name, formatted_date);
		if (gfit.sitelat == 0.0 && gfit.sitelong == 0.0) {
			siril_log_color_message(_("No topocentric data available. Set to geocentric\n"), "salmon");
		} else {
			siril_log_message(_("at lat: %f, long: %f, alt: %f\n"), gfit.sitelat,
				gfit.sitelong, gfit.siteelev);
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

	g_free(name);
	gchar *url = g_string_free(string_url, FALSE);
	gchar *cleaned_url = url_cleanup(url);
	g_free(url);
	siril_debug_print("URL: %s\n", cleaned_url);

	char *result = fetch_url(cleaned_url);

	g_free(cleaned_url);
	return result;
#endif
}
