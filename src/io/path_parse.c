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
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/siril_world_cs.h"
#include "io/fits_sequence.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"


// static const char *path_parse_error_to_string(path_parse_errors err) {
// 	switch (err) {
// 	case PATH_PARSE_KEY_NOT_FOUND:
// 		return "Key not found\n";
// 	case PATH_PARSE_WRONG_CALL:
// 		return "Call to read_key_from_header_text is malformed\n";
// 	case PATH_PARSE_HEADER_NULL:
// 		return "Header is empty\n";
// 	case PATH_PARSE_WRONG_FORMAT:
// 		return "Wrong format\n";
// 	case PATH_PARSE_UNSUPPORTED_FORMAT:
// 		return "Unsupported format\n";
// 	default:
// 		return NULL;
// 	}
// }

static path_parse_errors read_key_from_header_text(gchar **headers, gchar *key, double *numvalue, gchar *strvalue) {
	path_parse_errors status = PATH_PARSE_OK;
	gboolean keyfound = FALSE;
	char searchstr[10];
	g_sprintf(searchstr, "%-8s=", key);
	for (int i = 0; i < g_strv_length(headers); i++) {
		if (g_str_has_prefix(headers[i], searchstr)) {
			keyfound = TRUE;
			gchar **subs = g_strsplit(headers[i], "=", 2);
			gchar **valsubs = g_strsplit(subs[1], "/", 2);
			if (numvalue) {
				*numvalue = g_ascii_strtod(valsubs[0], NULL); 
			} else if (strvalue) {
				gchar *currstr = g_strdup(valsubs[0]);
				currstr = g_shell_unquote(currstr, NULL);
				strncpy(strvalue, currstr, FLEN_VALUE - 1);
				g_free(currstr);
			} else {
				g_free(valsubs);
				g_free(subs);
				return PATH_PARSE_WRONG_CALL;
			}
			g_strfreev(subs);
			g_strfreev(valsubs);
			break;
		}
	}
	if (!keyfound) return PATH_PARSE_KEY_NOT_FOUND;
	return status;
}

/*
This function takes an expression with potentially a wildcard in it
Searches the directory for files matching the pattern
and returns the first occurence of the match if any
*/
static gchar *wildcard_check(gchar *expression, int *status) {
	gchar *out = NULL, *dirname = NULL, *basename = NULL;
	const gchar *file = NULL;
	GDir *dir;
	GError *error = NULL;
	GPatternSpec *patternspec;
	gint count = 0;

	if (!g_utf8_strchr(expression, -1, '*')) {
		return g_strdup(expression);
	}

	dirname = g_path_get_dirname(expression);
	basename = g_path_get_basename(expression);

	if ((dir = g_dir_open(dirname, 0, &error)) == NULL) {
		siril_debug_print("wildcard dircheck: %s\n", error->message);
		siril_log_color_message(_("Problem parsing expression %s - Error code %d - aborting\n"), "red", dirname, *status);
		g_clear_error(&error);
		g_free(dirname);
		g_free(basename);
		return out;
	}
	patternspec = g_pattern_spec_new(basename);
	while ((file = g_dir_read_name(dir)) != NULL) {
#if GLIB_CHECK_VERSION(2,70,0)
		if (g_pattern_spec_match(patternspec, strlen(file), file, NULL)) {
#else
		if (g_pattern_match(patternspec, strlen(file), file, NULL)) {
#endif
			if (!count) {
				out = g_build_filename(dirname, file, NULL);
			}
			count++;
		}
	}
	if (count > 1) {
		*status = PATH_PARSE_MORE_THAN_ONE_HIT;
		siril_log_color_message(_("Found more than one file matching %s - Warning code %d\nUsing first:%s\n"), "salmon", basename, *status, out);

	} else if (!count) {
		*status = PATH_PARSE_MORE_THAN_ONE_HIT;
		siril_log_color_message(_("Problem parsing expression %s - Error code %d - aborting\n"), "red", expression, *status);
	}
	g_free(dirname);
	g_free(basename);
	g_pattern_spec_free(patternspec);
	return out;
}

/*
This is the main function. It takes as argument a fit, an expresion to be parsed
and a mode (read or write).
It searches for reserved keywords (starting with "lib")
and fetches the libs set in preferences if required.
It then parses all the tokens in between $ signs in the form $KEY:fmt$ where KEY
is a valid HEADER key and fmt either a %d %f or %s specifier or special formatters
such as "dm12" (date minus 12 hrs), "dm0" (date), "ra", "dec", "ran", "decn" which format
RA and DEC values (the suffix "n" indicates the input key should be read as numerical)
It returns the expression with all the tokens replaced by the formatted values
In read mode, it also replaces * by searching the directory for a file matching the pattern
*/
gchar *path_parse(fits *fit, gchar *expression, path_parse_mode mode, int *status) {
	gchar *out = NULL, *localexpression = NULL;
	if (!fit->header) {
		*status = PATH_PARSE_HEADER_NULL;
		return out;
	}
	if (g_str_has_prefix(expression, "lib")) { // using reserved keywords libbias, libdark, libflat
		if (!g_strcmp0(expression + 3, "bias")) {
			localexpression = g_strdup(com.pref.prepro.bias_lib);
		} else if (!g_strcmp0(expression + 3, "dark")) {
			localexpression = g_strdup(com.pref.prepro.dark_lib);
		} else if (!g_strcmp0(expression + 3, "flat")) {
			localexpression = g_strdup(com.pref.prepro.flat_lib);
		} else {
			*status = PATH_PARSE_WRONG_RESERVED_KEYWORD;
			siril_log_color_message(_("Unknown reserved keyword %s - Error code %d - aborting\n"), "red", expression, *status);
			return out;
		}
		if (strlen(localexpression) == 0) {
			*status = PATH_PARSE_LIBRARY_NOTDEFINED;
			siril_log_color_message(_("Library %s is not defined in preferences - Error code %d - aborting\n"), "red", expression, *status);
			g_free(localexpression);
			return out;
		}
	} else {
		localexpression = g_strdup(expression);
	}
	gchar *pattern = "\\$(.+?)\\$";
	gchar **tokens = g_regex_split_simple(pattern, localexpression, G_REGEX_RAW, 0);
	gchar **headerkeys = g_strsplit(fit->header, "\n", 0);
	for (int i = 0; i < g_strv_length(tokens); i++) {
		gchar **subs = g_strsplit(tokens[i], ":", 2);
		if (g_strv_length(subs) == 1) {
			g_strfreev(subs);
			continue;
		}
		if (strlen(subs[0]) == 1 && subs[1][0] == '\\') { // dealing with Windows drive letter "C:""
			g_strfreev(subs);
			continue;
		}
		gchar buf[FLEN_VALUE];
		gchar key[9];
		buf[0] = '\0';
		// Check if expression starts with a * wildcard
		// Behavior will depend if we are in read or write mode
		if (subs[0][0] == '*') {
			if (mode == PATH_PARSE_READ) {
				strncpy(buf, "*", FLEN_VALUE - 1);
			} else {
				strncpy(key, subs[0] + 1, 9);
			}
		} else {
			strncpy(key, subs[0], 9);
		}
		if (buf[0] == '*') {
			printf("Wildcard in READ mode\n");
		} else if (g_str_has_suffix(subs[1], "d") || g_str_has_suffix(subs[1], "f")) { // case %d or %f
			gboolean isint = g_str_has_suffix(subs[1], "d");
			double val;
			*status = read_key_from_header_text(headerkeys, key, &val, NULL);
			if (*status) {
				siril_log_color_message(_("Problem reading keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			(isint) ? sprintf(buf, subs[1], (int)val) : sprintf(buf, subs[1], val);
		} else if (g_str_has_suffix(subs[1], "s")) { // case %s
			char val[FLEN_VALUE];
			*status = read_key_from_header_text(headerkeys, key, NULL, val);
			if (*status) {
				siril_log_color_message(_("Problem reading keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			remove_spaces_from_str(val);
			sprintf(buf, subs[1], val); // just in case there is a fancy formatting directive
		} else if (g_str_has_prefix(subs[1],"dm")) { // case dm12 - date minus 12hrs or dm0
			double minus_hour = -1. * g_ascii_strtod(subs[1] + 2, NULL);
			char val[FLEN_VALUE];
			*status = read_key_from_header_text(headerkeys, key, NULL, val);
			if (*status) {
				siril_log_color_message(_("Problem reading keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			GDateTime *read_time = FITS_date_to_date_time(val);
			if (!read_time) {
				*status = PATH_PARSE_WRONG_DATE;
				siril_log_color_message(_("Could not read date from %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			GDateTime *read_time_corr = g_date_time_add_hours(read_time, minus_hour);
			gchar *fmtdate = date_time_to_date(read_time_corr);
			strncpy(buf, fmtdate, FLEN_VALUE - 1);
			g_date_time_unref(read_time);
			g_date_time_unref(read_time_corr);
			g_free(fmtdate);
		} else if (g_str_has_prefix(subs[1], "ra") || g_str_has_prefix(subs[1], "dec")) { // case ra and dec (str), ran and decn (num)
			gboolean is_ra = g_str_has_prefix(subs[1],"ra");
			gboolean is_float = g_str_has_suffix(subs[1],"n");
			SirilWorldCS *target_coords = NULL;
			char val[FLEN_VALUE];
			double valf;
			if (is_float)
				*status = read_key_from_header_text(headerkeys, key, &valf, NULL);
			else
				*status = read_key_from_header_text(headerkeys, key, NULL, val);
			if (*status) {
				siril_log_color_message(_("Problem reading keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			if (is_ra)
				if (is_float)
					target_coords = siril_world_cs_new_from_a_d(valf, 0.0);
				else
					target_coords = siril_world_cs_new_from_objct_ra_dec(val, "+00 00 00");
			else
				if (is_float)
					target_coords = siril_world_cs_new_from_a_d(0.0, valf);
				else
					target_coords = siril_world_cs_new_from_objct_ra_dec("00 00 00", val);
			if (!target_coords) {
				*status = PATH_PARSE_WRONG_WCS;
				siril_log_color_message(_("Problem reading keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			gchar *fmtcoord;
			if (is_ra)
				fmtcoord = siril_world_cs_alpha_format(target_coords, "%02dh%02dm%02ds");
			else
				fmtcoord = siril_world_cs_delta_format(target_coords, "%c%02dd%02dm%02ds");
			strncpy(buf, fmtcoord, FLEN_VALUE - 1);
			siril_world_cs_unref(target_coords);
			g_free(fmtcoord);
		} else {
			siril_log_color_message(_("Unsupported format %s - aborting\n"), "red", subs[1], *status);
			g_strfreev(subs);
			goto free_and_exit;
		}
		if (buf[0] == 0) {
			*status = PATH_PARSE_WRONG_FORMAT;
			siril_log_color_message(_("Problem parsing expression %s - Error code %d - aborting\n"), "red", tokens[i], *status);
			g_strfreev(subs);
			goto free_and_exit;
		}
		g_free(tokens[i]);
		tokens[i] = g_strdup(buf);
		g_strfreev(subs);
	}
	out = g_strjoinv("", tokens);
	siril_debug_print("String in: %s\n", expression);
	siril_debug_print("String out: %s\n", out);
	if (mode == PATH_PARSE_READ) {
		gchar *foundmatch = wildcard_check(out, status);
		g_free(out);
		out = NULL;
		if (*status <= 0) out = g_strdup(foundmatch);
		g_free(foundmatch);
	}
free_and_exit:
	if (tokens) g_strfreev(tokens);
	if (headerkeys) g_strfreev(headerkeys);
	if (localexpression) g_free(localexpression);
	return out;
}
