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
				remove_spaces_from_str(currstr);
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

gchar *path_parse(fits *fit, gchar *expression, int *status) {
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
		if (g_str_has_suffix(subs[1], "d") || g_str_has_suffix(subs[1], "f")) { // case %d or %f
			gboolean isint = g_str_has_suffix(subs[1], "d");
			double val;
			*status = read_key_from_header_text(headerkeys, subs[0], &val, NULL);
			if (*status) {
				siril_log_color_message(_("Problem reading keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			(isint) ? sprintf(buf, subs[1], (int)val) : sprintf(buf, subs[1], val);
		} else if (g_str_has_suffix(subs[1], "s")) { // case %s
			char val[FLEN_VALUE];
			*status = read_key_from_header_text(headerkeys, subs[0], NULL, val);
			if (*status) {
				siril_log_color_message(_("Problem reading keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
				g_strfreev(subs);
				goto free_and_exit;
			}
			sprintf(buf, subs[1], val); // just in case there is a fancy formatting directive like uppercase or else
		} else {
			siril_log_color_message(_("Unsupported format %s - aborting\n"), "red", subs[1], *status);
			g_strfreev(subs);
			goto free_and_exit;
		}
		if (buf[0] == 0) {
			*status = PATH_PARSE_WRONG_FORMAT;
			siril_log_color_message(_("Problem parsing keyword %s - Error code %d - aborting\n"), "red", subs[0], *status);
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
free_and_exit:
	if (tokens) g_strfreev(tokens);
	if (localexpression) g_free(localexpression);
	return out;
}
