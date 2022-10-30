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

static void display_path_parse_error(pathparse_errors err, gchar *addstr) {
	if (!err) return;
	gchar *endstr = (err < 0) ? _("going on") : _("aborting");
	gchar *msg = NULL;
	switch (err) {
		case PATHPARSE_ERR_HEADER_NULL:
		case PATHPARSE_ERR_HEADER_NULL_NOFAIL:
			msg = _("Header is empty");
			break;
		case PATHPARSE_ERR_WRONG_RESERVED_KEYWORD:
		case PATHPARSE_ERR_WRONG_RESERVED_KEYWORD_NOFAIL:
			msg = _("Wrong reserved keyword : ");
			break;
		case PATHPARSE_ERR_LIBRARY_NOTDEFINED:
		case PATHPARSE_ERR_LIBRARY_NOTDEFINED_NOFAIL:
			msg = _("Library not defined in preferences : ");
			break;
		case PATHPARSE_ERR_KEY_NOT_FOUND:
		case PATHPARSE_ERR_KEY_NOT_FOUND_NOFAIL:
			msg = _("Key not found : ");
			break;
		case PATHPARSE_ERR_WRONG_DATE:
		case PATHPARSE_ERR_WRONG_DATE_NOFAIL:
			msg = _("Wrong date : ");
			break;
		case PATHPARSE_ERR_WRONG_WCS:
		case PATHPARSE_ERR_WRONG_WCS_NOFAIL:
			msg = _("Wrong coordinates: ");
			break;
		case PATHPARSE_ERR_UNSUPPORTED_FORMAT:
		case PATHPARSE_ERR_UNSUPPORTED_FORMAT_NOFAIL:
			msg = _("Unsupported format : ");
			break;
		case PATHPARSE_ERR_MORE_THAN_ONE_HIT:
			msg = _("More than one match for : ");
			break;
		case PATHPARSE_ERR_NO_HIT_FOUND:
			msg = _("No match found for : ");
			break;
		case PATHPARSE_ERR_NO_DIR:
			msg = _("Problem with path : ");
			break;
		case PATHPARSE_ERR_WRONG_CALL:
		default:
			msg = _("Internal error");
			break;
	}
	siril_log_color_message(_("Err %3d - %s%s - %s\n"), "red", err, msg, addstr, endstr);
}

static pathparse_errors read_key_from_header_text(gchar **headers, gchar *key, double *numvalue, gchar *strvalue) {
	pathparse_errors status = PATHPARSE_ERR_OK;
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
				status = PATHPARSE_ERR_WRONG_CALL; // internal error, should not be thrown
				return status;
			}
			g_strfreev(subs);
			g_strfreev(valsubs);
			break;
		}
	}
	if (!keyfound) return PATHPARSE_ERR_KEY_NOT_FOUND;
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
		*status = PATHPARSE_ERR_NO_DIR;
		display_path_parse_error(*status, dirname);
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
		*status = PATHPARSE_ERR_MORE_THAN_ONE_HIT;
		display_path_parse_error(*status, basename);
		siril_log_color_message(_("Using first match: %s\n"), "salmon", out);
	} else if (!count) {
		*status = PATHPARSE_ERR_NO_HIT_FOUND;
		display_path_parse_error(*status, expression);
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
In write mode, the * is omitted and the token is parsed as per specifier
In write mode "nofail", it will try to return something no matter what
*/
gchar *path_parse(fits *fit, gchar *expression, pathparse_mode mode, int *status) {
	gchar *out = NULL, *localexpression = NULL;
	int nofail = (mode == PATHPARSE_MODE_WRITE_NOFAIL) ? -1 : 1; // statuses in nofail mode are transformed to warnings if multiplied by -1
	if (!fit->header) {
		*status = nofail * PATHPARSE_ERR_HEADER_NULL;
		if (*status > 0) {
			display_path_parse_error(*status, NULL);
			return out;
		}
	} 
	if (g_str_has_prefix(expression, "lib")) { // using reserved keywords libbias, libdark, libflat
		if (!g_strcmp0(expression + 3, "bias")) {
			localexpression = g_strdup(com.pref.prepro.bias_lib);
		} else if (!g_strcmp0(expression + 3, "dark")) {
			localexpression = g_strdup(com.pref.prepro.dark_lib);
		} else if (!g_strcmp0(expression + 3, "flat")) {
			localexpression = g_strdup(com.pref.prepro.flat_lib);
		} else {
			*status = nofail * PATHPARSE_ERR_WRONG_RESERVED_KEYWORD;
			display_path_parse_error(*status, expression);
			out = (*status > 0) ? NULL : g_strdup(expression); // using libsmthg as a fallback
			g_free(localexpression);
			return out;
		}
		if (strlen(localexpression) == 0) {
			*status = nofail * PATHPARSE_ERR_LIBRARY_NOTDEFINED;
			display_path_parse_error(*status, expression);
			out = (*status > 0) ? NULL : g_strdup(expression); // using libsmthg as a fallback
			g_free(localexpression);
			return out;
		}
	} else {
		localexpression = g_strdup(expression);
	}
	gchar *pattern = "\\$(.+?)\\$";
	gchar **tokens = g_regex_split_simple(pattern, localexpression, G_REGEX_RAW, 0);
	gchar **headerkeys = NULL;
	if (fit->header) { // avoiding split in case we are in nofail and header is empty
		headerkeys = g_strsplit(fit->header, "\n", 0);
	}
	for (int i = 0; i < g_strv_length(tokens); i++) {
		gchar **subs = g_strsplit(tokens[i], ":", 2);
		if (g_strv_length(subs) == 1) {
			g_strfreev(subs);
			continue;
		}
		if (strlen(subs[0]) == 1 && subs[1][0] == '\\') { // dealing with Windows drive letter "C:"
			g_strfreev(subs);
			continue;
		}
		gchar buf[FLEN_VALUE];
		gchar key[9];
		buf[0] = '\0';
		// Check if expression starts with a * wildcard
		// Behavior will depend if we are in read or write mode
		if (subs[0][0] == '*') {
			if (mode == PATHPARSE_MODE_READ) {
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
			*status = nofail * read_key_from_header_text(headerkeys, key,&val, NULL);
			display_path_parse_error(*status, key);
			if (*status > 0) {
				g_strfreev(subs);
				goto free_and_exit;
			}
			if (*status == 0) // can still be neg with write_nofail
				(isint) ? sprintf(buf, subs[1], (int)val) : sprintf(buf, subs[1], val);
		} else if (g_str_has_suffix(subs[1], "s")) { // case %s
			char val[FLEN_VALUE];
			*status = nofail * read_key_from_header_text(headerkeys, key, NULL, val);
			display_path_parse_error(*status, key);
			if (*status > 0) {
				g_strfreev(subs);
				goto free_and_exit;
			}
			if (*status == 0) { // can still be neg with write_nofail
				remove_spaces_from_str(val);
				sprintf(buf, subs[1], val); // just in case there is a fancy formatting directive
			}
		} else if (g_str_has_prefix(subs[1],"dm")) { // case dm12 - date minus 12hrs or dm0
			double minus_hour = -1. * g_ascii_strtod(subs[1] + 2, NULL);
			char val[FLEN_VALUE];
			*status = nofail * read_key_from_header_text(headerkeys, key, NULL, val);
			display_path_parse_error(*status, key);
			if (*status > 0) {
				g_strfreev(subs);
				goto free_and_exit;
			}
			if (*status == 0) {
				GDateTime *read_time = FITS_date_to_date_time(val);
				if (!read_time) {
					*status = nofail * PATHPARSE_ERR_WRONG_DATE;
					display_path_parse_error(*status, key);
					if (status > 0) {
						g_strfreev(subs);
						goto free_and_exit;
					}
				}
				if (*status == 0) {
					GDateTime *read_time_corr = g_date_time_add_hours(read_time, minus_hour);
					gchar *fmtdate = date_time_to_date(read_time_corr);
					strncpy(buf, fmtdate, FLEN_VALUE - 1);
					g_date_time_unref(read_time);
					g_date_time_unref(read_time_corr);
					g_free(fmtdate);
				}
			}
		} else if (g_str_has_prefix(subs[1], "ra") || g_str_has_prefix(subs[1], "dec")) { // case ra and dec (str), ran and decn (num)
			gboolean is_ra = g_str_has_prefix(subs[1],"ra");
			gboolean is_float = g_str_has_suffix(subs[1],"n");
			SirilWorldCS *target_coords = NULL;
			char val[FLEN_VALUE];
			double valf;
			if (is_float)
				*status = nofail * read_key_from_header_text(headerkeys, key,&valf, NULL);
			else
				*status = nofail * read_key_from_header_text(headerkeys, key, NULL, val);
			display_path_parse_error(*status, key);
			if (*status > 0) {
				g_strfreev(subs);
				goto free_and_exit;
			}
			if (*status == 0) {
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
					*status = nofail * PATHPARSE_ERR_WRONG_WCS;
					display_path_parse_error(*status, key);
					if (*status > 0) {
						g_strfreev(subs);
						goto free_and_exit;
					}
				}
				if (*status == 0) {
					gchar *fmtcoord;
					if (is_ra)
						fmtcoord = siril_world_cs_alpha_format(target_coords, "%02dh%02dm%02ds");
					else
						fmtcoord = siril_world_cs_delta_format(target_coords, "%c%02dd%02dm%02ds");
					strncpy(buf, fmtcoord, FLEN_VALUE - 1);
					siril_world_cs_unref(target_coords);
					g_free(fmtcoord);
				}
			}
		} else {
			*status = nofail * PATHPARSE_ERR_UNSUPPORTED_FORMAT;
			display_path_parse_error(*status, subs[1]);
			if (*status > 0) {
				g_strfreev(subs);
				goto free_and_exit;
			}
		}
		if (buf[0] == '\0' && mode == PATHPARSE_MODE_WRITE_NOFAIL) {
			strncpy(buf, key, 9);
		}
		g_free(tokens[i]);
		tokens[i] = g_strdup(buf);
		g_strfreev(subs);
	}
	out = g_strjoinv("", tokens);
	siril_debug_print("String in: %s\n", expression);
	siril_debug_print("String out: %s\n", out);
	if (mode == PATHPARSE_MODE_READ) {
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
