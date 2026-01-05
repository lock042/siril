/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/siril_world_cs.h"
#include "io/fits_sequence.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"

#define DM_PATTERN "\\d{4}-\\d{2}-\\d{2}"
#define DT_PATTERN  "\\d{4}-\\d{2}-\\d{2}T\\d{2}-\\d{2}-\\d{2}"

static void display_path_parse_error(pathparse_errors err, const gchar *addstr) {
	if (!err) return;
	gchar *startstr = (err < 0) ? _("Warning code:"): _("Error code:");
	gchar *endstr = (err < 0) ? _("going on") : _("aborting");
	gchar *msg = NULL;
	gchar *color = (err < 0) ? "salmon" : "red";
	gchar addbuf[256];
	g_snprintf(addbuf, 255, "%s", (!addstr) ? "" : addstr);

	switch (err) {
		case PATHPARSE_ERR_HEADER_NULL:
		case PATHPARSE_ERR_HEADER_NULL_NOFAIL:
			msg = _("Header is empty");
			break;
		case PATHPARSE_ERR_WRONG_RESERVED_KEYWORD:
		case PATHPARSE_ERR_WRONG_RESERVED_KEYWORD_NOFAIL:
			msg = _("Wrong reserved keyword: ");
			break;
		case PATHPARSE_ERR_LIBRARY_NOTDEFINED:
		case PATHPARSE_ERR_LIBRARY_NOTDEFINED_NOFAIL:
			msg = _("Library not defined in preferences: ");
			break;
		case PATHPARSE_ERR_KEY_NOT_FOUND:
		case PATHPARSE_ERR_KEY_NOT_FOUND_NOFAIL:
			msg = _("Key not found: ");
			break;
		case PATHPARSE_ERR_WRONG_DATE:
		case PATHPARSE_ERR_WRONG_DATE_NOFAIL:
			msg = _("Wrong date: ");
			break;
		case PATHPARSE_ERR_WRONG_WCS:
		case PATHPARSE_ERR_WRONG_WCS_NOFAIL:
			msg = _("Wrong coordinates: ");
			break;
		case PATHPARSE_ERR_UNSUPPORTED_FORMAT:
		case PATHPARSE_ERR_UNSUPPORTED_FORMAT_NOFAIL:
			msg = _("Unsupported format: ");
			break;
		case PATHPARSE_ERR_NOSEQLOADED:
		case PATHPARSE_ERR_NOSEQLOADED_NOFAIL:
			msg = _("No sequence loaded");
			break;
		case PATHPARSE_ERR_BADSTRING:
		case PATHPARSE_ERR_BADSTRING_NOFAIL:
			msg = _("Wrongly formatted string: ");
			break;
		case PATHPARSE_ERR_WILDCARD_SYNTAX:
		case PATHPARSE_ERR_WILDCARD_SYNTAX_NOFAIL:
			msg = _("Wrong wildcard usage: ");
			break;
		case PATHPARSE_ERR_MORE_THAN_ONE_HIT:
			msg = _("More than one match for: ");
			break;
		case PATHPARSE_ERR_HITS_ALL_NEWER:
			msg = _("All hits newer than: ");
			break;
		case PATHPARSE_ERR_NO_HIT_FOUND:
			msg = _("No match found for: ");
			break;
		case PATHPARSE_ERR_NO_DIR:
			msg = _("Problem with path: ");
			break;
		case PATHPARSE_ERR_TMPFIT:
			msg = _("Could not create temp FITS for header update: ");
			break;
		case PATHPARSE_ERR_WRONG_CALL:
		default:
			msg = _("Internal error");
			break;
	}
	siril_log_color_message("%s %d - %s%s - %s\n", color, startstr, err, msg, addbuf, endstr);
}

pathparse_errors read_key_from_header_text(gchar **headers, gchar *key, double *numvalue, gchar *strvalue) {
	pathparse_errors status = PATHPARSE_ERR_OK;
	if (!headers)
		return PATHPARSE_ERR_HEADER_NULL;
	gboolean keyfound = FALSE;
	char searchstr[10];
	g_sprintf(searchstr, "%-8s=", key);
	for (int i = 0; i < g_strv_length(headers); i++) {
		if (g_str_has_prefix(headers[i], searchstr)) {
			char val[FLEN_VALUE];
			int fstatus = 0;
			fits_parse_value(headers[i], val, NULL, &fstatus);
			if (fstatus) {
				return PATHPARSE_ERR_BADSTRING;
			}
			keyfound = TRUE;
			if (numvalue) {
				*numvalue = g_ascii_strtod(val, NULL);
			} else if (strvalue) {
				gchar *ucurrstr = g_shell_unquote(val, NULL);
				if (ucurrstr)
					strncpy(strvalue, ucurrstr, FLEN_VALUE - 1);
				else
					status = PATHPARSE_ERR_BADSTRING;
				g_free(ucurrstr);
			} else {
				status = PATHPARSE_ERR_WRONG_CALL; // internal error, should not be thrown
				return status;
			}
			break;
		}
	}
	if (!keyfound) return PATHPARSE_ERR_KEY_NOT_FOUND;
	return status;
}

typedef struct _file_date {
	GDateTime *date;
	gchar *filename;
} file_date;

static file_date *new_file_date(const gchar *date, const gchar *filename, gboolean isdatetime) {
	file_date *fd = g_slice_new(file_date);
	gboolean set = FALSE;
	gint year, month, day, hour, min, sec;
	GTimeZone *tz = g_time_zone_new_utc();
	if (isdatetime) {
		if (sscanf(date, "%d-%d-%dT%d-%d-%d", &year, &month, &day, &hour, &min,	&sec) == 6) {
			fd->date = g_date_time_new(tz, year, month, day, hour, min, (double)sec);
			set = TRUE;
		}
	} else {
		if (sscanf(date, "%d-%d-%d", &year, &month, &day) == 3) {
			fd->date = g_date_time_new(tz, year, month, day, 0, 0, 0.);
			set = TRUE;
		}
	}
	g_time_zone_unref(tz);
	if (!set) {
		g_slice_free(file_date, fd);
		return NULL;
	}
	fd->filename = g_strdup(filename);
	return fd;
}

static gint file_date_compare(const file_date *a, const file_date *b) {
	int res = g_date_time_compare(a->date, b->date);
	if (!res)
		return (strlen(a->filename) == 1) ? 1 : -1; // we want filename "." (the ref file), to be after
	return res;
}

static gint file_date_match_file(const file_date *a, gchar *data) {
	return g_strcmp0(a->filename, data);
}

static void file_date_free(gpointer data) {
	file_date *fd = (file_date *)data;
	g_free(fd->filename);
	g_date_time_unref(fd->date);
	g_slice_free(file_date, fd);
}

/*
This function takes an expression with potentially a wildcard in it
Searches the directory for files matching the pattern
and returns the first occurence of the match if any
*/
static gchar *wildcard_check(gchar *expression, int *status, gchar *target_date, gboolean isdatetime) {
	gchar *out = NULL, *dirname = NULL, *basename = NULL, *currfile = NULL;
	const gchar *file = NULL;
	GDir *dir;
	GError *error = NULL;
	gint count = 0;
	struct stat fileInfo, currfileInfo;
	GList *fds = NULL;

	//we need to check that the folder path does not contain wildcards
	GString *newexp = g_string_new(expression);
	g_string_replace(newexp, (isdatetime) ? DT_PATTERN : DM_PATTERN, "[DATE]", 0);
	gchar *tmpexp = g_string_free(newexp, FALSE);
	gchar *tmpdirname = g_path_get_dirname(tmpexp);
	gchar *tmpbasename = g_path_get_basename(tmpexp);
	if (g_strstr_len(tmpdirname, -1, "[DATE]") || g_strstr_len(tmpdirname, -1, "*")) {
		*status = PATHPARSE_ERR_WILDCARD_SYNTAX;
		display_path_parse_error(*status, tmpdirname);
		g_free(tmpexp);
		g_free(tmpdirname);
		g_free(tmpbasename);
		return NULL;
	}
	GString *newdir = g_string_new(tmpdirname);
	g_string_replace(newdir, "[DATE]", (isdatetime) ? DT_PATTERN : DM_PATTERN, 0);
	dirname =  g_string_free(newdir, FALSE);
	GString *newbase = g_string_new(tmpbasename);
	g_string_replace(newbase, "[DATE]", (isdatetime) ? DT_PATTERN : DM_PATTERN, 0);
	basename =  g_string_free(newbase, FALSE);

	if ((dir = g_dir_open(dirname, 0, &error)) == NULL) {
		siril_debug_print("wildcard dircheck: %s\n", error->message);
		*status = PATHPARSE_ERR_NO_DIR;
		display_path_parse_error(*status, dirname);
		g_clear_error(&error);
		g_free(dirname);
		g_free(basename);
		g_free(tmpdirname);
		g_free(tmpbasename);
		return out;
	}

	if (target_date) { // Init with target_date
		fds = g_list_append(fds, new_file_date(target_date, ".", isdatetime));
	}

	while ((file = g_dir_read_name(dir)) != NULL) {
		if (g_regex_match_simple(basename, file, G_REGEX_RAW, 0)) {
			if (!target_date) { // No date_obs we fetch the most recent file based on stat
				if (!count) {
					out = g_build_filename(dirname, file, NULL);
						if (stat(out, &fileInfo))
							siril_debug_print("stat() failed\n");
				} else {
					currfile = g_build_filename(dirname, file, NULL);
					if (stat(currfile, &currfileInfo))
							siril_debug_print("stat() failed\n");
					if (currfileInfo.st_ctime > fileInfo.st_ctime) { // currfile is more recent
						out = g_build_filename(dirname, file, NULL);
						if (stat(out, &fileInfo))
							siril_debug_print("stat() failed\n");
					}
				}
				count++;
			} else {
				out = g_build_filename(dirname, file, NULL);
				GRegex *regex = g_regex_new((isdatetime) ? DT_PATTERN : DM_PATTERN, 0, 0, NULL);
				GMatchInfo *match_info;
				g_regex_match(regex, out, 0, &match_info);
				if (g_match_info_matches(match_info)) {
					gchar *date = g_match_info_fetch(match_info, 0);
					siril_debug_print("Captured date: %s\n", date);
					fds = g_list_append(fds, new_file_date(date, out, isdatetime));
					g_free(date);
					count++;
				}
				g_free(out);
				out = NULL;
				g_match_info_free(match_info);
				g_regex_unref(regex);
			}
		}
	}

	if (!target_date && count > 1) {
		*status = PATHPARSE_ERR_MORE_THAN_ONE_HIT;
		display_path_parse_error(*status, basename);
		siril_log_color_message(_("Using most recent matching file: %s\n"), "salmon", out);
	} else if (target_date && count >= 1) { // we sort the date_file list
		fds = g_list_sort(fds, (GCompareFunc)file_date_compare);
		GList *ref = g_list_find_custom(fds, ".", (GCompareFunc)file_date_match_file);
		gint refpos = g_list_position(fds, ref);
		if (refpos == 0) { // all files matching have newer DATE-OBS
			out = g_strdup(((file_date *)g_list_nth_data(fds, 1))->filename);
			*status = PATHPARSE_ERR_HITS_ALL_NEWER;
			display_path_parse_error(*status, basename);
			siril_log_color_message(_("Using the closest matching file: %s\n"), "salmon", out);
		} else { // we return the file which has DATE-OBS lteq to reference
			out = g_strdup(((file_date *)g_list_nth_data(fds, refpos - 1))->filename);
		}
	} else if (!count) {
		*status = PATHPARSE_ERR_NO_HIT_FOUND;
		display_path_parse_error(*status, expression);
	}
	g_dir_close(dir);
	g_free(dirname);
	g_free(basename);
	g_free(currfile);
	g_free(tmpdirname);
	g_free(tmpbasename);
	g_list_free_full(fds, (GDestroyNotify)file_date_free);
	return out;
}

/*
This is the main function. It takes as argument a fit, an expresion to be parsed
and a mode (read or write). It returns a newly allocated string.
It first searches for reserved keywords (starting with "$lib")
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
gchar *path_parse(fits *fit, const gchar *expression, pathparse_mode mode, int *status) {
	*status = PATHPARSE_ERR_OK;
	gboolean has_wildcard = FALSE;
	gchar *target_date = NULL; // for date wildcard in READMODE
	gboolean isdatetime = FALSE; // for date wildcard in READMODE
	if (!g_utf8_strchr(expression, -1, '$')) { // nothing to parse, return original string
		return g_strdup(expression);
	}
	gchar *out = NULL, *localexpression = NULL;
	int nofail = (mode == PATHPARSE_MODE_WRITE_NOFAIL) ? -1 : 1; // statuses in nofail mode are transformed to warnings if multiplied by -1
	if (!fit->header) {
		*status = nofail * PATHPARSE_ERR_HEADER_NULL;
		if (*status > 0) {
			display_path_parse_error(*status, NULL);
			return out;
		}
	}

	// Checking that we don't have 2 date with wildcards, this is not allowed
	if (count_pattern_occurence(expression, "\\$\\*DATE") > 1) {
		*status = nofail * PATHPARSE_ERR_WILDCARD_SYNTAX;
		if (*status > 0) {
			display_path_parse_error(*status, expression);
			return out;
		}
	}

	if (g_str_has_prefix(expression, "$def")) { // using reserved keywords $defbias, $defdark, $defflat, $defstack
		if (!g_strcmp0(expression + 4, "bias")) {
			localexpression = g_strdup(com.pref.prepro.bias_lib);
		} else if (!g_strcmp0(expression + 4, "dark")) {
			localexpression = g_strdup(com.pref.prepro.dark_lib);
		} else if (!g_strcmp0(expression + 4, "flat")) {
			localexpression = g_strdup(com.pref.prepro.flat_lib);
		} else if (!g_strcmp0(expression + 4, "stack")) {
			localexpression = g_strdup(com.pref.prepro.stack_default);
		} else if (!g_strcmp0(expression + 4, "disto")) {
			localexpression = g_strdup(com.pref.prepro.disto_lib);
		} else {
			*status = nofail * PATHPARSE_ERR_WRONG_RESERVED_KEYWORD;
			display_path_parse_error(*status, expression);
			out = (*status > 0) ? NULL : g_strdup(expression); // using libsmthg as a fallback
			g_free(localexpression);
			return out;
		}
		if (!localexpression || strlen(localexpression) == 0) {
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
		gchar *fullheader = g_strconcat(fit->header, "\n", fit->unknown_keys, NULL);
		headerkeys = g_strsplit(fullheader, "\n", 0);
		g_free(fullheader);
	}
	for (int i = 0; i < g_strv_length(tokens); i++) {
		if (!tokens[i] || tokens[i][0] == '\0')
			continue;
		gchar **subs = g_strsplit(tokens[i], ":", 2);
		gchar buf[FLEN_VALUE];
		gchar key[9] = "";
		buf[0] = '\0';
		if (g_strv_length(subs) == 1) {
			if (!g_strcmp0(subs[0], "seqname")) { // reserved keyword $seqname$
				if (sequence_is_loaded() && com.seq.seqname) {
					g_snprintf(buf, FLEN_VALUE - 1, "%s%s", com.seq.seqname,
					g_str_has_suffix(com.seq.seqname, "_") ? "" :
					(g_str_has_suffix(com.seq.seqname, "-") ? "" : "_"));
					g_free(tokens[i]);
					tokens[i] = g_strdup(buf);
					g_strfreev(subs);
					continue;
				} else {
					g_snprintf(key, 9, "%s", subs[0]); // to be used if no fail
					*status = nofail * PATHPARSE_ERR_NOSEQLOADED;
					display_path_parse_error(*status, NULL);
					if (*status > 0) {
						g_strfreev(subs);
						goto free_and_exit;
					}
				}
			}
			g_strfreev(subs);
			continue;
		}
		if (strlen(subs[0]) == 1 && (subs[1][0] == '\\' || subs[1][0] == '/')) {// dealing with Windows drive letter "C:"
			g_strfreev(subs);
			continue;
		}
		// Check if expression starts with a * wildcard
		// Behavior will depend if we are in read or write mode
		gboolean has_wildcard_date = FALSE;
		if (subs[0][0] == '*') {
			if (mode == PATHPARSE_MODE_READ) {
				has_wildcard = TRUE;
				if ((g_str_has_prefix(subs[1],"dm"))) {
					g_snprintf(buf, strlen(DM_PATTERN) + 1, "%s", DM_PATTERN);
					has_wildcard_date = TRUE; // we still need to parse the target date without copying to buf
					g_snprintf(key, 9, "%s", subs[0] + 1);
				} else if ((g_str_has_prefix(subs[1],"dt"))) {
					g_snprintf(buf, strlen(DT_PATTERN) + 1, "%s", DT_PATTERN);
					has_wildcard_date = TRUE; // we still need to parse the target date without copying to buf
					g_snprintf(key, 9, "%s", subs[0] + 1);
					isdatetime = TRUE;
				} else {
					g_snprintf(buf, 2, "%s", "*");
				}
			} else {
				g_snprintf(key, 9, "%s", subs[0] + 1);
			}
		} else {
			g_snprintf(key, 9, "%s", subs[0]);
		}
		if (buf[0] == '*') {
			printf("Wildcard in READ mode\n");
		} else if (subs[1][0] == '%' && (g_str_has_suffix(subs[1], "d") || g_str_has_suffix(subs[1], "f"))) { // case %d or %f
			gboolean isint = g_str_has_suffix(subs[1], "d");
			double val = 0.;
			*status = nofail * read_key_from_header_text(headerkeys, key,&val, NULL);
			display_path_parse_error(*status, key);
			if (*status > 0) {
				g_strfreev(subs);
				goto free_and_exit;
			}
			if (*status == 0) {// can still be neg with write_nofail
				int success = (isint) ? sprintf(buf, subs[1], (int)val) : sprintf(buf, subs[1], val);
				if (success < 0) {// parsing failed probably due to format
					*status = nofail * PATHPARSE_ERR_UNSUPPORTED_FORMAT;
					display_path_parse_error(*status, subs[1]);
					if (*status > 0) {
						g_strfreev(subs);
						goto free_and_exit;
					}
				}
				remove_spaces_from_str(buf); // just in case the formatter introduced spaces
			}
		} else if (subs[1][0] == '%' && g_str_has_suffix(subs[1], "s")) { // case %s
			char val[FLEN_VALUE];
			if (!headerkeys) // ensure null pointer isn't passed to read_key_from_header_text()
				*status = 1;
			else {
				*status = nofail * read_key_from_header_text(headerkeys, key, NULL, val);
				display_path_parse_error(*status, key);
			}
			if (*status > 0) {
				g_strfreev(subs);
				goto free_and_exit;
			}
			if (*status == 0) { // can still be neg with write_nofail
				int success = sprintf(buf, subs[1], val); // just in case there is a fancy formatting directive
				if (success < 0) {// parsing failed probably due to format
					*status = nofail * PATHPARSE_ERR_UNSUPPORTED_FORMAT;
					display_path_parse_error(*status, subs[1]);
					if (*status > 0) {
						g_strfreev(subs);
						goto free_and_exit;
					}
				}
				g_strstrip(buf);
				replace_spaces_from_str(buf, '_');
				replace_invalid_chars(buf ,'_');
			}
		} else if (g_str_has_prefix(subs[1],"dm")) { // case dm12 - date minus 12hrs or dm0
			double minus_hour = -1. * g_ascii_strtod(subs[1] + 2, NULL);
			char val[FLEN_VALUE];
			if (!headerkeys) {
				*status = nofail * PATHPARSE_ERR_HEADER_NULL;
			} else {
				*status = nofail * read_key_from_header_text(headerkeys, key, NULL, val);
				display_path_parse_error(*status, key);
			}
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
					if (has_wildcard_date) {
						target_date = g_strdup(fmtdate);
					} else {
						strncpy(buf, fmtdate, FLEN_VALUE - 1);
					}
					g_date_time_unref(read_time);
					g_date_time_unref(read_time_corr);
					g_free(fmtdate);
				}
			}
		} else if (g_str_has_prefix(subs[1],"dt")) { // case dt (datetime)
			char val[FLEN_VALUE];
			if (!headerkeys) {
				*status = 1;
			} else {
				*status = nofail * read_key_from_header_text(headerkeys, key, NULL, val);
				display_path_parse_error(*status, key);
			}
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
					gchar *fmtdate = date_time_to_date_time(read_time);
					if (has_wildcard_date) {
						target_date = g_strdup(fmtdate);
					} else {
						strncpy(buf, fmtdate, FLEN_VALUE - 1);
					}
					g_date_time_unref(read_time);
					g_free(fmtdate);
				}
			}
		} else if (g_str_has_prefix(subs[1], "ra") || g_str_has_prefix(subs[1], "dec")) { // case ra and dec (str), ran and decn (num)
			gboolean is_ra = g_str_has_prefix(subs[1],"ra");
			gboolean is_float = g_str_has_suffix(subs[1],"n");
			SirilWorldCS *target_coords = NULL;
			char val[FLEN_VALUE];
			double valf = 0.f;
			if (is_float) {
				if (!headerkeys) {
					*status = 1;
				} else {
					*status = nofail * read_key_from_header_text(headerkeys, key,&valf, NULL);
					display_path_parse_error(*status, key);
				}
			} else {
				if (!headerkeys) {
					*status = 1;
				} else {
					*status = nofail * read_key_from_header_text(headerkeys, key, NULL, val);
					display_path_parse_error(*status, key);
				}
			}
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
	if (has_wildcard) {
		gchar *foundmatch = wildcard_check(out, status, target_date, isdatetime);
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

/*
Same as path_parse but makes sure the header is updated before parsing
A copy of the fits metadata is made into a temporary fit that gets updated
and its header string generated.
This temporary fit is then called by path_parse.
If creatdir is TRUE, creates directory if required
*/
gchar *update_header_and_parse(fits *fit, gchar *expression, pathparse_mode mode, gboolean createdir, int *status) {
	*status = PATHPARSE_ERR_OK;
	if (!g_utf8_strchr(expression, -1, '$')) { // nothing to parse, return original string
		return g_strdup(expression);
	}
	gchar *parsedname = NULL, *dirname = NULL;
	update_fits_header(fit);
	parsedname = path_parse(fit, expression, mode, status);
	if (parsedname && createdir) {
		dirname = g_path_get_dirname(parsedname);
		if (siril_mkdir_with_parents(dirname, 0755) < 0) {
			g_free(parsedname);
			parsedname = NULL;
		}
	}
	g_free(dirname);
	return parsedname;
}
