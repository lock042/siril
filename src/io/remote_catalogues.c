/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <unistd.h>
#endif
#endif

#if defined(HAVE_JSON_GLIB) && defined(HAVE_NETWORKING)
#include <json-glib/json-glib.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "algos/search_objects.h"
#include "algos/siril_wcs.h"
#include "io/annotation_catalogues.h"
#include "algos/astrometry_solver.h"
#include "algos/comparison_stars.h"
#include "io/siril_catalogues.h"
#include "io/remote_catalogues.h"
#include "registration/matching/misc.h"

// These statics define the formatting for some fields used when writing catalog names
static const gchar *catcodefmt = "%02d", *rafmt = "%08.4f", *decfmt = "%+08.4f", 
					*radiusfmt = "%3.2f", *limitmagfmt = "%3.1f";

// this function returns a heap-allocated structure cat_tap_query_fields
// cointaining the fields to be queried through TAP (only!)
// to be freed by the caller
// When adding a new catalog that can be queried through TAP:
// - the required query fields, catalogue code and server need to be added
// - it needs to be added to the siril_cat_index enum
// Warning: For Vizier, the catcode needs to be enclosed between ", hence the %22 chars
static cat_tap_query_fields *catalog_to_tap_fields(siril_cat_index cat) {
	cat_tap_query_fields *tap = calloc(1, sizeof(cat_tap_query_fields));
	switch (cat) {
		case CAT_TYCHO2:
			tap->catcode = g_strdup("%22I/259/tyc2%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAmdeg");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEmdeg");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("VTmag");
			tap->tap_columns[CAT_FIELD_BMAG] = g_strdup("BTmag");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("pmRA");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("pmDE");
			break;
		case CAT_NOMAD:
			tap->catcode = g_strdup("%22I/297/out%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("pmRA");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("pmDE");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("Vmag");
			tap->tap_columns[CAT_FIELD_BMAG] = g_strdup("Bmag");
			break;
		case CAT_GAIADR3:
			tap->catcode = g_strdup("%22I/355/gaiadr3%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("pmRA");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("pmDE");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("Gmag");
			tap->tap_columns[CAT_FIELD_BMAG] = g_strdup("BPmag");
			break;
		case CAT_PPMXL:
			tap->catcode = g_strdup("%22I/317/sample%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("pmRA");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("pmDE");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("Jmag");
			break;
		case CAT_BSC:
			tap->catcode = g_strdup("%22V/50/catalog%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("pmRA");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("pmDE");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("Vmag");
			tap->tap_columns[CAT_FIELD_NAME] = g_strdup("Name");
			break;
		case CAT_APASS:
			tap->catcode = g_strdup("%22II/336/apass9%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("Vmag");
			tap->tap_columns[CAT_FIELD_BMAG] = g_strdup("Bmag");
			tap->tap_columns[CAT_FIELD_E_MAG] = g_strdup("e_Vmag");
			tap->tap_columns[CAT_FIELD_E_BMAG] = g_strdup("e_Bmag");
			break;
		case CAT_GCVS:
			tap->catcode = g_strdup("%22B/gcvs/gcvs_cat%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("magMax");
			tap->tap_columns[CAT_FIELD_NAME] = g_strdup("VarName");
			break;
		case CAT_VSX:
			tap->catcode = g_strdup("%22B/vsx/vsx%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("max");
			tap->tap_columns[CAT_FIELD_NAME] = g_strdup("Name");
			break;
		case CAT_SIMBAD:
			tap->catcode = g_strdup("basic+JOIN+allfluxes+ON+oidref+=+oid");
			tap->tap_server = g_strdup(SIMBAD_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("ra");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("dec");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("V");
			tap->tap_columns[CAT_FIELD_BMAG] = g_strdup("B");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("pmra");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("pmdec");
			tap->tap_columns[CAT_FIELD_NAME] = g_strdup("main_id");
			break;
		case CAT_PGC:
			tap->catcode = g_strdup("%22VII/237/pgc%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RAJ2000");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DEJ2000");
			tap->tap_columns[CAT_FIELD_NAME] = g_strdup("PGC");
			tap->tap_columns[CAT_FIELD_DIAMETER] = g_strdup("EXP(logD25)/10"); // diameter in arcmin
			break;
		case CAT_EXOPLANETARCHIVE:
			tap->catcode = g_strdup("pscomppars"); // we query pscomppars instead of ps to get a single entry per planet
			tap->tap_server = g_strdup(EXOPLANETARCHIVE_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("ra");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("dec");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("sy_pmra");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("sy_pmdec");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("sy_vmag");
			tap->tap_columns[CAT_FIELD_NAME] = g_strdup("pl_name");
			break;
		default:
			free(tap);
			return NULL;
	}
	return tap;
}

static void free_cat_tap_query_fields(cat_tap_query_fields *tap) {
	if (!tap)
		return;
	g_free(tap->tap_server);
	g_free(tap->catcode);
	for (int i = 0; i < MAX_TAP_QUERY_COLUMNS; i++) {
		g_free(tap->tap_columns[i]);
	}
	free(tap);
}

/*              _ _                                    _
 *   ___  _ __ | (_)_ __   ___    __ _ _   _  ___ _ __(_) ___  ___
 *  / _ \| '_ \| | | '_ \ / _ \  / _` | | | |/ _ \ '__| |/ _ \/ __|
 * | (_) | | | | | | | | |  __/ | (_| | |_| |  __/ |  | |  __/\__ \
 *  \___/|_| |_|_|_|_| |_|\___|  \__, |\__,_|\___|_|  |_|\___||___/
 *                                  |_|
 * querying servers
 */

/* TODO: fetch_url is also defined in siril update checking, can we merge them? */

#ifdef HAVE_LIBCURL // the alternative is glib-networking, see the else below

static CURL *curl;
static const int DEFAULT_FETCH_RETRIES = 3;

struct ucontent {
	char *data;
	size_t len;
};

static void init() {
	if (!curl) {
		siril_debug_print("initializing CURL\n");
		curl_global_init(CURL_GLOBAL_ALL);
		curl = curl_easy_init();
		if (g_getenv("CURL_CA_BUNDLE"))
			if (curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE")))
				siril_debug_print("Error in curl_easy_setopt()\n");
	}

	if (!curl) {
		fprintf(stderr, "CURL won't initialize\n");
		exit(EXIT_FAILURE);
	}
}

static size_t cbk_curl(void *buffer, size_t size, size_t nmemb, void *userp) {
	size_t realsize = size * nmemb;
	struct ucontent *mem = (struct ucontent *) userp;

	mem->data = realloc(mem->data, mem->len + realsize + 1);

	memcpy(&(mem->data[mem->len]), buffer, realsize);
	mem->len += realsize;
	mem->data[mem->len] = 0;

	return realsize;
}

char *fetch_url(const gchar *url) {
	struct ucontent *content = malloc(sizeof(struct ucontent));
	char *result = NULL;
	long code;
	int retries;
	unsigned int s;

	init();
	retries = DEFAULT_FETCH_RETRIES;

retrieve:
	content->data = malloc(1);
	content->data[0] = '\0';
	content->len = 0;

	CURLcode ret = curl_easy_setopt(curl, CURLOPT_URL, url);
	ret |= curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);
	ret |= curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	ret |= curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	ret |= curl_easy_setopt(curl, CURLOPT_USERAGENT, PACKAGE_STRING);
	if (ret)
		siril_debug_print("Error in curl_easy_setopt()\n");

	siril_debug_print("fetch_url(): %s\n", url);
	if (curl_easy_perform(curl) == CURLE_OK) {
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);

		switch (code) {
		case 200:
			result = content->data;
			siril_debug_print("downloaded %zd bytes\n", content->len);
			break;
		case 500:
		case 502:
		case 503:
		case 504:
			siril_debug_print("Failed to download page %s (error %ld)\n", url, code);
			if (retries) {
				s = 2 * (DEFAULT_FETCH_RETRIES - retries) + 2;
				siril_debug_print("Wait %us before retry\n", s);
				g_usleep(s * 1E6);
				free(content->data);
				retries--;
				goto retrieve;
			} else {
				siril_log_color_message(_("After %ld tries, Server unreachable or unresponsive. (%s)\n"), "salmon", DEFAULT_FETCH_RETRIES, content->data);
			}
			break;
		default:
			break;
		}
	}
	else siril_log_color_message(_("Internet Connection failure.\n"), "red");

	curl_easy_cleanup(curl);
	curl = NULL;
	if (!result || content->len == 0) {
		free(content->data);
		result = NULL;
	}
	free(content);

	return result;
}

void free_fetch_result(char *result) {
	free(result);
}

#elif defined HAVE_NETWORKING

gchar *fetch_url(const gchar *url) {
	GFile *file = g_file_new_for_uri(url);
	GError *error = NULL;
	gchar *content = NULL;

	siril_debug_print("fetch_url(): %s\n", url);

	if (!g_file_load_contents(file, NULL, &content, NULL, NULL, &error)) {
		siril_log_color_message(_("Server unreachable or unresponsive. (%s)\n"), "salmon", error->message);
		g_clear_error(&error);
	}
	g_object_unref(file);
	return content;
}

void free_fetch_result(gchar *result) {
	g_free(result);
}
#endif

// Returns url to be queried based on catalog type and query
static gchar *siril_catalog_conesearch_get_url(siril_catalogue *siril_cat) {
#ifndef HAVE_NETWORKING
	siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
	return NULL;
#endif
	GString *url;
	gchar *fmtstr, *dt;
	const gchar **cat_columns = get_cat_colums_names();
	switch (siril_cat->cat_index){
		/////////////////////////////////////////////////////////////
		// TAP QUERY to csv - preferred way as it requires no parsing
		/////////////////////////////////////////////////////////////
		case CAT_TYCHO2 ... CAT_EXOPLANETARCHIVE:;
			cat_tap_query_fields *fields = catalog_to_tap_fields(siril_cat->cat_index);
			uint32_t catcols = siril_catalog_columns(siril_cat->cat_index);
			url = g_string_new(fields->tap_server);
			gboolean first = TRUE;
			for (int i = 0; i < MAX_TAP_QUERY_COLUMNS; i++) {
				if (fields->tap_columns[i]) {
					g_string_append_printf(url, "%s%s+as+%s", (first) ? "" : ",", fields->tap_columns[i], cat_columns[i]);
					if (first)
						first = FALSE;
				}
			}
			g_string_append_printf(url,"+FROM+%s", fields->catcode);
			g_string_append_printf(url,"+WHERE+CONTAINS(POINT('ICRS',%s,%s),", fields->tap_columns[CAT_FIELD_RA], fields->tap_columns[CAT_FIELD_DEC]);
			fmtstr = g_strdup_printf("CIRCLE('ICRS',%s,%s,%s))=1", rafmt, decfmt, radiusfmt);
			g_string_append_printf(url, fmtstr, siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius / 60.);
			g_free(fmtstr);
			if (siril_cat->limitmag > 0 && catcols & (1 << CAT_FIELD_MAG)) {
				fmtstr = g_strdup_printf("+AND+(%%s<=%s)", limitmagfmt);
				g_string_append_printf(url, fmtstr,  fields->tap_columns[CAT_FIELD_MAG], siril_cat->limitmag);
				g_free(fmtstr);
			}
			free_cat_tap_query_fields(fields);
			return g_string_free(url, FALSE);
		//////////////////////////////////
		// AAVSO chart of comparison stars
		//////////////////////////////////
		case CAT_AAVSO_CHART:
			url = g_string_new(AAVSOCHART_QUERY);
			fmtstr = g_strdup_printf("&ra=%s&dec=%s&fov=%s&maglimit=%s", rafmt, decfmt, radiusfmt, limitmagfmt);
			g_string_append_printf(url, fmtstr, siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius, siril_cat->limitmag);
			g_free(fmtstr);
			return g_string_free(url, FALSE);
		////////////////////////////
		// IMCCE - skybot conesearch
		////////////////////////////
		case CAT_IMCCE:
			if (!siril_cat->dateobs) { // should have been caught before... just in case
				siril_log_color_message(_("This command only works on images that have observation date information\n"), "red");
				return NULL;
			}
			url = g_string_new(IMCCE_QUERY);
			dt = date_time_to_FITS_date(siril_cat->dateobs);
			g_string_append_printf(url,"&-ep=%s", dt);
			g_free(dt);
			fmtstr = g_strdup_printf("&-ra=%s&-dec=%s&-rd=%s", rafmt, decfmt, radiusfmt);
			g_string_append_printf(url, fmtstr, siril_cat->center_ra, siril_cat->center_dec, 2. * siril_cat->radius / 60.); // request uses diameter not radius (despite what's written in the doc)
			g_free(fmtstr);
			g_string_append_printf(url,"&-loc=%s", (siril_cat->IAUcode) ? siril_cat->IAUcode : "500");
			return g_string_free(url, FALSE);
		default:
			break;
	}
	return NULL;
}

// Parsers for non-TAP queries
// IMCCE skybot - txt input delimited with " | " parsed to csv
static gboolean parse_IMCCE_buffer(gchar *buffer, GOutputStream *output_stream) {
	if (!buffer || buffer[0] == '\0' || !g_str_has_prefix(buffer, "# Flag:"))
		return FALSE;
	if (!g_str_has_prefix(buffer, "# Flag: 1") && !g_str_has_prefix(buffer, "# Flag: 0")) {
		siril_log_color_message(_("IMCCE server returned:\n"), "red");
		gchar **err_lines = g_strsplit(buffer, "\n", -1);
		int n_err = min(g_strv_length(err_lines), 3); // displaying max 3 lines
		for (int i = 0; i < n_err; i++) {
			if (err_lines[i][0] != '\0')
				siril_log_color_message("%s\n", "red", err_lines[i]);
		}
		g_strfreev(err_lines);
		return FALSE;
	}
	// we fill a siril_catalogue struct to use the generic catalogue writer
	// may be a bit slower than direct write to the output_stream but will ease maintenance
	// anyway, those catalogs are usually small so little impact on performance is to be expected
	siril_catalogue *siril_cat = calloc(1, sizeof(siril_catalogue));
	siril_cat->cat_index = CAT_IMCCE;
	siril_cat->columns = siril_catalog_columns(CAT_IMCCE);
	gchar **token = g_strsplit(buffer, "\n", -1);
	int nb_lines = g_strv_length(token);
	int nstars = nb_lines - 3;
	if (g_str_has_prefix(buffer, "# Flag: 0") || !nstars) {
		siril_debug_print("No items in the IMCCE catalog, will just write header\n");
		nstars = 0;
	}
	cat_item *cat_items = NULL;
	int n = 0;
	for (int i = 3; i < nb_lines; i++) {
		if (i == 3) {
			cat_items = calloc(nstars, sizeof(cat_item));
			if (!cat_items) {
				PRINT_ALLOC_ERR;
				siril_catalog_free(siril_cat);
				return FALSE;
			}
		}
		// format is '# Num | Name | RA(h) | DE(deg) | Class | Mv | Err(arcsec) | d(arcsec) | dRA(arcsec/h) | dDEC(arcsec/h) | Dg(ua) | Dh(ua)'
		gchar **vals = g_strsplit(token[i], " | ", -1);
		if (g_strv_length(vals) != 12) {
			g_strfreev(vals);
			continue;
		}
		double ra = parse_hms(vals[2]);	// in hours
		double dec = parse_dms(vals[3]);
		if (!isnan(ra) && !isnan(dec)) {
			cat_items[n].ra = ra;
			cat_items[n].dec = dec;
			cat_items[n].mag = (float)g_strtod(vals[5], NULL);
			cat_items[n].name = g_strdup(vals[1]);
			cat_items[n].vra = (float)g_strtod(vals[8], NULL);
			cat_items[n].vdec = (float)g_strtod(vals[9], NULL);
			cat_items[n].type = g_strdup(vals[4]);
			n++;
		}
		g_strfreev(vals);
	}
	g_strfreev(token);
	if (nstars && n < nstars) {
		if (!n) {
			free(cat_items);
			cat_items = NULL;
		} else {
			cat_item *new_array = realloc(cat_items, n * sizeof(cat_item));
			if (!new_array) {
				PRINT_ALLOC_ERR;
				siril_catalog_free(siril_cat);
				return FALSE;
			}
			cat_items = new_array;
		}
	}
	siril_cat->cat_items = cat_items;
	siril_cat->nbitems = n;
	gboolean ret = siril_catalog_write_to_output_stream(siril_cat, output_stream);
	siril_catalog_free(siril_cat);
	return ret;
}

// AAVSO chart - json input read with json-glib parsed to csv
static gboolean parse_AAVSO_Chart_buffer(gchar *buffer, GOutputStream *output_stream) {
#ifndef HAVE_JSON_GLIB
	siril_log_color_message(_("json-glib was not found at build time, cannot proceed. Install and rebuild.\n"), "red");
	return FALSE;
#else
	GError *error = NULL;
	JsonParser *parser = json_parser_new();
	if (!json_parser_load_from_data(parser, buffer, -1, &error)) {
		siril_log_color_message(_("Could not parse AAVSO chart buffer: %s\n"), "red", error->message);
		g_clear_object(&parser);
		g_clear_error(&error);
		return FALSE;
	}
	// we fill a siril_catalogue struct to use the generic catalogue writer
	// may be a bit slower than direct write to the output_stream but will ease maintenance
	// anyway, those catalogs are usually small so little impact on performance is to be expected
	siril_catalogue *siril_cat = calloc(1, sizeof(siril_catalogue));
	siril_cat->cat_index = CAT_AAVSO_CHART;
	siril_cat->columns = siril_catalog_columns(CAT_AAVSO_CHART);

	// parsing the AAVSO chart id in the comments section
	JsonReader *reader = json_reader_new(json_parser_get_root(parser));
	json_reader_read_member(reader, "chartid");
	const gchar *id = json_reader_get_string_value(reader);
	json_reader_end_member(reader);
	siril_cat->header = g_strdup_printf("#ChartID:%s", id);

	json_reader_read_member(reader, "photometry");
	int nstars = json_reader_count_elements(reader);
	int n = 0;
	cat_item *cat_items = NULL;
	for (int i = 0; i < nstars; i++) {
		if (i == 0) {
			cat_items = calloc(nstars, sizeof(cat_item));
		}
		json_reader_read_element(reader, i);
		// auid
		json_reader_read_member(reader, "auid");
		const gchar *name = json_reader_get_string_value(reader);
		json_reader_end_member(reader);
		// reading the data or V and B bands, bands being an array
		json_reader_read_member(reader, "bands");
		int nbands = json_reader_count_elements(reader);
		const gchar *band = NULL;
		double mag = 0., e_mag = 0., bmag = 0., e_bmag = 0.;
		for (int j = 0; j < nbands; j++) {
			json_reader_read_element(reader, j);
			json_reader_read_member(reader, "band");
			band = json_reader_get_string_value(reader);
			json_reader_end_member(reader);
			if (band && !strcmp(band, "V")) {
				json_reader_read_member(reader, "error");
				e_mag = json_reader_get_double_value(reader);
				json_reader_end_member(reader);
				json_reader_read_member(reader, "mag");
				mag = json_reader_get_double_value(reader);
				json_reader_end_member(reader);
			} else if (band && !strcmp(band, "B")) {
				json_reader_read_member(reader, "error");
				e_bmag = json_reader_get_double_value(reader);
				json_reader_end_member(reader);
				json_reader_read_member(reader, "mag");
				bmag = json_reader_get_double_value(reader);
				json_reader_end_member(reader);
			}
			json_reader_end_element(reader);
		}
		json_reader_end_member(reader);
		//dec
		json_reader_read_member(reader, "dec");
		const gchar *decstr = json_reader_get_string_value(reader);
		double dec = parse_dms(decstr);
		json_reader_end_member(reader);
		//ra
		json_reader_read_member(reader, "ra");
		const gchar *rastr = json_reader_get_string_value(reader);
		double ra = parse_hms(rastr);	// in hours
		json_reader_end_member(reader);
		if (!isnan(ra) && !isnan(dec)) {
			cat_items[n].ra = ra;
			cat_items[n].dec = dec;
			cat_items[n].mag = mag;
			cat_items[n].bmag = bmag;
			cat_items[n].e_mag = e_mag;
			cat_items[n].e_bmag = e_bmag;
			cat_items[n].name = g_strdup(name);
			n++;
		}
		json_reader_end_element(reader);
	}
	g_object_unref(reader);
	g_object_unref(parser);
	if (nstars && n < nstars) {
		if (!n) {
			free(cat_items);
			cat_items = NULL;
		} else {
			cat_item *new_array = realloc(cat_items, n * sizeof(cat_item));
			if (!new_array) {
				PRINT_ALLOC_ERR;
				siril_catalog_free(siril_cat);
				return FALSE;
			}
			cat_items = new_array;
		}
	}
	siril_cat->cat_items = cat_items;
	siril_cat->nbitems = n;
	gboolean ret = siril_catalog_write_to_output_stream(siril_cat, output_stream);
	siril_catalog_free(siril_cat);
	return ret;
#endif
}

static gchar *parse_remote_catalogue_filename(siril_catalogue *siril_cat) {
	gchar *filename = NULL;
	gchar *dt = NULL, *fmtstr = NULL;
	switch (siril_cat->cat_index) {
		case CAT_TYCHO2 ... CAT_AAVSO_CHART:
			fmtstr = g_strdup_printf("cat_%s_%s_%s_%s_%s.csv", catcodefmt, rafmt, decfmt, radiusfmt, limitmagfmt);
			filename = g_strdup_printf(fmtstr,
				(int)siril_cat->cat_index,
				siril_cat->center_ra,
				siril_cat->center_dec,
				siril_cat->radius,
				siril_cat->limitmag);
			g_free(fmtstr);
			return filename;
			break;
		case CAT_IMCCE:
			if (!siril_cat->IAUcode || !siril_cat->dateobs) {
				siril_debug_print("Queries for solar system should pass date and location code\n");
				return NULL;
			}
			dt = date_time_to_date_time(siril_cat->dateobs);
			fmtstr = g_strdup_printf("cat_%s_%s_%s_%s_%%s_%%s.csv", catcodefmt, rafmt, decfmt, radiusfmt);
			filename = g_strdup_printf(fmtstr,
				(int)siril_cat->cat_index,
				siril_cat->center_ra,
				siril_cat->center_dec,
				siril_cat->radius,
				dt,
				siril_cat->IAUcode);
			g_free(fmtstr);
			g_free(dt);
			return filename;
			break;
		default:
			siril_debug_print("should not happen, download_catalog should only be called for online calls\n");
			return NULL;
	}
	siril_debug_print("catalog type not handled, should not happen\n");
	return NULL;
}

// Parses the catalogue name and checks if it exists in cache
// Returns the path to write to (if in_cache is FALSE)
// or the path to read directly (if in_cache is TRUE)
static gchar *get_remote_catalogue_cached_path(siril_catalogue *siril_cat, gboolean *in_cache) {
	GError *error = NULL;
	GFile *file = NULL;
	*in_cache = FALSE;
	gchar *filename = parse_remote_catalogue_filename(siril_cat);
	if (!filename)
		return NULL;
	siril_debug_print("Catalogue file: %s\n", filename);

	// check if download_cache folder exists, create it otherwise
	gchar *root = g_build_filename(siril_get_config_dir(), PACKAGE, "download_cache", NULL);
	gchar *filepath = g_build_filename(root, filename, NULL);

	if (!g_file_test(root, G_FILE_TEST_EXISTS)) {
		if (g_mkdir_with_parents(root, 0755) < 0) {
			siril_log_color_message(_("Cannot create output folder: %s\n"), "red", root);
			g_free(filepath);
			g_free(root);
			return NULL; // we won't be able to write to the file
		}
	}
	g_free(root);

	if (g_file_test(filepath, G_FILE_TEST_EXISTS)) { // file already exists in cache
		file = g_file_new_for_path(filepath);
		GFileInfo *info = g_file_query_info(file, G_FILE_ATTRIBUTE_TIME_MODIFIED "," G_FILE_ATTRIBUTE_STANDARD_SIZE, 0, NULL, NULL);
		if ((g_file_info_get_size(info)) == 0) { // test if not empty and delete in case it is
			if (!g_file_delete(file, NULL, &error)) {
				siril_log_color_message(_("A corrupted version of %s was found in cache but cannot be deleted, aborting\n"), "red", filepath);
				g_free(filepath);
				return NULL;
			}
			return filepath; // the corrupted version was successfully deleted, passing the filepath
		} else {
			*in_cache = TRUE;
			return filepath;
		}
	}
	// TODO: expand cache search to find catalogs matching query within certain criteria:
	// - same catalog type
	// - sufficient overlap of the sky region
	// - large enough mag
	// - valid date and sitelocation (for sso)
	return filepath;
}

/* Downloads and writes to download_cache (if required) the online catalogue
   as per given catalogue type, center, radius, limit mag (optionnaly obscode and date obs for sso)
   Returns the path to the file (whether already cached or downloaded)
*/
static gchar *download_catalog(siril_catalogue *siril_cat) {
#ifndef HAVE_NETWORKING
	siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
#else
	gchar *str = NULL, *filepath = NULL, *url = NULL, *buffer = NULL;
	GError *error = NULL;
	GOutputStream *output_stream = NULL;
	GFile *file = NULL;
	gboolean remove_file = FALSE, catalog_is_in_cache = FALSE;

	/* check if catalogue already exists in cache */
	filepath = get_remote_catalogue_cached_path(siril_cat, &catalog_is_in_cache);
	g_free(str);

	if (catalog_is_in_cache) {
		siril_log_message(_("Using already downloaded catalogue %s\n"), catalog_to_str(siril_cat->cat_index));
		return filepath;
	}
	if (!filepath) { // if the path is NULL, an error was caught earlier, just free and abort
		g_free(str);
		return NULL;
	}

	// the catalog needs to be downloaded, prepare the output stream
	file = g_file_new_for_path(filepath);
	output_stream = (GOutputStream*) g_file_create(file, G_FILE_CREATE_NONE, NULL, &error);
	if (!output_stream) {
		goto download_error;
	}

	/* download */
	url = siril_catalog_conesearch_get_url(siril_cat);
	if (!url) {
		remove_file = TRUE;
		goto download_error;
	}
	siril_log_message(_("Contacting server\n"));
	buffer = fetch_url(url);
	g_free(url);

	/* save (and parse if required)*/
	if (buffer) {
		switch (siril_cat->cat_index) {
			case CAT_TYCHO2 ... CAT_EXOPLANETARCHIVE: // TAP query, no parsing, we just write the whole buffer to the output stream
				if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
					g_warning("%s\n", error->message);
					remove_file = TRUE;
					goto download_error;
				}
				break;
			case CAT_IMCCE:
				if (!parse_IMCCE_buffer(buffer, output_stream)) {
					remove_file = TRUE;
					goto download_error;
				}
				break;
			case CAT_AAVSO_CHART:
				if (!parse_AAVSO_Chart_buffer(buffer, output_stream)) {
					remove_file = TRUE;
					goto download_error;
				}
				break;
			default:
				break;
		}
		g_object_unref(output_stream);
		g_free(buffer);
	} else { // remove the file from cache
		remove_file = TRUE;
		goto download_error;
	}
	return filepath;

download_error:
	if (error) {
		siril_log_color_message(_("Cannot create catalogue file %s (%s)\n"), "red", filepath, error->message);
		g_clear_error(&error);
		}
	g_free(buffer);
	if (output_stream)
		g_object_unref(output_stream);
	if (file) {
		if (remove_file)
			if (g_unlink(g_file_peek_path(file)))
				siril_debug_print(("Cannot delete catalogue file %s\n"), filepath);
		g_object_unref(file);
	}
	if (filepath)
		g_free(filepath);
#endif
	return NULL;
}

/* This function is the main interface to collect an online catalogue
   It sends a conesearch around given center, within given radius and for stars below limit_mag
   Internally, it uses download_catalog to search cache and download catalog as required
   The cache folder (named download_cache) stores all downloaded queries in the form 
   of csv files, named as 'cat-cat_index-ra-dec-radius[-mag].csv' or
   'cat-cat_index-ra-dec-radius-date-obscode.csv' for solar syatem queries (IMCCE)
   It fills the siril_catalogue given in input
   Returns the number of stars fetched, -1 if successful but empty, 0 otherwise
*/
int siril_catalog_get_stars_from_online_catalogues(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return 0;
	if (siril_cat->cat_index >= CAT_AN_MESSIER) {
		siril_debug_print("Online cat query - Should not happen\n");
		return 0;
	}
	gchar *catfile = download_catalog(siril_cat);
	if (!catfile)
		return 0;
	int retval = siril_catalog_load_from_file(siril_cat, catfile);
	if (!retval)
		return siril_cat->nbitems;
	if (retval == -1)
		return -1; // empty but not failed
	return 0;
}

