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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* useful if no libcurl */
#include <stdlib.h>

#include "yyjson.h"

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_networking.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "algos/astrometry_solver.h"
#include "io/siril_catalogues.h"
#include "io/remote_catalogues.h"

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
			tap->tap_columns[CAT_FIELD_TEFF] = g_strdup("Teff");
			break;
		case CAT_GAIADR3_DIRECT:
			tap->catcode = g_strdup("gaiadr3.gaia_source");
			tap->tap_server = g_strdup(GAIA_DR3_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("ra");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("dec");
			tap->tap_columns[CAT_FIELD_PMRA] = g_strdup("pmra");
			tap->tap_columns[CAT_FIELD_PMDEC] = g_strdup("pmdec");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("phot_g_mean_mag");
			tap->tap_columns[CAT_FIELD_BMAG] = g_strdup("phot_bp_mean_mag");
			tap->tap_columns[CAT_FIELD_TEFF] = g_strdup("teff_gspphot");
			tap->tap_columns[CAT_FIELD_GAIASOURCEID] = g_strdup("source_id");
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
		case CAT_VARISUM:
			tap->catcode = g_strdup("%22I/358/varisum%22");
			tap->tap_server = g_strdup(VIZIER_TAP_QUERY);
			tap->tap_columns[CAT_FIELD_RA] = g_strdup("RA_ICRS");
			tap->tap_columns[CAT_FIELD_DEC] = g_strdup("DE_ICRS");
			tap->tap_columns[CAT_FIELD_MAG] = g_strdup("Gmagmax");
			tap->tap_columns[CAT_FIELD_NAME] = g_strdup("Source");
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

/// Returns url to be queried based on catalog type and query
static gchar *siril_catalog_conesearch_get_url(siril_catalogue *siril_cat) {
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
			g_string_append_printf(url, fmtstr, siril_cat->center_ra, siril_cat->center_dec, 2. * siril_cat->radius, siril_cat->limitmag);
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
			g_string_append_printf(url,"&-observer=%s", (siril_cat->IAUcode) ? siril_cat->IAUcode : "500");
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
	siril_catalogue *siril_cat = siril_catalog_new(CAT_IMCCE);
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

static gboolean parse_AAVSO_Chart_buffer(gchar *buffer, GOutputStream *output_stream) {
	// Parse JSON data
	yyjson_read_err err = { 0 };
	yyjson_doc *doc = yyjson_read(buffer, strlen(buffer), YYJSON_READ_NOFLAG);
	if (!doc) {
		siril_log_color_message(_("Could not parse AAVSO chart buffer: %s\n"), "red", err.msg);
		return FALSE;
	}

	// we fill a siril_catalogue struct to use the generic catalogue writer
	// may be a bit slower than direct write to the output_stream but will ease maintenance
	// anyway, those catalogs are usually small so little impact on performance is to be expected
	siril_catalogue *siril_cat = siril_catalog_new(CAT_AAVSO_CHART);

	// Get root object and parse chart ID
	yyjson_val *root = yyjson_doc_get_root(doc);
	yyjson_val *chartid = yyjson_obj_get(root, "chartid");
	const char *id = yyjson_get_str(chartid);
	siril_cat->header = g_strdup_printf("#ChartID:%s", id);

	// Get photometry array
	yyjson_val *photometry = yyjson_obj_get(root, "photometry");
	size_t nstars = yyjson_arr_size(photometry);
	int n = 0;
	cat_item *cat_items = NULL;
	if (nstars > 0) {
		cat_items = calloc(nstars, sizeof(cat_item));
	}

	// Iterate through stars
	yyjson_val *star;
	size_t idx, max;
	yyjson_arr_foreach(photometry, idx, max, star) {
		const char *name = NULL;
		double ra = NAN, dec = NAN;
		double mag = 0., e_mag = 0., bmag = 0., e_bmag = 0.;

		// Get basic star information
		yyjson_val *auid = yyjson_obj_get(star, "auid");
		if (auid) name = yyjson_get_str(auid);
		yyjson_val *ra_val = yyjson_obj_get(star, "ra");
		if (ra_val) ra = parse_hms(yyjson_get_str(ra_val));  // in hours
		yyjson_val *dec_val = yyjson_obj_get(star, "dec");
		if (dec_val) dec = parse_dms(yyjson_get_str(dec_val));

		// Process bands array
		yyjson_val *bands = yyjson_obj_get(star, "bands");
		if (bands) {
			yyjson_val *band_obj;
			size_t bidx, bmax;  // Separate iteration variables for bands
			yyjson_arr_foreach(bands, bidx, bmax, band_obj) {
				yyjson_val *band_val = yyjson_obj_get(band_obj, "band");
				const char *band = yyjson_get_str(band_val);

				if (band && !strcmp(band, "V")) {
					yyjson_val *mag_val = yyjson_obj_get(band_obj, "mag");
					if (mag_val) {
						mag = (double) yyjson_get_num(mag_val);
					}
					yyjson_val *err_val = yyjson_obj_get(band_obj, "error");
					if (err_val && !yyjson_is_null(err_val)) {
						e_mag = (double) yyjson_get_num(err_val);
					}
				}
				else if (band && !strcmp(band, "B")) {
					yyjson_val *mag_val = yyjson_obj_get(band_obj, "mag");
					if (mag_val) {
						bmag = (double) yyjson_get_num(mag_val);
					}
					yyjson_val *err_val = yyjson_obj_get(band_obj, "error");
					if (err_val && !yyjson_is_null(err_val)) {
						e_bmag = (double) yyjson_get_num(err_val);
					}
				}
			}
		}

		// Add star to catalog if we have valid coordinates
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
	}

	// Free JSON document
	yyjson_doc_free(doc);

	// Resize array if needed
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

	// Set catalog items and write
	siril_cat->cat_items = cat_items;
	siril_cat->nbitems = n;
	gboolean ret = siril_catalog_write_to_output_stream(siril_cat, output_stream);
	siril_catalog_free(siril_cat);
	return ret;
}

static gchar *parse_remote_catalogue_filename(siril_catalogue *siril_cat, retrieval_type datalink_product) {
	gchar *filename = NULL;
	gchar *dt = NULL, *fmtstr = NULL;
	gchar *ext = NULL;
	if (datalink_product == NO_DATALINK_RETRIEVAL || siril_cat->cat_index != CAT_GAIADR3_DIRECT) {
		ext = g_strdup(".csv");
	} else {
		ext = g_strdup_printf("_%d.fit", datalink_product);
	}
	switch (siril_cat->cat_index) {
		case CAT_TYCHO2 ... CAT_AAVSO_CHART:
			fmtstr = g_strdup_printf("cat_%s_%s_%s_%s_%s%s", catcodefmt, rafmt, decfmt, radiusfmt, limitmagfmt, ext);
			filename = g_strdup_printf(fmtstr,
				(int)siril_cat->cat_index,
				siril_cat->center_ra,
				siril_cat->center_dec,
				siril_cat->radius,
				siril_cat->limitmag);
			g_free(fmtstr);
			g_free(ext);
			return filename;
			break;
		case CAT_IMCCE:
			if (!siril_cat->IAUcode || !siril_cat->dateobs) {
				siril_debug_print("Queries for solar system should pass date and location code\n");
				g_free(ext);
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
			g_free(ext);
			return filename;
			break;
		default:
			siril_debug_print("should not happen, download_catalog should only be called for online calls\n");
			g_free(ext);
			break;
	}
	return NULL;
}

// Parses the catalogue name and checks if it exists in cache
// Returns the path to write to (if in_cache is FALSE)
// or the path to read directly (if in_cache is TRUE)
static gchar *get_remote_catalogue_cached_path(siril_catalogue *siril_cat, gboolean *in_cache, retrieval_type datalink_product) {
	GError *error = NULL;
	GFile *file = NULL;
	*in_cache = FALSE;
	gchar *filename = parse_remote_catalogue_filename(siril_cat, datalink_product);
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
	gchar *str = NULL, *filepath = NULL, *url = NULL, *buffer = NULL;
	GError *error = NULL;
	GOutputStream *output_stream = NULL;
	GFile *file = NULL;
	gboolean remove_file = FALSE, catalog_is_in_cache = FALSE;
	int fetch_url_error = 0;

	/* check if catalogue already exists in cache */
	filepath = get_remote_catalogue_cached_path(siril_cat, &catalog_is_in_cache, NO_DATALINK_RETRIEVAL);
	g_free(str);

	if (catalog_is_in_cache) {
		siril_log_message(_("Using already downloaded catalogue %s\n"), catalog_to_str(siril_cat->cat_index));
		return filepath;
	}
	if (!filepath) { // if the path is NULL, an error was caught earlier, just free and abort
		g_free(str);
		return NULL;
	}
	if (!is_online()) {
		siril_log_message(_("Offline: cannot download catalog\n"));
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
	siril_debug_print("URL: %s\n", url);
	siril_log_message(_("Contacting server\n"));
	gsize length;
	buffer = fetch_url(url, &length, &fetch_url_error, FALSE);

	/* save (and parse if required)*/
	if (buffer && !error) {
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

	g_free(url);
	return filepath;

download_error:
	if (fetch_url_error) {
		siril_log_color_message(_("Error: unable to retrieve the remote catalogue from the server at %s. "
			"This is not a Siril bug. It means the catalogue server is unavailable. This is usually a "
			"short-lived problem however this server is not affiliated with Siril and the Siril team do "
			"not control it. We highly recommend using a local server such as the optimized Siril extract "
			"from Gaia DR3, installable using the Catalog_Installer.py script. This is faster and does not "
			"rely on third party servers.\n"), "red", url);
	}
	g_free(url);

	if (error) {
		siril_log_color_message(_("Cannot create catalogue file %s (%s)\n"), "red", filepath, error->message);
		g_clear_error(&error);
	} else {
		siril_log_color_message(_("Cannot create catalogue file %s (generic error)\n"), "red", filepath);
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

static int submit_async_request(const char *url, const char *post_data, char **job_id) {
	if (!is_online()) {
		siril_log_message(_("Offline, cannot retrieve URL.\n"));
		return 1;
	}
	gchar *post_response = NULL;
	int res = submit_post_request(url, post_data, &post_response);
//	siril_debug_print(post_response);
	if (res) return 1;
	gchar *location_header = strstr(post_response, "jobId");
	if (location_header != NULL) {
		const gchar *start = find_first_numeric(location_header);
		gchar *end = strchr(start, ']');
		if (end != NULL) {
			*end = '\0';
			*job_id = g_strdup(start);
			siril_debug_print("Job ID: %s\n", *job_id);
		}
	}
	g_free(post_response);
	return 0;
}

// Returns url to be submitted as a Gaia DR3 job. This allows the job ID to be used as
// the source of source_ids for a subsequent datalink query. This is intended for use
// with SPCC but may be useful in the future for other types of datalink query.
// The path to the datalink FITS is set and the catalogue is populated
int siril_gaiadr3_datalink_query(siril_catalogue *siril_cat, retrieval_type type, gchar** datalink_path, int max_datalink_sources) {
	if (!(siril_compiled_with_networking())) {
		siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
		return -1;
	} else if (!(is_online())) {
		siril_log_color_message(_("Offline: cannot fetch catalog\n"), "red");
		return -1;
	}
	GError *error = NULL;
	gsize length;
	const gchar *host = "https://gea.esac.esa.int";
	GString *querystring = NULL;
	const gchar *pathinfo = "/tap-server/tap/async/";
	gchar *str = NULL, *url = NULL, *data = NULL, *buffer = NULL;
	GString *datalink_url = NULL;
	GOutputStream *output_stream = NULL;
	GFile *file = NULL;
	gboolean remove_file = FALSE;
	gchar *job_id = NULL;

	gboolean catalog_is_in_cache, retrieval_product_is_in_cache;
	gchar *csvfilepath = get_remote_catalogue_cached_path(siril_cat, &catalog_is_in_cache, NO_DATALINK_RETRIEVAL);
	if (!csvfilepath) { // if the path is NULL, an error was caught earlier, just free and abort
		return -1;
	}
	gchar *filepath = get_remote_catalogue_cached_path(siril_cat, &retrieval_product_is_in_cache, type);
	if (!filepath) { // if the path is NULL, an error was caught earlier, just free and abort
		return -1;
	}

	if (!(catalog_is_in_cache && retrieval_product_is_in_cache)) {

		// Set up query
		gchar *fmtstr;
		max_datalink_sources = max(min(max_datalink_sources,5000), 1); // Limit the maximum number of sources to retrieve from Gaia to be between 1-5000.
		const gchar **cat_columns = get_cat_colums_names();
		cat_tap_query_fields *fields = catalog_to_tap_fields(siril_cat->cat_index);
		uint32_t catcols = siril_catalog_columns(siril_cat->cat_index);
		querystring = g_string_new("LANG=ADQL&FORMAT=csv&QUERY=SELECT+TOP+"); // We ignore the cat_server as the URL is dealt with elsewhere
		g_string_append_printf(querystring, "%d+", max_datalink_sources);
		gboolean first = TRUE;
		for (int i = 0; i < MAX_TAP_QUERY_COLUMNS; i++) {
			if (fields->tap_columns[i]) {
				g_string_append_printf(querystring, "%s%s+as+%s", (first) ? "" : ",", fields->tap_columns[i], cat_columns[i]);
				if (first)
					first = FALSE;
			}
		}
		g_string_append_printf(querystring,"+FROM+%s", fields->catcode);
		g_string_append_printf(querystring,"+WHERE+has_xp_sampled+=+'True'+AND+CONTAINS(POINT('ICRS',%s,%s),", fields->tap_columns[CAT_FIELD_RA], fields->tap_columns[CAT_FIELD_DEC]);
		fmtstr = g_strdup_printf("CIRCLE('ICRS',%s,%s,%s))=1", rafmt, decfmt, radiusfmt);
		g_string_append_printf(querystring, fmtstr, siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius / 60.);
		g_free(fmtstr);
		if (siril_cat->limitmag > 0 && catcols & (1 << CAT_FIELD_MAG)) {
			fmtstr = g_strdup_printf("+AND+(%%s<=%s)", limitmagfmt);
			g_string_append_printf(querystring, fmtstr,  fields->tap_columns[CAT_FIELD_MAG], siril_cat->limitmag);
			g_free(fmtstr);
		}
		g_string_append_printf(querystring, "+ORDER+BY+random_index"); // Avoids bias in the results by ordering by random_index
		free_cat_tap_query_fields(fields);

		// Create job
		url = g_strdup_printf("%s%s", host, pathinfo);
		data = g_strdup_printf(
			"PHASE=run&"
			"LANG=ADQL&"
			"FORMAT=csv&"
			"REQUEST=doQuery&"
			"%s", querystring->str
		);
		siril_debug_print("Query data: %s\n", data);
		siril_log_message(_("Submitting conesearch request to ESA Gaia DR3 catalog. This may take a few seconds to complete...\n"));
		if (submit_async_request(url, data, &job_id)) {
			siril_log_color_message(_("Error submitting conesearch request.\n"), "red");
			goto tap_error_and_cleanup;
		}

		// Print job id
		siril_debug_print("Gaia DR3 Job ID: %s\n", job_id);
		if (job_id == NULL) {
			siril_log_color_message(_("Job id not found in the response.\n"), "red");
			goto tap_error_and_cleanup;
		}

		// Wait until the job is finished
		// Possible IVOA UWS statuses are: PENDING, QUEUED, EXECUTING, COMPLETED, ERROR, ABORTED
		//Timeout 1 minute in us
		gchar* job_check = g_strdup_printf("https://gea.esac.esa.int/tap-server/tap/async/%s", job_id);
		uint64_t timer = 0;
		gboolean success = FALSE;
		while (1) {
			if (timer > ASYNC_JOB_TIMEOUT) { // Avoid infinite loop
				siril_log_color_message(_("Timeout on Gaia DR3 query\n"), "red");
				break;
			}
			int fetch_url_error;
			buffer = fetch_url(job_check, &length, &fetch_url_error, FALSE);
			if (fetch_url_error) {
				g_free(buffer);
				break;
			}
			gboolean error = (g_strrstr(buffer, "ERROR") != NULL);
			if (!error)
				error = (g_strrstr(buffer, "ABORTED") != NULL);
			if (error) {
				siril_log_color_message(_("Gaia DR3 async query failed, unable to continue.\n"), "red");
				break;
			}
			gboolean completed = (g_strrstr(buffer,"COMPLETED") != NULL);
			g_free(buffer);
			if (completed) {
				success = TRUE;
				break;
			}
			g_usleep(500000);
			timer += 500000;
		}
		g_free(job_check);

		if (!success) // ERROR, ABORTED or timed out
			goto tap_error_and_cleanup;

		// Retrieve the TAP+ query result
		gchar *job_retrieval = g_strdup_printf("https://gea.esac.esa.int/tap-server/tap/async/%s/results/result", job_id);
		int fetch_url_error;
		gchar *buffer = fetch_url(job_retrieval, &length, &fetch_url_error, FALSE);
		siril_debug_print("buffer length: %lu\n", length);
		g_free(job_retrieval);

		// buffer is the CSV data for the standard Gaia DR3 TAP+ query, it gets saved to the usual catalogue location
		/* check if catalogue already exists in cache */
		GOutputStream *csvoutput_stream = NULL;
		GFile *csvfile = NULL;

		// the catalog needs to be downloaded, prepare the output stream
		g_unlink(csvfilepath); // If no file exists this call may fail, that's fine
		csvfile = g_file_new_for_path(csvfilepath);
		csvoutput_stream = (GOutputStream*) g_file_create(csvfile, G_FILE_CREATE_NONE, NULL, &error);
		if (!csvoutput_stream || !buffer || length == 0) {
			goto tap_error_and_cleanup;
		}
		if (!g_output_stream_write_all(csvoutput_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			remove_file = TRUE;
			goto tap_error_and_cleanup;
		}
		siril_log_message(_("Gaia DR3 conesearch query succeeded: cached as %s\n"), csvfilepath);

		// Populate the siril_catalog with the data from the initial query
		int retval = siril_catalog_load_from_file(siril_cat, csvfilepath);
		if (retval) {
			goto tap_error_and_cleanup;
		}
		// Finished with the CSV filepath
		g_free(csvfilepath);

		// Now we can use the job ID to provide the source_ids for a Datalink query
		datalink_url = g_string_new("https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=");
		// Append the retrieval type string
		switch (type) {
			case EPOCH_PHOTOMETRY:
				g_string_append(datalink_url, "EPOCH_PHOTOMETRY");
				break;
			case XP_SAMPLED:
				g_string_append(datalink_url, "XP_SAMPLED");
				break;
			case XP_CONTINUOUS:
				g_string_append(datalink_url, "XP_CONTINUOUS");
				break;
			case MCMC_GSPPHOT:
				g_string_append(datalink_url, "MCMC_GSPPHOT");
				break;
			case MCMC_MSC:
				g_string_append(datalink_url, "MCMC_MSC");
				break;
			case RVS:
				g_string_append(datalink_url, "RVS");
				break;
			case ALL:
				g_string_append(datalink_url, "ALL");
				break;
			default:
				siril_debug_print("Error: siril_gaiadr3_datalink_query called with unsupported retrieval type, cannot proceed\n");
				goto datalink_download_error;
		}

		// Set the data structure and format. Siril will always consume the data as FITS,
		// as we already use libcfitsio and it's easier than adding support for VOTables
		g_string_append(datalink_url, "&DATA_STRUCTURE=RAW&FORMAT=FITS");

		// Append the source IDs (using the job number)
		g_string_append(datalink_url, "&ID=job:");
		g_string_append(datalink_url,job_id);
		g_string_append(datalink_url, ".source_id");

		siril_debug_print("Datalink url: %s\n", datalink_url->str);
		siril_log_message(_("Submitting spectral data request to ESA Gaia DR3 catalog. This may take several seconds to complete...\n"));
		gchar *datalink_buffer = fetch_url(datalink_url->str, &length, &fetch_url_error, FALSE);

		siril_debug_print("datalink_buffer length: %lu\n", length);
		g_string_free(datalink_url, TRUE);
		datalink_url = NULL;
		if (fetch_url_error)
			goto datalink_download_error;
		g_free(job_id);
		if (retrieval_product_is_in_cache) {
			siril_log_message(_("Using already downloaded datalink product\n"));
			return 0;
		}

		// the catalog needs to be downloaded, prepare the output stream
		g_unlink(filepath); // If no file exists this call may fail, that's fine
		file = g_file_new_for_path(filepath);
		output_stream = (GOutputStream*) g_file_create(file, G_FILE_CREATE_NONE, NULL, &error);
		if (!output_stream || !datalink_buffer || length == 0) {
			goto datalink_download_error;
		}
		if (!g_output_stream_write_all(output_stream, datalink_buffer, length, NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			remove_file = TRUE;
			goto datalink_download_error;
		}
		siril_log_message(_("Gaia DR3 datalink query succeeded: cached as %s\n"), filepath);
		*datalink_path = g_strdup(filepath);
		g_free(filepath);
		g_free(str);
		return 0;
	} else {
		siril_log_message(_("Using already downloaded catalogue %s\n"), catalog_to_str(siril_cat->cat_index));
		// Populate the siril_catalog with the data from the initial query
		int retval = siril_catalog_load_from_file(siril_cat, csvfilepath);
		*datalink_path = g_strdup(filepath);
		if (retval) {
			goto tap_error_and_cleanup;
		}
		g_free(csvfilepath);
		g_free(filepath);
		return 0;
	}

tap_error_and_cleanup:
	// Cleanup
    g_free(url);
    g_free(data);
	g_free(job_id);
	g_free(csvfilepath);
	g_free(filepath);
	g_string_free(querystring, TRUE);
	return -1;

datalink_download_error:
	if (datalink_url)
		g_string_free(datalink_url, TRUE);
	if (error) {
		siril_log_color_message(_("Cannot create catalogue file %s (%s)\n"), "red", filepath, error->message);
		g_clear_error(&error);
		}
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
	return -1;
}
