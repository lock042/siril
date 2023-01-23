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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h> // for waitpid(2)
#include <sys/wait.h> // for waitpid(2)

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_date.h"
#include "core/siril_world_cs.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/search_objects.h"
#include "algos/annotate.h"
#include "algos/photometry.h"
#include "algos/photometric_cc.h"
#include "algos/siril_wcs.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/catalogues.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/match.h"
#include "registration/matching/apply_match.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/project_coords.h"
#include "gui/message_dialog.h"

#include "astrometry_solver.h"

#define DOWNSAMPLE_FACTOR 0.25
#define CONV_TOLERANCE 1E-8

#undef DEBUG		/* get some of diagnostic output */

typedef struct {
	point size;
	SirilWorldCS *px_cat_center;	// the original target first, but can get refined
	SirilWorldCS *image_center;
	double crpix[2];
	double pixel_size;		// in µm
	double focal_length;		// in mm
	Homography H;			// for matching results printing
	gboolean image_is_flipped;
} solve_results;

struct sky_object platedObject[RESOLVER_NUMBER];

static void fov_in_DHMS(double var, gchar *fov) {
	int deg, decM;
	double decS;

	if (var < 0) {
		fprintf(stdout, "fov_in_DHMS: negative value, should not happened\n");
		return;
	}
	deg = (int) var;
	decM = abs((int) ((var - deg) * 60));
	decS = (fabs((var - deg) * 60) - decM) * 60;
	if (deg > 0)
		g_snprintf(fov, 256, "%02dd %02dm %.2lfs", deg, decM, decS);
	else if (decM > 0)
		g_snprintf(fov, 256, "%02d\' %.2lf\"", decM, decS);
	else if (decS > 0.0)
		g_snprintf(fov, 256, "%.2lf\"", decS);
}

void free_Platedobject() {
	for (int i = 0; i < RESOLVER_NUMBER; i++) {
		if (platedObject[i].name) {
			siril_world_cs_unref(platedObject[i].world_cs);
			free(platedObject[i].name);
			platedObject[i].name = NULL;
		}
	}
}

/* get resolution in arcsec per pixel */
double get_resolution(double focal, double pixel) {
	return RADCONV / focal * pixel;
}

/* get diagonal field of view in arcmin, resolution in arcsec/px */
double get_fov_arcmin(double resolution, int rx, int ry) {
	uint64_t sqr_radius = rx * rx + ry * ry;
	double radius = resolution * sqrt((double)sqr_radius);	// in arcsec
	return radius / 60.0;	// in arcminutes
}

/* get half field of view in arcmin, or angle from image centre, resolution in arcsec/px */
double get_radius_deg(double resolution, int rx, int ry) {
	uint64_t sqr_radius = (rx * rx + ry * ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in arcsec
	return radius / 3600.0;	// in degrees
}

double compute_mag_limit_from_fov(double fov_degrees) {
	// Empiric formula for 1000 stars at 20 deg of galactic latitude
	double autoLimitMagnitudeFactor = 14.5;
	double m = autoLimitMagnitudeFactor * pow(fov_degrees, -0.179);
	// for astrometry, it can be useful to go down to mag 20, for
	// photometry the catalog's limit is 17 for APASS and 18 for NOMAD
	return round(100.0 * min(20.0, max(7.0, m))) / 100;
}

static void compute_limit_mag(struct astrometry_data *args) {
	if (args->mag_mode == LIMIT_MAG_ABSOLUTE)
		args->limit_mag = args->magnitude_arg;
	else {
		args->limit_mag = compute_mag_limit_from_fov(args->used_fov / 60.0);
		if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
			args->limit_mag += args->magnitude_arg;
	}
	siril_debug_print("using limit magnitude %f\n", args->limit_mag);
}

gchar *get_catalog_url(SirilWorldCS *center, double mag_limit, double radius, int type) {
	GString *url;
	gchar *coordinates;
	gchar *mag;
	gchar *fov;

	coordinates = g_strdup_printf("%f+%f", siril_world_cs_get_alpha(center), siril_world_cs_get_delta(center));
	mag = g_strdup_printf("%2.2lf", mag_limit);
	fov = g_strdup_printf("%2.1lf", radius);

	url = g_string_new("http://vizier.u-strasbg.fr/viz-bin/asu-tsv?-source=");
	switch (type) {
	case CAT_NOMAD:
		url = g_string_append(url, "NOMAD&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Vmag%20Bmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	default:
	case CAT_TYCHO2:
		url = g_string_append(url, "I/259/tyc2&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAmdeg%20DEmdeg%20VTmag%20BTmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&VTmag=<");
		url = g_string_append(url, mag);
		break;
	case CAT_GAIADR3:
		url = g_string_append(url, "I/355/gaiadr3&-out.meta=-h-u-D&-out.add=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Gmag%20BPmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Gmag=<");
		url = g_string_append(url, mag);
		break;
	case CAT_PPMXL:
		url = g_string_append(url, "I/317&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Jmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Jmag=<");
		url = g_string_append(url, mag);
		break;
	case CAT_BRIGHT_STARS:
		url = g_string_append(url, "V/50/catalog&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out.add=_RAJ,_DEJ&-out=Vmag&-out=B-V");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	case CAT_APASS: // for photometry only
		url = g_string_append(url, "APASS&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Vmag%20Bmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	}

	g_free(coordinates);
	g_free(mag);
	g_free(fov);

	return g_string_free(url, FALSE);
}

#if defined HAVE_LIBCURL
/*****
 * HTTP functions
 ****/

static CURL *curl;
static const int DEFAULT_FETCH_RETRIES = 10;

struct ucontent {
	char *data;
	size_t len;
};

static void init() {
	if (!curl) {
		printf("initializing CURL\n");
		curl_global_init(CURL_GLOBAL_ALL);
		curl = curl_easy_init();
		if (g_getenv("CURL_CA_BUNDLE"))
			curl_easy_setopt(curl, CURLOPT_CAINFO, g_getenv("CURL_CA_BUNDLE"));
	}

	if (!curl)
		exit(EXIT_FAILURE);
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

static char *fetch_url(const char *url) {
	struct ucontent *content = malloc(sizeof(struct ucontent));
	char *result = NULL, *error = NULL;
	long code;
	int retries;
	unsigned int s;

	init();
	retries = DEFAULT_FETCH_RETRIES;

retrieve:
	content->data = malloc(1);
	content->data[0] = '\0';
	content->len = 0;

	curl_easy_setopt(curl, CURLOPT_URL, url);
	curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	curl_easy_setopt(curl, CURLOPT_USERAGENT, PACKAGE_STRING);

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
			siril_log_message(_("Failed to download page %s (error %ld)\n"), url, code);

			if (retries) {
				s = 2 * (DEFAULT_FETCH_RETRIES - retries) + 2;
				siril_debug_print("Wait %us before retry\n", s);
				sleep(s);

				free(content->data);
				retries--;
				goto retrieve;
			}

			break;
		default:
			if (content->data[0] == '#') {
				// special case of ephemcc where we need to parse the output
				gchar **token = g_strsplit(content->data, "\n", -1);
				int nlines = g_strv_length(token);
				int line;
				for (line = 0; line < nlines; line++) {
					if (token[line][0] != '\0' && token[line][0] != '#') {
						error = siril_log_message(_("Fetch failed with code %ld\n%s\n"), code, token[line]);
						break;
					}
				}
				g_strfreev(token);
			}
			if (!error)
				error = siril_log_message(_("Fetch failed with code %ld\n%s"), code, content->data);
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Online service error"), error);
		}
	}

	curl_easy_cleanup(curl);
	curl = NULL;
	if (!result || content->len == 0) {
		free(content->data);
		result = NULL;
	}
	free(content);

	return result;
}
#else
static gchar *fetch_url(const gchar *url) {
	GFile *file = g_file_new_for_uri(url);
	GError *error = NULL;
	gchar *content = NULL;

	siril_debug_print("fetch_url(): %s\n", url);

	if (!g_file_load_contents(file, NULL, &content, NULL, NULL, &error)) {
		siril_log_message(_("Error loading url: [%s] - %s\n"), url, error->message);
		g_clear_error(&error);
	}
	g_object_unref(file);
	return content;
}
#endif

gchar *search_in_online_conesearch(struct astrometry_data *args) {
	if (!args->fit->date_obs)
		return NULL;
	double ra, dec;
	center2wcs(args->fit, &ra, &dec);

	// https://vo.imcce.fr/webservices/skybot/?conesearch
	GString *string_url = g_string_new(SKYBOT);
	string_url = g_string_append(string_url, "&-ep=");
	gchar *formatted_date = date_time_to_FITS_date(args->fit->date_obs);
	string_url = g_string_append(string_url, formatted_date);
	string_url = g_string_append(string_url, "&-ra=");		// RA
	g_string_append_printf(string_url, "%lf", ra);
	string_url = g_string_append(string_url, "&-dec=");		// DEC
	g_string_append_printf(string_url, "%lf", dec);
	string_url = g_string_append(string_url, "&-rm=");		// FOV
	g_string_append_printf(string_url, "%lf", get_fov_arcmin(args->scale, args->fit->rx, args->fit->ry));
	string_url = g_string_append(string_url, "&-mime=text");
	string_url = g_string_append(string_url, "&-output=object");
	string_url = g_string_append(string_url, "&-loc=500");
	string_url = g_string_append(string_url, "&-filter=0");
	string_url = g_string_append(string_url, "&-objFilter=111");
	string_url = g_string_append(string_url, "&-refsys=EQJ2000");
	string_url = g_string_append(string_url, "&-from=Siril;");

	gchar *url = g_string_free(string_url, FALSE);
	gchar *cleaned_url = url_cleanup(url);
	gchar *result = fetch_url(cleaned_url);
	siril_debug_print(_("URL: %s \n"), cleaned_url);

	g_free(cleaned_url);
	g_free(url);
	g_free(formatted_date);

	return result;
}

gchar *search_in_online_catalogs(const gchar *object, query_server server) {
	GString *string_url;
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
		siril_log_message(_("Searching for solar system object %s on observation date %s\n"), name, formatted_date);

		if (!gfit.sitelat || !gfit.sitelong) {
			siril_log_color_message(_("No topocentric data available. Set to geocentric\n"), "salmon");
		} else {
			siril_log_message(_("at LAT: %f, LONG: %f\n"), gfit.sitelat, gfit.sitelong);
		}
		g_free(formatted_site);
		g_free(formatted_date);
		break;
	case QUERY_SERVER_SIMBAD_PHOTO:  // SIMBAD request to get the magnitudes (BVRIJ) for a particular star
		string_url = g_string_new(name);
		g_string_replace(string_url, "+", "%2B", 0);
		g_string_replace(string_url, "-", "%2D", 0);
		string_url = g_string_prepend(string_url, SIMBADPHOTO);
		string_url = g_string_append(string_url, "';");
		siril_log_message(_("Searching %s in SIMBAD(photo)...\n"), name);
		break;
	}

	gchar *url = g_string_free(string_url, FALSE);
	gchar *cleaned_url = url_cleanup(url);
	gchar *result = fetch_url(cleaned_url);
	siril_debug_print(_("URL: %s \n"), cleaned_url);

	g_free(cleaned_url);
	g_free(url);
	g_free(name);

	return result;
}

/* parse the result from search_in_catalogs(), for object name to coordinates conversion */
int parse_content_buffer(char *buffer, struct sky_object *obj) {
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

GFile *download_catalog(online_catalog onlineCatalog, SirilWorldCS *catalog_center, double radius_arcmin, double mag) {
	gchar *buffer = NULL;
	GError *error = NULL;
	/* ------------------- get Vizier catalog in catalog.dat -------------------------- */

	/* check if catalogue already exist in cache */
	gchar *str = g_strdup_printf("cat-%d-%lf-%lf-%lf-%lf.cat",
			(int) onlineCatalog,
			siril_world_cs_get_alpha(catalog_center),
			siril_world_cs_get_delta(catalog_center),
			radius_arcmin, mag);
	siril_debug_print("Catalogue file: %s\n", str);
	GFile *file = g_file_new_build_filename(g_get_tmp_dir(), str, NULL);
	g_free(str);

	GOutputStream *output_stream = (GOutputStream*) g_file_create(file, G_FILE_CREATE_NONE, NULL, &error);

	if (!output_stream) {
		if (error) {
			/* if file already exists */
			if (error->code == G_IO_ERROR_EXISTS) {
				GFileInfo *info = g_file_query_info(file, G_FILE_ATTRIBUTE_TIME_MODIFIED "," G_FILE_ATTRIBUTE_STANDARD_SIZE, 0, NULL, NULL);
				/* test if size is not 0 */
				if ((g_file_info_get_size(info)) == 0) {
					if (g_file_delete(file, NULL, &error)) {
						output_stream = (GOutputStream*) g_file_create(file, G_FILE_CREATE_NONE, NULL, &error);
						if (!output_stream) {
							goto download_error;
						}
					} else {
						goto download_error;
					}
				} else {
					siril_log_message(_("Using already downloaded star catalogue\n"));
				}
				g_clear_error(&error);
			} else {
				goto download_error;
			}
		} else {
			siril_log_color_message(_("Cannot create catalogue file %s for plate solving (%s)\n"), "red", g_file_peek_path(file), "unknown error");
			g_object_unref(file);
			return NULL;
		}
	}

	if (output_stream) {
		/* download and save */
		gchar *url = get_catalog_url(catalog_center, mag, radius_arcmin, onlineCatalog);
		buffer = fetch_url(url);
		g_free(url);

		if (buffer) {
			if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
				g_warning("%s\n", error->message);
				g_clear_error(&error);
				g_free(buffer);
				g_object_unref(output_stream);
				g_object_unref(file);
				return NULL;
			}
			g_object_unref(output_stream);
			g_free(buffer);
		}
	}
	return file;

download_error:
	g_warning("%s\n", error->message);
	siril_log_color_message(_("Cannot create catalogue file %s for plate solving (%s)\n"), "red", g_file_peek_path(file), error->message);
	g_clear_error(&error);
	g_object_unref(file);
	return NULL;
}

static gchar *project_catalog(GFile *catalogue_name, SirilWorldCS *catalog_center) {
	GError *error = NULL;
	gchar *foutput = NULL;
	/* --------- Project coords of Vizier catalog and save it into catalog.proj ------- */

	GFile *fproj = g_file_new_build_filename(g_get_tmp_dir(), "catalog.proj", NULL);

	/* We want to remove the file if already exisit */
	if (!g_file_delete(fproj, NULL, &error)
			&& !g_error_matches(error, G_IO_ERROR, G_IO_ERROR_NOT_FOUND)) {
		// deletion failed for some reason other than the file not existing:
		// so report the error
		g_warning("Failed to delete %s: %s", g_file_peek_path(fproj),
				error->message);
	}

	convert_catalog_coords(catalogue_name, catalog_center, fproj);
	foutput = g_file_get_path(fproj);
	g_object_unref(fproj);
	return foutput;
}

gboolean has_any_keywords() {
	return (gfit.focal_length > 0.0 ||
			gfit.pixel_size_x > 0.f ||
			gfit.pixel_size_y > 0.f ||
			(gfit.wcsdata.crval[0] > 0.0 && gfit.wcsdata.crval[1] != 0.0) ||
			(gfit.wcsdata.objctra[0] != '\0' && gfit.wcsdata.objctdec[0] != '\0') ||
			(gfit.wcsdata.ra != 0.0 && gfit.wcsdata.dec != 0.0));
}

SirilWorldCS *get_eqs_from_header(fits *fit) {
	if (fit->wcsdata.ra != 0.0 && fit->wcsdata.dec != 0.0)
		return siril_world_cs_new_from_a_d(fit->wcsdata.ra, fit->wcsdata.dec);

	else if (fit->wcsdata.objctra[0] != '\0' && fit->wcsdata.objctdec[0] != '\0')
		return siril_world_cs_new_from_objct_ra_dec(fit->wcsdata.objctra, fit->wcsdata.objctdec);

	else if (fit->wcsdata.crval[0] != 0.0 && fit->wcsdata.crval[1] != 0.0)
		return siril_world_cs_new_from_a_d(fit->wcsdata.crval[0], fit->wcsdata.crval[1]);
	return NULL;
}

/* Extract CDELT from CD matrix.*/
static void extract_cdelt_from_cd(double cd1_1, double cd1_2, double cd2_1,
		double cd2_2, double *cdelt1, double *cdelt2) {
	int sign;
	if ((cd1_1 * cd2_2 - cd1_2 * cd2_1) >= 0)
		sign = +1;
	else
		sign = -1;

	*cdelt1 = sqrt((cd1_1 * cd1_1) + (cd2_1 * cd2_1)) * sign;
	*cdelt2 = sqrt((cd1_2 * cd1_2) + (cd2_2 * cd2_2));
}

static void print_platesolving_results(solve_results *image, gboolean downsample) {
	double rotation, det, scaleX, scaleY, resolution;
	double inliers;
	char field_x[256] = { 0 };
	char field_y[256] = { 0 };

	float factor = (downsample) ? DOWNSAMPLE_FACTOR : 1.0;
	Homography H = image->H;

	/* Matching information */
	gchar *str = ngettext("%d pair match.\n", "%d pair matches.\n", H.pair_matched);
	str = g_strdup_printf(str, H.pair_matched);
	siril_log_message(str);
	g_free(str);
	inliers = 1.0 - (((double) H.pair_matched - (double) H.Inliers) / (double) H.pair_matched);
	siril_log_message(_("Inliers:%*.3f\n"), 14, inliers);

	/* Plate Solving */
	scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	resolution = (scaleX + scaleY) * 0.5 * factor; // we assume square pixels
	siril_log_message(_("Resolution:%*.3lf arcsec/px\n"), 11, resolution);

	/* rotation */
	rotation = atan2(H.h00 + H.h01, H.h10 + H.h11) * RADTODEG + 135.0;
	det = (H.h00 * H.h11 - H.h01 * H.h10); // determinant of rotation matrix (ad - bc)
	/* If the determinant of the top-left 2x2 rotation matrix is > 0
	 * the transformation is orientation-preserving. */

	if (det < 0)
		rotation = -90.0 - rotation;
	if (rotation < -180.0)
		rotation += 360.0;
	if (rotation > 180.0)
		rotation -= 360.0;
	siril_log_message(_("Rotation:%+*.2lf deg %s\n"), 12, rotation, det < 0.0 ? _("(flipped)") : "");

	siril_log_message(_("Focal length:%*.2lf mm\n"), 8, RADCONV * image->pixel_size / resolution);
	siril_log_message(_("Pixel size:%*.2lf µm\n"), 10, image->pixel_size);
	fov_in_DHMS(resolution * image->size.x / 3600.0, field_x);
	fov_in_DHMS(resolution * image->size.y / 3600.0, field_y);
	siril_log_message(_("Field of view:    %s x %s\n"), field_x, field_y);
}

static void print_image_center(solve_results *image) {
	gchar *alpha = siril_world_cs_alpha_format(image->image_center, " %02dh%02dm%02ds");
	gchar *delta = siril_world_cs_delta_format(image->image_center, "%c%02d°%02d\'%02d\"");
	siril_log_message(_("Image center: alpha: %s, delta: %s\n"), alpha, delta);
	g_free(alpha);
	g_free(delta);
}

static int read_NOMAD_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {

		if (line[0] == COMMENT_CHAR || is_blank(line) || g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0;
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);
		g_free(line);
		if (Vmag >= 30.0)
			continue;
		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = n < 5 || Bmag >= 30.0 ? -99.9 : Bmag - Vmag;
		star->phot = NULL;
		cstars[i++] = star;
		cstars[i] = NULL;
	}
	g_object_unref(data_input);
	sort_stars_by_mag(cstars, i);
	siril_log_message(_("Catalog NOMAD size: %d objects\n"), i);
	return i;
}

static int read_LOCAL_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {

		if (line[0] == COMMENT_CHAR || is_blank(line) || g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0;
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);
		g_free(line);
		if (Vmag >= 30.0)
			continue;
		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = n < 5 || Bmag >= 30.0 ? -99.9 : Bmag - Vmag;
		star->phot = NULL;
		cstars[i++] = star;
		cstars[i] = NULL;
	}
	g_object_unref(data_input);
	sort_stars_by_mag(cstars, i);
	siril_log_message(_("Local catalogs size: %d objects\n"), i);
	return i;
}

static int read_TYCHO2_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			continue;
		}
		if (is_blank(line)) {
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			continue;
		}
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = n < 5 ? -99.9 : Bmag - Vmag;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
	}
	g_object_unref(data_input);
	sort_stars_by_mag(cstars, i);
	siril_log_message(_("Catalog TYCHO-2 size: %d objects\n"), i);
	return i;
}

static int read_GAIA_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Gmag = 0.0, BPmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Gmag, &BPmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Gmag;
		star->BV = -99.9;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars_by_mag(cstars, i);
	siril_log_message(_("Catalog Gaia DR3 size: %d objects\n"), i);
	return i;
}

static int read_PPMXL_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Jmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		sscanf(line, "%lf %lf %lf %lf", &r, &x, &y, &Jmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Jmag;
		star->BV = -99.9;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars_by_mag(cstars, i);
	siril_log_message(_("Catalog PPMXL size: %d objects\n"), i);
	return i;
}

static int read_BRIGHT_STARS_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, BV = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &BV);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = BV;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars_by_mag(cstars, i);
	siril_log_message(_("Catalog Bright stars size: %d objects\n"), i);
	return i;
}

static int read_APASS_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = n < 5 ? -99.9 : Bmag - Vmag;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars_by_mag(cstars, i);
	siril_log_message(_("Catalog APASS size: %d objects\n"), i);
	return i;
}

static int read_catalog(GInputStream *stream, psf_star **cstars, int type) {
	switch (type) {
	case CAT_TYCHO2:
		return read_TYCHO2_catalog(stream, cstars);
	default:
	case CAT_LOCAL:
		return read_LOCAL_catalog(stream, cstars);
	case CAT_NOMAD:
		return read_NOMAD_catalog(stream, cstars);
	case CAT_GAIADR3:
		return read_GAIA_catalog(stream, cstars);
	case CAT_PPMXL:
		return read_PPMXL_catalog(stream, cstars);
	case CAT_BRIGHT_STARS:
		return read_BRIGHT_STARS_catalog(stream, cstars);
	case CAT_APASS:
		return read_APASS_catalog(stream, cstars);
	}
}

static TRANS H_to_linear_TRANS(Homography H) {
	TRANS trans = { 0 };

	trans.order = AT_TRANS_LINEAR;

	trans.a = H.h02;
	trans.b = H.h00;
	trans.c = H.h01;
	trans.d = H.h12;
	trans.e = H.h10;
	trans.f = H.h11;

	return trans;
}

static gboolean check_affine_TRANS_sanity(TRANS trans) {
	double var1 = fabs(trans.b) - fabs(trans.f);
	double var2 = fabs(trans.c) - fabs(trans.e);
	siril_debug_print("abs(b+f)=%f et abs(c+e)=%f\n", var1, var2);

	return ((fabs(var1) < 0.3) && fabs(var2) < 0.3);
}

static gboolean image_is_flipped(Homography H) {
	double det = (H.h00 * H.h11 - H.h01 * H.h10); // determinant of rotation matrix (ad - bc)
	return det < 0;
}

gboolean has_nonzero_coords() {
	for (int i = 0; i < RESOLVER_NUMBER; i++){
		if (fabs(platedObject[i].imageCenter.x) > 0.000001) return TRUE;
		if (fabs(platedObject[i].imageCenter.y) > 0.000001) return TRUE;
	}
	return FALSE;
}

// From projected starlist and center (ra,dec), go back to original ra and dec
// All formulas from AIPS memo #27 III.A.ii
// https://library.nrao.edu/public/memos/aips/memos/AIPSM_027.pdf

static void deproject_starlist(int num_stars, s_star *star_list, double ra0, double dec0, int doASEC) {
	ra0 *= DEGTORAD;
	dec0 *= DEGTORAD;
	s_star *currstar;
	currstar = star_list;
	for (int i = 0; i < num_stars; i++) {
		double xi = currstar->x;
		double eta = currstar->y;
		if (doASEC > 0) {
			xi /= RADtoASEC;
			eta /= RADtoASEC;
		}
		double delta_ra = atan(xi / (cos(dec0) - eta * sin(dec0)));
		double ra = ra0 + delta_ra;
		double dec = atan(cos(delta_ra) * (eta * cos(dec0) + sin(dec0)) / (cos(dec0) - eta * sin(dec0)));
		currstar->x = ra / DEGTORAD;
		currstar->y = dec / DEGTORAD;
		currstar = currstar->next;
	}
}

// From starlist in (ra,dec) and center (ra,dec), project in "pixels" (in arcsec)
// All formulas from AIPS memo #27 III.A.i
// https://library.nrao.edu/public/memos/aips/memos/AIPSM_027.pdf

static void project_starlist(int num_stars, s_star *star_list, double ra0, double dec0, int doASEC) {
	double delta_ra;
	dec0 *= DEGTORAD;
	s_star *currstar;
	currstar = star_list;
	for (int i = 0; i < num_stars; i++) {
		double ra = currstar->x;
		double dec = currstar->y;
		if ((ra < 10) && (ra0 > 350)) {
			delta_ra = (ra + 360) - ra0;
		} else if ((ra > 350) && (ra0 < 10)) {
			delta_ra = (ra - 360) - ra0;
		} else {
			delta_ra = ra - ra0;
		}
		delta_ra *= DEGTORAD;
		dec *= DEGTORAD;

		/*
		 * let's transform from (delta_RA, delta_Dec) to (xi, eta),
		 */
		double xx = cos(dec) * sin(delta_ra);
		double yy = sin(dec0) * sin(dec) + cos(dec0) * cos(dec) * cos(delta_ra);
		double xi = (xx / yy);
		xx = cos(dec0) * sin(dec) - sin(dec0) * cos(dec) * cos(delta_ra);
		double eta = (xx / yy);

		if (doASEC > 0) {
			xi *= RADtoASEC;
			eta *= RADtoASEC;
		}
		currstar->x = xi;
		currstar->y = eta;
		currstar = currstar->next;
	}
}

void print_updated_wcs_data(fits *fit) {
	/* debug output */
	siril_debug_print("****Updated WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, fit->wcsdata.crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, fit->wcsdata.crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, fit->wcsdata.crval[0]);
	siril_debug_print("crval2 = %*.12e\n", 20, fit->wcsdata.crval[1]);
	siril_debug_print("cdelt1 = %*.12e\n", 20, fit->wcsdata.cdelt[0]);
	siril_debug_print("cdelt2 = %*.12e\n", 20, fit->wcsdata.cdelt[1]);
	siril_debug_print("pc1_1  = %*.12e\n", 20, fit->wcsdata.pc[0][0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, fit->wcsdata.pc[0][1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, fit->wcsdata.pc[1][0]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, fit->wcsdata.pc[1][1]);
	siril_debug_print("******************************************\n");
}

/******
 *
 * Public functions
 */

void flip_bottom_up_astrometry_data(fits *fit) {
	/* flip pc matrix */
	fit->wcsdata.pc[0][1] = -fit->wcsdata.pc[0][1];
	fit->wcsdata.pc[1][1] = -fit->wcsdata.pc[1][1];

	/* update crpix */
	fit->wcsdata.crpix[1] = fit->ry - fit->wcsdata.crpix[1];

	print_updated_wcs_data(fit);
}

void reframe_astrometry_data(fits *fit, Homography H) {
	double pc1_1, pc1_2, pc2_1, pc2_2;
	point refpointout;

	pc1_1 = H.h00 * fit->wcsdata.pc[0][0] + H.h01 * fit->wcsdata.pc[0][1];
	pc1_2 = H.h10 * fit->wcsdata.pc[0][0] + H.h11 * fit->wcsdata.pc[0][1];
	pc2_1 = H.h00 * fit->wcsdata.pc[1][0] + H.h01 * fit->wcsdata.pc[1][1];
	pc2_2 = H.h10 * fit->wcsdata.pc[1][0] + H.h11 * fit->wcsdata.pc[1][1];

	point refpointin = {fit->wcsdata.crpix[0], fit->wcsdata.crpix[1]};
	cvTransformImageRefPoint(H, refpointin, &refpointout);

	fit->wcsdata.pc[0][0] = pc1_1;
	fit->wcsdata.pc[0][1] = pc1_2;
	fit->wcsdata.pc[1][0] = pc2_1;
	fit->wcsdata.pc[1][1] = pc2_2;
	fit->wcsdata.crpix[0] = refpointout.x;
	fit->wcsdata.crpix[1] = refpointout.y;

	print_updated_wcs_data(fit);
}

void wcs_cd_to_pc(double cd[][2], double pc[][2], double cdelt[2]) {
	extract_cdelt_from_cd(cd[0][0], cd[0][1], cd[1][0], cd[1][1], &cdelt[0], &cdelt[1]);

	pc[0][0] = cd[0][0] / cdelt[0];
	pc[0][1] = cd[0][1] / cdelt[0];
	pc[1][0] = cd[1][0] / cdelt[1];
	pc[1][1] = cd[1][1] / cdelt[1];
}

void wcs_pc_to_cd(double pc[][2], const double cdelt[2], double cd[][2]) {
	cd[0][0] = pc[0][0] * cdelt[0];
	cd[0][1] = pc[0][1] * cdelt[0];
	cd[1][0] = pc[1][0] * cdelt[1];
	cd[1][1] = pc[1][1] * cdelt[1];
}

static int match_catalog(psf_star **stars, int n_fit, psf_star **cstars, int n_cat, struct astrometry_data *args, solve_results *solution);
static int local_asnet_platesolve(psf_star **stars, int n_fit, struct astrometry_data *args, solve_results *solution);

#define CHECK_FOR_CANCELLATION if (!get_thread_run()) { args->message = g_strdup(_("Cancelled")); args->ret = 1; goto clearup; }

/* entry point for plate solving */
gpointer plate_solver(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	psf_star **stars = NULL;	// image stars
	psf_star **cstars = NULL;	// catalogue stars
	int n_fit = 0, n_cat = 0;	// number of image and catalogue stars

	args->ret = 1;
	args->message = NULL;

	if (args->verbose) {
		if (args->onlineCatalog == CAT_ASNET) {
			siril_log_message(_("Plate solving image from local catalogues for a field of view of %.2f degrees\n"), args->used_fov / 60.0);
		} else if (args->use_local_cat) {
			siril_log_message(_("Plate solving image from local catalogues for a field of view of %.2f"
						" degrees%s, using a limit magnitude of %.2f\n"),
					args->used_fov / 60.0,
					args->uncentered ? _(" (uncentered)") : "", args->limit_mag);
		} else {
			siril_log_message(_("Plate solving image from an online catalogue for a field of view of %.2f"
						" degrees%s, using a limit magnitude of %.2f\n"),
					args->used_fov / 60.0,
					args->uncentered ? _(" (uncentered)") : "", args->limit_mag);
		}
	}

	/* 1. Get catalogue stars for the field of view */
	if (args->onlineCatalog != CAT_ASNET) {
		/* obtaining a star catalogue */
		if (!args->catalog_file && !args->use_local_cat) {
			args->catalog_file = download_catalog(args->onlineCatalog, args->cat_center,
					args->used_fov * 0.5, args->limit_mag);
			if (!args->catalog_file) {
				args->message = g_strdup(_("Could not download the online star catalogue."));
				goto clearup;
			}
		}
		CHECK_FOR_CANCELLATION;

		cstars = new_fitted_stars(MAX_STARS);
		if (!cstars) {
			PRINT_ALLOC_ERR;
			goto clearup;
		}

		/* project and open the file */
		if (args->use_local_cat) {
			args->catalogStars = get_and_project_local_catalog(args->cat_center,
					args->used_fov / 120.0, args->limit_mag, FALSE);
		} else {
			args->catalogStars = project_catalog(args->catalog_file, args->cat_center);
			if (!args->catalogStars) {
				args->message = g_strdup(_("Cannot project the star catalog."));
				goto clearup;
			}
		}
		CHECK_FOR_CANCELLATION;
		GFile *catalog = g_file_new_for_path(args->catalogStars);
		GError *error = NULL;
		GInputStream *input_stream = (GInputStream*) g_file_read(catalog, NULL, &error);
		if (!input_stream) {
			if (error != NULL) {
				args->message = g_strdup_printf(_("Could not load the star catalog (%s)."), error->message);
				g_clear_error(&error);
			}
			args->message = g_strdup_printf(_("Could not load the star catalog (%s)."), "generic error");
			goto clearup;
		}

		n_cat = read_catalog(input_stream, cstars, args->onlineCatalog);
		g_object_unref(input_stream);
		g_object_unref(catalog);
		CHECK_FOR_CANCELLATION;
	}

	/* 2. Get image stars */
	if (!args->manual) {
		fits fit_backup = { 0 };	// original image in case of downscale
		gchar *header_backup = NULL;
		if (args->downsample) {
			siril_log_message(_("Down-sampling image for faster star detection by a factor %.2f\n"), DOWNSAMPLE_FACTOR);
			copyfits(args->fit, &fit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
			copy_fits_metadata(args->fit, &fit_backup);
			header_backup = g_strdup(args->fit->header);
			cvResizeGaussian(args->fit, DOWNSAMPLE_FACTOR * args->fit->rx, DOWNSAMPLE_FACTOR * args->fit->ry, OPENCV_AREA, FALSE);
			//args->fit->pixel_size_x /= DOWNSAMPLE_FACTOR;
			//args->fit->pixel_size_y /= DOWNSAMPLE_FACTOR;
			args->scale /= DOWNSAMPLE_FACTOR;
		}

		image im = { .fit = args->fit, .from_seq = NULL, .index_in_seq = -1 };
		// capping the detection to max usable number of stars
		int max_stars = args->for_photometry_cc ? n_cat : min(n_cat, BRIGHTEST_STARS);

		// TODO: use good layer
		stars = peaker(&im, 0, &com.pref.starfinder_conf, &n_fit, &(args->solvearea), FALSE,
				args->onlineCatalog != CAT_ASNET, max_stars,
				com.pref.starfinder_conf.profile, com.max_thread);

		if (args->downsample) {
			clearfits(args->fit);
			memcpy(args->fit, &fit_backup, sizeof(fits));
			args->fit->header = header_backup;
			memset(&fit_backup, 0, sizeof(fits));
		}
	} else {
		stars = com.stars;
		if (com.stars)
			while (com.stars[n_fit])
				n_fit++;
	}
	CHECK_FOR_CANCELLATION;

	if (!stars || n_fit < AT_MATCH_STARTN_LINEAR) {
		args->message = g_strdup_printf(_("There are not enough stars picked in the image. "
				"At least %d are needed."), AT_MATCH_STARTN_LINEAR);
		goto clearup;
	}
	if (args->verbose)
		siril_log_message(_("Using %d detected stars from image.\n"), n_fit);

	/* 3. Plate solving */
	solve_results solution = { 0 };
	solution.size.x = args->fit->rx;
	solution.size.y = args->fit->ry;
	solution.pixel_size = args->pixel_size;

	if (args->onlineCatalog == CAT_ASNET)
		args->ret = local_asnet_platesolve(stars, n_fit, args, &solution);
	else match_catalog(stars, n_fit, cstars, n_cat, args, &solution);
	if (args->ret)
		goto clearup;

	/* 4. Print results */
	if (args->verbose && args->onlineCatalog != CAT_ASNET)
		print_platesolving_results(&solution, args->downsample);
	print_image_center(&solution);

	/* 5. Run photometric color correction, if enabled */
	if (args->for_photometry_cc) {
		pcc_star *pcc_stars = NULL;
		int nb_pcc_stars;
#ifndef HAVE_WCSLIB
		siril_log_color_message(_("This operation (PCC) relies on the missing WCSLIB software, cannot continue.\n"), "red");
		args->ret = 1;
		goto clearup;
#endif
		if (args->use_local_cat) {
			double tra = siril_world_cs_get_alpha(solution.image_center);
			double tdec = siril_world_cs_get_delta(solution.image_center);
			double res = get_resolution(solution.focal_length, args->pixel_size);
			double radius = get_radius_deg(res, args->fit->rx, args->fit->ry);
			// for photometry, we can use fainter stars, 1.5 seems ok above instead of 2.0
			if (args->verbose)
				siril_log_message(_("Getting stars from local catalogues for PCC, limit magnitude %.2f\n"), args->limit_mag);
			if (get_stars_from_local_catalogues(tra, tdec, radius, args->fit, args->limit_mag, &pcc_stars, &nb_pcc_stars)) {
				siril_log_color_message(_("Failed to get data from the local catalogue, is it installed?\n"), "red");
				args->ret = 1;
			}
		} else {
			args->ret = project_catalog_with_WCS(args->catalog_file, args->fit,
					&pcc_stars, &nb_pcc_stars);
		}
		if (args->ret) {
			args->message = g_strdup(_("Using plate solving to identify catalogue stars in the image failed, is plate solving wrong?\n"));
			goto clearup;
		}
		args->pcc->stars = pcc_stars;
		args->pcc->nb_stars = nb_pcc_stars;
		args->pcc->fwhm = filtered_FWHM_average(stars, n_fit);
		if (args->downsample)
			args->pcc->fwhm /= DOWNSAMPLE_FACTOR;

		args->ret = photometric_cc(args->pcc);

		args->pcc = NULL; // freed in PCC code
		free(pcc_stars);
		pcc_stars = NULL;
		if (args->ret) {
			args->message = g_strdup_printf(_("An astrometric solution was found but photometry analysis of the %d stars failed. This generally happens if they are saturated in the image or if they are too faint to have B-V index information (mag > 18)\n"), nb_pcc_stars);
			//goto clearup; // still flip
		}
	}

	/* 6. Flip image if needed */
	if (args->flip_image && solution.image_is_flipped) {
		if (args->verbose)
			siril_log_color_message(_("Flipping image and updating astrometry data.\n"), "salmon");
		fits_flip_top_to_bottom(args->fit);
		flip_bottom_up_astrometry_data(args->fit);
		load_WCS_from_memory(args->fit);
		args->image_flipped = TRUE;
	}

	/* 7. Clean-up */
	if (solution.px_cat_center)
		siril_world_cs_unref(solution.px_cat_center);
	args->new_center = solution.image_center;

clearup:
	if (stars && !args->manual) {
		for (int i = 0; i < n_fit; i++)
			free_psf(stars[i]);
		free(stars);
	}
	if (!args->for_sequence)
		siril_world_cs_unref(args->cat_center);
	if (cstars)
		free_fitted_stars(cstars);
	if (args->catalog_file)
		g_object_unref(args->catalog_file);
	g_free(args->catalogStars);

	if (!args->for_sequence)
		siril_add_idle(end_plate_solver, args);
	int retval = args->ret;
	if (com.script) {
		if (args->ret)
			siril_log_message(_("Plate solving failed: %s\n"), args->message);
		if (!args->for_sequence)
			free(args);
	}
	return GINT_TO_POINTER(retval);
}

/* entry point for siril's plate solver based on catalogue matching */
static int match_catalog(psf_star **stars, int n_fit, psf_star **cstars, int n_cat, struct astrometry_data *args, solve_results *solution) {
	Homography H = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int max_trials = 0;
	s_star *star_list_A = NULL, *star_list_B = NULL;

	if (args->uncentered)
		max_trials = 20; //retry to converge if solve is done at an offset from the center

	/* make sure that arrays are not too small
	 * make sure that the max of stars is BRIGHTEST_STARS */
	int n = min(min(n_fit, n_cat), BRIGHTEST_STARS);

	double scale_min = args->scale - 0.2;
	double scale_max = args->scale + 0.2;
	int attempt = 1;
	while (args->ret && attempt <= 3) {
		free_stars(&star_list_A);
		free_stars(&star_list_B);
		args->ret = new_star_match(stars, cstars, n, nobj,
				scale_min, scale_max, &H,
				FALSE, NULL, NULL,
				AFFINE_TRANSFORMATION, &star_list_A, &star_list_B);
		if (attempt == 1) {
			scale_min = -1.0;
			scale_max = -1.0;
			nobj += 10;
		} else {
			nobj += 30;
		}
		attempt++;
		CHECK_FOR_CANCELLATION;
	}
	if (args->ret)	// give generic error message
		goto clearup;

	double conv = DBL_MAX;
	solution->px_cat_center = siril_world_cs_ref(args->cat_center);
	/* we only want to compare with linear function
	 * Maybe one day we will apply match with homography matrix
	 */
	TRANS trans = H_to_linear_TRANS(H);
	if (!check_affine_TRANS_sanity(trans)) {
		args->message = g_strdup(_("Transformation matrix is invalid, solve failed"));
		args->ret = 1;
		goto clearup;
	}

	double ra0, dec0;
	solution->crpix[0] = ((args->fit->rx - 1) / 2.0);
	solution->crpix[1] = ((args->fit->ry - 1) / 2.0);

	apply_match(solution->px_cat_center, solution->crpix, trans, &ra0, &dec0);
	int num_matched = H.pair_matched;
	int trial = 0;

	/* try to get a better solution in case of uncentered selection */
	while (conv > CONV_TOLERANCE && trial < max_trials){
		double rainit = siril_world_cs_get_alpha(args->cat_center);
		double decinit = siril_world_cs_get_delta(args->cat_center);
		double orig_ra0 = ra0;
		double orig_dec0 = dec0;

		deproject_starlist(num_matched, star_list_B, rainit, decinit, 1);
		siril_debug_print("Deprojecting from: alpha: %s, delta: %s\n",
				siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
				siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));
		args->cat_center = siril_world_cs_new_from_a_d(ra0, dec0);
		siril_world_cs_unref(solution->px_cat_center);
		solution->px_cat_center = siril_world_cs_new_from_a_d(ra0, dec0);

		project_starlist(num_matched, star_list_B, ra0, dec0, 1);
		siril_debug_print("Reprojecting to: alpha: %s, delta: %s\n",
				siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
				siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));

		double scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
		double scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
		double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels

		double focal = RADCONV * solution->pixel_size / resolution;
		siril_debug_print("Current focal: %0.2fmm\n", focal);

		if (atPrepareHomography(num_matched, star_list_A, num_matched, star_list_B, &H, FALSE, NULL, NULL, AFFINE_TRANSFORMATION)){
			args->message = g_strdup(_("Updating homography failed."));
			args->ret = 1;
			break;
		}
		trans = H_to_linear_TRANS(H);
		apply_match(solution->px_cat_center, solution->crpix, trans, &ra0, &dec0);

		conv = fabs((dec0 - orig_dec0) / orig_dec0) + fabs((ra0 - orig_ra0) / orig_ra0);

		trial++;
		CHECK_FOR_CANCELLATION;
	}
	if (args->ret)	// after the break
		goto clearup;

	memcpy(&solution->H, &H, sizeof(Homography));
	double scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	double scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels
	solution->focal_length = RADCONV * solution->pixel_size / resolution;
	solution->image_center = siril_world_cs_new_from_a_d(ra0, dec0);
	if (max_trials == 0) {
		siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
	} else if (trial == max_trials) {
		siril_debug_print("No convergence found: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
	} else {
		siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f at iteration #%d\n", ra0, dec0, trial);
	}

	solution->image_is_flipped = image_is_flipped(H);

	/* compute cd matrix */
	double ra7, dec7, delta_ra;

	/* first, convert center coordinates from deg to rad: */
	dec0 *= DEGTORAD;
	ra0 *= DEGTORAD;

	/* make 1 step in direction crpix1 */
	double crpix1[] = { solution->crpix[0] + DOWNSAMPLE_FACTOR, solution->crpix[1] };
	apply_match(solution->px_cat_center, crpix1, trans, &ra7, &dec7);

	dec7 *= DEGTORAD;
	ra7 *= DEGTORAD;

	delta_ra = ra7 - ra0;
	if (delta_ra > +M_PI)
		delta_ra = 2.0 * M_PI - delta_ra;
	if (delta_ra < -M_PI)
		delta_ra = delta_ra - 2.0 * M_PI;
	double cd1_1 = (delta_ra) * cos(dec0) * RADTODEG;
	double cd2_1 = (dec7 - dec0) * RADTODEG;

	/* make 1 step in direction crpix2
	 * WARNING: we use -1 because of the Y axis reversing */
	double crpix2[] = { solution->crpix[0], solution->crpix[1] - DOWNSAMPLE_FACTOR };
	apply_match(solution->px_cat_center, crpix2, trans, &ra7, &dec7);

	dec7 *= DEGTORAD;
	ra7 *= DEGTORAD;

	delta_ra = ra7 - ra0;
	if (delta_ra > +M_PI)
		delta_ra = 2.0 * M_PI - delta_ra;
	if (delta_ra < -M_PI)
		delta_ra = delta_ra - 2.0 * M_PI;
	double cd1_2 = (delta_ra) * cos(dec0) * RADTODEG;
	double cd2_2 = (dec7 - dec0) * RADTODEG;

	if (args->downsample)
		solution->focal_length /= DOWNSAMPLE_FACTOR;

	CHECK_FOR_CANCELLATION;

	// saving state for undo before modifying fit structure
	if (!com.script) {
		const char *undo_str = args->for_photometry_cc ? _("Photometric CC") : _("Plate Solve");
		undo_save_state(args->fit, undo_str);
	}

	/**** Fill wcsdata fit structure ***/
	args->fit->wcsdata.equinox = 2000.0;
	args->fit->focal_length = solution->focal_length;
	args->fit->pixel_size_x = args->fit->pixel_size_y = solution->pixel_size;
	solution->crpix[0] /= DOWNSAMPLE_FACTOR;
	solution->crpix[1] /= DOWNSAMPLE_FACTOR;

	args->fit->wcsdata.ra = siril_world_cs_get_alpha(solution->image_center);
	args->fit->wcsdata.dec = siril_world_cs_get_delta(solution->image_center);

	args->fit->wcsdata.crpix[0] = solution->crpix[0];
	args->fit->wcsdata.crpix[1] = solution->crpix[1];
	args->fit->wcsdata.crval[0] = args->fit->wcsdata.ra;
	args->fit->wcsdata.crval[1] = args->fit->wcsdata.dec;

	args->fit->wcsdata.pltsolvd = TRUE;
	g_snprintf(args->fit->wcsdata.pltsolvd_comment, FLEN_COMMENT, "Siril internal solver");

	gchar *ra = siril_world_cs_alpha_format(solution->image_center, "%02d %02d %.3lf");
	gchar *dec = siril_world_cs_delta_format(solution->image_center, "%c%02d %02d %.3lf");

	g_sprintf(args->fit->wcsdata.objctra, "%s", ra);
	g_sprintf(args->fit->wcsdata.objctdec, "%s", dec);

	g_free(ra);
	g_free(dec);

	CHECK_FOR_CANCELLATION;
	double cdelt1, cdelt2;

	extract_cdelt_from_cd(cd1_1, cd1_2, cd2_1, cd2_2, &cdelt1, &cdelt2);

	args->fit->wcsdata.cdelt[0] = cdelt1;
	args->fit->wcsdata.cdelt[1] = cdelt2;

	/* PC + CDELT seems to be the preferred approach
	 * according to Calabretta private discussion
	 *
	 *    |cd11 cd12|  = |cdelt1      0| * |pc11 pc12|
	 *    |cd21 cd22|    |0      cdelt2|   |pc21 pc22|
	 */

	args->fit->wcsdata.pc[0][0] = cd1_1 / cdelt1;
	args->fit->wcsdata.pc[0][1] = cd1_2 / cdelt1;
	args->fit->wcsdata.pc[1][0] = cd2_1 / cdelt2;
	args->fit->wcsdata.pc[1][1] = cd2_2 / cdelt2;

	siril_debug_print("****Solution found: WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, solution->crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, solution->crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, args->fit->wcsdata.ra);
	siril_debug_print("crval2 = %*.12e\n", 20, args->fit->wcsdata.dec);
	siril_debug_print("cdelt1 = %*.12e\n", 20, cdelt1);
	siril_debug_print("cdelt2 = %*.12e\n", 20, cdelt2);
	siril_debug_print("pc1_1  = %*.12e\n", 20, args->fit->wcsdata.pc[0][0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, args->fit->wcsdata.pc[0][1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, args->fit->wcsdata.pc[1][0]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, args->fit->wcsdata.pc[1][1]);
	siril_debug_print("******************************************\n");

	load_WCS_from_memory(args->fit);
clearup:
	free_stars(&star_list_A);
	free_stars(&star_list_B);
	return args->ret;
}

static int local_asnet_platesolve(psf_star **stars, int n_fit, struct astrometry_data *args, solve_results *solution) {
	// command to run:
	// solve-field -p -N none -R none -M none -B none -U none --temp-axy -S none --crpix-center -T -H 0.2 -u app -X XIMAGE -Y YIMAGE sources.xyls
	const char *table_filename = "temp_sources.xyls";	// not reentrant because of this
	if (save_list_as_FITS_table(table_filename, stars, n_fit, args->fit->rx, args->fit->ry)) {
		siril_log_message(_("Failed to create the input data for solve-field\n"));
		return 1;
	}

	char low_scale[16], high_scale[16];
	sprintf(low_scale, "%f", args->scale * 0.8);
	sprintf(high_scale, "%f", args->scale * 1.2);

	char *sfargs[50] = {
#ifdef _WIN32
		"solve-field.exe",
#else
		"solve-field",
#endif
		"-p", "-O", "-N", "none", "-R", "none", "-M", "none", "-B", "none",
		"-U", "none", "--temp-axy", "-S", "none", "--crpix-center", "-T",
		"-u", "arcsecperpix", "-L", low_scale, "-H", high_scale, NULL };

	if (args->cat_center) {
		char start_ra[16], start_dec[16];
		sprintf(start_ra, "%f", siril_world_cs_get_alpha(args->cat_center));
		sprintf(start_dec, "%f", siril_world_cs_get_delta(args->cat_center));
		char *additional_args[] = { "--ra", start_ra, "--dec", start_dec, "--radius", "10", (char*)table_filename, NULL };
		append_elements_to_array(sfargs, additional_args);
	} else {
		char *additional_args[] = { (char*)table_filename, NULL };
		append_elements_to_array(sfargs, additional_args);
	}
	gchar *command = build_string_from_words(sfargs);
	siril_debug_print("Calling solve-field:\n%s\n", command);
	g_free(command);

	/* call solve-field */
	int pipefd[2];
	if (pipe(pipefd)) {
		perror("failed to create pipe for solve-field");
		return -1;
	}
	int pid = fork();
	if (pid == -1) {
		perror("failed to fork for solve-field");
		return -1;
	}
	if (pid == 0) {
		close(pipefd[0]);
		dup2(pipefd[1], 1);
		dup2(pipefd[1], 2);
		close(pipefd[1]);
		if (execvp(sfargs[0], sfargs)) {
			perror("could not execute solve-field");
			return -1;
		}
	}
	else {
		char buffer[1024];
		close(pipefd[1]);
		FILE *file = fdopen(pipefd[0], "r");
		if (!file) {
			perror("fdopen");
			waitpid(pid, NULL, 0);
			close(pipefd[0]);
			return -1;
		}
		gboolean success = FALSE;
		while (fgets(buffer, 1024, file)) {
			siril_debug_print("solver: %s", buffer);
			if (!strncmp(buffer, "Did not solve", strlen("Did not solve"))) {
				siril_log_color_message(_("No astrometric solution found\n"), "red");
				break;
			}
			if (!strncmp(buffer, "Field center: (RA,Dec)", 22)) {
				siril_debug_print("Found a solution, waiting for EOF and exit\n");
				success = TRUE;
			}
		}
		waitpid(pid, NULL, 0);
		siril_debug_print("solver exited\n");
		fclose(file);
		close(pipefd[0]);
		if (!success)
			return 1;

		/* get the results from the .wcs file */
		const char *wcs_filename = "temp_sources.wcs";
		fits result = { 0 };
		if (read_fits_metadata_from_path(wcs_filename, &result)) {
			siril_log_color_message(_("Could not read the solution from solve-field (expected in file %s)\n"), "red", wcs_filename);
			return 1;
		}

		memcpy(&args->fit->wcsdata, &result.wcsdata, sizeof(wcs_info));
		memset(&result.wcsdata, 0, sizeof(wcs_info));
#ifdef HAVE_WCSLIB
		args->fit->wcslib = result.wcslib;
		result.wcslib = NULL;
#endif
		clearfits(&result);
		g_unlink(table_filename);
		g_unlink(wcs_filename);

		args->fit->wcsdata.pltsolvd = TRUE;
		strcpy(args->fit->wcsdata.pltsolvd_comment, "This is a WCS header was created by Astrometry.net.");
		// asnet puts more info in the HISTORY and the console log in COMMENT fields
		/* populate the solution structure */
		//solution->px_cat_center = siril_world_cs_ref(args->cat_center);
		//solution->crpix[0] = args->fit->wcsdata.crpix[0];
		//solution->crpix[1] = args->fit->wcsdata.crpix[1];
		//solution->H ?
		solution->image_center = siril_world_cs_new_from_a_d(
				args->fit->wcsdata.crval[0],
				args->fit->wcsdata.crval[1]);
		// TODO: handle args->downsample
		return 0;
	}
	return 0;
}

void process_plate_solver_input(struct astrometry_data *args) {
	args->scale = get_resolution(args->focal_length, args->pixel_size);

	rectangle croparea = { 0 };
	if (!args->manual) {
		// first checking if there is a selection or if the full field is to be used
		if (com.selection.w != 0 && com.selection.h != 0) {
			memcpy(&croparea, &com.selection, sizeof(rectangle));
			siril_log_color_message(_("Warning: using the current selection to detect stars\n"), "salmon");
		} else {
			croparea.x = 0;
			croparea.y = 0;
			croparea.w = args->fit->rx;
			croparea.h = args->fit->ry;
		}
		double fov_arcmin = get_fov_arcmin(args->scale, croparea.w, croparea.h);
		siril_debug_print("image fov for given sampling: %f arcmin\n", fov_arcmin);

		// then apply or not autocropping to 5deg (300 arcmin)
		args->used_fov = args->autocrop ? min(fov_arcmin, 300.) : fov_arcmin;
		double cropfactor = (args->used_fov < fov_arcmin) ? args->used_fov / fov_arcmin : 1.0;
		if (cropfactor != 1.0) {
			croparea.x += (int) ((croparea.w - croparea.w * cropfactor) / 2);
			croparea.y += (int) ((croparea.h - croparea.h * cropfactor) / 2);
			croparea.w = (int) (cropfactor * croparea.w);
			croparea.h = (int) (cropfactor * croparea.h);
			siril_debug_print("Auto-crop factor: %.2f\n", cropfactor);
		}

		if (com.selection.w != 0 && com.selection.h != 0) {
			// detect if the selection is not centered enough that it matters
			double thr = max(args->fit->rx, args->fit->ry) / 10.0;
			args->uncentered =
				fabs(croparea.x + 0.5 * croparea.w - 0.5 * args->fit->rx) > thr ||
				fabs(croparea.y + 0.5 * croparea.h - 0.5 * args->fit->ry) > thr;
			if (args->uncentered)
				siril_debug_print("detected uncentered selection\n");
			else siril_debug_print("selection considered centered\n");
		} else {
			args->uncentered = FALSE;
		}

		if (args->downsample) {
			croparea.w *= DOWNSAMPLE_FACTOR;
			croparea.h *= DOWNSAMPLE_FACTOR;
			croparea.x *= DOWNSAMPLE_FACTOR;
			croparea.y *= DOWNSAMPLE_FACTOR;
		}
	} else { //stars manual selection - use full field centered
		args->used_fov = get_fov_arcmin(args->scale, args->fit->rx, args->fit->ry);
		args->uncentered = FALSE;
		if (com.selection.w != 0 && com.selection.h != 0)
			siril_log_message(_("Selection is not used in manual star selection mode\n"));
		// TODO: we could actually check if stars are in the selection
	}

	if (croparea.w == args->fit->rx && croparea.h == args->fit->ry)
		memset(&croparea, 0, sizeof(rectangle));
	else siril_debug_print("reduced area for the solve: %d, %d, %d x %d%s\n", croparea.x, croparea.y, croparea.w, croparea.h, args->downsample ? " (down-sampled)" : "");
	memcpy(&(args->solvearea), &croparea, sizeof(rectangle));

	compute_limit_mag(args); // to call after having set args->used_fov
}

static int astrometry_prepare_hook(struct generic_seq_args *arg) {
	struct astrometry_data *args = (struct astrometry_data *)arg->user;
	fits fit = { 0 };
	// load ref metadata in fit
	if (seq_read_frame_metadata(arg->seq, sequence_find_refimage(arg->seq), &fit))
		return 1;
	if (!args->cat_center)
		args->cat_center = get_eqs_from_header(&fit);
	if (args->pixel_size <= 0.0) {
		args->pixel_size = max(fit.pixel_size_x, fit.pixel_size_y);
		if (args->pixel_size <= 0.0) {
			args->pixel_size = com.pref.starfinder_conf.pixel_size_x;
			if (args->pixel_size <= 0.0) {
				siril_log_color_message(_("Pixel size not found in image or in settings, cannot proceed\n"), "red");
				return 1;
			}
		}
	}
	if (args->focal_length <= 0.0) {
		args->focal_length = fit.focal_length;
		if (args->focal_length <= 0.0) {
			// TODO: which one should we use here?
			args->focal_length = com.pref.starfinder_conf.focal_length;
			//args->focal_length = com.pref.focal;
			if (args->focal_length <= 0.0) {
				siril_log_color_message(_("Focal length not found in image or in settings, cannot proceed\n"), "red");
				return 1;
			}
		}
	}

	clearfits(&fit);
	return 0;
}

static int astrometry_image_hook(struct generic_seq_args *arg, int o, int i, fits *fit, rectangle *area, int threads) {
	struct astrometry_data *args = (struct astrometry_data *)arg->user;
	args->fit = fit;	// not reentrant because of that and other fields of the struct
	process_plate_solver_input(args);
	return GPOINTER_TO_INT(plate_solver(args));
}

/* TODO:
 * improvements:
 * - reuse the same catalog data, which will also allow parallelism
 * - write only image header
 */
void start_sequence_astrometry(sequence *seq, struct astrometry_data *args) {
	struct generic_seq_args *seqargs = create_default_seqargs(seq);
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seq->selnum;
	seqargs->stop_on_error = FALSE;
	seqargs->parallel = FALSE;		// I don't think local catalogues are reentrant
	seqargs->prepare_hook = astrometry_prepare_hook;
	seqargs->image_hook = astrometry_image_hook;
	seqargs->has_output = TRUE;
	seqargs->new_seq_prefix = "ps_";
	seqargs->description = "plate solving";
	seqargs->user = args;

	start_in_new_thread(generic_sequence_worker, seqargs);
}

