/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/sleef.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_world_cs.h"
#include "core/undo.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/annotate.h"
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
#include "gui/photometric_cc.h"

#include "astrometry_solver.h"

#define DOWNSAMPLE_FACTOR 0.25
#define CONV_TOLERANCE 1E-8

#undef DEBUG		/* get some of diagnostic output */

typedef struct {
	point size;
	SirilWorldCS *px_cat_center;
	SirilWorldCS *image_center;
	double crpix[2];
	double pixel_size, focal;
	Homography H;
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
static double get_fov_arcmin(double resolution, int rx, int ry) {
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

static void compute_mag_limit(struct astrometry_data *args) {
	if (args->auto_magnitude)
		args->limit_mag = compute_mag_limit_from_fov(args->used_fov / 60.0);
	else args->limit_mag = args->forced_magnitude;
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
	case NOMAD:
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
	case TYCHO2:
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
	case GAIA:
		url = g_string_append(url, "I/345/gaia2&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Gmag%20BPmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Gmag=<");
		url = g_string_append(url, mag);
		break;
	case GAIAEDR3:
		url = g_string_append(url, "I/350/gaiaedr3&-out.meta=-h-u-D&-out.add=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Gmag%20BPmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Gmag=<");
		url = g_string_append(url, mag);
		break;
	case PPMXL:
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
	case BRIGHT_STARS:
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
	case APASS: // for photometry only
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
	char *result = NULL, *error;
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
			error = siril_log_message(_("Fetch failed with code %ld for URL %s\n"), code, url);
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), error);
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
		siril_log_message(_("Error loading url: %s: %s\n"), url, error->message);
		g_clear_error(&error);
	}
	g_object_unref(file);
	return content;
}
#endif

gchar *search_in_catalogs(const gchar *object, query_server server) {
	GString *string_url;
	gchar *name = g_utf8_strup(object, -1);

	switch(server) {
	case 0:
		string_url = g_string_new(CDSSESAME);
		string_url = g_string_append(string_url, "/-oI/A?");
		string_url = g_string_append(string_url, name);
		break;
	case 1:
		string_url = g_string_new(VIZIERSESAME);
		string_url = g_string_append(string_url, "/-oI/A?");
		string_url = g_string_append(string_url, name);
		break;
	default:
	case 2:
		string_url = g_string_new(SIMBADSESAME);
		string_url = g_string_append(string_url, name);
		string_url = g_string_append(string_url, "';");
		break;

	}

	gchar *url = g_string_free(string_url, FALSE);
	gchar *cleaned_url = url_cleanup(url);
	gchar *result = fetch_url(cleaned_url);

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
			/* if file exist and user uses cache */
			if (error->code == G_IO_ERROR_EXISTS) {
				siril_log_message(_("Using already downloaded star catalogue\n"));
				g_clear_error(&error);
			} else {
				g_warning("%s\n", error->message);
				siril_log_color_message(_("Cannot create catalogue file %s for plate solving (%s)\n"), "red", g_file_peek_path(file), error->message);
				g_clear_error(&error);
				g_object_unref(file);
				return NULL;
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
	gchar *alpha, *delta;
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

	alpha = siril_world_cs_alpha_format(image->image_center, " %02dh%02dm%02ds");
	delta = siril_world_cs_delta_format(image->image_center, "%c%02d°%02d\'%02d\"");
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

static int read_GAIA_catalog(GInputStream *stream, psf_star **cstars, const gchar *version) {
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
	siril_log_message(_("Catalog Gaia %s size: %d objects\n"), version, i);
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
	default:
	case TYCHO2:
		return read_TYCHO2_catalog(stream, cstars);
	case NOMAD:
		return read_NOMAD_catalog(stream, cstars);
	case GAIA:
		return read_GAIA_catalog(stream, cstars, "DR2");
	case GAIAEDR3:
		return read_GAIA_catalog(stream, cstars, "EDR3");
	case PPMXL:
		return read_PPMXL_catalog(stream, cstars);
	case BRIGHT_STARS:
		return read_BRIGHT_STARS_catalog(stream, cstars);
	case APASS:
		return read_APASS_catalog(stream, cstars);
	}
}

static TRANS H_to_linear_TRANS(Homography H) {
	TRANS trans;

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

/* Determine if the image need to be flipped
 * 1 - is the image flipped? Need to compute the determinant of H
 * 2 - did the user asked for a flip?
 */
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

#define CHECK_FOR_CANCELLATION if (!get_thread_run()) { args->message = g_strdup(_("Cancelled")); args->ret = 1; goto clearup; }

/* entry point for plate solving */
gpointer match_catalog(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	GError *error = NULL;
	psf_star **cstars = NULL;
	int n_fit = 0, n_cat = 0;
	Homography H = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int max_trials = 0;
	GFile *catalog = NULL;
	GInputStream *input_stream = NULL;
	s_star *star_list_A = NULL, *star_list_B = NULL;
	fits fit_backup = { 0 };
	solve_results solution = { 0 };
	psf_star **stars = NULL;

	args->ret = 1;
	args->message = NULL;

	if (args->use_local_cat) {
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

	if (!args->catalog_file && !args->use_local_cat) {
		args->catalog_file = download_catalog(args->onlineCatalog, args->cat_center,
				args->used_fov * 0.5, args->limit_mag);
		if (!args->catalog_file) {
			args->message = g_strdup(_("Could not download the online star catalogue."));
			goto clearup;
		}
	}
	CHECK_FOR_CANCELLATION;

	if (args->downsample) {
		copyfits(args->fit, &fit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		cvResizeGaussian(args->fit, DOWNSAMPLE_FACTOR * args->fit->rx, DOWNSAMPLE_FACTOR * args->fit->ry, OPENCV_AREA);
	}

	if (!args->manual) {
		com.pref.starfinder_conf.pixel_size_x = com.pref.pitch;
		com.pref.starfinder_conf.focal_length = com.pref.focal;

		image im = { .fit = args->fit, .from_seq = NULL, .index_in_seq = -1 };

		stars = peaker(&im, 0, &com.pref.starfinder_conf, &n_fit, &(args->solvearea), FALSE, FALSE, MAX_STARS_FITTED, com.max_thread); // TODO: use good layer
		com.pref.starfinder_conf.pixel_size_x = 0.;
		com.pref.starfinder_conf.focal_length = 0.;
	} else {
		stars = com.stars;
		if (com.stars)
			while (com.stars[n_fit++]);
	}
	CHECK_FOR_CANCELLATION;

	if (!stars || n_fit < AT_MATCH_STARTN_LINEAR) {
		args->message = g_strdup_printf(_("There are not enough stars picked in the image. "
				"At least %d stars are needed."), AT_MATCH_STARTN_LINEAR);
		goto clearup;
	}
	siril_debug_print("using %d detected stars from image\n", n_fit);
	if (args->uncentered)
		max_trials = 20; //retry to converge if solve is done at an offset from the center

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
	catalog = g_file_new_for_path(args->catalogStars);
	input_stream = (GInputStream*) g_file_read(catalog, NULL, &error);
	if (!input_stream) {
		if (error != NULL) {
			args->message = g_strdup_printf(_("Could not load the star catalog (%s)."), error->message);
			g_clear_error(&error);
		}
		args->message = g_strdup_printf(_("Could not load the star catalog (%s)."), "generic error");
		goto clearup;
	}

	n_cat = read_catalog(input_stream, cstars, args->onlineCatalog);
	CHECK_FOR_CANCELLATION;

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
	solution.px_cat_center = siril_world_cs_ref(args->cat_center);
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
	solution.crpix[0] = ((args->fit->rx - 1) / 2.0);
	solution.crpix[1] = ((args->fit->ry - 1) / 2.0);
	if (args->downsample) {
		clearfits(args->fit);
		memcpy(args->fit, &fit_backup, sizeof(fits));
		memset(&fit_backup, 0, sizeof(fits));
	}
	solution.size.x = args->fit->rx;
	solution.size.y = args->fit->ry;
	solution.pixel_size = args->pixel_size;

	apply_match(solution.px_cat_center, solution.crpix, trans, &ra0, &dec0);
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
		solution.px_cat_center = siril_world_cs_new_from_a_d(ra0, dec0);

		project_starlist(num_matched, star_list_B, ra0, dec0, 1);
		siril_debug_print("Reprojecting to: alpha: %s, delta: %s\n",
				siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
				siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));

		double scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
		double scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
		double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels

		double focal = RADCONV * solution.pixel_size / resolution;
		siril_debug_print("Current focal: %0.2fmm\n", focal);

		if (atPrepareHomography(num_matched, star_list_A, num_matched, star_list_B, &H, FALSE, NULL, NULL, AFFINE_TRANSFORMATION)){
			args->message = g_strdup(_("Updating homography failed."));
			args->ret = 1;
			break;
		}
		trans = H_to_linear_TRANS(H);
		apply_match(solution.px_cat_center, solution.crpix, trans, &ra0, &dec0);

		conv = fabs((dec0 - orig_dec0) / orig_dec0) + fabs((ra0 - orig_ra0) / orig_ra0);

		trial++;
		CHECK_FOR_CANCELLATION;
	}
	if (args->ret)	// after the break
		goto clearup;

	memcpy(&solution.H, &H, sizeof(Homography));
	double scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	double scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels
	solution.focal = RADCONV * solution.pixel_size / resolution;

	solution.image_center = siril_world_cs_new_from_a_d(ra0, dec0);
	if (max_trials == 0) {
		siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
	} else if (trial == max_trials) {
		siril_debug_print("No convergence found: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
	} else {
		siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f at iteration #%d\n", ra0, dec0, trial);
	}

	double scalefactor = (args->downsample) ? 1.0 / DOWNSAMPLE_FACTOR : 1.0;
	if (args->downsample)
		solution.focal *= scalefactor;

	CHECK_FOR_CANCELLATION;
	/* compute cd matrix */
	double ra7, dec7, delta_ra;

	/* first, convert center coordinates from deg to rad: */
	dec0 *= DEGTORAD;
	ra0 *= DEGTORAD;

	/* make 1 step in direction crpix1 */
	double crpix1[] = { solution.crpix[0] + 1.0 / scalefactor, solution.crpix[1] };
	apply_match(solution.px_cat_center, crpix1, trans, &ra7, &dec7);

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
	double crpix2[] = { solution.crpix[0], solution.crpix[1] - 1.0 / scalefactor };
	apply_match(solution.px_cat_center, crpix2, trans, &ra7, &dec7);

	dec7 *= DEGTORAD;
	ra7 *= DEGTORAD;

	delta_ra = ra7 - ra0;
	if (delta_ra > +M_PI)
		delta_ra = 2.0 * M_PI - delta_ra;
	if (delta_ra < -M_PI)
		delta_ra = delta_ra - 2.0 * M_PI;
	double cd1_2 = (delta_ra) * cos(dec0) * RADTODEG;
	double cd2_2 = (dec7 - dec0) * RADTODEG;

	// saving state for undo before modifying fit structure
	if (!com.script) {
		const char *undo_str = args->for_photometry_cc ? _("Photometric CC") : _("Plate Solve");
		undo_save_state(args->fit, undo_str);
	}

	/**** Fill wcsdata fit structure ***/

	args->fit->wcsdata.equinox = 2000.0;
	args->fit->focal_length = solution.focal;
	args->fit->pixel_size_x = args->fit->pixel_size_y = solution.pixel_size;
	solution.crpix[0] *= scalefactor;
	solution.crpix[1] *= scalefactor;

	args->fit->wcsdata.crpix[0] = solution.crpix[0];
	args->fit->wcsdata.crpix[1] = solution.crpix[1];
	args->fit->wcsdata.crval[0] = ra0 * RADTODEG;
	args->fit->wcsdata.crval[1] = dec0 * RADTODEG;

	args->fit->wcsdata.ra = siril_world_cs_get_alpha(solution.image_center);
	args->fit->wcsdata.dec = siril_world_cs_get_delta(solution.image_center);

	args->fit->wcsdata.pltsolvd = TRUE;
	g_snprintf(args->fit->wcsdata.pltsolvd_comment, 21, "Siril internal solver");

	gchar *ra = siril_world_cs_alpha_format(solution.image_center, "%02d %02d %.3lf");
	gchar *dec = siril_world_cs_delta_format(solution.image_center, "%c%02d %02d %.3lf");

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
	siril_debug_print("crpix1 = %*.12e\n", 20, solution.crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, solution.crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, ra0 * RADTODEG);
	siril_debug_print("crval2 = %*.12e\n", 20, dec0 * RADTODEG);
	siril_debug_print("cdelt1 = %*.12e\n", 20, cdelt1);
	siril_debug_print("cdelt2 = %*.12e\n", 20, cdelt2);
	siril_debug_print("pc1_1  = %*.12e\n", 20, args->fit->wcsdata.pc[0][0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, args->fit->wcsdata.pc[0][1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, args->fit->wcsdata.pc[1][0]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, args->fit->wcsdata.pc[1][1]);
	siril_debug_print("******************************************\n");

	load_WCS_from_memory(args->fit);
	print_platesolving_results(&solution, args->downsample);
	if (args->for_photometry_cc) {
		pcc_star *pcc_stars = NULL;
		int nb_pcc_stars;
		if (args->use_local_cat) {
			double tra = siril_world_cs_get_alpha(solution.image_center);
			double tdec = siril_world_cs_get_delta(solution.image_center);
			double res = get_resolution(solution.focal, args->pixel_size);
			double radius = get_radius_deg(res, args->fit->rx, args->fit->ry);
			double mag = compute_mag_limit_from_fov(radius * 1.5);
			// for photometry, we can use fainter stars, 1.5 seems ok above instead of 2.0
			siril_log_message(_("Getting stars from local catalogues for PCC, limit magnitude %.2f\n"), mag);
			if (get_stars_from_local_catalogues(tra, tdec, radius, args->fit, mag, &pcc_stars, &nb_pcc_stars)) {
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
	if (args->flip_image && image_is_flipped(H)) {
		siril_log_color_message(_("Flipping image and updating astrometry data.\n"), "salmon");
		fits_flip_top_to_bottom(args->fit);
		flip_bottom_up_astrometry_data(args->fit);
		load_WCS_from_memory(args->fit);
		args->image_flipped = TRUE;
	}
	/* free solution */
	if (solution.px_cat_center)
		siril_world_cs_unref(solution.px_cat_center);
	args->new_center = solution.image_center;

clearup:
	if (fit_backup.rx != 0) {
		clearfits(args->fit);
		memcpy(args->fit, &fit_backup, sizeof(fits));
	}
	if (stars && !args->manual) {
		for (int i = 0; i < n_fit; i++)
			free_psf(stars[i]);
		free(stars);
	}
	free_stars(&star_list_A);
	free_stars(&star_list_B);
	siril_world_cs_unref(args->cat_center);
	if (cstars)
		free_fitted_stars(cstars);
	if (input_stream)
		g_object_unref(input_stream);
	if (catalog)
		g_object_unref(catalog);
	if (args->catalog_file)
		g_object_unref(args->catalog_file);
	g_free(args->catalogStars);

	siril_add_idle(end_plate_solver, args);
	int retval = args->ret;
	if (com.script) {
		if (args->ret)
			siril_log_message(_("Plate solving failed: %s\n"), args->message);
		free(args);
	}
	return GINT_TO_POINTER(retval);
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

	compute_mag_limit(args); // to call after having set args->used_fov
}
