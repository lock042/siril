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

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <unistd.h>
#endif
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "algos/search_objects.h"
#include "algos/siril_wcs.h"
#include "algos/annotate.h"
#include "algos/astrometry_solver.h"
#include "algos/comparison_stars.h"
#include "io/remote_catalogues.h"
#include "registration/matching/misc.h"
#include "gui/message_dialog.h"	// why a siril_message_dialog() call in here?


/* There are several ways of obtaining star catalogue data in Siril.
 * The raw catalogue contain in general star name, RA and dec coords, V and B
 * magnitudes. The download_catalog function makes a request to an online
 * source with a centre, a radius and a limit magnitude and stores that in a
 * raw cache file named 'cat-type-ra-dec-radius-mag.cat'.
 *
 * Then these raw catalogues can be used in different ways:
 * - The astrometric solver reads them and projects stars on a 2D image
 *   centered on the catalogue center and saves that to another file (the
 *   function project_catalog does that) called catalog.proj. Then the function
 *   read_projected_catalog below reads this intermediary projection catalogue
 *   and stores that in an array of psf_star objects. The obtained stars can be
 *   used for registration, but do not correspond to image coordinates.
 *
 * - The PCC reads them and projects stars on a plate-solved image using WCS
 *   and stores them in condensed form (pcc_star struct containing only
 *   x,y,b,v), done in the function project_catalog_with_WCS
 *
 * - Comparison star list creation needs equatorial coordinates and B-V
 *   magnitudes, projection is also used but only to check if a star is inside
 *   the image and far from its borders. The required object is psf_star.
 */

const char *catalog_to_str(online_catalog cat) {
	switch (cat) {
		case CAT_TYCHO2:
			return _("Tycho-2");
		case CAT_NOMAD:
			return _("NOMAD");
		case CAT_GAIADR3:
			return _("Gaia DR3");
		case CAT_PPMXL:
			return _("PPMXL");
		case CAT_BRIGHT_STARS:
			return _("bright stars");
		case CAT_APASS:
			return _("APASS");
		case CAT_AAVSO:
			return _("AAVSO");
		case CAT_LOCAL:
			return _("local Tycho-2+NOMAD");
		case CAT_ASNET:
			return _("local astrometry.net");
		default:
			return _("unknown");
	}
}

/*
 *                     _ _                                   _ _
 *  _ __ ___  __ _  __| (_)_ __   __ _   _ __ ___  ___ _   _| | |_ ___
 * | '__/ _ \/ _` |/ _` | | '_ \ / _` | | '__/ _ \/ __| | | | | __/ __|
 * | | |  __/ (_| | (_| | | | | | (_| | | | |  __/\__ \ |_| | | |_\__ \
 * |_|  \___|\__,_|\__,_|_|_| |_|\__, | |_|  \___||___/\__,_|_|\__|___/
 *                               |___/
 *
 * responses of online catalogue requests are saved to files for caching then
 * projected to 2D coordinates, this section is about reading those files in
 * which equatorial coordinates are replaced by approximative pixel coordinate
 * TODO: maybe only one function can read them all...
 */

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
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0, e_Vmag = 0.0, e_Bmag = 0.0;

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

		int n = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag, &e_Vmag, &e_Bmag);
		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->Bmag = Bmag;
		star->s_mag = e_Vmag;
		star->s_Bmag = e_Bmag;
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

int read_projected_catalog(GInputStream *stream, psf_star **cstars, online_catalog cat) {
	switch (cat) {
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
static const int DEFAULT_FETCH_RETRIES = 10;

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
			siril_log_message(_("Failed to download page %s (error %ld)\n"), url, code);

			if (retries) {
				s = 2 * (DEFAULT_FETCH_RETRIES - retries) + 2;
				siril_debug_print("Wait %us before retry\n", s);
				g_usleep(s * 1E6);

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
			// TODO: this GUI code happens only for curl, not the alternative, to remove?
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
		siril_log_message(_("Error loading url: [%s] - %s\n"), url, error->message);
		g_clear_error(&error);
	}
	g_object_unref(file);
	return content;
}

void free_fetch_result(gchar *result) {
	g_free(result);
}
#endif

gchar *get_catalog_url(SirilWorldCS *center, double mag_limit, double radius, int type) {
	GString *url;
	gchar *coordinates = g_strdup_printf("%f+%f", siril_world_cs_get_alpha(center), siril_world_cs_get_delta(center));
	gchar *mag = g_strdup_printf("%2.2lf", mag_limit);
	gchar *fov = g_strdup_printf("%2.1lf", radius);

	switch (type) {
	case CAT_NOMAD:
		url = g_string_new(VIZIER_QUERY);
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
		url = g_string_new(VIZIER_QUERY);
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
		url = g_string_new(VIZIER_QUERY);
		url = g_string_append(url, "I/355/gaiadr3&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Gmag%20BPmag%20e_Gmag%20e_BPmag%20DR3Name");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		/*if (com.target_star) {
			g_string_append_printf(url, "&Gmag=<%2.3f", com.target_star->mag + com.delta_vmag);
			g_string_append_printf(url, "&Gmag=>%2.3f", com.target_star->mag - com.delta_vmag);
		}*/
		break;
	case CAT_PPMXL:
		url = g_string_new(VIZIER_QUERY);
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
		url = g_string_new(VIZIER_QUERY);
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
		url = g_string_new(VIZIER_QUERY);
		// TODO: modify to get the correct response
		url = g_string_append(url, "APASS&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Vmag%20Bmag%20e_Vmag%20e_Bmag%20recno");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	case CAT_AAVSO: // for photometry only
		{
			int ra_h, ra_m;
			double ra_s;
			siril_world_cs_get_ra_hour_min_sec(center, &ra_h, &ra_m, &ra_s);
			int dec_deg, dec_m;
			double dec_s;
			siril_world_cs_get_dec_deg_min_sec(center, &dec_deg, &dec_m, &dec_s);
			url = g_string_new(AAVSO_QUERY);
			g_string_append_printf(url, "ra=%02d:%02d:%02.lf", ra_h, ra_m, ra_s);
			g_string_append_printf(url, "&dec=%02d:%02d:%02.lf", dec_deg, dec_m, dec_s);
			url = g_string_append(url, "&fov=");
			url = g_string_append(url, fov);
			g_string_append_printf(url, "&maglimit=%s&format=json", mag);
		}
		break;
	}

	g_free(coordinates);
	g_free(mag);
	g_free(fov);

	return g_string_free(url, FALSE);
}

gpointer search_in_online_conesearch(gpointer p) {
#ifndef HAVE_NETWORKING
	siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
	return GINT_TO_POINTER(1);
#else
	struct astrometry_data *args = (struct astrometry_data *) p;
	if (!args->fit->date_obs) {
		free(args);
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(-1);
	}
	double ra, dec;
	center2wcs(args->fit, &ra, &dec);
	int retval = 0;

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

	if (!gfit.date_obs) {
		siril_log_color_message(_("This command only works on images that have observation date information\n"), "red");
		return NULL;
	}
	siril_log_message(_("Solar System Objects search on observation date %s\n"), formatted_date);

	gchar *url = g_string_free(string_url, FALSE);
	gchar *cleaned_url = url_cleanup(url);
	gchar *result = fetch_url(cleaned_url);
	siril_debug_print(_("URL: %s\n"), cleaned_url);

	g_free(cleaned_url);
	g_free(url);
	g_free(formatted_date);

	if (result) {
		retval = parse_conesearch_buffer(result, args->limit_mag);
	}
#if defined HAVE_LIBCURL
	free(result);
#else
	g_free(result);
#endif
	if (!retval) {
		siril_add_idle(end_process_sso, args);
	} else {
		free(args);
		siril_add_idle(end_generic, NULL);
	}

	return GINT_TO_POINTER(retval);
#endif
}


GFile *download_catalog(online_catalog onlineCatalog, SirilWorldCS *catalog_center, double radius_arcmin, double mag) {
#ifndef HAVE_NETWORKING
	siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
#else
	gchar *buffer = NULL;
	GError *error = NULL;
	/* ---------------- get Vizier catalog in a cache catalogue file --------------------- */

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
#endif
	return NULL;
}


/*
 *      _ _               _                        _
 *   __| (_)_ __ ___  ___| |_   _ __ ___  __ _  __| |
 *  / _` | | '__/ _ \/ __| __| | '__/ _ \/ _` |/ _` |
 * | (_| | | | |  __/ (__| |_  | | |  __/ (_| | (_| |
 *  \__,_|_|_|  \___|\___|\__| |_|  \___|\__,_|\__,_|
 *
 * of catalogue files to memory representations (pcc_star or psf_star)
 */

int project_catalog_with_WCS(GFile *catalog_file, fits *fit, gboolean phot, pcc_star **ret_stars, int *ret_nb_stars) {
	GError *error = NULL;
	GInputStream *input_stream = NULL;
	/* catalog format should be 5 columns: distance from centre, RA, Dec, V, B */
	if (!(input_stream = (GInputStream*) g_file_read(catalog_file, NULL, &error))) {
		if (error) {
			siril_log_message(_("Can't open catalog file %s for PCC: %s\n"), g_file_peek_path(catalog_file), error->message);
			g_clear_error(&error);
		}
		*ret_stars = NULL;
		*ret_nb_stars = 0;
		return 1;
	}

	int nb_alloc = 1200, nb_stars = 0;
	pcc_star *stars = malloc(nb_alloc * sizeof(pcc_star));

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	gchar *line;
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		if (line[0] == COMMENT_CHAR || is_blank(line) || g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		double r = 0.0, ra = 0.0, dec = 0.0, Vmag = 0.0, Bmag = 0.0;
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &ra, &dec, &Vmag, &Bmag);
		g_free(line);
		if ((phot && n == 5 && Bmag < 30.0) ||	// 30 sometimes means not available in NOMAD
				(!phot && n >= 4)) {
			if (nb_stars >= nb_alloc) {
				nb_alloc *= 2;
				pcc_star *new_array = realloc(stars, nb_alloc * sizeof(pcc_star));
				if (!new_array) {
					PRINT_ALLOC_ERR;
					g_object_unref(data_input);
					g_object_unref(input_stream);
					free(stars);
					*ret_stars = NULL;
					*ret_nb_stars = 0;
					return 1;
				}
				stars = new_array;
			}

			double x, y;
			if (!wcs2pix(fit, ra, dec, &x, &y)) {
				stars[nb_stars].x = x;
				stars[nb_stars].y = y;
				stars[nb_stars].mag = Vmag;
				if (phot)
					stars[nb_stars].BV = Bmag - Vmag;
				else stars[nb_stars].BV = -99.9;
				nb_stars++;
			}
		}
	}

	g_object_unref(data_input);
	siril_debug_print("projected %d stars from the provided catalogue\n", nb_stars);
	*ret_stars = stars;
	*ret_nb_stars = nb_stars;
	return 0;
}

// simply load a catalogue from its downloaded cache file to a list of psf_star *
// if phot is true, only stars with correct B and V values are added to the result
int load_catalog(GFile *catalog_file, gboolean phot, psf_star **ret_stars, int *ret_nb_stars) {
	GError *error = NULL;
	GInputStream *input_stream = (GInputStream*) g_file_read(catalog_file, NULL, &error);
	if (!input_stream) {
		if (error != NULL) {
			siril_log_message(_("Could not load the star catalog (%s)."), error->message);
			g_clear_error(&error);
		}
		siril_log_message(_("Could not load the star catalog (%s)."), "generic error");
		return 1;
	}

	int nb_alloc = 1200, nb_stars = 0;
	psf_star *stars = malloc(nb_alloc * sizeof(psf_star));

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	gchar *line;
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		if (line[0] == COMMENT_CHAR || is_blank(line) || g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		// APASS fields: RAJ2000 DEJ2000 Vmag Bmag e_Vmag e_Bmag recno
		// NOMAD fields: RAJ2000 DEJ2000 Vmag Bmag
		// GAIA3 fields: RAJ2000 DEJ2000 Gmag BPmag e_Gmag e_BPmag DR3Name
		int size;
		double r = 0.0, ra = 0.0, dec = 0.0, Vmag = 0.0, Bmag = 0.0, e_Vmag = 0.0, e_Bmag = 0.0;
		int n = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %n",
				&r, &ra, &dec, &Vmag, &Bmag, &e_Vmag, &e_Bmag, &size);
		gchar *name = NULL;
		if (n == 7)
			name = line + size;
		// magnitudes above 30 are a code for 'undefined'
		if ((!phot && n >= 4) || (phot && n >= 5 && Vmag < 30.0 && Bmag < 30.0)) {
			// 30 sometimes means not available (NOMAD)
			if (nb_stars >= nb_alloc) {
				nb_alloc *= 2;
				psf_star *new_array = realloc(stars, nb_alloc * sizeof(psf_star));
				if (!new_array) {
					PRINT_ALLOC_ERR;
					g_object_unref(data_input);
					g_object_unref(input_stream);
					free(stars);
					*ret_stars = NULL;
					*ret_nb_stars = 0;
					return 1;
				}
				stars = new_array;
			}

			psf_star_init(&stars[nb_stars]);
			stars[nb_stars].star_name = g_strdup(name);
			stars[nb_stars].ra = ra;
			stars[nb_stars].dec = dec;
			stars[nb_stars].mag = Vmag;
			stars[nb_stars].Bmag = Bmag;
			stars[nb_stars].s_mag = e_Vmag;
			stars[nb_stars].s_Bmag = e_Bmag;
			nb_stars++;
		}
		g_free(line);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);
	*ret_stars = stars;	// warning: what's after nb_stars is not initialized
	*ret_nb_stars = nb_stars;
	siril_debug_print("read %d%s stars from catalogue\n", nb_stars, phot ? " photometric" : "");
	return 0;
}

// TODO to be reviewed
// uses gfit for an is_inside() check
int read_photo_aavso_buffer(const char *buffer, struct compstars_arg *args) {
	gchar **token = g_strsplit(buffer, "auid", -1);
	int nargs = g_strv_length(token);
	args->cat_stars = malloc((min(MAX_STARS, nargs) + 1) * sizeof(psf_star *));
	int nb_stars = 0;
	char chartid[28];

	if (nargs == 0) return 0;

	if (g_str_has_prefix(token[0], "{\"chartid\":\"")) {
		gchar **fields = g_strsplit(token[0], "\":\"", -1);
		sscanf(fields[0], "%s", chartid);
		g_strfreev (fields);
	}

	for (int i = 0; i < nargs; i++) {
		double ra = 0.0, dec = 0.0, Vmag = 0.0, Bmag = 0.0, e_Vmag = 0.0, e_Bmag = 0.0;
		double hours = 0.0, min = 0.0, seconds = 0.0, degres = 0.0;
		char sname[28];
		char chartid[28];
		char star_uri[128];

		// Fields containing the "chart id" and the chart URL according to the AAVSO
		// TODO: why here and above?
		if (g_str_has_prefix(token[i], "{\"chartid\":\"")) {
			gchar **fields = g_strsplit(token[i], "\":\"", -1);
			gchar **subfields = g_strsplit(fields[1], "\",\"", -1);
			// TODO: string length checks to not overflow buffer
			sscanf(subfields[0], "%s", chartid);		// chartid is a fixed 8 caracters string
			// TODO: g_strfreev required
			subfields = g_strsplit(fields[2], "\",\"", -1);
			sscanf(subfields[0], "%s", star_uri);		// star_uri is a fixed 56 caracters string
			g_strfreev (subfields);
			g_strfreev (fields);
		}

		gchar **fields = g_strsplit(token[i], "\":\"", -1);

		if (g_str_has_prefix(token[i], "\":\"")) {
			int ind = 1;
			while (fields[ind]) {
				gchar **part = g_strsplit(fields[ind], "\",\"", -1);

				// Field containing the name
				if (ind == 1)	sscanf(part[0], "%s", sname);
				// Field containing the RA
				if (ind == 2) {
					gchar **slice = g_strsplit(part[0], ":", -1);
					sscanf(slice[0], "%lf", &hours);
					sscanf(slice[1], "%lf", &min);
					sscanf(slice[2], "%lf", &seconds);
					ra = 360.0*hours/24.0 + 360.0*min/(24.0*60.0) + 360.0*seconds/(24.0*60.0*60.0);
					g_free(slice);
				}
				// Field containing the DEC
				if (ind == 3) {
					gchar **slice = g_strsplit(part[0], ":", -1);
					sscanf(slice[0], "%lf", &degres);
					sscanf(slice[1], "%lf", &min);
					sscanf(slice[2], "%lf", &seconds);
					if (degres < 0.0)
						dec = degres - min/60.0 - seconds/3600.0;
					else dec = degres + min/60.0 + seconds/3600.0;
					g_free(slice);
				}
				// Field containing Vmag/Bmag and e_Vmag/e_Bmag
				if (ind >= 4 && ind <= 5) {
					gchar *mag = NULL;
					char *start = strchr(fields[ind], ':');
					char *end = strstr(start, ",\"e");
					mag = g_strndup(start + 1, end - start - 3);
					if (ind == 4) sscanf(mag, "%lf", &Vmag);
					if (ind == 5) sscanf(mag, "%lf", &Bmag);
					g_free(mag);

					gchar *e_mag = NULL;
					start = strstr(fields[ind], "r\":");
					end = strrchr(start, '}');
					e_mag = g_strndup(start + 3, end - start - 3);
					if (ind == 4) sscanf(e_mag, "%lf", &e_Vmag);
					if (ind == 5) sscanf(e_mag, "%lf", &e_Bmag);
					g_free(e_mag);
				}
				g_free(part);
				ind++;
			}
		}
		g_strfreev (fields);
		double x, y;	// for display purposes only
		if (!is_inside2(&gfit, ra, dec, &x, &y))
			continue;
		args->cat_stars[nb_stars].star_name = g_strdup(sname);
		args->cat_stars[nb_stars].xpos = x;
		args->cat_stars[nb_stars].ypos = gfit.ry - y - 1;	// why the hell do we need that here?
		args->cat_stars[nb_stars].ra = ra;
		args->cat_stars[nb_stars].dec = dec;
		args->cat_stars[nb_stars].mag = Vmag;
		args->cat_stars[nb_stars].Bmag = Bmag;
		args->cat_stars[nb_stars].s_mag = e_Vmag;
		args->cat_stars[nb_stars].s_Bmag = e_Bmag;
		args->cat_stars[nb_stars].BV = Bmag - Vmag;
		args->cat_stars[nb_stars].phot = NULL;
		args->AAVSO_uri = g_strdup(star_uri);
		args->AAVSO_chartid = g_strdup(chartid);
		siril_debug_print(_("name: %s, ra: %lf, dec: %lf, Vmag: %lf, e_Vmag, %lf, Bmag: %lf, e_Bmag: %lf\n"),
				g_strdup(sname), ra, dec, Vmag, e_Vmag, Bmag, e_Bmag);
		if (nb_stars == MAX_STARS)
			break;
	}
	g_strfreev(token);
	return nb_stars;
}

