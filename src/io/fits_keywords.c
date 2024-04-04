/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "core/siril_world_cs.h"
#include "algos/siril_wcs.h"
#include "io/image_format_fits.h"

#include "fits_keywords.h"

#define DEFAULT_DOUBLE_VALUE -999.0
#define DEFAULT_INT_VALUE -INT_MAX //FIXME: 0 ?
#define DEFAULT_UINT_VALUE 0
#define DEFAULT_USHORT_VALUE DEFAULT_UINT_VALUE

#define KEYWORD(group, key, type, comment, data, handler) { group, key, type, comment, data, handler, TRUE, FALSE }
#define KEYWORD_FIXED(group, key, type, comment, data, handler) { group, key, type, comment, data, handler, TRUE, TRUE }
#define KEYWORD_WCS(group, key, type) { group, key, type, NULL, NULL, NULL, FALSE, TRUE }

static gboolean should_use_keyword(const fits *fit, const gchar *group, const gchar *keyword) {

	if (g_strcmp0(keyword, "ROWORDER") == 0) {
		return ((g_strcmp0(fit->keywords.row_order, "BOTTOM-UP") == 0)
				|| (g_strcmp0(fit->keywords.row_order, "TOP-DOWN") == 0));
	} else if (g_strcmp0(keyword, "FLENGTH") == 0) {
		return FALSE;
	} else if (g_strcmp0(keyword, "XBAYROFF") == 0) {
        return fit->keywords.bayer_pattern[0] != '\0';
    } else if (g_strcmp0(keyword, "YBAYROFF") == 0) {
        return fit->keywords.bayer_pattern[0] != '\0';
    } else if (g_strcmp0(keyword, "DFTNORM2") == 0) {
        return fit->naxes[2] > 1;
    } else if (g_strcmp0(keyword, "DFTNORM3") == 0) {
        return fit->naxes[2] > 1;
    } else if (g_strcmp0(keyword, "RA_D") == 0) {
        return FALSE;
    } else if (g_strcmp0(keyword, "DEC_D") == 0) {
        return FALSE;
    }
    return TRUE;
}

/****************************** handlers ******************************/

static void pixel_x_handler(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.pixel_size_x > 0.0) {
		fit->pixelkey = TRUE;
	}
}

static void binning_x_handler(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.binning_x <= 0)
		fit->keywords.binning_x = 1;
}

static void binning_y_handler(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.binning_y <= 0)
		fit->keywords.binning_y = 1;
}

static void roworder_handler(fits *fit, const char *comment, KeywordInfo *info) {
	if (!strcasecmp(fit->keywords.bayer_pattern, "NONE")) {
		memset(fit->keywords.bayer_pattern, 0, sizeof(char) * FLEN_VALUE);
	}
}

static void focal_length_handler(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.focal_length > 0.0)
		fit->focalkey = TRUE;
}

static void flength_handler(fits *fit, const char *comment, KeywordInfo *info) {
	fit->keywords.focal_length = fit->keywords.flength * 1000.0; // convert m to mm
	fit->focalkey = TRUE;
}

static void sitelong_handler(fits *fit, const char *comment, KeywordInfo *info) {
	char sitelong_dump[FLEN_VALUE] = { 0 };
	char sitelong_dump_tmp[FLEN_VALUE] = { 0 };
	double d_sitelong_dump = 0.0;

	gchar **token = g_strsplit(sitelong_dump, ":", -1); // Handles PRISM special parsing for SITELONG
	gsize token_size = g_strv_length(token);
	if (token_size > 1 && token[1])	{
		for (int i = 0; i < token_size; ++i) {
			g_strlcat(sitelong_dump_tmp, token[i], sizeof(sitelong_dump_tmp));
			if (i < 3) strncat(sitelong_dump_tmp, i < 2 ? ":" : ".", 2);
			d_sitelong_dump = parse_dms(sitelong_dump_tmp);
		}
	} else d_sitelong_dump = parse_dms(sitelong_dump);

	g_strfreev(token);

	if (isnan(d_sitelong_dump)) {	// Cases SITELONG and SITELAT keyword are numbers (only NINA and Seq. Generator, for now)
		gchar *end;
		fit->keywords.sitelong = g_ascii_strtod(sitelong_dump, &end);
		if (sitelong_dump == end) {
			siril_debug_print("Cannot read SITELONG\n");
		}
	} else {
		fit->keywords.sitelong = d_sitelong_dump;
	}
}

static void sitelat_handler(fits *fit, const char *comment, KeywordInfo *info) {
	char sitelat_dump[FLEN_VALUE] = { 0 };
	char sitelat_dump_tmp[FLEN_VALUE] = { 0 };
	double d_sitelat_dump = 0.0;

	gchar **token = g_strsplit(sitelat_dump, ":", -1); // Handles PRISM special parsing for SITELAT
	gsize token_size = g_strv_length(token);
	if (token_size > 1 && token[1])	{	// Denotes presence of ":"
		for (int i = 0; i < token_size; ++i) {
			g_strlcat(sitelat_dump_tmp, token[i], sizeof(sitelat_dump_tmp));
			if (i < 3) strncat(sitelat_dump_tmp, i < 2 ? ":" : ".", 2);
			d_sitelat_dump = parse_dms(sitelat_dump_tmp);
		}
	} else d_sitelat_dump = parse_dms(sitelat_dump);

	g_strfreev(token);

	if (isnan(d_sitelat_dump)) {	// Cases SITELONG and SITELAT keyword are numbers (only NINA and Seq. Generator, for now)
		gchar *end;
		fit->keywords.sitelat = g_ascii_strtod(sitelat_dump, &end);
		if (sitelat_dump == end) {
			siril_debug_print("Cannot read SITELAT\n");
		}
	} else {
		fit->keywords.sitelat = d_sitelat_dump;
	}
}

void pltsolvd_handler(fits *fit, const char *comment, KeywordInfo *info) {
	strncpy(fit->keywords.wcsdata.pltsolvd_comment, comment, FLEN_COMMENT);
}

void datamax_handler(fits *fit, const char *comment, KeywordInfo *info) {
	gboolean not_from_siril = (strstr(fit->keywords.program, "Siril") == NULL);
	if ((fit->bitpix == FLOAT_IMG && not_from_siril) || fit->bitpix == DOUBLE_IMG) {
		float mini, maxi;
		fit_stats(fit->fptr, &mini, &maxi);
		// override data_max if needed. In some images there are differences between max and data_max
		fit->keywords.data_max = (double) maxi;
		fit->keywords.data_min = (double) mini;
	}
}

/*****************************************************************************/

KeywordInfo *initialize_keywords(fits *fit) {
	KeywordInfo keyword_list[] = {
		// FIXME: add MIPS keywords
        KEYWORD( "image", "MIPS-HI;CWHITE", KTYPE_USHORT, "Upper visualization cutoff", &(fit->keywords.hi), NULL),
        KEYWORD( "image", "MIPS-LO;CBLACK", KTYPE_USHORT, "Lower visualization cutoff", &(fit->keywords.lo), NULL),
		/* IMPORTANT: PROGRAM MUST BE BEFORE DATAMAX */
        KEYWORD( "image", "PROGRAM", KTYPE_STR, "Software that created this HDU", &(fit->keywords.program), NULL),
        KEYWORD( "image", "DATAMAX", KTYPE_DOUBLE, "Order of the rows in image array", &(fit->keywords.data_max), datamax_handler),
        KEYWORD( "image", "ROWORDER", KTYPE_STR, "Order of the rows in image array", &(fit->keywords.row_order), roworder_handler),
        KEYWORD( "setup", "INSTRUME", KTYPE_STR, "Instrument name", &(fit->keywords.instrume), NULL),
        KEYWORD( "setup", "TELESCOP", KTYPE_STR, "Telescope used to acquire this image", &(fit->keywords.telescop), NULL),
        KEYWORD( "setup", "OBSERVER", KTYPE_STR, "Observer name", &(fit->keywords.observer), NULL),
        KEYWORD( "date",  "DATE", KTYPE_DATE, "UTC date that FITS file was created", &(fit->keywords.date), NULL),
        KEYWORD( "date",  "DATE-OBS", KTYPE_DATE, "YYYY-MM-DDThh:mm:ss observation start, UT", &(fit->keywords.date_obs), NULL),
        KEYWORD( "image", "STACKCNT;NCOMBINE", KTYPE_UINT, "Stack frames", &(fit->keywords.stackcnt), NULL),
        KEYWORD( "image", "EXPTIME;EXPOSURE", KTYPE_DOUBLE, "Exposure time [s]", &(fit->keywords.exposure), NULL),
        KEYWORD( "image", "LIVETIME", KTYPE_DOUBLE, "Exposure time after deadtime correction", &(fit->keywords.livetime), NULL),
        KEYWORD( "image", "EXPSTART", KTYPE_DOUBLE, "Exposure start time (standard Julian date)", &(fit->keywords.expstart), NULL),
        KEYWORD( "image", "EXPEND", KTYPE_DOUBLE, "Exposure end time (standard Julian date)", &(fit->keywords.expend), NULL),
        KEYWORD( "image", "XPIXSZ;XPIXELSZ;PIXSIZE1;PIXSIZEX;XPIXSIZE", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler),
        KEYWORD( "image", "YPIXSZ;YPIXELSZ;PIXSIZE2;PIXSIZEY;YPIXSIZE", KTYPE_DOUBLE, "Y pixel size microns", &(fit->keywords.pixel_size_y), NULL),
        KEYWORD( "image", "XBINNING;BINX", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_x), binning_x_handler),
        KEYWORD( "image", "YBINNING;BINY", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_y), binning_y_handler),
        KEYWORD( "image", "FOCALLEN;FOCAL", KTYPE_DOUBLE, "Camera focal length", &(fit->keywords.focal_length), focal_length_handler),
        KEYWORD( "image", "FLENGTH", KTYPE_DOUBLE, "Camera focal length [m]", &(fit->keywords.flength), flength_handler),
        KEYWORD( "image", "CCD-TEMP;CCD_TEMP;CCDTEMP;TEMPERAT;CAMTCCD", KTYPE_DOUBLE, "CCD temp in C", &(fit->keywords.ccd_temp), NULL),
        KEYWORD( "image", "SET-TEMP", KTYPE_DOUBLE, "Temperature setting in C", &(fit->keywords.set_temp), NULL),
        KEYWORD( "image", "FILTER;FILT-1", KTYPE_STR, "Active filter name", &(fit->keywords.filter), NULL),
        KEYWORD( "image", "IMAGETYP;FRAMETYP", KTYPE_STR, "Type of image", &(fit->keywords.image_type), NULL),
        KEYWORD( "image", "OBJECT", KTYPE_STR, "Name of the object of interest", &(fit->keywords.object), NULL),
        KEYWORD( "image", "APERTURE", KTYPE_DOUBLE, "Aperture of the instrument", &(fit->keywords.aperture), NULL),
        KEYWORD( "image", "ISOSPEED", KTYPE_DOUBLE, "ISO camera setting", &(fit->keywords.iso_speed), NULL),
        KEYWORD( "image", "BAYERPAT", KTYPE_STR, "Bayer color pattern", &(fit->keywords.bayer_pattern), NULL),
        KEYWORD( "image", "XBAYROFF", KTYPE_INT, "X offset of Bayer array", &(fit->keywords.bayer_xoffset), NULL),
        KEYWORD( "image", "YBAYROFF", KTYPE_INT, "Y offset of Bayer array", &(fit->keywords.bayer_yoffset), NULL),
        KEYWORD( "image", "GAIN", KTYPE_USHORT, "Camera gain", &(fit->keywords.key_gain), NULL),
        KEYWORD( "image", "OFFSET;BLKLEVEL", KTYPE_USHORT, "Camera offset", &(fit->keywords.key_offset), NULL),
        KEYWORD( "image", "CVF;EGAIN", KTYPE_DOUBLE, "Conversion factor (e-/adu)", &(fit->keywords.cvf), NULL),
        KEYWORD( "image", "AIRMASS", KTYPE_DOUBLE, "Airmass", &(fit->keywords.airmass), NULL),
        KEYWORD( "image", "SITELAT;SITE-LAT;OBSLAT", KTYPE_DOUBLE, "[deg] Observation site latitude", &(fit->keywords.sitelat), sitelat_handler),
        KEYWORD( "image", "SITELONG;SITE-LONOBSLONG", KTYPE_DOUBLE, "[deg] Observation site longitude", &(fit->keywords.sitelong), sitelong_handler),
        KEYWORD( "image", "SITEELEV", KTYPE_DOUBLE, "[m] Observation site elevation", &(fit->keywords.siteelev), NULL),
        KEYWORD( "dft",   "DFTTYPE", KTYPE_STR, "Module/Phase of a Discrete Fourier Transform", &(fit->keywords.dft.type), NULL),
        KEYWORD( "dft",   "DFTORD", KTYPE_STR, "Low/High spatial freq. are located at image center", &(fit->keywords.dft.ord), NULL),
        KEYWORD( "dft",   "DFTNORM1", KTYPE_DOUBLE, "Normalisation value for channel #1", &(fit->keywords.dft.norm[0]), NULL),
        KEYWORD( "dft",   "DFTNORM2", KTYPE_DOUBLE, "Normalisation value for channel #2", &(fit->keywords.dft.norm[1]), NULL),
        KEYWORD( "dft",   "DFTNORM3", KTYPE_DOUBLE, "Normalisation value for channel #3", &(fit->keywords.dft.norm[2]), NULL),

        KEYWORD_FIXED( "wcsdata", "CTYPE3", KTYPE_STR, "RGB image", "RGB", NULL),
        KEYWORD( "wcsdata", "OBJCTRA", KTYPE_STR, "Image center Right Ascension (hms)", &(fit->keywords.wcsdata.objctra), NULL),
        KEYWORD( "wcsdata", "OBJCTDEC", KTYPE_STR, "Image center Declination (dms)", &(fit->keywords.wcsdata.objctdec), NULL),
        KEYWORD( "wcsdata", "RA", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra), NULL),
        KEYWORD( "wcsdata", "RA_D", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra), NULL),
        KEYWORD( "wcsdata", "DEC", KTYPE_DOUBLE, "Image center Declination (deg)", &(fit->keywords.wcsdata.dec), NULL),
        KEYWORD( "wcsdata", "DEC_D", KTYPE_DOUBLE, "Image center Declination (deg)", &(fit->keywords.wcsdata.dec), NULL),
        KEYWORD( "wcsdata", "PLTSOLVD", KTYPE_BOOL, NULL, &(fit->keywords.wcsdata.pltsolvd), pltsolvd_handler),

		/* This group must be the last one !!
		 * It is not used. We write keywords just so that Siril knows about them
		 */
		KEYWORD_WCS( "wcslib", "CTYPE1", KTYPE_STR),
		KEYWORD_WCS( "wcslib", "CTYPE2", KTYPE_STR),
		KEYWORD_WCS( "wcslib", "CUNIT1", KTYPE_STR),
		KEYWORD_WCS( "wcslib", "CUNIT2", KTYPE_STR),
		KEYWORD_WCS( "wcslib", "EQUINOX", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CRPIX1", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CRPIX2", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CRVAL1", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CRVAL2", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "LONPOLE", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CDELT1", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CDELT2", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "PC1_1", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "PC1_2", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "PC2_1", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "PC2_2", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CD1_1", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CD1_2", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CD2_1", KTYPE_DOUBLE),
		KEYWORD_WCS( "wcslib", "CD2_2", KTYPE_DOUBLE),

		{NULL, NULL, KTYPE_BOOL, NULL, NULL, FALSE, TRUE}
    };

	// Count the number of keywords in the list
	int num_keywords = 0;
	while (keyword_list[num_keywords].group != NULL) {
		num_keywords++;
	}

    // Allocate memory dynamically for the keyword array
	KeywordInfo *all_keywords = (KeywordInfo*) malloc((num_keywords + 1) * sizeof(KeywordInfo));

    // Copy keyword information from the list to the dynamic array and set if keyword must be used
	for (int i = 0; i < num_keywords; i++) {
		all_keywords[i] = keyword_list[i];
		all_keywords[i].is_used = should_use_keyword(fit, keyword_list[i].group, keyword_list[i].key);

		// Set default values based on keyword type
		if (g_strcmp0(all_keywords[i].group, "wcslib")) { // This group is initialized somewhere
			switch (all_keywords[i].type) {
			case KTYPE_INT:
				if (*((int*) all_keywords[i].data) == 0)
					*((int*) all_keywords[i].data) = DEFAULT_INT_VALUE;
				break;
			case KTYPE_UINT:
				if (*((guint*)all_keywords[i].data) == 0)
					*((guint*) all_keywords[i].data) = DEFAULT_UINT_VALUE;
				break;
			case KTYPE_USHORT:
				if (*((gushort*)all_keywords[i].data) == 0)
					*((gushort*) all_keywords[i].data) = DEFAULT_USHORT_VALUE;
				break;
			case KTYPE_DOUBLE:
				if (*((double*)all_keywords[i].data) == 0)
					*((double*) all_keywords[i].data) = DEFAULT_DOUBLE_VALUE;
				break;
			default:
				break;
			}
		}
	}

    // Mark the end of the list
	all_keywords[num_keywords] = keyword_list[num_keywords];

	return all_keywords;
}

int save_fits_keywords(fits *fit) {
	KeywordInfo *keys = initialize_keywords(fit);
	KeywordInfo *keys_start = keys;
	int status;
	gchar *str;
	gushort us;
	guint ui;
	int ii;
	double dbl, zero, scale;
	GDateTime *date;

	/* Let's start by most important keywords */
	switch (fit->bitpix) {
	case BYTE_IMG:
	case SHORT_IMG:
		zero = 0.0;
		scale = 1.0;
		break;
	case FLOAT_IMG:
		zero = 0.0;
		scale = 1.0;
		break;
	default:
	case USHORT_IMG:
		zero = 32768.0;
		scale = 1.0;
		break;
	}
	status = 0;
	fits_update_key(fit->fptr, TDOUBLE, "BZERO", &zero, "Offset data range to that of unsigned short", &status);

	status = 0;
	fits_update_key(fit->fptr, TDOUBLE, "BSCALE", &scale, "Default scaling factor",	&status);

	strncpy(fit->keywords.program, "Siril "PACKAGE_VERSION, FLEN_VALUE - 1);

	/* Let's save all other keywords */
	while (keys->group) {
		if (!keys->is_used || g_strcmp0(keys->group, "wcslib") == 0) {
			keys++;
			continue;
		}
		gchar** tokens = g_strsplit(keys->key, ";", -1);
		switch (keys->type) {
			case KTYPE_INT:
				status = 0;
				ii = (*((int*)keys->data));
				if (ii) {
					fits_update_key(fit->fptr, TINT, tokens[0], &ii, keys->comment, &status);
				}
				break;
			case KTYPE_UINT:
				status = 0;
				ui = (*((guint*)keys->data));
				if (ui) {
					fits_update_key(fit->fptr, TUINT, tokens[0], &ui, keys->comment, &status);
				}
				break;
			case KTYPE_USHORT:
				status = 0;
				us = (*((int*)keys->data));
				if (us) {
					fits_update_key(fit->fptr, TUSHORT, tokens[0], &us, keys->comment, &status);
				}
				break;
			case KTYPE_DOUBLE:
				status = 0;
				dbl = *((double*)keys->data);
				if (dbl > DEFAULT_DOUBLE_VALUE) {
					fits_update_key(fit->fptr, TDOUBLE, tokens[0], &dbl, keys->comment, &status);
				}
				break;
			case KTYPE_STR:
				status = 0;
				str = ((gchar*)keys->data);
				if (str && str[0] != '\0') {
					fits_update_key(fit->fptr, TSTRING, tokens[0], str, keys->comment, &status);
				}
				break;
			case KTYPE_DATE:
				status = 0;
				date = *((GDateTime**)keys->data);
				if (date) {
					if (!g_strcmp0("DATE", tokens[0])) {
						int itmp;
						char fit_date[40];
						fits_get_system_time(fit_date, &itmp, &status);
						fits_update_key(fit->fptr, TSTRING, tokens[0], fit_date, keys->comment, &status);
					} else {
						gchar *formatted_date = date_time_to_FITS_date(date);
						fits_update_key(fit->fptr, TSTRING, tokens[0], formatted_date, keys->comment, &status);
						g_free(formatted_date);
					}
				}
				break;
			case KTYPE_BOOL:
				status = 0;
				fits_update_key(fit->fptr, TLOGICAL, tokens[0], &(*((gboolean*)keys->data)), keys->comment, &status);
				break;
			default:
				siril_debug_print("Save_fits_keywords: Error. Type is not handled.\n");
		}
		keys++;
	}

	free(keys_start);

	return 0;
}

int save_fits_unknown_keywords(fits *fit) {
	int status = 0;
	/*** Save list of unknown keys ***/
	if (fit->unknown_keys) {
		for (int i = 0; i < 2; i++) {
			status = 0;
			fits_write_comment(fit->fptr, "************************************************************", &status);
		}
		status = 0;
		fits_write_comment(fit->fptr, "**********************Unknown keywords**********************", &status);
		for (int i = 0; i < 2; i++) {
			status = 0;
			fits_write_comment(fit->fptr, "************************************************************", &status);
		}
		status = associate_header_to_memfile(fit->unknown_keys, fit->fptr);
	}
	return status;
}

int save_history_keywords(fits *fit) {
	int status = 0;

	if (fit->history) {
		GSList *list;
		for (list = fit->history; list; list = list->next) {
			fits_write_history(fit->fptr, (char *)list->data, &status);
		}
	}

	status = 0;
	if (com.history) {
		for (int i = 0; i < com.hist_display; i++) {
			if (com.history[i].history[0] != '\0')
				fits_write_history(fit->fptr, com.history[i].history, &status);
		}
	}

	return status;
}

/* FIXME: DATE-OBS should be used in the new structure, in the handler */
static void read_fits_date_obs_header(fits *fit) {
	int status = 0;
	char ut_start[FLEN_VALUE] = { 0 };
	char date_obs[FLEN_VALUE] = { 0 };

	fits_read_key(fit->fptr, TSTRING, "DATE-OBS", &date_obs, NULL, &status);

	/* In some cases, date is divided in two:
	 * - DATE-OBS
	 * - TIME-OBS
	 * We need to check if we find the "T" inside DATE-OBS.
	 * If not, then try to check for TIME-OBS to get the time
	 */
	if (!g_strstr_len(date_obs, -1, "T")) {
		status = 0;
		char time_obs[FLEN_VALUE] = { 0 };
		fits_read_key(fit->fptr, TSTRING, "TIME-OBS", &time_obs, NULL, &status);
		if (!status) {
			strcat(date_obs, "T");
			strcat(date_obs, time_obs);
		}
	}

	/** Case seen in some FITS files. Needed to get date back in SER conversion **/
	status = 0;
	fits_read_key(fit->fptr, TSTRING, "UT-START", &ut_start, NULL, &status);
	if (status == 0 && ut_start[0] != '\0' && date_obs[2] == '/') {
		int year, month, day;
		if (sscanf(date_obs, "%02d/%02d/%04d", &day, &month, &year) == 3) {
			g_snprintf(date_obs, sizeof(date_obs), "%04d-%02d-%02dT%s", year, month, day, ut_start);
		}
	}

	fit->keywords.date_obs = FITS_date_to_date_time(date_obs);
}

int read_fits_keywords(fits *fit) {
	KeywordInfo *keys = initialize_keywords(fit);
	KeywordInfo *keys_start = keys;
	int status = 0;
	int key_number = 1;
	GString *UnknownKeys = g_string_new(NULL);

	/***** Special cases *****/

	double scale = 0.0, zero = 0.0;
	fits_read_key(fit->fptr, TDOUBLE, "BSCALE", &scale, NULL, &status);
	if (!status && 1.0 != scale) {
		siril_log_message(_("Loaded FITS file has a BSCALE different than 1 (%f)\n"), scale);
		status = 0;
		/* We reset the scaling factors as we don't use it */
		fits_set_bscale(fit->fptr, 1.0, 0.0, &status);
	}

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "BZERO", &zero, NULL, &status);
	if (!status && 0.0 != zero && fit->bitpix == FLOAT_IMG) {
		fprintf(stdout, "ignoring BZERO\n");
		fits_set_bscale(fit->fptr, 1.0, 0.0, &status);
	}

//	read_fits_date_obs_header(fit);

	/***** End of special cases *****/

	fits_get_hdrspace(fit->fptr, &key_number, NULL, &status); /* get # of keywords */

	for (int ii = 1; ii <= key_number; ii++) {
		char card[FLEN_CARD];
		status = 0;
		if (fits_read_record(fit->fptr, ii, card, &status)) {
			fits_report_error(stderr, status);
			break;
		}
		char keyname[FLEN_KEYWORD];
		char value[FLEN_VALUE] = { 0 };
		char comment[FLEN_COMMENT];
		int length = 0;
		char type;

		fits_get_keyname(card, keyname, &length, &status);
		fits_parse_value(card, value, comment, &status);
		fits_get_keytype(value, &type, &status);
		status = 0;

		/* FIXME: need to handle MIPS, ... */

	    /* These have been already processed */
	    if (g_strcmp0(keyname, "BSCALE") == 0 || g_strcmp0(keyname, "BZERO") == 0) {
	        continue;
	    }

	    KeywordInfo *current_key = keys_start;
		gboolean value_set = FALSE; // Flag to indicate if a value was set
		while (current_key->group && !value_set) {
			if (fits_get_keyclass(card) == TYP_STRUC_KEY) {
				value_set = TRUE;
				current_key++;
				continue;
			}
			gchar** tokens = g_strsplit(current_key->key, ";", -1);
			int n = g_strv_length(tokens);
			for (int i = 0; i < n && !value_set; i++) {
				if (g_strcmp0(tokens[i], keyname) == 0) {
					if (g_strcmp0("EQUINOX", keyname) == 0) {
						printf("stop\n");
					}
					if (current_key->fixed_value) {
						value_set = TRUE;
						break;
					}
					int int_value;
					guint uint_value;
					gushort ushort_value;
					double double_value;
					gchar *str_value;
					GDateTime *date;
					gboolean bool;
					char *end;

					switch (current_key->type) {
					case KTYPE_INT:
						int_value = g_ascii_strtoll(value, &end, 10);
						if (value != end) {
							*((int*) current_key->data) = int_value;
							value_set = TRUE;
						}
						break;
					case KTYPE_UINT:
						uint_value = g_ascii_strtoll(value, &end, 10);
						if (value != end) {
							*((guint*) current_key->data) = uint_value;
							value_set = TRUE;
						}
						break;
					case KTYPE_USHORT:
						ushort_value = g_ascii_strtoll(value, &end, 10);
						if (value != end) {
							*((gushort*) current_key->data) = ushort_value;
							value_set = TRUE;
						}
						break;
					case KTYPE_DOUBLE:
						double_value = g_ascii_strtod(value, &end);
						if (value != end) {
							*((double*) current_key->data) = double_value;
							value_set = TRUE;
						}
						break;
					case KTYPE_STR:
						str_value = g_shell_unquote(value, NULL);
						strncpy((char*) current_key->data, str_value, FLEN_VALUE - 1);
						value_set = TRUE;
						break;
					case KTYPE_DATE:
						str_value = g_shell_unquote(value, NULL);
						date = FITS_date_to_date_time(str_value);
						if (date) {
							*((GDateTime**) current_key->data) = date;
							value_set = TRUE;
						}
						break;
					case KTYPE_BOOL:
						bool = value[0] == 'T' ? TRUE : FALSE;
						*((gboolean*) current_key->data) = bool;
						value_set = TRUE;
						break;
					default:
						break;
					}
					/* Handle special cases */
					if (current_key->special_handler) {
						current_key->special_handler(fit, comment, current_key);
					}
				}
			}

			g_strfreev(tokens);
			current_key++;
		}
		/* output unknown keywords */
		/** FIXME: need to exclude all keywords not handled but known. Like WCS and SIP keywords */
		if (!value_set) {
			UnknownKeys = g_string_append(UnknownKeys, card);
			UnknownKeys = g_string_append(UnknownKeys, "\n");
		}
	}
	free(keys_start);

	/** Retrieve unknown keys **/
	if (fit->unknown_keys)
		g_free(fit->unknown_keys);
	fit->unknown_keys = g_string_free(UnknownKeys, FALSE);

	return 0;
}



