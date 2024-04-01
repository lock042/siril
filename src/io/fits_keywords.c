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

#include "fits_keywords.h"

#define DBL_FLAG -999.0

#define KEYWORD(group, key, type, comment, data) { group, key, type, comment, data, TRUE }

static gboolean should_use_keyword(const fits *fit, const gchar *keyword) {
    if (g_strcmp0(keyword, "XBAYROFF") == 0) {
        return fit->keywords.bayer_pattern[0] != '\0';
    } else if (g_strcmp0(keyword, "YBAYROFF") == 0) {
        return fit->keywords.bayer_pattern[0] != '\0';
    } else if (g_strcmp0(keyword, "DFTNORM2") == 0) {
        return fit->naxes[2] > 1;
    } else if (g_strcmp0(keyword, "DFTNORM3") == 0) {
        return fit->naxes[2] > 1;
    } else if (g_strcmp0(keyword, "CTYPE3") == 0) {
        return (fit->naxes[2] > 1  && com.pref.rgb_aladin);
    }

    return TRUE;
}

KeywordInfo *initialize_keywords(fits *fit) {
	KeywordInfo keyword_list[] = {
//        KEYWORD( "image", "BZERO", KTYPE_DOUBLE, "Offset data range to that of unsigned short", &(fit->keywords.bzero)),
//        KEYWORD( "image", "BSCALE", KTYPE_DOUBLE, "Default scaling factor", &(fit->keywords.bscale)),
        KEYWORD( "image", "ROWORDER", KTYPE_STR, "Order of the rows in image array", &(fit->keywords.row_order)),
        KEYWORD( "setup", "INSTRUME", KTYPE_STR, "Instrument name", &(fit->keywords.instrume)),
        KEYWORD( "setup", "TELESCOP", KTYPE_STR, "Telescope used to acquire this image", &(fit->keywords.telescop)),
        KEYWORD( "setup", "OBSERVER", KTYPE_STR, "Observer name", &(fit->keywords.observer)),
        KEYWORD( "date",  "DATE", KTYPE_DATE, "UTC date that FITS file was created", &(fit->keywords.date)),
        KEYWORD( "date",  "DATE-OBS", KTYPE_DATE, "YYYY-MM-DDThh:mm:ss observation start, UT", &(fit->keywords.date_obs)),
        KEYWORD( "image", "STACKCNT;NCOMBINE", KTYPE_UINT, "Stack frames", &(fit->keywords.stackcnt)),
        KEYWORD( "image", "EXPTIME;EXPOSURE", KTYPE_DOUBLE, "Exposure time [s]", &(fit->keywords.exposure)),
        KEYWORD( "image", "LIVETIME", KTYPE_DOUBLE, "Exposure time after deadtime correction", &(fit->keywords.livetime)),
        KEYWORD( "image", "EXPSTART", KTYPE_DOUBLE, "Exposure start time (standard Julian date)", &(fit->keywords.expstart)),
        KEYWORD( "image", "EXPEND", KTYPE_DOUBLE, "Exposure end time (standard Julian date)", &(fit->keywords.expend)),
        KEYWORD( "image", "XPIXSZ;XPIXELSZ;PIXSIZE1;PIXSIZEX;XPIXSIZE", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x)),
        KEYWORD( "image", "YPIXSZ;YPIXELSZ;PIXSIZE2;PIXSIZEY;YPIXSIZE", KTYPE_DOUBLE, "Y pixel size microns", &(fit->keywords.pixel_size_y)),
        KEYWORD( "image", "XBINNING;BINX", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_x)),
        KEYWORD( "image", "YBINNING;BINY", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_y)),
        KEYWORD( "image", "FOCALLEN;FOCAL", KTYPE_DOUBLE, "Camera focal length", &(fit->keywords.focal_length)),
        KEYWORD( "image", "CCD-TEMP;CCD_TEMP;CCDTEMP;TEMPERAT;CAMTCCD", KTYPE_DOUBLE, "CCD temp in C", &(fit->keywords.ccd_temp)),
        KEYWORD( "image", "SET-TEMP", KTYPE_DOUBLE, "Temperature setting in C", &(fit->keywords.set_temp)),
        KEYWORD( "image", "FILTER;FILT-1", KTYPE_STR, "Active filter name", &(fit->keywords.filter)),
        KEYWORD( "image", "IMAGETYP;FRAMETYP", KTYPE_STR, "Type of image", &(fit->keywords.image_type)),
        KEYWORD( "image", "OBJECT", KTYPE_STR, "Name of the object of interest", &(fit->keywords.object)),
        KEYWORD( "image", "APERTURE", KTYPE_DOUBLE, "Aperture of the instrument", &(fit->keywords.aperture)),
        KEYWORD( "image", "ISOSPEED", KTYPE_DOUBLE, "ISO camera setting", &(fit->keywords.iso_speed)),
        KEYWORD( "image", "BAYERPAT", KTYPE_STR, "Bayer color pattern", &(fit->keywords.bayer_pattern)),
        KEYWORD( "image", "XBAYROFF", KTYPE_INT, "X offset of Bayer array", &(fit->keywords.bayer_xoffset)),
        KEYWORD( "image", "YBAYROFF", KTYPE_INT, "Y offset of Bayer array", &(fit->keywords.bayer_yoffset)),
        KEYWORD( "image", "GAIN", KTYPE_USHORT, "Camera gain", &(fit->keywords.key_gain)),
        KEYWORD( "image", "OFFSET;BLKLEVEL", KTYPE_USHORT, "Camera offset", &(fit->keywords.key_offset)),
        KEYWORD( "image", "CVF;EGAIN", KTYPE_DOUBLE, "Conversion factor (e-/adu)", &(fit->keywords.cvf)),
        KEYWORD( "image", "AIRMASS", KTYPE_DOUBLE, "Airmass", &(fit->keywords.airmass)),
        KEYWORD( "image", "SITELAT;SITE-LAT;OBSLAT", KTYPE_DOUBLE, "[deg] Observation site latitude", &(fit->keywords.sitelat)),
        KEYWORD( "image", "SITELONG;SITE-LONOBSLONG", KTYPE_DOUBLE, "[deg] Observation site longitude", &(fit->keywords.sitelong)),
        KEYWORD( "image", "SITEELEV", KTYPE_DOUBLE, "[m] Observation site elevation", &(fit->keywords.siteelev)),
        KEYWORD( "dft",   "DFTTYPE", KTYPE_STR, "Module/Phase of a Discrete Fourier Transform", &(fit->keywords.dft.type)),
        KEYWORD( "dft",   "DFTORD", KTYPE_STR, "Low/High spatial freq. are located at image center", &(fit->keywords.dft.ord)),
        KEYWORD( "dft",   "DFTNORM1", KTYPE_DOUBLE, "Normalisation value for channel #1", &(fit->keywords.dft.norm[0])),
        KEYWORD( "dft",   "DFTNORM2", KTYPE_DOUBLE, "Normalisation value for channel #2", &(fit->keywords.dft.norm[1])),
        KEYWORD( "dft",   "DFTNORM3", KTYPE_DOUBLE, "Normalisation value for channel #3", &(fit->keywords.dft.norm[2])),
        KEYWORD( "image", "PROGRAMM", KTYPE_STR, "Software that created this HDU", "Siril "PACKAGE_VERSION),

        KEYWORD( "wcs",   "CTYPE3", KTYPE_STR, "RGB image", "RGB"),
        KEYWORD( "wcs",   "OBJCTRA", KTYPE_STR, "Image center Right Ascension (hms)", &(fit->keywords.wcsdata.objctra)),
        KEYWORD( "wcs",   "OBJCTDEC", KTYPE_STR, "Image center Declination (dms)", &(fit->keywords.wcsdata.objctdec)),
        KEYWORD( "wcs",   "RA", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra)),
        KEYWORD( "wcs",   "DEC", KTYPE_DOUBLE, "Image center Declination (deg)", &(fit->keywords.wcsdata.dec)),

		{NULL, NULL, KTYPE_BOOL, NULL, NULL, FALSE}
    };

	// Count the number of keywords in the list
	int num_keywords = 0;
	while (keyword_list[num_keywords].group != NULL) {
		num_keywords++;
	}

    // Allocate memory dynamically for the keyword array
	KeywordInfo *all_keywords = (KeywordInfo*) malloc((num_keywords + 1) * sizeof(KeywordInfo));

    // Copy keyword information from the list to the dynamic array
	for (int i = 0; i < num_keywords; i++) {
		all_keywords[i] = keyword_list[i];
		all_keywords[i].is_used = should_use_keyword(fit, keyword_list[i].key);
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
	double dbl;
	GDateTime *date;
	while (keys->group) {
		if (!keys->is_used) {
			keys++;
			continue;
		}
		gchar** tokens = g_strsplit(keys->key, ";", -1);
		switch (keys->type) {
			case KTYPE_INT:
				status = 0;
				fits_update_key(fit->fptr, TINT, tokens[0], &(*((int*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_UINT:
				status = 0;
				fits_update_key(fit->fptr, TUINT, tokens[0], &(*((guint*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_USHORT:
				status = 0;
				fits_update_key(fit->fptr, TUSHORT, tokens[0], &(*((int*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_DOUBLE:
				status = 0;
				dbl = *((double*)keys->data);
				if (dbl > DBL_FLAG) {
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
			default:
				siril_debug_print("Save_fits_keywords: Error. Type is not handled.\n");
		}
		keys++;
	}

	free(keys_start);

	/*******************************************************************
	 * ********************* HISTORY KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
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

	return 0;
}
