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
#include "algos/siril_wcs.h"


#include "fits_keywords.h"

#define DBL_FLAG -999.0

#define KEYWORD(group, key, type, comment, data) { group, key, type, comment, data, TRUE, FALSE }
#define KEYWORD_FIXED(group, key, type, comment, data) { group, key, type, comment, data, TRUE, TRUE }

static gboolean should_use_keyword(const fits *fit, const gchar *group, const gchar *keyword) {
	gboolean use_keyword = TRUE;
	if (g_strcmp0(group, "wcslib") == 0) {
		use_keyword = (fit->keywords.wcslib != NULL);
	}

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
    } else if (g_strcmp0(keyword, "CDELT1") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_1);
    } else if (g_strcmp0(keyword, "CDELT2") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_1);
    } else if (g_strcmp0(keyword, "PC1_1") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_1);
    } else if (g_strcmp0(keyword, "PC1_2") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_1);
    } else if (g_strcmp0(keyword, "PC2_1") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_1);
    } else if (g_strcmp0(keyword, "PC2_2") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_1);
    } else if (g_strcmp0(keyword, "CD1_1") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_2);
    } else if (g_strcmp0(keyword, "CD1_2") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_2);
    } else if (g_strcmp0(keyword, "CD2_1") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_2);
    } else if (g_strcmp0(keyword, "CD2_2") == 0) {
        return (use_keyword && com.pref.wcs_formalism == WCS_FORMALISM_2);
    }

    return use_keyword;
}

KeywordInfo *initialize_keywords(fits *fit) {
	KeywordInfo keyword_list[] = {
			// FIXME: add MIPS keywords
	    KEYWORD( "image", "MIPS-HI", KTYPE_USHORT, "Upper visualization cutoff", &(fit->keywords.hi)),
	    KEYWORD( "image", "MIPS-LO", KTYPE_USHORT, "Lower visualization cutoff", &(fit->keywords.lo)),
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
		KEYWORD_FIXED( "image", "PROGRAMM", KTYPE_STR, "Software that created this HDU", "Siril "PACKAGE_VERSION),

		KEYWORD_FIXED( "wcsdata",   "CTYPE3", KTYPE_STR, "RGB image", "RGB"),
        KEYWORD( "wcsdata",   "OBJCTRA", KTYPE_STR, "Image center Right Ascension (hms)", &(fit->keywords.wcsdata.objctra)),
        KEYWORD( "wcsdata",   "OBJCTDEC", KTYPE_STR, "Image center Declination (dms)", &(fit->keywords.wcsdata.objctdec)),
        KEYWORD( "wcsdata",   "RA", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra)),
        KEYWORD( "wcsdata",   "DEC", KTYPE_DOUBLE, "Image center Declination (deg)", &(fit->keywords.wcsdata.dec)),
//        KEYWORD( "wcslib",   "CTYPE1", KTYPE_STR, "TAN (gnomic) projection", "RA---TAN"), // FIXME: handle both version of comments
//        KEYWORD( "wcslib",   "CTYPE2", KTYPE_STR, "TAN (gnomic) projection", "DEC---TAN"), // FIXME: handle both version of comments
		KEYWORD_FIXED( "wcslib", "CUNIT1", KTYPE_STR, "Unit of coordinates", "deg"),
		KEYWORD_FIXED( "wcslib", "CUNIT1", KTYPE_STR, "Unit of coordinates", "deg"),
        KEYWORD( "wcslib", "EQUINOX", KTYPE_DOUBLE, "Equatorial equinox", &(fit->keywords.wcslib->equinox)),
        KEYWORD( "wcslib", "CRPIX1", KTYPE_DOUBLE, "Axis1 reference pixel", &(fit->keywords.wcslib->crpix[0])),
        KEYWORD( "wcslib", "CRPIX2", KTYPE_DOUBLE, "Axis2 reference pixel", &(fit->keywords.wcslib->crpix[1])),
        KEYWORD( "wcslib", "CRVAL1", KTYPE_DOUBLE, "Axis1 reference value (deg)", &(fit->keywords.wcslib->crval[0])),
        KEYWORD( "wcslib", "CRVAL2", KTYPE_DOUBLE, "Axis2 reference value (deg)", &(fit->keywords.wcslib->crval[1])),
        KEYWORD( "wcslib", "LONPOLE", KTYPE_DOUBLE, "Native longitude of celestial pole", &(fit->keywords.wcslib->lonpole)),
        KEYWORD( "wcslib", "CDELT1", KTYPE_DOUBLE, "X pixel size (deg)", &(fit->keywords.wcslib->cdelt[0])),
        KEYWORD( "wcslib", "CDELT2", KTYPE_DOUBLE, "X pixel size (deg)", &(fit->keywords.wcslib->cdelt[1])),
        KEYWORD( "wcslib", "PC1_1", KTYPE_DOUBLE, "Linear transformation matrix (1, 1)", &(fit->keywords.wcslib->pc[0])),
        KEYWORD( "wcslib", "PC1_2", KTYPE_DOUBLE, "Linear transformation matrix (1, 2)", &(fit->keywords.wcslib->pc[1])),
        KEYWORD( "wcslib", "PC2_1", KTYPE_DOUBLE, "Linear transformation matrix (2, 1)", &(fit->keywords.wcslib->pc[2])),
        KEYWORD( "wcslib", "PC2_2", KTYPE_DOUBLE, "Linear transformation matrix (2, 2)", &(fit->keywords.wcslib->pc[3])),
        KEYWORD( "wcslib", "CD1_1", KTYPE_DOUBLE, "Scale matrix (1, 1)", &(fit->keywords.wcslib->cd[0])),
        KEYWORD( "wcslib", "CD1_2", KTYPE_DOUBLE, "Scale matrix (1, 2)", &(fit->keywords.wcslib->cd[1])),
        KEYWORD( "wcslib", "CD2_1", KTYPE_DOUBLE, "Scale matrix (2, 1)", &(fit->keywords.wcslib->cd[2])),
        KEYWORD( "wcslib", "CD2_2", KTYPE_DOUBLE, "Scale matrix (2, 2)", &(fit->keywords.wcslib->cd[3])),

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
	ushort us;
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

	/* Let's save all other keywords */
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
				us = (*((int*)keys->data));
				if (us) {
					fits_update_key(fit->fptr, TUSHORT, tokens[0], &us, keys->comment, &status);
				}
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

int read_fits_keywords(fits *fit) {
	KeywordInfo *keys = initialize_keywords(fit);
	KeywordInfo *keys_start = keys;
	int status = 0;
	int key_number = 1;
	while (1) {
		char card[FLEN_CARD];
		if (fits_read_record(fit->fptr, key_number++, card, &status)) {
			break;
		}
		char keyname[FLEN_KEYWORD];
		char value[FLEN_VALUE] = { 0 };
		int length = 0;
		char type;

		fits_get_keyname(card, keyname, &length, &status);
		fits_parse_value(card, value, NULL, &status);
		fits_get_keytype(value, &type, &status);

		/* FIXME: need to handle MIPS, BZERO, BSCALE, ... */
//		status = 0;
//		fits_read_key(fit->fptr, TDOUBLE, "BSCALE", &scale, NULL, &status);
//		if (!status && 1.0 != scale) {
//			siril_log_message(_("Loaded FITS file has a BSCALE different than 1 (%f)\n"), scale);
//			status = 0;
//			/* We reset the scaling factors as we don't use it */
//			fits_set_bscale(fit->fptr, 1.0, 0.0, &status);
//		}
//
//		status = 0;
//		fits_read_key(fit->fptr, TDOUBLE, "BZERO", &zero, NULL, &status);
//		if (!status && 0.0 != zero && fit->bitpix == FLOAT_IMG) {
//			fprintf(stdout, "ignoring BZERO\n");
//			fits_set_bscale(fit->fptr, 1.0, 0.0, &status);
//		}

		while (keys->group) {
			if (keys->fixed_value) {
				/* if fixed value, we do not read it, we do not store it */
				keys++;
				continue;
			}
			gchar** tokens = g_strsplit(keys->key, ";", -1);
			int n = g_strv_length(tokens);
			for (int i = 0; i < n; i++) {
				if (g_strcmp0(tokens[i], keyname) == 0) {
					int int_value = 0;
					guint uint_value = 0;
					gushort ushort_value = 0;
					double double_value = DBL_FLAG;

					switch (keys->type) {
					case KTYPE_INT:
						sscanf(value, "%d", &int_value);
						*((int*) keys->data) = int_value;
						break;
					case KTYPE_UINT:
						sscanf(value, "%u", &uint_value);
						*((guint*) keys->data) = uint_value;
						break;
					case KTYPE_USHORT:
						sscanf(value, "%hu", &ushort_value);
						*((gushort*) keys->data) = ushort_value;
						break;
					case KTYPE_DOUBLE:
						sscanf(value, "%lf", &double_value);
						*((double*) keys->data) = double_value;
						break;
					case KTYPE_STR:
						strcpy((char*) keys->data, value);
						break;
					case KTYPE_DATE:
						status = 0;
						// FIXME: convert to GDateTime
						break;
					default:
						break;
					}
					break;
				}
			}
			keys++;
		}
		free(keys_start);
	}
	return 0;
}



