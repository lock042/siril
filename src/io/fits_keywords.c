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
#define DEFAULT_FLOAT_VALUE -999.f
#define DEFAULT_INT_VALUE -INT_MAX
#define DEFAULT_UINT_VALUE 0
#define DEFAULT_USHORT_VALUE DEFAULT_UINT_VALUE

#define KEYWORD_PRIMARY(group, key, type, comment, data, handler_read, handler_save) { group, key, type, comment, data, handler_read, handler_save, TRUE, FALSE }
#define KEYWORD_SECONDA(group, key, type, comment, data, handler_read, handler_save) { group, key, type, comment, data, handler_read, handler_save, FALSE, FALSE }
#define KEYWORD_FIXED(group, key, type, comment, data, handler_read, handler_save) { group, key, type, comment, data, handler_read, handler_save, TRUE, TRUE }
#define KEYWORD_WCS(group, key, type) { group, key, type, NULL, NULL, NULL, NULL, FALSE, TRUE }
#define UNKNOWN_KEYWORDS "**********************Unknown keywords**********************"

static gboolean should_use_keyword(const fits *fit, KeywordInfo keyword) {

	if (g_strcmp0(keyword.key, "MIPS-HI") == 0) {
		return (fit->type == DATA_USHORT);
	} else if (g_strcmp0(keyword.key, "MIPS-LO") == 0) {
		return (fit->type == DATA_USHORT);
	} else if (g_strcmp0(keyword.key, "MIPS-FLO") == 0) {
		return (fit->type == DATA_FLOAT);
	} else if (g_strcmp0(keyword.key, "MIPS-FHI") == 0) {
		return (fit->type == DATA_FLOAT);
	} else if (g_strcmp0(keyword.key, "ROWORDER") == 0) {
		return ((g_strcmp0(fit->keywords.row_order, "BOTTOM-UP") == 0)
				|| (g_strcmp0(fit->keywords.row_order, "TOP-DOWN") == 0));
	} else if (g_strcmp0(keyword.key, "FLENGTH") == 0) {
		return FALSE;
	} else if (g_strcmp0(keyword.key, "XBAYROFF") == 0) {
		return fit->keywords.bayer_pattern[0] != '\0';
	} else if (g_strcmp0(keyword.key, "YBAYROFF") == 0) {
		return fit->keywords.bayer_pattern[0] != '\0';
	} else if (g_strcmp0(keyword.key, "DFTNORM2") == 0) {
		return fit->naxes[2] > 1;
	} else if (g_strcmp0(keyword.key, "DFTNORM3") == 0) {
		return fit->naxes[2] > 1;
	} else if (g_strcmp0(keyword.key, "RA_D") == 0) {
		return FALSE;
	} else if (g_strcmp0(keyword.key, "DEC_D") == 0) {
		return FALSE;
	}
	return keyword.is_used;
}

/****************************** handlers ******************************/

static void pixel_x_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.pixel_size_x > 0.0) {
		fit->pixelkey = TRUE;
	}
}

static void binning_x_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.binning_x <= 0)
		fit->keywords.binning_x = 1;
}

static void binning_y_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.binning_y <= 0)
		fit->keywords.binning_y = 1;
}

static void roworder_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (!strcasecmp(fit->keywords.bayer_pattern, "NONE")) {
		memset(fit->keywords.bayer_pattern, 0, sizeof(char) * FLEN_VALUE);
	}
}

static void focal_length_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (fit->keywords.focal_length > 0.0)
		fit->focalkey = TRUE;
}

static void flength_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	fit->keywords.focal_length = fit->keywords.flength * 1000.0; // convert m to mm
	fit->focalkey = TRUE;
}

static void sitelong_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
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

static void sitelat_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
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

static void datamax_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	gboolean not_from_siril = (strstr(fit->keywords.program, "Siril") == NULL);
	if ((fit->bitpix == FLOAT_IMG && not_from_siril) || fit->bitpix == DOUBLE_IMG) {
		float mini, maxi;
		fit_stats(fit->fptr, &mini, &maxi);
		// override data_max if needed. In some images there are differences between max and data_max
		fit->keywords.data_max = (double) maxi;
		fit->keywords.data_min = (double) mini;
	}
}

static void flo_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (!fit->keywords.hi && (fit->orig_bitpix == FLOAT_IMG || fit->orig_bitpix == DOUBLE_IMG)) {
		fit->keywords.lo = float_to_ushort_range(fit->keywords.flo);
	}
}

static void fhi_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (!fit->keywords.hi && (fit->orig_bitpix == FLOAT_IMG || fit->orig_bitpix == DOUBLE_IMG)) {
		fit->keywords.hi = float_to_ushort_range(fit->keywords.fhi);
	}
}


static void flo_handler_save(fits *fit, KeywordInfo *info) {
	if (!fit->keywords.hi && (fit->orig_bitpix == FLOAT_IMG || fit->orig_bitpix == DOUBLE_IMG)) {
		fit->keywords.flo = ushort_to_float_range(fit->keywords.lo);
	}
}

static void fhi_handler_save(fits *fit, KeywordInfo *info) {
	if (!fit->keywords.hi && (fit->orig_bitpix == FLOAT_IMG || fit->orig_bitpix == DOUBLE_IMG)) {
		fit->keywords.fhi = ushort_to_float_range(fit->keywords.hi);
	}
	if (fit->keywords.fhi == 0.0) {
		fit->keywords.fhi = DEFAULT_FLOAT_VALUE;
		fit->keywords.flo = DEFAULT_FLOAT_VALUE;
	}
}



/*****************************************************************************/

KeywordInfo *initialize_keywords(fits *fit, GHashTable **hash) {
	KeywordInfo keyword_list[] = {
        KEYWORD_PRIMARY( "image", "MIPS-HI", KTYPE_USHORT, "Lower visualization cutoff", &(fit->keywords.hi), NULL, NULL),
        KEYWORD_SECONDA( "image", "CWHITE", KTYPE_USHORT, "Lower visualization cutoff", &(fit->keywords.hi), NULL, NULL),
        KEYWORD_PRIMARY( "image", "MIPS-LO", KTYPE_USHORT, "Upper visualization cutoff", &(fit->keywords.lo), NULL, NULL),
        KEYWORD_SECONDA( "image", "CBLACK", KTYPE_USHORT, "Upper visualization cutoff", &(fit->keywords.lo), NULL, NULL),
        KEYWORD_PRIMARY( "image", "MIPS-FLO", KTYPE_FLOAT, "Lower visualization cutoff", &(fit->keywords.flo), flo_handler_read, flo_handler_save),
        KEYWORD_PRIMARY( "image", "MIPS-FHI", KTYPE_FLOAT, "Upper visualization cutoff", &(fit->keywords.fhi), fhi_handler_read, fhi_handler_save),
        /* ATTENTION: PROGRAM MUST BE BEFORE DATAMAX */
        KEYWORD_PRIMARY( "image", "PROGRAM", KTYPE_STR, "Software that created this HDU", &(fit->keywords.program), NULL, NULL),
        KEYWORD_PRIMARY( "image", "DATAMAX", KTYPE_DOUBLE, "Order of the rows in image array", &(fit->keywords.data_max), datamax_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "ROWORDER", KTYPE_STR, "Order of the rows in image array", &(fit->keywords.row_order), roworder_handler_read, NULL),
        KEYWORD_PRIMARY( "setup", "INSTRUME", KTYPE_STR, "Instrument name", &(fit->keywords.instrume), NULL, NULL),
        KEYWORD_PRIMARY( "setup", "TELESCOP", KTYPE_STR, "Telescope used to acquire this image", &(fit->keywords.telescop), NULL, NULL),
        KEYWORD_PRIMARY( "setup", "OBSERVER", KTYPE_STR, "Observer name", &(fit->keywords.observer), NULL, NULL),
        KEYWORD_PRIMARY( "date",  "DATE", KTYPE_DATE, "UTC date that FITS file was created", &(fit->keywords.date), NULL, NULL),
        KEYWORD_PRIMARY( "date",  "DATE-OBS", KTYPE_DATE, "YYYY-MM-DDThh:mm:ss observation start, UT", &(fit->keywords.date_obs), NULL, NULL),
        KEYWORD_PRIMARY( "image", "STACKCNT", KTYPE_UINT, "Stack frames", &(fit->keywords.stackcnt), NULL, NULL),
        KEYWORD_SECONDA( "image", "NCOMBINE", KTYPE_UINT, "Stack frames", &(fit->keywords.stackcnt), NULL, NULL),
        KEYWORD_PRIMARY( "image", "EXPTIME", KTYPE_DOUBLE, "Exposure time [s]", &(fit->keywords.exposure), NULL, NULL),
        KEYWORD_SECONDA( "image", "EXPOSURE", KTYPE_DOUBLE, "Exposure time [s]", &(fit->keywords.exposure), NULL, NULL),
        KEYWORD_PRIMARY( "image", "LIVETIME", KTYPE_DOUBLE, "Exposure time after deadtime correction", &(fit->keywords.livetime), NULL, NULL),
        KEYWORD_PRIMARY( "image", "EXPSTART", KTYPE_DOUBLE, "Exposure start time (standard Julian date)", &(fit->keywords.expstart), NULL, NULL),
        KEYWORD_PRIMARY( "image", "EXPEND", KTYPE_DOUBLE, "Exposure end time (standard Julian date)", &(fit->keywords.expend), NULL, NULL),
        KEYWORD_PRIMARY( "image", "XPIXSZ", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "XPIXELSZ", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "PIXSIZE1", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "PIXSIZEX", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "XPIXSIZE", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "YPIXSZ", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "YPIXELSZ", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "PIXSIZE2", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "PIXSIZEY", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "YPIXSIZE", KTYPE_DOUBLE, "X pixel size microns", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "XBINNING", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_x), binning_x_handler_read, NULL),
        KEYWORD_SECONDA( "image", "BINX", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_x), binning_x_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "YBINNING", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_y), binning_y_handler_read, NULL),
        KEYWORD_SECONDA( "image", "BINY", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_y), binning_y_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "FOCALLEN", KTYPE_DOUBLE, "Camera focal length", &(fit->keywords.focal_length), focal_length_handler_read, NULL),
        KEYWORD_SECONDA( "image", "FOCAL", KTYPE_DOUBLE, "Camera focal length", &(fit->keywords.focal_length), focal_length_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "FLENGTH", KTYPE_DOUBLE, "Camera focal length [m]", &(fit->keywords.flength), flength_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "CCD-TEMP", KTYPE_DOUBLE, "CCD temp in C", &(fit->keywords.ccd_temp), NULL, NULL),
        KEYWORD_SECONDA( "image", "CCD_TEMP", KTYPE_DOUBLE, "CCD temp in C", &(fit->keywords.ccd_temp), NULL, NULL),
        KEYWORD_SECONDA( "image", "CCDTEMP", KTYPE_DOUBLE, "CCD temp in C", &(fit->keywords.ccd_temp), NULL, NULL),
        KEYWORD_SECONDA( "image", "TEMPERAT", KTYPE_DOUBLE, "CCD temp in C", &(fit->keywords.ccd_temp), NULL, NULL),
        KEYWORD_SECONDA( "image", "CAMTCCD", KTYPE_DOUBLE, "CCD temp in C", &(fit->keywords.ccd_temp), NULL, NULL),
        KEYWORD_PRIMARY( "image", "SET-TEMP", KTYPE_DOUBLE, "Temperature setting in C", &(fit->keywords.set_temp), NULL, NULL),
        KEYWORD_PRIMARY( "image", "FILTER", KTYPE_STR, "Active filter name", &(fit->keywords.filter), NULL, NULL),
        KEYWORD_SECONDA( "image", "FILT-1", KTYPE_STR, "Active filter name", &(fit->keywords.filter), NULL, NULL),
        KEYWORD_PRIMARY( "image", "IMAGETYP", KTYPE_STR, "Type of image", &(fit->keywords.image_type), NULL, NULL),
        KEYWORD_SECONDA( "image", "FRAMETYP", KTYPE_STR, "Type of image", &(fit->keywords.image_type), NULL, NULL),
        KEYWORD_PRIMARY( "image", "OBJECT", KTYPE_STR, "Name of the object of interest", &(fit->keywords.object), NULL, NULL),
        KEYWORD_PRIMARY( "image", "APERTURE", KTYPE_DOUBLE, "Aperture of the instrument", &(fit->keywords.aperture), NULL, NULL),
        KEYWORD_PRIMARY( "image", "ISOSPEED", KTYPE_DOUBLE, "ISO camera setting", &(fit->keywords.iso_speed), NULL, NULL),
        KEYWORD_PRIMARY( "image", "BAYERPAT", KTYPE_STR, "Bayer color pattern", &(fit->keywords.bayer_pattern), NULL, NULL),
        KEYWORD_PRIMARY( "image", "XBAYROFF", KTYPE_INT, "X offset of Bayer array", &(fit->keywords.bayer_xoffset), NULL, NULL),
        KEYWORD_PRIMARY( "image", "YBAYROFF", KTYPE_INT, "Y offset of Bayer array", &(fit->keywords.bayer_yoffset), NULL, NULL),
        KEYWORD_PRIMARY( "image", "GAIN", KTYPE_USHORT, "Camera gain", &(fit->keywords.key_gain), NULL, NULL),
        KEYWORD_PRIMARY( "image", "OFFSET", KTYPE_USHORT, "Camera offset", &(fit->keywords.key_offset), NULL, NULL),
        KEYWORD_SECONDA( "image", "BLKLEVEL", KTYPE_USHORT, "Camera offset", &(fit->keywords.key_offset), NULL, NULL),
        KEYWORD_PRIMARY( "image", "CVF", KTYPE_DOUBLE, "Conversion factor (e-/adu)", &(fit->keywords.cvf), NULL, NULL),
        KEYWORD_SECONDA( "image", "EGAIN", KTYPE_DOUBLE, "Conversion factor (e-/adu)", &(fit->keywords.cvf), NULL, NULL),
        KEYWORD_PRIMARY( "image", "AIRMASS", KTYPE_DOUBLE, "Airmass", &(fit->keywords.airmass), NULL, NULL),
        KEYWORD_PRIMARY( "image", "SITELAT", KTYPE_DOUBLE, "[deg] Observation site latitude", &(fit->keywords.sitelat), sitelat_handler_read, NULL),
        KEYWORD_SECONDA( "image", "SITE-LAT", KTYPE_DOUBLE, "[deg] Observation site latitude", &(fit->keywords.sitelat), sitelat_handler_read, NULL),
        KEYWORD_SECONDA( "image", "OBSLAT", KTYPE_DOUBLE, "[deg] Observation site latitude", &(fit->keywords.sitelat), sitelat_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "SITELONG", KTYPE_DOUBLE, "[deg] Observation site longitude", &(fit->keywords.sitelong), sitelong_handler_read, NULL),
        KEYWORD_SECONDA( "image", "SITE-LON", KTYPE_DOUBLE, "[deg] Observation site longitude", &(fit->keywords.sitelong), sitelong_handler_read, NULL),
        KEYWORD_SECONDA( "image", "OBSLONG", KTYPE_DOUBLE, "[deg] Observation site longitude", &(fit->keywords.sitelong), sitelong_handler_read, NULL),
        KEYWORD_PRIMARY( "image", "SITEELEV", KTYPE_DOUBLE, "[m] Observation site elevation", &(fit->keywords.siteelev), NULL, NULL),
        KEYWORD_PRIMARY( "dft",   "DFTTYPE", KTYPE_STR, "Module/Phase of a Discrete Fourier Transform", &(fit->keywords.dft.type), NULL, NULL),
        KEYWORD_PRIMARY( "dft",   "DFTORD", KTYPE_STR, "Low/High spatial freq. are located at image center", &(fit->keywords.dft.ord), NULL, NULL),
        KEYWORD_PRIMARY( "dft",   "DFTNORM1", KTYPE_DOUBLE, "Normalisation value for channel #1", &(fit->keywords.dft.norm[0]), NULL, NULL),
        KEYWORD_PRIMARY( "dft",   "DFTNORM2", KTYPE_DOUBLE, "Normalisation value for channel #2", &(fit->keywords.dft.norm[1]), NULL, NULL),
        KEYWORD_PRIMARY( "dft",   "DFTNORM3", KTYPE_DOUBLE, "Normalisation value for channel #3", &(fit->keywords.dft.norm[2]), NULL, NULL),

        KEYWORD_FIXED(   "wcsdata", "CTYPE3", KTYPE_STR, "RGB image", "RGB", NULL, NULL),
        KEYWORD_PRIMARY( "wcsdata", "OBJCTRA", KTYPE_STR, "Image center Right Ascension (hms)", &(fit->keywords.wcsdata.objctra), NULL, NULL),
        KEYWORD_PRIMARY( "wcsdata", "OBJCTDEC", KTYPE_STR, "Image center Declination (dms)", &(fit->keywords.wcsdata.objctdec), NULL, NULL),
        KEYWORD_PRIMARY( "wcsdata", "RA", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra), NULL, NULL),
        KEYWORD_SECONDA( "wcsdata", "RA_D", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra), NULL, NULL),
        KEYWORD_PRIMARY( "wcsdata", "DEC", KTYPE_DOUBLE, "Image center Declination (deg)", &(fit->keywords.wcsdata.dec), NULL, NULL),
        KEYWORD_SECONDA( "wcsdata", "DEC_D", KTYPE_DOUBLE, "Image center Declination (deg)", &(fit->keywords.wcsdata.dec), NULL, NULL),

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

        KEYWORD_WCS( "wcslib", "A_ORDER", KTYPE_INT),
        KEYWORD_WCS( "wcslib", "A_0_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_1_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_0_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_2_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_1_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_0_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_3_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_2_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_1_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_0_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_4_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_3_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_2_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_1_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_0_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_5_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_4_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_3_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_2_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_1_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "A_0_5", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_ORDER", KTYPE_INT),
        KEYWORD_WCS( "wcslib", "B_0_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_1_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_0_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_2_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_1_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_0_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_3_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_2_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_1_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_0_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_4_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_3_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_2_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_1_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_0_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_5_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_4_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_3_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_2_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_1_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "B_0_5", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_ORDER", KTYPE_INT),
        KEYWORD_WCS( "wcslib", "AP_0_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_1_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_0_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_2_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_1_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_0_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_3_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_2_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_1_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_0_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_4_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_3_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_2_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_1_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_0_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_5_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_4_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_3_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_2_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_1_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "AP_0_5", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_ORDER", KTYPE_INT),
        KEYWORD_WCS( "wcslib", "BP_0_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_1_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_0_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_2_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_1_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_0_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_3_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_2_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_1_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_0_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_4_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_3_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_2_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_1_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_0_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_5_0", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_4_1", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_3_2", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_2_3", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_1_4", KTYPE_DOUBLE),
        KEYWORD_WCS( "wcslib", "BP_0_5", KTYPE_DOUBLE),

		KEYWORD_WCS( "wcsdata", "PLTSOLVD", KTYPE_BOOL),
        {NULL, NULL, KTYPE_BOOL, NULL, NULL, NULL, FALSE, TRUE}
    };

	// Count the number of keywords in the list
	int num_keywords = 0;
	while (keyword_list[num_keywords].group != NULL) {
		num_keywords++;
	}

    // Allocate memory dynamically for the keyword array
    KeywordInfo *all_keywords = (KeywordInfo*) malloc((num_keywords + 1) * sizeof(KeywordInfo));
    if (all_keywords == NULL) {
        // Handle memory allocation failure
        return NULL;
    }

	GHashTable *hash_table = g_hash_table_new_full(g_str_hash, g_str_equal, (GDestroyNotify) g_free, NULL);

    // Copy keyword information from the list to the dynamic array and set if keyword must be used
    for (int i = 0; i < num_keywords; i++) {
        all_keywords[i] = keyword_list[i];
        all_keywords[i].is_used = should_use_keyword(fit, keyword_list[i]);

        // Set default values based on keyword type
        if (g_strcmp0(all_keywords[i].group, "wcslib")) { // This group is initialized somewhere
            switch (all_keywords[i].type) {
                case KTYPE_INT:
                    if (*((int*) all_keywords[i].data) == 0)
                        *((int*) all_keywords[i].data) = DEFAULT_INT_VALUE;
                    break;
                case KTYPE_UINT:
                    if (*((guint*) all_keywords[i].data) == 0)
                        *((guint*) all_keywords[i].data) = DEFAULT_UINT_VALUE;
                    break;
                case KTYPE_USHORT:
                    if (*((gushort*) all_keywords[i].data) == 0)
                        *((gushort*) all_keywords[i].data) = DEFAULT_USHORT_VALUE;
                    break;
                case KTYPE_DOUBLE:
                    if (*((double*) all_keywords[i].data) == 0)
                        *((double*) all_keywords[i].data) = DEFAULT_DOUBLE_VALUE;
                    break;
                case KTYPE_FLOAT:
                    if (*((float*) all_keywords[i].data) == 0)
                        *((float*) all_keywords[i].data) = DEFAULT_FLOAT_VALUE;
                    break;
                default:
                    break;
            }
        }

        g_hash_table_insert(hash_table, g_strdup(keyword_list[i].key), &(all_keywords[i]));
    }

    // Mark the end of the list
    all_keywords[num_keywords] = keyword_list[num_keywords];

    if (hash)
        *hash = hash_table;

    return all_keywords;
}

int save_fits_keywords(fits *fit) {
	KeywordInfo *keys = initialize_keywords(fit, NULL);
	KeywordInfo *keys_start = keys;
	int status;
	gchar *str;
	gushort us;
	guint ui;
	int ii;
	double dbl, zero, scale;
	float flt;
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

		/* Handle special cases */
		if (keys->special_handler_save) {
			keys->special_handler_save(fit, keys);
		}

		switch (keys->type) {
			case KTYPE_INT:
				status = 0;
				ii = (*((int*)keys->data));
				if (ii > DEFAULT_INT_VALUE) {
					fits_update_key(fit->fptr, TINT, keys->key, &ii, keys->comment, &status);
				}
				break;
			case KTYPE_UINT:
				status = 0;
				ui = (*((guint*)keys->data));
				if (ui) {
					fits_update_key(fit->fptr, TUINT, keys->key, &ui, keys->comment, &status);
				}
				break;
			case KTYPE_USHORT:
				status = 0;
				us = (*((int*)keys->data));
				if (us) {
					fits_update_key(fit->fptr, TUSHORT, keys->key, &us, keys->comment, &status);
				}
				break;
			case KTYPE_DOUBLE:
				status = 0;
				dbl = *((double*)keys->data);
				if (dbl > DEFAULT_DOUBLE_VALUE) {
					fits_update_key(fit->fptr, TDOUBLE, keys->key, &dbl, keys->comment, &status);
				}
				break;
			case KTYPE_FLOAT:
				status = 0;
				flt = *((float*)keys->data);
				if (flt > DEFAULT_FLOAT_VALUE) {
					fits_update_key(fit->fptr, TFLOAT, keys->key, &flt, keys->comment, &status);
				}
				break;
			case KTYPE_STR:
				status = 0;
				str = ((gchar*)keys->data);
				if (str && str[0] != '\0') {
					fits_update_key(fit->fptr, TSTRING, keys->key, str, keys->comment, &status);
				}
				break;
			case KTYPE_DATE:
				status = 0;
				date = *((GDateTime**)keys->data);
				if (date) {
					/* special case */
					if (!g_strcmp0("DATE", keys->key)) {
						int itmp;
						char fit_date[40];
						fits_get_system_time(fit_date, &itmp, &status);
						fits_update_key(fit->fptr, TSTRING, keys->key, fit_date, keys->comment, &status);
					} else {
						gchar *formatted_date = date_time_to_FITS_date(date);
						fits_update_key(fit->fptr, TSTRING, keys->key, formatted_date, keys->comment, &status);
						g_free(formatted_date);
					}
				}
				break;
			case KTYPE_BOOL:
				status = 0;
				fits_update_key(fit->fptr, TLOGICAL, keys->key, &(*((gboolean*)keys->data)), keys->comment, &status);
				break;
			default:
				siril_debug_print("Save_fits_keywords: Error. Type is not handled.\n");
		}
		keys++;
	}

	free(keys_start);

	return 0;
}

int save_wcs_keywords(fits *fit) {
	int status = 0;
	/* Needed for Aladin compatibility */
	if (fit->naxes[2] == 3 && com.pref.rgb_aladin) {
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CTYPE3", "RGB", "RGB image", &status);
	}

	if (fit->keywords.wcsdata.objctra[0] != '\0') {
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "OBJCTRA", &(fit->keywords.wcsdata.objctra),	"Image center Right Ascension (hms)", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "OBJCTDEC", &(fit->keywords.wcsdata.objctdec), "Image center Declination (dms)", &status);
	}
	if (fit->keywords.wcsdata.ra > 0) {
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "RA", &(fit->keywords.wcsdata.ra), "Image center Right Ascension (deg)", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "DEC", &(fit->keywords.wcsdata.dec), "Image center Declination (deg)", &status);
	}
	status = 0;

	if (fit->keywords.wcslib) {
		gboolean has_sip = fit->keywords.wcslib->lin.dispre != NULL; // we don't handle the disseq terms for now
		if (!has_sip) {// no distortions
			fits_update_key(fit->fptr, TSTRING, "CTYPE1", "RA---TAN", "TAN (gnomic) projection", &status);
			status = 0;
			fits_update_key(fit->fptr, TSTRING, "CTYPE2", "DEC--TAN", "TAN (gnomic) projection", &status);
			status = 0;
		} else {
			fits_update_key(fit->fptr, TSTRING, "CTYPE1", "RA---TAN-SIP", "TAN (gnomic) projection + SIP distortions", &status);
			status = 0;
			fits_update_key(fit->fptr, TSTRING, "CTYPE2", "DEC--TAN-SIP", "TAN (gnomic) projection + SIP distortions", &status);
			status = 0;
		}
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CUNIT1", "deg","Unit of coordinates", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CUNIT2", "deg","Unit of coordinates", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "EQUINOX", &(fit->keywords.wcslib->equinox),	"Equatorial equinox", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRPIX1", &(fit->keywords.wcslib->crpix[0]), "Axis1 reference pixel", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRPIX2", &(fit->keywords.wcslib->crpix[1]), "Axis2 reference pixel", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL1", &(fit->keywords.wcslib->crval[0]), "Axis1 reference value (deg)", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL2", &(fit->keywords.wcslib->crval[1]), "Axis2 reference value (deg)", &status);
		if (fit->keywords.wcslib->lonpole) {
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "LONPOLE", &(fit->keywords.wcslib->lonpole), "Native longitude of celestial pole", &status);
		}
		if (com.pref.wcs_formalism == WCS_FORMALISM_1) {
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CDELT1", &(fit->keywords.wcslib->cdelt[0]), "X pixel size (deg)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CDELT2", &(fit->keywords.wcslib->cdelt[1]), "Y pixel size (deg)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC1_1", &(fit->keywords.wcslib->pc[0]), "Linear transformation matrix (1, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC1_2", &(fit->keywords.wcslib->pc[1]), "Linear transformation matrix (1, 2)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC2_1", &(fit->keywords.wcslib->pc[2]), "Linear transformation matrix (2, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC2_2", &(fit->keywords.wcslib->pc[3]), "Linear transformation matrix (2, 2)", &status);
			status = 0;
		} else {
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD1_1", &(fit->keywords.wcslib->cd[0]), "Scale matrix (1, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD1_2", &(fit->keywords.wcslib->cd[1]), "Scale matrix (1, 2)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD2_1", &(fit->keywords.wcslib->cd[2]), "Scale matrix (2, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD2_2", &(fit->keywords.wcslib->cd[3]), "Scale matrix (2, 2)", &status);
			status = 0;
		}
		if (has_sip) {
			// we deal with images up to order 6, we need 7 to hold 0_6 terms
			double A[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
			double B[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
			double AP[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
			double BP[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
			struct disprm *dis = fit->keywords.wcslib->lin.dispre;
			int order = extract_SIP_order_and_matrices(dis, A, B, AP, BP);
			// we know the order of the distortions, we can now write them
			// A terms
			fits_update_key(fit->fptr, TINT, "A_ORDER", &order, "SIP polynomial degree, axis 1, pixel-to-sky", &status);
			for (int i = 0; i <= order; i++) {
				for (int j = i; j >= 0; j--) {
					int k = i - j;
					char key[6];
					g_snprintf(key, 6, "A_%d_%d", j, k);
					fits_update_key(fit->fptr, TDOUBLE, key, &A[j][k], NULL, &status);
				}
			}
			// B terms
			fits_update_key(fit->fptr, TINT, "B_ORDER", &order, "SIP polynomial degree, axis 2, pixel-to-sky", &status);
			for (int i = 0; i <= order; i++) {
				for (int j = i; j >= 0; j--) {
					int k = i - j;
					char key[6];
					g_snprintf(key, 6, "B_%d_%d", j, k);
					fits_update_key(fit->fptr, TDOUBLE, key, &B[j][k], NULL, &status);
				}
			}
			// AP terms
			fits_update_key(fit->fptr, TINT, "AP_ORDER", &order, "SIP polynomial degree, axis 1, sky-to-pixel", &status);
			for (int i = 0; i <= order; i++) {
				for (int j = i; j >= 0; j--) {
					int k = i - j;
					char key[7];
					g_snprintf(key, 7, "AP_%d_%d", j, k);
					fits_update_key(fit->fptr, TDOUBLE, key, &AP[j][k], NULL, &status);
				}
			}
			// BP terms
			fits_update_key(fit->fptr, TINT, "BP_ORDER", &order, "SIP polynomial degree, axis 2, sky-to-pixel", &status);
			for (int i = 0; i <= order; i++) {
				for (int j = i; j >= 0; j--) {
					int k = i - j;
					char key[7];
					g_snprintf(key, 7, "BP_%d_%d", j, k);
					fits_update_key(fit->fptr, TDOUBLE, key, &BP[j][k], NULL, &status);
				}
			}
		}
	}
	if (fit->keywords.wcsdata.pltsolvd) {
		fits_update_key(fit->fptr, TLOGICAL, "PLTSOLVD", &(fit->keywords.wcsdata.pltsolvd), fit->keywords.wcsdata.pltsolvd_comment, &status);
	}

	return 0;
}


int save_fits_unknown_keywords(fits *fit) {
	int status = 0;
	/*** Save list of unknown keys ***/
	if (fit->unknown_keys) {
		status = 0;
		fits_write_comment(fit->fptr, UNKNOWN_KEYWORDS, &status);
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

void read_fits_date_obs_header(fits *fit) {
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
	// Initialize keywords and get hash table
	GHashTable *keys_hash;
	KeywordInfo *keys = initialize_keywords(fit, &keys_hash);
	int status = 0;
	int key_number = 1;
	gboolean end_of_header = FALSE;
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

	read_fits_date_obs_header(fit);

	fits_get_hdrspace(fit->fptr, &key_number, NULL, &status); /* get # of keywords */

	// Loop through each keyword
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

		if (g_strcmp0(comment, UNKNOWN_KEYWORDS) == 0) {
			end_of_header = TRUE;
			continue;
		}

		/* These have been already processed */
		if (g_strcmp0(keyname, "BSCALE") == 0 || g_strcmp0(keyname, "BZERO") || g_strcmp0(keyname, "DATE-OBS") == 0) {
			continue;
		}

		// Retrieve KeywordInfo from the hash table
		KeywordInfo *current_key = g_hash_table_lookup(keys_hash, keyname);
		if (((current_key == NULL && fits_get_keyclass(card) != TYP_STRUC_KEY) || end_of_header)) {
			if (strncmp(card, "HISTORY", 7) == 0) continue;
			// Handle unknown keys
			UnknownKeys = g_string_append(UnknownKeys, card);
			UnknownKeys = g_string_append(UnknownKeys, "\n");
			continue;
		}

		gboolean value_set = FALSE; // Flag indicating whether a value has been set.
		if (fits_get_keyclass(card) == TYP_STRUC_KEY) {
			value_set = TRUE;
			continue;
		} else if (current_key->fixed_value) {
			value_set = TRUE;
			continue;
		}
		int int_value;
		guint uint_value;
		gushort ushort_value;
		double double_value;
		float float_value;
		gchar *str_value;
		GDateTime *date;
		gboolean bool_value;
		char *end;

		// Process the value based on the type of the keyword
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
		case KTYPE_FLOAT:
			float_value = g_ascii_strtod(value, &end);
			if (value != end) {
				*((float*) current_key->data) = float_value;
				value_set = TRUE;
			}
			break;
		case KTYPE_STR:
			str_value = g_strstrip(g_shell_unquote(value, NULL));
			strncpy((char*) current_key->data, str_value, FLEN_VALUE - 1);
			value_set = TRUE;
			break;
		case KTYPE_DATE:
			str_value = g_strstrip(g_shell_unquote(value, NULL));
			date = FITS_date_to_date_time(str_value);
			if (date) {
				*((GDateTime**) current_key->data) = date;
				value_set = TRUE;
			}
			break;
		case KTYPE_BOOL:
			bool_value = value[0] == 'T' ? TRUE : FALSE;
			*((gboolean*) current_key->data) = bool_value;
			value_set = TRUE;
			break;
		default:
			break;
		}

		// Handle special cases
		if (current_key->special_handler_read != NULL) {
			current_key->special_handler_read(fit, comment, current_key);
		}

		/* output unknown keywords */
		if (!value_set) {
			UnknownKeys = g_string_append(UnknownKeys, card);
			UnknownKeys = g_string_append(UnknownKeys, "\n");
		}
	}

	if (fit->unknown_keys != NULL) {
		g_free(fit->unknown_keys);
	}
	fit->unknown_keys = g_string_free(UnknownKeys, FALSE);

	// Free the hash table and unknown keys
	g_hash_table_destroy(keys_hash);
	free(keys);
	return 0;
}

