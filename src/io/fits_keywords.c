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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "core/siril_world_cs.h"
#include "algos/siril_wcs.h"
#include "gui/keywords_tree.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/path_parse.h"

#include "fits_keywords.h"

// Uncomment to print debug verbose
// #define DEBUG_PRINT_HEADER

#define PRINT_PARSING_ERROR \
		do { \
				siril_debug_print("Error parsing value: %s\n", value); \
				siril_debug_print("Parsing stopped at: %s\n", end); \
		} while (0)


#define KEYWORD_PRIMARY(group, key, type, comment, data, handler_read, handler_save) { group, key, type, comment, data, handler_read, handler_save, TRUE, FALSE }
#define KEYWORD_SECONDA(group, key, type, comment, data, handler_read, handler_save) { group, key, type, comment, data, handler_read, handler_save, FALSE, FALSE }
#define KEYWORD_FIXED(group, key, type, comment, data, handler_read, handler_save) { group, key, type, comment, data, handler_read, handler_save, TRUE, TRUE }
#define KEYWORD_WCS(group, key, type) { group, key, type, NULL, NULL, NULL, NULL, FALSE, TRUE }

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
	} else if (g_strcmp0(keyword.key, "XBAYROFF") == 0) {
		return fit->keywords.bayer_pattern[0] != '\0';
	} else if (g_strcmp0(keyword.key, "YBAYROFF") == 0) {
		return fit->keywords.bayer_pattern[0] != '\0';
	} else if (g_strcmp0(keyword.key, "CTYPE3") == 0) {
		return (fit->naxes[2] == 3 && com.pref.rgb_aladin);
	}
	return keyword.is_saved;
}

/****************************** handlers ******************************/

static void bscale_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (1.0 != fit->keywords.bscale) {
		siril_log_message(_("Loaded FITS file has a BSCALE different than 1 (%f)\n"), fit->keywords.bscale);
		int status = 0;
		/* We reset the scaling factors as we don't use it */
		fits_set_bscale(fit->fptr, 1.0, 0.0, &status);
	}
}

static void bzero_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	if (0.0 != fit->keywords.bzero && fit->bitpix == FLOAT_IMG) {
		fprintf(stdout, "ignoring BZERO\n");
		int status = 0;
		fits_set_bscale(fit->fptr, 1.0, 0.0, &status);
	}
}

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
	char sitelong_dump_tmp[FLEN_VALUE] = { 0 };
	double d_sitelong_dump = 0.0;

	gchar **token = g_strsplit(fit->keywords.sitelong_str, ":", -1); // Handles PRISM special parsing for SITELONG
	gsize token_size = g_strv_length(token);
	if (token_size > 1 && token[1])	{
		for (int i = 0; i < token_size; ++i) {
			if (g_strlcat(sitelong_dump_tmp, token[i], sizeof(sitelong_dump_tmp)) >= sizeof(sitelong_dump_tmp))
				siril_debug_print("Truncation occurred in g_strlcat\n");
			if (i < 3) strncat(sitelong_dump_tmp, i < 2 ? ":" : ".", 2);
			d_sitelong_dump = parse_dms(sitelong_dump_tmp);
		}
	} else d_sitelong_dump = parse_dms(fit->keywords.sitelong_str);

	g_strfreev(token);

	if (isnan(d_sitelong_dump)) {	// Cases SITELONG and SITELAT keyword are numbers (only NINA and Seq. Generator, for now)
		gchar *end;
		fit->keywords.sitelong = g_ascii_strtod(fit->keywords.sitelong_str, &end);
		if (fit->keywords.sitelong_str == end) {
			fit->keywords.sitelong = DEFAULT_DOUBLE_VALUE;
			info->used = FALSE;
			siril_debug_print("Cannot read SITELONG\n");
		}
	} else {
		fit->keywords.sitelong = d_sitelong_dump;
	}
}

static void sitelat_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	char sitelat_dump_tmp[FLEN_VALUE] = { 0 };
	double d_sitelat_dump = 0.0;

	gchar **token = g_strsplit(fit->keywords.sitelat_str, ":", -1); // Handles PRISM special parsing for SITELAT
	gsize token_size = g_strv_length(token);
	if (token_size > 1 && token[1])	{	// Denotes presence of ":"
		for (int i = 0; i < token_size; ++i) {
			if (g_strlcat(sitelat_dump_tmp, token[i], sizeof(sitelat_dump_tmp)) >= sizeof(sitelat_dump_tmp))
				siril_debug_print("Truncation occurred in g_strlcat\n");
			if (i < 3) strncat(sitelat_dump_tmp, i < 2 ? ":" : ".", 2);
			d_sitelat_dump = parse_dms(sitelat_dump_tmp);
		}
	} else d_sitelat_dump = parse_dms(fit->keywords.sitelat_str);

	g_strfreev(token);

	if (isnan(d_sitelat_dump)) {	// Cases SITELONG and SITELAT keyword are numbers (only NINA and Seq. Generator, for now)
		gchar *end;
		fit->keywords.sitelat = g_ascii_strtod(fit->keywords.sitelat_str, &end);
		if (fit->keywords.sitelat_str == end) {
			fit->keywords.sitelat = DEFAULT_DOUBLE_VALUE;
			info->used = FALSE;
			siril_debug_print("Cannot read SITELAT\n");
		}
	} else {
		fit->keywords.sitelat = d_sitelat_dump;
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

static void ra_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	fit->keywords.wcsdata.ra = parse_hms(fit->keywords.wcsdata.objctra);
	if (isnan(fit->keywords.wcsdata.ra)) {
		fit->keywords.wcsdata.ra = DEFAULT_DOUBLE_VALUE;
		info->used = FALSE;
	}
	else
		siril_debug_print("read RA as HMS\n");
}

static void dec_handler_read(fits *fit, const char *comment, KeywordInfo *info) {
	fit->keywords.wcsdata.dec = parse_dms(fit->keywords.wcsdata.objctdec);
	if (isnan(fit->keywords.wcsdata.dec)) {
		fit->keywords.wcsdata.dec = DEFAULT_DOUBLE_VALUE;
		info->used = FALSE;
	}
	else
		siril_debug_print("read DEC as DMS\n");

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

static void bzero_handler_save(fits *fit, KeywordInfo *info) {
	switch (fit->bitpix) {
	default:
	case BYTE_IMG:
	case SHORT_IMG:
	case FLOAT_IMG:
		fit->keywords.bzero = 0.0;
		break;
	case USHORT_IMG:
		fit->keywords.bzero = 32768.0;
		break;
	}
}

static void bscale_handler_save(fits *fit, KeywordInfo *info) {
	fit->keywords.bscale = 1.0;
}

static void program_handler_save(fits *fit, KeywordInfo *info) {
	strncpy(fit->keywords.program, "Siril "PACKAGE_VERSION, FLEN_VALUE - 1);
}

static void focal_length_handler_save(fits *fit, KeywordInfo *info) {
	info->is_saved = fit->focalkey;
}

static void pixel_x_handler_save(fits *fit, KeywordInfo *info) {
	info->is_saved = fit->pixelkey;
}

static void pixel_y_handler_save(fits *fit, KeywordInfo *info) {
	info->is_saved = fit->pixelkey;
}

/*****************************************************************************/

KeywordInfo *initialize_keywords(fits *fit, GHashTable **hash) {
	/*
	 * KEYWORD_PRIMARY represents keywords read and saved by Siril
	 * KEYWORD_SECONDA are those read but not saved by Siril
	 * KEYWORD_FIXED are keywords whose value is fixed and does not change.
	 * KEYWORD_WCS Used for keywords in the wcslib group. They are recognized as known
	 * keywords but are read by another routine. They are also saved in a special function.
	 *
	 * A series of keywords representing the same data usually begins with KEYWORD_PRIMARY and
	 * is followed by a list of KEYWORD_SECONDA. They should normally be linked to the same
	 * variable in the keywords structure.
	 * The only exception is FLENGHT, where the units are not the same. We could have managed
	 * the units to avoid this, but the software that manages these keywords doesn't save
	 * the units in compliance with the standard. So we abandoned the idea. It's something we
	 * can do if we have the data to do it. All we have to do is retrieve the units and add a
	 * variable to the read handle.
	 */
	KeywordInfo keyword_list[] = {
			KEYWORD_PRIMARY( "image", "BZERO", KTYPE_DOUBLE, "Offset data range to that of unsigned short", &(fit->keywords.bzero), bzero_handler_read, bzero_handler_save),
			KEYWORD_PRIMARY( "image", "BSCALE", KTYPE_DOUBLE, "Default scaling factor", &(fit->keywords.bscale), bscale_handler_read, bscale_handler_save),
			KEYWORD_PRIMARY( "image", "MIPS-HI", KTYPE_USHORT, "Lower visualization cutoff", &(fit->keywords.hi), NULL, NULL),
			KEYWORD_SECONDA( "image", "CWHITE", KTYPE_USHORT, "Lower visualization cutoff", &(fit->keywords.hi), NULL, NULL),
			KEYWORD_PRIMARY( "image", "MIPS-LO", KTYPE_USHORT, "Upper visualization cutoff", &(fit->keywords.lo), NULL, NULL),
			KEYWORD_SECONDA( "image", "CBLACK", KTYPE_USHORT, "Upper visualization cutoff", &(fit->keywords.lo), NULL, NULL),
			KEYWORD_PRIMARY( "image", "MIPS-FLO", KTYPE_FLOAT, "Lower visualization cutoff", &(fit->keywords.flo), flo_handler_read, flo_handler_save),
			KEYWORD_PRIMARY( "image", "MIPS-FHI", KTYPE_FLOAT, "Upper visualization cutoff", &(fit->keywords.fhi), fhi_handler_read, fhi_handler_save),
			/* ATTENTION: PROGRAM MUST BE BEFORE DATAMAX */
			KEYWORD_PRIMARY( "image", "PROGRAM", KTYPE_STR, "Software that created this HDU", &(fit->keywords.program), NULL, program_handler_save),
			KEYWORD_SECONDA( "image", "SWCREATE", KTYPE_STR, "Software that created this HDU", &(fit->keywords.program), NULL, NULL),
			/* We do not want to save datamax */
			KEYWORD_PRIMARY( "image", "FILENAME", KTYPE_STR, "", &(fit->keywords.filename), NULL, NULL), // no comment to save space for filename
			KEYWORD_SECONDA( "image", "DATAMAX", KTYPE_DOUBLE, "Maximum pixel value", &(fit->keywords.data_max), NULL, NULL),
			KEYWORD_PRIMARY( "date",  "DATE", KTYPE_DATE, "UTC date that FITS file was created", &(fit->keywords.date), NULL, NULL),
			KEYWORD_PRIMARY( "date",  "DATE-OBS", KTYPE_DATE, "YYYY-MM-DDThh:mm:ss observation start, UT", &(fit->keywords.date_obs), NULL, NULL),
			KEYWORD_PRIMARY( "image", "IMAGETYP", KTYPE_STR, "Type of image", &(fit->keywords.image_type), NULL, NULL),
			KEYWORD_SECONDA( "image", "FRAMETYP", KTYPE_STR, "Type of image", &(fit->keywords.image_type), NULL, NULL),
			KEYWORD_PRIMARY( "image", "ROWORDER", KTYPE_STR, "Order of the rows in image array", &(fit->keywords.row_order), roworder_handler_read, NULL),
			KEYWORD_PRIMARY( "image", "EXPTIME", KTYPE_DOUBLE, "[s] Exposure time duration", &(fit->keywords.exposure), NULL, NULL),
			KEYWORD_SECONDA( "image", "EXPOSURE", KTYPE_DOUBLE, "[s] Exposure time duration", &(fit->keywords.exposure), NULL, NULL),
			KEYWORD_PRIMARY( "telescope", "TELESCOP", KTYPE_STR, "Telescope used to acquire this image", &(fit->keywords.telescop), NULL, NULL),
			KEYWORD_PRIMARY( "telescope", "OBSERVER", KTYPE_STR, "Observer name", &(fit->keywords.observer), NULL, NULL),
			KEYWORD_PRIMARY( "telescope", "FILTER", KTYPE_STR, "Active filter name", &(fit->keywords.filter), NULL, NULL),
			KEYWORD_SECONDA( "telescope", "FILT-1", KTYPE_STR, "Active filter name", &(fit->keywords.filter), NULL, NULL),
			KEYWORD_PRIMARY( "telescope", "APERTURE", KTYPE_DOUBLE, "Aperture of the instrument", &(fit->keywords.aperture), NULL, NULL),
			KEYWORD_PRIMARY( "telescope", "ISOSPEED", KTYPE_DOUBLE, "ISO camera setting", &(fit->keywords.iso_speed), NULL, NULL),
			KEYWORD_PRIMARY( "telescope", "FOCALLEN", KTYPE_DOUBLE, "[mm] Focal length", &(fit->keywords.focal_length), focal_length_handler_read, focal_length_handler_save),
			KEYWORD_SECONDA( "telescope", "FOCAL", KTYPE_DOUBLE, "[mm] Focal length", &(fit->keywords.focal_length), focal_length_handler_read, NULL),
			KEYWORD_SECONDA( "telescope", "FLENGTH", KTYPE_DOUBLE, "[m] Focal length", &(fit->keywords.flength), flength_handler_read, NULL),
			KEYWORD_PRIMARY( "telescope", "CENTALT", KTYPE_DOUBLE, "[deg] Altitude of telescope", &(fit->keywords.centalt), NULL, NULL),
			KEYWORD_SECONDA( "telescope", "OBJCTALT", KTYPE_DOUBLE, "[deg] Altitude of telescope", &(fit->keywords.centalt), NULL, NULL),
			KEYWORD_PRIMARY( "telescope", "CENTAZ", KTYPE_DOUBLE, "[deg] Azimuth of telescope", &(fit->keywords.centaz), NULL, NULL),
			KEYWORD_SECONDA( "telescope", "OBJCTAZ", KTYPE_DOUBLE, "[deg] Azimuth of telescope", &(fit->keywords.centaz), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "XBINNING", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_x), binning_x_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "BINX", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_x), binning_x_handler_read, NULL),
			KEYWORD_PRIMARY( "camera", "YBINNING", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_y), binning_y_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "BINY", KTYPE_UINT, "Camera binning mode", &(fit->keywords.binning_y), binning_y_handler_read, NULL),
			KEYWORD_PRIMARY( "camera", "XPIXSZ", KTYPE_DOUBLE,   "[um] Pixel X axis size", &(fit->keywords.pixel_size_x), pixel_x_handler_read, pixel_x_handler_save),
			KEYWORD_SECONDA( "camera", "XPIXELSZ", KTYPE_DOUBLE, "[um] Pixel X axis size", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "PIXSIZE1", KTYPE_DOUBLE, "[um] Pixel X axis size", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "PIXSIZEX", KTYPE_DOUBLE, "[um] Pixel X axis size", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "XPIXSIZE", KTYPE_DOUBLE, "[um] Pixel X axis size", &(fit->keywords.pixel_size_x), pixel_x_handler_read, NULL),
			KEYWORD_PRIMARY( "camera", "YPIXSZ", KTYPE_DOUBLE,   "[um] Pixel Y axis size", &(fit->keywords.pixel_size_y), pixel_x_handler_read, pixel_y_handler_save),
			KEYWORD_SECONDA( "camera", "YPIXELSZ", KTYPE_DOUBLE, "[um] Pixel Y axis size", &(fit->keywords.pixel_size_y), pixel_x_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "PIXSIZE2", KTYPE_DOUBLE, "[um] Pixel Y axis size", &(fit->keywords.pixel_size_y), pixel_x_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "PIXSIZEY", KTYPE_DOUBLE, "[um] Pixel Y axis size", &(fit->keywords.pixel_size_y), pixel_x_handler_read, NULL),
			KEYWORD_SECONDA( "camera", "YPIXSIZE", KTYPE_DOUBLE, "[um] Pixel Y axis size", &(fit->keywords.pixel_size_y), pixel_x_handler_read, NULL),
			KEYWORD_PRIMARY( "camera", "INSTRUME", KTYPE_STR, "Instrument name", &(fit->keywords.instrume), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "CCD-TEMP", KTYPE_DOUBLE, "[degC] CCD temperature", &(fit->keywords.ccd_temp), NULL, NULL),
			KEYWORD_SECONDA( "camera", "CCD_TEMP", KTYPE_DOUBLE, "[degC] CCD temperature", &(fit->keywords.ccd_temp), NULL, NULL),
			KEYWORD_SECONDA( "camera", "CCDTEMP", KTYPE_DOUBLE, "[degC] CCD temperature", &(fit->keywords.ccd_temp), NULL, NULL),
			KEYWORD_SECONDA( "camera", "TEMPERAT", KTYPE_DOUBLE, "[degC] CCD temperature", &(fit->keywords.ccd_temp), NULL, NULL),
			KEYWORD_SECONDA( "camera", "CAMTCCD", KTYPE_DOUBLE, "[degC] CCD temperature", &(fit->keywords.ccd_temp), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "SET-TEMP", KTYPE_DOUBLE, "[degC] CCD temperature setpoint", &(fit->keywords.set_temp), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "GAIN", KTYPE_UINT, "Sensor gain", &(fit->keywords.key_gain), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "OFFSET", KTYPE_UINT, "Sensor gain offset", &(fit->keywords.key_offset), NULL, NULL),
			KEYWORD_SECONDA( "camera", "BLKLEVEL", KTYPE_UINT, "Sensor gain offset", &(fit->keywords.key_offset), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "CVF", KTYPE_DOUBLE, "[e-/ADU] Electrons per A/D unit", &(fit->keywords.cvf), NULL, NULL),
			KEYWORD_SECONDA( "camera", "EGAIN", KTYPE_DOUBLE, "[e-/ADU] Electrons per A/D unit", &(fit->keywords.cvf), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "BAYERPAT", KTYPE_STR, "Bayer color pattern", &(fit->keywords.bayer_pattern), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "XBAYROFF", KTYPE_INT, "X offset of Bayer array", &(fit->keywords.bayer_xoffset), NULL, NULL),
			KEYWORD_PRIMARY( "camera", "YBAYROFF", KTYPE_INT, "Y offset of Bayer array", &(fit->keywords.bayer_yoffset), NULL, NULL),
			KEYWORD_PRIMARY( "focuser", "FOCNAME", KTYPE_STR, "Focusing equipment name", &(fit->keywords.focname), NULL, NULL),
			KEYWORD_PRIMARY( "focuser", "FOCPOS", KTYPE_INT, "[step] Focuser position", &(fit->keywords.focuspos), NULL, NULL),
			KEYWORD_SECONDA( "focuser", "FOCUSPOS", KTYPE_INT, "[step] Focuser position", &(fit->keywords.focuspos), NULL, NULL),
			KEYWORD_PRIMARY( "focuser", "FOCUSSZ", KTYPE_INT, "[um] Focuser step size", &(fit->keywords.focussz), NULL, NULL),
			KEYWORD_PRIMARY( "focuser", "FOCTEMP", KTYPE_DOUBLE, "[degC] Focuser temperature", &(fit->keywords.foctemp), NULL, NULL),
			KEYWORD_SECONDA( "focuser", "FOCUSTEM", KTYPE_DOUBLE, "[degC] Focuser temperature", &(fit->keywords.foctemp), NULL, NULL),
			KEYWORD_PRIMARY( "stack", "STACKCNT", KTYPE_UINT, "Stack frames", &(fit->keywords.stackcnt), NULL, NULL),
			KEYWORD_SECONDA( "stack", "NCOMBINE", KTYPE_UINT, "Stack frames", &(fit->keywords.stackcnt), NULL, NULL),
			KEYWORD_PRIMARY( "stack", "LIVETIME", KTYPE_DOUBLE, "[s] Exposure time after deadtime correction", &(fit->keywords.livetime), NULL, NULL),
			KEYWORD_PRIMARY( "stack", "EXPSTART", KTYPE_DOUBLE, "[JD] Exposure start time (standard Julian date)", &(fit->keywords.expstart), NULL, NULL),
			KEYWORD_PRIMARY( "stack", "EXPEND", KTYPE_DOUBLE, "[JD] Exposure end time (standard Julian date)", &(fit->keywords.expend), NULL, NULL),
			KEYWORD_PRIMARY( "target", "OBJECT", KTYPE_STR, "Name of the object of interest", &(fit->keywords.object), NULL, NULL),
			KEYWORD_PRIMARY( "target", "AIRMASS", KTYPE_DOUBLE, "Airmass at frame center (Gueymard 1993)", &(fit->keywords.airmass), NULL, NULL),
			/* SITELAT and SITELONG can be provided in double or string. We need to check both. We save as double */
			KEYWORD_SECONDA( "geo", "SITELAT", KTYPE_STR, "[deg] Observation site latitude", &(fit->keywords.sitelat_str), sitelat_handler_read, NULL),
			KEYWORD_SECONDA( "geo", "SITE-LAT", KTYPE_STR, "[deg] Observation site latitude", &(fit->keywords.sitelat_str), sitelat_handler_read, NULL),
			KEYWORD_SECONDA( "geo", "OBSLAT", KTYPE_STR, "[deg] Observation site latitude", &(fit->keywords.sitelat_str), sitelat_handler_read, NULL),
			KEYWORD_PRIMARY( "geo", "SITELAT", KTYPE_DOUBLE, "[deg] Observation site latitude", &(fit->keywords.sitelat), NULL, NULL),
			KEYWORD_SECONDA( "geo", "SITELONG", KTYPE_STR, "[deg] Observation site longitude", &(fit->keywords.sitelong_str), sitelong_handler_read, NULL),
			KEYWORD_SECONDA( "geo", "SITE-LON", KTYPE_STR, "[deg] Observation site longitude", &(fit->keywords.sitelong_str), sitelong_handler_read, NULL),
			KEYWORD_SECONDA( "geo", "OBSLONG", KTYPE_STR, "[deg] Observation site longitude", &(fit->keywords.sitelong_str), sitelong_handler_read, NULL),
			KEYWORD_PRIMARY( "geo", "SITELONG", KTYPE_DOUBLE, "[deg] Observation site longitude", &(fit->keywords.sitelong), NULL, NULL),
			KEYWORD_PRIMARY( "geo", "SITEELEV", KTYPE_DOUBLE, "[m] Observation site elevation", &(fit->keywords.siteelev), NULL, NULL),
			KEYWORD_PRIMARY( "dft", "DFTTYPE", KTYPE_STR, "Module/Phase of a Discrete Fourier Transform", &(fit->keywords.dft.type), NULL, NULL),
			KEYWORD_PRIMARY( "dft", "DFTORD", KTYPE_STR, "Low/High spatial freq. are located at image center", &(fit->keywords.dft.ord), NULL, NULL),
			KEYWORD_PRIMARY( "dft", "DFTNORM1", KTYPE_DOUBLE, "Normalisation value for channel #1", &(fit->keywords.dft.norm[0]), NULL, NULL),
			KEYWORD_PRIMARY( "dft", "DFTNORM2", KTYPE_DOUBLE, "Normalisation value for channel #2", &(fit->keywords.dft.norm[1]), NULL, NULL),
			KEYWORD_PRIMARY( "dft", "DFTNORM3", KTYPE_DOUBLE, "Normalisation value for channel #3", &(fit->keywords.dft.norm[2]), NULL, NULL),

			KEYWORD_FIXED(   "wcsdata", "CTYPE3", KTYPE_STR, "RGB image", "RGB", NULL, NULL),
			KEYWORD_PRIMARY( "wcsdata", "OBJCTRA", KTYPE_STR, "Image center Right Ascension (hms)", &(fit->keywords.wcsdata.objctra), NULL, NULL),
			KEYWORD_PRIMARY( "wcsdata", "OBJCTDEC", KTYPE_STR, "Image center Declination (dms)", &(fit->keywords.wcsdata.objctdec), NULL, NULL),
			KEYWORD_SECONDA( "wcsdata", "RA", KTYPE_STR, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.objctra), ra_handler_read, NULL),
			KEYWORD_PRIMARY( "wcsdata", "RA", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra), NULL, NULL),
			KEYWORD_SECONDA( "wcsdata", "RA_D", KTYPE_DOUBLE, "Image center Right Ascension (deg)", &(fit->keywords.wcsdata.ra), NULL, NULL),
			KEYWORD_SECONDA( "wcsdata", "DEC", KTYPE_STR, "Image center Declination (deg)", &(fit->keywords.wcsdata.objctdec), dec_handler_read, NULL),
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
			{NULL, NULL, KTYPE_BOOL, NULL, NULL, NULL, FALSE, TRUE }
	};

	int num_keywords = G_N_ELEMENTS(keyword_list) - 1; // remove last line

	// Allocate memory dynamically for the keyword array
	KeywordInfo *all_keywords = (KeywordInfo*) malloc((num_keywords + 1) * sizeof(KeywordInfo));
	if (!all_keywords) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	GHashTable *hash_table = NULL;

	if (hash)
		hash_table = g_hash_table_new_full(g_str_hash, g_str_equal, (GDestroyNotify) g_free, NULL);

	// Copy keyword information from the list to the dynamic array and set if keyword must be used
	for (int i = 0; i < num_keywords; i++) {
		all_keywords[i] = keyword_list[i];
		all_keywords[i].is_saved = should_use_keyword(fit, keyword_list[i]);

		if (hash_table)
			g_hash_table_insert(hash_table, g_strdup(keyword_list[i].key), &(all_keywords[i]));
	}

	// Mark the end of the list
	all_keywords[num_keywords] = keyword_list[num_keywords];

	if (hash_table)
		*hash = hash_table;

	return all_keywords;
}

int save_fits_keywords(fits *fit) {
	KeywordInfo *keys = initialize_keywords(fit, NULL);
	KeywordInfo *keys_start = keys;
	super_bool sbool;
	gboolean bool;
	int status;
	gchar *str;
	gushort us;
	guint ui;
	int ii;
	double dbl;
	float flt;
	GDateTime *date = NULL;

	while (keys->group) {
		/* Handle special cases */
		if (keys->special_handler_save) {
			keys->special_handler_save(fit, keys);
		}

		if (!keys->is_saved || g_strcmp0(keys->group, "wcslib") == 0) {
			keys++;
			continue;
		}

		switch (keys->type) {
		case KTYPE_INT:
			status = 0;
			ii = (*((int*) keys->data));
			if (ii > DEFAULT_INT_VALUE) {
				fits_update_key(fit->fptr, TINT, keys->key, &ii, keys->comment, &status);
			}
			break;
		case KTYPE_UINT:
			status = 0;
			ui = (*((guint*) keys->data));
			if (ui != DEFAULT_UINT_VALUE) {
				fits_update_key(fit->fptr, TUINT, keys->key, &ui, keys->comment, &status);
			}
			break;
		case KTYPE_USHORT:
			status = 0;
			us = (*((int*) keys->data));
			if (us) {
				fits_update_key(fit->fptr, TUSHORT, keys->key, &us, keys->comment, &status);
			}
			break;
		case KTYPE_DOUBLE:
			status = 0;
			dbl = *((double*) keys->data);
			if (dbl > DEFAULT_DOUBLE_VALUE) {
				fits_update_key(fit->fptr, TDOUBLE, keys->key, &dbl, keys->comment, &status);
			}
			break;
		case KTYPE_FLOAT:
			status = 0;
			flt = *((float*) keys->data);
			if (flt > DEFAULT_FLOAT_VALUE) {
				fits_update_key(fit->fptr, TFLOAT, keys->key, &flt, keys->comment, &status);
			}
			break;
		case KTYPE_STR:
			status = 0;
			str = ((gchar*) keys->data);
			if (str && str[0] != '\0') {
				fits_update_key(fit->fptr, TSTRING, keys->key, str, keys->comment, &status);
			}
			break;
		case KTYPE_DATE:
			status = 0;
			if (g_strcmp0("DATE", keys->key) == 0) {
				int itmp;
				char fit_date[40];
				fits_get_system_time(fit_date, &itmp, &status);
				fits_update_key(fit->fptr, TSTRING, keys->key, fit_date, keys->comment, &status);
			} else {
				date = *((GDateTime**) keys->data);
				if (date) {
					gchar *formatted_date = date_time_to_FITS_date(date);
					fits_update_key(fit->fptr, TSTRING, keys->key, formatted_date, keys->comment, &status);
					g_free(formatted_date);
				}
			}
			break;
		case KTYPE_BOOL:
			status = 0;
			sbool = *((super_bool*) keys->data);
			if (sbool != BOOL_NOT_SET) {
				bool = (sbool == BOOL_FALSE) ? FALSE : TRUE;
				fits_update_key(fit->fptr, TLOGICAL, keys->key, &(bool), keys->comment, &status);
			}
			break;
		default:
			siril_debug_print("Save_fits_keywords: Error. Type is not handled: %s.\n", keys->key);
		}
		keys++;
	}

	free(keys_start);

	return 0;
}

int remove_all_fits_keywords(fits *fit) {
	int status = 0, nkeys = 0;
	char keyname[FLEN_KEYWORD];
	char value[FLEN_VALUE];

	fits_get_hdrspace(fit->fptr, &nkeys, NULL, &status); /* get # of keywords */
	for (int i = nkeys; i > 0; i--) { // we start from the end because the keys are removed in place
		fits_read_keyn(fit->fptr, i, keyname, value, NULL, &status);
		siril_debug_print("%3d:%s=%s\n", i, keyname, value);
		if (!keyword_is_protected(keyname)) {
			fits_delete_record(fit->fptr, i, &status);
			siril_debug_print("%s removed\n", keyname);
		}
	}
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
		fits_update_key(fit->fptr, TSTRING, "OBJCTRA", &(fit->keywords.wcsdata.objctra), "[H M S] Image center Right Ascension", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "OBJCTDEC", &(fit->keywords.wcsdata.objctdec), "[D M S] Image center Declination", &status);
	}
	if (fit->keywords.wcsdata.ra > 0) {
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "RA", &(fit->keywords.wcsdata.ra), "[deg] Image center Right Ascension", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "DEC", &(fit->keywords.wcsdata.dec), "[deg] Image center Declination", &status);
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
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL1", &(fit->keywords.wcslib->crval[0]), "[deg] Axis1 reference value", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL2", &(fit->keywords.wcslib->crval[1]), "[deg] Axis2 reference value", &status);
		if (fit->keywords.wcslib->lonpole) {
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "LONPOLE", &(fit->keywords.wcslib->lonpole), "Native longitude of celestial pole", &status);
		}
		if (com.pref.wcs_formalism == WCS_FORMALISM_1) {
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CDELT1", &(fit->keywords.wcslib->cdelt[0]), "[deg] X pixel size", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CDELT2", &(fit->keywords.wcslib->cdelt[1]), "[deg] Y pixel size", &status);
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

static void set_to_default_not_used(fits *fit, GHashTable *keys_hash) {
	GHashTableIter iter;
	gpointer key, value;
	g_hash_table_iter_init(&iter, keys_hash);
	while (g_hash_table_iter_next(&iter, &key, &value)) {
		KeywordInfo *keyword_info = (KeywordInfo*) value;

		if (!keyword_info->used && should_use_keyword(fit, *keyword_info)) {
			switch (keyword_info->type) {
			case KTYPE_INT:
				if (keyword_info->data && *((int*) keyword_info->data) == 0)
					*((int*) keyword_info->data) = DEFAULT_INT_VALUE;
				break;
			case KTYPE_UINT:
				if (keyword_info->data && *((guint*) keyword_info->data) == 0)
					*((guint*) keyword_info->data) = DEFAULT_UINT_VALUE;
				break;
			case KTYPE_USHORT:
				if (keyword_info->data && *((gushort*) keyword_info->data) == 0)
					*((gushort*) keyword_info->data) = DEFAULT_USHORT_VALUE;
				break;
			case KTYPE_DOUBLE:
				if (keyword_info->data && *((double*) keyword_info->data) == 0.0)
					*((double*) keyword_info->data) = DEFAULT_DOUBLE_VALUE;
				break;
			case KTYPE_FLOAT:
				if (keyword_info->data && *((float*) keyword_info->data) == 0.f)
					*((float*) keyword_info->data) = DEFAULT_FLOAT_VALUE;
				break;
			default:
				break;
			}
		}
	}
}

void set_all_keywords_default(fits *fit) {
	GHashTable *keys_hash;
	KeywordInfo *keys = initialize_keywords(fit, &keys_hash);

	set_to_default_not_used(fit, keys_hash);

	/* Special cases */
	if (fit->keywords.stackcnt == DEFAULT_UINT_VALUE) {
		// DEFAULT_UINT_VALUE doesn't make any sense for a default stack count.
		// Change it to 1.
		fit->keywords.stackcnt = 1;
	}
	// Add any other special cases here...

	// Free the hash table and unknown keys
	g_hash_table_destroy(keys_hash);
	free(keys);
}

#ifdef DEBUG_PRINT_HEADER
const char* fits_type_to_string(const char type) {
	switch (type) {
		case 'C': return "string";
		case 'L': return "logical";
		case 'I': return "integer";
		case 'F': return "float";
		case 'X': return "complex";
		default: return "unknown";
	}
}
#endif

int read_fits_keywords(fits *fit) {
	// Initialize keywords and get hash table
	GHashTable *keys_hash;
	KeywordInfo *keys = initialize_keywords(fit, &keys_hash);
	int status = 0;
	int key_number = 1;
	GString *unknown_keys = g_string_new(NULL);

	read_fits_date_obs_header(fit); // handle very special case

	fits_get_hdrspace(fit->fptr, &key_number, NULL, &status); /* get # of keywords */

	// Loop through each keyword
#ifdef DEBUG_PRINT_HEADER
	printf("FITS Header:\n");
#endif
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

#ifdef DEBUG_PRINT_HEADER
		printf("%2d: Keyname: %8s, Type: %7s, Value: %30s, Comment: %s\n",
					ii, keyname, fits_type_to_string(type), value, comment);
#endif
		status = 0;

		// Skip the FITS structure keys and the DATE-OBS key since they have already been processed.
		if (g_strcmp0(keyname, "DATE-OBS") == 0 || keyword_is_protected(card)) {
			continue;
		}

		// Retrieve KeywordInfo from the hash table
		KeywordInfo *current_key = g_hash_table_lookup(keys_hash, keyname);

		// If the keyword is not found in the hash table, it is either an unknown or HISTORY keyword.
		// we don't want to load checksum keywords neither
		if (current_key == NULL) {
			GRegex *regex = g_regex_new("TR[0-9]+_[0-9]+|CROTA[0-9]", 0, 0, NULL);

			if (strncmp(card, "HISTORY", 7) == 0
					|| strncmp(card, "CHECKSUM", 8) == 0
					|| strncmp(card, "DATASUM", 7) == 0) {
				continue;
			}

			GMatchInfo *match_info = NULL;
			if (g_regex_match(regex, card, 0, &match_info)) {
				g_match_info_free(match_info);
				g_regex_unref(regex);
				continue;
			}

			if (match_info) {
				g_match_info_free(match_info);
			}

			g_regex_unref(regex);

			unknown_keys = g_string_append(unknown_keys, card);
			unknown_keys = g_string_append(unknown_keys, "\n");
			continue;
		}

		// At this point, the keyword is known and we can process it via the KeywordInfo list.

		if (current_key->fixed_value) {
			continue;
		}
		int int_value;
		guint uint_value;
		gushort ushort_value;
		double double_value;
		float float_value;
		gchar *str_value = NULL, *unquoted = NULL;
		GDateTime *date;
		super_bool bool_value;
		char *end = NULL;

		// Process the value based on the type of the keyword
		// Expected integer types may be included in the FITS file in scientific notation, e.g., 5.600E+01.
		// Hence we read integer values as double and then cast them to the respective integer type.
		switch (current_key->type) {
		case KTYPE_INT:
			double_value = g_ascii_strtod(value, &end);
			if (double_value < G_MININT || double_value > G_MAXINT) {
				siril_log_color_message("Warning: FITS value for keyname '%s' out of range for INT: %s\n",
						"salmon", keyname, value);
			}
			int_value = (int) double_value;
			if (value != end) {
				*((int*) current_key->data) = int_value;
				current_key->used = TRUE;
			} else {
				PRINT_PARSING_ERROR;
			}
			break;
		case KTYPE_UINT:
			double_value = g_ascii_strtod(value, &end);
			if (double_value < 0 || double_value > G_MAXUINT) {
				siril_log_color_message("Warning: FITS value for keyname '%s' out of range for UINT: %s\n",
						"salmon", keyname, value);
			}
			uint_value = (guint) double_value;
			if (value != end) {
				*((guint*) current_key->data) = uint_value;
				current_key->used = TRUE;
			} else {
				PRINT_PARSING_ERROR;
			}
			break;
		case KTYPE_USHORT:
			double_value = g_ascii_strtod(value, &end);
			if (double_value < 0 || double_value > G_MAXUSHORT) {
				siril_log_color_message("Warning: FITS value for keyname '%s' out of range for USHORT: %s\n",
						"salmon", keyname, value);
			}
			ushort_value = (gushort) double_value;
			if (value != end) {
				*((gushort*) current_key->data) = ushort_value;
				current_key->used = TRUE;
			} else {
				PRINT_PARSING_ERROR;
			}
			break;
		case KTYPE_DOUBLE:
			double_value = g_ascii_strtod(value, &end);
			if (value != end) {
				*((double*) current_key->data) = double_value;
				current_key->used = TRUE;
			} else {
				PRINT_PARSING_ERROR;
			}
			break;
		case KTYPE_FLOAT:
			float_value = g_ascii_strtod(value, &end);
			if (value != end) {
				*((float*) current_key->data) = float_value;
				current_key->used = TRUE;
			} else {
				PRINT_PARSING_ERROR;
			}
			break;
		case KTYPE_STR:
			unquoted = g_shell_unquote(value, NULL);
			str_value = g_strstrip(unquoted);
			strncpy((char*) current_key->data, str_value, FLEN_VALUE - 1);
			g_free(unquoted);
			current_key->used = TRUE;
			break;
		case KTYPE_DATE:
			unquoted = g_shell_unquote(value, NULL);
			str_value = g_strstrip(unquoted);
			date = FITS_date_to_date_time(str_value);
			if (date) {
				*((GDateTime**) current_key->data) = date;
				current_key->used = TRUE;
			} else {
				PRINT_PARSING_ERROR;
			}
			g_free(unquoted);
			break;
		case KTYPE_BOOL:
			bool_value = value[0] == 'T' ? BOOL_TRUE : BOOL_FALSE;
			*((super_bool*) current_key->data) = bool_value;
			current_key->used = TRUE;
			break;
		default:
			break;
		}

		// Handle special cases
		if (current_key->special_handler_read != NULL) {
			current_key->special_handler_read(fit, comment, current_key);
		}
	}

	gboolean not_from_siril = (strstr(fit->keywords.program, "Siril") == NULL);
	if ((fit->bitpix == FLOAT_IMG && not_from_siril) || fit->bitpix == DOUBLE_IMG) {
		float mini, maxi;
		fit_stats(fit->fptr, &mini, &maxi);
		// override data_max if needed. In some images there are differences between max and data_max
		fit->keywords.data_max = (double) maxi;
		fit->keywords.data_min = (double) mini;
	}

	if (fit->unknown_keys != NULL) {
		g_free(fit->unknown_keys);
	}
	fit->unknown_keys = g_string_free(unknown_keys, FALSE);

	set_to_default_not_used(fit, keys_hash);

	// Free the hash table and unknown keys
	g_hash_table_destroy(keys_hash);
	free(keys);
	return 0;
}


static void remove_keyword(const gchar *keyword, fits *fit, GHashTable *keys_hash) {
	GHashTableIter iter;
	gpointer key, value;
	g_hash_table_iter_init(&iter, keys_hash);
	while (g_hash_table_iter_next(&iter, &key, &value)) {
		KeywordInfo *keyword_info = (KeywordInfo*) value;
		if (keyword_info->data && !g_strcmp0(keyword_info->key, keyword)) {
			switch (keyword_info->type) {
			case KTYPE_BOOL:
				*((super_bool*) keyword_info->data) = BOOL_NOT_SET;
				break;
			case KTYPE_INT:
				*((int*) keyword_info->data) = DEFAULT_INT_VALUE;
				break;
			case KTYPE_UINT:
				*((guint*) keyword_info->data) = DEFAULT_UINT_VALUE;
				break;
			case KTYPE_USHORT:
				*((gushort*) keyword_info->data) = DEFAULT_USHORT_VALUE;
				break;
			case KTYPE_DOUBLE:
				*((double*) keyword_info->data) = DEFAULT_DOUBLE_VALUE;
				break;
			case KTYPE_FLOAT:
				*((float*) keyword_info->data) = DEFAULT_FLOAT_VALUE;
				break;
			case KTYPE_STR:
				memset((char*) keyword_info->data, 0, FLEN_VALUE);
				break;
			case KTYPE_DATE:
				if ((GDateTime*) keyword_info->data) {
					g_date_time_unref((GDateTime*) keyword_info->data);
					keyword_info->data = NULL;
				}
				break;
			default:
				break;
			}
		}
	}
}

void remove_keyword_in_fit_keywords(const gchar *keyword, fits *fit) {
	GHashTable *keys_hash;
	KeywordInfo *keys = initialize_keywords(fit, &keys_hash);
	remove_keyword(keyword, fit, keys_hash);

	g_hash_table_destroy(keys_hash);
	free(keys);
}

static int keywords_prepare_hook(struct generic_seq_args *arg) {
	if (arg->seq->type == SEQ_FITSEQ) {
		gchar *filename = g_strdup(arg->seq->fitseq_file->filename);
		// it was opened in READONLY mode, we close it
		if (fitseq_close_file(arg->seq->fitseq_file)) {
			siril_log_color_message(_("Error when closing fitseq\n"), "red");
			g_free(filename);
			return 1;
		}
		arg->seq->fitseq_file->fptr = NULL;
		arg->seq->fitseq_file->filename = g_strdup(filename); // freed in fitseq_destroy
		// and we reopen in READWRITE mode to update it
		if (fitseq_open(filename, arg->seq->fitseq_file, READWRITE)) {
			siril_log_color_message(_("Error when reopening fitseq\n"), "red");
			g_free(filename);
			return 1;
		}
		g_free(filename);
	}
	return 0;
}

static int keywords_image_hook(struct generic_seq_args *arg, int o, int i, fits *fit, rectangle *area, int threads) {
	struct keywords_data *kargs = (struct keywords_data *)arg->user;

	int retval = updateFITSKeyword(fit, kargs->FITS_key, kargs->newkey, kargs->value, kargs->comment, FALSE, arg->seq->type == SEQ_FITSEQ);
	if (arg->seq->type == SEQ_REGULAR) {
		char root[256];
		if (!retval) {
			if (!fit_sequence_get_image_filename(arg->seq, i, root, FALSE)) {
				return 1;
			}
			savefits(root, fit);
		}
	} else if (arg->seq->type == SEQ_FITSEQ) { // case SEQ_FITSEQ, fit already holds its fptr, we just update
		save_fits_header(fit);
		siril_log_color_message(_("FITS header of image %d updated\n"), "salmon", i + 1);
	}
	return 0;
}

static int keywords_finalize_hook(struct generic_seq_args *arg) {
	struct keywords_data *kargs = (struct keywords_data *)arg->user;
	int retval = 0;
	if (arg->seq->type == SEQ_FITSEQ) {
		gchar *filename = g_strdup(arg->seq->fitseq_file->filename);
		// it was opened in READWRITE mode, we close it to save everything
		if (fitseq_close_file(arg->seq->fitseq_file)) {
			siril_debug_print("error when closing again fitseq\n");
			g_free(filename);
			retval = 1;
			goto finish;
		}
		arg->seq->fitseq_file->fptr = NULL;
		siril_log_color_message(_("File %s updated\n"), "salmon", filename);
		arg->seq->fitseq_file->filename = filename; // we may need to reopen in the idle so we save it here
		arg->seq->fitseq_file->hdu_index = NULL;
	}
	if (!arg->retval)
		writeseqfile(arg->seq);
finish:
	free(kargs->FITS_key);
	free(kargs->value);
	free(kargs->comment);
	free(kargs);
	return retval;
}

gboolean end_keywords_sequence(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	if (check_seq_is_comseq(args->seq)) {
		if (args->seq->type == SEQ_FITSEQ) { // if FITSEQ, we need to repoen in READONLY mode
			if (fitseq_open(args->seq->fitseq_file->filename, args->seq->fitseq_file, READONLY)) {
				siril_debug_print("error when finally re-opening fitseq\n");
			}
		}
		update_sequences_list(args->seq->seqname);
		gui_function(refresh_keywords_dialog, NULL);
	}
	if (!check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);
	free(p);
	return end_generic(NULL);
}

void start_sequence_keywords(sequence *seq, struct keywords_data *args) {
	struct generic_seq_args *seqargs = create_default_seqargs(seq);
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seq->selnum;
	seqargs->stop_on_error = FALSE;
	seqargs->parallel = seq->type != SEQ_FITSEQ;
	seqargs->prepare_hook = keywords_prepare_hook;
	seqargs->image_hook = keywords_image_hook;
	seqargs->finalize_hook =  keywords_finalize_hook;
	seqargs->idle_function = end_keywords_sequence;
	seqargs->has_output = (seq->type == SEQ_SER); // we don't update keywords for SER file
	seqargs->output_type = get_data_type(seq->bitpix);
	seqargs->description = "keywords update";
	if (seq->type == SEQ_SER) {
		siril_log_color_message(_("This command won't work for SER sequence.\n"), "red");
		free(seqargs);
		return;
	}
	seqargs->user = args;
	if (!start_in_new_thread(generic_sequence_worker, seqargs)) {
		free_generic_seq_args(seqargs, TRUE);
	}
}


int parse_wcs_image_dimensions(fits *fit, int *rx, int *ry) {
	*rx = 0;
	*ry = 0;
	if (!fit)
		return 1;
	if (strlen(fit->unknown_keys) > 0) {
		gchar **headerkeys = g_strsplit(fit->unknown_keys, "\n", 0);
		double drx = 0., dry = 0.;
		int status_w = read_key_from_header_text(headerkeys, "IMAGEW", &drx, NULL);
		int status_h = read_key_from_header_text(headerkeys, "IMAGEH", &dry, NULL);
		g_strfreev(headerkeys);
		if (status_w || status_h)
			return 1;
		*rx = (int)drx;
		*ry = (int)dry;
		siril_debug_print("IMAGEW: %d, IMAGEH: %d\n", *rx, *ry);
		return 0;
	}
	return 1;
}

void clear_Bayer_information(fits *fit) {
	memset(fit->keywords.bayer_pattern, 0, FLEN_VALUE);
}