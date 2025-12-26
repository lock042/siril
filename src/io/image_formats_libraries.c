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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lcms2.h>

#ifdef HAVE_LIBTIFF
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#endif
#ifdef HAVE_LIBJPEG
#include <jpeglib.h>
#include <jconfig.h>
#include <jerror.h>
#endif
#ifdef HAVE_LIBPNG
#include <png.h>
#include <setjmp.h>
#endif
#ifdef HAVE_LIBRAW
#include <libraw/libraw.h>
#include <libraw/libraw_version.h>
#endif
#ifdef HAVE_LIBHEIF
#include <libheif/heif.h>
#endif
#ifdef HAVE_LIBXISF
#include "io/SirilXISFWraper.h"
#endif
#ifdef HAVE_LIBJXL
#include "io/SirilJpegXLWrapper.h"
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "core/siril_log.h"
#include "core/exif.h"
#include "io/fits_keywords.h"
#include "algos/geometry.h"
#include "algos/demosaicing.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "image_format_fits.h"

static void fill_date_obs_if_any(fits *fit, const char *file) {
	gchar *date_time = NULL;
	int year = 1, month = 1, day = 1, h = 0, m = 0, s = 0;
	date_time = siril_get_date_from_exif(file);
	if (date_time) {
		int n = sscanf(date_time, "%04d:%02d:%02d %02d:%02d:%02d", &year, &month, &day, &h, &m, &s);
		if (n == 6) {
			GTimeZone *tz = g_time_zone_new_utc();
			fit->keywords.date_obs = g_date_time_new(tz, year, month, day, h, m, (double) s);
			g_time_zone_unref(tz);
		}
		g_free(date_time);
	}
}

/********************* TIFF IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBTIFF

static int readtifstrip(TIFF* tif, uint32_t width, uint32_t height, uint16_t nsamples, uint16_t color, WORD **data) {
	uint32_t rowsperstrip;
	uint16_t config;
	int retval = nsamples;

	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &config);
	TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

	size_t npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	if (!*data) {
		PRINT_ALLOC_ERR;
		return OPEN_IMAGE_ERROR;
	}
	WORD *gbuf[3] = {*data, *data, *data};
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[GLAYER] = *data + npixels;
		gbuf[BLAYER] = *data + npixels * 2;
	}

	const tmsize_t scanline = TIFFScanlineSize(tif);
	WORD *buf = (WORD *)_TIFFmalloc(TIFFStripSize(tif));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return OPEN_IMAGE_ERROR;
	}
	for (uint32_t row = 0; row < height; row += rowsperstrip){
		uint32_t nrow = (row + rowsperstrip > height ? height - row : rowsperstrip);
		switch (config) {
		case PLANARCONFIG_CONTIG:
			if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, 0), buf, nrow * scanline) < 0) {
				siril_log_color_message(_("An unexpected error was encountered while trying to read the file.\n"), "red");
				retval = OPEN_IMAGE_ERROR;
				break;
			}
			for (size_t i = 0; i < width * nrow; i++) {
				if (color == PHOTOMETRIC_MINISWHITE) {
					*gbuf[RLAYER]++ = USHRT_MAX - buf[i * nsamples + 0];
				} else {
					*gbuf[RLAYER]++ = buf[i * nsamples + 0];
				}

				if ((nsamples == 3) || (nsamples == 4)) {
					if (color == PHOTOMETRIC_MINISWHITE) {
						*gbuf[GLAYER]++ = USHRT_MAX - buf[i * nsamples + 1];
						*gbuf[BLAYER]++ = USHRT_MAX - buf[i * nsamples + 2];
					} else {
						*gbuf[GLAYER]++ = buf[i * nsamples + 1];
						*gbuf[BLAYER]++ = buf[i * nsamples + 2];
					}
				}
			}
			break;
		case PLANARCONFIG_SEPARATE:
			if (nsamples >= 3)		//don't need to read the alpha
				nsamples = 3;
			for (int j = 0; j < nsamples; j++) {	//loop on the layer
				if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, j), buf, nrow * scanline) < 0) {
					siril_log_color_message(_("An unexpected error was encountered while trying to read the file.\n"), "red");
					retval = OPEN_IMAGE_ERROR;
					break;
				}
				for (size_t i = 0; i < width * nrow; i++) {
					if (color == PHOTOMETRIC_MINISWHITE) {
						*gbuf[j]++ = USHRT_MAX - buf[i];
					} else {
						*gbuf[j]++ = buf[i];
					}
				}
			}
			break;
		default:
			siril_log_color_message(_("Unknown TIFF file.\n"), "red");
			retval = OPEN_IMAGE_ERROR;
		}
	}
	_TIFFfree(buf);
	return retval;
}

static int readtifstrip32(TIFF* tif, uint32_t width, uint32_t height, uint16_t nsamples, uint16_t color, float **data) {
	uint32_t rowsperstrip;
	uint16_t config;
	int retval = nsamples;

	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &config);
	TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

	size_t npixels = width * height;
	*data = malloc(npixels * sizeof(float) * nsamples);
	if (!*data) {
		PRINT_ALLOC_ERR;
		return OPEN_IMAGE_ERROR;
	}
	float *gbuf[3] = { *data, *data, *data };
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	const tmsize_t scanline = TIFFScanlineSize(tif);
	float *buf = (float *)_TIFFmalloc(TIFFStripSize(tif));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return OPEN_IMAGE_ERROR;
	}
	for (uint32_t row = 0; row < height; row += rowsperstrip) {
		uint32_t nrow = (row + rowsperstrip > height ? height - row : rowsperstrip);
		switch (config) {
		case PLANARCONFIG_CONTIG:
			if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, 0), buf, nrow * scanline) < 0) {
				siril_log_color_message(_("An unexpected error was encountered while trying to read the file.\n"), "red");
				retval = OPEN_IMAGE_ERROR;
				break;
			}
			for (size_t i = 0; i < width * nrow; i++) {
				if (color == PHOTOMETRIC_MINISWHITE) {
					*gbuf[RLAYER]++ = USHRT_MAX_SINGLE - buf[i * nsamples + 0];
				} else {
					*gbuf[RLAYER]++ = buf[i * nsamples + 0];
				}
				if ((nsamples == 3) || (nsamples == 4)) {
					if (color == PHOTOMETRIC_MINISWHITE) {
						*gbuf[GLAYER]++ = USHRT_MAX_SINGLE - buf[i * nsamples + 1];
						*gbuf[BLAYER]++ = USHRT_MAX_SINGLE - buf[i * nsamples + 2];
					} else {
						*gbuf[GLAYER]++ = buf[i * nsamples + 1];
						*gbuf[BLAYER]++ = buf[i * nsamples + 2];
					}
				}
			}
			break;
		case PLANARCONFIG_SEPARATE:
			if (nsamples >= 3)		//don't need to read the alpha
				nsamples = 3;
			for (int j = 0; j < nsamples; j++) {	//loop on the layer
				if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, j),
						buf, nrow * scanline) < 0) {
					siril_log_color_message(_("An unexpected error was encountered while trying to read the file.\n"), "red");
					retval = OPEN_IMAGE_ERROR;
					break;
				}
				for (size_t i = 0; i < width * nrow; i++) {
					if (color == PHOTOMETRIC_MINISWHITE) {
						*gbuf[j]++ = USHRT_MAX_SINGLE - buf[i];
					} else {
						*gbuf[j]++ = buf[i];
					}
				}
			}
			break;
		default:
			siril_log_color_message(_("Unknown TIFF file.\n"), "red");
			retval = OPEN_IMAGE_ERROR;
		}
	}
	_TIFFfree(buf);
	return retval;
}

static int readtifstrip32uint(TIFF* tif, uint32_t width, uint32_t height, uint16_t nsamples, uint16_t color, float **data) {
	uint32_t rowsperstrip;
	uint16_t config;
	int retval = nsamples;

	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &config);
	TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

	size_t npixels = width * height;
	*data = malloc(npixels * sizeof(float) * nsamples);
	if (!*data) {
		PRINT_ALLOC_ERR;
		return OPEN_IMAGE_ERROR;
	}
	float *gbuf[3] = { *data, *data, *data };
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	const tmsize_t scanline = TIFFScanlineSize(tif);
	uint32_t *buf = (uint32_t *)_TIFFmalloc(TIFFStripSize(tif));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return OPEN_IMAGE_ERROR;
	}
	for (uint32_t row = 0; row < height; row += rowsperstrip) {
		uint32_t nrow = (row + rowsperstrip > height ? height - row : rowsperstrip);
		switch (config) {
		case PLANARCONFIG_CONTIG:
			if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, 0), buf, nrow * scanline) < 0) {
				siril_log_color_message(_("An unexpected error was encountered while trying to read the file.\n"), "red");
				retval = OPEN_IMAGE_ERROR;
				break;
			}
			for (size_t i = 0; i < width * nrow; i++) {
				if (color == PHOTOMETRIC_MINISWHITE) {
					*gbuf[RLAYER]++ = USHRT_MAX_SINGLE - (buf[i * nsamples + 0] / (float) UINT32_MAX);
				} else {
					*gbuf[RLAYER]++ = buf[i * nsamples + 0] / (float) UINT32_MAX;
				}
				if ((nsamples == 3) || (nsamples == 4)) {
					if (color == PHOTOMETRIC_MINISWHITE) {
						*gbuf[GLAYER]++ = USHRT_MAX_SINGLE - (buf[i * nsamples + 1] / (float) UINT32_MAX);
						*gbuf[BLAYER]++ = USHRT_MAX_SINGLE - (buf[i * nsamples + 2] / (float) UINT32_MAX);
					} else {
						*gbuf[GLAYER]++ = buf[i * nsamples + 1] / (float) UINT32_MAX;
						*gbuf[BLAYER]++ = buf[i * nsamples + 2] / (float) UINT32_MAX;
					}
				}
			}
			break;
		case PLANARCONFIG_SEPARATE:
			if (nsamples >= 3)		//don't need to read the alpha
				nsamples = 3;
			for (int j = 0; j < nsamples; j++) {	//loop on the layer
				if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, j),
						buf, nrow * scanline) < 0) {
					siril_log_color_message(_("An unexpected error was encountered while trying to read the file.\n"), "red");
					retval = OPEN_IMAGE_ERROR;
					break;
				}
				for (size_t i = 0; i < width * nrow; i++) {
					if (color == PHOTOMETRIC_MINISWHITE) {
						*gbuf[j]++ = USHRT_MAX_SINGLE - (buf[i] / (float) UINT32_MAX);
					} else {
						*gbuf[j]++ = buf[i] / (float) UINT32_MAX;
					}
				}
			}
			break;
		default:
			siril_log_color_message(_("Unknown TIFF file.\n"), "red");
			retval = OPEN_IMAGE_ERROR;
		}
	}
	_TIFFfree(buf);
	return retval;
}

static int readtif8bits(TIFF* tif, uint32_t width, uint32_t height, uint16_t nsamples, uint16_t color, WORD **data) {
	int retval = nsamples;

	size_t npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	if (!*data) {
		PRINT_ALLOC_ERR;
		return OPEN_IMAGE_ERROR;
	}
	WORD *gbuf[3] = { *data, *data, *data };
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	/* get the data */
	uint32_t *raster = (uint32_t*) _TIFFmalloc(npixels * sizeof(uint32_t));
	if (raster != NULL) {
		if (TIFFReadRGBAImage(tif, width, height, raster, 0)) {
			for (int j = 0; j < height; j++) {
				int istart = j * width;
				for (int i = 0; i < width; i++) {
					if (color == PHOTOMETRIC_MINISWHITE) {
						*gbuf[RLAYER]++ = UCHAR_MAX - (WORD)TIFFGetR(raster[istart + i]);
					} else {
						*gbuf[RLAYER]++ = (WORD)TIFFGetR(raster[istart + i]);
					}
					if ((nsamples == 3) || (nsamples == 4)) {
						if (color == PHOTOMETRIC_MINISWHITE) {
							*gbuf[GLAYER]++ = UCHAR_MAX - (WORD)TIFFGetG(raster[istart + i]);
							*gbuf[BLAYER]++ = UCHAR_MAX - (WORD)TIFFGetB(raster[istart + i]);
						} else {
							*gbuf[GLAYER]++ = (WORD)TIFFGetG(raster[istart + i]);
							*gbuf[BLAYER]++ = (WORD)TIFFGetB(raster[istart + i]);
						}
					}
				}
			}
		}
		else {
			siril_log_color_message(_("An unexpected error was encountered while trying to read the file.\n"), "red");
			retval = OPEN_IMAGE_ERROR;
		}
		_TIFFfree(raster);
	}
	else retval = OPEN_IMAGE_ERROR;
	return retval;
}

gboolean get_tiff_compression() {
	if (!com.headless) {
		GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("radiobuttonCompDeflate"));
		if (gtk_toggle_button_get_active(button))
			return TRUE;
	}
	return FALSE;
}

static TIFF* Siril_TIFFOpen(const char *name, const char *mode) {
#ifdef _WIN32
	wchar_t *wname;

	wname = g_utf8_to_utf16(name, -1, NULL, NULL, NULL);
	if (wname == NULL) {
		return NULL;
	}

	TIFF* tif = TIFFOpenW(wname, mode);
	g_free(wname);
	return tif;
#else
	return(TIFFOpen(name, mode));
#endif
}

/* reads a TIFF file and stores it in the fits argument.
 * If file loading fails, the argument is untouched.
 */
int readtif(const char *name, fits *fit, gboolean force_float, gboolean verbose) {
	int retval = 0;
	uint32_t height, width;
	uint16_t nbits, nsamples, color, orientation;
	WORD *data = NULL;
	float *fdata = NULL;
	uint16_t sampleformat = 0;
	/* EXIFS */
	gchar *description = NULL;
	double exposure = 0.0;
	double aperture = 0.0;
	double focal_length = 0.0;

	TIFF* tif = Siril_TIFFOpen(name, "r");
	if (!tif) {
		siril_log_message(_("Could not open the TIFF file %s\n"), name);
		return OPEN_IMAGE_ERROR;
	}

	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH, &height);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &sampleformat);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &nbits);
	TIFFGetFieldDefaulted(tif, TIFFTAG_PHOTOMETRIC, &color);
	TIFFGetFieldDefaulted(tif, TIFFTAG_ORIENTATION, &orientation);

	// Retrieve the Date/Time as in the TIFF TAG
	gchar *date_time = NULL;
	int year = 1, month = 1, day = 1, h = 0, m = 0, s = 0;

	if (TIFFGetField(tif, TIFFTAG_DATETIME, &date_time)) {
		sscanf(date_time, "%04d:%02d:%02d %02d:%02d:%02d", &year, &month, &day, &h, &m, &s);
	}

	// Try to read EXIF data
	toff_t ExifID;

	if (TIFFGetField(tif, TIFFTAG_EXIFIFD, &ExifID)) {
		uint16_t currentIFD = TIFFCurrentDirectory(tif);

		if (TIFFReadEXIFDirectory(tif, ExifID)) {
			TIFFGetField(tif, EXIFTAG_EXPOSURETIME, &exposure);
			TIFFGetField(tif, EXIFTAG_FNUMBER, &aperture);
			TIFFGetField(tif, EXIFTAG_FOCALLENGTH, &focal_length);
		}
		// Revert IFD to status quo ante TIFFReadEXIFDirectory
		TIFFSetDirectory(tif, currentIFD);
	}

	// Try to read embedded ICC profile data
	cmsUInt32Number EmbedLen = 0;
	cmsUInt8Number* EmbedBuffer = NULL;
	(void) TIFFGetField(tif, TIFFTAG_ICCPROFILE, &EmbedLen, &EmbedBuffer);

	// Retrieve Description field
	char *desc = NULL;
	if (TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &desc)) {
		description = g_strdup(desc);
	}

	size_t npixels = width * height;

	switch(nbits){
		case 8:
			/* High level functions in readtif8bits: should read every 8-bit TIFF file */
			retval = readtif8bits(tif, width, height, nsamples, color, &data);
			break;

		case 16:
			retval = readtifstrip(tif, width, height, nsamples, color, &data);
			break;

		case 32:
			if (sampleformat == SAMPLEFORMAT_IEEEFP) {
				retval = readtifstrip32(tif, width, height, nsamples, color, &fdata);
			} else if (sampleformat == SAMPLEFORMAT_UINT) {
				retval = readtifstrip32uint(tif, width, height, nsamples, color, &fdata);
			} else {
				siril_log_color_message(_("Siril cannot read this TIFF format.\n"), "red");
				retval = OPEN_IMAGE_ERROR;
			}
			break;

		default :
			siril_log_color_message(_("Siril cannot read this TIFF format.\n"), "red");
			retval = OPEN_IMAGE_ERROR;
	}

	/* We clear fits. Everything written above is erased */
	/* note: this has been moved slightly, it has to happen before we initialize the ICC profile
	 * which in turn has to happen before we close the TIFF */
	clearfits(fit);

	set_all_keywords_default(fit);

	if (date_time) {
		GTimeZone *tz = g_time_zone_new_utc();
		fit->keywords.date_obs = g_date_time_new(tz, year, month, day, h, m, (double) s);
		g_time_zone_unref(tz);
	}
	cmsUInt8Number *embed = NULL;
	if (EmbedLen) {
		embed = malloc(EmbedLen * sizeof(cmsUInt8Number));
		memcpy(embed, EmbedBuffer, EmbedLen * sizeof(cmsUInt8Number));
	}
	cmsUInt32Number len = EmbedLen;

	TIFFClose(tif);
	if (retval < 0) {
		free(data);
		free(fdata);
		g_free(description);
		free(embed);
		return OPEN_IMAGE_ERROR;
	}

	fit->rx = width;
	fit->ry = height;
	fit->naxes[0] = width;
	fit->naxes[1] = height;
	fit->data = data;
	fit->fdata = fdata;
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	if (nsamples == 1 || nsamples == 2) {
		fit->naxes[2] = 1;
		fit->naxis = 2;
		if (data) {
			fit->pdata[RLAYER] = fit->data;
			fit->pdata[GLAYER] = fit->data;
			fit->pdata[BLAYER] = fit->data;
		} else {
			fit->fpdata[RLAYER] = fit->fdata;
			fit->fpdata[GLAYER] = fit->fdata;
			fit->fpdata[BLAYER] = fit->fdata;
		}
	} else {
		fit->naxes[2] = 3;
		fit->naxis = 3;
		if (data) {
			fit->pdata[RLAYER] = fit->data;
			fit->pdata[GLAYER] = fit->data + npixels;
			fit->pdata[BLAYER] = fit->data + npixels * 2;
		} else {
			fit->fpdata[RLAYER] = fit->fdata;
			fit->fpdata[GLAYER] = fit->fdata + npixels;
			fit->fpdata[BLAYER] = fit->fdata + npixels * 2;
		}
	}
	switch (nbits) {
	case 8:
		fit->bitpix = BYTE_IMG;
		fit->type = DATA_USHORT;
		if (force_float) {
			size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
			fit_replace_buffer(fit, ushort8_buffer_to_float(fit->data, ndata), DATA_FLOAT);
		}
		break;
	case 16:
		fit->bitpix = USHORT_IMG;
		fit->type = DATA_USHORT;
		if (force_float) {
			size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
			fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
		}

		if (orientation == ORIENTATION_TOPLEFT) {
			mirrorx(fit, FALSE);
		} else if (orientation == ORIENTATION_TOPRIGHT) {
			mirrorx(fit, FALSE);
			mirrory(fit, FALSE);
		} else if (orientation == ORIENTATION_BOTRIGHT) {
			mirrory(fit, FALSE);
		} else if (orientation == ORIENTATION_BOTLEFT) {
			; // do nothing
		} else {
			siril_debug_print(_("TIFFTAG Orientation not handled.\n"));
		}
		break;
	case 32:
		fit->bitpix = FLOAT_IMG;
		fit->type = DATA_FLOAT;
		if (orientation == ORIENTATION_TOPLEFT) {
			mirrorx(fit, FALSE);
		} else if (orientation == ORIENTATION_TOPRIGHT) {
			mirrorx(fit, FALSE);
			mirrory(fit, FALSE);
		} else if (orientation == ORIENTATION_BOTRIGHT) {
			mirrory(fit, FALSE);
		} else if (orientation == ORIENTATION_BOTLEFT) {
			; // do nothing
		} else {
			siril_debug_print(_("TIFFTAG Orientation not handled.\n"));
		}
	}
	fit->orig_bitpix = fit->bitpix;
	g_snprintf(fit->keywords.row_order, FLEN_VALUE, "%s", "TOP-DOWN");

	/* fill exifs is exist */
	if (exposure > 0.0)
		fit->keywords.exposure = exposure;
	if (aperture > 0.0)
		fit->keywords.aperture = aperture;
	if (focal_length > 0.0) {
		fit->keywords.focal_length = focal_length;
		fit->focalkey = TRUE;
	}
	if (description) {
		if (g_str_has_prefix(description, "SIMPLE  =")) {
			// It is FITS header, copy it
			siril_debug_print("ASTRO-TIFF detected.\n");
			if (fit->header) free(fit->header);
			fit->header = description;
			int ret = fits_parse_header_str(fit, description);
			if (ret) {
				siril_debug_print("ASTRO-TIFF is not well formed.\n");
			}
		} else {
			free(description);
		}
	}

	fits_initialize_icc(fit, embed, len);
	free(embed);

	retval = nsamples;

	gchar *basename = g_path_get_basename(name);
	if (verbose)
		siril_log_message(_("Reading TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
						nbits, basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return retval;
}

void get_tif_data_from_ui(fits *fit, gchar **description, gchar **copyright) {
	if (!com.script && !com.headless) {
		/*******************************************************************
		 * If the user saves a tif from the graphical menu, he can set
		 * the Description and the Copyright of the Image
		 ******************************************************************/

		GtkTextIter itDebut;
		GtkTextIter itFin;

		GtkTextView *description_txt_view = GTK_TEXT_VIEW(lookup_widget("Description_txt"));
		GtkTextBuffer *desbuf = gtk_text_view_get_buffer(description_txt_view);
		gtk_text_buffer_get_start_iter(desbuf, &itDebut);
		gtk_text_buffer_get_end_iter(desbuf, &itFin);
		*description = gtk_text_buffer_get_text(desbuf, &itDebut, &itFin, TRUE);

		GtkTextView *copyright_txt_view = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));
		GtkTextBuffer *copybuf = gtk_text_view_get_buffer(copyright_txt_view);
		gtk_text_buffer_get_start_iter(copybuf, &itDebut);
		gtk_text_buffer_get_end_iter(copybuf, &itFin);
		*copyright = gtk_text_buffer_get_text(copybuf, &itDebut, &itFin, TRUE);

	}
}

/*** This function save the current image into a uncompressed 8- or 16-bit file *************/
int savetif(const char *name, fits *fit, uint16_t bitspersample,
		const gchar *description, const gchar *copyright,
		gboolean tiff_compression, gboolean embeded_icc, gboolean verbose) {
	int retval = 0;
	float norm;
	gchar *filename = g_strdup(name);
	uint32_t profile_len = 0;
	unsigned char *profile = NULL;
	gboolean write_ok = TRUE;

	if (!g_str_has_suffix(filename, ".tif") && (!g_str_has_suffix(filename, ".tiff"))) {
		filename = str_append(&filename, ".tif");
	}

	TIFF* tif = Siril_TIFFOpen(filename, "w");
	if (!tif) {
		siril_log_color_message(_("Siril cannot create TIFF file.\n"), "red");
		free(filename);
		return 1;
	}
	const uint16_t nsamples = (uint16_t) fit->naxes[2];
	const uint32_t width = (uint32_t) fit->rx;
	const uint32_t height = (uint32_t) fit->ry;


	/*******************************************************************/

	/* TIFF TAG FIELD */
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, bitspersample == 32 ? SAMPLEFORMAT_IEEEFP : SAMPLEFORMAT_UINT);
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, -1));
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, nsamples);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, tiff_compression ? COMPRESSION_ADOBE_DEFLATE : COMPRESSION_NONE);
	if (description) {
		TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, description);
	}
	if (copyright) {
		TIFFSetField(tif, TIFFTAG_COPYRIGHT, copyright);
	}
	TIFFSetField(tif, TIFFTAG_SOFTWARE, PACKAGE " v" VERSION);

	gboolean src_is_float = (fit->type == DATA_FLOAT);
	if (nsamples == 1) {
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	} else if (nsamples == 3) {
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	} else {
		TIFFClose(tif);
		siril_log_color_message(_("TIFF file has unexpected number of channels (not 1 or 3).\n"), "red");
		free(filename);
		return 1;
	}

	if (fit->keywords.date_obs) {
		gchar *date_time = g_date_time_format(fit->keywords.date_obs, "%Y:%m:%d %H:%M:%S");

		TIFFSetField(tif, TIFFTAG_DATETIME, date_time);
		g_free(date_time);
	}

	WORD *gbuf[3] =	{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	float *gbuff[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	// Color management (only if the image is color managed and has a profile)
	void *buf = NULL;
	void *dest = NULL;
	cmsHTRANSFORM save_transform = NULL;
	gboolean threaded = !get_thread_run();
	if (fit->color_managed && fit->icc_profile) {
		// Transform the data
		buf = src_is_float ? (void *) fit->fdata : (void *) fit->data;
		size_t npixels = fit->rx * fit->ry;
		cmsColorSpaceSignature sig = cmsGetColorSpace(fit->icc_profile);
		cmsUInt32Number trans_type = get_planar_formatter_type(sig, fit->type, FALSE);
		if (src_is_float) {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(float));
		} else {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(WORD));
		}
		// Check what is the appropriate color space to save in
		// Covers 8- and high-bitdepth files
		// bitspersample is *destination* bits per sample (because that affects the export
		// ICC profile) but transforms are done source bit depth to source bit depth
		if (bitspersample == 8) { // 8-bit
			if (nsamples == 1) { // mono
				cmsHPROFILE srgb_mono_out = NULL;
				switch (com.pref.icc.export_8bit_method) {
					case EXPORT_SRGB:
						srgb_mono_out = gray_srgbtrc();
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(srgb_mono_out, &profile_len);
						cmsCloseProfile(srgb_mono_out);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.mono_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(dest);
						free(filename);
						return 1;
				}
			} else { // rgb
				switch (com.pref.icc.export_8bit_method) {
					case EXPORT_SRGB:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.srgb_out, &profile_len);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.working_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.working_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(dest);
						free(filename);
						return 1;
				}
			}
		} else  if (bitspersample == 16) { // 16-bit or 32-bit save from < 32-bit input
			if (nsamples == 1) { // mono
				cmsHPROFILE srgb_mono_out = NULL;
				switch (com.pref.icc.export_16bit_method) {
					case EXPORT_SRGB:
						srgb_mono_out = gray_srgbtrc();
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(srgb_mono_out, &profile_len);
						cmsCloseProfile(srgb_mono_out);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.mono_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(dest);
						free(filename);
						return 1;
				}
			} else { // rgb
				switch (com.pref.icc.export_16bit_method) {
					case EXPORT_SRGB:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.srgb_out, &profile_len);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.working_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.working_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(dest);
						free(filename);
						return 1;
				}
			}
		} else {
			// 32-bit images are always saved with the image profile
			profile = get_icc_profile_data(fit->icc_profile, &profile_len);
		}
		cmsUInt32Number datasize = fit->type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
		cmsUInt32Number bytesperline = width * datasize;
		cmsUInt32Number bytesperplane = npixels * datasize;
		if (save_transform) {
			// Do the transform
			cmsDoTransformLineStride(save_transform, buf, dest, width, height, bytesperline, bytesperline, bytesperplane, bytesperplane);
			cmsDeleteTransform(save_transform);
		} else {
			// For "use image ICC profile" save_transform will be NULL, no need to transform the data
			memcpy(dest, buf, bytesperplane * nsamples);
		}
		gbuf[0] = (WORD *) dest;
		gbuf[1] = (WORD *) dest + (fit->rx * fit->ry);
		gbuf[2] = (WORD *) dest + (fit->rx * fit->ry * 2);
		gbuff[0] = (float *) dest;
		gbuff[1] = (float *) dest + (fit->rx * fit->ry);
		gbuff[2] = (float *) dest + (fit->rx * fit->ry * 2);
	} else if (com.pref.icc.default_to_srgb) { // Default for non-color managed files is to assume sRGB but there is a preference for this
		profile = get_icc_profile_data(fit->naxes[2] == 1 ? com.icc.mono_out : com.icc.srgb_out, &profile_len);
	}
	// If there is an ICC profile, embed it in the file.
	if (profile) {
		TIFFSetField(tif, TIFFTAG_ICCPROFILE, profile_len, profile);
	}

	switch (bitspersample) {
	case 8:
		siril_debug_print("Saving 8-bit TIFF file.\n");
		BYTE *buf8 = _TIFFmalloc(width * sizeof(unsigned char) * nsamples);
		if (!buf8) {
			PRINT_ALLOC_ERR;
			retval = OPEN_IMAGE_ERROR;
			write_ok = FALSE;
			break;
		}

		norm = fit->orig_bitpix != BYTE_IMG ? UCHAR_MAX_SINGLE / USHRT_MAX_SINGLE : 1.f;

		for (uint32_t row = height; row-- > 0;) {
			for (uint32_t col = 0; col < width; col++) {
				for (uint16_t n = 0; n < nsamples; n++) {
					buf8[col * nsamples + n] =
							(fit->type == DATA_USHORT) ?
									gbuf[n][col + row * width] * norm :
									float_to_uchar_range(gbuff[n][col + row * width]);
				}
			}
			if (TIFFWriteScanline(tif, buf8, height - 1 - row, 0) < 0) {
				siril_debug_print("Error while writing in TIFF File.\n");
				retval = OPEN_IMAGE_ERROR;
				write_ok = FALSE;
				break;
			}
		}
		_TIFFfree(buf8);
		break;
	case 16:
		siril_debug_print("Saving 16-bit TIFF file.\n");
		WORD *buf16 = _TIFFmalloc(width * sizeof(WORD) * nsamples);
		if (!buf16) {
			PRINT_ALLOC_ERR;
			retval = OPEN_IMAGE_ERROR;
			write_ok = FALSE;
			break;
		}

		norm = fit->orig_bitpix == BYTE_IMG ? USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE : 1.f;

		for (uint32_t row = height; row-- > 0;) {
			for (uint32_t col = 0; col < width; col++) {
				for (uint16_t n = 0; n < nsamples; n++) {
					buf16[col * nsamples + n] =
							(fit->type == DATA_USHORT) ?
									gbuf[n][(col + row * width)] * norm :
									float_to_ushort_range(gbuff[n][col + row * width]);
				}
			}
			if (TIFFWriteScanline(tif, buf16, height - 1 - row, 0) < 0) {
				siril_debug_print("Error while writing in TIFF File.\n");
				retval = OPEN_IMAGE_ERROR;
				write_ok = FALSE;
				break;
			}
		}
		_TIFFfree(buf16);
		break;
	case 32:
		siril_debug_print("Saving 32-bit TIFF file.\n");
		float *buf32 = _TIFFmalloc(width * sizeof(float) * nsamples);
		if (!buf32) {
			PRINT_ALLOC_ERR;
			retval = OPEN_IMAGE_ERROR;
			write_ok = FALSE;
			break;
		}

		for (uint32_t row = height; row-- > 0;) {
			for (uint32_t col = 0; col < width; col++) {
				for (uint16_t n = 0; n < nsamples; n++) {
					buf32[col * nsamples + n] =
							(fit->type == DATA_USHORT) ?
									(fit->orig_bitpix == BYTE_IMG ?
											gbuf[n][col + row * width] / UCHAR_MAX_SINGLE :
											gbuf[n][col + row * width] / USHRT_MAX_SINGLE) : gbuff[n][col + row * width];
				}
			}
			if (TIFFWriteScanline(tif, buf32, height - 1 - row, 0) < 0) {
				siril_debug_print("Error while writing in TIFF File.\n");
				retval = OPEN_IMAGE_ERROR;
				write_ok = FALSE;
				break;
			}
		}
		_TIFFfree(buf32);
		break;
	default:		// Should not happen
		retval = OPEN_IMAGE_ERROR;
		write_ok = FALSE;
	}

	if (TIFFFlush(tif) != 1) {
		write_ok = FALSE;
	}

	TIFFClose(tif);

	if (!write_ok) {
		siril_log_color_message(_("Saving TIFF: Cannot write TIFF file.\n"), "red");
		retval = OPEN_IMAGE_ERROR;
		if (g_remove(filename))
			fprintf(stderr, "Error removing file\n");
	} else {
		if (verbose)
			siril_log_message(_("Saving TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
				bitspersample, filename, nsamples, width, height);
	}
	free(dest);
	free(profile);
	g_free(filename);
	return retval;
}
#endif	// HAVE_LIBTIFF

/********************* XISF IMPORT *********************/

#ifdef HAVE_LIBXISF

/**
 * Reformats a FITS header to ensure it contains all necessary information
 * in the correct order (SIMPLE, BITPIX, NAXIS, etc.)
 *
 * @param fit Fits structure containing image information
 * @param original_header The original header that might be malformed
 * @return A newly formatted header (to be freed by the caller) or NULL if error
 */
char* format_fits_header_for_xisf(fits *fit, const char *original_header) {
	if (!fit || !original_header) {
		return NULL;
	}

	// Check if the header starts with SIMPLE
	if (strncmp(original_header, "SIMPLE", 6) != 0) {
		// If the header doesn't start with SIMPLE, create a minimal header with essential information
		char *basic_header = malloc(1024); // Arbitrary size for a minimalist header
		if (!basic_header) {
			return NULL;
		}

		// Create minimal header with essential information
		int offset = 0;
		offset += snprintf(basic_header + offset, 1024 - offset,
				"SIMPLE  =                    T / File conforms to FITS standard\n");

		// Add the rest of the original header
		strcat(basic_header, original_header);
		original_header = basic_header;
	} else {
		// If the header starts with SIMPLE, we can duplicate it for modification
		original_header = strdup(original_header);
		if (!original_header) {
			return NULL;
		}
	}

	// Allocate buffer to build the new header (generous estimation)
	size_t header_size = strlen(original_header) + 1024; // Extra for new lines
	char *formatted_header = malloc(header_size);
	if (!formatted_header) {
		free((char*)original_header);
		return NULL;
	}

	formatted_header[0] = '\0'; // Initialize string as empty

	// Extract the first line (SIMPLE)
	char simple_line[81] = {0};
	const char *first_line_end = strchr(original_header, '\n');
	if (first_line_end) {
		size_t first_line_len = first_line_end - original_header;
		if (first_line_len > 80) first_line_len = 80;
		strncpy(simple_line, original_header, first_line_len);
		simple_line[first_line_len] = '\0';
		// Add newline if needed
		strcat(simple_line, "\n");
	} else {
		// If there's no newline, use the entire string
		strncpy(simple_line, original_header, 80);
		simple_line[80] = '\0';
		strcat(simple_line, "\n");
	}

	// Start the formatted header with SIMPLE
	strcat(formatted_header, simple_line);

	// Add BITPIX (position 2)
	char bitpix_line[81] = {0};
	switch (fit->bitpix) {
		case BYTE_IMG:
			snprintf(bitpix_line, 80, "BITPIX  =                    8 / Bits per data value");
			break;
		case USHORT_IMG:
			snprintf(bitpix_line, 80, "BITPIX  =                   16 / Bits per data value");
			break;
		case LONG_IMG:
			snprintf(bitpix_line, 80, "BITPIX  =                   32 / Bits per data value");
			break;
		case FLOAT_IMG:
			snprintf(bitpix_line, 80, "BITPIX  =                  -32 / Bits per data value");
			break;
		case DOUBLE_IMG:
			snprintf(bitpix_line, 80, "BITPIX  =                  -64 / Bits per data value");
			break;
		default:
			snprintf(bitpix_line, 80, "BITPIX  =                   16 / Bits per data value");
	}
	strcat(bitpix_line, "\n");
	strcat(formatted_header, bitpix_line);

	// Add NAXIS (position 3)
	char naxis_line[81] = {0};
	snprintf(naxis_line, 80, "NAXIS   =                    %d / Number of axes", fit->naxis);
	strcat(naxis_line, "\n");
	strcat(formatted_header, naxis_line);

	// Add NAXIS1 (position 4)
	char naxis1_line[81] = {0};
	snprintf(naxis1_line, 80, "NAXIS1  =                 %5ld / Size of the first axis", fit->naxes[0]);
	strcat(naxis1_line, "\n");
	strcat(formatted_header, naxis1_line);

	// Add NAXIS2 (position 5)
	char naxis2_line[81] = {0};
	snprintf(naxis2_line, 80, "NAXIS2  =                 %5ld / Size of the second axis", fit->naxes[1]);
	strcat(naxis2_line, "\n");
	strcat(formatted_header, naxis2_line);

	// Add NAXIS3 if needed (position 6)
	if (fit->naxis > 2) {
		char naxis3_line[81] = {0};
		snprintf(naxis3_line, 80, "NAXIS3  =                 %5ld / Size of the third axis", fit->naxes[2]);
		strcat(naxis3_line, "\n");
		strcat(formatted_header, naxis3_line);
	}

	// Add the other lines from the original header ignoring those already added
	const char *line_start = original_header;
	const char *line_end;

	while ((line_end = strchr(line_start, '\n')) != NULL) {
		char line[81] = {0};
		size_t line_len = line_end - line_start;
		if (line_len > 80) line_len = 80;

		strncpy(line, line_start, line_len);
		line[line_len] = '\0';

		// Skip lines that we've already added
		if (strncmp(line, "SIMPLE", 6) != 0 &&
			strncmp(line, "BITPIX", 6) != 0 &&
			strncmp(line, "NAXIS ", 6) != 0 &&
			strncmp(line, "NAXIS1", 6) != 0 &&
			strncmp(line, "NAXIS2", 6) != 0 &&
			strncmp(line, "NAXIS3", 6) != 0) {

			strcat(formatted_header, line);
			strcat(formatted_header, "\n");
		}

		line_start = line_end + 1;
	}

	// Check if the header ends with END
	if (strstr(formatted_header, "END     ") == NULL &&
		strstr(formatted_header, "END\n") == NULL) {
		strcat(formatted_header, "END\n");
	}

	// Free memory allocated for the duplicated original header
	free((char*)original_header);

	return formatted_header;
}

int readxisf(const char* name, fits *fit, gboolean force_float) {
	struct xisf_data *xdata = (struct xisf_data *) calloc(1, sizeof(struct xisf_data));

	siril_get_xisf_buffer(name, xdata);
	size_t npixels = xdata->width * xdata->height;

	clearfits(fit);

	set_all_keywords_default(fit);

	if (xdata->channelCount == 1)
		fit->naxis = 2;
	else
		fit->naxis = 3;
	fit->rx = xdata->width;
	fit->ry = xdata->height;
	fit->naxes[0] = xdata->width;
	fit->naxes[1] = xdata->height;
	fit->naxes[2] = xdata->channelCount;

	uint32_t *buffer32;
	double *buffer64;
	unsigned char *buffer8;

	switch (xdata->sampleFormat) {
	case BYTE_IMG:
			buffer8 = (unsigned char *)xdata->data;
			fit->data = (WORD *)malloc(npixels * fit->naxes[2] * sizeof(WORD));
			if (!fit->data) {
				siril_log_message(_("Memory allocation error for image data.\n"));
				free(xdata->fitsHeader);
				free(xdata->icc_buffer);
				free(xdata->data);
				free(xdata);
				return -1;
			}

			for (size_t i = 0; i < npixels * fit->naxes[2]; i++) {
				fit->data[i] = (WORD)buffer8[i];
			}

			free(xdata->data);
			xdata->data = NULL;

			fit->pdata[RLAYER] = fit->data;
			fit->pdata[GLAYER] = fit->naxes[2] == 3 ? fit->data + npixels : fit->data;
			fit->pdata[BLAYER] = fit->naxes[2] == 3 ? fit->data + npixels * 2 : fit->data;
			fit->bitpix = fit->orig_bitpix = BYTE_IMG;
			fit->type = DATA_USHORT;

			if (force_float) {
				fit_replace_buffer(fit, ushort8_buffer_to_float(fit->data, npixels * fit->naxes[2]), DATA_FLOAT);
			}
			break;
	case USHORT_IMG:
		fit->data = (WORD *)xdata->data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->naxes[2] == 3 ? fit->data + npixels : fit->data;
		fit->pdata[BLAYER] = fit->naxes[2] == 3 ? fit->data + npixels * 2 : fit->data;
		fit->bitpix = fit->orig_bitpix = USHORT_IMG;
		fit->type = DATA_USHORT;
		if (force_float) {
			fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, npixels * fit->naxes[2]), DATA_FLOAT);
		}
		break;
	case LONG_IMG:
		buffer32 = (uint32_t *)xdata->data;
		fit->fdata = (float *)xdata->data;
		for (int i = 0; i < npixels * fit->naxes[2]; i++)
			fit->fdata[i] = (float)buffer32[i] / 4294967295.f;

		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = fit->naxes[2] == 3 ? fit->fdata + npixels : fit->fdata;
		fit->fpdata[BLAYER] = fit->naxes[2] == 3 ? fit->fdata + npixels * 2 : fit->fdata;
		fit->bitpix = fit->orig_bitpix = FLOAT_IMG;
		fit->type = DATA_FLOAT;
		break;
	case FLOAT_IMG:
		fit->fdata = (float *)xdata->data;
		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = fit->naxes[2] == 3 ? fit->fdata + npixels : fit->fdata;
		fit->fpdata[BLAYER] = fit->naxes[2] == 3 ? fit->fdata + npixels * 2 : fit->fdata;
		fit->bitpix = fit->orig_bitpix = FLOAT_IMG;
		fit->type = DATA_FLOAT;
		break;
	case DOUBLE_IMG:
		buffer64 = (double *)xdata->data;
		fit->fdata = (float *)xdata->data;
		for (int i = 0; i < npixels * fit->naxes[2]; i++)
			fit->fdata[i] = (float)buffer64[i];

		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = fit->naxes[2] == 3 ? fit->fdata + npixels : fit->fdata;
		fit->fpdata[BLAYER] = fit->naxes[2] == 3 ? fit->fdata + npixels * 2 : fit->fdata;
		fit->bitpix = fit->orig_bitpix = FLOAT_IMG;
		fit->type = DATA_FLOAT;
		break;
	default:
		siril_log_message(_("This image type is not handled.\n"));
		free(xdata->fitsHeader);
		free(xdata->icc_buffer);
		free(xdata->data);
		free(xdata);
		return -1;
	}

	/* Assign the ICC profile, if there is one */
	if (xdata->icc_buffer && xdata->icc_length > 0) {
		fit->icc_profile = cmsOpenProfileFromMem(xdata->icc_buffer, xdata->icc_length);
		color_manage(fit, TRUE);
	} else {
		color_manage(fit, FALSE);
	}
	free(xdata->icc_buffer);

	/* let's do it before header parsing. */
	g_snprintf(fit->keywords.row_order, FLEN_VALUE, "%s", "TOP-DOWN");

	// Format the header to ensure it's properly structured
	char *formatted_header = format_fits_header_for_xisf(fit, xdata->fitsHeader);

	if (formatted_header) {
		fit->header = formatted_header;
		int ret = fits_parse_header_str(fit, formatted_header);
		if (ret) {
			siril_debug_print("XISF Header cannot be read despite formatting.\n");
		}
	} else {
		// If formatting fails, use the original header
		fit->header = strdup(xdata->fitsHeader);
		siril_debug_print("Failed to format XISF header, using original.\n");

		int ret = fits_parse_header_str(fit, fit->header);
		if (ret) {
			siril_debug_print("XISF Header cannot be read.\n");
		}
	}

	fits_flip_top_to_bottom(fit);
	siril_log_message(_("Reading XISF: file %s, %ld layer(s), %ux%u pixels\n"),
			name, fit->naxes[2], fit->rx, fit->ry);

	/* free data */
	free(xdata->fitsHeader);
	free(xdata);

	return 0;
}

#endif

/********************* JPEG IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBJPEG
int readjpg(const char* name, fits *fit){
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;

	FILE *f = g_fopen(name, "rb");
	if (f == NULL) {
		siril_log_color_message(_("Sorry but Siril cannot open the file: %s.\n"), "red", name);
		return OPEN_IMAGE_ERROR;
	}
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, f);

#if LIBJPEG_TURBO_VERSION_NUMBER >= 2000000
	/* tell the lib to save APP2 data (ICC profiles) */
	jpeg_save_markers (&cinfo, JPEG_APP0 + 2, 0xffff);
#endif

	(void) jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	size_t npixels = cinfo.output_width * cinfo.output_height;
	WORD *data = malloc(npixels * sizeof(WORD) * 3);
	if (!data) {
		PRINT_ALLOC_ERR;
		fclose(f);
		return OPEN_IMAGE_ERROR;
	}
	WORD *buf[3] = { data, data + npixels, data + npixels * 2 };
	int row_stride = cinfo.output_width * cinfo.output_components;
	JSAMPARRAY pJpegBuffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE,	row_stride, 1);

	while (cinfo.output_scanline < cinfo.output_height) {
		jpeg_read_scanlines(&cinfo, pJpegBuffer, 1);
		for (int i = 0; i < cinfo.output_width; i++) {
			*buf[RLAYER]++ = pJpegBuffer[0][cinfo.output_components * i + 0];
			*buf[GLAYER]++ = pJpegBuffer[0][cinfo.output_components * i + 1];
			*buf[BLAYER]++ = pJpegBuffer[0][cinfo.output_components * i + 2];
		}
	}
	// TODO: this doesn't work despite being the same as the reference djpeg.c
	// Check for an ICC profile
	JOCTET *EmbedBuffer = NULL;
	unsigned int EmbedLen = 0;
#if LIBJPEG_TURBO_VERSION_NUMBER >= 2000000
	if (jpeg_read_icc_profile(&cinfo, &EmbedBuffer, &EmbedLen)) {
		siril_log_message(_("Read ICC profile from JPEG.\n"));
	}
	else if (cinfo.err->msg_code != JWRN_BOGUS_ICC) {
		siril_log_message(_("Cannot read an ICC profile from this JPEG, assuming sRGB.\n"));
	}
#else
	siril_log_message(_("JPEG ICC profile support unavailable, assuming sRGB.\n"));
#endif

	fclose(f);
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);

	clearfits(fit);

	set_all_keywords_default(fit);

	fit->bitpix = fit->orig_bitpix = BYTE_IMG;
	if (cinfo.output_components == 1)
		fit->naxis = 2;
	else
		fit->naxis = 3;
	fit->rx = cinfo.output_width;
	fit->ry = cinfo.output_height;
	fit->naxes[0] = cinfo.output_width;
	fit->naxes[1] = cinfo.output_height;
	fit->naxes[2] = cinfo.output_components;
	fit->data = data;
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data + npixels;
	fit->pdata[BLAYER] = fit->data + npixels * 2;
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	fit->type = DATA_USHORT;
	mirrorx(fit, FALSE);
	fill_date_obs_if_any(fit, name);

	// Initialize ICC profile and display transform
	fits_initialize_icc(fit, (cmsUInt8Number*) EmbedBuffer,
							 (cmsUInt32Number) EmbedLen);
	free(EmbedBuffer);

	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading JPG: file %s, %ld layer(s), %ux%u pixels\n"),
			basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return cinfo.output_components;
}

int savejpg(const char *name, fits *fit, int quality, gboolean verbose) {
	struct jpeg_compress_struct cinfo;    // Basic info for JPEG properties.
	struct jpeg_error_mgr jerr;           // In case of error.

	//## ALLOCATE AND INITIALIZE JPEG COMPRESSION OBJECT
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	char *filename = strdup(name);
	if (!g_str_has_suffix(filename, ".jpg") && (!g_str_has_suffix(filename, ".jpeg"))) {
		filename = str_append(&filename, ".jpg");
	}

	//## OPEN FILE FOR DATA DESTINATION:
	FILE *f = g_fopen(filename, "wb");
	if (f == NULL) {
		siril_log_color_message(_("Siril cannot create JPG file.\n"), "red");
		free(filename);
		return 1;
	}
	jpeg_stdio_dest(&cinfo, f);

	//## SET PARAMETERS FOR COMPRESSION:
	cinfo.image_width  = fit->rx;   // |-- Image width and height in pixels.
	cinfo.image_height = fit->ry;   // |
	cinfo.input_components = fit->naxes[2];     // Number of color components per pixel.
	cinfo.in_color_space = (fit->naxes[2] == 3) ? JCS_RGB : JCS_GRAYSCALE; // Colorspace of input image as RGB.

	WORD *gbuf[3] =	{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	float *gbuff[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	// Apply the colorspace transform
	void *dest = NULL;
	JOCTET *profile = NULL;
	unsigned int profile_len = 0;
	cmsHTRANSFORM save_transform = NULL;
	if (fit->color_managed && fit->icc_profile) {
		void *buf = NULL;
		gboolean threaded = !get_thread_run();
		gboolean src_is_float = (fit->type == DATA_FLOAT);
		buf = src_is_float ? (void *) fit->fdata : (void *) fit->data;
		size_t npixels = fit->rx * fit->ry;
		size_t nchans = fit->naxes[2];
		cmsColorSpaceSignature sig = cmsGetColorSpace(fit->icc_profile);
		cmsUInt32Number trans_type = get_planar_formatter_type(sig, fit->type, FALSE);
		if (src_is_float) {
			dest = (float*) malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(float));
		} else {
			dest = (WORD*) malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(WORD));
		}

#if LIBJPEG_TURBO_VERSION_NUMBER >= 2000000
		if (nchans == 1) { // mono
			cmsHPROFILE srgb_mono_out = NULL;
			switch (com.pref.icc.export_8bit_method) {
				case EXPORT_SRGB:
					srgb_mono_out = gray_srgbtrc();
					save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
					profile = get_icc_profile_data(srgb_mono_out, &profile_len);
					cmsCloseProfile(srgb_mono_out);
					break;
				case EXPORT_WORKING:
					save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_out, trans_type, com.pref.icc.export_intent, 0);
					profile = get_icc_profile_data(com.icc.mono_out, &profile_len);
					break;
				case EXPORT_IMAGE_ICC:
					profile = get_icc_profile_data(fit->icc_profile, &profile_len);
					break;
				default:
					free(dest);
					free(filename);
					return 1;
			}
		} else { // rgb
			switch (com.pref.icc.export_8bit_method) {
				case EXPORT_SRGB:
					save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_out, trans_type, com.pref.icc.export_intent, 0);
					profile = get_icc_profile_data(com.icc.srgb_out, &profile_len);
					break;
				case EXPORT_WORKING:
					save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.working_out, trans_type, com.pref.icc.export_intent, 0);
					profile = get_icc_profile_data(com.icc.working_out, &profile_len);
					break;
				case EXPORT_IMAGE_ICC:
					profile = get_icc_profile_data(fit->icc_profile, &profile_len);
					break;
				default:
					free(dest);
					free(filename);
					return 1;
			}
		}
#else
		if (nchans == 1) {
			cmsHPROFILE srgb_mono_out = gray_srgbtrc();
			save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
			cmsCloseProfile(srgb_mono_out);
		} else {
			save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_out, trans_type, com.pref.icc.export_intent, 0);
		}
#endif
		cmsUInt32Number datasize = fit->type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
		cmsUInt32Number bytesperline = fit->rx * datasize;
		cmsUInt32Number bytesperplane = npixels * datasize;
		if (save_transform) { // save_transform will be NULL if saving in the current image colorspace
			cmsDoTransformLineStride(save_transform, buf, dest, fit->rx, fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
			cmsDeleteTransform(save_transform);
		}
		gbuf[0] = (WORD*) dest;
		gbuf[1] = (WORD*) dest + npixels;
		gbuf[2] = (WORD*) dest + 2 * npixels;
		gbuff[0] = (float*) dest;
		gbuff[1] = (float*) dest + npixels;
		gbuff[2] = (float*) dest + 2 * npixels;
	} else if (com.pref.icc.default_to_srgb) { // Default for non-color managed files is to assume sRGB but there is a preference for this
		profile = get_icc_profile_data(fit->naxes[2] == 1 ? com.icc.mono_out : com.icc.srgb_out, &profile_len);
	}

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	//## CREATE IMAGE BUFFER TO WRITE FROM AND MODIFY THE IMAGE TO LOOK LIKE CHECKERBOARD:
	unsigned char *image_buffer = (unsigned char*) malloc(
			cinfo.image_width * cinfo.image_height * cinfo.num_components);
	if (!image_buffer) {
		PRINT_ALLOC_ERR;
		free(filename);
		fclose(f);
		return 1;
	}

	float norm = (fit->orig_bitpix != BYTE_IMG ?
			UCHAR_MAX_SINGLE / USHRT_MAX_SINGLE : 1.f);

	for (int i = (cinfo.image_height - 1); i >= 0; i--) {
		for (int j = 0; j < cinfo.image_width; j++) {
			int pixelIdx = ((i * cinfo.image_width) + j) * cinfo.input_components;
			if (fit->type == DATA_USHORT) {
				WORD red = *gbuf[RLAYER]++;
				image_buffer[pixelIdx + 0] = round_to_BYTE(red * norm); // r |-- Set r,g,b components to
				if (cinfo.input_components == 3) {
					WORD green = *gbuf[GLAYER]++;
					WORD blue = *gbuf[BLAYER]++;
					image_buffer[pixelIdx + 1] = round_to_BYTE(green * norm); // g |   make this pixel
					image_buffer[pixelIdx + 2] = round_to_BYTE(blue * norm); // b |
				}
			} else {
				float red = *gbuff[RLAYER]++;
				image_buffer[pixelIdx + 0] = float_to_uchar_range(red); // r |-- Set r,g,b components to
				if (cinfo.input_components == 3) {
					float green = *gbuff[GLAYER]++;
					float blue = *gbuff[BLAYER]++;
					image_buffer[pixelIdx + 1] = float_to_uchar_range(green); // g |   make this pixel
					image_buffer[pixelIdx + 2] = float_to_uchar_range(blue); // b |
				}
			}
		}
	}
	//## START COMPRESSION:
	jpeg_start_compress(&cinfo, TRUE);

	// Write the ICC profile, if the image is color managed
#if LIBJPEG_TURBO_VERSION_NUMBER >= 2000000
	if (profile)
		jpeg_write_icc_profile(&cinfo, (const JOCTET*) profile, profile_len);
	else
		siril_log_color_message(_("Error: failed to write ICC profile to JPG\n"), "red");
#endif

	int row_stride = cinfo.image_width * cinfo.input_components;        // JSAMPLEs per row in image_buffer

	JSAMPROW row_pointer[1];
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &image_buffer[cinfo.next_scanline * row_stride];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	// NOTE: jpeg_write_scanlines expects an array of pointers to scanlines.
	//       Here the array is only one element long, but you could pass
	//       more than one scanline at a time if that's more convenient.

	//## FINISH COMPRESSION AND CLOSE FILE:
	jpeg_finish_compress(&cinfo);

	fclose(f);
	jpeg_destroy_compress(&cinfo);
	free(image_buffer);
	if (verbose)
		siril_log_message(_("Saving JPG: file %s, quality=%d%%, %ld layer(s), %ux%u pixels\n"),
						filename, quality, fit->naxes[2], fit->rx, fit->ry);
	free(filename);
	free(profile);
	free(dest);
	return OPEN_IMAGE_OK;
}

#endif	// HAVE_LIBJPEG

/********************* PNG IMPORT *********************/

#ifdef HAVE_LIBPNG
/* reads a PNG file and stores it in the fits argument.
 */
int readpng(const char *name, fits* fit) {
	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL,	NULL);
	if (!png) {
		siril_log_color_message(_("Sorry but Siril cannot open the file: %s.\n"), "red", name);
		return OPEN_IMAGE_ERROR;
	}

	png_infop info = png_create_info_struct(png);
	png_infop end_info = png_create_info_struct(png);
	if (!info || !end_info)
		return OPEN_IMAGE_ERROR;

	if (setjmp(png_jmpbuf(png)))
		return OPEN_IMAGE_ERROR;

	FILE *f = g_fopen(name, "rb");
	if (!f) {
		siril_log_color_message(_("Error opening the file %s\n"), "red", name);
		return OPEN_IMAGE_ERROR;
	}
	png_init_io(png, f);

	png_read_info(png, info);

	const int width = png_get_image_width(png, info);
	const int height = png_get_image_height(png, info);
	size_t npixels = width * height;
	png_byte color_type = png_get_color_type(png, info);
	png_byte bit_depth = png_get_bit_depth(png, info);

	WORD *data = malloc(npixels * sizeof(WORD) * 3);
	if (!data) {
		PRINT_ALLOC_ERR;
		fclose(f);
		png_destroy_read_struct(&png, &info, &end_info);
		return OPEN_IMAGE_ERROR;
	}
	WORD *buf[3] = { data, data + npixels, data + npixels * 2 };

	if (color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_palette_to_rgb(png);

	// PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
	if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
		png_set_expand_gray_1_2_4_to_8(png);

	if (png_get_valid(png, info, PNG_INFO_tRNS))
		png_set_tRNS_to_alpha(png);

	// These color_type don't have an alpha channel then fill it with 0xff.
	if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY
			|| color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

	if (color_type == PNG_COLOR_TYPE_GRAY
			|| color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png);

	png_read_update_info(png, info);

	png_bytep *row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
	for (int y = 0; y < height; y++) {
		row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png, info));
	}

	cmsUInt8Number *embed = NULL;
	cmsUInt32Number len = 0;
	{
		png_charp name_str;
		int comp_type;
#if ((PNG_LIBPNG_VER_MAJOR << 8) | PNG_LIBPNG_VER_MINOR << 0) < \
    ((1 << 8) | (5 << 0))
		png_charp profile;
#else  // >= libpng 1.5.0
		png_bytep profile;
#endif
		if (png_get_iCCP(png, info,
						&name_str, &comp_type, &profile, &len) ==
						PNG_INFO_iCCP) {
			embed = malloc(len * sizeof(cmsUInt8Number));
			memcpy(embed, profile, len * sizeof(cmsUInt8Number));
		}
    }

	png_read_image(png, row_pointers);

	fclose(f);
	png_destroy_read_struct(&png, &info, &end_info);

	if (bit_depth == 16) {			//in 16-bit: it is stored as RRGGBB
		for (int y = height - 1; y > -1; y--) {
			png_byte* row = row_pointers[y];
			for (int x = 0; x < width; x++) {
				const png_byte* ptr = &(row[x * 8]);
				*buf[RLAYER]++ = (ptr[0] << 8) + ptr[1];
				*buf[GLAYER]++ = (ptr[2] << 8) + ptr[3];
				*buf[BLAYER]++ = (ptr[4] << 8) + ptr[5];
			}
		}
	} else {
		for (int y = height - 1; y > -1; y--) {
			png_byte* row = row_pointers[y];
			for (int x = 0; x < width; x++) {
				const png_byte* ptr = &(row[x * 4]);
				*buf[RLAYER]++ = ptr[0];
				*buf[GLAYER]++ = ptr[1];
				*buf[BLAYER]++ = ptr[2];
			}

		}
	}
	// We define the number of channel we have
	int nbplanes;
	switch (color_type) {
	case PNG_COLOR_TYPE_RGB_ALPHA:
	case PNG_COLOR_TYPE_RGB:
		nbplanes = 3;
		break;
	case PNG_COLOR_TYPE_GRAY_ALPHA:
	case PNG_COLOR_TYPE_GRAY:
		nbplanes = 1;
		break;
	default:
		nbplanes = 0;
	}

	// free allocated memory
	for (int y = 0; y < height; y++)
		free(row_pointers[y]);
	free(row_pointers);

	if (data != NULL) {
		clearfits(fit);

		set_all_keywords_default(fit);

		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->naxes[2] = nbplanes;
		if (nbplanes == 1)
			fit->naxis = 2;
		else
			fit->naxis = 3;
		fit->bitpix = (bit_depth == 8 ? BYTE_IMG : USHORT_IMG);
		fit->type = DATA_USHORT;
		fit->orig_bitpix = fit->bitpix;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->naxes[2] == 3 ? fit->data + npixels : fit->data;
		fit->pdata[BLAYER] = fit->naxes[2] == 3 ? fit->data + npixels * 2 : fit->data;

		fit->keywords.binning_x = fit->keywords.binning_y = 1;
		g_snprintf(fit->keywords.row_order, FLEN_VALUE, "%s", "TOP-DOWN");
		fill_date_obs_if_any(fit, name);
		// Initialize ICC profile and display transform
		fits_initialize_icc(fit, embed, len);
		free(embed);
	}

	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading PNG: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
			bit_depth, basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return nbplanes;
}

static WORD *convert_data(fits *image) {
	size_t ndata = image->rx * image->ry;
	int ch = image->naxes[2];

	WORD *buffer = malloc(ndata * ch * sizeof(WORD));
	if (!buffer) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	for (size_t i = 0, j = 0; i < ndata * ch; i += ch, j++) {
		if (image->type == DATA_USHORT) {
			buffer[i + 0] = image->pdata[RLAYER][j];
			if (ch > 1) {
				buffer[i + 1] = image->pdata[GLAYER][j];
				buffer[i + 2] = image->pdata[BLAYER][j];
			}
		} else if (image->type == DATA_FLOAT) {
			buffer[i + 0] = float_to_ushort_range(image->fpdata[RLAYER][j]);
			if (ch > 1) {
				buffer[i + 1] = float_to_ushort_range(image->fpdata[GLAYER][j]);
				buffer[i + 2] = float_to_ushort_range(image->fpdata[BLAYER][j]);
			}
		}
		else {
			free(buffer);
			return NULL;
		}
	}
	return buffer;
}

static uint8_t *convert_data8(fits *image) {
	size_t ndata = image->rx * image->ry;
	const long ch = image->naxes[2];

	uint8_t *buffer = malloc(ndata * ch * sizeof(uint8_t));
	if (!buffer) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	for (size_t i = 0, j = 0; i < ndata * ch; i += ch, j++) {
		if (image->type == DATA_USHORT) {
			buffer[i + 0] = (uint8_t) image->pdata[RLAYER][j];
			if (ch > 1) {
				buffer[i + 1] = (uint8_t) image->pdata[GLAYER][j];
				buffer[i + 2] = (uint8_t) image->pdata[BLAYER][j];
			}
		} else if (image->type == DATA_FLOAT) {
			buffer[i + 0] = float_to_uchar_range(image->fpdata[RLAYER][j]);
			if (ch > 1) {
				buffer[i + 1] = float_to_uchar_range(image->fpdata[GLAYER][j]);
				buffer[i + 2] = float_to_uchar_range(image->fpdata[BLAYER][j]);
			}
		}
		else {
			free(buffer);
			return NULL;
		}
	}
	return buffer;
}

int savepng(const char *name, fits *fit, uint32_t bytes_per_sample,
		gboolean is_colour) {
	int32_t ret = -1;
	png_structp png_ptr;
	png_infop info_ptr;
	const uint32_t width = fit->rx;
	const uint32_t height = fit->ry;

	char *filename = strdup(name);
	if (!g_str_has_suffix(filename, ".png")) {
		filename = str_append(&filename, ".png");
	}

	FILE *p_png_file = g_fopen(filename, "wb");
	if (p_png_file == NULL) {
		free(filename);
		return ret;
	}

	/* Create and initialize the png_struct with the desired error handler
	 * functions.  If you want to use the default stderr and longjump method,
	 * you can supply NULL for the last three parameters.  We also check that
	 * the library version is compatible with the one used at compile time,
	 * in case we are using dynamically linked libraries.  REQUIRED.
	 */
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fclose(p_png_file);
		free(filename);
		return ret;
	}

	/* Allocate/initialize the image information data.  REQUIRED */
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fclose(p_png_file);
		png_destroy_write_struct(&png_ptr, NULL);
		free(filename);
		return ret;
	}

	/* Set error handling.  REQUIRED if you aren't supplying your own
	 * error handling functions in the png_create_write_struct() call.
	 */
	if (setjmp(png_jmpbuf(png_ptr))) {
		/* If we get here, we had a problem writing the file */
		fclose(p_png_file);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		free(filename);
		return ret;
	}

	/* Set up the output control if you are using standard C streams */
	png_init_io(png_ptr, p_png_file);

	/* Set the image information here.  Width and height are up to 2^31,
	 * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
	 * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
	 * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
	 * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
	 * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
	 * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
	 */
	uint32_t profile_len = 0;
	unsigned char *profile = NULL;

	if (is_colour) {
		png_set_IHDR(png_ptr, info_ptr, width, height, bytes_per_sample * 8,
				PNG_COLOR_TYPE_RGB,
				PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
				PNG_FILTER_TYPE_DEFAULT);
	} else {
		png_set_IHDR(png_ptr, info_ptr, width, height, bytes_per_sample * 8,
				PNG_COLOR_TYPE_GRAY,
				PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
				PNG_FILTER_TYPE_DEFAULT);
	}

	int samples_per_pixel;
	if (is_colour) {
		samples_per_pixel = 3;
	} else {
		samples_per_pixel = 1;
	}

	gboolean threaded = get_thread_run();
	cmsHTRANSFORM save_transform = NULL;
	cmsUInt32Number trans_type;
	if (fit->color_managed && fit->icc_profile) {
		// Check what is the appropriate color space to save in
		// Covers 8- and high-bitdepth files
		if (bytes_per_sample == 1) { // 8-bit
			trans_type = (fit->naxes[2] == 1 ? TYPE_GRAY_8 : TYPE_RGB_8);
			if (samples_per_pixel == 1) { // mono
				cmsHPROFILE srgb_mono_out = NULL;
				switch (com.pref.icc.export_8bit_method) {
					case EXPORT_SRGB:
						srgb_mono_out = gray_srgbtrc();
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(srgb_mono_out, &profile_len);
						cmsCloseProfile(srgb_mono_out);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.mono_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			} else { // rgb
				switch (com.pref.icc.export_8bit_method) {
					case EXPORT_SRGB:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.srgb_out, &profile_len);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.working_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.working_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			}
		} else { // high bitdepth
			trans_type = (fit->naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16);
			if (samples_per_pixel == 1) { // mono
				cmsHPROFILE srgb_mono_out = NULL;
				switch (com.pref.icc.export_16bit_method) {
					case EXPORT_SRGB:
						srgb_mono_out = gray_srgbtrc();
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(srgb_mono_out, &profile_len);
						cmsCloseProfile(srgb_mono_out);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.mono_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			} else { // rgb
				switch (com.pref.icc.export_16bit_method) {
					case EXPORT_SRGB:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.srgb_out, &profile_len);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.working_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.working_out, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			}
		}
	} else if(com.pref.icc.default_to_srgb) {
		profile = get_icc_profile_data((samples_per_pixel == 1 ? com.icc.mono_out : com.icc.srgb_out), &profile_len);
	}

	if (profile && profile_len > 0) {
			png_set_iCCP(png_ptr, info_ptr, "icc", 0, (png_const_bytep) profile, profile_len);
	}

	/* Write the file header information.  REQUIRED */
	png_write_info(png_ptr, info_ptr);

	png_bytep *row_pointers = malloc((size_t) height * sizeof(png_bytep));

	WORD *data = NULL;
	uint8_t *data8 = NULL;
	if (bytes_per_sample == 2) {
		/* swap bytes of 16 bit files to most significant bit first */
		png_set_swap(png_ptr);
		data = convert_data(fit);
		// Apply ICC transform (only for color managed images)
		if (fit->color_managed && fit->icc_profile && save_transform) {
			cmsUInt32Number datasize = sizeof(WORD);
			cmsUInt32Number bytesperline = fit->rx * datasize * fit->naxes[2];
			cmsUInt32Number bytesperplane = fit->rx * fit->ry * datasize;
			cmsDoTransformLineStride(save_transform, data, data, fit->rx, fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
			cmsDeleteTransform(save_transform);
		}
		for (unsigned i = 0, j = height - 1; i < height; i++)
			row_pointers[j--] = (png_bytep) ((uint16_t*) data + (size_t) samples_per_pixel * i * width);
	} else {
		data8 = convert_data8(fit);
		// Apply ICC transform
		if (fit->color_managed && fit->icc_profile && save_transform) {
			cmsUInt32Number datasize = sizeof(BYTE);
			cmsUInt32Number bytesperline = fit->rx * datasize;
			cmsUInt32Number bytesperplane = fit->rx * fit->ry * datasize;
			cmsDoTransformLineStride(save_transform, data8, data8, fit->rx, fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
			cmsDeleteTransform(save_transform);
		}
		for (unsigned i = 0, j = height - 1; i < height; i++)
			row_pointers[j--] = (uint8_t*) data8 + (size_t) samples_per_pixel * i * width;
	}

	png_write_image(png_ptr, row_pointers);

	/* Clean up after the write, and free any memory allocated */
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);

	siril_log_message(_("Saving PNG: file %s, %ld layer(s), %ux%u pixels\n"),
				filename, fit->naxes[2], fit->rx, fit->ry);

	/* Close the file */
	fclose(p_png_file);
	if (data) free(data);
	if (data8) free(data8);
	free(row_pointers);
	free(profile);
	free(filename);
	return 0;
}
#endif	// HAVE_LIBPNG

/********************* RAW IMPORT *********************/
#ifdef HAVE_LIBRAW

#if LIBRAW_VERSION < LIBRAW_MAKE_VERSION(0, 18, 0)
#define LIBRAW_FORMAT_1INCH 5
#endif

/* this is an estimation of the pixel size. Indeed, we cannot know
 * the real width resolution with libraw.
 * However, this approximation should be good enough.
 */
static float estimate_pixel_pitch(libraw_data_t *raw) {
	float s_width;

	switch (raw->lens.makernotes.CameraFormat) {
	case LIBRAW_FORMAT_APSC:
		if (!g_ascii_strncasecmp("Canon", raw->idata.make, 5))
			s_width = 22.3f;
		else
			s_width = 23.6f;
		break;
	case LIBRAW_FORMAT_FF:
		if (!g_ascii_strncasecmp("Sony", raw->idata.make, 4))
			s_width = 35.6f;
		else
			s_width = 36.0f;
		break;
	case LIBRAW_FORMAT_FT:
		s_width = 17.3f;
		break;
	case LIBRAW_FORMAT_APSH:
		s_width = 28.7f;
		break;
	case LIBRAW_FORMAT_1INCH:
		s_width = 13.2f;
		break;
	case LIBRAW_FORMAT_MF:
		s_width = 44.0f;
		break;
	default:
		s_width = 0.0f;
		break;
	}
//	printf("s_width=%f\n", s_width);
	float pitch = s_width / (float) raw->sizes.width * 1000.f;
	return roundf(pitch * 100.f) / 100.f;
}

static int siril_libraw_open_file(libraw_data_t* rawdata, const char *name) {
/* libraw_open_wfile is not defined for all windows compilers in previous LibRaw versions */
#if (defined(_WIN32) && !defined(__MINGW32__) && defined(_MSC_VER) && (_MSC_VER > 1310)) || (defined(_WIN32) && LIBRAW_VERSION >= LIBRAW_MAKE_VERSION(0, 21, 0))
	wchar_t *wname;

	wname = g_utf8_to_utf16(name, -1, NULL, NULL, NULL);
	if (wname == NULL) {
		return 1;
	}

	int ret = libraw_open_wfile(rawdata, wname);
	g_free(wname);
	return ret;
#elif defined(_WIN32)
	gchar *localefilename = g_win32_locale_filename_from_utf8(name);
	int ret = libraw_open_file(rawdata, localefilename);
	g_free(localefilename);
	return ret;
#else
	return(libraw_open_file(rawdata, name));
#endif
}

#define FC(filters, row, col) \
	(filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)

static const char filter[16][16] =
{ { 2,1,1,3,2,3,2,0,3,2,3,0,1,2,1,0 },
  { 0,3,0,2,0,1,3,1,0,1,1,2,0,3,3,2 },
  { 2,3,3,2,3,1,1,3,3,1,2,1,2,0,0,3 },
  { 0,1,0,1,0,2,0,2,2,0,3,0,1,3,2,1 },
  { 3,1,1,2,0,1,0,2,1,3,1,3,0,1,3,0 },
  { 2,0,0,3,3,2,3,1,2,0,2,0,3,2,2,1 },
  { 2,3,3,1,2,1,2,1,2,1,1,2,3,0,0,1 },
  { 1,0,0,2,3,0,0,3,0,3,0,3,2,1,2,3 },
  { 2,3,3,1,1,2,1,0,3,2,3,0,2,3,1,3 },
  { 1,0,2,0,3,0,3,2,0,1,1,2,0,1,0,2 },
  { 0,1,1,3,3,2,2,1,1,3,3,0,2,1,3,2 },
  { 2,3,2,0,0,1,3,0,2,0,1,2,3,0,1,0 },
  { 1,3,1,2,3,2,3,2,0,2,0,1,1,0,3,0 },
  { 0,2,0,3,1,0,0,1,1,3,3,2,3,2,2,1 },
  { 2,1,3,2,3,1,2,1,0,3,0,2,0,2,0,2 },
  { 0,3,1,0,0,2,0,3,2,1,3,1,1,3,1,3 } };

static int fcol(libraw_data_t *raw, int row, int col) {
	if (raw->idata.filters == 1)
		return filter[(row + raw->rawdata.sizes.top_margin) & 15][(col
				+ raw->rawdata.sizes.left_margin) & 15];
	if (raw->idata.filters == 9)
		return raw->idata.xtrans[(row + 6) % 6][(col + 6) % 6];
	return FC(raw->idata.filters, row, col);
}

static int readraw_in_cfa(const char *name, fits *fit) {
	libraw_data_t *raw = libraw_init(0);
	char pattern[FLEN_VALUE];

	int ret = siril_libraw_open_file(raw, name);
	if (ret) {
		siril_log_color_message("Error in libraw %s\n", "red", libraw_strerror(ret));
		return OPEN_IMAGE_ERROR;
	}

	ret = libraw_unpack(raw);
	if (ret) {
		siril_log_color_message("Error in libraw %s\n", "red", libraw_strerror(ret));
		return OPEN_IMAGE_ERROR;
	}

	/* This test checks if raw data exist. Sometimes it doesn't. This is
	 * the case for DNG built from lightroom for example */
	if (raw->rawdata.raw_image == NULL
			&& (raw->rawdata.color3_image || raw->rawdata.color4_image)) {
		siril_log_color_message(_("Siril cannot open this file in CFA mode: "
				"no RAW data available.\n"), "red");
		return OPEN_IMAGE_ERROR;
	}

	raw->params.user_flip = 0;				/* no flip                                 */
	raw->params.output_color = 0;			/* output colorspace, 0=raw, 1=sRGB, 2=Adobe, 3=Wide, 4=ProPhoto, 5=XYZ*/

	const ushort raw_width = raw->sizes.raw_width;
	const ushort raw_height = raw->sizes.raw_height;
	const ushort left_margin = raw->rawdata.sizes.left_margin;
	const ushort top_margin = raw->rawdata.sizes.top_margin;

	ushort width, height;

	if (raw->rawdata.ioparams.fuji_width) {
		const ushort right_margin = raw_width - raw->rawdata.ioparams.fuji_width
				- left_margin;
		width = raw_width - right_margin;
		height = raw_height;
	} else {
		width = raw->sizes.iwidth;
		height = raw->sizes.iheight;
	}

	float pitch = estimate_pixel_pitch(raw);
	size_t npixels = (size_t) width * (size_t) height;

	if (raw->other.shutter > 0 && raw->other.shutter < 1)
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=1/%0.1f sec)\n"),
						raw->idata.make, raw->idata.model, raw->other.iso_speed, 1/raw->other.shutter);
	else
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=%0.1f sec)\n"),
						raw->idata.make, raw->idata.model, raw->other.iso_speed, raw->other.shutter);

	unsigned filters = raw->idata.filters;

	if (filters) {
		int fhigh = 2, fwide = 2;
		if ((filters ^ (filters >> 8)) & 0xff)
			fhigh = 4;
		if ((filters ^ (filters >> 16)) & 0xffff)
			fhigh = 8;
		if (filters == 1) /* Leaf Catchlight with 16x16 bayer matrix */
			fhigh = fwide = 16;
		if (filters == 9) /* Fuji X-Trans (6x6 matrix) */
			fhigh = fwide = 6;

		int j = 0;
		for (int i = 0; i < fhigh; i++) {
			for (int c = i /*&& (pattern[j++] = '/')*/ && 0; c < fwide; c++) {
				pattern[j++] = raw->idata.cdesc[fcol(raw, i, c)];
			}
		}
		pattern[j++] = '\0';
		siril_log_message(_("Filter pattern: %s\n"), pattern);
	}

	WORD *data = (WORD*) calloc(1, npixels * sizeof(WORD));
	if (!data) {
		PRINT_ALLOC_ERR;
		libraw_recycle(raw);
		libraw_close(raw);
		return OPEN_IMAGE_ERROR;
	}

	WORD *buf = data;

	int offset = raw_width * top_margin + left_margin;

	if (!raw->rawdata.raw_image) {
		libraw_recycle(raw);
		libraw_close(raw);
		free(buf);
		return OPEN_IMAGE_ERROR;
	}
	int i = 0;
	for (int row = height - 1; row > -1; row--) {
		for (int col = 0; col < width; col++) {
			buf[i++] = raw->rawdata.raw_image[offset + col + (raw_width * row)];
		}
	}

	clearfits(fit);

	set_all_keywords_default(fit);

	fit->bitpix = fit->orig_bitpix = USHORT_IMG;
	fit->type = DATA_USHORT;
	fit->rx = (unsigned int) (width);
	fit->ry = (unsigned int) (height);
	fit->naxes[0] = (long) (width);
	fit->naxes[1] = (long) (height);
	fit->naxes[2] = 1;
	fit->naxis = 2;
	fit->data = data;
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data;
	fit->pdata[BLAYER] = fit->data;
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	// RAW files are always mono, and pretty certainly always linear: they
	// should have the mono_linear ICC profile. However for consistency with
	// other straight-from-the-camera formats we do not set a profile: the user
	// may assign one if they wish.

	color_manage(fit, FALSE);
	if (pitch > 0.f) {
		fit->keywords.pixel_size_x = fit->keywords.pixel_size_y = pitch;
		fit->pixelkey = TRUE;
	}
	if (raw->other.focal_len > 0.f) {
		fit->keywords.focal_length = raw->other.focal_len;
		fit->focalkey = TRUE;
	}
	if (raw->other.iso_speed > 0.f)
		fit->keywords.iso_speed = raw->other.iso_speed;
	if (raw->other.shutter > 0.f)
		fit->keywords.exposure = raw->other.shutter;
	if (raw->other.aperture > 0.f)
		fit->keywords.aperture = raw->other.aperture;
	g_snprintf(fit->keywords.instrume, FLEN_VALUE, "%s %s", raw->idata.make,
			raw->idata.model);
	fit->keywords.date_obs = g_date_time_new_from_unix_utc(raw->other.timestamp);
	if (filters)
		g_snprintf(fit->keywords.bayer_pattern, FLEN_VALUE, "%s", pattern);

	g_snprintf(fit->keywords.row_order, FLEN_VALUE, "%s", "BOTTOM-UP");

	libraw_recycle(raw);
	libraw_close(raw);
	return 1;
}

int open_raw_files(const char *name, fits *fit) {
	int retval = readraw_in_cfa(name, fit);

	if (retval >= 0) {
		gchar *basename = g_path_get_basename(name);
		siril_log_message(_("Reading RAW: file %s, %ld layer(s), %ux%u pixels\n"),
				basename, fit->naxes[2], fit->rx, fit->ry);
		g_free(basename);
	}
	return retval;
}
#endif

#ifdef HAVE_LIBHEIF
#define MAX_THUMBNAIL_SIZE com.pref.gui.thumbnail_size

struct HeifImage {
	uint32_t ID;
	char caption[100]; // image text (filled with resolution description)
	struct heif_image *thumbnail;
	int width, height;
};

static gboolean load_thumbnails(struct heif_context *heif, struct HeifImage *images) {
	int numImages = heif_context_get_number_of_top_level_images(heif);

	// get list of all (top level) image IDs

	uint32_t *IDs = malloc(numImages * sizeof(uint32_t));
	heif_context_get_list_of_top_level_image_IDs(heif, IDs, numImages);

	// --- Load a thumbnail for each image.

	for (int i = 0; i < numImages; i++) {

		images[i].ID = IDs[i];
		images[i].caption[0] = 0;
		images[i].thumbnail = NULL;

		// get image handle

		struct heif_image_handle *handle;
		struct heif_error err = heif_context_get_image_handle(heif, IDs[i],
				&handle);
		if (err.code) {
			g_printf("%s\n", err.message);
			continue;
		}

		// generate image caption

		int width = heif_image_handle_get_width(handle);
		int height = heif_image_handle_get_height(handle);

		if (heif_image_handle_is_primary_image(handle)) {
			sprintf(images[i].caption, "%dx%d (%s)", width, height,
					_("primary"));
		} else {
			sprintf(images[i].caption, "%dx%d", width, height);
		}

		// get handle to thumbnail image
		// if there is no thumbnail image, just the the image itself (will be scaled down later)

		struct heif_image_handle *thumbnail_handle;
		heif_item_id thumbnail_ID;

		int nThumbnails = heif_image_handle_get_list_of_thumbnail_IDs(handle,
				&thumbnail_ID, 1);

		if (nThumbnails > 0) {
			err = heif_image_handle_get_thumbnail(handle, thumbnail_ID,
					&thumbnail_handle);
			if (err.code) {
				g_printf("%s\n", err.message);
				continue;
			}
		} else {
			err = heif_context_get_image_handle(heif, IDs[i],
					&thumbnail_handle);
			if (err.code) {
				g_printf("%s\n", err.message);
				continue;
			}
		}

		// decode the thumbnail image

		struct heif_image *thumbnail_img;
		err = heif_decode_image(thumbnail_handle, &thumbnail_img,
				heif_colorspace_RGB, heif_chroma_interleaved_RGB,
				NULL);
		if (err.code) {
			g_printf("%s\n", err.message);
			continue;
		}

		// if thumbnail image size exceeds the maximum, scale it down

		int thumbnail_width = heif_image_handle_get_width(thumbnail_handle);
		int thumbnail_height = heif_image_handle_get_height(thumbnail_handle);

		if (thumbnail_width > MAX_THUMBNAIL_SIZE
				|| thumbnail_height > MAX_THUMBNAIL_SIZE) {

			// compute scaling factor to fit into a max sized box

			float factor_h = thumbnail_width / (float) MAX_THUMBNAIL_SIZE;
			float factor_v = thumbnail_height / (float) MAX_THUMBNAIL_SIZE;

			int new_width, new_height;

			if (factor_v > factor_h) {
				new_height = MAX_THUMBNAIL_SIZE;
				new_width = thumbnail_width / factor_v;
			} else {
				new_height = thumbnail_height / factor_h;
				new_width = MAX_THUMBNAIL_SIZE;
			}

			// scale the image

			struct heif_image *scaled_img = NULL;

			err = heif_image_scale_image(thumbnail_img,
					&scaled_img, new_width, new_height,
					NULL);
			if (err.code) {
				g_printf("%s\n", err.message);
				continue;
			}

			// release the old image and only keep the scaled down version

			heif_image_release(thumbnail_img);
			thumbnail_img = scaled_img;

			thumbnail_width = new_width;
			thumbnail_height = new_height;
		}

		heif_image_handle_release(thumbnail_handle);
		heif_image_handle_release(handle);

		// remember the HEIF thumbnail image (we need it for the GdkPixbuf)

		images[i].thumbnail = thumbnail_img;

		images[i].width = thumbnail_width;
		images[i].height = thumbnail_height;
	}

	return TRUE;
}

static gboolean heif_dialog(struct heif_context *heif, uint32_t *selected_image) {
	int numImages = heif_context_get_number_of_top_level_images(heif);

	struct HeifImage *heif_images = malloc(numImages * sizeof(struct HeifImage));
	gboolean success = load_thumbnails(heif, heif_images);
	if (!success) {
		free(heif_images);
		return FALSE;
	}

	GtkWidget *dlg = gtk_dialog_new_with_buttons(_("Load HEIF image content"),
			GTK_WINDOW(lookup_widget("control_window")), GTK_DIALOG_MODAL,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_OK"), GTK_RESPONSE_OK, NULL);
	gtk_dialog_set_default_response(GTK_DIALOG(dlg), GTK_RESPONSE_OK);

	GtkContainer *content_area = GTK_CONTAINER(gtk_dialog_get_content_area(GTK_DIALOG(dlg)));
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 12);

	GtkWidget *frame = gtk_frame_new(_("Select image"));
	gtk_container_add(content_area, GTK_WIDGET(frame));
	gtk_widget_show(frame);

// prepare list store with all thumbnails and caption

	GtkListStore *liststore;
	GtkTreeIter iter;

	liststore = gtk_list_store_new(2, G_TYPE_STRING, GDK_TYPE_PIXBUF);

	for (int i = 0; i < numImages; i++) {
		gtk_list_store_append(liststore, &iter);
		gtk_list_store_set(liststore, &iter, 0, heif_images[i].caption, -1);

		int stride;
		const uint8_t *data = heif_image_get_plane_readonly(
				heif_images[i].thumbnail, heif_channel_interleaved, &stride);

		GdkPixbuf *pixbuf = gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB,
				FALSE, 8, heif_images[i].width, heif_images[i].height, stride,
				NULL, NULL);
		gtk_list_store_set(liststore, &iter, 1, pixbuf, -1);
	}

	GtkWidget *iconview = gtk_icon_view_new();
	gtk_icon_view_set_model((GtkIconView*) iconview, (GtkTreeModel*) liststore);
	gtk_icon_view_set_text_column((GtkIconView*) iconview, 0);
	gtk_icon_view_set_pixbuf_column((GtkIconView*) iconview, 1);
	gtk_icon_view_set_item_width((GtkIconView*) iconview, MAX_THUMBNAIL_SIZE);

	GtkWidget *scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_widget_set_size_request(scroll, -1, 400);
	g_object_set(scroll, "expand", TRUE, NULL);

	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),	GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);

	gtk_container_add(GTK_CONTAINER(frame), scroll);
	gtk_container_add(GTK_CONTAINER(scroll), iconview);

	gtk_widget_show(scroll);
	gtk_widget_show(iconview);

// pre-select the primary image

	int selected_idx = -1;
	for (int i = 0; i < numImages; i++) {
		if (heif_images[i].ID == *selected_image) {
			selected_idx = i;
			break;
		}
	}

	if (selected_idx != -1) {
		GtkTreePath *path = gtk_tree_path_new_from_indices(selected_idx, -1);
		gtk_icon_view_select_path((GtkIconView*) iconview, path);
		gtk_tree_path_free(path);
	}

	gtk_widget_show(dlg);

	gboolean run = (gtk_dialog_run(GTK_DIALOG(dlg)) == GTK_RESPONSE_OK);

	if (run) {
		GList *selected_items = gtk_icon_view_get_selected_items(
				(GtkIconView*) iconview);

		if (selected_items) {
			GtkTreePath *path = (GtkTreePath*) (selected_items->data);
			const gint *indices = gtk_tree_path_get_indices(path);

			*selected_image = heif_images[indices[0]].ID;

			g_list_free_full(selected_items,
					(GDestroyNotify) gtk_tree_path_free);
		}
	}

	gtk_widget_destroy(dlg);

// release thumbnail images

	for (int i = 0; i < numImages; i++) {
		heif_image_release(heif_images[i].thumbnail);
	}

	free(heif_images);

	return run;
}

static cmsHPROFILE nclx_to_icc_profile (const struct heif_color_profile_nclx *nclx) {
	const gchar *primaries_name = "";
	const gchar *trc_name = "";
	cmsHPROFILE profile = NULL;
	cmsCIExyY whitepoint;
	cmsCIExyYTRIPLE primaries;
	cmsToneCurve *curve[3];

	cmsFloat64Number srgb_parameters[5] =
	{ 2.4, 1.0 / 1.055,  0.055 / 1.055, 1.0 / 12.92, 0.04045 };

	cmsFloat64Number rec709_parameters[5] =
	{ 2.2, 1.0 / 1.099,  0.099 / 1.099, 1.0 / 4.5, 0.081 };

	if (nclx == NULL) {
		return NULL;
	}

	if (nclx->color_primaries == heif_color_primaries_unspecified) {
		return NULL;
	}

	if (nclx->color_primaries == heif_color_primaries_ITU_R_BT_709_5) {
		if (nclx->transfer_characteristics == heif_transfer_characteristic_IEC_61966_2_1) {
			return srgb_trc();
		}

		if (nclx->transfer_characteristics == heif_transfer_characteristic_linear) {
			return srgb_linear();
		}
	}

	whitepoint.x = nclx->color_primary_white_x;
	whitepoint.y = nclx->color_primary_white_y;
	whitepoint.Y = 1.0f;

	primaries.Red.x = nclx->color_primary_red_x;
	primaries.Red.y = nclx->color_primary_red_y;
	primaries.Red.Y = 1.0f;

	primaries.Green.x = nclx->color_primary_green_x;
	primaries.Green.y = nclx->color_primary_green_y;
	primaries.Green.Y = 1.0f;

	primaries.Blue.x = nclx->color_primary_blue_x;
	primaries.Blue.y = nclx->color_primary_blue_y;
	primaries.Blue.Y = 1.0f;

	switch (nclx->color_primaries) {
		case heif_color_primaries_ITU_R_BT_709_5:
			primaries_name = "BT.709";
			break;
		case   heif_color_primaries_ITU_R_BT_470_6_System_M:
			primaries_name = "BT.470-6 System M";
			break;
		case heif_color_primaries_ITU_R_BT_470_6_System_B_G:
			primaries_name = "BT.470-6 System BG";
			break;
		case heif_color_primaries_ITU_R_BT_601_6:
			primaries_name = "BT.601";
			break;
		case heif_color_primaries_SMPTE_240M:
			primaries_name = "SMPTE 240M";
			break;
		case 8:
			primaries_name = "Generic film";
			break;
		case 9:
			primaries_name = "BT.2020";
			break;
		case 10:
			primaries_name = "XYZ";
			break;
		case 11:
			primaries_name = "SMPTE RP 431-2";
			break;
		case 12:
			primaries_name = "SMPTE EG 432-1 (DCI P3)";
			break;
		case 22:
			primaries_name = "EBU Tech. 3213-E";
			break;
		default:
			g_warning ("%s: Unsupported color_primaries value %d.",
					G_STRFUNC, nclx->color_primaries);
			return NULL;
		break;
	}

	switch (nclx->transfer_characteristics) {
		case heif_transfer_characteristic_ITU_R_BT_709_5:
			curve[0] = curve[1] = curve[2] = cmsBuildParametricToneCurve (NULL, 4,
											rec709_parameters);
			profile = cmsCreateRGBProfile (&whitepoint, &primaries, curve);
			cmsFreeToneCurve (curve[0]);
			trc_name = "Rec709 RGB";
			break;
		case heif_transfer_characteristic_ITU_R_BT_470_6_System_M:
			curve[0] = curve[1] = curve[2] = cmsBuildGamma (NULL, 2.2f);
			profile = cmsCreateRGBProfile (&whitepoint, &primaries, curve);
			cmsFreeToneCurve (curve[0]);
			trc_name = "Gamma2.2 RGB";
			break;
		case heif_transfer_characteristic_ITU_R_BT_470_6_System_B_G:
			curve[0] = curve[1] = curve[2] = cmsBuildGamma (NULL, 2.8f);
			profile = cmsCreateRGBProfile (&whitepoint, &primaries, curve);
			cmsFreeToneCurve (curve[0]);
			trc_name = "Gamma2.8 RGB";
			break;
		case heif_transfer_characteristic_linear:
			curve[0] = curve[1] = curve[2] = cmsBuildGamma (NULL, 1.0f);
			profile = cmsCreateRGBProfile (&whitepoint, &primaries, curve);
			cmsFreeToneCurve (curve[0]);
			trc_name = "linear RGB";
			break;
		case heif_transfer_characteristic_IEC_61966_2_1:
			/* same as default */
			curve[0] = curve[1] = curve[2] = cmsBuildParametricToneCurve (NULL, 4,
											srgb_parameters);
			profile = cmsCreateRGBProfile (&whitepoint, &primaries, curve);
			cmsFreeToneCurve (curve[0]);
			trc_name = "sRGB-TRC RGB";
			break;
		default:
			siril_log_color_message(_("Error: the specified NCLX TRC is not yet supported in Siril. "
										"A linear TRC will be used: you may need to fix this file "
										"up manually.\n"), "red");
			curve[0] = curve[1] = curve[2] = cmsBuildGamma (NULL, 1.0f);
			profile = cmsCreateRGBProfile (&whitepoint, &primaries, curve);
			cmsFreeToneCurve (curve[0]);
			trc_name = "linear RGB";
			break;
	}

	if (profile) {
		gchar *description = g_strdup_printf ("%s %s", primaries_name, trc_name);

		icc_profile_set_tag (profile, cmsSigProfileDescriptionTag,
											description);
		icc_profile_set_tag (profile, cmsSigDeviceMfgDescTag,
											"GIMP");
		icc_profile_set_tag (profile, cmsSigDeviceModelDescTag,
											description);
		icc_profile_set_tag (profile, cmsSigCopyrightTag,
											"Public Domain");
		g_free (description);
		return profile;
	}

	return NULL;
}

int readheif(const char* name, fits *fit, gboolean interactive){
	struct heif_error err;
#if LIBHEIF_HAVE_VERSION(1,13,0)
	heif_init(NULL);
#endif

	struct heif_context *ctx = heif_context_alloc();
	err = heif_context_read_from_file(ctx, name, NULL);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_context_free(ctx);
#if LIBHEIF_HAVE_VERSION(1,13,0)
		heif_deinit();
#endif
		return OPEN_IMAGE_ERROR;
	}

	// analyze image content
	int num = heif_context_get_number_of_top_level_images(ctx);
	if (num == 0) {
		siril_log_color_message(_("Input file contains no readable images.\n"), "red");
		heif_context_free(ctx);
#if LIBHEIF_HAVE_VERSION(1,13,0)
		heif_deinit();
#endif
		return OPEN_IMAGE_ERROR;
	}

	  // get the primary image

	heif_item_id primary;

	err = heif_context_get_primary_image_ID(ctx, &primary);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_context_free(ctx);
#if LIBHEIF_HAVE_VERSION(1,13,0)
		heif_deinit();
#endif
		return OPEN_IMAGE_ERROR;
	}

	// if primary image is no top level image or not present (invalid file), just take the first image

	if (!heif_context_is_top_level_image_ID(ctx, primary)) {
		int n = heif_context_get_list_of_top_level_image_IDs(ctx, &primary, 1);
		g_assert(n == 1);
	}

	heif_item_id selected_image = primary;

	if (num > 1) {
		if (!interactive) {
			siril_log_message(_("This is a sequence of %d images: "
					"loading the primary one.\n"), num);
		} else {
			if (!heif_dialog(ctx, &selected_image)) {
				heif_context_free(ctx);
#if LIBHEIF_HAVE_VERSION(1,13,0)
				heif_deinit();
#endif
				return OPEN_IMAGE_CANCEL;
			}
		}
	}

	// get a handle to the primary image
	struct heif_image_handle *handle;
	err = heif_context_get_image_handle(ctx, selected_image, &handle);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_context_free(ctx);
#if LIBHEIF_HAVE_VERSION(1,13,0)
		heif_deinit();
#endif
		return OPEN_IMAGE_ERROR;
	}
	uint8_t* icc_buffer = NULL;
	uint32_t icc_length = 0;
	// Get the ICC profile, if there is one
	enum heif_color_profile_type cp_type = heif_image_handle_get_color_profile_type(handle);
	const char* fourcc = cp_type == heif_color_profile_type_rICC ? "rICC" : cp_type == heif_color_profile_type_prof ? "prof" : cp_type == heif_color_profile_type_nclx ? "nclx" : "none";
	if (cp_type == heif_color_profile_type_rICC || cp_type == heif_color_profile_type_prof) {
		icc_length = heif_image_handle_get_raw_color_profile_size(handle);
		if (icc_length > 0) {
			icc_buffer = malloc(icc_length);
			err = heif_image_handle_get_raw_color_profile(handle, icc_buffer);
			if (err.code) {
				siril_log_color_message(_("Error getting ICC profile from HEIF file. Continuing: you will need to manually assign an ICC profile\n"), "red");
			}
		}
	} else if (cp_type == heif_color_profile_type_nclx) {
		struct heif_color_profile_nclx *nclx = NULL;
		err = heif_image_handle_get_nclx_color_profile(handle, &nclx);
		if (!err.code) {
			cmsHPROFILE profile = nclx_to_icc_profile(nclx);
			heif_nclx_color_profile_free(nclx);
			icc_buffer = get_icc_profile_data(profile, &icc_length);
			cmsCloseProfile(profile);
		}
	} else if (cp_type == heif_color_profile_type_not_present) {
		siril_debug_print("HEIF does not contain any color profile. Assuming sRGB.\n");
		icc_buffer = get_icc_profile_data(com.icc.srgb_profile, &icc_length);
	}
	siril_debug_print("ICC profile type %s, length %u read from HEIF file.\n", fourcc, icc_length);

	int has_alpha = heif_image_handle_has_alpha_channel(handle);
	int bit_depth;

	enum heif_chroma chroma = heif_chroma_interleaved_RGB;

	bit_depth = heif_image_handle_get_luma_bits_per_pixel (handle);
	if (bit_depth < 0) {
		siril_log_color_message(_("Input image has undefined bit-depth.\n"), "red");
		heif_image_handle_release (handle);
		heif_context_free (ctx);

		return OPEN_IMAGE_ERROR;
	} else {
		siril_debug_print("HEIF reports bit depth: %d has_alpha: %d\n", bit_depth, has_alpha);
	}

	if (bit_depth == 8) {
		if (has_alpha) {
			chroma = heif_chroma_interleaved_RGBA;
		} else {
			chroma = heif_chroma_interleaved_RGB;
		}
	} else { /* high bit depth */
#if LIBHEIF_HAVE_VERSION(1,8,0)
#if ( G_BYTE_ORDER == G_LITTLE_ENDIAN )
		if (has_alpha) {
			chroma = heif_chroma_interleaved_RRGGBBAA_LE;
		} else {
			chroma = heif_chroma_interleaved_RRGGBB_LE;
		}
#else
		if (has_alpha) {
			chroma = heif_chroma_interleaved_RRGGBBAA_BE;
		} else {
			chroma = heif_chroma_interleaved_RRGGBB_BE;
		}
#endif
#endif
	}

	// Note: the HEIF library does not appear to provide good grayscale support at present
	// so imported HEIF images will always be 3-channel
	struct heif_image *img = 0;
	err = heif_decode_image(handle, &img, heif_colorspace_RGB, chroma, NULL);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_image_handle_release(handle);
		heif_context_free(ctx);
#if LIBHEIF_HAVE_VERSION(1,13,0)
		heif_deinit();
#endif
		return OPEN_IMAGE_ERROR;
	}
	int stride;
	const uint8_t* udata = heif_image_get_plane_readonly(img, heif_channel_interleaved, &stride);
	const int width = heif_image_get_width(img, heif_channel_interleaved);
	const int height = heif_image_get_height(img, heif_channel_interleaved);

	size_t npixels = width * height;
	gboolean mono = TRUE;
	WORD *data = NULL;
	unsigned int nchannels = has_alpha ? 4 : 3;

	data = malloc(height * width * 3 * sizeof(WORD));
	if (!data) {
		PRINT_ALLOC_ERR;
		heif_image_handle_release(handle);
		heif_context_free(ctx);
#if LIBHEIF_HAVE_VERSION(1,13,0)
	heif_deinit();
#endif
		return OPEN_IMAGE_ERROR;
	}

	gboolean threaded = !get_thread_run();

	if (bit_depth > 8) { /* high bit depth */
		const float scale = bit_depth == 10 ? 1023.f : bit_depth == 12 ? 4095.f : 65535.f;
		WORD *buf[3] = { data, data + npixels, data + npixels * 2 };
		const uint16_t *src = (const uint16_t*) udata;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (threaded)
#endif
		for (int row = 0; row < height; row += stride) {
			int nrow = (row + stride > height ? height - row : stride);
			WORD r, g, b;
#ifdef _OPENMP
#pragma omp simd
#endif
			for (int i = 0; i < width * nrow; i++) {
				r = (int) ( ( (float) (0x0fff & (src[i * nchannels + RLAYER]))  / scale) * 65535.0f + 0.5f);
				g = (int) ( ( (float) (0x0fff & (src[i * nchannels + GLAYER]))  / scale) * 65535.0f + 0.5f);
				b = (int) ( ( (float) (0x0fff & (src[i * nchannels + BLAYER]))  / scale) * 65535.0f + 0.5f);
				*buf[RLAYER]++ = r;
				*buf[GLAYER]++ = g;
				*buf[BLAYER]++ = b;
				if (mono && !(r == g && r == b && g == b))
					mono = FALSE;
			}
		}
	} else {
		WORD *buf[3] = { data, data + npixels, data + npixels * 2 };
		unsigned int nchannels = has_alpha ? 4 : 3;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (threaded)
#endif
		for (int row = 0; row < height; row += stride) {
			int nrow = (row + stride > height ? height - row : stride);
			WORD r, g, b;
#ifdef _OPENMP
#pragma omp simd
#endif
			for (int i = 0; i < width * nrow; i++) {
				r = udata[i * nchannels + RLAYER];
				g = udata[i * nchannels + GLAYER];
				b = udata[i * nchannels + BLAYER];
				*buf[RLAYER]++ = r;
				*buf[GLAYER]++ = g;
				*buf[BLAYER]++ = b;
				if (mono && !(r == g && r == b && g == b))
					mono = FALSE;
			}
		}
	}

	if (mono) {
		WORD *tmp = realloc(data, npixels* sizeof(WORD));
		if (tmp) {
			data = tmp;
		} else {
			PRINT_ALLOC_ERR;
			return OPEN_IMAGE_ERROR;
		}
		nchannels = 1;
	} else {
		nchannels = 3;
	}

	clearfits(fit);

	set_all_keywords_default(fit);

	fit->bitpix = fit->orig_bitpix = bit_depth == 8 ? BYTE_IMG : USHORT_IMG;
	fit->type = DATA_USHORT;
	fit->naxis = nchannels == 1 ? 2 : 3;
	fit->rx = width;
	fit->ry = height;
	fit->naxes[0] = fit->rx;
	fit->naxes[1] = fit->ry;
	fit->naxes[2] = nchannels;
	fit->data = data;
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data + npixels * (nchannels == 3);
	fit->pdata[BLAYER] = fit->data + npixels * 2 * (nchannels == 3);
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	mirrorx(fit, FALSE);

	fits_initialize_icc(fit, icc_buffer, icc_length);
	color_manage(fit, (fit->icc_profile != NULL));
	free(icc_buffer);

	heif_image_handle_release(handle);
	heif_context_free(ctx);
	heif_image_release(img);
#if LIBHEIF_HAVE_VERSION(1,13,0)
	heif_deinit();
#endif
	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading HEIF: file %s, %ld layer(s), %ux%u pixels, bitdepth %d\n"),
			basename, fit->naxes[2], fit->rx, fit->ry, bit_depth);
	g_free(basename);

	return OPEN_IMAGE_OK;
}
#endif

#ifdef HAVE_LIBJXL

int readjxl(const char* name, fits *fit) {
	GError* error = NULL;
	gsize jxl_size;
	uint8_t* jxl_data = NULL;
	gboolean success = g_file_get_contents(name, (gchar**) &jxl_data, &jxl_size, &error);
	if (!success) {
		siril_log_color_message(_("Sorry but Siril cannot open the file: %s.\n"), "red", name);
		g_error_free(error);
		return OPEN_IMAGE_ERROR;
	}
	uint8_t* icc_profile = NULL;
	size_t icc_profile_length = 0;
	uint8_t* internal_icc_profile = NULL;
	size_t internal_icc_profile_length = 0;
	size_t xsize = 0, ysize = 0, zsize = 0, extra_channels = 0;
	uint8_t bitdepth = 0;
	float* pixels = NULL;
	if (DecodeJpegXlOneShotWrapper(jxl_data, jxl_size, &pixels,
			&xsize, &ysize, &zsize, &extra_channels, &bitdepth,
			&icc_profile, &icc_profile_length, &internal_icc_profile,
			&internal_icc_profile_length)) {
		siril_debug_print("Error while decoding the jxl file\n");
		return 1;
	}
	siril_debug_print("Image decoded as %d bits per pixel\n", bitdepth);
	clearfits(fit);

	set_all_keywords_default(fit);

	fit->bitpix = (bitdepth == 16 || com.pref.force_16bit) ? USHORT_IMG : FLOAT_IMG;
	fit->type = fit->bitpix == FLOAT_IMG ? DATA_FLOAT : DATA_USHORT;
	if (zsize == 1)
		fit->naxis = 2;
	else
		fit->naxis = 3;
	fit->rx = xsize;
	fit->ry = ysize;
	fit->naxes[0] = xsize;
	fit->naxes[1] = ysize;
	fit->naxes[2] = zsize;
	size_t npixels = xsize * ysize;
	if (fit->bitpix == FLOAT_IMG) {
		fit->fdata = malloc(xsize * ysize * zsize * sizeof(float));
		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = zsize == 3 ? fit->fdata + npixels : fit->fdata;
		fit->fpdata[BLAYER] = zsize == 3 ? fit->fdata + npixels * 2 : fit->fdata;
	} else {
		fit->data = malloc(xsize * ysize * zsize * sizeof(WORD));
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = zsize == 3 ? fit->data + npixels : fit->data;
		fit->pdata[BLAYER] = zsize == 3 ? fit->data + npixels * 2 : fit->data;
	}
	fit->keywords.binning_x = fit->keywords.binning_y = 1;

	if (fit->naxes[2] == 1) {
		if (fit->type == DATA_FLOAT) {
			memcpy(fit->fdata, pixels, xsize * ysize * zsize * sizeof(float));
		} else {
			for (int i = 0 ; i < xsize * ysize ; i++) {
				fit->data[i] = roundf_to_WORD(pixels[i] * USHRT_MAX);
			}
		}
	} else {
		const uint32_t totchans = zsize + extra_channels;
		siril_debug_print("Channels: total %u, extra %lu\n", totchans, extra_channels);
		if (fit->type == DATA_FLOAT) {
			for (size_t i = 0 ; i < xsize * ysize ; i++) {
				size_t pixel = i * totchans;
				for (int j = 0 ; j < 3 ; j++) {
					fit->fpdata[j][i] = pixels[pixel + j];
				}
			}
		} else {
			for (size_t i = 0 ; i < xsize * ysize ; i++) {
				size_t pixel = i * totchans;
				for (int j = 0 ; j < 3 ; j++) {
					fit->pdata[j][i] = roundf_to_WORD(pixels[pixel + j] * USHRT_MAX);
				}
			}
		}
	}
	free(pixels);
	cmsHPROFILE internal = cmsOpenProfileFromMem(internal_icc_profile, internal_icc_profile_length);
	cmsHPROFILE original = cmsOpenProfileFromMem(icc_profile, icc_profile_length);
	if (internal && original) {
		fit->icc_profile = copyICCProfile(internal);
		fit->color_managed = TRUE; // Don't use color_manage() here as we don't want the GUI updated yet
		gchar* orig_desc = siril_color_profile_get_description(original);
		gchar* int_desc = siril_color_profile_get_description(internal);
		siril_debug_print("Transforming from %s to %s\n", int_desc, orig_desc);
		g_free(orig_desc);
		g_free(int_desc);
		siril_colorspace_transform(fit, original);
		if (fit->icc_profile) cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = copyICCProfile(original);
	} else if (internal) {
		fit->icc_profile = copyICCProfile(internal);
	}
	color_manage(fit, (fit->icc_profile != NULL));
	if (original) cmsCloseProfile(original);
	if (internal) cmsCloseProfile(internal);
	free(icc_profile);
	free(internal_icc_profile);

	mirrorx(fit, FALSE);
//	fill_date_obs_if_any(fit, name);
	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading JPEG XL: file %s, %ld layer(s), %ux%u pixels\n"),
						basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return zsize;
}

int savejxl(const char *name, fits *fit, int effort, double quality, gboolean force_8bit) {
	gboolean threaded = !get_thread_run();

	char *filename = strdup(name);
	if (!g_str_has_suffix(filename, ".jxl")) {
		filename = str_append(&filename, ".jxl");
	}

	int max_bitdepth = force_8bit ? 8 : fit->type == DATA_USHORT ? 16 : 32;

	void *buffer = NULL;
	int bitdepth;
	if (interleave(fit, max_bitdepth, &buffer, &bitdepth, FALSE)) {
		siril_log_color_message(_("Error interleaving image\n"), "red");
	}

	uint32_t profile_len = 0;
	uint8_t *profile = NULL;
	cmsUInt32Number trans_type;
	cmsHTRANSFORM save_transform = NULL;

	if (fit->color_managed && fit->icc_profile) {
		// Check what is the appropriate color space to save in
		// Covers 8- and high-bitdepth files
		if (max_bitdepth == 8) { // 8-bit
			trans_type = (fit->naxes[2] == 1 ? TYPE_GRAY_8 : TYPE_RGB_8);
			if (fit->naxes[2] == 1) { // mono
				cmsHPROFILE srgb_mono_out = NULL;
				switch (com.pref.icc.export_8bit_method) {
					case EXPORT_SRGB:
						srgb_mono_out = gray_srgbtrc();
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(srgb_mono_out, &profile_len);
						cmsCloseProfile(srgb_mono_out);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_standard, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.mono_standard, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			} else { // rgb
				switch (com.pref.icc.export_8bit_method) {
					case EXPORT_SRGB:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_profile, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.srgb_profile, &profile_len);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.working_standard, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.working_standard, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			}
		} else { // high bitdepth
			if (fit->type == DATA_USHORT)
				trans_type = (fit->naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16);
			else // fit->type == DATA_FLOAT
				trans_type = (fit->naxes[2] == 1 ? TYPE_GRAY_FLT : TYPE_RGB_FLT);
			if (fit->naxes[2] == 1) { // mono
				cmsHPROFILE srgb_mono_out = NULL;
				switch (com.pref.icc.export_16bit_method) {
					case EXPORT_SRGB:
						srgb_mono_out = gray_srgbtrc();
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, srgb_mono_out, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(srgb_mono_out, &profile_len);
						cmsCloseProfile(srgb_mono_out);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_standard, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.mono_standard, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			} else { // rgb
				switch (com.pref.icc.export_16bit_method) {
					case EXPORT_SRGB:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_profile, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.srgb_profile, &profile_len);
						break;
					case EXPORT_WORKING:
						save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.working_standard, trans_type, com.pref.icc.export_intent, 0);
						profile = get_icc_profile_data(com.icc.working_standard, &profile_len);
						break;
					case EXPORT_IMAGE_ICC:
						profile = get_icc_profile_data(fit->icc_profile, &profile_len);
						break;
					default:
						free(filename);
						return 1;
				}
			}
		}
	} else if(com.pref.icc.default_to_srgb) {
		profile = get_icc_profile_data((fit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.srgb_profile), &profile_len);
	}

	cmsUInt32Number datasize = max_bitdepth == 8 ? 1 : fit->type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
	cmsUInt32Number bytesperline = fit->rx * datasize * fit->naxes[2];
	cmsUInt32Number bytesperplane = fit->rx * fit->ry * datasize * fit->naxes[2];
	if (save_transform) { // For "use image ICC profile" save_transform will be NULL, no need to transform the data
		cmsDoTransformLineStride(save_transform, buffer, buffer, fit->rx, fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
		cmsDeleteTransform(save_transform);
	}
	uint8_t *compressed = NULL;
	size_t compressed_length;

	if (quality == 100.0)
		siril_log_message(_("Saving JPEG XL: file %s, lossless, effort=%d %ld layer(s), %ux%u pixels, bit depth: %d\n"),
						filename, effort, fit->naxes[2], fit->rx, fit->ry, bitdepth);
	else
		siril_log_message(_("Saving JPEG XL: file %s, quality=%.3f, effort=%d %ld layer(s), %ux%u pixels, bit depth: %d\n"),
						filename, quality, effort, fit->naxes[2], fit->rx, fit->ry, bitdepth);

	EncodeJpegXlOneshotWrapper(buffer, fit->rx,
					fit->ry, fit->naxes[2], bitdepth,
					&compressed, &compressed_length, effort,
					quality, profile, profile_len);

	siril_log_color_message(_("Save complete.\n"), "green");
	GError *error = NULL;
	g_file_set_contents(name, (const gchar *) compressed, compressed_length, &error);
	free(buffer);
	free(compressed);
	free(filename);
	g_error_free(error);
	return OPEN_IMAGE_OK;
}

#endif
