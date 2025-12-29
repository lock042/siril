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

/* Management of Siril's internal image formats */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_statistics.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "core/masks.h"
#include "filters/mtf.h"
#include "io/sequence.h"
#include "io/fits_sequence.h"
#include "gui/progress_and_log.h"
#include "gui/siril_preview.h"
#include "algos/statistics.h"
#include "algos/demosaicing.h"
#include "algos/spcc.h"
#include "algos/siril_wcs.h"
#include "io/fits_keywords.h"
#include "image_format_fits.h"

#define RECIPSQRT2 0.70710678f // 1/sqrt(2) as float

const char *fit_extension[] = {
		".fit",
		".fits",
		".fts"
};

static char *EXPOSURE[] = { "EXPTIME", "EXPOSURE", NULL };
static char *NB_STACKED[] = { "STACKCNT", "NCOMBINE", NULL };

static int CompressionMethods[] = { RICE_1, GZIP_1, GZIP_2, HCOMPRESS_1};

#define __tryToFindKeywords(fptr, type, keywords, value, status) \
{ \
	int __iter__ = 0; \
	do { \
		*status = 0; \
		fits_read_key(fptr, type, keywords[__iter__], value, NULL, status); \
		__iter__++; \
	} while ((keywords[__iter__]) && (*status > 0)); \
}


void fit_get_photometry_data(fits *fit) {
	int status = 0;
	read_fits_date_obs_header(fit);
	__tryToFindKeywords(fit->fptr, TDOUBLE, EXPOSURE, &fit->keywords.exposure, &status);
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "AIRMASS", &fit->keywords.airmass, NULL, &status);
}

int fit_stats(fitsfile *fptr, float *mini, float *maxi) {
	int status = 0;
	int ii;
	long npixels = 1L;
	long anaxes[3] = { 1L, 1L, 1L }, firstpix[3] = { 1L, 1L, 1L };
	float *pix;
	float minval = 1.E33, maxval = -1.E33;

	/* initialize value in case where it does not work */
	*mini = 0;
	*maxi = 0;

	fits_get_img_size(fptr, 3, anaxes, &status);

	if (status) {
		report_fits_error(status); /* print error message */
		return(status);
	}

	npixels = anaxes[0];  /* no. of pixels to read in each row */
	pix = malloc(npixels * sizeof(float)); /* memory for 1 row */
	if (pix == NULL) {
		PRINT_ALLOC_ERR;
		return (1);
	}

	/* loop over all planes of the cube (2D images have 1 plane) */
	for (firstpix[2] = 1; firstpix[2] <= anaxes[2]; firstpix[2]++) {
		/* loop over all rows of the plane */
		for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
			/* give starting pixel coordinate and number of pixels to read */
			if (fits_read_pix(fptr, TFLOAT, firstpix, npixels, NULL, pix,
						NULL, &status))
				break; /* jump out of loop on error */

			for (ii = 0; ii < npixels; ii++) {
				if (pix[ii] < minval)
					minval = pix[ii]; /* find min and  */
				if (pix[ii] > maxval)
					maxval = pix[ii]; /* max values    */
			}
		}
	}    /* end of loop over planes */
	free(pix);

	if (status) {
		report_fits_error(status); /* print any error message */
	} else {
		siril_debug_print("  minimum value = %f\n", minval);
		siril_debug_print("  maximum value = %f\n", maxval);
		*maxi = maxval;
		*mini = minval;
	}
	return status;
}

static void read_history_in_hdu(fitsfile *fptr, GSList **list) {
	int nkeys, status = 0;
	char card[FLEN_CARD];
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	for (int i = 1; i <= nkeys; i++) {
		if (fits_read_record(fptr, i, card, &status))
			break;
		if (!strncmp(card, "HISTORY", 7)) {
			*list = g_slist_prepend(*list, g_strdup(card + 8));
		}
	}
}

/* copy the history into a list of strings
 * header is read from current HDU and the following HDU as long as they don't contain an image.
 * Original active HDU is restored */
static void fits_read_history(fitsfile *fptr, GSList **history) {
	GSList *list = NULL;
	read_history_in_hdu(fptr, &list); // read from current HDU

	// browse following HDU
	int orig_hdu, type, status;
	gboolean hdu_changed = FALSE;
	fits_get_hdu_num(fptr, &orig_hdu);
	do {
		status = 0;
		fits_movrel_hdu(fptr, 1, &type, &status);
		if (status) {
			//fits_report_error(stderr, status);
			break;
		}
		hdu_changed = TRUE;
		if (type == IMAGE_HDU)
			break;
		siril_debug_print("history read from another HDU (CHDU changed)\n");
		read_history_in_hdu(fptr, &list);
	} while (1);

	// restore
	if (hdu_changed) {
		status = 0;
		fits_movabs_hdu(fptr, orig_hdu, NULL, &status);
	}

	if (*history)
		g_slist_free_full(*history, g_free);
	list = g_slist_reverse(list);
	*history = list;
}

/* reading the FITS header to get useful information
 * stored in the fit, requires an opened file descriptor */
void read_fits_header(fits *fit) {
	/* use new keywords structure */
	read_fits_keywords(fit);

	/* so now, fill the wcslib structure. */
	load_WCS_from_fits(fit);

	fits_read_history(fit->fptr, &(fit->history));
}

GSList *read_header_keyvals_strings(fitsfile *fptr) {
	int nkeys, status = 0;
	if (fits_get_hdrspace(fptr, &nkeys, NULL, &status)) {
		report_fits_error(status);
		return NULL;
	}
	if (fits_read_record(fptr, 0, NULL, &status))
		return NULL;

	GSList *entries = NULL;
	for (int n = 1; n <= nkeys; n++) {
		char key[FLEN_KEYWORD], value[FLEN_VALUE], comment[FLEN_COMMENT];
		status = 0;
		if (fits_read_keyn(fptr, n, key, value, comment, &status)) {
			report_fits_error(status);
			break;
		}
		if (!strcmp(key, "COMMENT"))
			continue;
		int len = strlen(value);
		// pretty-print strings: remove quotes and trailing spaces
		if (len > 1 && value[0] == '\'' && value[len - 1] == '\'') {
			len -= 2;
			for (int i = 0; i < len; i++)
				value[i] = value[i + 1];
			value[len] = '\0';
			g_strchomp(value);
		}
		if (value[0] == '\0' && comment[0] == '\0')
			continue;
		header_record *r = malloc(sizeof(header_record));
		if (!r) {
			PRINT_ALLOC_ERR;
			break;
		}
		r->key = strdup(key);
		r->value = strdup(value[0] == '\0' ? comment : value);
		entries = g_slist_prepend(entries, r);
	}
	entries = g_slist_reverse(entries);
	return entries;
}

/* copy the header for the current HDU in a heap-allocated string */
static int copy_header_from_hdu(fitsfile *fptr, char **header, int *strsize, int *strlength) {
	int nkeys, status = 0;
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	if (status || nkeys < 0) {
		free(*header);
		*header = NULL;
		return 1;
	}
	for (int i = 1; i <= nkeys; i++) {
		int cardlen;
		char *newstr;
		char card[FLEN_CARD];
		if (fits_read_record(fptr, i, card, &status))
			break;
		cardlen = strlen(card);
		if (*strlength + cardlen + 1 >= *strsize) {
			*strsize += 567;
			newstr = realloc(*header, *strsize);
			if (!newstr) {
				PRINT_ALLOC_ERR;
				free(*header);
				*header = NULL;
				return 1;
			}
			*header = newstr;
		}
		strcpy(*header + *strlength, card);
		(*strlength) += cardlen;
		strcpy(*header + *strlength, "\n");
		(*strlength)++;
	}
	if (*strlength + 3 + 1 >= *strsize) {
		*strsize += 4;
		char *newstr = realloc(*header, *strsize);
		if (!newstr) {
			PRINT_ALLOC_ERR;
			free(*header);
			*header = NULL;
			return 1;
		}
		*header = newstr;
	}
	strcpy(*header + *strlength, "END");
	(*strlength) += 3;
	return 0;
}

/* copy the complete header in a heap-allocated string */
char *copy_header(fits *fit) {
	int strsize, strlength;
	char *header;

	/* each line in the FITS header is 80 character wide
	 * in our string, we also keep new lines, so that's 81.
	 * initial allocation is 20 times 81 = 1620
	 * reallocations are 7 lines, 567 */
	strsize = 1620;
	if (!(header = malloc(strsize))) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	strlength = 0;
	if (copy_header_from_hdu(fit->fptr, &header, &strsize, &strlength))
		return NULL;

	int orig_hdu, type, status;
	gboolean hdu_changed = FALSE;
	fits_get_hdu_num(fit->fptr, &orig_hdu);
	do {
		status = 0;
		fits_movrel_hdu(fit->fptr, 1, &type, &status);
		if (status)
			break;
		hdu_changed = TRUE;
		if (type == IMAGE_HDU)
			break;
		siril_debug_print("header read from another HDU (CHDU changed)\n");
		if (copy_header_from_hdu(fit->fptr, &header, &strsize, &strlength))
			break;
	} while (1);

	// restore
	if (hdu_changed) {
		status = 0;
		fits_movabs_hdu(fit->fptr, orig_hdu, NULL, &status);
	}

	if (header[0] == '\0') {
		free(header);
		header = NULL;
	}
	if (!header)
		return NULL;

	/* we need to test if text is utf8
	 * indeed some header are not.
	 *
	 */
	if (!g_utf8_validate(header, -1, NULL)) {
		gchar *str = g_utf8_make_valid(header, -1);
		free(header);
		header = strdup(str);
		g_free(str);
	}
	return header;
}

/***** data reading and transformation ******/

data_type get_data_type(int bitpix) {
	if (bitpix == BYTE_IMG || bitpix == SHORT_IMG || bitpix == USHORT_IMG)
		return DATA_USHORT;
	if (bitpix == LONGLONG_IMG)
		return DATA_UNSUPPORTED;
	return DATA_FLOAT;
}


void report_fits_error(int status) {
	if (status) {
		char errmsg[FLEN_ERRMSG];
		while (fits_read_errmsg(errmsg)) {
			siril_log_message(_("FITS error: %s\n"),
					errmsg[0] != '\0' ? errmsg : "unknown" );
		}
	}
}

static void conv_8_to_16(WORD *data, size_t nbdata) {
	for (size_t i = 0; i < nbdata; i++) {
		double tmp = (double) data[i] / UCHAR_MAX_DOUBLE * USHRT_MAX_DOUBLE;
		data[i] = round_to_WORD(tmp);
	}
}

static void conv_16_to_32(WORD *udata, float *fdata, size_t nbdata) {
	for (size_t i = 0; i < nbdata; i++) {
		fdata[i] = (double) udata[i] / USHRT_MAX_DOUBLE;
	}
}

/* convert FITS data formats to siril native.
 * nbdata is the number of pixels, w * h.
 * from is not freed, to must be allocated and can be the same as from */
static void convert_data_ushort(int bitpix, const void *from, WORD *to, size_t nbdata, gboolean values_above_1) {
	BYTE *data8;
	int16_t *data16;
	size_t i;

	switch (bitpix) {
		case BYTE_IMG:
			data8 = (BYTE *)from;
			for (i = 0; i < nbdata; i++)
				to[i] = (WORD)data8[i];
			break;
		case USHORT_IMG:	// siril 0.9 native
			// nothing to do
			break;
		case SHORT_IMG:
			// add 2^15 to the read data to obtain unsigned
			data16 = (int16_t *)from;
			for (i = 0; i < nbdata; i++) {
				int sum = 32768 + (int)data16[i];
				to[i] = (WORD)sum;
			}
			break;
		case LONGLONG_IMG:	// 64-bit integer pixels
		default:
			siril_log_message(_("Unknown FITS data format in internal conversion\n"));
	}
}

/* convert FITS data formats to siril native.
 * nbdata is the number of pixels, w * h.
 * from is not freed, to must be allocated and can be the same as from */
static void convert_data_float(int bitpix, const void *from, float *to, size_t nbdata, double data_max) {
	size_t i;
	BYTE *data8;
	WORD *ushort;
	int16_t *data16;
	double *pixels_double;
	long *sdata32;	// TO BE TESTED on 32-bit arch, seems to be a cfitsio bug
	unsigned long *data32;
	double mini = DBL_MAX;
	double maxi = -DBL_MAX;
	const float *data32f;

	switch (bitpix) {
		case BYTE_IMG:
			data8 = (BYTE *)from;
			for (i = 0; i < nbdata; i++)
				to[i] = (float)data8[i] * INV_UCHAR_MAX_SINGLE;
			break;
		case USHORT_IMG:	// siril 0.9 native
			ushort = (WORD *)from;
			for (i = 0; i < nbdata; i++) {
				to[i] = (float)ushort[i] * INV_USHRT_MAX_SINGLE;
			}
			break;
		case SHORT_IMG:
			// add 2^15 to the read data to obtain unsigned
			data16 = (int16_t *)from;
			for (i = 0; i < nbdata; i++) {
				int sum = 32768 + (int)data16[i];
				to[i] = (float)sum * INV_USHRT_MAX_SINGLE;
			}
			break;
		case ULONG_IMG:		// 32-bit unsigned integer pixels
			data32 = (unsigned long *)from;
			if (data_max > 0.0) {
				siril_debug_print("Normalizing image data with DATA_MAX of %d\n", round_to_int(data_max));
				for (i = 0; i < nbdata; i++)
					to[i] = (float)(((double)(data32[i])) / data_max);
			} else {
				// N.I.N.A. gives the normalization value in DATAMAX, if it's not
				// there, compute from data. This is not suitable for sequence work
				unsigned long min = ULONG_MAX, max = 0;
				for (i = 0; i < nbdata; i++) {
					min = min(data32[i], min);
					max = max(data32[i], max);
				}
				siril_log_message(_("Normalizing the image data with [%ld, %ld], not suitable for sequence operations\n"), min, max);
				mini = (double)min;
				maxi = (double)max;
				for (i = 0; i < nbdata; i++)
					to[i] = (float)((data32[i] - mini)) / (maxi - mini);
			}
			break;
		case LONG_IMG:		// 32-bit signed integer pixels
			sdata32 = (long *)from;
			if (data_max > 0.0) {
				siril_debug_print("Normalizing image data with DATA_MAX of %d\n", round_to_int(data_max));
				for (i = 0; i < nbdata; i++)
					to[i] = (float)(((double)(sdata32[i])) / data_max);
			} else {
				long min = LONG_MAX, max = LONG_MIN;
				for (i = 0; i < nbdata; i++) {
					min = min(sdata32[i], min);
					max = max(sdata32[i], max);
				}
				siril_log_message(_("Normalizing the image data with [%ld, %ld], not suitable for sequence operations\n"), min, max);
				mini = (double)min;
				maxi = (double)max;
				for (i = 0; i < nbdata; i++)
					to[i] = (float)((sdata32[i] - mini)) / (maxi - mini);
			}
			break;
		case FLOAT_IMG:		// 32-bit floating point pixels, we use it only if float is not in the [0, 1] range
			data32f = (float *)from;
			for (i = 0; i < nbdata; i++)
				to[i] = (data32f[i] / USHRT_MAX_SINGLE);
			break;
		case DOUBLE_IMG:	// 64-bit floating point pixels
			pixels_double = (double *)from;
			for (i = 0; i < nbdata; i++) {
				to[i] = (float)pixels_double[i];
			}
			break;
		case LONGLONG_IMG:	// 64-bit integer pixels
		default:
			siril_log_message(_("Unknown FITS data format in internal conversion\n"));
	}
}

static void convert_floats(int bitpix, float *data, size_t nbdata, gboolean verbose) {
	size_t i;
	switch (bitpix) {
		case BYTE_IMG:
			for (i = 0; i < nbdata; i++)
				data[i] = data[i] * INV_UCHAR_MAX_SINGLE;
			break;
		case FLOAT_IMG:
			if (verbose)
				siril_log_message(_("Normalizing input data to our float range [0, 1]\n")); // @suppress("No break at end of case")
				// Fallthrough is deliberate, this handles FITS with floating point data
				// where MAX is 65565 (e.g. JWST)
		default:
		case USHORT_IMG:	// siril 0.9 native
			for (i = 0; i < nbdata; i++)
				data[i] = data[i] * INV_USHRT_MAX_SINGLE;
			break;
		case SHORT_IMG:
			// add 2^15 to the read data to obtain unsigned
			for (i = 0; i < nbdata; i++) {
				data[i] = (32768.f + data[i]) * INV_USHRT_MAX_SINGLE;
			}
			break;
	}
}

static int get_compression_type(int siril_compression_fits_method) {
	if (siril_compression_fits_method >= G_N_ELEMENTS(CompressionMethods)) {
		return -1;
	} else {
		return CompressionMethods[siril_compression_fits_method];
	}
}

/* Move to the first HDU of an opened FIT file
 * with IMAGE type and with dimension > 0
 */
static int siril_fits_move_first_image(fitsfile* fp) {
	int status = 0;
	int naxis = 0;

	/* Move to the first HDU of type image */
	fits_movabs_hdu(fp, 1, IMAGE_HDU, &status);
	if (status) {
		siril_log_message(_("Cannot move to first image in FITS file\n"));
		report_fits_error(status); /* print error message */
		return status;
	}

	/* Find the first HDU of type image with dimension > 0 */
	do {
		//Check naxis for current image HDU
		fits_get_img_dim(fp, &naxis, &status);
		if (status) {
			siril_log_message(_("Cannot get dimension of FITS file\n"));
			report_fits_error(status); /* print error message */
			return status;
		}
		if (naxis > 0) {
			break;
		}
		//Check next image HDU
		fits_movrel_hdu(fp, 1, IMAGE_HDU, &status);
		if (status) {
			siril_log_message(_("Cannot move to next image in FITS file\n"));
			report_fits_error(status); /* print error message */
			return status;
		}
	} while (!status);

	//siril_debug_print("Found image HDU (changed CHDU) with naxis %d (status %d)\n", naxis, status);
	return status;
}

/* read buffer from an already open FITS file, fit should have all metadata
 * correct, and convert the buffer to fit->data with the given type, which
 * currently should be TBYTE or TUSHORT because fit doesn't contain other data.
 * filename is for error reporting
 */
int read_fits_with_convert(fits* fit, const char* filename, gboolean force_float) {
	int status = 0, zero = 0, datatype;
	BYTE *data8;
	unsigned long *pixels_long;
	// orig ^ gives the coordinate in each dimension of the first pixel to be read
	size_t nbpix = fit->naxes[0] * fit->naxes[1];
	size_t nbdata = nbpix * fit->naxes[2];
	// with force_float, image is read as float data, type is stored as DATA_FLOAT
	int fake_bitpix = force_float ? FLOAT_IMG : fit->bitpix;
	double data_max;

	switch (fake_bitpix) {
		case BYTE_IMG:
		case SHORT_IMG:
		case USHORT_IMG:
			/* we store these types as unsigned short */
			if ((fit->data = malloc(nbdata * sizeof(WORD))) == NULL) {
				PRINT_ALLOC_ERR;
				return -1;
			}
			fit->pdata[RLAYER] = fit->data;
			if (fit->naxis == 3) {
				fit->pdata[GLAYER] = fit->data + nbpix;
				fit->pdata[BLAYER] = fit->data + nbpix * 2;
			} else {
				fit->pdata[GLAYER] = fit->data;
				fit->pdata[BLAYER] = fit->data;
			}
			fit->type = DATA_USHORT;
			break;

		case ULONG_IMG:		// 32-bit unsigned integer pixels
		case LONG_IMG:		// 32-bit signed integer pixels
		case DOUBLE_IMG:	// 64-bit floating point pixels
		case FLOAT_IMG:		// 32-bit floating point pixels
			/* we store these types as float */
			if ((fit->fdata = malloc(nbdata * sizeof(float))) == NULL) {
				PRINT_ALLOC_ERR;
				return -1;
			}
			fit->fpdata[RLAYER] = fit->fdata;
			if (fit->naxis == 3) {
				fit->fpdata[GLAYER] = fit->fdata + nbpix;
				fit->fpdata[BLAYER] = fit->fdata + nbpix * 2;
			} else {
				fit->fpdata[GLAYER] = fit->fdata;
				fit->fpdata[BLAYER] = fit->fdata;
			}
			fit->type = DATA_FLOAT;
			break;

		case LONGLONG_IMG:	// 64-bit integer pixels
		default:
			siril_log_message(_("FITS image format %d is not supported by Siril.\n"), fit->bitpix);
			return -1;
	}

	status = 0;
	switch (fake_bitpix) {
	case BYTE_IMG:
		data8 = malloc(nbdata * sizeof(BYTE));
		datatype = fit->bitpix == BYTE_IMG ? TBYTE : TSBYTE;
		fits_read_img(fit->fptr, datatype, 1, nbdata, &zero, data8, &zero, &status);
		if (status) break;
		convert_data_ushort(fit->bitpix, data8, fit->data, nbdata, FALSE);
		free(data8);
		break;
	case SHORT_IMG:
		fits_read_img(fit->fptr, TSHORT, 1, nbdata, &zero, fit->data, &zero, &status);
		if (status) break;
		convert_data_ushort(fit->bitpix, fit->data, fit->data, nbdata, FALSE);
		fit->bitpix = USHORT_IMG;
		break;
	case USHORT_IMG:
		// siril 0.9 native, no conversion required
		fits_read_img(fit->fptr, TUSHORT, 1, nbdata, &zero, fit->data, &zero, &status);
		if (status == NUM_OVERFLOW) {
			// in case there are errors, we try short data
			status = 0;
			fits_read_img(fit->fptr, TSHORT, 1, nbdata, &zero, fit->data,
					&zero, &status);
			if (status)
				break;
			convert_data_ushort(SHORT_IMG, fit->data, fit->data, nbdata, FALSE);
			if (fit->keywords.lo)
				fit->keywords.lo += 32768;
			if (fit->keywords.hi)
				fit->keywords.hi += 32768;
			fit->bitpix = USHORT_IMG;
		}
		break;

	case ULONG_IMG:		// 32-bit unsigned integer pixels
	case LONG_IMG:		// 32-bit signed integer pixels
		pixels_long = malloc(nbdata * sizeof(unsigned long));
		datatype = fit->bitpix == LONG_IMG ? TLONG : TULONG;
		fits_read_img(fit->fptr, datatype, 1, nbdata, &zero, pixels_long, &zero, &status);
		if (status) break;
		fits_read_key(fit->fptr, TDOUBLE, "DATAMAX", &data_max, NULL, &status);
		if (status) {
			data_max = 0.0;
			status = 0;
		}
		convert_data_float(fit->bitpix, pixels_long, fit->fdata, nbdata, data_max);
		free(pixels_long);
		fit->bitpix = FLOAT_IMG;
		break;
	case FLOAT_IMG:		// 32-bit floating point pixels
		// siril 1.0 native, no conversion required
	case DOUBLE_IMG:	// 64-bit floating point pixels
		// let cfitsio do the conversion
		/* we assume we are in the range [0, 1]. But, for some images
		 * some values can be negative
		 */
		fits_read_img(fit->fptr, TFLOAT, 1, nbdata, &zero, fit->fdata, &zero, &status);
		if ((fit->bitpix == USHORT_IMG || fit->bitpix == SHORT_IMG
				// needed for some FLOAT_IMG. 10.0 is probably a good number to represent the limit at which we judge that these are not clip-on pixels.
				|| fit->bitpix == BYTE_IMG) || fit->keywords.data_max > 10.0) {
			convert_floats(fit->bitpix, fit->fdata, nbdata, TRUE);
		}
		fit->bitpix = FLOAT_IMG;
		fit->orig_bitpix = FLOAT_IMG; // force this, to avoid problems saving the FITS if needed
		break;
	}

	if (status) {
		siril_log_message(_("Fitsio error reading data, file: %s.\n"), filename);
		report_fits_error(status);
		return -1;
	}

	return 0;
}

/* This function reads partial data on one layer from the opened FITS and
 * convert it to siril's format (USHORT) */
int internal_read_partial_fits(fitsfile *fptr, unsigned int ry,
		int bitpix, void *dest, int layer, const rectangle *area) {
	double data_max = -1.0;
	int datatype;
	BYTE *data8;
	long *pixels_long;
	long fpixel[3], lpixel[3], inc[3] = { 1L, 1L, 1L };
	int zero = 0, status = 0;

	/* fpixel is first pixel, lpixel is last pixel, starts with value 1 */
	fpixel[0] = area->x + 1;        // in siril, it starts with 0
	fpixel[1] = ry - area->y - area->h + 1;
	fpixel[2] = layer + 1;
	lpixel[0] = area->x + area->w;  // with w and h at least 1, we're ok
	lpixel[1] = ry - area->y;
	lpixel[2] = layer + 1;

	size_t nbdata = area->w * area->h;

	switch (bitpix) {
		case BYTE_IMG:
			data8 = malloc(nbdata * sizeof(BYTE));
			datatype = TBYTE;
			fits_read_subset(fptr, datatype, fpixel, lpixel, inc, &zero, data8,
					&zero, &status);
			if (status) break;
			convert_data_ushort(bitpix, data8, dest, nbdata, FALSE);
			free(data8);
			break;
		case SHORT_IMG:
			fits_read_subset(fptr, TSHORT, fpixel, lpixel, inc, &zero, dest,
					&zero, &status);
			convert_data_ushort(bitpix, dest, dest, nbdata, FALSE);
			break;
		case USHORT_IMG:
			fits_read_subset(fptr, TUSHORT, fpixel, lpixel, inc, &zero, dest,
					&zero, &status);
			break;

		/* types below are read as a 32-bit float image: this breaks compatibility
		 * with functions working only with 16-bit int data */
		case ULONG_IMG:		// 32-bit unsigned integer pixels
		case LONG_IMG:		// 32-bit signed integer pixels
			pixels_long = malloc(nbdata * sizeof(long));
			status = 0;
			datatype = bitpix == LONG_IMG ? TLONG : TULONG;
			fits_read_subset(fptr, datatype, fpixel, lpixel, inc, &zero,
					pixels_long, &zero, &status);
			if (status) break;
			fits_read_key(fptr, TDOUBLE, "DATAMAX", &data_max, NULL, &status);
			if (status) {
				data_max = 0.0;
				status = 0;
			}
			convert_data_float(bitpix, pixels_long, dest, nbdata, data_max);
			free(pixels_long);
			break;
		case DOUBLE_IMG:	// 64-bit floating point pixels
		case FLOAT_IMG:		// 32-bit floating point pixels
			fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &zero, dest, &zero, &status);
			if (status) break;
			int status2 = 0;
			fits_read_key(fptr, TDOUBLE, "DATAMAX", &data_max, NULL, &status2);
			if (status2 == KEY_NO_EXIST && nbdata > 3) {
				data_max = 0.;
				for (size_t i = 0; i < nbdata; i += nbdata / 3) { // we check the 3 pixels in diagonal and hope for the best not all will be null...
					data_max = max(data_max, ((float *)dest)[i]);
				}
				status2 = 0;
			}
			if (status2 == 0 && data_max > 10.0) { // needed for some FLOAT_IMG
				convert_floats(bitpix, dest, nbdata, FALSE);
			}
			break;
		case LONGLONG_IMG:	// 64-bit integer pixels
		default:
			siril_log_message(_("FITS image format %d is not supported by Siril.\n"), bitpix);
			return -1;
	}
	return status;
}

int siril_fits_create_diskfile(fitsfile **fptr, const char *filename, int *status) {
	*status = 0;
	gchar *dirname = g_path_get_dirname(filename);
	GDir *dir = g_dir_open(dirname, 0, NULL);
	if (!dir && siril_mkdir_with_parents(dirname, 0755) < 0) {
		*status = 1;
		g_free(dirname);
		return *status;
	}
	if (dir)
		g_dir_close(dir);
	g_free(dirname);
	fits_create_diskfile(fptr, filename, status);
	return *status;
}

// updates the header string from fit->header
// by creating an in-memory fits file (as oppsed to on-disk file)
void update_fits_header(fits *fit) {
	void *memptr;
	size_t memsize = IOBUFLEN;
	int status = 0;
	fitsfile *fptr = NULL;
	memptr = malloc(memsize);
	if (!memptr) {
		PRINT_ALLOC_ERR;
		return;
	}
	fits_create_memfile(&fptr, &memptr, &memsize, IOBUFLEN, realloc, &status);
	if (status) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return;
	}
	if (fits_create_img(fptr, fit->bitpix, fit->naxis, fit->naxes, &status)) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return;
	}
	fits tmpfit = { 0 };
	copy_fits_metadata(fit, &tmpfit);
	tmpfit.fptr = fptr;
	save_fits_header(&tmpfit);
	if (fit->header)
		free(fit->header);
	fit->header = copy_header(&tmpfit);
	fits_close_file(fptr, &status);
	clearfits(&tmpfit);
	free(memptr);
}

void save_fits_header(fits *fit) {
	save_fits_keywords(fit);
	save_wcs_keywords(fit);
	save_history_keywords(fit);
	save_fits_unknown_keywords(fit);
}

/********************** public functions ************************************/

/* gets specific data from the header (for stacking).
 * if exposure is not found, it defaults to 0
 * if livetime is not found, exposure is assigned to it
 * if stack_count is not found, it defaults to 1
 */
void get_date_data_from_fitsfile(fitsfile *fptr, GDateTime **dt, double *exposure, double *livetime, unsigned int *stack_count) {
	*exposure = 0.0;
	*stack_count = 1;
	*dt = NULL;
	int status = 0;

	__tryToFindKeywords(fptr, TDOUBLE, EXPOSURE, exposure, &status);
	status = 0;
	__tryToFindKeywords(fptr, TUINT, NB_STACKED, stack_count, &status);
	status = 0;
	if (fits_read_key(fptr, TDOUBLE, "LIVETIME", livetime, NULL, &status))
		*livetime = *exposure;

	char date_obs[FLEN_VALUE];
	status = 0;
	if (!fits_read_key(fptr, TSTRING, "DATE-OBS", &date_obs, NULL, &status))
		*dt = FITS_date_to_date_time(date_obs);
}

int import_metadata_from_fitsfile(fitsfile *fptr, fits *to) {
	fits from = { 0 };
	from.fptr = fptr;
	read_fits_header(&from);
	copy_fits_metadata(&from, to);
	clearfits(&from);
	return 0;
}

int write_icc_profile_to_fptr(fitsfile *fptr, cmsHPROFILE icc_profile) {
	int status = 0;     // CFITSIO status variable
	cmsUInt32Number profile_length;
	cmsUInt8Number *profile = NULL;
	status = cmsSaveProfileToMem(icc_profile, NULL, &profile_length);
	if (profile_length > 0) {
		profile = malloc(profile_length * sizeof(cmsUInt8Number));
		cmsBool ret = cmsSaveProfileToMem(icc_profile, (void*) profile, &profile_length);
		status = !ret;
	}

	// Move to the last HDU
	int nhdus;
	fits_get_num_hdus(fptr, &nhdus, &status);

	if (nhdus && fits_movabs_hdu(fptr, nhdus, NULL, &status)) {
		fits_report_error(stderr, status);
		goto ERROR_MESSAGE_AND_RETURN;
	}

	// Create the image extension
	long naxis = 1, naxes = profile_length;  // Image dimensions
	long fpixel = 1;       // First pixel to write
	if (fits_create_img(fptr, BYTE_IMG, naxis, &naxes, &status)) {
		fits_report_error(stderr, status);
		goto ERROR_MESSAGE_AND_RETURN;
	}
	int bitpix = 8, pcount = 0, gcount = 1;
	// Populate the header with the field values
	status = (fits_update_key(fptr, TSTRING, "XTENSION", "IMAGE", NULL, &status));
	if (!status) status = (fits_update_key(fptr, TINT, "BITPIX", &bitpix, NULL, &status));
	if (!status) status = (fits_update_key(fptr, TINT, "NAXIS", &naxis, NULL, &status));
	if (!status) status = (fits_update_key(fptr, TINT, "NAXIS1", &profile_length, NULL, &status));
	if (!status) status = (fits_update_key(fptr, TINT, "PCOUNT", &pcount, NULL, &status));
	if (!status) status = (fits_update_key(fptr, TINT, "GCOUNT", &gcount, NULL, &status));
	if (!status) status = (fits_update_key(fptr, TSTRING, "EXTNAME", "ICCProfile", NULL, &status));
	if (status) {
		fits_report_error(stderr, status);
		goto ERROR_MESSAGE_AND_RETURN;
	}

	// Write the ICC profile data to the image extension
	if (fits_write_img(fptr, TBYTE, fpixel, profile_length, profile, &status)) {
		fits_report_error(stderr, status);
		goto ERROR_MESSAGE_AND_RETURN;
	}

	// Success! Clean up and return
	g_free(profile);
	siril_debug_print("ICC profile embedded in FITS file\n");
	return 0;

ERROR_MESSAGE_AND_RETURN:
	free(profile);
	siril_log_color_message(_("Warning: error encountered writing ICC profile to FITS.\n"), "salmon");
	return status;
}

int write_icc_profile_to_fits(fits *fit) {
	int retval = write_icc_profile_to_fptr(fit->fptr, fit->icc_profile);
	return retval;
}

/* Look for a HDU containing an ICC profile; if one is found, open it */
cmsHPROFILE read_icc_profile_from_fptr(fitsfile *fptr) {
	cmsHPROFILE icc_profile;
	int status = 0;
	char extname[FLEN_VALUE], comment[FLEN_COMMENT];
	int ihdu, nhdus, hdutype, orig_hdu = 1;
	fits_get_hdu_num(fptr, &orig_hdu);
	fits_get_num_hdus(fptr, &nhdus, &status);
	for (ihdu = 2 ; ihdu <= nhdus ; ihdu++) {
		fits_movabs_hdu(fptr,ihdu, &hdutype, &status);
		fits_read_key(fptr, TSTRING, "EXTNAME", &extname, comment, &status);
		if (status) {
			status = 0;
			continue; /* next HDU */
		}
		if (!g_str_has_prefix(extname, "ICCProfile"))
			continue; /* next HDU */
		break; /* current HDU matches */
	}
	if (ihdu > nhdus) {
		/* no matching HDU */
		status = BAD_HDU_NUM;
		fits_movabs_hdu(fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return NULL;
	}
	int strsize = 1620;
	int strlength = 0;
	char *header = NULL;
	if (!(header = malloc(strsize))) {
		PRINT_ALLOC_ERR;
		fits_movabs_hdu(fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return NULL;
	}
	status = copy_header_from_hdu(fptr, &header, &strsize, &strlength);
	if (status) {
		free(header);
		fits_movabs_hdu(fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return NULL;
	}
	// Get the ICC Profile length
	uint32_t profile_length, bitpix;
	fits_read_key(fptr, TUINT, "NAXIS1", &profile_length, comment, &status);
	fits_read_key(fptr, TUINT, "BITPIX", &bitpix, comment, &status);
	if (bitpix != 8 || status != 0) {
		free(header);
		fits_movabs_hdu(fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return NULL;
	}
	int zero = 0;
	BYTE *profile = NULL;
	if (!(profile = malloc(profile_length * sizeof(BYTE)))) {
		PRINT_ALLOC_ERR;
		free(header);
		fits_movabs_hdu(fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return NULL;
	}
	fits_read_img(fptr, TBYTE, 1, profile_length, &zero, profile, &zero, &status);
	if (status) {
		free(profile);
		free(header);
		fits_movabs_hdu(fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return NULL;
	}
	icc_profile = cmsOpenProfileFromMem(profile, profile_length);
	if (icc_profile)
		siril_debug_print("Embedded ICC profile read from FITS\n");
	free(profile);
	free(header);
	fits_movabs_hdu(fptr, orig_hdu, &hdutype, &status);
	if (status)
		siril_debug_print("Error returning to original HDU!\n");
	return icc_profile;
}

int read_icc_profile_from_fits(fits *fit) {
	int status = 0;
	char extname[FLEN_VALUE], comment[FLEN_COMMENT];
	int ihdu, nhdus, hdutype, orig_hdu = 1;
	fits_get_hdu_num(fit->fptr, &orig_hdu);
	// siril_debug_print("Original HDU before looking for ICC profile: %d\n", orig_hdu);
	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = NULL;
	fits_get_num_hdus(fit->fptr, &nhdus, &status);
	for (ihdu = 2 ; ihdu <= nhdus ; ihdu++) {
		fits_movabs_hdu(fit->fptr,ihdu, &hdutype, &status);
		fits_read_key(fit->fptr, TSTRING, "EXTNAME", &extname, comment, &status);
		if (status) {
			status = 0;
			continue; /* next HDU */
		}
		if (!g_str_has_prefix(extname, "ICCProfile"))
			continue; /* next HDU */
		break; /* current HDU matches */
	}
	if (ihdu > nhdus) {
		/* no matching HDU */
		fits_movabs_hdu(fit->fptr, orig_hdu, &hdutype, &status);
		status = BAD_HDU_NUM;
		return 1;
	}
	int strsize = 1620;
	int strlength = 0;
	char *header = NULL;
	if (!(header = malloc(strsize))) {
		PRINT_ALLOC_ERR;
		fits_movabs_hdu(fit->fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return 1;
	}
	status = copy_header_from_hdu(fit->fptr, &header, &strsize, &strlength);
	if (status) {
		free(header);
		fits_movabs_hdu(fit->fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return 1;
	}
	// Get the ICC Profile length
	uint32_t profile_length, bitpix;
	fits_read_key(fit->fptr, TUINT, "NAXIS1", &profile_length, comment, &status);
	fits_read_key(fit->fptr, TUINT, "BITPIX", &bitpix, comment, &status);
	if (bitpix != 8 || status != 0) {
		free(header);
		fits_movabs_hdu(fit->fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return 1;
	}
	int zero = 0;
	BYTE *profile = NULL;
	if (!(profile = malloc(profile_length * sizeof(BYTE)))) {
		PRINT_ALLOC_ERR;
		free(header);
		fits_movabs_hdu(fit->fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return 1;
	}
	fits_read_img(fit->fptr, TBYTE, 1, profile_length, &zero, profile, &zero, &status);
	if (status) {
		free(profile);
		free(header);
		fits_movabs_hdu(fit->fptr, orig_hdu, &hdutype, &status);
		if (status)
			siril_debug_print("Error returning to original HDU!\n");
		return 1;
	}
	fit->icc_profile = cmsOpenProfileFromMem(profile, profile_length);
	if (fit->icc_profile) {
		siril_debug_print("Embedded ICC profile read from FITS\n");
		color_manage(fit, TRUE);
	} else {
		color_manage(fit, FALSE);
	}
	free(profile);
	free(header);
	fits_movabs_hdu(fit->fptr, orig_hdu, &hdutype, &status);
	if (status)
		siril_debug_print("Error returning to original HDU!\n");
	fits_get_hdu_num(fit->fptr, &orig_hdu);
	return 0;
}

/* from bitpix, depending on BZERO, bitpix and orig_bitpix are set.
 *
 * since USHORT_IMG is a cfitsio trick and doesn't really exist in the
 * file, when we read it, the bitpix is given as SHORT_IMG. If BZERO is
 * 2^15, it means that it is USHORT data and we force the bitpix to it
 * in order to read it properly. Same thing for LONG and ULONG.
 * https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node23.html
 */
void manage_bitpix(fitsfile *fptr, int *bitpix, int *orig_bitpix) {
	double offset;
	int status = 0;
	fits_read_key(fptr, TDOUBLE, "BZERO", &offset, NULL, &status);
	if (!status) {
		if (*bitpix == SHORT_IMG && offset != 0.0) {
			*bitpix = USHORT_IMG;
		}
		else if (*bitpix == LONG_IMG && offset != 0.0)
			*bitpix = ULONG_IMG;
	} else {
		/* but some software just put unsigned 16-bit data in the file
		 * and don't set the BZERO keyword... */
		if (status == KEY_NO_EXIST && *bitpix == SHORT_IMG)
			*bitpix = USHORT_IMG;
	}
	// and we store the original bitpix to reuse it during later partial
	// reads when we have no access to the header or a fits * struct.
	*orig_bitpix = *bitpix;
}

// return 0 on success, fills realname if not NULL with the opened file's name
int readfits(const char *filename, fits *fit, char *realname, gboolean force_float) {
	int status, retval = 1, dataok = 0, hduok = 0;
	char *name = NULL;
	unsigned long datasum, hdusum;
	gchar *basename;
	image_type imagetype;

	if (stat_file(filename, &imagetype, &name)) {
		siril_log_message(_("%s.[any_allowed_extension] not found.\n"),
				filename);
		free(name);
		return 1;
	}
	if (imagetype != TYPEFITS) {
		siril_log_message(_("The file %s is not a FITS file or doesn't exists with FITS extensions.\n"), filename);
		free(name);
		return 1;
	}

	if (realname)
		strcpy(realname, name);

	status = 0;
	siril_fits_open_diskfile_img(&(fit->fptr), name, READONLY, &status);
	if (status) {
		report_fits_error(status);
		free(name);
		return status;
	}
	free(name);

	status = read_fits_metadata(fit);
	if (status)
		goto close_readfits;

	if (com.pref.use_checksum) {
		fits_verify_chksum(fit->fptr, &dataok, &hduok, &status);
		if (hduok == -1 || dataok == -1) {
			status = 0;
			fits_get_chksum(fit->fptr, &datasum, &hdusum, &status);
			if (hduok == -1) {
				char checksum[FLEN_VALUE], ascii[FLEN_VALUE];
				fits_read_key(fit->fptr, TSTRING, "CHECKSUM", &checksum, NULL, &status);
				fits_encode_chksum(hdusum, TRUE, ascii);
				siril_log_color_message(_("Error: HDU checksum mismatch. Expected %s, got %s.\n"), "red", ascii, checksum);
			}
			if (dataok == -1) {
				char checksum[FLEN_VALUE];
				fits_read_key(fit->fptr, TSTRING, "DATASUM", &checksum, NULL, &status);
				siril_log_color_message(_("Error: Data checksum mismatch. Expected %lu, got %s.\n"), "red", datasum, checksum);
			}
		}
	}

	retval = read_fits_with_convert(fit, filename, force_float);
	fit->top_down = FALSE;

	if (!retval) {
		basename = g_path_get_basename(filename);
		siril_log_message(_("Reading FITS: file %s, %ld layer(s), %ux%u pixels, %d bits\n"),
				basename, fit->naxes[2], fit->rx, fit->ry, fit->type == DATA_USHORT ? 16 : 32) ;
		g_free(basename);
	}
	check_profile_correct(fit);

close_readfits:
	status = 0;
	fits_close_file(fit->fptr, &status);
	return retval;
}

static int siril_fits_open_diskfile(fitsfile **fptr, const char *filename, int iomode, int *status) {
	fits_open_diskfile(fptr, filename, iomode, status);
	return *status;
}

int siril_fits_open_diskfile_img(fitsfile **fptr, const char *filename, int iomode, int *status) {
	fits_open_diskfile(fptr, filename, iomode, status);
	if (!(*status)) {
		*status = siril_fits_move_first_image(*fptr);
	}
	return *status;
}

GDateTime* get_date_from_fits(const gchar *filename) {
	fitsfile *fptr = NULL;
	GDateTime *date = NULL;
	int status = 0;
	fits_open_diskfile(&fptr, filename, READONLY, &status);
	if (!status) {
		status = siril_fits_move_first_image(fptr);
		char date_obs[FLEN_VALUE] = { 0 };
		fits_read_key(fptr, TSTRING, "DATE-OBS", &date_obs, NULL, &status);
		if (!status)
			date = FITS_date_to_date_time(date_obs);
	}
	status = 0;
	fits_close_file(fptr, &status);
	return date;
}

// reset a fit data structure, deallocates everything in it but keep the data:
// useful in processing internal_fits in SEQ_INTERNAL sequences
void clearfits_header(fits *fit) {
	if (fit == NULL)
		return;
	if (fit->header) {
		free(fit->header);
		fit->header = NULL;
	}
	if (fit->unknown_keys) {
		g_free(fit->unknown_keys);
		fit->unknown_keys = NULL;
	}
	if (fit->history) {
		g_slist_free_full(fit->history, g_free);
		fit->history = NULL;
	}
	if (fit->keywords.date_obs) {
		g_date_time_unref(fit->keywords.date_obs);
		fit->keywords.date_obs = NULL;
	}
	if (fit->keywords.date) {
		g_date_time_unref(fit->keywords.date);
		fit->keywords.date = NULL;
	}
	if (fit->stats) {
		for (int i = 0; i < fit->naxes[2]; i++)
			free_stats(fit->stats[i]);
		free(fit->stats);
		fit->stats = NULL;
	}
	color_manage(fit, FALSE);
	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = NULL;
	free_wcs(fit);
	reset_wcsdata(fit);
	if (fit == gfit && is_preview_active())
		clear_backup();
	memset(fit, 0, sizeof(fits));
}

// reset a fit data structure, deallocates everything in it and zero the data
void clearfits(fits *fit) {
	if (fit == NULL)
		return;
	if (fit->data) {
		free(fit->data);
		fit->data = NULL;
	}
	if (fit->fdata) {
		free(fit->fdata);
		fit->fdata = NULL;
	}
	if (fit->mask) {
		free_mask(fit->mask);
		fit->mask = NULL;
	}
	clearfits_header(fit);
}

/* Read a rectangular section of a FITS image in Siril's format, pointed by its
 * exact filename. Only layer layer is read.
 * Returned fit->data is upside-down. */
int readfits_partial(const char *filename, int layer, fits *fit,
		const rectangle *area, gboolean do_photometry) {
	int status;
	size_t nbdata;
	double data_max = 0.0;

	status = 0;
	if (siril_fits_open_diskfile_img(&(fit->fptr), filename, READONLY, &status)) {
		report_fits_error(status);
		return status;
	}

	status = 0;
	fits_get_img_param(fit->fptr, 3, &(fit->bitpix), &(fit->naxis), fit->naxes,
			&status);
	if (status) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return status;
	}

	manage_bitpix(fit->fptr, &(fit->bitpix), &(fit->orig_bitpix));

	if (do_photometry)
		fit_get_photometry_data(fit);

	if (fit->naxis == 2 && fit->naxes[2] == 0)
		fit->naxes[2] = 1;	// see readfits for the explanation
	if (layer > fit->naxes[2] - 1) {
		siril_log_message(_("FITS read partial: there is no layer %d in the image %s\n"),
				layer + 1, filename);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	if (fit->naxis == 3 && fit->naxes[2] != 3) {
		siril_log_message(_("Unsupported FITS image format (%ld axes).\n"),
				fit->naxes[2]);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	nbdata = area->w * area->h;
	fit->type = get_data_type(fit->bitpix);
	if (fit->type == DATA_UNSUPPORTED) {
		siril_log_message(_("Unknown FITS data format in internal conversion\n"));
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	if (fit->type == DATA_USHORT) {
		/* realloc fit->data to the image size */
		WORD *olddata = fit->data;
		if ((fit->data = realloc(fit->data, nbdata * sizeof(WORD))) == NULL) {
			PRINT_ALLOC_ERR;
			status = 0;
			fits_close_file(fit->fptr, &status);
			if (olddata)
				free(olddata);
			return -1;
		}
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data;
		fit->pdata[BLAYER] = fit->data;

		status = internal_read_partial_fits(fit->fptr, fit->naxes[1],
				fit->bitpix, fit->data, layer, area);
	} else {
		if (fit->bitpix == FLOAT_IMG || fit->bitpix == DOUBLE_IMG) {
			char str[FLEN_VALUE] = { 0 };
			status = 0;
			fits_read_key(fit->fptr, TSTRING, "PROGRAM", &str, NULL, &status);
			gboolean not_from_siril = status || g_ascii_strncasecmp(str, PACKAGE, strlen(PACKAGE));

			status = 0;
			fits_read_key(fit->fptr, TDOUBLE, "DATAMAX", &data_max, NULL, &status);
			if ((fit->bitpix == FLOAT_IMG && not_from_siril) || fit->bitpix == DOUBLE_IMG) {
				// override data_max if needed. In some images there are differences between max and data_max
				float mini, maxi;
				fit_stats(fit->fptr, &mini, &maxi);
				fit->keywords.data_max = (double) maxi;
				fit->keywords.data_min = (double) mini;
			}
		}

		/* realloc fit->fdata to the image size */
		float *olddata = fit->fdata;
		if ((fit->fdata = realloc(fit->fdata, nbdata * sizeof(float))) == NULL) {
			PRINT_ALLOC_ERR;
			status = 0;
			fits_close_file(fit->fptr, &status);
			if (olddata)
				free(olddata);
			return -1;
		}
		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = fit->fdata;
		fit->fpdata[BLAYER] = fit->fdata;

		status = internal_read_partial_fits(fit->fptr, fit->naxes[1],
				fit->bitpix, fit->fdata, layer, area);
	}

	if (status) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	read_fits_metadata(fit);
	fit->orig_ry = fit->naxes[1];
	fit->x_offset = area->x;
	fit->y_offset = fit->orig_ry - area->y - area->h;
	fit->naxes[0] = area->w;
	fit->naxes[1] = area->h;
	fit->rx = fit->naxes[0];
	fit->ry = fit->naxes[1];
	fit->naxes[2] = 1;
	fit->naxis = 2;


	/* Note: reading one channel of a multi-channel FITS poses a color management challenge:
	 * the original FITS would have had a 3-channel ICC profile (eg linear RGB) which would
	 * expect 3-channel data at the transform endpoint. Now we only have 1 channel of the 3.
	 * We set the icc_profile to NULL and fit->color_mananged to FALSE to indicate that this
	 * is raw data and no longer associated with a color managed image. Care must be taken
	 * regarding subsequent reintegration of this data into a color managed workflow.
	 */
	fit->icc_profile = NULL;
	color_manage(fit, FALSE);

	status = 0;
	fits_close_file(fit->fptr, &status);
	siril_debug_print("Loaded partial FITS file %s\n", filename);
	return 0;
}

/* Read a rectangular section of a FITS image in Siril's format, pointed by its
 * exact filename. All layers are read
 * Returned fit->data is upside-down. */
int readfits_partial_all_layers(const char *filename, fits *fit, const rectangle *area) {
	int status;
	size_t nbdata = 0;
	int nblayer = 1;

	status = 0;
	if (siril_fits_open_diskfile_img(&(fit->fptr), filename, READONLY, &status)) {
		report_fits_error(status);
		return status;
	}

	status = 0;
	fits_get_img_param(fit->fptr, 3, &(fit->bitpix), &(fit->naxis), fit->naxes,	&status);
	if (status) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return status;
	}

	if (fit->naxis == 2 && fit->naxes[2] == 0)
		fit->naxes[2] = 1;	// see readfits for the explanation

	if (fit->naxis == 3 && fit->naxes[2] != 3) {
		siril_log_message(_("Unsupported FITS image format (%ld axes).\n"),
				fit->naxes[2]);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	nbdata = area->w * area->h;
	nblayer = fit->naxes[2];
	fit->type = get_data_type(fit->bitpix);
	if (fit->type == DATA_UNSUPPORTED) {
		siril_log_message(_("Unknown FITS data format in internal conversion\n"));
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	if (fit->type == DATA_USHORT) {
		/* realloc fit->data to the image size */
		WORD *olddata = fit->data;
		if ((fit->data = realloc(fit->data, nblayer * nbdata * sizeof(WORD))) == NULL) {
			PRINT_ALLOC_ERR;
			status = 0;
			fits_close_file(fit->fptr, &status);
			if (olddata)
				free(olddata);
			return -1;
		}
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + nbdata;
		fit->pdata[BLAYER] = fit->data + nbdata * 2;

		for (int layer = 0; layer < nblayer; layer++) {
			status = internal_read_partial_fits(fit->fptr, fit->naxes[1],
					fit->bitpix, fit->pdata[layer], layer, area);
			if (status) {
				report_fits_error(status);
				status = 0;
				fits_close_file(fit->fptr, &status);
				if (olddata)
					free(olddata);
				return status;
			}
		}
	} else {
		/* realloc fit->fdata to the image size */
		float *olddata = fit->fdata;
		if ((fit->fdata = realloc(fit->fdata, nblayer * nbdata * sizeof(float))) == NULL) {
			PRINT_ALLOC_ERR;
			status = 0;
			fits_close_file(fit->fptr, &status);
			if (olddata)
				free(olddata);
			return -1;
		}
		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = fit->fdata + nbdata;
		fit->fpdata[BLAYER] = fit->fdata + nbdata * 2;

		for (int layer = 0; layer < nblayer; layer++) {
			status = internal_read_partial_fits(fit->fptr, fit->naxes[1],
				fit->bitpix, fit->fpdata[layer], layer, area);
			if (status) {
				report_fits_error(status);
				status = 0;
				fits_close_file(fit->fptr, &status);
				if (olddata)
					free(olddata);
				return status;
			}
		}
	}

	read_fits_metadata(fit);
	fit->naxes[0] = area->w;
	fit->naxes[1] = area->h;
	fit->rx = fit->naxes[0];
	fit->ry = fit->naxes[1];
	fit->naxes[2] = nblayer;
	fit->naxis = 3;

	status = 0;
	fits_close_file(fit->fptr, &status);
	siril_debug_print("Loaded partial FITS file %s\n", filename);
	return 0;
}

int read_fits_metadata(fits *fit) {
	int status = 0;
	fit->naxes[2] = 1;
	fits_get_img_param(fit->fptr, 3, &(fit->bitpix), &(fit->naxis), fit->naxes, &status);
	if (status) {
		siril_log_message(_("FITSIO error getting image parameters.\n"));
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return 1;
	}

	manage_bitpix(fit->fptr, &(fit->bitpix), &(fit->orig_bitpix));

	fit->rx = fit->naxes[0];
	fit->ry = fit->naxes[1];

	if (fit->naxis == 3 && fit->naxes[2] != 3) {
		siril_log_color_message(_("The FITS image contains more than 3 channels (%ld). Opening only the three first.\n"), "salmon", fit->naxes[2]);
		if (fit->naxis == 3) fit->naxes[2] = 3;
	}

	if (fit->naxis == 2 && fit->naxes[2] == 0) {
		fit->naxes[2] = 1;
		/* naxes[2] is set to 1 because:
		 * - it doesn't matter, since naxis is 2, it's not used
		 * - it's very convenient to use it in multiplications as the number of layers
		 */
	}
	if (fit->bitpix == LONGLONG_IMG) {
		siril_log_message(
				_("FITS images with 64 bits signed integer per pixel.channel are not supported.\n"));
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	read_icc_profile_from_fits(fit); // read color management data

	read_fits_header(fit);	// stores useful header data in fit
	fit->header = copy_header(fit);

	return 0;
}

static int read_fits_metadata_from_path_internal(const char *filename, fits *fit, gboolean image_file) {
	int status = 0;
	if (image_file)
		siril_fits_open_diskfile_img(&(fit->fptr), filename, READONLY, &status);
	else siril_fits_open_diskfile(&(fit->fptr), filename, READONLY, &status);
	if (status) {
		report_fits_error(status);
		return status;
	}

	read_fits_metadata(fit);
	status = 0;
	fits_close_file(fit->fptr, &status);

	return status;
}

// for an image
int read_fits_metadata_from_path(const char *filename, fits *fit) {
	return read_fits_metadata_from_path_internal(filename, fit, TRUE);
}

// for any type of FITS
int read_fits_metadata_from_path_first_HDU(const char *filename, fits *fit) {
	return read_fits_metadata_from_path_internal(filename, fit, FALSE);
}

gboolean flip_buffer(int bitpix, void *buffer, const rectangle *area) {
	/* reverse the read data, because it's stored upside-down */
	if (get_data_type(bitpix) == DATA_FLOAT) {
		int line_size = area->w * sizeof(float);
		void *swap = malloc(line_size);
		if (!swap) {
			return FALSE;
		}
		float *buf = (float *)buffer;
		int i;
		for (i = 0; i < area->h/2 ; i++) {
			memcpy(swap, buf + i*area->w, line_size);
			memcpy(buf + i*area->w, buf + (area->h - i - 1)*area->w, line_size);
			memcpy(buf + (area->h - i - 1)*area->w, swap, line_size);
		}
		free(swap);
	} else {
		int line_size = area->w * sizeof(WORD);
		void *swap = malloc(line_size);
		if (!swap) {
			return FALSE;
		}
		WORD *buf = (WORD *)buffer;
		int i;
		for (i = 0; i < area->h/2 ; i++) {
			memcpy(swap, buf + i*area->w, line_size);
			memcpy(buf + i*area->w, buf + (area->h - i - 1)*area->w, line_size);
			memcpy(buf + (area->h - i - 1)*area->w, swap, line_size);
		}
		free(swap);
	}
	return TRUE;
}


/* read subset of an opened fits file.
 * The rectangle's coordinates x,y start at 0,0 for first pixel in the image.
 * layer and index also start at 0.
 * buffer has to be allocated with enough space to store the area.
 */
int read_opened_fits_partial(sequence *seq, int layer, int index, void *buffer,
		const rectangle *area) {
	int status;

	if (!seq || !seq->fptr || !seq->fptr[index]) {
		printf("data initialization error in read fits partial\n");
		return 1;
	}
	int rx = (seq->is_variable) ? seq->imgparam[index].rx : seq->rx;
	int ry = (seq->is_variable) ? seq->imgparam[index].ry : seq->ry;
	if (area->x < 0 || area->y < 0 || area->x >= rx || area->y >= ry
			|| area->w <= 0 || area->h <= 0 || area->x + area->w > rx
			|| area->y + area->h > ry) {
		fprintf(stderr, "partial read from FITS file has been requested outside image bounds or with invalid size (%d: %d %d %d %d)\n", index + 1,
					area->x, area->y, area->w, area->h);
		return 1;
	}

#ifdef _OPENMP
	g_assert(seq->fd_lock);
	omp_set_lock(&seq->fd_lock[index]);
#endif

	status = internal_read_partial_fits(seq->fptr[index], ry, seq->bitpix, buffer, layer, area);

#ifdef _OPENMP
	omp_unset_lock(&seq->fd_lock[index]);
#endif
	if (status)
		return 1;

	if (!flip_buffer(seq->bitpix, buffer, area)) {
		return 1;
	}
	return 0;
}

int siril_fits_compress(fits *f) {
	int status = 0;
	int comp_type = -1;
	siril_debug_print("Compressing FIT file with method %d and quantization %f\n",
				com.pref.comp.fits_method,
				com.pref.comp.fits_quantization);
	comp_type = get_compression_type(com.pref.comp.fits_method);
	siril_debug_print("cfitsio compression type %d\n",
				comp_type);
	if (comp_type < 0) {
		siril_log_message(_("Unknown FITS compression method in internal conversion\n"));
		return 1;
	}

	/** In the doc: The default algorithm is called ``SUBTRACTIVE_DITHER_1''.
	 * A second variation called ``SUBTRACTIVE_DITHER_2'' is also available,
	 * which does the same thing except that any pixels with a value of 0.0 are not
	 * dithered and instead the zero values are exactly preserved in the compressed image
	 *
	 * Here we need to use SUBTRACTIVE_DITHER in order to not have border artifacts at the
	 * end of stacking with compressed images.
	 */
	if (fits_set_quantize_dither(f->fptr, SUBTRACTIVE_DITHER_2, &status)) {
		report_fits_error(status);
		return 1;
	}
	status = 0;

	if (fits_set_compression_type(f->fptr, comp_type, &status)) {
		report_fits_error(status);
		return 1;
	}
	status = 0;

	if (fits_set_quantize_level(f->fptr, com.pref.comp.fits_quantization, &status)) {
		report_fits_error(status);
		return 1;
	}

	status = 0;

	/* Set the Hcompress scale factor if relevant */
	if (comp_type == HCOMPRESS_1) {
		if (fits_set_hcomp_scale(f->fptr, com.pref.comp.fits_hcompress_scale, &status)) {
			report_fits_error(status);
			return 1;
		}
		siril_debug_print("FITS HCompress scale factor %f\n",
				com.pref.comp.fits_hcompress_scale);
		status = 0;
	}
	return status;
}

gchar *set_right_extension(const char *name) {
	gchar *filename = NULL;

	gboolean comp_flag = FALSE;
	/* first check if there is fz extension */
	if (g_str_has_suffix(name, ".fz")) {
		comp_flag = TRUE;
	}

	gboolean right_extension = FALSE;
	for (int i = 0; i < G_N_ELEMENTS(fit_extension); i++) {
		gchar *extension;
		if (comp_flag) {
			extension = g_strdup_printf("%s.fz", fit_extension[i]);
		} else {
			extension = g_strdup(fit_extension[i]);
		}
		if (g_str_has_suffix(name, extension)) {
			right_extension = TRUE;
			g_free(extension);
			break;
		}
		g_free(extension);
	}

	if (!right_extension) {
		if (com.pref.comp.fits_enabled) {
			filename = g_strdup_printf("%s%s.fz", name, com.pref.ext);
		} else {
			filename = g_strdup_printf("%s%s", name, com.pref.ext);
		}
	} else {
		if (comp_flag && !com.pref.comp.fits_enabled) {
			/* we remove .fz */
			gchar *tmp = g_strdup(name);
			tmp[strlen(tmp) - 3] = '\0';
			filename = g_strdup_printf("%s", tmp);

			g_free(tmp);
		} else if (!comp_flag && com.pref.comp.fits_enabled) {
			filename = g_strdup_printf("%s.fz", name);

		} else {
			filename = g_strdup_printf("%s", name);
		}
	}
	return filename;
}

/* creates, saves and closes the file associated to f, overwriting previous  */
int savefits(const char *name, fits *f) {
	int status;

	f->naxes[0] = f->rx;
	f->naxes[1] = f->ry;

	if (f->naxis == 3 && f->naxes[2] != 3) {
		printf("Trying to save a FITS color file with more than 3 channels?\n");
		return 1;
	}

	gchar *filename = set_right_extension(name);
	if (!filename) return 1;

	if (g_unlink(filename))
		siril_debug_print("g_unlink() failed\n"); /* Delete old file if it already exists */

	status = 0;
	if (siril_fits_create_diskfile(&(f->fptr), filename, &status)) { /* create new FITS file */
		report_fits_error(status);
		g_free(filename);
		return 1;
	}

	if (com.pref.comp.fits_enabled) {
		status = siril_fits_compress(f);
		if (status) {
			report_fits_error(status);
			g_free(filename);
			return 1;
		}
	}

	if (fits_create_img(f->fptr, f->bitpix, f->naxis, f->naxes, &status)) {
		report_fits_error(status);
		g_free(filename);
		return 1;
	}

	if (save_opened_fits(f)) {
		status = 0;
		fits_close_file(f->fptr, &status);
		f->fptr = NULL;
		g_free(filename);
		return 1;
	}

	if (com.pref.fits_save_icc && f->color_managed) {
		/* Only write the ICC profile for color managed FITS. This avoids writing
		 * ICC profiles to things like extracted channels where it doesn't really
		 * make sense. */
		if (f->icc_profile) {
			write_icc_profile_to_fits(f);
		} else {
			siril_debug_print("Info: FITS has no assigned ICC profile, saving without one.\n");
		}
	}

	if (f->checksum) {
		status = 0;
		fits_write_chksum(f->fptr, &status);
		if (status) {
			fits_report_error(stderr, status);
			return 1;
		}
	}

	status = 0;
	fits_close_file(f->fptr, &status);
	if (!status) {
		siril_log_message(_("Saving FITS: file %s, %ld layer(s), %ux%u pixels, %d bits\n"),
				filename, f->naxes[2], f->rx, f->ry,
				f->bitpix == FLOAT_IMG ? 32 : (f->bitpix == SHORT_IMG || f->bitpix == USHORT_IMG ? 16 : 8));
	}
	g_free(filename);
	return 0;
}

int save_opened_fits(fits *f) {
	BYTE *data8;
	long orig[3] = { 1L, 1L, 1L };
	size_t i, pixel_count;
	int status = 0;
	signed short *data;

	save_fits_header(f);

	pixel_count = f->naxes[0] * f->naxes[1] * f->naxes[2];

	status = 0;
	switch (f->bitpix) {
	case BYTE_IMG:
		data8 = malloc(pixel_count * sizeof(BYTE));
		if (f->type == DATA_FLOAT) {
			for (i = 0; i < pixel_count; i++) {
				data8[i] = float_to_uchar_range(f->fdata[i]);
			}
		} else {
			double norm = get_normalized_value(f);
			for (i = 0; i < pixel_count; i++) {
				if (norm == USHRT_MAX_DOUBLE)
					data8[i] = conv_to_BYTE((double)f->data[i]);
				else
					data8[i] = truncate_to_BYTE(f->data[i]);
			}
		}
		if (fits_write_pix(f->fptr, TBYTE, orig, pixel_count, data8, &status)) {
			report_fits_error(status);
			free(data8);
			return 1;
		}
		f->keywords.lo >>= 8;
		f->keywords.hi >>= 8;
		g_free(data8);
		break;
	case SHORT_IMG:
		if (f->type == DATA_FLOAT) {
			data = float_buffer_to_short(f->fdata, f->naxes[0] * f->naxes[1] * f->naxes[2]);
		} else {
			if (f->orig_bitpix == BYTE_IMG) {
				conv_8_to_16(f->data, pixel_count);
			}
			data = ushort_buffer_to_short(f->data, f->naxes[0] * f->naxes[1] * f->naxes[2]);
		}
		if (fits_write_pix(f->fptr, TSHORT, orig, pixel_count, data, &status)) {
			report_fits_error(status);
			free(data);
			return 1;
		}
		free(data);
		break;
	case USHORT_IMG:
		if (f->type == DATA_FLOAT) {
			WORD *datau = float_buffer_to_ushort(f->fdata, f->naxes[0] * f->naxes[1] * f->naxes[2]);
			if (fits_write_pix(f->fptr, TUSHORT, orig, pixel_count, datau, &status)) {
				report_fits_error(status);
				g_free(datau);
				return 1;
			}
			free(datau);
		} else {
			if (f->orig_bitpix == BYTE_IMG) {
				conv_8_to_16(f->data, pixel_count);
			}
			if (fits_write_pix(f->fptr, TUSHORT, orig, pixel_count, f->data, &status)) {
				report_fits_error(status);
				return 1;
			}
		}
		break;
	case FLOAT_IMG:
		if (f->type == DATA_USHORT) {
			if (f->orig_bitpix == BYTE_IMG) {
				conv_8_to_16(f->data, pixel_count);
			}
			f->fdata = malloc(pixel_count * sizeof(float));
			conv_16_to_32(f->data, f->fdata, pixel_count);
			fit_replace_buffer(f, f->fdata, DATA_FLOAT);
		}
		if (fits_write_pix(f->fptr, TFLOAT, orig, pixel_count, f->fdata, &status)) {
			report_fits_error(status);
			return 1;
		}
		break;
	case LONG_IMG:
	case LONGLONG_IMG:
	case DOUBLE_IMG:
	default:
		siril_log_message(_("ERROR: trying to save a FITS image "
				"with an unsupported format (%d).\n"), f->bitpix);
		fits_close_file(f->fptr, &status);
		return 1;
	}

	if (!status) {
		// copy the entire header in memory
		if (f->header)
			free(f->header);
		f->header = copy_header(f);
	}

	return 0;
}

/* Duplicates some of a fits data into another, with various options; the third
 * parameter, oper, indicates with bits what operations will be done:
 *
 * - CP_ALLOC: allocates the to->data pointer to the size of from->data and
 *   sets to->pdata; required if data is not already allocated with the
 *   correct size or at all. No data is copied
 * - CP_INIT: initialize to->data with zeros, same size of the image in from,
 *   but no other data is modified. Ignored if not used with CP_ALLOC.
 * - CP_COPYA: copies the actual data, from->data to to->data on all layers,
 *   but no other information from the source. Should not be used with CP_INIT
 * - CP_FORMAT: copy all metadata and leaves data to null
 * - CP_COPYMASK: copy the image mask
 * - CP_EXPAND: forces the destination number of layers to be taken as 3, but
 *   the other operations have no modifications, meaning that if the source
 *   image has one layer, the output image will have only one actual layer
 *   filled, and two filled with random data unless CP_INIT is used to fill it
 *   with zeros.
 *
 * Example: to duplicate a fits from one to an unknown-allocated other, those
 * flags should be used:	CP_ALLOC | CP_COPYA | CP_FORMAT
 *
 */
int copyfits(fits *from, fits *to, unsigned char oper, int layer) {
	int depth, i;
	size_t nbdata = from->naxes[0] * from->naxes[1];

	if ((oper & CP_EXPAND))
		depth = 3;
	else depth = from->naxes[2];

	if ((oper & CP_FORMAT)) {
		// free anything that might need deallocating in to
		clearfits(to);
		// copying metadata, not data or stats which are kept null
		memcpy(to, from, sizeof(fits));
		to->naxis = depth == 3 ? 3 : from->naxis;
		to->naxes[2] = depth;
		if (depth != from->naxes[2]) {
			to->maxi = -1.0;
		}
		to->stats = NULL;
		to->fptr = NULL;
		to->data = NULL;
		to->pdata[0] = NULL;
		to->pdata[1] = NULL;
		to->pdata[2] = NULL;
		to->fdata = NULL;
		to->fpdata[0] = NULL;
		to->fpdata[1] = NULL;
		to->fpdata[2] = NULL;
		to->mask = NULL;
		to->header = NULL;
		to->unknown_keys = NULL;
		to->history = NULL;
		to->keywords.date = NULL;
		to->keywords.date_obs = NULL;
		to->icc_profile = NULL;
		to->color_managed = FALSE;
		to->keywords.wcslib = NULL;
	}

	if ((oper & CP_ALLOC)) {
		// allocating to->data and assigning to->pdata
		if (from->type == DATA_USHORT) {
			WORD *olddata = to->data;
			if (!(to->data = realloc(to->data, nbdata * depth * sizeof(WORD)))) {
				PRINT_ALLOC_ERR;
				if (olddata)
					free(olddata);
				return -1;
			}
			to->type = DATA_USHORT;
			to->pdata[RLAYER] = to->data;
			if (depth == 3) {
				to->pdata[GLAYER] = to->data + nbdata;
				to->pdata[BLAYER] = to->data + 2 * nbdata;
			} else {
				to->pdata[GLAYER] = to->data;
				to->pdata[BLAYER] = to->data;
			}

			if ((oper & CP_INIT)) {
				// clearing to->data allocated above
				memset(to->data, 0, nbdata * depth * sizeof(WORD));
			}
		}
		else if (from->type == DATA_FLOAT) {
			float *olddata = to->fdata;
			if (!(to->fdata = realloc(to->fdata, nbdata * depth * sizeof(float)))) {
				PRINT_ALLOC_ERR;
				if (olddata)
					free(olddata);
				return -1;
			}
			to->type = DATA_FLOAT;
			to->fpdata[RLAYER] = to->fdata;
			if (depth == 3) {
				to->fpdata[GLAYER] = to->fdata + nbdata;
				to->fpdata[BLAYER] = to->fdata + 2 * nbdata;
			} else {
				to->fpdata[GLAYER] = to->fdata;
				to->fpdata[BLAYER] = to->fdata;
			}

			if ((oper & CP_INIT)) {
				// clearing to->fdata allocated above
				memset(to->fdata, 0, nbdata * depth * sizeof(float));
			}
		}
		else {
			fprintf(stderr, "unsupported copy\n");
			return -1;
		}
	}

	if ((oper & CP_COPYA)) {
		// copying data
		if (to->type == DATA_USHORT) {
			if (!(to->data)) {
				fprintf(stderr, "error: data ptr unallocated\n");
				return -1;
			}
			if (!from->data) {
				fprintf(stderr, "error: no suitable data in src fits\n");
				return -1;
			}
			memcpy(to->data, from->data, nbdata * depth * sizeof(WORD));
		} else if (to->type == DATA_FLOAT) {
			if (!(to->fdata)) {
				fprintf(stderr, "error: fdata ptr unallocated\n");
				return -1;
			}
			if (!from->fdata) {
				fprintf(stderr, "error: no suitable data in src fits\n");
				return -1;
			}
			memcpy(to->fdata, from->fdata, nbdata * depth * sizeof(float));
		} else {
			fprintf(stderr, "unsupported copy\n");
			return -1;
		}

		// copying stats
		if (from->stats) {
			for (i = 0; i < from->naxes[2]; i++) {
				if (from->stats[i])
					add_stats_to_fit(to, i, from->stats[i]);
			}
		} else {
			invalidate_stats_from_fit(to);
		}
	}

	if (oper & CP_COPYMASK) {
		if (from->mask) {
			to->mask = malloc(nbdata * sizeof(uint8_t));
			memcpy(to->mask, from->mask, nbdata * sizeof(uint8_t));
		} else {
			to->mask = NULL;
		}
	}

	if ((oper & CP_ALLOC) || (oper & CP_COPYA)) {
		// copy color management data
		to->color_managed = from->color_managed;
		if (to->color_managed) {
			to->icc_profile = copyICCProfile(from->icc_profile);
		} else {
			to->icc_profile = NULL;
		}
		color_manage(to, to->color_managed);
	}

	return 0;
}

/* Changes a FITS to [layers] channels. Data in the GLAYER and BLAYER
 * is not zeroed so will initially be filled with junk, it is assumed that the
 * calling function will subsequently populate the data in these channels.
 */

int fits_change_depth(fits *fit, int layers) {
	if (layers == fit->naxes[2]) return 0;
	g_assert (layers == 1 || layers == 3); // Can only change depth to mono or 3-channel

	size_t nbdata = fit->naxes[0] * fit->naxes[1];
	if (fit->type == DATA_USHORT) {
		WORD *tmp;
		if (!(tmp = realloc(fit->data, nbdata * layers * sizeof(WORD)))) {
			PRINT_ALLOC_ERR;
			return -1;
		}
		fit->data = tmp;
		fit->naxes[2] = layers;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = layers == 3 ? fit->data + nbdata : fit->data;
		fit->pdata[BLAYER] = layers == 3 ? fit->data + 2 * nbdata : fit->data;
	}
	else if (fit->type == DATA_FLOAT) {
		float *tmp;
		if (!(tmp = realloc(fit->fdata, nbdata * layers * sizeof(float)))) {
			PRINT_ALLOC_ERR;
			return -1;
		}
		fit->fdata = tmp;
		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = layers == 3 ? fit->fdata + nbdata : fit->fdata;
		fit->fpdata[BLAYER] = layers == 3 ? fit->fdata + 2 * nbdata: fit->fdata;
	}
	fit->naxes[2] = layers;
	fit->naxis = layers == 1 ? 2 : 3;
	return 0;
}

int extract_fits(fits *from, fits *to, int channel, gboolean to_float) {
	size_t nbdata = from->naxes[0] * from->naxes[1];
	// copying metadata, not data or stats which are kept null
	memcpy(to, from, sizeof(fits));
	to->naxis = 2;
	to->naxes[2] = 1L;
	to->maxi = -1.0;
	to->stats = NULL;
	to->fptr = NULL;
	to->header = NULL;
	to->unknown_keys = NULL;
	to->history = NULL;
	to->keywords.date = NULL;
	to->keywords.date_obs = NULL;
	/* Note: extracting one channel of a multi-channel FITS poses a color management challenge:
	 * the original FITS would have had a 3-channel ICC profile (eg linear RGB) which would
	 * expect 3-channel data at the transform endpoint. Now we only have 1 channel of the 3.
	 * We therefore set a NULL icc_profile and color_managed = FALSE to indicate that this
	 * extracted channel is considered non-color managed raw data. Care must be taken if
	 * subsequently reintegrating it into a color managed workflow.
	 */
	color_manage(to, FALSE);
	to->icc_profile = NULL;
	to->keywords.wcslib = NULL;

	if (from->type == DATA_USHORT)
		if (to_float) {
			/*if (from->bitpix == BYTE_IMG)
				to->fdata = ushort8_buffer_to_float(from->pdata[channel], nbdata);
			else to->fdata = ushort_buffer_to_float(from->pdata[channel], nbdata); */
			to->fdata = malloc(nbdata * sizeof(float));
			if (!to->fdata) {
				PRINT_ALLOC_ERR;
				return -1;
			}
			for (int i = 0; i < nbdata; i++)
				to->fdata[i] = (float)from->pdata[channel][i];
			to->type = DATA_FLOAT;
			to->bitpix = FLOAT_IMG;
			to->data = NULL;
			to->pdata[0] = NULL;
			to->pdata[1] = NULL;
			to->pdata[2] = NULL;
			to->fpdata[0] = to->fdata;
			to->fpdata[1] = to->fdata;
			to->fpdata[2] = to->fdata;
		}
		else {
			to->data = malloc(nbdata * sizeof(WORD));
			if (!to->data) {
				PRINT_ALLOC_ERR;
				return -1;
			}
			memcpy(to->data, from->pdata[channel], nbdata * sizeof(WORD));
			to->pdata[0] = to->data;
			to->pdata[1] = to->data;
			to->pdata[2] = to->data;
		}
	else if (from->type == DATA_FLOAT) {
		to->fdata = malloc(nbdata * sizeof(float));
		if (!to->fdata) {
			PRINT_ALLOC_ERR;
			return -1;
		}
		memcpy(to->fdata, from->fpdata[channel], nbdata * sizeof(float));
		to->fpdata[0] = to->fdata;
		to->fpdata[1] = to->fdata;
		to->fpdata[2] = to->fdata;
	}
	return 0;
}

void keep_only_first_channel(fits *fit) {
	fit->naxis = 2;
	fit->naxes[2] = 1L;
	fit->maxi = -1.0;
	fit->history = NULL;
	invalidate_stats_from_fit(fit);
	size_t nbdata = fit->naxes[0] * fit->naxes[1];
	if (fit->type == DATA_USHORT) {
		fit->data = realloc(fit->data, nbdata * sizeof(WORD));
		fit->pdata[1] = fit->data;
		fit->pdata[2] = fit->data;
	} else if (fit->type == DATA_FLOAT) {
		fit->fdata = realloc(fit->fdata, nbdata * sizeof(float));
		fit->fpdata[1] = fit->fdata;
		fit->fpdata[2] = fit->fdata;
	}
	/* Note: keeping one channel of a multi-channel FITS poses a color management challenge:
	 * the original FITS would have had a 3-channel ICC profile (eg linear RGB) which would
	 * expect 3-channel data at the transform endpoint. Now we only have 1 channel of the 3.
	 * As above, we set icc_profile to NULL and color_managed to FALSE.
	 */
	fit->icc_profile = NULL;
	color_manage(fit, FALSE);
}

void copy_fits_metadata(fits *from, fits *to) {
	// Copy simple fields
	memcpy(&to->keywords, &from->keywords, sizeof(fkeywords));

	// Copy other structures
	memcpy(&to->keywords.dft, &from->keywords.dft, sizeof(dft_info));
	memcpy(&to->keywords.wcsdata, &from->keywords.wcsdata, sizeof(wcs_info));

	to->keywords.date = NULL; // will be set at save
	// Copy date_obs
	to->keywords.date_obs = NULL;
	if (from->keywords.date_obs) {
		to->keywords.date_obs = g_date_time_ref(from->keywords.date_obs);
	}

	if (from->keywords.wcslib) {
		int status = -1;
		to->keywords.wcslib = wcs_deepcopy(from->keywords.wcslib, &status);
		if (status) {
			wcsfree(to->keywords.wcslib);
			siril_debug_print("could not copy wcslib struct\n");
		}
	}

	if (from->unknown_keys) {
		to->unknown_keys = g_strdup(from->unknown_keys);
	}
	// Set boolean flags
//	to->pixelkey = (from->keywords.pixel_size_x > 0.);
//	to->focalkey = (from->keywords.focal_length > 0.);
//	/* override these two keys. */
	to->pixelkey = from->pixelkey;
	to->focalkey = from->focalkey;

	// copy from->history?

}


int copy_fits_from_file(const char *source, const char *destination) {
	fitsfile *infptr, *outfptr; /* FITS file pointers defined in fitsio.h */
	int status = 0; /* status must always be initialized = 0  */

	/* Open the input file */
	if (!siril_fits_open_diskfile_img(&infptr, source, READONLY, &status)) {
		/* Create the output file */
		if (!siril_fits_create_diskfile(&outfptr, destination, &status)) {

			/* copy the previous, current, and following HDUs */
			fits_copy_file(infptr, outfptr, 1, 1, 1, &status);

			fits_close_file(outfptr, &status);
		}
		fits_close_file(infptr, &status);
	}

	/* if error occured, print out error message */
	if (status)
		report_fits_error(status);
	return (status);
}


int save1fits16(const char *filename, fits *fit, int layer) {
	if (layer != RLAYER) {
		size_t nbdata = fit->naxes[0] * fit->naxes[1];
		memcpy(fit->data, fit->data + layer * nbdata, nbdata * sizeof(WORD));
	}
	fit->naxis = 2;
	fit->naxes[2] = 1;
	if (fit->icc_profile) {
		cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = NULL;
	}
	color_manage(fit, FALSE);
	int retval = savefits(filename, fit);
	return retval;
}

int save1fits32(const char *filename, fits *fit, int layer) {
	if (layer != RLAYER) {
		size_t nbdata = fit->naxes[0] * fit->naxes[1];
		memcpy(fit->fdata, fit->fdata + layer * nbdata, nbdata * sizeof(float));
	}
	fit->naxis = 2;
	fit->naxes[2] = 1;
	if (fit->icc_profile) {
		cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = NULL;
	}
	color_manage(fit, FALSE);
	int retval = savefits(filename, fit);
	return retval;
}

/* this method converts 24-bit RGB or BGR data (no padding) to 48-bit FITS data.
 * order is RGB when inverted is FALSE, BGR when inverted is TRUE
 * fit->data has to be already allocated and fit->rx and fit->ry must be correct */
void rgb24bit_to_fits48bit(unsigned char *rgbbuf, fits *fit, gboolean inverted) {
	size_t nbdata = fit->naxes[0] * fit->naxes[1];
	WORD *rdata, *gdata, *bdata;
	rdata = fit->pdata[RLAYER] = fit->data;
	gdata = fit->pdata[GLAYER] = fit->data + nbdata;
	bdata = fit->pdata[BLAYER] = fit->data + 2 * nbdata;
	for (size_t i = 0; i < nbdata; ++i) {
		if (inverted)
			*bdata++ = (WORD) *rgbbuf++;
		else
			*rdata++ = (WORD) *rgbbuf++;
		*gdata++ = (WORD) *rgbbuf++;
		if (inverted)
			*rdata++ = (WORD) *rgbbuf++;
		else
			*bdata++ = (WORD) *rgbbuf++;
	}
}

/* this method converts 8-bit gray data to 16-bit FITS data.
 * fit->data has to be already allocated and fit->rx and fit->ry must be correct */
void rgb8bit_to_fits16bit(const unsigned char *graybuf, fits *fit) {
	WORD *data;
	size_t i, nbdata = fit->naxes[0] * fit->naxes[1];
	fit->pdata[0] = fit->data;
	fit->pdata[1] = fit->data;
	fit->pdata[2] = fit->data;
	data = fit->data;
	for (i = 0; i < nbdata; ++i) {
		*data++ = *graybuf++;
	}
}

/* this method converts 48-bit RGB or BGR data (no padding) to 48-bit FITS data.
 * order is RGB when inverted is FALSE, BGR when inverted is TRUE
 * the endianness of the data, since we have two byte per value, may not match the endianness
 * of our FITS files, so the change_endian parameter allows to flip the endian.
 * fit->data has to be already allocated and fit->rx and fit->ry must be correct */
void rgb48bit_to_fits48bit(const WORD *rgbbuf, fits *fit, gboolean inverted,
		gboolean change_endian) {
	size_t i, nbdata = fit->naxes[0] * fit->naxes[1];
	WORD *rdata, *gdata, *bdata, curval;
	rdata = fit->pdata[0] = fit->data;
	gdata = fit->pdata[1] = fit->data + nbdata;
	bdata = fit->pdata[2] = fit->data + 2 * nbdata;
	for (i = 0; i < nbdata; ++i) {
		curval = *rgbbuf++;
		if (change_endian)
			curval = change_endianness16(curval);
		if (inverted)
			*bdata++ = curval;
		else
			*rdata++ = curval;

		curval = *rgbbuf++;
		if (change_endian)
			curval = change_endianness16(curval);
		*gdata++ = curval;

		curval = *rgbbuf++;
		if (change_endian)
			curval = change_endianness16(curval);
		if (inverted)
			*rdata++ = curval;
		else
			*bdata++ = curval;
	}
}

/* this method flips top-bottom of fit data.
 * fit->rx, fit->ry, fit->naxes[2] and fit->pdata[*] are required to be assigned correctly */
static void fits_flip_top_to_bottom_ushort(fits *fit) {
	int line, axis, line_size;
	WORD *swapline, *src, *dst;

	line_size = fit->rx * sizeof(WORD);
	swapline = malloc(line_size);

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
}

static void fits_flip_top_to_bottom_float(fits *fit) {
	int line, axis, line_size;
	float *swapline, *src, *dst;

	line_size = fit->rx * sizeof(float);
	swapline = malloc(line_size);

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->fpdata[axis] + line * fit->rx;
			dst = fit->fpdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
}

void fits_flip_top_to_bottom(fits *fit) {
	if (fit->type == DATA_USHORT)
		fits_flip_top_to_bottom_ushort(fit);
	else if (fit->type == DATA_FLOAT)
		fits_flip_top_to_bottom_float(fit);
}

/* This function copies an area from the fits 'from' on layer 'layer' into
 * another and initializes all relevant data */
/* the crop function does the same but in place and for all channels without
 * reallocating */
static void extract_region_from_fits_ushort(fits *from, int layer, fits *to,
		const rectangle *area) {
	int x, y, d, ystart, yend;
	clearfits(to);
	to->data = malloc(area->w * area->h * sizeof(WORD));

	d = 0;
	ystart = from->ry - area->y - area->h;
	yend = from->ry - area->y;
	for (y = ystart; y < yend; y++) {
		for (x = area->x; x < area->x + area->w; x++) {
			to->data[d++] = from->pdata[layer][x + y * from->rx];
		}
	}

	to->rx = area->w;
	to->ry = area->h;
	to->naxes[0] = area->w;
	to->naxes[1] = area->h;
	to->naxes[2] = 1;
	to->naxis = 2;
	to->pdata[0] = to->data;
	to->pdata[1] = to->data;
	to->pdata[2] = to->data;
	to->bitpix = from->bitpix;
	to->type = DATA_USHORT;
	to->icc_profile = NULL;
	color_manage(to, FALSE);
}

static void extract_region_from_fits_float(fits *from, int layer, fits *to,
		const rectangle *area) {
	int x, y, d, ystart, yend;
	clearfits(to);
	to->fdata = malloc(area->w * area->h * sizeof(float));

	d = 0;
	ystart = from->ry - area->y - area->h;
	yend = from->ry - area->y;
	for (y = ystart; y < yend; y++) {
		for (x = area->x; x < area->x + area->w; x++) {
			to->fdata[d++] = from->fpdata[layer][x + y * from->rx];
		}
	}

	to->rx = area->w;
	to->ry = area->h;
	to->naxes[0] = area->w;
	to->naxes[1] = area->h;
	to->naxes[2] = 1;
	to->naxis = 2;
	to->fpdata[0] = to->fdata;
	to->fpdata[1] = to->fdata;
	to->fpdata[2] = to->fdata;
	to->bitpix = from->bitpix;
	to->type = DATA_FLOAT;
	to->icc_profile = NULL;
	color_manage(to, FALSE);
}

void extract_region_from_fits(fits *from, int layer, fits *to,
		const rectangle *area) {
	if (from->type == DATA_USHORT)
		extract_region_from_fits_ushort(from, layer, to, area);
	else if (from->type == DATA_FLOAT)
		extract_region_from_fits_float(from, layer, to, area);
	to->orig_ry = from->ry;
	to->x_offset = area->x;
	to->y_offset = to->orig_ry - area->y - area->h;
}

int new_fit_image(fits **fit, int width, int height, int nblayer, data_type type) {
	return new_fit_image_with_data(fit, width, height, nblayer, type, NULL);
}

/* creates a new fit image from scratch (NULL fit) or into a fits * previously
 * allocated, with provided or non-cleared (random) data. */
int new_fit_image_with_data(fits **fit, int width, int height, int nblayer, data_type type, void *data) {
	size_t npixels, data_size;
	gboolean data_is_local = FALSE;
	g_assert(width > 0);
	g_assert(height > 0);
	g_assert(nblayer == 1 || nblayer == 3);

	npixels = width * height;
	data_size = type == DATA_USHORT ? sizeof(WORD) : sizeof(float);

	if (!data) {
		data = calloc(npixels * nblayer, data_size);
		if (!data) {
			PRINT_ALLOC_ERR;
			return -1;
		}
		data_is_local = TRUE;
	}

	if (*fit)
		clearfits(*fit);
	else {
		*fit = calloc(1, sizeof(fits));
		if (!*fit) {
			PRINT_ALLOC_ERR;
			if (data_is_local)
				free(data);
			return -1;
		}
	}

	if (nblayer > 1)
		(*fit)->naxis = 3;
	else (*fit)->naxis = 2;
	(*fit)->rx = width;
	(*fit)->ry = height;
	(*fit)->naxes[0] = width;
	(*fit)->naxes[1] = height;
	(*fit)->naxes[2] = nblayer;
	(*fit)->type = type;

	if (type == DATA_USHORT) {
		(*fit)->bitpix = USHORT_IMG;
		(*fit)->orig_bitpix = USHORT_IMG;
		(*fit)->data = (WORD *)data;
		(*fit)->pdata[RLAYER] = (*fit)->data;
		if ((*fit)->naxis == 3) {
			(*fit)->pdata[GLAYER] = (*fit)->data + npixels;
			(*fit)->pdata[BLAYER] = (*fit)->data + npixels* 2;
		} else {
			(*fit)->pdata[GLAYER] = (*fit)->data;
			(*fit)->pdata[BLAYER] = (*fit)->data;
		}
	} else if (type == DATA_FLOAT) {
		(*fit)->bitpix = FLOAT_IMG;
		(*fit)->orig_bitpix = FLOAT_IMG;
		(*fit)->fdata = (float *)data;
		(*fit)->fpdata[RLAYER] = (*fit)->fdata;
		if ((*fit)->naxis == 3) {
			(*fit)->fpdata[GLAYER] = (*fit)->fdata + npixels;
			(*fit)->fpdata[BLAYER] = (*fit)->fdata + npixels * 2;
		} else {
			(*fit)->fpdata[GLAYER] = (*fit)->fdata;
			(*fit)->fpdata[BLAYER] = (*fit)->fdata;
		}
	} else { // DATA_UNSUPPORTED
		free(data);
	}
	// Note: FITS created in this way will initially not be color managed. The
	// user, or the calling function, must assign an ICC profile if required.
	return 0;
}

/* use for type change, it doesn't free old data if it's the same type */
void fit_replace_buffer(fits *fit, void *newbuf, data_type newtype) {
	fit->type = newtype;
	invalidate_stats_from_fit(fit);
	size_t nbdata = fit->naxes[0] * fit->naxes[1];
	if (newtype == DATA_USHORT) {
		fit->bitpix = USHORT_IMG;
		fit->orig_bitpix = USHORT_IMG;
		fit->data = (WORD *)newbuf;
		fit->pdata[RLAYER] = fit->data;
		if (fit->naxis == 3) {
			fit->pdata[GLAYER] = fit->data + nbdata;
			fit->pdata[BLAYER] = fit->data + nbdata * 2;
		} else {
			fit->pdata[GLAYER] = fit->data;
			fit->pdata[BLAYER] = fit->data;
		}
		if (fit->fdata) {
			free(fit->fdata);
			fit->fdata = NULL;
		}
		fit->fpdata[0] = NULL;
		fit->fpdata[1] = NULL;
		fit->fpdata[2] = NULL;

		siril_debug_print("Changed a fit data (WORD)\n");
	} else if (newtype == DATA_FLOAT) {
		fit->bitpix = FLOAT_IMG;
		fit->orig_bitpix = FLOAT_IMG;
		fit->fdata = (float *)newbuf;
		fit->fpdata[RLAYER] = fit->fdata;
		if (fit->naxis == 3) {
			fit->fpdata[GLAYER] = fit->fdata + nbdata;
			fit->fpdata[BLAYER] = fit->fdata + nbdata * 2;
		} else {
			fit->fpdata[GLAYER] = fit->fdata;
			fit->fpdata[BLAYER] = fit->fdata;
		}
		if (fit->data) {
			free(fit->data);
			fit->data = NULL;
		}
		fit->pdata[0] = NULL;
		fit->pdata[1] = NULL;
		fit->pdata[2] = NULL;
		siril_debug_print("Changed a fit data (FLOAT)\n");
	}
}

void fit_debayer_buffer(fits *fit, void *newbuf) {
	size_t nbdata = fit->naxes[0] * fit->naxes[1];

	/* before changing naxis, we clear fit->stats that was allocated for
	 * one channel */
	full_stats_invalidation_from_fit(fit);

	fit->naxis = 3;
	fit->naxes[2] = 3;
	if (fit->type == DATA_USHORT) {
		if (fit->data)
			free(fit->data);
		fit->data = (WORD *)newbuf;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + nbdata;
		fit->pdata[BLAYER] = fit->data + nbdata * 2;
	}
	else if (fit->type == DATA_FLOAT) {
		if (fit->fdata)
			free(fit->fdata);
		fit->fdata = (float *)newbuf;
		fit->fpdata[RLAYER] = fit->fdata;
		fit->fpdata[GLAYER] = fit->fdata + nbdata;
		fit->fpdata[BLAYER] = fit->fdata + nbdata * 2;
	}
}

static inline void gray2rgb(float gray, guchar *rgb) {
	guchar val = (guchar) roundf_to_BYTE(255.f * gray);
	*rgb++ = val;
	*rgb++ = val;
	*rgb++ = val;
}

static inline void set_rgb(float r, float g, float b, guchar *rgb) {
	*rgb++ = (guchar) roundf_to_BYTE(255.f * r);
	*rgb++ = (guchar) roundf_to_BYTE(255.f * g);
	*rgb++ = (guchar) roundf_to_BYTE(255.f * b);
}

// Choose thread count based on workload
static inline int choose_num_threads(int W, int H, int max_threads) {
	int pixels = W * H;
	if (pixels < 65536) { // < 64k pixels → single-thread
		return 1;
	}
	int threads = pixels / 16384; // aim ~16k pixels per thread
	if (threads > max_threads) threads = max_threads;
	if (threads < 1) threads = 1;
	return threads;
}

static GdkPixbufDestroyNotify free_preview_data(guchar *pixels, gpointer data) {
	free(pixels);
	free(data);
	return FALSE;
}

#define TRYFITS(f, ...) \
	do{ \
		status = FALSE; \
		f(__VA_ARGS__, &status); \
		if(status){ \
			if(ima_data) free(ima_data); \
			fits_close_file(fp, &status); \
			return NULL; \
		} \
	} while(0)

/**
 * Create a preview of a FITS file in a GdkPixbuf (color if available)
 * @param filename
 * @return a GdkPixbuf containing the preview or NULL
 */
GdkPixbuf* get_thumbnail_from_fits(char *filename, gchar **descr) {
	fitsfile *fp;
	gchar *description = NULL;
	const int MAX_SIZE = com.pref.gui.thumbnail_size;
	float nullval = 0.;
	int naxis, dtype, stat, status, frames;

	long naxes[4];
	float *ima_data = NULL;

	TRYFITS(siril_fits_open_diskfile, &fp, filename, READONLY);

	if (siril_fits_move_first_image(fp)) {
		siril_log_message(_("Selecting the primary header failed, is the FITS file '%s' malformed?\n"), filename);
		return NULL;
	}

	TRYFITS(fits_get_img_param, fp, 4, &dtype, &naxis, naxes);

	const int w = naxes[0];
	const int h = naxes[1];
	const int n_channels = (naxis >= 3 && naxes[2] >= 3) ? 3 : 1;
	const gboolean is_color = (n_channels == 3);

	if (w <= 0 || h <= 0)
		return NULL;

	size_t sz = (size_t)w * h * n_channels;
	ima_data = malloc(sz * sizeof(float));
	if (!ima_data) {
		fits_close_file(fp, &status);
		return NULL;
	}

	TRYFITS(fits_read_img, fp, TFLOAT, 1, sz, &nullval, ima_data, &stat);

	const int x = (int) ceil((float) w / MAX_SIZE);
	const int y = (int) ceil((float) h / MAX_SIZE);
	const int pixScale = (x > y) ? x : y;   // picture scale factor
	const int Ws = w / pixScale;            // preview width
	const int Hs = h / pixScale;            // preview height

	if (fitseq_is_fitseq(filename, &frames)) {
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)\n%d %s", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), abs(dtype), frames,
				ngettext("frame", "frames", frames));
	} else {
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), abs(dtype));
	}

	/* Allocate preview_data */
	size_t prev_size = (size_t)Ws * Hs;
	float *preview_data = malloc(prev_size * n_channels * sizeof(float));
	if (!preview_data) {
		free(ima_data);
		fits_close_file(fp, &status);
		g_free(description);
		return NULL;
	}
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
#ifdef _OPENMP
#pragma omp for
#endif
		for (int i = 0; i < Hs; i++) { // cycle through blocks by lines
			int M = i * pixScale;

			for (int j = 0; j < Ws; j++) { // cycle through blocks by columns
				int N = j * pixScale;

				for (int ch = 0; ch < n_channels; ch++) {
					float sum = 0.f;
					unsigned int count = 0;

					for (int l = 0; l < pixScale && (M + l) < h; l++) {
						for (int k = 0; k < pixScale && (N + k) < w; k++) {
							int idx;
							if (is_color) {
								idx = ch * w * h + (M + l) * w + (N + k);
							} else {
								idx = (M + l) * w + (N + k);
							}
							sum += ima_data[idx];
							count++;
						}
					}

					int preview_idx = ch * Ws * Hs + i * Ws + j;
					preview_data[preview_idx] = (count > 0) ? sum / count : 0.0f;
				}
			}
		}
	}

	// Set min, max and scale, avoiding caclulations where possible
	float maxmax = 1.f;
	float minmin = 0.f;
	float scale = 1.f;
	switch (dtype) {
		case BYTE_IMG:;
			scale = INV_UCHAR_MAX_SINGLE;
			break;
		case SHORT_IMG:;
			scale = INV_USHRT_MAX_SINGLE;
			break;
		case USHORT_IMG:;
			scale = INV_USHRT_MAX_SINGLE;
			minmin = -32768.f; // min value for 16-bit signed
			break;
		default:; // FLOAT_IMG, LONG_IMG, ULONG_IMG
			/* Find per-channel min/max */
			float min_vals[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
			float max_vals[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

			for (int ch = 0; ch < n_channels; ch++) {
				for (size_t i = 0; i < prev_size; i++) {
					int idx = ch * prev_size + i;
					float val = preview_data[idx];
					if (val < min_vals[ch]) min_vals[ch] = val;
					if (val > max_vals[ch]) max_vals[ch] = val;
				}
			}
			maxmax = is_color ? fmaxf(fmaxf(max_vals[0], max_vals[1]), max_vals[2]) : max_vals[0];
			if (maxmax < 10.f) maxmax = 1.f;	// Allow maxmax to handle integer-range and JWST images but clamp
												// to typical 0-1 range otherwise, for consistent preview
			minmin = is_color ? fminf(fminf(min_vals[0], min_vals[1]), min_vals[2]) : min_vals[0];
			if (minmin > -1.f) minmin = 0.f; // Allow minmin to handle SHORT_IMG but clamp to 0.f otherwise

			if (dtype == FLOAT_IMG) {
				scale = (maxmax > 10.f) ? INV_USHRT_MAX_SINGLE : 1.f;
				break;
			}
			maxmax = is_color ? fmaxf(fmaxf(max_vals[0], max_vals[1]), max_vals[2]) : max_vals[0];
			minmin = is_color ? fminf(fminf(min_vals[0], min_vals[1]), min_vals[2]) : min_vals[0];
			scale = 1.f / (maxmax - minmin);
			break;
	}

	int num_threads = choose_num_threads(Ws, Hs, com.max_thread);
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) if (num_threads > 1)
#endif
	for (int idx = 0 ; idx < (int)(prev_size * n_channels); idx++) {
		preview_data[idx] = (preview_data[idx] - minmin) * scale;
	}

	fits *tmp = NULL;
	new_fit_image_with_data(&tmp, Ws, Hs, n_channels, DATA_FLOAT, preview_data);
	struct mtf_params mtfp[3] = {
		{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE },
		{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE },
		{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE }
	};

	find_unlinked_midtones_balance_default(tmp, mtfp);
	apply_unlinked_mtf_to_fits(tmp, tmp, mtfp);
	tmp->fdata = NULL;
	tmp->fpdata[0] = NULL;
	tmp->fpdata[1] = NULL;
	tmp->fpdata[2] = NULL;
	clearfits(tmp);
	free(tmp);

	guchar *pixbuf_data = malloc(3 * prev_size * sizeof(guchar));

	// Move this outside the loop to avoid unnecessary multiplications
	// in the loop
	int twice_prev_size = prev_size * 2;

	// Recalculate num_threads as we rely on simd in the inner loop
	num_threads = choose_num_threads(1, Hs, com.max_thread);
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) if(num_threads > 1)
#endif
	for (int i = 0; i < Hs; i++) {
		int src_row_offset  = i * Ws;
		int dest_row_offset = (Hs - 1 - i) * Ws * 3;
		for (int j = 0; j < Ws; j++) {
			int src_idx  = src_row_offset + j;
			int dest_idx = dest_row_offset + j * 3;

			if (is_color) {
				float r = preview_data[src_idx];
				float g = preview_data[prev_size + src_idx];
				float b = preview_data[twice_prev_size + src_idx];
				set_rgb(r, g, b, &pixbuf_data[dest_idx]);
			} else {
				float gray = preview_data[src_idx];
				gray2rgb(gray, &pixbuf_data[dest_idx]);
			}
		}
	}

	fits_close_file(fp, &status);
	free(ima_data);
	free(preview_data);

	GdkPixbuf *pixbuf = gdk_pixbuf_new_from_data(pixbuf_data,
			GDK_COLORSPACE_RGB,
			FALSE,
			8,
			Ws, Hs,
			Ws * 3,
			(GdkPixbufDestroyNotify) free_preview_data,
			NULL);
	*descr = description;
	return pixbuf;
}

/* verify that the parameters of the image pointed by fptr are the same as some reference values */
int check_fits_params(fitsfile *fptr, int *oldbitpix, int *oldnaxis, long *oldnaxes, gboolean relax_dimcheck) {
	int status = 0;
	long naxes[3] = { 0L };
	int bitpix = 0, naxis = -1;
	fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status);
	if (status) {
		siril_log_message(_("Opening image failed\n"));
		fits_report_error(stderr, status); /* print error message */
		return -1;
	}
	if (naxis > 3) {
		siril_log_message(_("Stacking error: images with > 3 dimensions "
					"are not supported\n"));
		return -1;
	}

	if (*oldnaxis > 0) {
		if (!relax_dimcheck && (naxis != *oldnaxis ||
				oldnaxes[0] != naxes[0] ||
				oldnaxes[1] != naxes[1] ||
				oldnaxes[2] != naxes[2])) {
			siril_log_message(_("Stacking error: input images have "
						"different sizes\n"));
			return -1;
		}
	} else {
		*oldnaxis = naxis;
		oldnaxes[0] = naxes[0];
		oldnaxes[1] = naxes[1];
		oldnaxes[2] = naxes[2];
	}

	if (*oldbitpix != 0) {
		if (bitpix != *oldbitpix) {
			siril_log_message(_("Stacking error: input images have "
						"different precision\n"));
			return -1;
		}
	} else {
		*oldbitpix = bitpix;
	}

	return 0;
}

/* NULL-terminated list of fits */
int check_loaded_fits_params(fits *ref, ...) {
	int retval = 0;
	va_list args;
	va_start(args, ref);

	const fits *test;
	while (!retval && (test = va_arg(args, fits *))) {
		if (test->bitpix != ref->bitpix ||
				test->naxis != ref->naxis) {
			retval = 1;
			break;
		}
		for (int chan = 0; chan < ref->naxis; chan++) {
			if (test->naxes[chan] != ref->naxes[chan]) {
				retval = 1;
				break;
			}
		}
	}

	va_end(args);
	return retval;
}

// f is NULL-terminated and not empty
void merge_fits_headers_to_result2(fits *result, fits **f, gboolean do_sum) {
	// input validation
	if (!f || !f[0] || !result) {
		siril_debug_print("merge_fits_headers_to_result2: No headers to merge\n");
		return;
	}
	/* copy all from the first */
	copy_fits_metadata(f[0], result);

	/* then refine the variable fields */
	gboolean found_WCS = has_wcs(f[0]);
	GDateTime *date_obs = result->keywords.date_obs;	// already referenced
	double expstart = f[0]->keywords.expstart;
	double expend = f[0]->keywords.expend;
	int image_count = 1;
	double exposure = f[0]->keywords.exposure;
	result->keywords.stackcnt = f[0]->keywords.stackcnt != DEFAULT_UINT_VALUE ? max(1, f[0]->keywords.stackcnt) : 1;
	result->keywords.livetime = f[0]->keywords.livetime > 0 ? f[0]->keywords.livetime : exposure;

	fits *current;
	while ((current = f[image_count])) {
		// take the first WCS information we find
		if (!found_WCS && has_wcs(current)) {
			int status = -1;
			result->keywords.wcslib = wcs_deepcopy(current->keywords.wcslib, &status);
			if (status)
				siril_debug_print("could not copy wcslib struct\n");
			else {
				result->keywords.wcsdata = current->keywords.wcsdata;
				found_WCS = TRUE;
			}
		}
		// set date_obs, the date of obs start, to the earliest found
		if (date_obs && current->keywords.date_obs && g_date_time_compare(date_obs, current->keywords.date_obs) == 1) {
			g_date_time_unref(date_obs);
			date_obs = g_date_time_ref(current->keywords.date_obs);
		}
		// set exposure start to the earliest found
		if (expstart > current->keywords.expstart)
			expstart = current->keywords.expstart;
		// set exposure end to the latest found
		if (expend < current->keywords.expend)
			expend = current->keywords.expend;
		// do not store conflicting filter information
		if (strcmp(result->keywords.filter, current->keywords.filter))
			strcpy(result->keywords.filter, "mixed");

		if (do_sum) {
			// add the exposure times and number of stacked images
			result->keywords.stackcnt += current->keywords.stackcnt != DEFAULT_UINT_VALUE ? max(1, current->keywords.stackcnt) : 1;
			result->keywords.livetime += current->keywords.livetime > 0 ? current->keywords.livetime : current->keywords.exposure;

			/* to add if one day we keep FITS comments: discrepancies in
			 * various fields like exposure, instrument, observer,
			 * telescope, ... */
		}
		// average exposure
		exposure += current->keywords.exposure;

		image_count++;
	}
	result->keywords.exposure = exposure / (double)image_count;
	result->keywords.date_obs = date_obs;
	result->keywords.expstart = expstart;
	result->keywords.expend = expend;
}

// NULL-terminated list of fits, given with decreasing importance
// HISTORY is not managed, neither is some conflicting information
void merge_fits_headers_to_result(fits *result, gboolean do_sum, fits *f1, ...) {
	if (!f1 || !result) return;
	// converting variadic to array of args
	va_list ap;
	va_start(ap, f1);
	int nb_fits = 1;
	fits *current;
	while ((current = va_arg(ap, fits *)))
		nb_fits++;
	va_end(ap);
	fits **array = malloc((nb_fits + 1) * sizeof(fits *));

	va_start(ap, f1);
	int i = 0;
	array[i++] = f1;
	while ((current = va_arg(ap, fits *)))
		array[i++] = current;
	va_end(ap);
	array[i] = NULL;

	merge_fits_headers_to_result2(result, array, do_sum);
	free(array);
}

typedef gsize (*StrlFunc)(char *dest, const char *src, gsize maxlen);

static void strl_with_check(char *dest, const char *src, gsize maxlen, StrlFunc strl_func) {
	gsize len = strl_func(dest, src, maxlen);
	if (len >= maxlen) {
		siril_debug_print("Exceeded FITS card length: %s, %lu, %lu\n", dest, len, maxlen);
	}
}

static int ffs2c(const char *instr, /* I - null terminated input string  */
char *outstr, /* O - null terminated quoted output string */
const int *status) /* IO - error status */
/*
 convert an input string to a quoted string. Leading spaces
 are significant.  FITS string keyword values must be at least
 8 chars long so pad out string with spaces if necessary.
 Example:   km/s ==> 'km/s    '
 Single quote characters in the input string will be replace by
 two single quote characters. e.g., o'brian ==> 'o''brian'
 */
{
	size_t len, ii, jj;

	if (*status > 0) /* inherit input status value if > 0 */
		return (*status);

	if (!instr) /* a null input pointer?? */
	{
		strcpy(outstr, "''"); /* a null FITS string */
		return (*status);
	}

	outstr[0] = '\''; /* start output string with a quote */

	len = strlen(instr);
	if (len > 68)
		len = 68; /* limit input string to 68 chars */

	for (ii = 0, jj = 1; ii < len && jj < 69; ii++, jj++) {
		outstr[jj] = instr[ii]; /* copy each char from input to output */
		if (instr[ii] == '\'') {
			jj++;
			outstr[jj] = '\''; /* duplicate any apostrophies in the input */
		}
	}

	for (; jj < 9; jj++) /* pad string so it is at least 8 chars long */
		outstr[jj] = ' ';

	if (jj == 70) /* only occurs if the last char of string was a quote */
		outstr[69] = '\0';
	else {
		outstr[jj] = '\''; /* append closing quote character */
		outstr[jj + 1] = '\0'; /* terminate the string */
	}

	return (*status);
}

void process_keyword_string_value(const char *input, char *output, gboolean condition) {
	int status = 0;
	if (!input || !output) return;

	if (condition) {
		ffs2c(input, output, &status);
	} else {
		strncpy(output, input, 70);
		output[70] = '\0';
	}
}

int updateFITSKeyword(fits *fit, const gchar *key, const gchar *newkey, const gchar *value, const gchar *comment, gboolean verbose, gboolean isfitseq) {
	char card[FLEN_CARD] = { 0 }, newcard[FLEN_CARD] = { 0 };
	char oldvalue[FLEN_VALUE] = { 0 }, oldcomment[FLEN_COMMENT] = { 0 };
	int keytype;
	void *memptr;
	size_t memsize = IOBUFLEN;
	int status = 0;
	fitsfile *fptr = NULL;

	memptr = malloc(memsize);
	if (!memptr) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fits_create_memfile(&fptr, &memptr, &memsize, IOBUFLEN, realloc, &status);
	if (status) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return 1;
	}

	if (fits_create_img(fptr, fit->bitpix, fit->naxis, fit->naxes, &status)) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return 1;
	}

	fits tmpfit = { 0 };
	copy_fits_metadata(fit, &tmpfit);

	tmpfit.fptr = fptr;
	save_fits_header(&tmpfit);

	if (key && fits_read_card(fptr, key, card, &status)) { // key is NULL if only comment
		siril_debug_print("Keyword does not exist\n");
		card[0] = '\0';
		oldcomment[0] = '\0';
		status = 0; /* reset status after error */
	} else
		siril_debug_print("%s\n", card);

	/* check if this is a protected keyword that must not be changed */
	if (*card && keyword_is_protected(card, fit)) {
		siril_log_color_message("Protected keyword cannot be modified.\n", "red");
		return 1;
	} else {
		/* Modifying keyname */
		if (newkey != NULL) {
			strl_with_check(newcard, "- ", FLEN_CARD, g_strlcpy);
			strl_with_check(newcard, key, FLEN_CARD, g_strlcat);
			strl_with_check(newcard, " ", FLEN_CARD, g_strlcat);
			strl_with_check(newcard, newkey, FLEN_CARD, g_strlcat);
			/* Deleting keyname */
		} else if (comment == NULL && value == NULL) {
			strl_with_check(newcard, "- ", FLEN_CARD, g_strlcpy);
			strl_with_check(newcard, key, FLEN_CARD, g_strlcat);
			/* Adding comment */
		} else if (key == NULL && value == NULL && comment != NULL) {
			strl_with_check(newcard, "        ", FLEN_CARD, g_strlcpy);
			strl_with_check(newcard, comment, FLEN_CARD, g_strlcat);
		} else {
			/* get the comment string */
			if (*card)
				fits_parse_value(card, oldvalue, oldcomment, &status);

			/* construct template for new keyword */
			strl_with_check(newcard, key, FLEN_CARD, g_strlcpy);
			strl_with_check(newcard, " = ", FLEN_CARD, g_strlcat);
			strl_with_check(newcard, value, FLEN_CARD, g_strlcat);

			if (*oldcomment || comment) { /* Restore comment if exist, or use new one */
				strl_with_check(newcard, " / ", FLEN_CARD, g_strlcat);
				if (comment) {
					strl_with_check(newcard, comment, FLEN_CARD, g_strlcat);
				} else {
					strl_with_check(newcard, oldcomment, FLEN_CARD, g_strlcat);
				}
			}
		}
		status = 0;
		fits_parse_template(newcard, card, &keytype, &status);

		switch (keytype) {
		case -2:
			// Rename the key: the old name is returned in the first 8 chars of card
			// and the new name is returned in characters 41-48 of card
			card[8] = '\0';
			char *new_name = card + 40;
			char *end = strchr(new_name, ' ');
			if (end)
				*end = '\0';
			card[47] = '\0';
			fits_modify_name(tmpfit.fptr, card, new_name, &status);
			remove_keyword_in_fit_keywords(key, &tmpfit); // needed to manage known keywords
			if (verbose) {
				siril_log_color_message("Keyword %s has been renamed to %s\n", "green", card, new_name);
			}
			break;
		case -1:
			// Delete the key
			fits_delete_key(tmpfit.fptr, key, &status);
			remove_keyword_in_fit_keywords(key, &tmpfit); // needed to manage known keywords
			if (verbose) {
				siril_log_color_message("Keyword %s has been removed\n", "green", key);
			}
			break;
		case 0:
			// Update the card if it already exists, otherwise append a new card
			fits_update_card(tmpfit.fptr, key, card, &status);
			if (verbose) {
				siril_log_color_message("Keyword has been changed to:\n", "green");
				siril_log_message("%s\n", card);
			}
			break;
		case 1:
			// Append the record (for HISTORY or COMMENT cards)
			//fits_write_record(tmpfit.fptr, card, &status);
			fits_write_comment(tmpfit.fptr, g_strstrip(card), &status); // here we use card for comment
			if (verbose) {
				siril_log_color_message("Comment \"%s\" has been added\n", "green", g_strstrip(card));
			}
			break;
		case 2:
		default:
			// This is for END records: do nothing
			break;
		}

		/* populate all structures */
		read_fits_header(&tmpfit);
		if (isfitseq)
			remove_all_fits_keywords(fit);
		copy_fits_metadata(&tmpfit, fit);

		if (fit->header)
			free(fit->header);
		fit->header = copy_header(&tmpfit);
	}

	if (status)
		fits_report_error(stderr, status);

	fits_close_file(tmpfit.fptr, &status);
	clearfits(&tmpfit);
	free(memptr);

	return status;
}

int associate_header_to_memfile(const char *header, fitsfile *fptr) {
	int status = 0;
	char *saveptr = NULL;

	char *header_copy = strdup(header);
	if (!header_copy) {
		return -1;
	}

	char *line = strtok_r(header_copy, "\n", &saveptr);

	while (line != NULL) {
		if (fits_write_record(fptr, line, &status)) {
			report_fits_error(status);
			free(header_copy);
			return status;
		}

		line = strtok_r(NULL, "\n", &saveptr);
	}

	free(header_copy);

	return 0;
}

int fits_parse_header_str(fits *fit, const char *header){
	void *memptr;
	size_t memsize = IOBUFLEN;
	int status = 0;
	fitsfile *fptr = NULL;

	memptr = malloc(memsize);
	if (!memptr) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	fits_create_memfile(&fptr, &memptr, &memsize, IOBUFLEN, realloc, &status);
	if (status) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return 1;
	}

	fits tmpfit = { 0 };
	tmpfit.fptr = fptr;

	associate_header_to_memfile(header, fptr);

	/* populate all structures */
	read_fits_header(&tmpfit);
	copy_fits_metadata(&tmpfit, fit);

	if (fit->header)
		free(fit->header);
	fit->header = copy_header(&tmpfit);

	fits_close_file(tmpfit.fptr, &status);
	clearfits(&tmpfit);
	free(memptr);

	return status;
}

int save_wcs_fits(fits *f, const gchar *name) {
	int status;

	if (!name)
		return 1;

	if (g_unlink(name))
		siril_debug_print("g_unlink() failed\n");
	

	status = 0;
	if (siril_fits_create_diskfile(&(f->fptr), name, &status)) {
		report_fits_error(status);
		return 1;
	}
	// Prepare the minimal header
	fits_write_key(f->fptr, TLOGICAL, "SIMPLE", &(int){1}, "conforms to FITS standard", &status);
	fits_write_key(f->fptr, TINT, "BITPIX", &(int){8}, "ASCII or bytes array", &status);
	fits_write_key(f->fptr, TINT, "NAXIS", &(int){0}, "Minimal header", &status);
	fits_write_key(f->fptr, TLOGICAL, "EXTEND", &(int){1}, "There may be FITS ext", &status);
	fits_write_key(f->fptr, TINT, "WCSAXES", &(int){2}, NULL, &status);
	fits_write_key(f->fptr, TINT, "IMAGEW", &f->rx, "Image width in pixels", &status);
	fits_write_key(f->fptr, TINT, "IMAGEH", &f->ry, "Image height in pixels", &status);
	if (f->keywords.date_obs) {
		gchar *formatted_date = date_time_to_FITS_date(f->keywords.date_obs);
		fits_update_key(f->fptr, TSTRING, "DATE-OBS", formatted_date, "YYYY-MM-DDThh:mm:ss observation start, UT", &status);
		g_free(formatted_date);
	}

	// Write the WCS data
	if (save_wcs_keywords(f)) {
		report_fits_error(status);
		fits_close_file(f->fptr, &status);
		f->fptr = NULL;
		return 1;
	}

	status = 0;
	fits_close_file(f->fptr, &status);
	if (!status) {
		siril_log_message(_("Saving WCS file %s\n"), name);
	} else {
		report_fits_error(status);
	}
	f->fptr = NULL;
	return status;
}

// writes an 32b bit distance mask to disk
int save_mask_fits(int rx, int ry, float *buffer, const gchar *name) {
	int status;
	fitsfile *fptr;
	long orig[2] = { 1L, 1L};
	long naxes[2] = { rx, ry };

	if (!name)
		return 1;

	if (g_unlink(name))
		siril_debug_print("g_unlink() failed\n");

	status = 0;
	if (siril_fits_create_diskfile(&fptr, name, &status)) {
		report_fits_error(status);
		return 1;
	}

	if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status)) {
		report_fits_error(status);
		return 1;
	}

	if (fits_write_pix(fptr, TFLOAT, orig, (size_t)(rx * ry), buffer, &status)) {
		report_fits_error(status);
		return 1;
	}

	fits_close_file(fptr, &status);
	if (!status) {
		siril_log_message(_("Saving mask file %s\n"), name);
	} else {
		report_fits_error(status);
	}
	return status;
}

// read an area form a 32b mask
// we have this dedicated function to avoid the automatic rescaling
// and dealing with all the types and headers
int read_mask_fits_area(const gchar *name, rectangle *area, int ry, float *mask) {
	int status = 0;
	fitsfile *fptr = NULL;
	int naxis = 2;
	long fpixel[2], lpixel[2], inc[2] = { 1L, 1L };
	long naxes[2] = { 1L, 1L};

	fpixel[0] = area->x + 1;        // in siril, it starts with 0
	fpixel[1] = ry - area->y - area->h + 1;
	lpixel[0] = area->x + area->w;  // with w and h at least 1, we're ok
	lpixel[1] = ry - area->y;

	if (!name)
		return 1;
	siril_fits_open_diskfile_img(&fptr, name, READONLY, &status);
	if (status) {
		report_fits_error(status);
		return 1;
	}
	fits_get_img_size(fptr, naxis, naxes, &status);
	if (status) {
		report_fits_error(status);
		return 1;
	}
	if (naxes[0] < area->w) {
		siril_debug_print("area too wide\n");
		fits_close_file(fptr, &status);
		return 1;
	}
	if (naxes[1] < area->h) {
		siril_debug_print("area too high\n");
		fits_close_file(fptr, &status);
		return 1;
	}
	fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, NULL, mask,	NULL, &status);
	if (status) {
		report_fits_error(status);
		return 1;
	}
	fits_close_file(fptr, &status);
	if (!status) {
		siril_debug_print("Read mask file %s\n", name);
	} else {
		report_fits_error(status);
	}
	return status;
}

// TODO: should not be different than read_mask_fits_area, to be improved
// read an area form a 32b drizz_weight
// we have this dedicated function to avoid the automatic rescaling
// and dealing with all the types and headers
// if we pass layer == 4, it reads the whole file
int read_drizz_fits_area(const gchar *name, int layer, rectangle *area, int ry, float *drizz) {
	if (layer < 0)
		return read_mask_fits_area(name, area, ry, drizz);
	int status = 0;
	fitsfile *fptr;
	int naxis = 3;
	long fpixel[3], lpixel[3], inc[3] = { 1L, 1L , 1L};
	long naxes[3] = { 1L, 1L, 1L};

	/* fpixel is first pixel, lpixel is last pixel, starts with value 1 */
	fpixel[0] = area->x + 1;        // in siril, it starts with 0
	fpixel[1] = ry - area->y - area->h + 1;
	fpixel[2] = layer + 1;
	lpixel[0] = area->x + area->w;  // with w and h at least 1, we're ok
	lpixel[1] = ry - area->y;
	lpixel[2] = layer + 1;

	if (!name)
		return 1;
	siril_fits_open_diskfile_img(&fptr, name, READONLY, &status);
	if (status) {
		report_fits_error(status);
		return 1;
	}
	fits_get_img_size(fptr, naxis, naxes, &status);
	if (status) {
		report_fits_error(status);
		return 1;
	}
	if (naxes[0] < area->w) {
		siril_debug_print("area too wide\n");
		fits_close_file(fptr, &status);
		return 1;
	}
	if (naxes[1] < area->h) {
		siril_debug_print("area too high\n");
		fits_close_file(fptr, &status);
		return 1;
	}
	if (layer == 4) { // read the whole file
		size_t nbdata = area->w * area->h * 3;
		fits_read_img(fptr, TFLOAT, 1, nbdata, NULL, drizz, NULL, &status);
	} else {
		fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, NULL, drizz, NULL, &status);
	}
	
	if (status) {
		report_fits_error(status);
		return 1;
	}
	fits_close_file(fptr, &status);
	if (!status) {
		siril_debug_print("Read drizz file %s\n", name);
	} else {
		report_fits_error(status);
	}
	return status;
}

// Swaps all image-data related elements of two FITS (the data pointers themselves,
// also dimensions and statistics, but not the header, keywords or fptr).

int fits_swap_image_data(fits *a, fits *b) {
	if (a == NULL || b == NULL)
		return 1;
	float *ftmp;
	WORD *tmp;
	tmp = a->data;
	a->data = b->data;
	b->data = tmp;
	for (int i = 0 ; i < 3 ; i++) {
		tmp = a->pdata[i];
		a->pdata[i] = b->pdata[i];
		b->pdata[i] = tmp;
	}
	ftmp = a->fdata;
	a->fdata = b->fdata;
	b->fdata = ftmp;
	for (int i = 0 ; i < 3 ; i++) {
		ftmp = a->fpdata[i];
		a->fpdata[i] = b->fpdata[i];
		b->fpdata[i] = ftmp;
	}
	int temp = a->bitpix;
	a->bitpix = b->bitpix;
	b->bitpix = temp;
	temp = a->orig_bitpix;
	a->orig_bitpix = b->orig_bitpix;
	b->orig_bitpix = temp;
	temp = a->naxis;
	a->naxis = b->naxis;
	b->naxis = temp;
	unsigned int uitmp = a->rx;
	a->rx = b->rx;
	b->rx = uitmp;
	uitmp = a->ry;
	a->ry = b->ry;
	b->ry = uitmp;
	long temp_naxes[3];
	memcpy(temp_naxes, a->naxes, sizeof(temp_naxes));
	memcpy(a->naxes, b->naxes, sizeof(temp_naxes));
	memcpy(b->naxes, temp_naxes, sizeof(temp_naxes));
	data_type dtmp = a->type;
	a->type = b->type;
	b->type = dtmp;;
	gboolean btmp = a->top_down;
	a->top_down = b->top_down;
	b->top_down = btmp;
	gboolean ctmp = a->debayer_checked;
	a->debayer_checked = b->debayer_checked;
	b->debayer_checked = ctmp;
	imstats **stmp = a->stats;
	a->stats = b->stats;
	b->stats = stmp;
	double fftmp = a->mini;
	a->mini = b->mini;
	b->mini = fftmp;
	fftmp = a->maxi;
	a->maxi = b->maxi;
	b->maxi = fftmp;
	float tmpf = a->neg_ratio;
	a->neg_ratio = b->neg_ratio;
	b->neg_ratio = tmpf;

	return 0;
}

/** Calculate the bayer pattern color from the row and column **/

// These interpolation routines will work for X-Trans as well as Bayer patterns

static void interpolate_nongreen_float(fits *fit, BYTE cfa[36], int cfadim) {
	uint32_t width = fit->rx;
	uint32_t height = fit->ry;
	for (int row = 0; row < height - 1; row++) {
		for (int col = 0; col < width - 1; col++) {
			if (FC_array(row, col, cfa, cfadim) == 1)
				continue;
			uint32_t index = col + row * width;
			float interp = 0.f;
			float weight = 0.f;
			for (int dy = -1; dy <= 1; dy++) {
				for (int dx = -1; dx <= 1; dx++) {
					if (dx != 0 || dy != 0) {
						// Check if the neighboring pixel is within the image bounds and green
						int nx = col + dx;
						int ny = row + dy;
						if (FC_array(nx, ny, cfa, cfadim) == 1 && nx >= 0 && nx < width && ny >= 0 && ny < height) {
							// Calculate distance from the current pixel
							float distance = dx + dy;
							// Calculate weight for the neighboring pixel
							float weight_contrib = (distance == 1) ? 1 : RECIPSQRT2;
							// Accumulate weighted green values and total weight
							interp += weight_contrib * fit->fdata[nx + ny * width];
							weight += weight_contrib;
						}
					}
				}
			}
			fit->fdata[index] = interp / weight;
		}
	}
	fit->keywords.bayer_pattern[0] = '\0'; // Mark this as no longer having a Bayer pattern
}

static void interpolate_nongreen_ushort(fits *fit, BYTE cfa[36], int cfadim) {
	uint32_t width = fit->rx;
	uint32_t height = fit->ry;
	for (int row = 0; row < height - 1; row++) {
		for (int col = 0; col < width - 1; col++) {
			if (FC_array(row, col, cfa, cfadim) == 1)
				continue;
			uint32_t index = col + row * width;
			float interp = 0.f;
			float weight = 0.f;
			for (int dy = -1; dy <= 1; dy++) {
				for (int dx = -1; dx <= 1; dx++) {
					if (dx != 0 || dy != 0) {
						// Check if the neighboring pixel is within the image bounds and green
						int nx = col + dx;
						int ny = row + dy;
						if (nx >= 0 && nx < width && ny >= 0 && ny < height && FC_array(nx, ny, cfa, cfadim) == 1) {
							// Calculate distance from the current pixel
							float distance = dx + dy;
							// Calculate weight for the neighboring pixel
							float weight_contrib = (distance == 1) ? 1 : RECIPSQRT2;
							// Accumulate weighted green values and total weight
							interp += weight_contrib * (float) fit->data[nx + ny * width];
							weight += weight_contrib;
						}
					}
				}
			}
			fit->data[index] = roundf_to_WORD(interp / weight);
		}
	}
	fit->keywords.bayer_pattern[0] = '\0'; // Mark this as no longer having a Bayer pattern
}

#undef RECIPSQRT2

void interpolate_nongreen(fits *fit) {
	int cfadim = 0;
	BYTE cfa[36];
	if (fit->naxes[2] != 1 || get_compiled_pattern(fit, cfa, &cfadim, FALSE))
		return;
	if (fit->type == DATA_FLOAT) {
		interpolate_nongreen_float(fit, cfa, cfadim);
	} else {
		interpolate_nongreen_ushort(fit, cfa, cfadim);
	}
	invalidate_stats_from_fit(fit);
	siril_debug_print("Interpolating non-green pixels\n");
}
