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
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "io/sequence.h"
#include "io/fits_sequence.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "algos/statistics.h"
#include "algos/astrometry_solver.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "image_format_fits.h"
#include "algos/siril_wcs.h"

const char *fit_extension[] = {
		".fit",
		".fits",
		".fts"
};

static char *MIPSHI[] = {"MIPS-HI", "CWHITE", "DATAMAX", NULL };
static char *MIPSLO[] = {"MIPS-LO", "CBLACK", "DATAMIN", NULL };
static char *PIXELSIZEX[] = { "XPIXSZ", "XPIXELSZ", "PIXSIZE1", "PIXSIZEX", "XPIXSIZE", NULL };
static char *PIXELSIZEY[] = { "YPIXSZ", "YPIXELSZ", "PIXSIZE2", "PIXSIZEY", "YPIXSIZE", NULL };
static char *BINX[] = { "XBINNING", "BINX", NULL };
static char *BINY[] = { "YBINNING", "BINY", NULL };
static char *FOCAL[] = { "FOCAL", "FOCALLEN", NULL };
static char *CCD_TEMP[] = { "CCD-TEMP", "CCD_TEMP", "CCDTEMP", "TEMPERAT", "CAMTCCD", NULL };
static char *EXPOSURE[] = { "EXPTIME", "EXPOSURE", NULL };
static char *FILTER[] = {"FILTER", "FILT-1", NULL };
static char *CVF[] = { "CVF", "EGAIN", NULL };
static char *IMAGETYP[] = { "IMAGETYP", "FRAMETYP", NULL };
static char *OFFSETLEVEL[] = { "OFFSET", "BLKLEVEL", NULL };  //Used for synthetic offset
static char *NB_STACKED[] = { "STACKCNT", "NCOMBINE", NULL };
static char *SITELAT[] = { "SITELAT", "SITE-LAT", "OBSLAT", NULL };
static char *SITELONG[] = { "SITELONG", "SITE-LON", "OBSLONG", NULL };


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

static void read_fits_locdata_header(fits *fit) {
	int status = 0;
	char sitelat_dump[FLEN_VALUE] = { 0 };
	char sitelong_dump[FLEN_VALUE] = { 0 };
	double d_sitelat_dump = 0.0,  d_sitelong_dump = 0.0;
	char sitelat_dump_tmp[FLEN_VALUE] = { 0 };
	char sitelong_dump_tmp[FLEN_VALUE] = { 0 };


	__tryToFindKeywords(fit->fptr, TSTRING, SITELAT, &sitelat_dump, &status);

	if (status == 0) {
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
	}

	__tryToFindKeywords(fit->fptr, TSTRING, SITELONG, &sitelong_dump, &status);

	if (status == 0) {
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
	}

	d_sitelat_dump = parse_dms(sitelat_dump);
	d_sitelong_dump = parse_dms(sitelong_dump);
	if (isnan(d_sitelat_dump) || isnan(d_sitelong_dump)) {	// Cases SITELONG and SITELAT keyword are numbers (only NINA and Seq. Generator, for now)
		fit->sitelat = strtod(sitelat_dump, NULL);
		fit->sitelong = strtod(sitelong_dump, NULL);
	} else {
		fit->sitelat = d_sitelat_dump;
		fit->sitelong = d_sitelong_dump;
	}

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "SITEELEV", &(fit->siteelev), NULL, &status);	// Handles SITEELEV keyword cases
	if (status) {
		fit->siteelev = 0.0;	// set to 0.0 if no elevation keyword (all except NINA and Seq. Generator, for now)
	}
}

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

	fit->date_obs = FITS_date_to_date_time(date_obs);
}

void fit_get_photometry_data(fits *fit) {
	int status = 0;
	read_fits_date_obs_header(fit);
	__tryToFindKeywords(fit->fptr, TDOUBLE, EXPOSURE, &fit->exposure, &status);
}

static int fit_stats(fitsfile *fptr, float *mini, float *maxi) {
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
		if (status)
			break;
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

static int try_read_float_lo_hi(fitsfile *fptr, WORD *lo, WORD *hi) {
	float fhi, flo;
	int status = 0;
	fits_read_key(fptr, TFLOAT, "MIPS-FHI", &fhi, NULL, &status);
	if (!status) {
		*hi = float_to_ushort_range(fhi);
		status = 0;
		fits_read_key(fptr, TFLOAT, "MIPS-FLO", &flo, NULL, &status);
		if (!status) {
			*lo = float_to_ushort_range(flo);
		}
	}
	return status;
}

static void new_wcs_from_old(fits *fit, const double *crota) {
	if ((fit->wcsdata.pc[0][0] * fit->wcsdata.pc[1][1] - fit->wcsdata.pc[1][0] * fit->wcsdata.pc[0][1]) == 0.0) {
		double cd[2][2];
		int sign;

		cd[0][0] = fit->wcsdata.cdelt[0] * cos(crota[1] * M_PI / 180.0);
		if (fit->wcsdata.cdelt[0] >= 0)
			sign = +1;
		else
			sign = -1;
		cd[0][1] = fabs(fit->wcsdata.cdelt[1]) * sign * sin(crota[1] * M_PI / 180.0);
		if (fit->wcsdata.cdelt[1] >= 0)
			sign = +1;
		else
			sign = -1;
		cd[1][0] = -fabs(fit->wcsdata.cdelt[0]) * sign * sin(crota[1] * M_PI / 180.0);
		cd[1][1] = fit->wcsdata.cdelt[1] * cos(crota[1] * M_PI / 180.0);

		wcs_cd_to_pc(cd, fit->wcsdata.pc, fit->wcsdata.cdelt);
	}
}

static void load_wcs_keywords(fits *fit) {
	int status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "EQUINOX", &(fit->wcsdata.equinox), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "OBJCTRA", &(fit->wcsdata.objctra), NULL, &status);

	// RA and DEC can have values either in double or in string
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "RA", &(fit->wcsdata.ra), NULL, &status);
	if (status) {
		status = 0;
		char hms[FLEN_VALUE];
		fits_read_key(fit->fptr, TSTRING, "RA", hms, NULL, &status);
		if (!status) {
			fit->wcsdata.ra = parse_hms(hms);
			if (isnan(fit->wcsdata.ra)) fit->wcsdata.ra = 0.0;
			else siril_debug_print("read RA as HMS\n");
		}
	}

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "OBJCTDEC", &(fit->wcsdata.objctdec), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DEC", &(fit->wcsdata.dec), NULL, &status);
	if (status) {
		status = 0;
		char dms[FLEN_VALUE];
		fits_read_key(fit->fptr, TSTRING, "DEC", dms, NULL, &status);
		if (!status) {
			fit->wcsdata.dec = parse_dms(dms);
			if (isnan(fit->wcsdata.dec)) fit->wcsdata.dec = 0.0;
			else siril_debug_print("read DEC as DMS\n");
		}
	}

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRPIX1", &(fit->wcsdata.crpix[0]), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRPIX2", &(fit->wcsdata.crpix[1]), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRVAL1", &(fit->wcsdata.crval[0]), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRVAL2", &(fit->wcsdata.crval[1]), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CDELT1", &(fit->wcsdata.cdelt[0]), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CDELT2", &(fit->wcsdata.cdelt[1]), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "PC1_1", &(fit->wcsdata.pc[0][0]), NULL, &status);

	/* IF PC does not exist, check for CD */
	if (status == KEY_NO_EXIST) {
		double cd[2][2];
		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "CD1_1", &cd[0][0], NULL, &status);

		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "CD1_2", &cd[0][1], NULL, &status);

		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "CD2_1", &cd[1][0], NULL, &status);

		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "CD2_2", &cd[1][1], NULL, &status);

		if (status == 0) {
			/* we compute PC and cdelt, but before check that CD is not singular */
			if ((cd[0][0] * cd[1][1] - cd[1][0] * cd[0][1]) != 0.0) {
				wcs_cd_to_pc(cd, fit->wcsdata.pc, fit->wcsdata.cdelt);
			}
		}
	} else {
		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "PC1_2", &(fit->wcsdata.pc[0][1]), NULL, &status);

		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "PC2_1", &(fit->wcsdata.pc[1][0]), NULL, &status);

		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "PC2_2", &(fit->wcsdata.pc[1][1]), NULL, &status);
	}

	double crota[2] = { 0 };
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CROTA1", &crota[0], NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CROTA2", &crota[1], NULL, &status);

	/* Test for pc matrix: in some cases (ekos?) data are written in old (deprecated) wcs formalism
	 * So we need to test if the pc matrix is singular or not. If yes we compute cd from cdelt and crota
	 * then we compute pc */
	if (crota[0] != 0.0 && crota[1] != 0.0 &&
			(fit->wcsdata.pc[0][0] * fit->wcsdata.pc[1][1] - fit->wcsdata.pc[1][0] * fit->wcsdata.pc[0][1]) == 0.0) {

		new_wcs_from_old(fit, crota);
	}

	status = 0;
	fits_read_key(fit->fptr, TLOGICAL, "PLTSOLVD", &(fit->wcsdata.pltsolvd), fit->wcsdata.pltsolvd_comment, &status);
}


/* reading the FITS header to get useful information
 * stored in the fit, requires an opened file descriptor */
void read_fits_header(fits *fit) {
	/* about the status argument: http://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node28.html */
	int status = 0;
	double scale, zero;
	char str[FLEN_VALUE] = { 0 };

	__tryToFindKeywords(fit->fptr, TUSHORT, MIPSLO, &fit->lo, &status);
	status = 0;
	__tryToFindKeywords(fit->fptr, TUSHORT, MIPSHI, &fit->hi, &status);
	if (fit->orig_bitpix == FLOAT_IMG || fit->orig_bitpix == DOUBLE_IMG) {
		try_read_float_lo_hi(fit->fptr, &fit->lo, &fit->hi);
	}

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "BSCALE", &scale, NULL, &status);
	if (!status && 1.0 != scale) {
		siril_log_message(_("Loaded FITS file has "
					"a BSCALE different than 1 (%f)\n"), scale);
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

	/* check if file is created by Siril. If so, and if the file is FLOAT_IMG
	 * Then we are confident enough that the range is [0, 1]. So, no need to
	 * compute fit_stats
	 */
	status = 0;
	fits_read_key(fit->fptr, TSTRING, "PROGRAM", &str, NULL, &status);
	gboolean not_from_siril = status || g_ascii_strncasecmp(str, PACKAGE, strlen(PACKAGE));

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DATAMAX", &(fit->data_max), NULL, &status);

	if ((fit->bitpix == FLOAT_IMG && not_from_siril) || fit->bitpix == DOUBLE_IMG) {
		float mini, maxi;
		fit_stats(fit->fptr, &mini, &maxi);
		// override data_max if needed. In some images there are differences between max and data_max
		fit->data_max = (double) maxi;
		fit->data_min = (double) mini;
	}

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "ROWORDER", &(fit->row_order), NULL,
			&status);

	/*******************************************************************
	 * ************* CAMERA AND INSTRUMENT KEYWORDS ********************
	 * ****************************************************************/

	status = 0;
	__tryToFindKeywords(fit->fptr, TFLOAT, PIXELSIZEX, &fit->pixel_size_x, &status);
	status = 0;
	__tryToFindKeywords(fit->fptr, TFLOAT, PIXELSIZEY, &fit->pixel_size_y, &status);
#ifdef _WIN32 //TODO: remove after cfitsio is fixed
	status = 0;
	__tryToFindKeywords(fit->fptr, TINT, BINX, &fit->binning_x, &status);
	status = 0;
	__tryToFindKeywords(fit->fptr, TINT, BINY, &fit->binning_y, &status);
#else
	status = 0;
	__tryToFindKeywords(fit->fptr, TUINT, BINX, &fit->binning_x, &status);
	status = 0;
	__tryToFindKeywords(fit->fptr, TUINT, BINY, &fit->binning_y, &status);
#endif
	if (fit->binning_x <= 0) fit->binning_x = 1;
	if (fit->binning_y <= 0) fit->binning_y = 1;

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "INSTRUME", &(fit->instrume), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "TELESCOP", &(fit->telescop), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "OBSERVER", &(fit->observer), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "BAYERPAT", &(fit->bayer_pattern), NULL, &status);
	/* if value is NONE we do not want to have the keyword */
	if (status == 0 && !strcasecmp(fit->bayer_pattern, "NONE")) {
		memset(fit->bayer_pattern, 0, sizeof(char) * FLEN_VALUE);
	}

	status = 0;
	fits_read_key(fit->fptr, TINT, "XBAYROFF", &(fit->bayer_xoffset), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TINT, "YBAYROFF", &(fit->bayer_yoffset), NULL, &status);

	read_fits_date_obs_header(fit);

	status = 0;
	char date[FLEN_VALUE];
	fits_read_key(fit->fptr, TSTRING, "DATE", &date, NULL, &status);
	fit->date = FITS_date_to_date_time(date);

	status = 0;
	__tryToFindKeywords(fit->fptr, TDOUBLE, FOCAL, &fit->focal_length, &status);
	if (fit->focal_length <= 0.0) {
		/* this keyword is seen in some professional images, FLENGTH is in m. */
		double flength;
		status = 0;
		fits_read_key(fit->fptr, TDOUBLE, "FLENGTH", &flength, NULL, &status);
		if (!status) {
			fit->focal_length = flength * 1000.0; // convert m to mm
		}
	}

	status = 0;
	__tryToFindKeywords(fit->fptr, TDOUBLE, CCD_TEMP, &fit->ccd_temp, &status);
	if (status == KEY_NO_EXIST) {
		fit->ccd_temp = -999.9;
	}
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "SET-TEMP", &(fit->set_temp), NULL, &status);
	if (status == KEY_NO_EXIST) {
		fit->set_temp = -999.9;
	}
	status = 0;
	__tryToFindKeywords(fit->fptr, TDOUBLE, EXPOSURE, &fit->exposure, &status);
#ifdef _WIN32
	status = 0;
	__tryToFindKeywords(fit->fptr, TINT, NB_STACKED, &fit->stackcnt, &status);
#else
	status = 0;
	__tryToFindKeywords(fit->fptr, TUINT, NB_STACKED, &fit->stackcnt, &status);
#endif

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "LIVETIME", &(fit->livetime), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "EXPSTART", &(fit->expstart), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "EXPEND", &(fit->expend), NULL, &status);

	status = 0;
	__tryToFindKeywords(fit->fptr, TSTRING, FILTER, &fit->filter, &status);

	status = 0;
	__tryToFindKeywords(fit->fptr, TSTRING, IMAGETYP, &fit->image_type, &status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "OBJECT", &(fit->object), NULL, &status);

	read_fits_locdata_header(fit);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "APERTURE", &(fit->aperture), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "ISOSPEED", &(fit->iso_speed), NULL, &status);	// Non-standard keywords used in MaxIm DL

	status = 0;
	__tryToFindKeywords(fit->fptr, TDOUBLE, CVF, &fit->cvf, &status); // conversion gain in e-/ADU

	status = 0;
	fits_read_key(fit->fptr, TUSHORT, "GAIN", &(fit->key_gain), NULL, &status);  // Gain setting from camera

	status = 0;
	__tryToFindKeywords(fit->fptr, TUSHORT, OFFSETLEVEL, &fit->key_offset, &status); // Offset setting from camera


	/* first fill wcsdata FITS structure */
	load_wcs_keywords(fit);

	/* so now, fill the wcslib structure.
	 * Not the same because wcslib is not mandatory */
	load_WCS_from_file(fit);

	/**
	 * Others
	 */

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "AIRMASS", &(fit->airmass), NULL, &status);

	/*******************************************************************
	 * ************************* DFT KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFTNORM0", &(fit->dft.norm[0]), NULL, &status);
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFTNORM1", &(fit->dft.norm[1]), NULL, &status);
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFTNORM2", &(fit->dft.norm[2]), NULL, &status);
	status = 0;
	fits_read_key(fit->fptr, TSTRING, "DFTORD", &(fit->dft.ord), NULL, &status);
	status = 0;
	fits_read_key(fit->fptr, TSTRING, "DFTTYPE", &(fit->dft.type), NULL, &status);

	fits_read_history(fit->fptr, &(fit->history));
}

static gboolean siril_str_has_prefix(char *card, char **key) {
	int iter = 0;
	do {
		gchar *new_key = g_strdup_printf("%-8s=", key[iter]);
		if (g_str_has_prefix(card, new_key)) {
			g_free(new_key);
			return TRUE;
		}
		g_free(new_key);
		iter++;
	} while ((key[iter]));
	return FALSE;
}

/* remove leading and trailing ' and whitespaces */
static void copy_string_key(char *to, char *from) {
	if (!from || from[0]=='\0') return;
	g_strlcpy(to, from + 1, FLEN_VALUE);
	to[strlen(to) - 1] = '\0'; // remove the last ' in the string
	g_strstrip(to);
}

int fits_parse_header_string(fits *fit, gchar *header) {
	int retval = 1;
	double cd[2][2] = { { 0.0 } };
	gchar **token = g_strsplit(header, "\n", -1);
	guint nargs = g_strv_length(token);

	// we need to init these 2 fields
	fit->ccd_temp = -999.9;
	fit->set_temp = -999.9;

	for (int i = 0; i < nargs; i++) {
		int status = 0;
		char card[FLEN_CARD] = { 0 };
		char value[FLEN_VALUE] = { 0 };
		char comment[FLEN_COMMENT] = { 0 };
		int keytype;

		char *tok = token[i];
		/* in case of CR + LF, we remove \r */
		if (tok[strlen(tok) - 1] == '\r') tok[strlen(tok) - 1] = '\0';
		fits_parse_template(tok, card, &keytype, &status);
		if (status) return status;
		fits_parse_value(card, value, comment, &status);

		if (g_str_has_prefix(card, "END")) {
			retval = 0;
			break;
		} else if (siril_str_has_prefix(card, PIXELSIZEX)) {
			fit->pixel_size_x = g_ascii_strtod(value, NULL);
		} else if (siril_str_has_prefix(card, PIXELSIZEY)) {
			fit->pixel_size_y = g_ascii_strtod(value, NULL);
		} else if (siril_str_has_prefix(card, BINX)) {
			fit->binning_x = g_ascii_strtoull(value, NULL, 10);
		} else if (siril_str_has_prefix(card, BINY)) {
			fit->binning_y = g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(card, "ROWORDER=")) {
			copy_string_key(fit->row_order, value);
		} else if (g_str_has_prefix(card, "INSTRUME=")) {
			copy_string_key(fit->instrume, value);
		} else if (g_str_has_prefix(card, "TELESCOP=")) {
			copy_string_key(fit->telescop, value);
		} else if (g_str_has_prefix(card, "OBSERVER=")) {
			copy_string_key(fit->observer, value);
		} else if (g_str_has_prefix(card, "BAYERPAT=")) {
			copy_string_key(fit->bayer_pattern, value);
			if (!strcasecmp(fit->bayer_pattern, "NONE")) {
				memset(fit->bayer_pattern, 0, sizeof(char) * FLEN_VALUE);
			}
		} else if (g_str_has_prefix(card, "XBAYROFF=")) {
			fit->bayer_xoffset = g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(card, "YBAYROFF=")) {
			fit->bayer_yoffset = g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(card, "DATE    =")) {
			char tmp[FLEN_VALUE];
			copy_string_key(tmp, value);
			fit->date = FITS_date_to_date_time(tmp);
		} else if (g_str_has_prefix(card, "DATE-OBS=")) {
			char tmp[FLEN_VALUE];
			copy_string_key(tmp, value);
			fit->date_obs = FITS_date_to_date_time(tmp);
		} else if (siril_str_has_prefix(card, FOCAL)) {
			fit->focal_length = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "FLENGTH =")) {
			fit->focal_length = g_ascii_strtod(value, NULL) * 1000.0;
		} else if (siril_str_has_prefix(card, CCD_TEMP)) {
			fit->ccd_temp = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "SET-TEMP=")) {
			fit->set_temp = g_ascii_strtod(value, NULL);
		} else if (siril_str_has_prefix(card, EXPOSURE)) {
			fit->exposure = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "STACKCNT=")) {
			fit->stackcnt = g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(card, "LIVETIME=")) {
			fit->livetime = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "EXPSTART=")) {
			fit->expstart = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "EXPEND  =")) {
			fit->expend = g_ascii_strtod(value, NULL);
		} else if (siril_str_has_prefix(card, FILTER)) {
			copy_string_key(fit->filter, value);
		} else if (siril_str_has_prefix(card, IMAGETYP)) {
			copy_string_key(fit->image_type, value);
		} else if (g_str_has_prefix(card, "OBJECT  =")) {
			copy_string_key(fit->object, value);
		} else if (g_str_has_prefix(card, "APERTURE=")) {
			fit->aperture = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "ISOSPEED=")) {
			fit->iso_speed = g_ascii_strtod(value, NULL);
		} else if (siril_str_has_prefix(card, CVF)) {
			fit->cvf = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "GAIN    =")) {
			fit->key_gain = g_ascii_strtoull(value, NULL, 10);
		} else if (siril_str_has_prefix(card, OFFSETLEVEL)) {
			fit->key_offset = g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(card, "EQUINOX =")) {
			fit->wcsdata.equinox = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "OBJCTRA =")) {
			copy_string_key(fit->wcsdata.objctra, value);
		} else if (g_str_has_prefix(card, "RA      =")) {
			fit->wcsdata.ra = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "SITELAT =")) {	// It seems either SITELATE or SITE-LATE can be seen, depending on the software...
			fit->sitelat = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "SITE-LAT=")) {
			fit->sitelat = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "SITELONG=")) {	// It seems either SITELONG or SITE-LON can be seen, depending on the software...
			fit->sitelong = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "SITE-LON=")) {
			fit->sitelong = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "SITEELEV=")) {	// It seems either SITELEV or TELALT can be seen, depending on the software...
			fit->sitelong = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "OBJCTDEC=")) {
			copy_string_key(fit->wcsdata.objctdec, value);
		} else if (g_str_has_prefix(card, "DEC     =")) {
			fit->wcsdata.dec = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CRPIX1  =")) {
			fit->wcsdata.crpix[0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CRPIX2  =")) {
			fit->wcsdata.crpix[1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CRVAL1  =")) {
			fit->wcsdata.crval[0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CRVAL2  =")) {
			fit->wcsdata.crval[1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CDELT1  =")) {
			fit->wcsdata.cdelt[0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CDELT2  =")) {
			fit->wcsdata.cdelt[1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "PC1_1   =")) {
			fit->wcsdata.pc[0][0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "PC1_2   =")) {
			fit->wcsdata.pc[0][1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "PC2_1   =")) {
			fit->wcsdata.pc[1][0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "PC2_2   =")) {
			fit->wcsdata.pc[1][1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CD1_1   =")) {
			cd[0][0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CD1_2   =")) {
			cd[0][1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CD2_1   =")) {
			cd[1][0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "CD2_2   =")) {
			cd[1][1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "PLTSOLVD=")) {
			fit->wcsdata.pltsolvd = !g_strcmp0(value, "T") ? TRUE : FALSE;
			strncpy(fit->wcsdata.pltsolvd_comment, comment, FLEN_COMMENT);
		} else if (g_str_has_prefix(card, "AIRMASS =")) {
			fit->airmass = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "DFTNORM0=")) {
			fit->dft.norm[0] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "DFTNORM1=")) {
			fit->dft.norm[1] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "DFTNORM2=")) {
			fit->dft.norm[2] = g_ascii_strtod(value, NULL);
		} else if (g_str_has_prefix(card, "DFTORD  =")) {
			copy_string_key(fit->dft.ord, value);
		} else if (g_str_has_prefix(card, "DFTTYPE =")) {
			copy_string_key(fit->dft.type, value);
		}
	}
	/* we compute PC and cdelt, but before check that CD is not singular */
	if ((cd[0][0] * cd[1][1] - cd[1][0] * cd[0][1]) != 0.0) {
		wcs_cd_to_pc(cd, fit->wcsdata.pc, fit->wcsdata.cdelt);
	}
	load_WCS_from_memory(fit);

	return retval;
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
	float *data32f;

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

static void convert_floats(int bitpix, float *data, size_t nbdata) {
	size_t i;
	switch (bitpix) {
		case BYTE_IMG:
			for (i = 0; i < nbdata; i++)
				data[i] = data[i] * INV_UCHAR_MAX_SINGLE;
			break;
		case FLOAT_IMG:
			siril_log_message(_("Normalizing input data to our float range [0, 1]\n"));
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
			if (fit->lo)
				fit->lo += 32768;
			if (fit->hi)
				fit->hi += 32768;
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
				|| fit->bitpix == BYTE_IMG) || fit->data_max > 2.0) { // needed for some FLOAT_IMG
			convert_floats(fit->bitpix, fit->fdata, nbdata);
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
			if (status2 == 0 && data_max > 2.0) { // needed for some FLOAT_IMG
				convert_floats(bitpix, dest, nbdata);
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
	gchar *localefilename = g_strdup(filename);
	fits_create_diskfile(fptr, localefilename, status);
	g_free(localefilename);
	return *status;
}

static void save_wcs_keywords(fits *fit) {
	int status = 0;

	if (fit->wcsdata.equinox > 0.0) {
		fits_update_key(fit->fptr, TSTRING, "CTYPE1", "RA---TAN", "Coordinate type for the first axis", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CTYPE2", "DEC--TAN", "Coordinate type for the second axis", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CUNIT1", "deg","Unit of coordinates", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CUNIT2", "deg","Unit of coordinates", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "EQUINOX", &(fit->wcsdata.equinox),	"Equatorial equinox", &status);
	}
	status = 0;

	/* Needed for Aladin compatibility */
	if (fit->naxes[2] == 3 && com.pref.rgb_aladin) {
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CTYPE3", "RGB", "RGB image", &status);
	}

	if (fit->wcsdata.objctra[0] != '\0') {
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "OBJCTRA", &(fit->wcsdata.objctra),	"Image center Right Ascension (hms)", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "OBJCTDEC", &(fit->wcsdata.objctdec), "Image center Declination (dms)", &status);
	}
	if (fit->wcsdata.ra > 0) {
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "RA", &(fit->wcsdata.ra), "Image center Right Ascension (deg)", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "DEC", &(fit->wcsdata.dec), "Image center Declination (deg)", &status);
	}
	if (fit->wcsdata.crpix[0] != '\0') {
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRPIX1", &(fit->wcsdata.crpix[0]), "Axis1 reference pixel", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRPIX2", &(fit->wcsdata.crpix[1]), "Axis2 reference pixel", &status);
	}
	if (fit->wcsdata.crval[0] != '\0') {
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL1", &(fit->wcsdata.crval[0]), "Axis1 reference value (deg)", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL2", &(fit->wcsdata.crval[1]), "Axis2 reference value (deg)", &status);
	}

	/* check if pc matrix exists */
	if ((fit->wcsdata.pc[0][0] * fit->wcsdata.pc[1][1] - fit->wcsdata.pc[1][0] * fit->wcsdata.pc[0][1]) != 0.0) {
		if (com.pref.wcs_formalism == WCS_FORMALISM_1) {
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CDELT1", &(fit->wcsdata.cdelt[0]), "X pixel size (deg)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CDELT2", &(fit->wcsdata.cdelt[1]), "Y pixel size (deg)", &status);

			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC1_1", &(fit->wcsdata.pc[0][0]), "Linear transformation matrix (1, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC1_2", &(fit->wcsdata.pc[0][1]), "Linear transformation matrix (1, 2)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC2_1", &(fit->wcsdata.pc[1][0]), "Linear transformation matrix (2, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "PC2_2", &(fit->wcsdata.pc[1][1]), "Linear transformation matrix (2, 2)", &status);
			status = 0;
		} else {
			double cd[2][2];

			wcs_pc_to_cd(fit->wcsdata.pc, fit->wcsdata.cdelt, cd);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD1_1", &cd[0][0], "Scale matrix (1, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD1_2", &cd[0][1], "Scale matrix (1, 2)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD2_1", &cd[1][0], "Scale matrix (2, 1)", &status);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, "CD2_2", &cd[1][1], "Scale matrix (2, 2)", &status);
			status = 0;
		}
	}
	if (fit->wcsdata.pltsolvd) {
		fits_update_key(fit->fptr, TLOGICAL, "PLTSOLVD", &(fit->wcsdata.pltsolvd), fit->wcsdata.pltsolvd_comment, &status);
	}
}

void save_fits_header(fits *fit) {
	int i, status = 0;
	double zero, scale;
	char comment[FLEN_COMMENT];

	if (fit->hi) { /* may not be initialized */
		if (fit->type == DATA_USHORT) {
			fits_update_key(fit->fptr, TUSHORT, "MIPS-HI", &(fit->hi),
					"Upper visualization cutoff", &status);
			fits_update_key(fit->fptr, TUSHORT, "MIPS-LO", &(fit->lo),
					"Lower visualization cutoff", &status);
		}
		else if (fit->type == DATA_FLOAT) {
			float fhi = ushort_to_float_range(fit->hi);
			fits_update_key(fit->fptr, TFLOAT, "MIPS-FHI", &(fhi),
					"Upper visualization cutoff", &status);
			float flo = ushort_to_float_range(fit->lo);
			fits_update_key(fit->fptr, TFLOAT, "MIPS-FLO", &(flo),
					"Lower visualization cutoff", &status);
		}
	}
	// physical_value = BZERO + BSCALE * array_value
	// https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node26.html
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
	fits_update_key(fit->fptr, TDOUBLE, "BZERO", &zero,
			"offset data range to that of unsigned short", &status);

	status = 0;
	fits_update_key(fit->fptr, TDOUBLE, "BSCALE", &scale, "default scaling factor",
			&status);

	status = 0;
	if (!g_strcmp0(fit->row_order, "BOTTOM-UP") || !g_strcmp0(fit->row_order, "TOP-DOWN")) {
		fits_update_key(fit->fptr, TSTRING, "ROWORDER", &fit->row_order,
				"Order of the rows in image array", &status);
	}

	/*******************************************************************
	 * ************* CAMERA AND INSTRUMENT KEYWORDS ********************
	 * ******************** AND DATES **********************************
	 * ****************************************************************/

	status = 0;
	if (fit->instrume[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "INSTRUME", &(fit->instrume),
				"instrument name", &status);
	status = 0;
	if (fit->telescop[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "TELESCOP", &(fit->telescop),
				"telescope used to acquire this image", &status);
	status = 0;
	if (fit->observer[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "OBSERVER", &(fit->observer),
				"observer name", &status);
	status = 0;
	int itmp;
	char fit_date[40];
	fits_get_system_time(fit_date, &itmp, &status);
	fits_update_key(fit->fptr, TSTRING, "DATE", fit_date,
			"UTC date that FITS file was created", &status);

	status = 0;
	if (fit->date_obs) {
		gchar *formatted_date = date_time_to_FITS_date(fit->date_obs);
		fits_update_key(fit->fptr, TSTRING, "DATE-OBS", formatted_date,
				"YYYY-MM-DDThh:mm:ss observation start, UT", &status);

		g_free(formatted_date);
	}

	status = 0;
	if (fit->stackcnt > 0)
		fits_update_key(fit->fptr, TUINT, "STACKCNT", &(fit->stackcnt),
				"Stack frames", &status);

	status = 0;
	if (fit->exposure > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "EXPTIME", &(fit->exposure),
				"Exposure time [s]", &status);

	status = 0;
	if (fit->livetime > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "LIVETIME", &(fit->livetime),
				"Exposure time after deadtime correction", &status);

	status = 0;
	if (fit->expstart > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "EXPSTART", &(fit->expstart),
				"Exposure start time (standard Julian date)", &status);

	status = 0;
	if (fit->expend > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "EXPEND", &(fit->expend),
				"Exposure end time (standard Julian date)", &status);

	status = 0;
	if (fit->pixel_size_x > 0.)
		fits_update_key(fit->fptr, TFLOAT, "XPIXSZ", &(fit->pixel_size_x),
				"X pixel size microns", &status);
	if (fit->pixel_size_y > 0.)
		fits_update_key(fit->fptr, TFLOAT, "YPIXSZ", &(fit->pixel_size_y),
				"Y pixel size microns", &status);

	status = 0;
	if (fit->binning_x)
		fits_update_key(fit->fptr, TUINT, "XBINNING", &(fit->binning_x),
				"Camera binning mode", &status);
	if (fit->binning_y)
		fits_update_key(fit->fptr, TUINT, "YBINNING", &(fit->binning_y),
				"Camera binning mode", &status);

	status = 0;
	if (fit->focal_length > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "FOCALLEN", &(fit->focal_length),
				"Camera focal length", &status);

	status = 0;
	if (fit->ccd_temp > -999.0)
		fits_update_key(fit->fptr, TDOUBLE, "CCD-TEMP", &(fit->ccd_temp),
				"CCD temp in C", &status);

	status = 0;
	if (fit->set_temp > -999.0)
		fits_update_key(fit->fptr, TDOUBLE, "SET-TEMP", &(fit->set_temp),
				"Temperature setting in C", &status);

	status = 0;
	if (fit->filter[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "FILTER", &(fit->filter),
				"Active filter name", &status);

	status = 0;
	if (fit->image_type[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "IMAGETYP", &(fit->image_type),
				"Type of image", &status);

	status = 0;
	if (fit->object[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "OBJECT", &(fit->object),
				"Name of the object of interest", &status);

	status = 0;
	if (fit->aperture > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "APERTURE", &(fit->aperture),
				"Aperture of the instrument", &status);

	status = 0;
	if (fit->iso_speed > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "ISOSPEED", &(fit->iso_speed),
				"ISO camera setting", &status);

	status = 0;
	if (fit->bayer_pattern[0] != '\0') {
		fits_update_key(fit->fptr, TSTRING, "BAYERPAT", &(fit->bayer_pattern),
				"Bayer color pattern", &status);

		status = 0;
		fits_update_key(fit->fptr, TINT, "XBAYROFF", &(fit->bayer_xoffset),
				"X offset of Bayer array", &status);

		status = 0;
		fits_update_key(fit->fptr, TINT, "YBAYROFF", &(fit->bayer_yoffset),
				"Y offset of Bayer array", &status);
	}

	status = 0;
	if (fit->key_gain > 0.)
		fits_update_key(fit->fptr, TUSHORT, "GAIN", &(fit->key_gain),
				"Camera gain", &status);

	status = 0;
	if (fit->key_offset > 0.)
		fits_update_key(fit->fptr, TUSHORT, "OFFSET", &(fit->key_offset),
				"Camera offset", &status);

	status = 0;
	if (fit->cvf > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "CVF", &(fit->cvf),
				"Conversion factor (e-/adu)", &status);

	/**
	 * others
	 */

	status = 0;
	if (fit->airmass > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "AIRMASS", &(fit->airmass),
				"Airmass", &status);

	status = 0;
	if (fit->sitelat)
		fits_update_key(fit->fptr, TDOUBLE, "SITELAT", &(fit->sitelat),
			"[deg] Observation site latitude", &status);

	status = 0;
	if (fit->sitelong)
		fits_update_key(fit->fptr, TDOUBLE, "SITELONG", &(fit->sitelong),
				"[deg] Observation site longitude", &status);

	status = 0;
	if (fit->siteelev)
		fits_update_key(fit->fptr, TDOUBLE, "SITEELEV", &(fit->siteelev),
				"[m] Observation site elevation", &status);

	/*******************************************************************
	 * ************************* DFT KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
	if (fit->dft.type[0] != '\0') {
		if (fit->dft.type[0] == 'S')
			strcpy(comment, "Module of a Discrete Fourier Transform");
		else if (fit->dft.type[0] == 'P')
			strcpy(comment, "Phase of a Discrete Fourier Transform");
		else
			status = 1;			// should not happen
		fits_update_key(fit->fptr, TSTRING, "DFTTYPE", &(fit->dft.type), comment, &status);
	}

	status = 0;
	if (fit->dft.ord[0] != '\0') {
		if (fit->dft.ord[0] == 'C')
			strcpy(comment, "Low spatial freq. are located at image center");
		else if (fit->dft.ord[0] == 'R')
			strcpy(comment, "High spatial freq. are located at image center");
		else
			status = 1;			// should not happen
		fits_update_key(fit->fptr, TSTRING, "DFTORD", &(fit->dft.ord), comment, &status);
	}

	for (i = 0; i < fit->naxes[2]; i++) {
		if (fit->dft.norm[i] > 0.) {
			const char *str1 = "DFTNORM";
			const char *str2 = "Normalisation value for channel #";
			char key_str[FLEN_KEYWORD], comment_str[FLEN_VALUE];
			g_sprintf(key_str, "%s%d", str1, i);
			g_sprintf(comment_str, "%s%d", str2, i);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, key_str, &(fit->dft.norm[i]),
					comment_str, &status);
		}
	}

	/*******************************************************************
	 * ********************** PROGRAM KEYWORD **************************
	 * ****************************************************************/

	status = 0;
	char programm[32];

	sprintf(programm, "%c%s v%s", toupper(PACKAGE[0]), (char*) PACKAGE + 1, VERSION);
	fits_update_key(fit->fptr, TSTRING, "PROGRAM", programm,
			"Software that created this HDU", &status);

	/*******************************************************************
	 * ******************* PLATE SOLVING KEYWORDS **********************
	 * ****************************************************************/

	save_wcs_keywords(fit);

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
		for (i = 0; i < com.hist_display; i++) {
			if (com.history[i].history[0] != '\0')
				fits_write_history(fit->fptr, com.history[i].history, &status);
		}
	}
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
#ifdef _WIN32 //TODO: remove after cfitsio is fixed
	status = 0;
	__tryToFindKeywords(fptr, TINT, NB_STACKED, stack_count, &status);
#else
	status = 0;
	__tryToFindKeywords(fptr, TUINT, NB_STACKED, stack_count, &status);
#endif
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
	int status, retval = 1;
	char *name = NULL;
	gchar *basename;
	image_type imagetype;

	if (stat_file(filename, &imagetype, &name)) {
		siril_log_message(_("%s.[any_allowed_extension] not found.\n"),
				filename);
		free(name);
		return 1;
	}
	if (imagetype != TYPEFITS) {
		siril_log_message(
				_("The file %s is not a FITS file or doesn't exists with FITS extensions.\n"),
						filename);
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

	retval = read_fits_with_convert(fit, filename, force_float);
	fit->top_down = FALSE;

	if (!retval) {
		basename = g_path_get_basename(filename);
		siril_log_message(_("Reading FITS: file %s, %ld layer(s), %ux%u pixels, %d bits\n"),
				basename, fit->naxes[2], fit->rx, fit->ry, fit->type == DATA_USHORT ? 16 : 32) ;
		g_free(basename);
	}

close_readfits:
	status = 0;
	fits_close_file(fit->fptr, &status);
	return retval;
}

static int siril_fits_open_diskfile(fitsfile **fptr, const char *filename, int iomode, int *status) {
	gchar *localefilename = g_strdup(filename);
	fits_open_diskfile(fptr, localefilename, iomode, status);
	g_free(localefilename);
	return *status;
}

int siril_fits_open_diskfile_img(fitsfile **fptr, const char *filename, int iomode, int *status) {
	gchar *localefilename = g_strdup(filename);
	fits_open_diskfile(fptr, localefilename, iomode, status);
	if (!(*status)) {
		*status = siril_fits_move_first_image(*fptr);
	}
	g_free(localefilename);
	return *status;
}

GDateTime* get_date_from_fits(const gchar *filename) {
	gchar *localefilename = g_strdup(filename);
	fitsfile *fptr = NULL;
	GDateTime *date = NULL;
	int status = 0;
	fits_open_diskfile(&fptr, localefilename, READONLY, &status);
	g_free(localefilename);
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

// reset a fit data structure, deallocates everything in it and zero the data
void clearfits(fits *fit) {
	if (fit == NULL)
		return;
	if (fit->data)
		free(fit->data);
	if (fit->fdata)
		free(fit->fdata);
	if (fit->header)
		free(fit->header);
	if (fit->history)
		g_slist_free_full(fit->history, g_free);
	if (fit->date_obs)
		g_date_time_unref(fit->date_obs);
	if (fit->date)
		g_date_time_unref(fit->date);
	if (fit->stats) {
		for (int i = 0; i < fit->naxes[2]; i++)
			free_stats(fit->stats[i]);
		free(fit->stats);
	}
	free_wcs(fit, FALSE);
	memset(fit, 0, sizeof(fits));
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
				fit->data_max = (double) maxi;
				fit->data_min = (double) mini;
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

	fit->naxes[0] = area->w;
	fit->naxes[1] = area->h;
	fit->rx = fit->naxes[0];
	fit->ry = fit->naxes[1];
	fit->naxes[2] = 1;
	fit->naxis = 2;

	status = 0;
	fits_close_file(fit->fptr, &status);
	fprintf(stdout, _("Loaded partial FITS file %s\n"), filename);
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

void flip_buffer(int bitpix, void *buffer, const rectangle *area) {
	/* reverse the read data, because it's stored upside-down */
	if (get_data_type(bitpix) == DATA_FLOAT) {
		int line_size = area->w * sizeof(float);
		void *swap = malloc(line_size);
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
		WORD *buf = (WORD *)buffer;
		int i;
		for (i = 0; i < area->h/2 ; i++) {
			memcpy(swap, buf + i*area->w, line_size);
			memcpy(buf + i*area->w, buf + (area->h - i - 1)*area->w, line_size);
			memcpy(buf + (area->h - i - 1)*area->w, swap, line_size);
		}
		free(swap);
	}
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
	if (area->x < 0 || area->y < 0 || area->x >= seq->rx || area->y >= seq->ry
			|| area->w <= 0 || area->h <= 0 || area->x + area->w > seq->rx
			|| area->y + area->h > seq->ry) {
		fprintf(stderr, "partial read from FITS file has been requested outside image bounds or with invalid size\n");
		return 1;
	}

#ifdef _OPENMP
	g_assert(seq->fd_lock);
	omp_set_lock(&seq->fd_lock[index]);
#endif

	status = internal_read_partial_fits(seq->fptr[index], seq->ry, seq->bitpix, buffer, layer, area);

#ifdef _OPENMP
	omp_unset_lock(&seq->fd_lock[index]);
#endif
	if (status)
		return 1;

	flip_buffer(seq->bitpix, buffer, area);
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
	 * end of sracking with compressed images.
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
		printf("Trying to save a FITS color file with more than 3 channels?");
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

	status = 0;
	fits_close_file(f->fptr, &status);
	if (!status) {
		siril_log_message(_("Saving FITS: file %s, %ld layer(s), %ux%u pixels, %d bits\n"),
				filename, f->naxes[2], f->rx, f->ry,
				f->type == DATA_USHORT ? 16 : 32);
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
		f->lo >>= 8;
		f->hi >>= 8;
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
			WORD *data = float_buffer_to_ushort(f->fdata, f->naxes[0] * f->naxes[1] * f->naxes[2]);
			if (fits_write_pix(f->fptr, TUSHORT, orig, pixel_count, data, &status)) {
				report_fits_error(status);
				g_free(data);
				return 1;
			}
			free(data);
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
		to->header = NULL;
		to->history = NULL;
		to->date = NULL;
		to->date_obs = NULL;
#ifdef HAVE_WCSLIB
		to->wcslib = NULL;
#endif
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
	to->history = NULL;
	to->date = NULL;
	to->date_obs = NULL;
#ifdef HAVE_WCSLIB
	to->wcslib = NULL;
#endif

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
}

/* copy non-mandatory keywords from 'from' to 'to' */
void copy_fits_metadata(fits *from, fits *to) {
	to->pixel_size_x = from->pixel_size_x;
	to->pixel_size_y = from->pixel_size_y;
	to->binning_x = from->binning_x;
	to->binning_y = from->binning_y;

	// 'date' is the file creation date, automatically set on save
	if (from->date_obs) {
		if (to->date_obs)
			g_date_time_unref(to->date_obs);
		to->date_obs = g_date_time_ref(from->date_obs);
	}
	strncpy(to->filter, from->filter, FLEN_VALUE);
	strncpy(to->image_type, from->image_type, FLEN_VALUE);
	strncpy(to->object, from->object, FLEN_VALUE);
	strncpy(to->instrume, from->instrume, FLEN_VALUE);
	strncpy(to->telescop, from->telescop, FLEN_VALUE);
	strncpy(to->observer, from->observer, FLEN_VALUE);
	strncpy(to->bayer_pattern, from->bayer_pattern, FLEN_VALUE);
	strncpy(to->row_order, from->row_order, FLEN_VALUE);

	to->bayer_xoffset = from->bayer_xoffset;
	to->bayer_yoffset = from->bayer_yoffset;
	to->focal_length = from->focal_length;
	to->iso_speed = from->iso_speed;
	to->exposure = from->exposure;
	to->expstart = from->expstart;
	to->expend = from->expend;
	to->livetime = from->livetime;
	to->stackcnt = from->stackcnt;
	to->aperture = from->aperture;
	to->ccd_temp = from->ccd_temp;
	to->set_temp = from->set_temp;
	to->cvf = from->cvf;
	to->key_gain = from->key_gain;
	to->key_offset = from->key_offset;
	to->airmass = from->airmass;
	to->sitelat = from->sitelat;
	to->sitelong = from->sitelong;
	to->siteelev = from->siteelev;

	memcpy(&to->dft, &from->dft, sizeof(dft_info));
	memcpy(&to->wcsdata, &from->wcsdata, sizeof(wcs_info));
#ifdef HAVE_WCSLIB
	//wcssub()?
#endif
	// copy from->history?
}

int copy_fits_from_file(char *source, char *destination) {
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
	return savefits(filename, fit);
}

int save1fits32(const char *filename, fits *fit, int layer) {
	if (layer != RLAYER) {
		size_t nbdata = fit->naxes[0] * fit->naxes[1];
		memcpy(fit->fdata, fit->fdata + layer * nbdata, nbdata * sizeof(float));
	}
	fit->naxis = 2;
	fit->naxes[2] = 1;
	return savefits(filename, fit);
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
void rgb8bit_to_fits16bit(unsigned char *graybuf, fits *fit) {
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
void rgb48bit_to_fits48bit(WORD *rgbbuf, fits *fit, gboolean inverted,
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
}

void extract_region_from_fits(fits *from, int layer, fits *to,
		const rectangle *area) {
	if (from->type == DATA_USHORT)
		extract_region_from_fits_ushort(from, layer, to, area);
	else if (from->type == DATA_FLOAT)
		extract_region_from_fits_float(from, layer, to, area);
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
		data = malloc(npixels * nblayer * data_size);
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
	}
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

static void gray2rgb(float gray, guchar *rgb) {
	*rgb++ = (guchar) roundf_to_BYTE(255.f * gray);
	*rgb++ = (guchar) roundf_to_BYTE(255.f * gray);
	*rgb++ = (guchar) roundf_to_BYTE(255.f * gray);
}

static GdkPixbufDestroyNotify free_preview_data(guchar *pixels, gpointer data) {
	free(pixels);
	free(data);
	return FALSE;
}

static double logviz(double arg) {
	return log1p(arg);
} // for PREVIEW_LOG

#define TRYFITS(f, ...) \
	do{ \
		status = FALSE; \
		f(__VA_ARGS__, &status); \
		if(status){ \
			free(ima_data); \
			fits_close_file(fp, &status); \
			return NULL; \
		} \
	} while(0)

/**
 * Create a monochrome preview of a FITS file in a GdkPixbuf
 * @param filename
 * @return a GdkPixbuf containing the preview or NULL
 */
GdkPixbuf* get_thumbnail_from_fits(char *filename, gchar **descr) {
	fitsfile *fp;
	gchar *description;
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
	if (w <= 0 || h <= 0)
		return(NULL);

	size_t sz = w * h;
	ima_data = malloc(sz * sizeof(float));

	TRYFITS(fits_read_img, fp, TFLOAT, 1, sz, &nullval, ima_data, &stat);

	const int x = (int) ceil((float) w / MAX_SIZE);
	const int y = (int) ceil((float) h / MAX_SIZE);
	const int pixScale = (x > y) ? x : y;	// picture scale factor
	const int Ws = w / pixScale; 			// picture width in pixScale blocks
	const int Hs = h / pixScale; 			// -//- height pixScale

	const int n_channels = naxis == 3 ? naxis : 1;

	if (fitseq_is_fitseq(filename, &frames)) { // FIXME: we reopen the file in this function
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)\n%d %s\n%s", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), abs(dtype), frames,
				ngettext("frame", "frames", frames), _("(Monochrome Preview)"));
	} else {
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)\n%s", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), abs(dtype),
				_("(Monochrome Preview)"));
	}

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		// array for preview picture line
		float pix[MAX_SIZE];
#ifdef _OPENMP
#pragma omp for
#endif
		for (int i = 0; i < Hs; i++) { // cycle through a blocks by lines
			int M = i * pixScale;
			for (int j = 0; j < MAX_SIZE; j++) { // zero line buffer
				pix[j] = 0;
			}
			unsigned int m = 0; // amount of strings read in block
			for (int l = 0; l < pixScale; l++, m++) { // cycle through a block lines
				float *ptr = &ima_data[M * w];
				int N = 0; // number of column
				for (int j = 0; j < Ws; j++) { // cycle through a blocks by columns
					unsigned int n = 0;	// amount of columns read in block
					float sum = 0.f; // average intensity in block
					for (int k = 0; k < pixScale; k++, n++) { // cycle through block pixels
						if (N++ < w) // row didn't end
							sum += *ptr++; // sum[(pix-min)/wd]/n = [sum(pix)/n-min]/wd
						else
							break;
					}
					pix[j] += sum / n; //(byte / n - min)/wd;
				}
				if (++M >= h)
					break;
			}
			// fill unused picture pixels
			float *ptr = &ima_data[i * Ws];
			for (int l = 0; l < Ws; l++)
				*ptr++ = pix[l] / m;
		}
	}

	float *ptr = ima_data;
	sz = Ws * Hs;
	float max = *ptr;
	float min = max;
	float avr = 0.f;
	for (size_t i = 0; i < sz; i++, ptr++) {
		const float val = *ptr;
		max = max(max, val);
		min = min(min, val);
		avr += val;
	}
	avr /= (float) sz;

	/* use FITS keyword if available for a better visualization */
	float lo = 0.f;
	float hi = 0.f;
	status = 0;
	__tryToFindKeywords(fp, TFLOAT, MIPSLO, &lo, &status);
	status = 0;
	__tryToFindKeywords(fp, TFLOAT, MIPSHI, &hi, &status);

	if (hi != lo && hi != 0.f && abs(dtype) <= USHORT_IMG) {
		min = lo;
		max = hi;
	} else if (dtype <= FLOAT_IMG) {	// means float or double image
		WORD wlo, whi;
		if (!try_read_float_lo_hi(fp, &wlo, &whi)) {
			min = (float) wlo / USHRT_MAX_SINGLE;
			max = (float) whi / USHRT_MAX_SINGLE;
		}
	}

	float wd = max - min;
	avr = (avr - min) / wd;	// normal average by preview
	avr = -logf(avr);		// scale factor
	if (avr > 1.) {
		wd /= avr;
	}

	guchar *pixbuf_data = malloc(3 * MAX_SIZE * MAX_SIZE * sizeof(guchar));

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
	for (int i = Hs - 1; i > -1; i--) {	// fill pixbuf mirroring image by vertical
		guchar *pptr = &pixbuf_data[Ws * i * 3];
		float *p = &ima_data[(Hs - i - 1) * Ws];
		for (int j = 0; j < Ws; j++) {
			gray2rgb(logviz((*p++ - min) / wd), pptr);
			pptr += 3;
		}
	}
	fits_close_file(fp, &status);
	free(ima_data);
	GdkPixbuf *pixbuf = gdk_pixbuf_new_from_data(pixbuf_data,	// guchar* data
			GDK_COLORSPACE_RGB,	// only this supported
			FALSE,				// no alpha
			8,				// number of bits
			Ws, Hs,				// size
			Ws * 3,				// line length in bytes
			(GdkPixbufDestroyNotify) free_preview_data, // function (*GdkPixbufDestroyNotify) (guchar *pixels, gpointer data);
			NULL);
	*descr = description;
	return pixbuf;
}

/* verify that the parameters of the image pointed by fptr are the same as some reference values */
int check_fits_params(fitsfile *fptr, int *oldbitpix, int *oldnaxis, long *oldnaxes) {
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
		if (naxis != *oldnaxis ||
				oldnaxes[0] != naxes[0] ||
				oldnaxes[1] != naxes[1] ||
				oldnaxes[2] != naxes[2]) {
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

	fits *test;
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
void merge_fits_headers_to_result2(fits *result, fits **f) {
	/* copy all from the first */
	copy_fits_metadata(f[0], result);

	/* then refine the variable fields */
	gboolean found_WCS = has_wcsdata(f[0]);
	GDateTime *date_obs = result->date_obs;	// already referenced
	double expstart = f[0]->expstart;
	double expend = f[0]->expend;
	int image_count = 1;
	double exposure = f[0]->exposure;

	fits *current;
	while ((current = f[image_count])) {
		// take the first WCS information we find
		if (!found_WCS && has_wcsdata(current)) {
			result->wcsdata = current->wcsdata;
			found_WCS = TRUE;
		}
		// set date_obs, the date of obs start, to the earliest found
		if (date_obs && current->date_obs && g_date_time_compare(date_obs, current->date_obs) == 1) {
			g_date_time_unref(date_obs);
			date_obs = g_date_time_ref(current->date_obs);
		}
		// set exposure start to the earliest found
		if (expstart > current->expstart)
			expstart = current->expstart;
		// set exposure end to the latest found
		if (expend < current->expend)
			expend = current->expend;
		// do not store conflicting filter information
		if (strcmp(result->filter, current->filter))
			strcpy(result->filter, "mixed");
		// add the exposure times and number of stacked images
		result->stackcnt += current->stackcnt;
		result->livetime += current->livetime;
		// average exposure
		exposure += current->exposure;

		/* to add if one day we keep FITS comments: discrepencies in
		 * various fields like exposure, instrument, observer,
		 * telescope, ... */
		image_count++;
	}
	result->exposure = exposure / image_count;
	result->date_obs = date_obs;
	result->expstart = expstart;
	result->expend = expend;
}

// NULL-terminated list of fits, given with decreasing importance
// HISTORY is not managed, neither is some conflicting information
void merge_fits_headers_to_result(fits *result, fits *f1, ...) {
	if (!f1) return;
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

	merge_fits_headers_to_result2(result, array);
	free(array);
}
