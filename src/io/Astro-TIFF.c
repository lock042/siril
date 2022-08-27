/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"

#include "algos/astrometry_solver.h"
#include "core/siril_date.h"

/*--------------------------------------------------------------------------*/
static int ffs2c(const char *instr, /* I - null terminated input string  */
          char *outstr,      /* O - null terminated quoted output string */
          int *status)       /* IO - error status */
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

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (!instr)            /* a null input pointer?? */
    {
       strcpy(outstr, "''");   /* a null FITS string */
       return(*status);
    }

    outstr[0] = '\'';      /* start output string with a quote */

    len = strlen(instr);
    if (len > 68)
        len = 68;    /* limit input string to 68 chars */

    for (ii=0, jj=1; ii < len && jj < 69; ii++, jj++)
    {
        outstr[jj] = instr[ii];  /* copy each char from input to output */
        if (instr[ii] == '\'')
        {
            jj++;
            outstr[jj]='\'';   /* duplicate any apostrophies in the input */
        }
    }

    for (; jj < 9; jj++)       /* pad string so it is at least 8 chars long */
        outstr[jj] = ' ';

    if (jj == 70)   /* only occurs if the last char of string was a quote */
        outstr[69] = '\0';
    else
    {
        outstr[jj] = '\'';         /* append closing quote character */
        outstr[jj+1] = '\0';          /* terminate the string */
    }

    return(*status);
}

static void siril_string_append_str(GString *str, char *value, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	ffs2c(value, valstring, &status);
	fits_make_key(key, valstring, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

static void siril_string_append_logical(GString *str, char *value, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	sprintf(valstring, "%s", value);
	fits_make_key(key, valstring, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

static void siril_string_append_comment(GString *str, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	fits_make_key(key, NULL, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

static void siril_string_append_float(GString *str, float value, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	sprintf(valstring, "%f", value);
	fits_make_key(key, valstring, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

static void siril_string_append_double(GString *str, float value, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	sprintf(valstring, "%f", value);
	fits_make_key(key, valstring, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

static void siril_string_append_int(GString *str, int value, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	sprintf(valstring, "%d", value);
	fits_make_key(key, valstring, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

static void siril_string_append_unsigned(GString *str, WORD value, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	sprintf(valstring, "%u", value);
	fits_make_key(key, valstring, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

static void siril_string_append_long(GString *str, long value, const char *key, const char *comment) {
	char valstring[FLEN_VALUE];
	char card[FLEN_CARD] = { 0 };
	int status = 0;
	sprintf(valstring, "%ld", value);
	fits_make_key(key, valstring, comment, card, &status);
	if (!status)
		g_string_append_printf(str, "%s\n", card);
}

/*
 *
    Mandatory FITS keywords are as follows:

    SIMPLE – always ”T”, indicating a FITS header.
    BITPIX – indicates array format. Options include unsigned 8-bit (8), signed 16 bit (16), signed 32 bit (32), 32-bit IEEE float (-32), and 64-bit IEEE float (-64). The standard format is 16; -64 can be read by MaxIm DL but is not written.
    NAXIS – number of axes in the data array. MaxIm DL uses 2 for monochrome images, and 3 for color images.
    NAXIS1 – corresponds to the X axis.
    NAXIS2 – corresponds to the Y axis.
    NAXIS3 – present only for color images; value is always 3 (red, green, blue color planes are present in that order).

    Optional keywords defined by the FITS standard and used in MaxIm DL:

    BSCALE – this value should be multiplied by the data array values when reading the FITS file. MaxIm DL always writes a value of 1 for this keyword.
    BZERO – this value should be added to the data array values when reading the FITS file. For 16-bit integer files, MaxIm DL writes 32768 (unless overridden by the Settings dialog).
    DATE-OBS – date of observation in the ISO standard 8601 format (Y2K compliant FITS): CCYY-MM-DDThh:mm:ss.sss. The Universal time at the start of the exposure is used. Note: the alternate format using DATE-OBS and TIME-OBS is not written, but MaxIm DL will correctly interpret it when read. The time is written to 10 ms resolution. The default behavior is to report the start of observation time, but individual camera drivers can change this.  As of Version 6.24 the DL Imaging driver sets the time to exposure midpoint.
    HISTORY – indicates the processing history of the image. This keyword may be repeated as many times as necessary.
    INSTRUME – camera information. Either user entered or obtained from the camera driver.
    OBJECT – name or catalog number of object being imaged, if available from Observatory Control window or specified by the user in Settings.
    OBSERVER – user-entered information; the observer’s name.
    TELESCOP – user-entered information about the telescope used.
 *
 *
 */

char *AstroTiff_build_header(fits *fit) {
	double bscale = 1.0, bzero = 0.0;
	int status = 0;
	GString *str = g_string_new(NULL);
	siril_string_append_logical(str, "T", "SIMPLE", "file does conform to FITS standard");
	siril_string_append_int(str, fit->bitpix, "BITPIX", "number of bits per data pixel");
	siril_string_append_int(str, fit->naxis, "NAXIS", "number of data axes");
	siril_string_append_long(str, fit->naxes[0], "NAXIS1", "length of data axis 1");
	siril_string_append_long(str, fit->naxes[1], "NAXIS2", "length of data axis 2");
	siril_string_append_long(str, fit->naxes[2], "NAXIS3", "length of data axis 3");

	siril_string_append_logical(str, "T", "EXTEND", "FITS dataset may contain extensions");
	str = g_string_append(str, "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy\n");
	str = g_string_append(str, "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H\n");

	if (fit->hi) { /* may not be initialized */
		if (fit->type == DATA_USHORT) {
			siril_string_append_unsigned(str, fit->hi, "MIPS-HI", "Upper visualization cutoff");
			siril_string_append_unsigned(str, fit->lo, "MIPS-LO", "Lower visualization cutoff");
		} else if (fit->type == DATA_FLOAT) {
			float fhi = ushort_to_float_range(fit->hi);
			float flo = ushort_to_float_range(fit->lo);
			siril_string_append_float(str, fhi, "MIPS-FHI", "Upper visualization cutoff");
			siril_string_append_float(str, flo, "MIPS-FLO", "Lower visualization cutoff");
		}
	}

	siril_string_append_double(str, bzero, "BZERO", "offset data range to that of unsigned short");
	siril_string_append_double(str, bscale, "BSCALE", "default scaling factor");

	int itmp;
	status = 0;
	char fit_date[40];
	fits_get_system_time(fit_date, &itmp, &status);
	siril_string_append_str(str, fit_date, "DATE", "UTC date that FITS file was created");

	if (fit->date_obs) {
		gchar *formatted_date = date_time_to_FITS_date(fit->date_obs);
		siril_string_append_str(str, formatted_date, "DATE-OBS", "YYYY-MM-DDThh:mm:ss observation start, UT");

		g_free(formatted_date);
	}

	siril_string_append_str(str, fit->instrume, "INSTRUME", "instrument name");
	siril_string_append_str(str, fit->observer, "OBSERVER", "observer name");
	siril_string_append_str(str, fit->telescop, "TELESCOP", "telescope used to acquire this image");

	if (!g_strcmp0(fit->row_order, "BOTTOM-UP") || !g_strcmp0(fit->row_order, "TOP-DOWN")) {
		siril_string_append_str(str, fit->row_order, "ROWORDER", "Order of the rows in image array");
	}

	if (fit->pixel_size_x > 0)
		siril_string_append_float(str, fit->pixel_size_x, "XPIXSZ", "X pixel size microns");
	if (fit->pixel_size_y > 0)
		siril_string_append_float(str, fit->pixel_size_x, "YPIXSZ", "Y pixel size microns");
	if (fit->binning_x > 0)
		siril_string_append_unsigned(str, fit->binning_x, "XBINNING", "Camera binning mode");
	if (fit->binning_y > 0)
		siril_string_append_unsigned(str, fit->binning_y, "YBINNING", "Camera binning mode");
	if (fit->bayer_pattern[0] != '\0')
		siril_string_append_str(str, fit->bayer_pattern, "BAYERPAT", "Bayer color pattern");
	if (fit->bayer_xoffset > 0)
		siril_string_append_int(str, fit->bayer_xoffset, "XBAYROFF", "X offset of Bayer array");
	if (fit->bayer_yoffset > 0)
		siril_string_append_int(str, fit->bayer_yoffset, "YBAYROFF", "Y offset of Bayer array");
	if (fit->focal_length > 0)
		siril_string_append_double(str, fit->focal_length, "FOCALLEN", "Camera focal length");
	if (fit->ccd_temp != -999.0)
		siril_string_append_double(str, fit->ccd_temp, "CCD-TEMP", "CCD temp in C");
	if (fit->exposure > 0.0)
		siril_string_append_double(str, fit->exposure, "EXPTIME", "Exposure time [s]");
	if (fit->stackcnt)
		siril_string_append_unsigned(str, fit->stackcnt, "STACKCNT", "Stack frames");
	if (fit->livetime > 0.0)
		siril_string_append_double(str, fit->livetime, "LIVETIME", "Exposure time after deadtime correction");
	if (fit->filter[0] !='\0')
		siril_string_append_str(str, fit->filter, "FILTER", "Active filter name");
	if (fit->image_type[0] !='\0')
		siril_string_append_str(str, fit->image_type, "IMAGETYP", "Type of image");
	if (fit->object[0] !='\0')
		siril_string_append_str(str, fit->object, "OBJECT", "Name of the object of interest");
	if (fit->aperture > 0.0)
		siril_string_append_double(str, fit->aperture, "APERTURE", "Aperture of the instrument");
	if (fit->iso_speed > 0.0)
		siril_string_append_double(str, fit->iso_speed, "ISOSPEED", "ISO camera setting");
	if (fit->cvf > 0.0)
		siril_string_append_double(str, fit->cvf, "CVF", "Conversion factor (e-/adu)");
	if (fit->key_gain > 0)
		siril_string_append_int(str, fit->key_gain, "GAIN", "Camera gain");
	if (fit->key_offset > 0)
		siril_string_append_int(str, fit->key_offset, "OFFSET", "Camera offset");

	if (fit->wcsdata.equinox > 0.0)
		siril_string_append_double(str, fit->wcsdata.equinox, "EQUINOX", "Equatorial equinox");
	if (fit->wcsdata.objctra[0] !='\0')
		siril_string_append_str(str, fit->wcsdata.objctra, "OBJCTRA", "Image center Right Ascension (hms)");
	if (fit->wcsdata.ra > 0.0)
		siril_string_append_double(str, fit->wcsdata.ra, "RA", "Image center Right Ascension (deg)");
	if (fit->wcsdata.objctdec[0] !='\0')
		siril_string_append_str(str, fit->wcsdata.objctdec, "OBJCTDEC", "Image center Declination (dms)");
	if (fit->wcsdata.dec > 0.0)
		siril_string_append_double(str, fit->wcsdata.dec, "DEC", "Image center Declination (deg)");
	if (fit->wcsdata.crpix[0] > 0.0)
		siril_string_append_double(str, fit->wcsdata.crpix[0], "CRPIX1", "Axis1 reference pixel");
	if (fit->wcsdata.crpix[1] > 0.0)
		siril_string_append_double(str, fit->wcsdata.crpix[1], "CRPIX2", "Axis2 reference pixel");
	if (fit->wcsdata.crval[0] > 0.0)
		siril_string_append_double(str, fit->wcsdata.crval[0], "CRVAL1", "Axis1 reference value (deg)");
	if (fit->wcsdata.crval[1] > 0.0)
		siril_string_append_double(str, fit->wcsdata.crval[1], "CRVAL2", "Axis2 reference value (deg)");

	/* check if pc matrix exists */
	if ((fit->wcsdata.pc[0][0] * fit->wcsdata.pc[1][1] - fit->wcsdata.pc[1][0] * fit->wcsdata.pc[0][1]) != 0.0) {
		if (com.pref.wcs_formalism == WCS_FORMALISM_1) {
			siril_string_append_double(str, fit->wcsdata.cdelt[0], "CDELT1", "X pixel size (deg)");
			siril_string_append_double(str, fit->wcsdata.cdelt[1], "CDELT2", "Y pixel size (deg)");
			siril_string_append_double(str, fit->wcsdata.pc[0][0], "PC1_1", "Linear transformation matrix (1, 1)");
			siril_string_append_double(str, fit->wcsdata.pc[0][1], "PC1_2", "Linear transformation matrix (1, 2)");
			siril_string_append_double(str, fit->wcsdata.pc[1][0], "PC2_1", "Linear transformation matrix (2, 1)");
			siril_string_append_double(str, fit->wcsdata.pc[1][1], "PC2_2", "Linear transformation matrix (2, 2)");

		} else {
			double cd[2][2];

			wcs_pc_to_cd(fit->wcsdata.pc, fit->wcsdata.cdelt, cd);
			siril_string_append_double(str, cd[0][0], "CD1_1", "Scale matrix (1, 1)");
			siril_string_append_double(str, cd[0][1], "CD1_2", "Scale matrix (1, 2)");
			siril_string_append_double(str, cd[1][0], "CD2_1", "Scale matrix (2, 1)");
			siril_string_append_double(str, cd[1][1], "CD2_2", "Scale matrix (2, 2)");
		}
	}
	if (fit->wcsdata.pltsolvd) {
		siril_string_append_logical(str, "T", "PLTSOLVD", "Siril internal solve");
	}

	if (fit->airmass > 0)
		siril_string_append_double(str, fit->airmass, "AIRMASS", "Airmass");
	if (fit->dft.norm[0] > 0)
		siril_string_append_double(str, fit->dft.norm[0], "DFTNORM0", "Normalisation value for channel #0");
	if (fit->dft.norm[1] > 0)
		siril_string_append_double(str, fit->dft.norm[1], "DFTNORM1", "Normalisation value for channel #1");
	if (fit->dft.norm[2] > 0)
		siril_string_append_double(str, fit->dft.norm[2], "DFTNORM2", "Normalisation value for channel 2");
	if (fit->dft.ord[0] !='\0') {
		char comment[FLEN_COMMENT];
		if (fit->dft.ord[0] == 'C')
			strcpy(comment, "Low spatial freq. are located at image center");
		else if (fit->dft.ord[0] == 'R')
			strcpy(comment, "High spatial freq. are located at image center");
		siril_string_append_str(str, fit->dft.ord, "DFTORD", comment);
	}
	if (fit->dft.type[0] !='\0') {
		char comment[FLEN_COMMENT];

		if (fit->dft.type[0] == 'S')
			strcpy(comment, "Module of a Discrete Fourier Transform");
		else if (fit->dft.type[0] == 'P')
			strcpy(comment, "Phase of a Discrete Fourier Transform");
		siril_string_append_str(str, fit->dft.type, "DFTTYPE", comment);
	}


	str = g_string_append(str, "END");
	return g_string_free(str, FALSE);
}
