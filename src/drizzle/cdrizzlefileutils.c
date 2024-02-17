/*
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
 *
 */

#include <stdio.h>
#include <fitsio.h>
#include <string.h>

int save_doubles_to_fits_bintable(const char* filename, double* data, int rx, int ry, int numChannels) {
    fitsfile *fptr;
    int status = 0;
    int bitpix = DOUBLE_IMG; // Updated to DOUBLE_IMG
    long naxis = 2;
    long naxes[2] = { rx, ry };
    char *ttype[numChannels];
    char *tform[numChannels];
    char *tunit[numChannels];
    double *array[numChannels]; // Updated to double
    long firstrow = 1; // First row number to write
    int i;

    // Open FITS file
    if (fits_create_file(&fptr, filename, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Create primary array
    if (fits_create_img(fptr, bitpix, naxis, naxes, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Prepare column information for each channel
    for (i = 0; i < numChannels; i++) {
        char colName[20];
        sprintf(colName, "CHANNEL_%d", i + 1);
        ttype[i] = strdup(colName);
        tform[i] = strdup("1D"); // Updated to 1D for double
        tunit[i] = strdup(" ");
    }

    // Write data to FITS file
    if (fits_create_tbl(fptr, BINARY_TBL, 0, numChannels, ttype, tform, tunit, "BINTABLE", &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Write channel data to each column
    for (i = 0; i < numChannels; i++) {
        array[i] = data + i * rx * ry; // Pointer arithmetic to access channel data
        if (fits_write_col(fptr, TDOUBLE, i + 1, firstrow, 1, rx * ry, array[i], &status)) { // Updated to TDOUBLE
            fits_report_error(stderr, status);
            return status;
        }
    }

    // Free column information
    for (i = 0; i < numChannels; i++) {
        free(ttype[i]);
        free(tform[i]);
        free(tunit[i]);
    }

    // Close FITS file
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    return 0;
}

int read_doubles_from_fits_bintable(const char* filename, double** data, int* rx, int* ry, int* nchans) {
    fitsfile *fptr;
    int status = 0, hdutype, ncols = 3;
    long nrows;
    char *ttype = NULL;
    char *tform = NULL;
    char *tunit = NULL;
    long naxes[3];

    // Open FITS file
    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Move to the primary extension to read size
    if (fits_movabs_hdu(fptr, 1, &hdutype, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Get image dimensions
    if (fits_get_img_size(fptr, ncols, naxes, &status)) { // Fixed typo ncols
        fits_report_error(stderr, status);
        return status;
    }

    // Move to the binary table extension
    if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) { // Assuming the binary table is the second HDU
        fits_report_error(stderr, status);
        return status;
    }

    // Get number of rows and columns in the binary table
    if (fits_get_num_rows(fptr, &nrows, &status) || fits_get_num_cols(fptr, &ncols, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Assuming each column corresponds to a channel, so number of channels = number of columns
    *nchans = ncols;

    // Allocate memory for reading the data
    *data = (double*)malloc((*nchans) * nrows * sizeof(double)); // Updated to double
    if (*data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return -1;
    }

    // Read data from each column
    for (int i = 0; i < ncols; i++) {
        // Get column information
        if (fits_get_bcolparms(fptr, i + 1, ttype, tform, tunit, NULL, NULL, NULL, NULL, NULL, &status)) {
            fits_report_error(stderr, status);
            return status;
        }
        // Read data from the column
        if (fits_read_col(fptr, TDOUBLE, i + 1, 1, 1, nrows, NULL, *data + i * nrows, NULL, &status)) { // Updated to TDOUBLE
            fits_report_error(stderr, status);
            return status;
        }
    }

    printf("Column data read.\n");

    // Close FITS file
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    return 0;
}
