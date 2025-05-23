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

/* Internal image formats import and export: BMP and PPM */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "gui/progress_and_log.h"
#include "io/image_format_fits.h"
#include "io/fits_keywords.h"

#ifndef O_BINARY
#define O_BINARY 0
#endif

static int bmp32tofits48(unsigned char *rvb, unsigned long rx, unsigned long ry, fits *fit) {
	size_t datasize, i;
	WORD *rdata, *gdata, *bdata, *olddata;

	datasize = rx * ry;

	olddata = fit->data;
	if ((fit->data = realloc(fit->data, 3 * datasize * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		if (olddata)
			free(fit->data);
		return 1;
	}

	rdata = fit->pdata[RLAYER] = fit->data;
	gdata = fit->pdata[GLAYER] = fit->data + datasize;
	bdata = fit->pdata[BLAYER] = fit->data + 2 * datasize;
	for (i = 0; i < datasize; i++) {
		*bdata++ = (WORD) *rvb++;
		*gdata++ = (WORD) *rvb++;
		*rdata++ = (WORD) *rvb++;
		rvb++;
	}
	set_all_keywords_default(fit);

	fit->bitpix = fit->orig_bitpix = BYTE_IMG;
	fit->naxis = 3;
	fit->rx = rx;
	fit->ry = ry;
	fit->naxes[0] = rx;
	fit->naxes[1] = ry;
	fit->naxes[2] = 3;
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	return 0;
}

static int bmp24tofits48(unsigned char *rvb, unsigned long rx, unsigned long ry, fits *fit) {
	int i, j;
	WORD *rdata, *gdata, *bdata, *olddata;

	int padsize = (4 - (rx * 3) % 4) % 4;
	size_t newdatasize = ry * rx;

	olddata = fit->data;
	if ((fit->data = realloc(fit->data, 3 * newdatasize * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		if (olddata)
			free(fit->data);
		return 1;
	}
	rdata = fit->pdata[RLAYER] = fit->data;
	gdata = fit->pdata[GLAYER] = fit->data + newdatasize;
	bdata = fit->pdata[BLAYER] = fit->data + 2 * newdatasize;
	for (i = 0; i < ry; i++) {
		for (j = 0; j < rx; j++) {
			*bdata++ = *rvb++;
			*gdata++ = *rvb++;
			*rdata++ = *rvb++;
		}
		rvb += padsize;
	}
	set_all_keywords_default(fit);

	fit->bitpix = fit->orig_bitpix = BYTE_IMG;
	fit->naxis = 3;
	fit->rx = rx;
	fit->ry = ry;
	fit->naxes[0] = rx;
	fit->naxes[1] = ry;
	fit->naxes[2] = 3;
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	return 0;
}

static int bmp16tofits48(unsigned char *rvb, unsigned long rx, unsigned long ry, fits *fit) {
	WORD *rdata, *gdata, *bdata, *olddata;
	size_t newdatasize = ry * rx;

	olddata = fit->data;
	if ((fit->data = realloc(fit->data, 3 * newdatasize * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		if (olddata)
			free(fit->data);
		return 1;
	}
	rdata = fit->pdata[RLAYER] = fit->data;
	gdata = fit->pdata[GLAYER] = fit->data + newdatasize;
	bdata = fit->pdata[BLAYER] = fit->data + 2 * newdatasize;
	for (size_t i = 0; i < newdatasize; i++) {
		unsigned char buf0 = *rvb++;
		unsigned char buf1 = *rvb++;
		unsigned pixel_data = buf0 | buf1 << 8;

		*rdata++ = ((pixel_data & 0x7c00) >> 10) * 255.0 / 31.0 + 0.5;
		*gdata++ = ((pixel_data & 0x03e0) >> 5) * 255.0 / 31.0 + 0.5;
		*bdata++ = ((pixel_data & 0x001f) >> 0) * 255.0 / 31.0 + 0.5;
	}
	set_all_keywords_default(fit);

	fit->bitpix = fit->orig_bitpix = BYTE_IMG;
	fit->naxis = 3;
	fit->rx = rx;
	fit->ry = ry;
	fit->naxes[0] = rx;
	fit->naxes[1] = ry;
	fit->naxes[2] = 3;
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	return 0;
}

static int bmp8tofits(unsigned char *rgb, unsigned long rx, unsigned long ry, fits *fit) {
	size_t nbdata, padsize;
	int i, j;
	WORD *data, *olddata;

	padsize = (4 - (rx % 4)) % 4;
	nbdata = rx * ry;

	olddata = fit->data;
	if ((fit->data = realloc(fit->data, nbdata * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		if (olddata)
			free(fit->data);
		return 1;
	}
	data = fit->pdata[BW_LAYER] = fit->data;
	for (i = 0; i < ry; i++) {
		for (j = 0; j < rx; j++) {
			*data++ = (WORD) *rgb++;
		}
		rgb += padsize;
	}
	set_all_keywords_default(fit);

	fit->bitpix = fit->orig_bitpix = BYTE_IMG;
	fit->rx = rx;
	fit->ry = ry;
	fit->naxes[0] = rx;
	fit->naxes[1] = ry;
	fit->naxes[2] = 1;
	fit->naxis = 2;
	fit->keywords.binning_x = fit->keywords.binning_y = 1;
	return 0;
}

static void get_image_size(BYTE *header, unsigned long *width,
		unsigned long *height) {
	unsigned long bitmapinfoheader = 0;
	unsigned long lx = 0, ly = 0;
	unsigned short sx = 0, sy = 0;

	memcpy(&bitmapinfoheader, header + 14, 4);
	if (bitmapinfoheader == 12) {
		memcpy(&sx, header + 18, 2);
		memcpy(&sy, header + 20, 2);
		*width = (unsigned long) sx;
		*height = (unsigned long) sy;
	} else {
		memcpy(&lx, header + 18, 4);
		memcpy(&ly, header + 22, 4);
		*width = lx;
		*height = ly;
	}
}

/* reads a BMP image at filename `name', and stores it into the fit argument */
int readbmp(const char *name, fits *fit) {
	BYTE header[256];
	FILE *file;
	long int count;
	unsigned char *buf;
	unsigned long data_offset = 0;
	unsigned long width = 0, height = 0;
	unsigned short nbplane = 0;

	if ((file = g_fopen(name, "rb")) == NULL) {
		siril_log_color_message(_("Error opening BMP.\n"), "red");
		return -1;
	}

	if ((count = fread(header, 1, 54, file)) != 54) {
		siril_log_color_message(_("readbmp: %ld header bytes read instead of 54\n"), "red", count);
		perror("readbmp");
		fclose(file);
		return -1;
	}

	/*	memcpy(&compression, header + 30, 4);*/

	get_image_size(header, &width, &height);
	if (width < 1 || height < 1 || width > MAX_IMAGE_DIM || height > MAX_IMAGE_DIM) {
		siril_log_color_message(_("readbmp: file reports negative, zero or excessive dimensions\n"), "red");
		perror("readbmp");
		fclose(file);
		return -1;
	}
	memcpy(&nbplane, header + 28, 2);
	nbplane = nbplane / 8;
	memcpy(&data_offset, header + 10, 4);

	unsigned int padsize = (4 - (width * nbplane) % 4) % 4;
	size_t nbdata = width * height * nbplane + height * padsize;

	if (fseek(file, data_offset, SEEK_SET) == -1) {
		siril_debug_print("BMP fseek for data");
		fclose(file);
		return -1;
	}

	buf = malloc(nbdata);
	if (!buf) {
		PRINT_ALLOC_ERR;
		fclose(file);
		return -1;
	}
	unsigned long f;
	if (nbdata != (f = fread(buf, 1, nbdata, file))) {
		siril_log_color_message(_("readbmp: could not read all data: (%zu, %lu)\n"), "red", nbdata, f);
		free(buf);
		fclose(file);
		return -1;
	}
	fclose(file);

	set_all_keywords_default(fit);

	switch (nbplane) {
		case 1:
			bmp8tofits(buf, width, height, fit);
			break;
		case 2:
			bmp16tofits48(buf, width, height, fit);
			break;
		case 3:
			bmp24tofits48(buf, width, height, fit);
			break;
		case 4:
			bmp32tofits48(buf, width, height, fit);
			break;
		default:
			siril_log_color_message(_("Sorry but Siril cannot "
						"open this kind of BMP. Try to convert it before.\n"), "red");
	}
	fit->type = DATA_USHORT;
	free(buf);
	// Initialize ICC profile. As the buffer is set to NULL, this sets the
	// profile as sRGB (or Gray g22) which is what we want for 8-bit BMPs
	fits_initialize_icc(fit, NULL, 0);

	char *basename = g_path_get_basename(name);
	siril_log_message(_("Reading BMP: file %s, %ld layer(s), %ux%u pixels\n"),
			basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);
	return (int) nbplane;
}

int savebmp(const char *name, fits *fit) {
	unsigned char bmpfileheader[14] = { 'B', 'M', 	//Magic Number
		0, 0, 0, 0, 	//Size in bytes, see below
		0, 0, 0, 0, 54, 0, 0, 0	//offset
	};
	unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, //info of the header size
		0, 0, 0, 0, 	//width, see below
		0, 0, 0, 0, 	//height, see below
		1, 0, 		//number color planes
		24, 0,		//bits per pixel
		0, 0, 0, 0, 	//no compression
		0, 0, 0, 0, 	//image bits size
		0, 0, 0, 0, 	//horizontal resolution, we don't care
		0, 0, 0, 0, 	//vertical resolution, we don't care neither
		0, 0, 0, 0, 	//colors in pallete
		0, 0, 0, 0, 	//important colors
	};
	unsigned int width = fit->rx, height = fit->ry;
	double norm;

	FILE *f;

	WORD *gbuf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	float *gbuff[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	// Apply the colorspace transform
	void *buf = NULL;
	void *dest = NULL;
	gboolean src_is_float = (fit->type == DATA_FLOAT);

	// Transform the data
	if (fit->icc_profile && fit->color_managed) {
		buf = src_is_float ? (void *) fit->fdata : (void *) fit->data;
		size_t npixels = fit->rx * fit->ry;
		size_t nchans = fit->naxes[2];
		cmsUInt32Number trans_type;
		if (src_is_float) {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(float));
			trans_type = nchans == 1 ? TYPE_GRAY_FLT : TYPE_RGB_FLT_PLANAR;
		} else {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(WORD));
			trans_type = nchans == 1 ? TYPE_GRAY_16 : TYPE_RGB_16_PLANAR;
		}
		gboolean threaded = !get_thread_run();
		cmsHTRANSFORM save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, (nchans == 1 ? com.icc.mono_out : com.icc.srgb_out), trans_type, com.pref.icc.export_intent, 0);
		cmsUInt32Number data_format_size = gfit.type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
		cmsUInt32Number bytesperline = gfit.rx * data_format_size;
		cmsUInt32Number bytesperplane = npixels * data_format_size;
		cmsDoTransformLineStride(save_transform, buf, dest, gfit.rx, gfit.ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
		cmsDeleteTransform(save_transform);
		gbuf[0] = (WORD *) dest;
		gbuf[1] = (WORD *) dest + (fit->rx * fit->ry);
		gbuf[2] = (WORD *) dest + (fit->rx * fit->ry * 2);
		gbuff[0] = (float *) dest;
		gbuff[1] = (float *) dest + (fit->rx * fit->ry);
		gbuff[2] = (float *) dest + (fit->rx * fit->ry * 2);
		// No profile will be embedded as the version of BITMAPINFOHEADER
		// we use does not support them: the assumption is that the data is
		// sRGB (or possibly "Windows color space" which seems to be more or
		// less the same thing).
	}

	unsigned int padsize = (4 - (width * 3) % 4) % 4;
	size_t datasize = width * height * 3 + padsize * height;
	size_t filesize = datasize + sizeof(bmpfileheader) + sizeof(bmpinfoheader);
	int i, j;
	WORD red, blue, green;
	float redf, bluef, greenf;
	unsigned char pixel[3];

	bmpfileheader[2] = (unsigned char) (filesize);
	bmpfileheader[3] = (unsigned char) (filesize >> 8);
	bmpfileheader[4] = (unsigned char) (filesize >> 16);
	bmpfileheader[5] = (unsigned char) (filesize >> 24);

	bmpinfoheader[4] = (unsigned char) (width);
	bmpinfoheader[5] = (unsigned char) (width >> 8);
	bmpinfoheader[6] = (unsigned char) (width >> 16);
	bmpinfoheader[7] = (unsigned char) (width >> 24);

	bmpinfoheader[8] = (unsigned char) (height);
	bmpinfoheader[9] = (unsigned char) (height >> 8);
	bmpinfoheader[10] = (unsigned char) (height >> 16);
	bmpinfoheader[11] = (unsigned char) (height >> 24);

	bmpinfoheader[24] = (unsigned char) (datasize);
	bmpinfoheader[25] = (unsigned char) (datasize >> 8);
	bmpinfoheader[26] = (unsigned char) (datasize >> 16);
	bmpinfoheader[27] = (unsigned char) (datasize >> 24);

	char *filename = strdup(name);
	if (!g_str_has_suffix(filename, ".bmp")) {
		filename = str_append(&filename, ".bmp");
	}

	f = g_fopen(filename, "wb");
	if (f == NULL) {
		siril_log_color_message(_("Can't create BMP file.\n"), "red");
		free(filename);
		return 1;
	}

	norm = (fit->orig_bitpix != BYTE_IMG ?
			UCHAR_MAX_DOUBLE / USHRT_MAX_DOUBLE : 1.0);

	fwrite(bmpfileheader, sizeof(bmpfileheader), 1, f);
	fwrite(bmpinfoheader, sizeof(bmpinfoheader), 1, f);

	if (fit->type == DATA_USHORT) {
		for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				red = *gbuf[RLAYER]++;
				if (fit->naxes[2] == 3) {
					green = *gbuf[GLAYER]++;
					blue = *gbuf[BLAYER]++;
				} else {
					green = red;
					blue = red;
				}

				pixel[0] = round_to_BYTE(blue * norm); /* swap Blue and Red */
				pixel[1] = round_to_BYTE(green * norm);
				pixel[2] = round_to_BYTE(red * norm);

				fwrite(pixel, sizeof(pixel), 1, f);
			}
			if (padsize != 0)
				fwrite("0", 1, padsize, f);		//We fill the end of width with 0
		}
	} else {
		for (i = 0; i < height; i++) {
			for (j = 0; j < width; j++) {
				redf = *gbuff[RLAYER]++;
				if (fit->naxes[2] == 3) {
					greenf = *gbuff[GLAYER]++;
					bluef = *gbuff[BLAYER]++;
				} else {
					greenf = redf;
					bluef = redf;
				}

				pixel[0] = float_to_uchar_range(bluef); /* swap Blue and Red */
				pixel[1] = float_to_uchar_range(greenf);
				pixel[2] = float_to_uchar_range(redf);

				fwrite(pixel, sizeof(pixel), 1, f);
			}
			if (padsize != 0)
				fwrite("0", 1, padsize, f);	//We fill the end of width with 0
		}
	}
	fclose(f);
	siril_log_message(_("Saving BMP: file %s, %ld layer(s), %ux%u pixels\n"), filename,
			fit->naxes[2], fit->rx, fit->ry);
	free(filename);
	free(dest);
	return 0;
}

/********************* NetPBM IMAGE LOADING **********************/
/* P1	Portable bitmap	ASCII
 * P2	Portable graymap	ASCII
 * P3	Portable pixmap	ASCII
 * P4	Portable bitmap	Binary
 * P5	Portable graymap	Binary
 * P6	Portable pixmap	Binary
 */
/* This method loads a pnm or pgm binary file into the fits image passed as argument. */
int import_pnm_to_fits(const char *filename, fits *fit) {
	FILE *file;
	char buf[256];
	size_t i, j;
	int max_val;
	size_t stride;

	if ((file = g_fopen(filename, "rb")) == NULL) {
		siril_log_color_message(_("Sorry but Siril cannot open this file.\n"), "red");
		return -1;
	}
	if (fgets(buf, 256, file) == NULL) {
		perror("reading pnm file");
		fclose(file);
		return -1;
	}
	if (buf[0] != 'P' || buf[1] < '5' || buf[1] > '6' || buf[2] != '\n') {
		siril_log_color_message(
				_("Wrong magic cookie in PNM file, ASCII types and"
					" b&w bitmaps are not supported.\n"), "red");
		fclose(file);
		return -1;
	}
	if (buf[1] == '6') {
		fit->naxis = 3;
		fit->naxes[2] = 3;
	} else {
		fit->naxes[2] = 1;
		fit->naxis = 2;
	}

	do {
		if (fgets(buf, 256, file) == NULL) {
			fclose(file);
			return -1;
		}
	} while (buf[0] == '#');
	i = 0;
	while (buf[i] >= '0' && buf[i] <= '9')
		i++;
	if (i == 0) {
		fclose(file);
		return -1;
	}
	buf[i] = '\0';
	fit->rx = fit->naxes[0] = g_ascii_strtoull(buf, NULL, 10);
	j = ++i;
	while (buf[j] >= '0' && buf[j] <= '9')
		j++;
	if (j == i) {
		fclose(file);
		return -1;
	}
	if (buf[j] != '\n') {
		fclose(file);
		return -1;
	}
	buf[j] = '\0';
	fit->ry = fit->naxes[1] = g_ascii_strtoull(buf + i, NULL, 10);

	do {
		if (fgets(buf, 256, file) == NULL) {
			fclose(file);
			return -1;
		}
	} while (buf[0] == '#');
	i = 0;
	while (buf[i] >= '0' && buf[i] <= '9')
		i++;
	if (buf[i] != '\n') {
		fclose(file);
		return -1;
	}
	buf[i] = '\0';
	max_val = g_ascii_strtoll(buf, NULL, 10);
	if (max_val < UCHAR_MAX) {
		fclose(file);
		return -1;
	}
	set_all_keywords_default(fit);

	if (max_val == UCHAR_MAX) {
		/* 8-bit file */
		unsigned char *tmpbuf = NULL;
		if (fit->naxes[2] == 1)
			stride = fit->rx;
		else
			stride = fit->rx * 3;
		tmpbuf = malloc(stride * fit->ry);
		WORD *tmp = realloc(fit->data, stride * fit->ry * sizeof(WORD));
		if (tmp == NULL || tmpbuf == NULL) {
			PRINT_ALLOC_ERR;
			fclose(file);
			if (fit->data && !tmp)
				free(fit->data);
			if (tmpbuf)
				free(tmpbuf);
			if (tmp)
				free(tmp);
			fit->data = NULL;
			return -1;
		}
		fit->data = tmp;
		if (fread(tmpbuf, stride, fit->ry, file) < fit->ry) {
			siril_log_color_message(_("Error reading 8-bit PPM image data.\n"), "red");
			fclose(file);
			free(tmpbuf);
			free(fit->data);
			fit->data = NULL;
			return -1;
		}
		if (fit->naxes[2] == 3)
			rgb24bit_to_fits48bit(tmpbuf, fit, FALSE);
		else
			rgb8bit_to_fits16bit(tmpbuf, fit);
		free(tmpbuf);
		fit->bitpix = BYTE_IMG;
		fit->keywords.binning_x = fit->keywords.binning_y = 1;
		fits_flip_top_to_bottom(fit);
	} else if (max_val == USHRT_MAX || max_val == SHRT_MAX) {
		/* 16-bit file */
		if (fit->naxes[2] == 1) {
			WORD *olddata = fit->data;
			stride = fit->rx * sizeof(WORD);
			fit->data = realloc(fit->data, stride * fit->ry * sizeof(WORD));
			if (fit->data == NULL) {
				PRINT_ALLOC_ERR;
				fclose(file);
				if (olddata)
					free(olddata);
				return -1;
			}
			if (fread(fit->data, stride, fit->ry, file) < fit->ry) {
				siril_log_color_message(
						_("Error reading 16-bit gray PPM image data.\n"), "red");
				fclose(file);
				free(fit->data);
				fit->data = NULL;
				return -1;
			}
			/* change endianness in place */
			size_t nbdata = fit->rx * fit->ry;
			for (i = 0; i < nbdata; i++)
				fit->data[i] = change_endianness16(fit->data[i]);
			fit->pdata[0] = fit->data;
			fit->pdata[1] = fit->data;
			fit->pdata[2] = fit->data;

		} else {
			/* RGB 16-bit image */
			WORD *tmpbuf;
			stride = fit->rx * 3 * sizeof(WORD);
			tmpbuf = malloc(stride * fit->ry);
			WORD *tmp = realloc(fit->data, stride * fit->ry * sizeof(WORD));
			if (tmp == NULL || tmpbuf == NULL) {
				PRINT_ALLOC_ERR;
				fclose(file);
				if (fit->data && !tmp)
					free(fit->data);
				if (tmpbuf)
					free(tmpbuf);
				if (tmp)
					free(tmp);
				return -1;
			}
			fit->data = tmp;
			if (fread(tmpbuf, stride, fit->ry, file) < fit->ry) {
				siril_log_color_message(
						_("Error reading 16-bit color PPM image data.\n"), "red");
				fclose(file);
				free(tmpbuf);
				free(fit->data);
				fit->data = NULL;
				return -1;
			}
			rgb48bit_to_fits48bit(tmpbuf, fit, FALSE, TRUE);
			free(tmpbuf);
		}
		fit->bitpix = USHORT_IMG;
		fit->keywords.binning_x = fit->keywords.binning_y = 1;
		fits_flip_top_to_bottom(fit);
	} else {
		siril_log_color_message(_("Not handled max value for PNM: %d.\n"), "red",
				max_val);
		fclose(file);
		return -1;
	}
	fit->type = DATA_USHORT;
	g_snprintf(fit->keywords.row_order, FLEN_VALUE, "%s", "TOP-DOWN");

	// Initialize ICC profile. As the buffer is set to NULL, this sets the
	// profile as sRGB (or Gray g22) which is what we want
	fits_initialize_icc(fit, NULL, 0);

	fclose(file);
	char *basename = g_path_get_basename(filename);
	siril_log_message(_("Reading NetPBM: file %s, %ld layer(s), %ux%u pixels\n"),
			basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);
	return fit->naxes[2];
}

static int saveppm(const char *name, fits *fit) {
	FILE *fp = g_fopen(name, "wb");
	if (!fp) {
		siril_log_color_message(_("Error opening file %s\n"), "red", name);
		return 1;
	}
	size_t i, ndata = fit->rx * fit->ry;
	double norm;
	const char *comment = "# CREATOR : SIRIL";

	fprintf(fp, "P6\n%s\n%u %u\n%u\n", comment, fit->rx, fit->ry, USHRT_MAX);
	WORD *gbuf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	float *gbuff[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	// Apply the colorspace transform
	gboolean src_is_float = (fit->type == DATA_FLOAT);
	void *buf = NULL;
	void *dest = NULL;

	// Colorspace transform the data
	if (fit->color_managed && fit->icc_profile) {
		buf = src_is_float ? (void *) fit->fdata : (void *) fit->data;
		const size_t npixels = fit->rx * fit->ry;
		if (src_is_float) {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(float));
		} else {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(WORD));
		}
		gboolean threaded = !get_thread_run();
		cmsColorSpaceSignature sig = cmsGetColorSpace(fit->icc_profile);
		cmsUInt32Number trans_type = get_planar_formatter_type(sig, fit->type, FALSE);
		cmsHTRANSFORM save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.srgb_out, trans_type, com.pref.icc.export_intent, 0);
		cmsUInt32Number datasize = gfit.type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
		cmsUInt32Number bytesperline = gfit.rx * datasize;
		cmsUInt32Number bytesperplane = npixels * datasize;
		cmsDoTransformLineStride(save_transform, buf, dest, gfit.rx, gfit.ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
		cmsDeleteTransform(save_transform);
		gbuf[0] = (WORD *) dest;
		gbuf[1] = (WORD *) dest + (fit->rx * fit->ry);
		gbuf[2] = (WORD *) dest + (fit->rx * fit->ry * 2);
		gbuff[0] = (float *) dest;
		gbuff[1] = (float *) dest + (fit->rx * fit->ry);
		gbuff[2] = (float *) dest + (fit->rx * fit->ry * 2);
	}

	fits_flip_top_to_bottom(fit);
	norm = (fit->orig_bitpix != BYTE_IMG) ? 1.0 : USHRT_MAX_DOUBLE / UCHAR_MAX_DOUBLE;
	if (fit->type == DATA_USHORT) {
		for (i = 0; i < ndata; i++) {
			WORD color[3];
			color[0] = *gbuf[RLAYER]++ * norm;
			color[1] = *gbuf[GLAYER]++ * norm;
			color[2] = *gbuf[BLAYER]++ * norm;

			color[0] = change_endianness16(color[0]);
			color[1] = change_endianness16(color[1]);
			color[2] = change_endianness16(color[2]);
			fwrite(color, sizeof(WORD), 3, fp);
		}
	} else {
		for (i = 0; i < ndata; i++) {
			WORD color[3];
			color[0] = float_to_ushort_range(*gbuff[RLAYER]++);
			color[1] = float_to_ushort_range(*gbuff[GLAYER]++);
			color[2] = float_to_ushort_range(*gbuff[BLAYER]++);

			color[0] = change_endianness16(color[0]);
			color[1] = change_endianness16(color[1]);
			color[2] = change_endianness16(color[2]);
			fwrite(color, sizeof(WORD), 3, fp);
		}
	}
	fclose(fp);
	fits_flip_top_to_bottom(fit);
	siril_log_message(_("Saving NetPBM: file %s, %ld layer(s), %ux%u pixels\n"),
			name, fit->naxes[2], fit->rx, fit->ry);
	free(dest);
	return 0;
}

static int savepgm(const char *name, fits *fit) {
	FILE *fp;
	size_t i, ndata = fit->rx * fit->ry;
	double norm;
	WORD *gbuf = fit->pdata[RLAYER];
	float *gbuff = fit->fpdata[RLAYER];

	// Apply the colorspace transform
	void *buf = NULL;
	void *dest = NULL;
	gboolean src_is_float = (fit->type == DATA_FLOAT);

	// Colorspace transform the data
	if (fit->color_managed && fit->icc_profile) {
		buf = src_is_float ? (void *) fit->fdata : (void *) fit->data;
		const size_t npixels = fit->rx * fit->ry;
		if (src_is_float) {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(float));
		} else {
			dest = malloc(fit->rx * fit->ry * fit->naxes[2] * sizeof(WORD));
		}
		gboolean threaded = get_thread_run();
		cmsColorSpaceSignature sig = cmsGetColorSpace(fit->icc_profile);
		cmsUInt32Number trans_type = get_planar_formatter_type(sig, fit->type, FALSE);
		cmsHTRANSFORM save_transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, trans_type, com.icc.mono_out, trans_type, com.pref.icc.export_intent, 0);
		cmsUInt32Number datasize = gfit.type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
		cmsUInt32Number bytesperline = gfit.rx * datasize;
		cmsUInt32Number bytesperplane = npixels * datasize;
		cmsDoTransformLineStride(save_transform, buf, dest, gfit.rx, gfit.ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
		cmsDeleteTransform(save_transform);
		gbuf = (WORD *) dest;
		gbuff = (float *) dest;
	}

	const char *comment = "# CREATOR : SIRIL";

	fp = g_fopen(name, "wb");
	if (!fp)
		return -1;
	fprintf(fp, "P5\n%s\n%u %u\n%u\n", comment, fit->rx, fit->ry, USHRT_MAX);

	fits_flip_top_to_bottom(fit);
	norm = (fit->orig_bitpix != BYTE_IMG) ? 1.0 : USHRT_MAX_DOUBLE / UCHAR_MAX_DOUBLE;
	if (fit->type == DATA_USHORT) {
		for (i = 0; i < ndata; i++) {
			WORD tmp = *gbuf++ * norm;
			/* change endianness in place */
			WORD data[1];
			data[0] = (tmp >> 8) | (tmp << 8);
			fwrite(data, sizeof(data), 1, fp);
		}
	} else {
		for (i = 0; i < ndata; i++) {
			WORD tmp = float_to_ushort_range(*gbuff++);
			/* change endianness in place */
			WORD data[1];
			data[0] = (tmp >> 8) | (tmp << 8);
			fwrite(data, sizeof(data), 1, fp);
		}
	}
	fclose(fp);
	fits_flip_top_to_bottom(fit);
	siril_log_message(_("Saving NetPBM: file %s, %ld layer(s), %ux%u pixels\n"),
			name, fit->naxes[2], fit->rx, fit->ry);
	free(dest);
	return 0;
}

int saveNetPBM(const char *name, fits *fit) {
	int retval;
	char *filename = strdup(name);

	if (fit->naxes[2] == 1) {
		if (!g_str_has_suffix(filename, ".pgm")) {
			filename = str_append(&filename, ".pgm");
		}
		retval = savepgm(filename, fit);
	} else {
		if (!g_str_has_suffix(filename, ".ppm") && !g_str_has_suffix(filename, ".pnm")) {
			filename = str_append(&filename, ".ppm");
		}
		retval = saveppm(filename, fit);
	}
	free(filename);
	return retval;
}


static int pictofit(const WORD *buf, fits *fit) {
	WORD *data, *olddata = fit->data;

	size_t i, nbdata = fit->rx * fit->ry;
	if ((fit->data = realloc(fit->data, nbdata * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		if (olddata)
			free(olddata);
		return -1;
	}
	data = fit->pdata[BW_LAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data;
	fit->pdata[BLAYER] = fit->data;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < nbdata; i++)
		data[i] = buf[i];
	fit->bitpix = fit->orig_bitpix = SHORT_IMG;
	fit->naxes[0] = fit->rx;
	fit->naxes[1] = fit->ry;
	fit->naxes[2] = 1;
	fit->naxis = 2;
	return 0;
}

static int pictofitrgb(const WORD *buf, fits *fit) {
	WORD *data[3], *olddata = fit->data;

	size_t i, nbdata = fit->rx * fit->ry;
	if ((fit->data = realloc(fit->data, nbdata * 3 * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		if (olddata)
			free(olddata);
		return -1;
	}
	data[RLAYER] = fit->pdata[RLAYER] = fit->data;
	data[GLAYER] = fit->pdata[GLAYER] = fit->data + nbdata;
	data[BLAYER] = fit->pdata[BLAYER] = fit->data + 2 * nbdata;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < nbdata; i++)
		data[RLAYER][i] = buf[i + (nbdata * RLAYER)];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < nbdata; i++)
		data[GLAYER][i] = buf[i + (nbdata * GLAYER)];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < nbdata; i++)
		data[BLAYER][i] = buf[i + (nbdata * BLAYER)];

	fit->bitpix = fit->orig_bitpix = SHORT_IMG;
	fit->naxis = 3;
	fit->naxes[0] = fit->rx;
	fit->naxes[1] = fit->ry;
	fit->naxes[2] = 3;
	return 0;
}

static int _pic_read_header(struct pic_struct *pic_file) {
	char header[290];
	if (!pic_file || !pic_file->file)
		return -1;
	if (sizeof(header) != fread(header, 1, sizeof(header), pic_file->file)) {
		perror("read");
		return -1;
	}

	memcpy(&pic_file->magic, header, 4);

	if (pic_file->magic != 0x12231fc) {
		siril_log_color_message(_("Wrong magic cookie in PIC file. "
					"This image is not supported.\n"), "red");
		return -1;
	}

	memcpy(&pic_file->width, header + 68, 2);
	memcpy(&pic_file->height, header + 70, 2);
	assert(pic_file->width > 0 && pic_file->height > 0);
	memcpy(pic_file->bin, header + 80, 12);
	memcpy(&pic_file->nbplane, header + 92, 2);
	assert(pic_file->nbplane != 0);
	memcpy(&pic_file->hi, header + 118, 2);
	memcpy(&pic_file->lo, header + 120, 2);
	pic_file->date = g_strndup(header + 94, 10);
	pic_file->time = g_strndup(header + 104, 12);
	return 0;
}

static int _pic_close_file(struct pic_struct *pic_file) {
	int retval = 0;
	if (!pic_file)
		return retval;
	if (pic_file->file) {
		retval = fclose(pic_file->file);
		pic_file->file = NULL;
	}
	g_free(pic_file->date);
	g_free(pic_file->time);
	free(pic_file);
	return retval;
}

int readpic(const char *name, fits *fit) {
	struct pic_struct *pic_file;
	WORD *buf;
	int retval = 0;

	pic_file = calloc(1, sizeof(struct pic_struct));

	if ((pic_file->file = g_fopen(name, "rb")) == NULL) {
		siril_log_color_message(
				_("Sorry but Siril cannot open the PIC file: %s.\n"), "red", name);
		free(pic_file);
		return -1;
	}

	set_all_keywords_default(fit);

	if (_pic_read_header(pic_file)) {
		_pic_close_file(pic_file);
		return -1;
	}

	fit->rx = (unsigned int) pic_file->width;
	fit->ry = (unsigned int) pic_file->height;
	fit->keywords.binning_x = (unsigned int) pic_file->bin[4];
	fit->keywords.binning_y = (unsigned int) pic_file->bin[5];
	fit->keywords.hi = pic_file->hi;
	fit->keywords.lo = pic_file->lo;
	fit->type = DATA_USHORT;

	size_t nbdata = fit->rx * fit->ry;

	if (fseek(pic_file->file, 290, SEEK_SET)) {
		siril_log_color_message(_("Error: seek failure in file.\n"), "red");
		_pic_close_file(pic_file);
		return -1;
	}
	buf = malloc(nbdata * pic_file->nbplane * sizeof(WORD));
	if (!buf) {
		siril_log_color_message(_("Error: memory allocation failure.\n"), "red");
		_pic_close_file(pic_file);
		return -1;
	}

	if ((fread(buf, 1, nbdata * pic_file->nbplane * sizeof(WORD), pic_file->file))
			!= nbdata * pic_file->nbplane * sizeof(WORD)) {
		siril_log_color_message(_("Error: Cannot read the data\n"), "red");
		free(buf);
		_pic_close_file(pic_file);
		return -1;
	}

	switch (pic_file->nbplane) {
		case 1:
			retval = pictofit(buf, fit);
			break;
		case 3:
			retval = pictofitrgb(buf, fit);
			break;
		default:
			retval = -1;
			siril_log_color_message(_("Sorry but Siril cannot open this file.\n"), "red");
	}
	free(buf);

	if (retval) {
		_pic_close_file(pic_file);
		return -1;
	}

	char *basename = g_path_get_basename(name);
	siril_log_message(_("Reading PIC: file %s, %ld layer(s), %ux%u pixels\n"),
			basename, fit->naxes[2], fit->rx, fit->ry);
	siril_log_message("(%d,%d)-(%d,%d) - Binning %dx%d\n", pic_file->bin[0],
			pic_file->bin[1], pic_file->bin[2], pic_file->bin[3],
			fit->keywords.binning_x, fit->keywords.binning_y);

	if (pic_file->date[0] != 0x00) {
		g_strchug(pic_file->date);	// removing left white spaces if exist
		siril_log_message(_("Date (of observation): %s\n"), pic_file->date);
	}
	if (pic_file->time[0] != 0x00) {
		g_strchug(pic_file->time);	// removing left white spaces if exist
		siril_log_message(_("Time (of observation): %s\n"), pic_file->time);
	}

	// Initialize ICC profile. As the buffer is set to NULL, this sets the
	// profile as sRGB (or Gray g22) which *I think* is what we want for
	// PIC files. If anyone is still using them, please correct me if I'm
	// wrong!
	fits_initialize_icc(fit, NULL, 0);

	_pic_close_file(pic_file);
	g_free(basename);
	return retval;
}
