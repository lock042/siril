/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/exif.h"
#include "gui-gtk4/histogram.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/utils.h"
#include "io/ser.h"
#include "io/image_format_fits.h"

#ifdef HAVE_LIBHEIF
#include <libheif/heif.h>
#endif
#ifdef HAVE_LIBXISF
#include "io/SirilXISFWraper.h"
#endif
#ifdef HAVE_LIBJXL
#include "io/SirilJpegXLWrapper.h"
#endif

#include "dialog_preview.h"
#include "core/proto.h"
#include "filters/mtf.h"

/* extract_thumbnail_from_fits is the raw-byte FITS thumbnail extractor;
 * its declaration is intentionally not in io/image_format_fits.h to
 * keep that core header free of gdk-pixbuf includes.  See the comment
 * in image_format_fits.h. */
extern guchar *extract_thumbnail_from_fits(const char *filename, gchar **descr,
                                            int *width_out, int *height_out);

/* Returns a malloc'd RGB888 byte buffer plus the thumbnail's dimensions
 * and a g_strdup'd description string.  Caller owns the buffer (free()).
 * Returns NULL on error.  Exposed (non-static) so the GTK4 file browser
 * can reuse it for SER previews without copy-pasting the MTF / scaling
 * logic. */
guchar *extract_thumbnail_from_ser(const char *filename, gchar **descr,
                                    int *width_out, int *height_out) {
	int MAX_SIZE = com.pref.gui.thumbnail_size;
	gchar *description = NULL;
	int i, j, k, l, N, M;
	int w, h, pixScale, Ws, Hs, n_channels, n_frames, bit;
	int sz;
	struct ser_struct ser;
	fits fit = { 0 };

	ser_init_struct(&ser);
	if (ser_open_file(filename, &ser)) {
		return NULL;
	}
	float *pix = malloc(MAX_SIZE * sizeof(float));
	float *ima_data = NULL, *ptr = NULL, byte, n, max, min, wd, avr;
	guchar *pixbuf_data = NULL;

	w = ser.image_width;
	h = ser.image_height;
	sz = w * h;

	switch (ser.color_id) {
	case SER_MONO:
		n_channels = 1;
		break;
	default:
		n_channels = 3;
	}

	ima_data = malloc(sz * n_channels * sizeof(float));
	pixbuf_data = malloc(3 * MAX_SIZE * MAX_SIZE * sizeof(guchar));

	if (ser_read_frame(&ser, 0, &fit, FALSE, FALSE)) {
		/* Without a frame, fit.pdata is NULL and the copy loops below would
		 * dereference it. Bail out cleanly. */
		ser_close_file(&ser);
		free(ima_data);
		free(pixbuf_data);
		free(pix);
		return NULL;
	}

	if (n_channels == 1) {
		for (i = 0; i < sz; i++) {
			ima_data[i] = (float)fit.pdata[RLAYER][i];
		}
	} else {
		for (i = 0; i < sz; i++) {
			ima_data[i * 3 + 0] = (float)fit.pdata[RLAYER][i];
			ima_data[i * 3 + 1] = (float)fit.pdata[GLAYER][i];
			ima_data[i * 3 + 2] = (float)fit.pdata[BLAYER][i];
		}
	}
	clearfits(&fit);

	i = (int) ceil((float) w / MAX_SIZE);
	j = (int) ceil((float) h / MAX_SIZE);
	pixScale = (i > j) ? i : j;
	if (pixScale == 0) {
		free(ima_data);
		free(pixbuf_data);
		free(pix);
		return NULL;
	}
	Ws = w / pixScale;
	Hs = h / pixScale;

	n_frames = ser.frame_count;
	bit = ser.bit_pixel_depth;

	if (n_channels == 1) {
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)\n%d %s", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), bit, n_frames,
				ngettext("frame", "frames", n_frames));
	} else {
		description = g_strdup_printf("%d x %d %s\n%d %s (%d bits)\n%d %s", w,
				h, ngettext("pixel", "pixels", h), n_channels,
				ngettext("channel", "channels", n_channels), bit, n_frames,
				ngettext("frame", "frames", n_frames));
	}

	float *pix_r = malloc(MAX_SIZE * sizeof(float));
	float *pix_g = malloc(MAX_SIZE * sizeof(float));
	float *pix_b = malloc(MAX_SIZE * sizeof(float));

	M = 0;
	for (i = 0; i < Hs; i++) {
		for (j = 0; j < MAX_SIZE; j++) {
			if (n_channels == 1) {
				pix[j] = 0;
			} else {
				pix_r[j] = 0;
				pix_g[j] = 0;
				pix_b[j] = 0;
			}
		}
		float m = 0.f;
		for (l = 0; l < pixScale; l++, m++) {
			if (n_channels == 1) {
				ptr = &ima_data[M * w];
			} else {
				ptr = &ima_data[M * w * 3];
			}
			N = 0;
			for (j = 0; j < Ws; j++) {
				n = 0.;
				if (n_channels == 1) {
					byte = 0.;
					for (k = 0; k < pixScale; k++, n++) {
						if (N++ < w)
							byte += *ptr++;
						else
							break;
					}
					if (n == 0)
						n = 1.f;
					pix[j] += byte / n;
				} else {
					float byte_r = 0., byte_g = 0., byte_b = 0.;
					for (k = 0; k < pixScale; k++, n++) {
						if (N++ < w) {
							byte_r += *ptr++;
							byte_g += *ptr++;
							byte_b += *ptr++;
						} else {
							break;
						}
					}
					if (n == 0)
						n = 1.f;
					pix_r[j] += byte_r / n;
					pix_g[j] += byte_g / n;
					pix_b[j] += byte_b / n;
				}
			}
			if (++M >= h)
				break;
		}
		/* m counts the rows actually accumulated above; if the inner loop
		 * broke on its first iteration (M reached h) the m++ never ran, so
		 * guard the division. */
		if (m == 0.f)
			m = 1.f;
		if (n_channels == 1) {
			ptr = &ima_data[i * Ws];
			for (l = 0; l < Ws; l++)
				*ptr++ = pix[l] / m;
		} else {
			for (l = 0; l < Ws; l++) {
				ima_data[(i * Ws + l) * 3 + 0] = pix_r[l] / m;
				ima_data[(i * Ws + l) * 3 + 1] = pix_g[l] / m;
				ima_data[(i * Ws + l) * 3 + 2] = pix_b[l] / m;
			}
		}
	}

	if (bit > 8) {
		int prev_size = Ws * Hs;
		float *ima_ch = NULL;
		if (n_channels > 1) {
			ima_ch = malloc(prev_size * n_channels * sizeof(float));
			if (!ima_ch) {
				ima_ch = NULL;
			} else {
				for (i = 0; i < prev_size; i++) {
					ima_ch[0 * prev_size + i] = ima_data[i * 3 + 0];
					ima_ch[1 * prev_size + i] = ima_data[i * 3 + 1];
					ima_ch[2 * prev_size + i] = ima_data[i * 3 + 2];
				}
			}
		} else {
			ima_ch = ima_data;
		}

		if (ima_ch != NULL) {
			float maxmax = 65535.f;
			float minmin = 0.f;
			float scale = 1.f / (maxmax - minmin);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
			for (int idx = 0 ; idx < prev_size * n_channels; idx++) {
				ima_ch[idx] = (ima_ch[idx] - minmin) * scale;
			}
			fits *tmp = NULL;
			new_fit_image_with_data(&tmp, Ws, Hs, n_channels, DATA_FLOAT, ima_ch);
			struct mtf_params mtfp[3] = {
				{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE },
				{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE },
				{ 0.f, 0.f, 0.f, TRUE, TRUE, TRUE }
			};
			find_unlinked_midtones_balance_default(tmp, mtfp);
			apply_unlinked_mtf_to_fits(tmp, tmp, mtfp);
			if (n_channels > 1) {
				for (i = 0; i < prev_size; i++) {
					ima_data[i * 3 + 0] = ima_ch[0 * prev_size + i];
					ima_data[i * 3 + 1] = ima_ch[1 * prev_size + i];
					ima_data[i * 3 + 2] = ima_ch[2 * prev_size + i];
				}
			}
			tmp->fdata = NULL;
			tmp->fpdata[0] = NULL;
			tmp->fpdata[1] = NULL;
			tmp->fpdata[2] = NULL;
			clearfits(tmp);
			free(tmp);
			if (n_channels > 1) {
				free(ima_ch);
			}
		}
	}

	if (n_channels == 1) {
		ptr = ima_data;
		sz = Ws * Hs;
		max = min = *ptr;
		avr = 0;
		for (i = 0; i < sz; i++, ptr++) {
			float tmp = *ptr;
			if (tmp > max) max = tmp;
			else if (tmp < min) min = tmp;
			avr += tmp;
		}
		avr /= (float) sz;
		wd = max - min;
		if (wd == 0.f) wd = 1.f;
		avr = (avr - min) / wd;
		if (avr > 1.)
			wd /= avr;
		ptr = ima_data;
		for (i = Hs - 1; i > -1; i--) {
			guchar *pptr = &pixbuf_data[Ws * i * 3];
			for (j = 0; j < Ws; j++) {
				guchar val = (guchar) roundf_to_BYTE(255.f * (*ptr - min) / wd);
				*pptr++ = val;
				*pptr++ = val;
				*pptr++ = val;
				ptr++;
			}
		}
	} else {
		float max_r = ima_data[0], min_r = ima_data[0];
		float max_g = ima_data[1], min_g = ima_data[1];
		float max_b = ima_data[2], min_b = ima_data[2];
		sz = Ws * Hs;
		for (i = 0; i < sz; i++) {
			float r = ima_data[i * 3 + 0];
			float g = ima_data[i * 3 + 1];
			float b = ima_data[i * 3 + 2];
			if (r > max_r) max_r = r; else if (r < min_r) min_r = r;
			if (g > max_g) max_g = g; else if (g < min_g) min_g = g;
			if (b > max_b) max_b = b; else if (b < min_b) min_b = b;
		}
		float wd_r = max_r - min_r;
		float wd_g = max_g - min_g;
		float wd_b = max_b - min_b;
		if (wd_r == 0.f) wd_r = 1.f;
		if (wd_g == 0.f) wd_g = 1.f;
		if (wd_b == 0.f) wd_b = 1.f;
		for (i = Hs - 1; i > -1; i--) {
			guchar *pptr = &pixbuf_data[Ws * i * 3];
			for (j = 0; j < Ws; j++) {
				int idx = ((Hs - 1 - i) * Ws + j) * 3;
				*pptr++ = (guchar) roundf_to_BYTE(255.f * (ima_data[idx + 0] - min_r) / wd_r);
				*pptr++ = (guchar) roundf_to_BYTE(255.f * (ima_data[idx + 1] - min_g) / wd_g);
				*pptr++ = (guchar) roundf_to_BYTE(255.f * (ima_data[idx + 2] - min_b) / wd_b);
			}
		}
	}

	ser_close_file(&ser);
	free(ima_data);
	free(pix);
	free(pix_r);
	free(pix_g);
	free(pix_b);
	*descr = description;
	if (width_out)  *width_out  = Ws;
	if (height_out) *height_out = Hs;
	return pixbuf_data;
}
