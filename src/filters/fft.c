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

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>
#include <fftw3.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "io/image_format_fits.h"
#include "algos/statistics.h"

#include "fft.h"

enum {
	TYPE_CENTERED,
	TYPE_REGULAR
};

static unsigned strategy = 0;

static void fft_to_spectra(fftwf_complex *frequency_repr, float *as, float *ps,
		float *maxi, size_t nbdata) {
	*maxi = 0.0;

	for (size_t i = 0; i < nbdata; i++) {
		float r = crealf(frequency_repr[i]);
		float im = cimagf(frequency_repr[i]);
		as[i] = hypotf(r, im);
		ps[i] = atan2f(im, r);
		*maxi = max(as[i], *maxi);
	}
}

static void fft_to_freq(fftwf_complex *frequency_repr, const float *as, const float *ps,
		size_t nbdata) {
	for (size_t i = 0; i < nbdata; i++) {
		frequency_repr[i] = as[i] * (cosf(ps[i]) + I * sinf(ps[i]));
	}
}

static void change_symmetry_forward(unsigned int width, unsigned int height,
		unsigned int i, unsigned int j, unsigned int *x, unsigned int *y) {

	*x = i + width / 2;
	if (*x >= width) {
		*x = *x - width;
	}
	*y = j + height / 2;
	if (*y >= height) {
		*y = *y - height;
	}
}

static void change_symmetry_backward(unsigned int width, unsigned int height,
		unsigned int i, unsigned int j, unsigned int *x, unsigned int *y) {

	*x = i + width - width / 2;
	if (*x >= width) {
		*x = *x - width;
	}
	*y = j + height - height / 2;
	if (*y >= height) {
		*y = *y - height;
	}
}

static void centered_ushort(WORD *buf, unsigned int width, unsigned int height,
		int type) {
	unsigned int i, j;

	WORD *temp = malloc(width * height * sizeof(WORD));
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			unsigned int x = i;
			unsigned int y = j;
			if (type == FFTW_FORWARD) {
				change_symmetry_forward(width, height, i, j, &x, &y);
			} else {
				change_symmetry_backward(width, height, i, j, &x, &y);
			}

			temp[j * width + i] = buf[y * width + x];
		}
	}

	memcpy(buf, temp, sizeof(WORD) * width * height);
	free(temp);
}

static void centered_float(float *buf, unsigned int width, unsigned int height,
		int type) {
	unsigned int i, j;

	float *temp = malloc(width * height * sizeof(float));
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			unsigned int x = i;
			unsigned int y = j;
			if (type == FFTW_FORWARD) {
				change_symmetry_forward(width, height, i, j, &x, &y);
			} else {
				change_symmetry_backward(width, height, i, j, &x, &y);
			}

			temp[j * width + i] = buf[y * width + x];
		}
	}

	memcpy(buf, temp, sizeof(float) * width * height);
	free(temp);
}

static void normalisation_spectra_ushort(unsigned int w, unsigned int h, const float *modul, const float *phase,
		WORD *abuf, WORD *pbuf, float maxi) {
	for (size_t i = 0; i < h * w; i++) {
		pbuf[i] = roundf_to_WORD(((phase[i] + (float)G_PI) * USHRT_MAX_SINGLE / (2.f * (float)G_PI)));
		abuf[i] = roundf_to_WORD((modul[i] * USHRT_MAX_SINGLE / maxi));
	}
}

static void normalisation_spectra_float(unsigned int w, unsigned int h, const float *modul, const float *phase,
		float *abuf, float *pbuf, float maxi) {
	for (size_t i = 0; i < h * w; i++) {
		pbuf[i] = (phase[i] + (float)G_PI) / (2.f * (float)G_PI);
		abuf[i] = (modul[i] / maxi);
	}
}

static void save_dft_information_in_gfit(fits *fit) {
	int i;

	strcpy(gfit->keywords.dft.ord, fit->keywords.dft.type);
	strcpy(gfit->keywords.dft.ord, fit->keywords.dft.ord);
	for (i = 0; i < fit->naxes[2]; i++)
		gfit->keywords.dft.norm[i] = fit->keywords.dft.norm[i];

}

static void FFTD_ushort(fits *fit, fits *x, fits *y, int type_order, int layer) {
	WORD *xbuf = x->pdata[layer];
	WORD *ybuf = y->pdata[layer];
	WORD *gbuf = fit->pdata[layer];
	unsigned int width = fit->rx, height = fit->ry;
	float maxi;
	size_t i, nbdata = width * height;

	fftwf_complex *spatial_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!spatial_repr) {
		PRINT_ALLOC_ERR;
		return;
	}
	fftwf_complex *frequency_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!frequency_repr) {
		PRINT_ALLOC_ERR;
		fftwf_free(spatial_repr);
		return;
	}


	/* we run the Fourier Transform */
	fftwf_plan p = fftwf_plan_dft_2d(height, width, spatial_repr, frequency_repr,
			FFTW_FORWARD, strategy);
	/* copying image selection into the fftw data */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if(nbdata > 15000)
#endif
	for (i = 0; i < nbdata; i++) {
		spatial_repr[i] = (float) gbuf[i];
	}
	fftwf_execute(p);

	/* we compute modulus and phase */
	float *modul = malloc(nbdata * sizeof(float));
	float *phase = malloc(nbdata * sizeof(float));

	fft_to_spectra(frequency_repr, modul, phase, &maxi, nbdata);

	//We normalize the modulus and the phase
	normalisation_spectra_ushort(width, height, modul, phase, xbuf, ybuf, maxi);
	if (type_order == TYPE_CENTERED) {
		strcpy(x->keywords.dft.ord, "CENTERED");
		centered_ushort(xbuf, width, height, FFTW_FORWARD);
		centered_ushort(ybuf, width, height, FFTW_FORWARD);
	} else {
		strcpy(x->keywords.dft.ord, "REGULAR");
	}
	strcpy(y->keywords.dft.ord, x->keywords.dft.ord);
	x->keywords.dft.norm[layer] = maxi / USHRT_MAX_SINGLE;

	free(modul);
	free(phase);
	fftwf_destroy_plan(p);
	fftwf_free(spatial_repr);
	fftwf_free(frequency_repr);
}

static void FFTD_float(fits *fit, fits *x, fits *y, int type_order, int layer) {
	float *xbuf = x->fpdata[layer];
	float *ybuf = y->fpdata[layer];
	float *gbuf = fit->fpdata[layer];
	unsigned int width = fit->rx, height = fit->ry;
	float maxi;
	size_t i, nbdata = width * height;

	fftwf_complex *spatial_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!spatial_repr) {
		PRINT_ALLOC_ERR;
		return;
	}
	fftwf_complex *frequency_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!frequency_repr) {
		PRINT_ALLOC_ERR;
		fftwf_free(spatial_repr);
		return;
	}


	/* we run the Fourier Transform */
	fftwf_plan p = fftwf_plan_dft_2d(height, width, spatial_repr, frequency_repr,
			FFTW_FORWARD, strategy);
	/* copying image selection into the fftw data */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if(nbdata > 15000)
#endif
	for (i = 0; i < nbdata; i++) {
		spatial_repr[i] = (float) gbuf[i];
	}
	fftwf_execute(p);

	/* we compute modulus and phase */
	float *modul = malloc(nbdata * sizeof(float));
	float *phase = malloc(nbdata * sizeof(float));

	fft_to_spectra(frequency_repr, modul, phase, &maxi, nbdata);

	//We normalize the modulus and the phase
	normalisation_spectra_float(width, height, modul, phase, xbuf, ybuf, maxi);
	if (type_order == TYPE_CENTERED) {
		strcpy(x->keywords.dft.ord, "CENTERED");
		centered_float(xbuf, width, height, FFTW_FORWARD);
		centered_float(ybuf, width, height, FFTW_FORWARD);
	} else {
		strcpy(x->keywords.dft.ord, "REGULAR");
	}
	strcpy(y->keywords.dft.ord, x->keywords.dft.ord);
	x->keywords.dft.norm[layer] = maxi;

	free(modul);
	free(phase);
	fftwf_destroy_plan(p);
	fftwf_free(spatial_repr);
	fftwf_free(frequency_repr);
}

static void FFTD(fits *fit, fits *xfit, fits *yfit, int type_order, int layer) {
	if (fit->type == DATA_USHORT) {
		FFTD_ushort(fit, xfit, yfit, type_order, layer);
	} else if (fit->type == DATA_FLOAT) {
		FFTD_float(fit, xfit, yfit, type_order, layer);
	}
}

static void FFTI_ushort(fits *fit, fits *xfit, fits *yfit, int type_order, int layer) {
	WORD *xbuf = xfit->pdata[layer];
	WORD *ybuf = yfit->pdata[layer];
	WORD *gbuf = fit->pdata[layer];
	unsigned int width = xfit->rx;
	unsigned int height = xfit->ry;
	size_t i, nbdata = width * height;

	float *modul = calloc(1, nbdata * sizeof(float));
	float *phase = calloc(1, nbdata * sizeof(float));

	if (type_order == TYPE_CENTERED) {
		centered_ushort(xbuf, width, height, FFTW_BACKWARD);
		centered_ushort(ybuf, width, height, FFTW_BACKWARD);
	}

	for (i = 0; i < height * width; i++) {
		modul[i] = (float) xbuf[i] * (xfit->keywords.dft.norm[layer]);
		phase[i] = (float) ybuf[i] * (2.f * (float)G_PI / USHRT_MAX_SINGLE);
		phase[i] -= (float)G_PI;
	}

	fftwf_complex* spatial_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!spatial_repr) {
		PRINT_ALLOC_ERR;
		free(modul);
		free(phase);
		return;
	}

	fftwf_complex* frequency_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!frequency_repr) {
		PRINT_ALLOC_ERR;
		fftwf_free(spatial_repr);
		free(modul);
		free(phase);
		return;
	}


	fftwf_plan p = fftwf_plan_dft_2d(height, width, frequency_repr, spatial_repr,
			FFTW_BACKWARD, strategy);
	fft_to_freq(frequency_repr, modul, phase, nbdata);
	fftwf_execute(p);

	for (i = 0; i < nbdata; i++) {
		float pxl = crealf(spatial_repr[i]) / nbdata;
		gbuf[i] = roundf_to_WORD(pxl);
	}
	gui_iface.delete_selection();
	invalidate_stats_from_fit(fit);

	free(modul);
	free(phase);
	fftwf_destroy_plan(p);
	fftwf_free(spatial_repr);
	fftwf_free(frequency_repr);
}

static void FFTI_float(fits *fit, fits *xfit, fits *yfit, int type_order, int layer) {
	float *xbuf = xfit->fpdata[layer];
	float *ybuf = yfit->fpdata[layer];
	float *gbuf = fit->fpdata[layer];
	unsigned int width = xfit->rx;
	unsigned int height = xfit->ry;
	size_t i, nbdata = width * height;

	float *modul = calloc(1, nbdata * sizeof(float));
	float *phase = calloc(1, nbdata * sizeof(float));

	if (type_order == TYPE_CENTERED) {
		centered_float(xbuf, width, height, FFTW_BACKWARD);
		centered_float(ybuf, width, height, FFTW_BACKWARD);
	}

	for (i = 0; i < height * width; i++) {
		modul[i] = xbuf[i] * (xfit->keywords.dft.norm[layer]);
		phase[i] = ybuf[i] * (2.f * (float)G_PI);
		phase[i] -= (float)G_PI;
	}

	fftwf_complex* spatial_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!spatial_repr) {
		PRINT_ALLOC_ERR;
		free(modul);
		free(phase);
		return;
	}

	fftwf_complex* frequency_repr = fftwf_malloc(sizeof(fftwf_complex) * nbdata);
	if (!frequency_repr) {
		PRINT_ALLOC_ERR;
		fftwf_free(spatial_repr);
		free(modul);
		free(phase);
		return;
	}


	fftwf_plan p = fftwf_plan_dft_2d(height, width, frequency_repr, spatial_repr,
			FFTW_BACKWARD, strategy);
	fft_to_freq(frequency_repr, modul, phase, nbdata);
	fftwf_execute(p);

	for (i = 0; i < nbdata; i++) {
		float pxl = crealf(spatial_repr[i]) / nbdata;
		gbuf[i] = pxl;
	}
	gui_iface.delete_selection();
	invalidate_stats_from_fit(fit);

	free(modul);
	free(phase);
	fftwf_destroy_plan(p);
	fftwf_free(spatial_repr);
	fftwf_free(frequency_repr);
}

static void FFTI(fits *fit, fits *xfit, fits *yfit, int type_order, int layer) {
	if (fit->type == DATA_USHORT) {
		FFTI_ushort(fit, xfit, yfit, type_order, layer);
	} else if (fit->type == DATA_FLOAT) {
		FFTI_float(fit, xfit, yfit, type_order, layer);
	}
}

/* end_fourier_transform and fft_idle moved to src/gui/fft.c */

/* Destructor for fft_data - used by free_generic_img_args */
void free_fft_data(void *p) {
	struct fft_data *args = (struct fft_data *)p;
	if (!args)
		return;
	g_free(args->type);
	g_free(args->modulus);
	g_free(args->phase);
	free(args);
}

/* Core FFT computation: reads from fit, writes result back into fit.
 * Returns 0 on success, non-zero on failure. */
int fft_compute_core(struct fft_data *args, fits *fit) {
	unsigned int width = fit->rx;
	unsigned int height = fit->ry;
	int chan;
	int retval = 0;
	fits *tmp = NULL, *tmp1 = NULL, *tmp2 = NULL;
	data_type type = fit->type;

#ifdef HAVE_FFTW3F_MULTITHREAD
	int n = com.pref.fftw_conf.multithreaded ? com.max_thread : 1;
	fftwf_plan_with_nthreads(n);
	fprintf(stdout, "fftwf initialized with %d threads\n", n);
#endif

	strategy = FFTW_ESTIMATE;
	set_wisdom_file();
	if (fftwf_import_wisdom_from_filename(com.pref.fftw_conf.wisdom_file) == 1) {
		siril_log_message(_("Siril FFT wisdom imported successfully...\n"));
	} else if (fftwf_import_system_wisdom() == 1) {
		siril_log_message(_("System FFT wisdom imported successfully...\n"));
	} else {
		siril_log_message(_("No FFT wisdom found to import...\n"));
	}

	/* "Fourier Transform: processing..." is logged by the framework before the hook runs */

	switch (args->type[3]) {
	default:
	case 'd':
	case 'D':
		if (new_fit_image(&tmp1, width, height, fit->naxes[2], type) ||
				new_fit_image(&tmp2, width, height, fit->naxes[2], type)) {
			retval = 1;
			goto end;
		}

		for (chan = 0; chan < fit->naxes[2]; chan++)
			FFTD(fit, tmp1, tmp2, args->type_order, chan);
		strcpy(tmp1->keywords.dft.type, "SPECTRUM");
		if (savefits(args->modulus, tmp1)) {
			retval = 1;
			goto end;
		}
		strcpy(tmp2->keywords.dft.type, "PHASE");
		if (savefits(args->phase, tmp2)) {
			retval = 1;
			goto end;
		}

		/* Copy modulus into fit for display */
		if (copyfits(tmp1, fit, CP_ALLOC | CP_FORMAT | CP_COPYA, -1)) {
			retval = 1;
			goto end;
		}

		/* copy the header informations */
		save_dft_information_in_gfit(tmp1);
		break;
	case 'i':
	case 'I':
		tmp = calloc(1, sizeof(fits));
		if (!tmp || readfits(args->modulus, tmp, NULL, FALSE)) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto end;
		}
		tmp1 = calloc(1, sizeof(fits));
		if (!tmp1 || readfits(args->phase, tmp1, NULL, FALSE)) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto end;
		}
		if ((tmp->rx != tmp1->rx) || (tmp->ry != tmp1->ry) ||
				(tmp->naxes[2] != tmp1->naxes[2]) || (tmp->bitpix != tmp1->bitpix)) {
			retval = 1;
			siril_log_error(_("Images must have same dimensions.\n"));
			goto end;
		}
		if (tmp->keywords.dft.ord[0] == 'C')
			args->type_order = TYPE_CENTERED;
		else if (tmp->keywords.dft.ord[0] == 'R')
			args->type_order = TYPE_REGULAR;
		else {
			retval = 1;
			siril_log_message(_("There is something wrong in your files\n"));
			goto end;
		}
		new_fit_image(&tmp2, width, height, tmp->naxes[2], type);
		for (chan = 0; chan < fit->naxes[2]; chan++)
			FFTI(tmp2, tmp, tmp1, args->type_order, chan);
		/* Copy inverse transform result into fit for display */
		if (copyfits(tmp2, fit, CP_ALLOC | CP_FORMAT | CP_COPYA, -1)) {
			retval = 1;
			goto end;
		}
	}

end:
	if (fftwf_export_wisdom_to_filename(com.pref.fftw_conf.wisdom_file) == 1) {
		siril_log_message(_("Siril FFT wisdom updated successfully...\n"));
	} else {
		siril_log_warning(_("Siril FFT wisdom update failed...\n"));
	}

	invalidate_stats_from_fit(fit);
	if (tmp)  { clearfits(tmp);  free(tmp);  }
	if (tmp1) { clearfits(tmp1); free(tmp1); }
	if (tmp2) { clearfits(tmp2); free(tmp2); }

	return retval;
}

/* generic_image_worker image hook */
int fft_image_hook(struct generic_img_args *gargs, fits *fit, int threads) {
	struct fft_data *args = (struct fft_data *)gargs->user;
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	int retval = fft_compute_core(args, fit);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	return retval;
}

gchar *fft_log_hook(gpointer p, log_hook_detail detail) {
	struct fft_data *args = (struct fft_data *)p;
	gboolean is_inverse = (args->type[3] == 'i' || args->type[3] == 'I');
	if (is_inverse)
		return g_strdup(_("Fourier Transform (inverse)"));
	return g_strdup(_("Fourier Transform (direct)"));
}

/* GUI idle functions and callbacks moved to src/gui/fft.c */
