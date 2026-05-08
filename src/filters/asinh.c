/*
 * Refactored asinh stretch using generic_image_worker — processing only.
 * GUI callbacks and dialog logic live in src/gui/asinh.c.
 */

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/processing.h"
#include "algos/statistics.h"

#include "asinh.h"

/* The actual asinh processing hook */
int asinh_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	asinh_params *params = (asinh_params *)args->user;
	if (!params)
		return 1;

	return asinhlut(fit, params->beta, params->offset, params->human_luminance, params->clip_mode);
}

gchar *asinh_log_hook(gpointer p, log_hook_detail detail) {
	asinh_params *params = (asinh_params*) p;
	gchar *message = g_strdup_printf(_("Asinh Transformation: (stretch=%6.1f, bp=%7.5f)"),
				params->beta, params->offset);
	return message;
}

/* Keep the existing asinhlut functions with clip_mode as parameter */
int asinhlut_ushort(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clip_mode) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	const gboolean do_channel[3] = { TRUE, TRUE, TRUE };
	const float m_CB = 1.f;
	float norm = get_normalized_value(fit);
	float invnorm = 1.0f / norm;
	float asinh_beta = asinh(beta);
	const float third = 1.f / 3.f;
	float factor_red = human_luminance ? 0.2126f : third;
	float factor_green = human_luminance ? 0.7152f : third;
	float factor_blue = human_luminance ? 0.0722f : third;

	float inv_1moffset = 1.f / (1.f - offset);
	size_t n = fit->naxes[0] * fit->naxes[1];
	float globalmax = -FLT_MAX;
	if (fit->naxes[2] > 1) {
		switch (clip_mode) {
			case CLIP:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0; i < n; i++) {
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (buf[chan][i] * invnorm - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					for (int chan = 0; chan < 3; chan++)
						buf[chan][i] = roundf_to_WORD(norm * fmaxf(0.f, fminf(prime[chan] * k, 1.f)));
				}
				break;
			case RESCALE:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0; i < n; i++) {
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (buf[chan][i] * invnorm - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					float sf[3] = { fmaxf(0.f, prime[0] * k), fmaxf(0.f, prime[1] * k), fmaxf(0.f, prime[2] * k) };
					float maxval = fmaxf(fmaxf(sf[0], sf[1]), sf[2]);
					if (maxval > 1.f) {
						float invmaxval = 1.f / maxval;
						for (int chan = 0; chan < 3; chan++)
							buf[chan][i] = roundf_to_WORD(norm * sf[chan] * invmaxval);
					} else {
						for (int chan = 0; chan < 3; chan++)
							buf[chan][i] = roundf_to_WORD(norm * sf[chan]);
					}
				}
				break;
			case RGBBLEND:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0; i < n; i++) {
					blend_data data = { .sf = { 0.f }, .tf = { 0.f }, .do_channel = do_channel };
					float val[3] = { buf[0][i] * invnorm, buf[1][i] * invnorm, buf[2][i] * invnorm };
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (val[chan] - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					for (int chan = 0; chan < 3; chan++)
						data.sf[chan] = fminf(1.f, fmaxf(0.f, prime[chan] * k));
					for (int chan = 0; chan < 3; chan++)
						data.tf[chan] = (prime[chan] == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * prime[chan]) / (prime[chan] * asinh_beta);
					float out[3] = { val[0], val[1], val[2] };
					rgbblend(&data, &out[0], &out[1], &out[2], m_CB);
					for (int chan = 0; chan < 3; chan++)
						buf[chan][i] = roundf_to_WORD(norm * out[chan]);
				}
				break;
			case RESCALEGLOBAL:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(max:globalmax)
#endif
				for (size_t i = 0; i < n; i++) {
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (buf[chan][i] * invnorm - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					float maxval = fmaxf(fmaxf(fminf(1.f, prime[0] * k),
					                           fminf(1.f, prime[1] * k)),
					                           fminf(1.f, prime[2] * k));
					if (maxval > globalmax)
						globalmax = maxval;
				}
				{
					float inv_globalmax = 1.f / fmaxf(globalmax, FLT_MIN);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
					for (size_t i = 0; i < n; i++) {
						float prime[3];
						for (int chan = 0; chan < 3; chan++)
							prime[chan] = fmaxf(0.f, (buf[chan][i] * invnorm - offset) * inv_1moffset);
						float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
						float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
						for (int chan = 0; chan < 3; chan++)
							buf[chan][i] = roundf_to_WORD(norm * fmaxf(0.f, fminf(1.f, prime[chan] * k) * inv_globalmax));
					}
				}
				break;
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			float x = buf[RLAYER][i] * invnorm;
			float xprime = fmaxf(0.f, (x - offset) * inv_1moffset);
			float k = (xprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = roundf_to_WORD(norm * fminf(1.f, fmaxf(0.f, xprime * k)));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int asinhlut_float(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clip_mode) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	const gboolean do_channel[3] = { TRUE, TRUE, TRUE };
	float m_CB = 1.f;
	float asinh_beta = asinhf(beta);
	float factor_red = human_luminance ? 0.2126f : 0.3333f;
	float factor_green = human_luminance ? 0.7152f : 0.3333f;
	float factor_blue = human_luminance ? 0.0722f : 0.3333f;
	float inv_1moffset = 1.f / (1.f - offset);

	size_t n = fit->naxes[0] * fit->naxes[1];
	float globalmax = -FLT_MAX;
	if (fit->naxes[2] > 1) {
		switch (clip_mode) {
			case CLIP:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0; i < n; i++) {
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (buf[chan][i] - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					for (int chan = 0; chan < 3; chan++)
						buf[chan][i] = fmaxf(0.f, fminf(prime[chan] * k, 1.f));
				}
				break;
			case RESCALE:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0; i < n; i++) {
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (buf[chan][i] - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					float sf[3] = { fmaxf(0.f, prime[0] * k), fmaxf(0.f, prime[1] * k), fmaxf(0.f, prime[2] * k) };
					float maxval = fmaxf(fmaxf(sf[0], sf[1]), sf[2]);
					if (maxval > 1.f) {
						float invmaxval = 1.f / maxval;
						for (int chan = 0; chan < 3; chan++)
							buf[chan][i] = sf[chan] * invmaxval;
					} else {
						for (int chan = 0; chan < 3; chan++)
							buf[chan][i] = sf[chan];
					}
				}
				break;
			case RGBBLEND:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0; i < n; i++) {
					blend_data data = { .sf = { 0.f }, .tf = { 0.f }, .do_channel = do_channel };
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (buf[chan][i] - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					for (int chan = 0; chan < 3; chan++)
						data.sf[chan] = fminf(1.f, fmaxf(0.f, prime[chan] * k));
					for (int chan = 0; chan < 3; chan++)
						data.tf[chan] = (prime[chan] == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * prime[chan]) / (prime[chan] * asinh_beta);
					rgbblend(&data, &buf[0][i], &buf[1][i], &buf[2][i], m_CB);
				}
				break;
			case RESCALEGLOBAL:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(max:globalmax)
#endif
				for (size_t i = 0; i < n; i++) {
					float prime[3];
					for (int chan = 0; chan < 3; chan++)
						prime[chan] = fmaxf(0.f, (buf[chan][i] - offset) * inv_1moffset);
					float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
					float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
					float maxval = fmaxf(fmaxf(fminf(1.f, prime[0] * k),
					                           fminf(1.f, prime[1] * k)),
					                           fminf(1.f, prime[2] * k));
					if (maxval > globalmax)
						globalmax = maxval;
				}
				{
					float inv_globalmax = 1.f / fmaxf(globalmax, FLT_MIN);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
					for (size_t i = 0; i < n; i++) {
						float prime[3];
						for (int chan = 0; chan < 3; chan++)
							prime[chan] = fmaxf(0.f, (buf[chan][i] - offset) * inv_1moffset);
						float x = factor_red * prime[0] + factor_green * prime[1] + factor_blue * prime[2];
						float k = (x == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * x) / (x * asinh_beta);
						for (int chan = 0; chan < 3; chan++)
							buf[chan][i] = fmaxf(0.f, fminf(1.f, prime[chan] * k) * inv_globalmax);
					}
				}
				break;
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			float xprime = fmaxf(0.f, (buf[RLAYER][i] - offset) * inv_1moffset);
			float k = (xprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = fminf(1.f, fmaxf(0.f, xprime * k));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int asinhlut(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clip_mode) {
	if (fit->type == DATA_USHORT)
		return asinhlut_ushort(fit, beta, offset, human_luminance, clip_mode);
	if (fit->type == DATA_FLOAT)
		return asinhlut_float(fit, beta, offset, human_luminance, clip_mode);
	return 1;
}

int command_asinh(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clipmode) {
	return asinhlut(fit, beta, offset, human_luminance, clipmode);
}
