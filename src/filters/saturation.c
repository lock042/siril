/*
 * Refactored saturation using generic_image_worker
 */

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "algos/colors.h"
#include "algos/statistics.h"

#include "saturation.h"

/* Helper to map hue types to degree ranges */
void satu_set_hues_from_types(saturation_params *args, int type) {
	switch (type) {
		case 0:		// Pink-Red to Red-Orange
			args->h_min = 346.0;
			args->h_max = 20.0;
			break;
		case 1:		// Orange-Brown to Yellow
			args->h_min = 21.0;
			args->h_max = 60.0;
			break;
		case 2:		// Yellow-Green to Green-Cyan
			args->h_min = 61.0;
			args->h_max = 200.0;
			break;
		case 3:		// Cyan
			args->h_min = 170.0;
			args->h_max = 200.0;
			break;
		case 4:		// Cyan-Blue to Blue-Magenta
			args->h_min = 201.0;
			args->h_max = 280.0;
			break;
		case 5:		// Magenta to Pink
			args->h_min = 281.0;
			args->h_max = 345.0;
			break;
		default:
		case 6:		// Global
			args->h_min = 0.0;
			args->h_max = 360.0;
	}
}

/* Core Algorithm: USHORT */
static int enhance_saturation_ushort(fits *fit, saturation_params *params) {
	double bg = 0.0;
	WORD *in[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	WORD *out[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };

	double h_min = params->h_min / 360.0;
	double h_max = params->h_max / 360.0;

	if (params->background_factor > 0.00) {
		imstats *stat = statistics(NULL, -1, fit, GLAYER, NULL, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg = (stat->median + stat->sigma) * params->background_factor;
		bg /= stat->normValue;
		free_stats(stat);
	}
	siril_debug_print("threshold for saturation: %f\n", bg);

	gboolean loop_range = h_min > h_max;
	double s_mult = 1.0 + params->coeff;
	double norm = fit->bitpix == BYTE_IMG ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
	double invnorm = 1.0 / norm;
	size_t i, n = fit->naxes[0] * fit->naxes[1];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < n; i++) {
		double h, s, l;
		double r = in[RLAYER][i] * invnorm;
		double g = in[GLAYER][i] * invnorm;
		double b = in[BLAYER][i] * invnorm;
		rgb_to_hsl(r, g, b, &h, &s, &l);
		if (l > bg) {
			if (loop_range) {
				if (h >= h_min || h <= h_max)
					s *= s_mult;
			} else {
				if (h >= h_min && h <= h_max)
					s *= s_mult;
			}
			s = (s < 0.0) ? 0.0 : (s > 1.0) ? 1.0 : s;
			hsl_to_rgb(h, s, l, &r, &g, &b);
		}
		out[RLAYER][i] = round_to_WORD(r * norm);
		out[GLAYER][i] = round_to_WORD(g * norm);
		out[BLAYER][i] = round_to_WORD(b * norm);
	}
	return 0;
}

/* Core Algorithm: FLOAT */
static int enhance_saturation_float(fits *fit, saturation_params *params) {
	float bg = 0.0f;
	float *in[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	float *out[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	float h_min = (float)params->h_min / 60.0f;
	float h_max = (float)params->h_max / 60.0f;

	if (params->background_factor > 0.00) {
		imstats *stat = statistics(NULL, -1, fit, GLAYER, NULL, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg = (float)((stat->median + stat->sigma) * params->background_factor);
		bg /= (float)stat->normValue;
		free_stats(stat);
	}
	siril_debug_print("threshold for saturation: %f\n", bg);

	gboolean loop_range = h_min > h_max;
	float s_mult = 1.f + (float)params->coeff;

	size_t i, n = fit->naxes[0] * fit->naxes[1];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(dynamic, fit->rx * 16)
#endif
	for (i = 0; i < n; i++) {
		float h, s, l;
		float r = in[RLAYER][i];
		float g = in[GLAYER][i];
		float b = in[BLAYER][i];
		rgb_to_hsl_float_sat(r, g, b, bg, &h, &s, &l);
		if (l > bg) {
			if (loop_range) {
				if (h >= h_min || h <= h_max)
					s *= s_mult;
			} else {
				if (h >= h_min && h <= h_max)
					s *= s_mult;
			}
			s = (s < 0.f) ? 0.f : (s > 1.f) ? 1.f : s;
			hsl_to_rgb_float_sat(h, s, l, &r, &g, &b);
		}
		out[RLAYER][i] = r;
		out[GLAYER][i] = g;
		out[BLAYER][i] = b;
	}
	return 0;
}

/* The Generic Processing Hook */
int saturation_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	saturation_params *params = (saturation_params *)args->user;
	if (!params)
		return 1;

	if (fit->type == DATA_USHORT) {
		return enhance_saturation_ushort(fit, params);
	} else if (fit->type == DATA_FLOAT) {
		return enhance_saturation_float(fit, params);
	}
	return 1;
}

gchar* satu_log_hook(gpointer p, log_hook_detail detail) {
	saturation_params *params = (saturation_params *) p;
	return g_strdup_printf(_("Color saturation %d%%, threshold %.2f"),
		round_to_int(params->coeff * 100.0), params->background_factor);
}

/* GUI callbacks moved to src/gui/saturation.c */
