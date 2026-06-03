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

/*
 * Expression-evaluation engine for PixelMath.  This file is GTK-free;
 * it contains only the computation worker, variable-image management,
 * and the post-operation idle that cleans up via gui_iface.
 *
 * The widget management (dialog callbacks, tree-view helpers, file
 * chooser, status labels) lives in gui/pixelmath.c.
 */

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/conversion.h"
#include "algos/demosaicing.h"

#include "tinyexpr.h"
#include "pixel_math_runner.h"

#define T_CURRENT "gfit"

static const gchar *variables[] = {
		"I1",
		"I2",
		"I3",
		"I4",
		"I5",
		"I6",
		"I7",
		"I8",
		"I9",
		"I10"
};

#define MAX_IMAGES G_N_ELEMENTS(variables)

fits var_fit[MAX_IMAGES] = { 0 };  /* exported: gui/pixelmath.c loads images here */
static gboolean var_fit_mask[MAX_IMAGES] = { 0 };

/* ── Post-operation cleanup idle ────────────────────────────────────────── */

static gboolean end_pixel_math_operation(gpointer p) {
	struct pixel_math_data *args = (struct pixel_math_data *)p;
	stop_processing_thread();

	if (!args->ret) {
		/* gfit was already written by the worker; just do GTK-side work */
		if (sequence_is_loaded())
			close_sequence(FALSE);
		gui_iface.on_image_loaded();
	}

	gui_iface.set_busy(FALSE);
	if (args->from_ui)
		gui_iface.update_pixel_math_status(args->ret);

	free(args->fit);
	free(args);
	return FALSE;
}

/* ── Processing helpers ─────────────────────────────────────────────────── */

static void update_metadata(fits *fit, gboolean do_sum) {
	fits **f = malloc((MAX_IMAGES + 1) * sizeof(fits *));
	int j = 0;
	for (int i = 0; i < MAX_IMAGES ; i++)
		if (var_fit[i].rx > 0 && var_fit_mask[i])
			f[j++] = &var_fit[i];
	f[j] = NULL;

	if (!f[0] && single_image_is_loaded() )
		/* if no fit used (only constants), copy metadata from gfit */
		copy_fits_metadata(gfit, fit);
	else
		merge_fits_headers_to_result2(fit, f, do_sum);
	update_fits_header(fit);
	free(f);
}

static gchar* parse_image_functions(gpointer p, int idx, int c) {
	struct pixel_math_data *args = (struct pixel_math_data*) p;
	gchar *expression;
	gchar **image = args->varname;
	int nb_images = args->nb_rows;
	switch (idx) {
	case 1:
		expression = args->expression1;
		break;
	case 2:
		expression = args->expression2;
		break;
	case 3:
		expression = args->expression3;
		break;
	default:
		return NULL;
	}
	if (!expression)
		return expression;

	gchar *result = g_strdup(expression);
	GRegex *regex = g_regex_new("(\\w+)\\((\\w+)\\)", 0, 0, NULL);

	gboolean replaced = TRUE;
	while (replaced) {
		replaced = FALSE;
		GMatchInfo *match_info;
		g_regex_match(regex, result, 0, &match_info);

		if (g_match_info_matches(match_info)) {
			gchar *function = g_match_info_fetch(match_info, 1);
			gchar *param = g_match_info_fetch(match_info, 2);
			gchar *full_match = g_match_info_fetch(match_info, 0);
			double median = 0.0, mean = 0.0, min = 0.0, max = 0.0, noise = 0.0, adev = 0.0, bwmv = 0.0, mad = 0.0, sdev = 0.0;
			double w = 0.0, h = 0.0;
			imstats *stats = NULL;

			if (g_strcmp0(param, T_CURRENT) == 0) {
				stats = statistics(NULL, -1, gfit, c, NULL, STATS_MAIN, MULTI_THREADED);
				if (!stats) {
					g_free(full_match);
					g_free(function);
					g_free(param);
					g_match_info_free(match_info);
					g_regex_unref(regex);
					return result;
				}
				median = stats->median;
				mean = stats->mean;
				min = stats->min;
				max = stats->max;
				noise = stats->bgnoise;
				adev = stats->avgDev;
				bwmv = stats->sqrtbwmv * stats->sqrtbwmv;
				mad = stats->mad;
				sdev = stats->sigma;
				w = (double) gfit->rx;
				h = (double) gfit->ry;
				free_stats(stats);
			} else {
				for (int j = 0; j < nb_images; j++) {
					if (g_strcmp0(param, image[j]) == 0) {
						stats = statistics(NULL, -1, &var_fit[j], c, NULL, STATS_MAIN, MULTI_THREADED);
						if (!stats) {
							g_free(full_match);
							g_free(function);
							g_free(param);
							g_match_info_free(match_info);
							g_regex_unref(regex);
							return result;
						}
						median = stats->median;
						mean = stats->mean;
						min = stats->min;
						max = stats->max;
						noise = stats->bgnoise;
						adev = stats->avgDev;
						bwmv = stats->sqrtbwmv * stats->sqrtbwmv;
						mad = stats->mad;
						sdev = stats->sigma;
						w = (double) var_fit[j].rx;
						h = (double) var_fit[j].ry;
						free_stats(stats);
						break;
					}
				}
			}

			gchar *replace = NULL;
			if (!g_strcmp0(function, "mean")) {
				replace = g_strdup_printf("%g", mean);
			} else if (!g_strcmp0(function, "med") || !g_strcmp0(function, "median")) {
				replace = g_strdup_printf("%g", median);
			} else if (!g_strcmp0(function, "min")) {
				replace = g_strdup_printf("%g", min);
			} else if (!g_strcmp0(function, "max")) {
				replace = g_strdup_printf("%g", max);
			} else if (!g_strcmp0(function, "noise")) {
				replace = g_strdup_printf("%g", noise);
			} else if (!g_strcmp0(function, "adev")) {
				replace = g_strdup_printf("%g", adev);
			} else if (!g_strcmp0(function, "bwmv")) {
				replace = g_strdup_printf("%g", bwmv);
			} else if (!g_strcmp0(function, "mad") || !g_strcmp0(function, "mdev")) {
				replace = g_strdup_printf("%g", mad);
			} else if (!g_strcmp0(function, "sdev")) {
				replace = g_strdup_printf("%g", sdev);
			} else if (!g_strcmp0(function, "width") || !g_strcmp0(function, "w")) {
				replace = g_strdup_printf("%g", w);
			} else if (!g_strcmp0(function, "height") || !g_strcmp0(function, "h")) {
				replace = g_strdup_printf("%g", h);
			}

			if (replace) {
				gchar *temp = result;
				gchar **split = g_strsplit(result, full_match, 2);
				if (split[0] && split[1]) {
					result = g_strconcat(split[0], replace, split[1], NULL);
				} else {
					result = g_strdup(result);
				}
				g_strfreev(split);
				g_free(temp);
				g_free(replace);
				replaced = TRUE;
				siril_log_debug("Expression%d: %s\n", c, result);
			}

			g_free(full_match);
			g_free(function);
			g_free(param);
		}

		g_match_info_free(match_info);
	}

	g_regex_unref(regex);

	for (int j = 0; j < nb_images; j++) {
		const gchar *test = g_strrstr(result, image[j]);
		if (test) {
			var_fit_mask[j] = TRUE;
			siril_log_debug("found image name %s in the expression %s\n", image[j], result);
		}
	}

	return result;
}

/* ── Main worker ────────────────────────────────────────────────────────── */

gpointer apply_pixel_math_operation(gpointer p) {
	struct pixel_math_data *args = (struct pixel_math_data *)p;

	te_expr *n1 = NULL, *n2 = NULL, *n3 = NULL;
	fits *fit = args->fit;
	int nb_rows = args->nb_rows;
	gboolean failed = FALSE;
	args->ret = 0;
	float maximum = -FLT_MAX;
	float minimum = +FLT_MAX;
	long width, height, nchan;

	/* Writer lock for the duration of OMP gfit reads (when args->has_gfit).
	 * Headless mode is excluded: the Python module does not run headless, and
	 * memcpy(gfit, args->fit, sizeof(fits)) in the headless path below would
	 * overwrite gfit->rwlock while it is held. */
	gboolean rwlocked = FALSE;
	if (!com.headless) {
		g_rw_lock_writer_lock(&gfit->rwlock);
		rwlocked = TRUE;
	}

	if (args->single_rgb && args->fit->naxes[2] > 1) {
		args->expression2 = g_strdup(args->expression1);
		args->expression3 = g_strdup(args->expression1);

		args->expression1 = parse_image_functions(args, 1, RLAYER);
		args->expression2 = parse_image_functions(args, 2, GLAYER);
		args->expression3 = parse_image_functions(args, 3, BLAYER);
	} else {
		args->expression1 = parse_image_functions(args, 1, RLAYER);
		args->expression2 = parse_image_functions(args, 2, RLAYER);
		args->expression3 = parse_image_functions(args, 3, RLAYER);
	}
	siril_log_message(_("Pixelmath expression 1: %s\n"), args->expression1);
	if (args->expression2 && strncmp(args->expression1, args->expression2, strlen(args->expression1)))
		siril_log_message(_("Pixelmath expression 2: %s\n"), args->expression2);
	if (args->expression3 && strncmp(args->expression1, args->expression3, strlen(args->expression1)))
		siril_log_message(_("Pixelmath expression 3: %s\n"), args->expression3);
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) firstprivate(n1,n2,n3)
#endif
	{
		int k = 0;
		if (args->has_gfit) k = 1;
		te_variable *vars = malloc((nb_rows + k) * sizeof(te_variable));
		double *x = malloc((nb_rows + k) * sizeof(double));
		if (!vars || !x) {
			failed = TRUE;
			goto failure;
		}
		for (int i = 0; i < nb_rows; i++) {
			vars[i].name = args->varname[i];
			vars[i].address = &x[i];
			vars[i].context = NULL;
			vars[i].type = 0;
		}
		if (args->has_gfit) {
			vars[nb_rows].name = g_strdup(T_CURRENT);
			vars[nb_rows].address = &x[nb_rows];
			vars[nb_rows].context = NULL;
			vars[nb_rows].type = 0;
		}
		int err = 0;
		n1 = te_compile(args->expression1, vars, nb_rows + k, &err);
		if (!n1) {
#ifdef _OPENMP
			if (omp_get_thread_num() == 0)
#endif
				siril_log_error(_("Error in pixel math expression '%s' at character %d\n"), args->expression1, err);
			failed = TRUE;
			goto failure;
		}

		if (args->expression2) {
			n2 = te_compile(args->expression2, vars, nb_rows + k, &err);
			if (!n2) {
#ifdef _OPENMP
				if (omp_get_thread_num() == 0)
#endif
					siril_log_error(_("Error in pixel math expression '%s' at character %d\n"), args->expression2, err);
				failed = TRUE;
				goto failure;
			}

			n3 = te_compile(args->expression3, vars, nb_rows + k, &err);
			if (!n3) {
#ifdef _OPENMP
				if (omp_get_thread_num() == 0)
#endif
					siril_log_error(_("Error in pixel math expression '%s' at character %d\n"), args->expression3, err);
				failed = TRUE;
				goto failure;
			}
		}
		if (args->has_gfit && nb_rows == 0) {
			width = gfit->naxes[0];
			height = gfit->naxes[1];
			nchan = gfit->naxes[2];
		} else {
			width = var_fit[0].naxes[0];
			height = var_fit[0].naxes[1];
			nchan = var_fit[0].naxes[2];
		}
		if (com.pref.force_16bit) {
#ifdef _OPENMP
#pragma omp for schedule(static) reduction(max:maximum) reduction(min:minimum)
#endif
			for (size_t px = 0; px < width * height * nchan; px++) {
				for (int i = 0; i < nb_rows; i++) {
					x[i] = (double) var_fit[i].fdata[px];
				}
				if (args->has_gfit) {
					if (gfit->type == DATA_USHORT) {
						x[nb_rows] = (double) gfit->data[px] / USHRT_MAX_DOUBLE;
					} else {
						x[nb_rows] = (double) gfit->fdata[px];
					}
				}

				if (!args->single_rgb) {
					fit->pdata[RLAYER][px] =  roundf_to_WORD((float) te_eval(n1) * USHRT_MAX_SINGLE);
					fit->pdata[GLAYER][px] =  roundf_to_WORD((float) te_eval(n2) * USHRT_MAX_SINGLE);
					fit->pdata[BLAYER][px] =  roundf_to_WORD((float) te_eval(n3) * USHRT_MAX_SINGLE);

					maximum = max(maximum, max(fit->pdata[RLAYER][px], max(fit->pdata[GLAYER][px], fit->pdata[BLAYER][px])));
					minimum = min(minimum, min(fit->pdata[RLAYER][px], min(fit->pdata[GLAYER][px], fit->pdata[BLAYER][px])));
				} else {
					if (px < (var_fit[0].naxes[0] * var_fit[0].naxes[1])) {
						fit->data[px] =  roundf_to_WORD((float) te_eval(n1) * USHRT_MAX_SINGLE);
					} else if (px < 2 * (var_fit[0].naxes[0] * var_fit[0].naxes[1])) {
						fit->data[px] =  roundf_to_WORD((float) te_eval(n2) * USHRT_MAX_SINGLE);
					} else {
						fit->data[px] =  roundf_to_WORD((float) te_eval(n3) * USHRT_MAX_SINGLE);
					}

					maximum = max(maximum, fit->data[px]);
					minimum = min(minimum, fit->data[px]);
				}
			}

		} else {
#ifdef _OPENMP
#pragma omp for schedule(static) reduction(max:maximum) reduction(min:minimum)
#endif
			for (size_t px = 0; px < width * height * nchan; px++) {
				for (int i = 0; i < nb_rows; i++) {
					x[i] = var_fit[i].fdata[px];
				}
				if (args->has_gfit) {
					if (gfit->type == DATA_USHORT) {
						x[nb_rows] = gfit->data[px] / USHRT_MAX_DOUBLE;
					} else {
						x[nb_rows] = gfit->fdata[px];
					}
				}

				if (!args->single_rgb) {
					fit->fpdata[RLAYER][px] = (float) te_eval(n1);
					fit->fpdata[GLAYER][px] = (float) te_eval(n2);
					fit->fpdata[BLAYER][px] = (float) te_eval(n3);

					maximum = max(maximum, max(fit->fpdata[RLAYER][px], max(fit->fpdata[GLAYER][px], fit->fpdata[BLAYER][px])));
					minimum = min(minimum, min(fit->fpdata[RLAYER][px], min(fit->fpdata[GLAYER][px], fit->fpdata[BLAYER][px])));
				} else {
					if (px < (width * height)) {
						fit->fdata[px] = (float) te_eval(n1);
					} else if (px < 2 * (width * height)) {
						fit->fdata[px] = (float) te_eval(n2);
					} else {
						fit->fdata[px] = (float) te_eval(n3);
					}

					maximum = max(maximum, fit->fdata[px]);
					minimum = min(minimum, fit->fdata[px]);
				}
			}
		}

failure:
		te_free(n1);
		if (args->expression2) {
			te_free(n2);
			te_free(n3);
		}
		free(vars);
		free(x);
	} // end of parallel block

	if (args->rescale) {
		if (com.pref.force_16bit) {
			args->min *= USHRT_MAX_SINGLE;
			args->max *= USHRT_MAX_SINGLE;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (int i = 0; i < fit->rx * fit->ry * fit->naxes[2]; i++) {
				fit->data[i] = roundf_to_WORD((float)(args->max - args->min) * (float)(fit->data[i] - minimum) / (float)(maximum - minimum) + args->min);
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (int i = 0; i < fit->rx * fit->ry * fit->naxes[2]; i++) {
				fit->fdata[i] = (((args->max - args->min) * (fit->fdata[i] - minimum)) / (maximum - minimum)) + args->min;
			}
		}
	}

	if (failed)
		args->ret = 1;
	else update_metadata(args->fit, args->do_sum);

	/* free memory */
	g_free(args->expression1);
	if (args->expression2) {
		g_free(args->expression2);
		g_free(args->expression3);
	}
	for (int i = 0; i < args->nb_rows; i++)
		g_free(args->varname[i]);
	free(args->varname);
	free_pm_var(args->nb_rows);

	/* manage result and display */
	if (com.headless) {
		if (!failed) {
			clearfits(gfit);
			memcpy(gfit, args->fit, sizeof(fits));
			com.seq.current = UNRELATED_IMAGE;
			create_uniq_from_gfit(strdup(_("Pixel Math result")), FALSE);
		}
		else clearfits(args->fit);
		free(args->fit);
		free(args);
	}
	else {
		/* Still holding the writer lock: write the result into gfit now so the
		 * idle (end_pixel_math_operation) is left with GTK-only tasks.
		 * Use offsetof to avoid overwriting gfit->rwlock itself. */
		if (!failed) {
			clearfits(gfit);
			gui_iface.invalidate_histogram();
			memcpy(gfit, args->fit, offsetof(fits, rwlock));
			icc_auto_assign(gfit, ICC_ASSIGN_ON_COMPOSITION);
			com.seq.current = UNRELATED_IMAGE;
			create_uniq_from_gfit(strdup(_("Pixel Math result")), FALSE);
		} else {
			clearfits(args->fit);
		}
		if (rwlocked)
			g_rw_lock_writer_unlock(&gfit->rwlock);
		gui_iface.execute_idle_sync(end_pixel_math_operation, args);
	}
	return GINT_TO_POINTER((gint)failed);
}

/* ── Variable image management (called from Python bridge / headless) ─── */

int load_pm_var(const gchar *var, int index, int *w, int *h, int *c) {
	if (index > MAX_IMAGES - 1) {
		siril_log_message(_("A maximum of %d images can be used in a single expression.\n"), MAX_IMAGES);
		return 1;
	}

	image_type imagetype;
	char *realname = NULL;
	if (stat_file(var, &imagetype, &realname)) {
		siril_log_error(_("File not found or not supported: %s\n"), var);
		*w = *h = *c = -1;
		return 1;
	}
	int load_retval = any_to_fits(imagetype, realname, &var_fit[index], FALSE, TRUE);
	free(realname);
	if (load_retval) {
		*w = *h = *c = -1;
		return 1;
	}
	debayer_if_needed(imagetype, &var_fit[index], FALSE);
	*w = var_fit[index].rx;
	*h = var_fit[index].ry;
	*c = var_fit[index].naxes[2];
	return 0;
}

void free_pm_var(int nb) {
	for (int i = 0; i < nb; i++) {
		clearfits(&var_fit[i]);
		var_fit_mask[i] = FALSE;
	}
}

int replace_t_with_gfit(struct pixel_math_data *args) {
	int retval = 0;
	if (args->single_rgb == FALSE &&
		((args->expression1 && g_strstr_len(args->expression1, -1, "$T") != NULL) ||
		 (args->expression2 && g_strstr_len(args->expression2, -1, "$T") != NULL) ||
		 (args->expression3 && g_strstr_len(args->expression3, -1, "$T") != NULL))) {
		retval = 1;
	}

	if ((args->expression1 && g_strstr_len(args->expression1, -1, "$T") != NULL) ||
		(args->expression2 && g_strstr_len(args->expression2, -1, "$T") != NULL) ||
		(args->expression3 && g_strstr_len(args->expression3, -1, "$T") != NULL)) {

		const gchar *pattern = "\\$T";
		const gchar *replacement = T_CURRENT;

		args->has_gfit = TRUE;

		if (args->expression1 != NULL) {
			gchar *new_expression1 = g_regex_replace(g_regex_new(pattern, 0, 0, NULL), args->expression1, -1, 0, replacement, 0, NULL);
			g_free(args->expression1);
			args->expression1 = new_expression1;
		}

		if (args->expression2 != NULL) {
			gchar *new_expression2 = g_regex_replace(g_regex_new(pattern, 0, 0, NULL), args->expression2, -1, 0, replacement, 0, NULL);
			g_free(args->expression2);
			args->expression2 = new_expression2;
		}

		if (args->expression3 != NULL) {
			gchar *new_expression3 = g_regex_replace(g_regex_new(pattern, 0, 0, NULL), args->expression3, -1, 0, replacement, 0, NULL);
			g_free(args->expression3);
			args->expression3 = new_expression3;
		}
	}

	return retval;
}

const gchar *pm_get_variable_name(int i) {
	if (i < 0 || i >= (int)MAX_IMAGES) return NULL;
	return variables[i];
}

int pm_get_max_images(void) {
	return (int)MAX_IMAGES;
}
