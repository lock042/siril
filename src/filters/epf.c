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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"

#include "core/fits_region.h"
#include "filters/epf.h"
#include "core/op_descriptors.h"
#include "core/nde/nde_history.h"

/* NDE serializers (phase 4.5 Convention 1 — file operands).
 *
 * The bilateral filter and self-guided filter take no external file: bilateral
 * (guidefit == NULL) serializes plain POD; the self-guided filter emits the
 * token guide=self and replay_pre points guidefit at the replay target.  A
 * guided filter with a SEPARATE guide image pins that file by path (+hash/size
 * at capture) and replay_pre readfits() it into a fresh heap guidefit (freed by
 * free_epf_args via guide_needs_freeing).  At capture a self-guide is detected
 * by guidefit == fit (both the command and GUI self-guide sites set that). */
static gchar *epf_serialize(gconstpointer user) {
	const struct epfargs *p = user;
	GString *kv = nde_kv_start();
	nde_kv_add_double(kv, "d", p->d);
	nde_kv_add_double(kv, "sigma_col", p->sigma_col);
	nde_kv_add_double(kv, "sigma_space", p->sigma_space);
	nde_kv_add_double(kv, "mod", p->mod);
	/* on-disk value: enum order is frozen by the NDE format — do not reorder */
	nde_kv_add_int(kv, "filter", p->filter);
	if (p->filter == EP_GUIDED && p->guidefit) {
		if (p->guidefit == p->fit) {
			nde_kv_add_str(kv, "guide", "self");
		} else {
			nde_kv_add_str(kv, "guide", "file");
			nde_kv_add_operand(kv, p->guide_path);
		}
	}
	return nde_kv_end(kv);
}

static gpointer epf_deserialize(const gchar *blob, int version) {
	if (version > op_desc_epf.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	double d, sigma_col, sigma_space, mod;
	gint64 filter;
	NDE_KV_REQUIRE(kv, nde_kv_get_double(kv, "d", &d) &&
	                   nde_kv_get_double(kv, "sigma_col", &sigma_col) &&
	                   nde_kv_get_double(kv, "sigma_space", &sigma_space) &&
	                   nde_kv_get_double(kv, "mod", &mod) &&
	                   nde_kv_get_int(kv, "filter", &filter));
	const char *guide = nde_kv_get_str(kv, "guide");
	/* A guided record must carry a guide token (self/file); a file guide must
	 * carry a path.  A record without these (pre-phase-4.5 guided capture) is
	 * not replayable. */
	if (filter == EP_GUIDED) {
		NDE_KV_REQUIRE(kv, guide);
		if (!strcmp(guide, "file")) {
			const char *path = nde_kv_get_str(kv, "operand_path");
			NDE_KV_REQUIRE(kv, path && *path);
		}
	}
	struct epfargs *p = new_epf_args();
	if (p) {
		p->d = d;
		p->sigma_col = sigma_col;
		p->sigma_space = sigma_space;
		p->mod = mod;
		p->filter = (ep_filter_t)filter;
		if (guide && !strcmp(guide, "file"))
			p->guide_path = strdup(nde_kv_get_str(kv, "operand_path"));
	}
	g_hash_table_unref(kv);
	return p;
}

/* replay_pre: bind guidefit for a guided filter.  self → point at the target;
 * file → verify + readfits into a fresh heap fits (guide_needs_freeing so
 * free_epf_args releases it).  Bilateral / no-guide records skip this. */
static int epf_replay_pre(gpointer user, GHashTable *kv, fits *target) {
	struct epfargs *p = user;
	if (p->filter != EP_GUIDED)
		return 0;
	const char *guide = nde_kv_get_str(kv, "guide");
	if (guide && !strcmp(guide, "self")) {
		p->guidefit = target;   /* not owned; guide_needs_freeing stays FALSE */
		return 0;
	}
	const char *path = nde_operand_get(kv);
	if (!path)
		return 1;
	fits *guidefit = calloc(1, sizeof(fits));
	if (!guidefit)
		return 1;
	/* both construction sites load the guide image with force_float == FALSE */
	if (readfits(path, guidefit, NULL, FALSE)) {
		free(guidefit);
		return 1;
	}
	p->guidefit = guidefit;
	p->guide_needs_freeing = TRUE;   /* freed by free_epf_args */
	return 0;
}

/* Op descriptor — single source of truth for this operation (op_descriptor.h) */
/* Input context a region preview needs (op_descriptor.h).
 *
 * Both filters are fixed-kernel, so the halo is their radius and it is exact:
 * with that much real context around the requested rectangle, OpenCV's border
 * extrapolation only ever touches pixels the region does not keep.
 *
 * The radii mirror what edge_preserving_filter() actually passes down, so this
 * has to track those adjustments rather than the raw GUI values:
 *
 *   bilateral — cv::bilateralFilter takes d as an int and uses radius = d/2,
 *               or cvRound(sigma_space * 1.5) when d <= 0.
 *   guided    — the filter is called with r = (d ? d : sigma_space) / 3, and
 *               its box filter is applied TWICE (once to build a and b, once
 *               to average them), so the support is 2r, not r.
 *
 * The 2r is the sort of thing that looks like an off-by-a-factor guess, which
 * is why roi_halo_test.c asserts the deviation is zero here and non-zero one
 * pixel below.
 */
static int epf_roi_halo(gconstpointer user) {
	const struct epfargs *p = user;
	if (!p)
		return 0;
	if (p->filter == EP_BILATERAL) {
		const int d = (int)p->d;
		const int radius = (d > 0) ? d / 2 : (int)lround(p->sigma_space * 1.5);
		return radius > 1 ? radius : 1;
	}
	/* guided */
	const double base = (p->d != 0.0) ? p->d : p->sigma_space;
	const int r = (int)(base / 3.0);
	return 2 * (r > 1 ? r : 1);
}

const op_descriptor op_desc_epf = {
	.id = "filters.epf", .version = 1,
	.image_hook = epf_image_hook,
	.log_hook = epf_log_hook,
	.description = N_("Edge Preserving Filter"),
	.mem_ratio = 3.0f,
	.flags = OP_MASK_CAPABLE | OP_ROI_CAPABLE,
	.roi_halo = epf_roi_halo, .roi_halo_exact = TRUE,
	.serialize = epf_serialize, .deserialize = epf_deserialize,
	.replay_pre = epf_replay_pre,
};

/*****************************************************************************
 *      E P F      A L L O C A T O R   A N D   D E S T R U C T O R          *
 ****************************************************************************/

/* Allocator for epfargs */
struct epfargs *new_epf_args() {
	struct epfargs *args = calloc(1, sizeof(struct epfargs));
	if (args) {
		args->destroy_fn = free_epf_args;
	}
	return args;
}

/* Destructor for epfargs */
void free_epf_args(void *ptr) {
	struct epfargs *args = (struct epfargs *)ptr;
	if (!args)
		return;

	if (args->guide_needs_freeing && args->guidefit) {
		clearfits(args->guidefit);
		free(args->guidefit);
		args->guidefit = NULL;
	}
	free(args->guide_path);
	free(ptr);
}

gchar *epf_log_hook(gpointer p, log_hook_detail detail) {
	struct epfargs *args = (struct epfargs*) p;
	gchar *message = NULL;
	if (detail == SUMMARY)
		message = g_strdup_printf(_("%s filter: d=%.2f, sig_col=%.2f, sig_spatial=%.2f, mod=%.2f"),
								args->filter == EP_BILATERAL ? _("Bilateral") : _("Guided"),
								args->d, args->sigma_col, args->sigma_space, args->mod);
		else
		message = g_strdup_printf(_("%s filter: d=%.3f, sigma_col=%.3f, sigma_spatial=%.3f, mod=%.3f"),
								args->filter == EP_BILATERAL ? _("Bilateral") : _("Guided"),
								args->d, args->sigma_col, args->sigma_space, args->mod);
	return message;
}

static int edge_preserving_filter(struct epfargs *args) {
	fits *fit = args->fit;
	fits *guide = args->guidefit;
	double d = args->d;
	double sigma_col = args->sigma_col;
	double sigma_space = args->sigma_space;
	double mod = args->mod;
	ep_filter_t filter_type = args->filter;

	if (sigma_col <= 0.0 || (sigma_space <= 0.0 && filter_type == EP_BILATERAL))
		return 1;
	sigma_col /= 100.0;

	if (fit->naxes[2] == 1) {
		// This makes the settings behave more consistently between color and mono images
		sigma_col /= 25.0;
		if (filter_type == EP_GUIDED) {
			sigma_col /= 10.0;
		}
	}
	if (filter_type == EP_GUIDED) {
		// This makes the settings behave more consistently between the two filter types
		sigma_col /= 5.0;
		if (d == 0)
			d = sigma_space;
		d /= 3.0;
	}
	if (fit->type == DATA_FLOAT) {
		// This makes the settings behave more consistently between 16-bit and 32-bit images
		sigma_col *= 2.0;
	}

	// cv::BilateralFilter() only works on 8u and 32f images, so we convert 16-bit to 32-bit
	size_t ndata = fit->rx * fit->ry * fit->naxes[2];
	data_type orig_type = fit->type;
	if (orig_type == DATA_USHORT) {
		fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
	}
	fits orig = { 0 }; // for use with modulation
	if (mod < (1.0 - DBL_EPSILON)) {
		copyfits(fit, &orig, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	}
	double eps = sigma_col * sigma_col;
	fits *guide_roi = NULL, *guidance = NULL;
	gboolean roi_fitting_needed;
	switch (filter_type) {
		case EP_BILATERAL:
			cvBilateralFilter(fit, d, eps, sigma_space);
			break;
		case EP_GUIDED:
			/* On a region preview `fit` is a crop of the image, so a
			 * separate guide image has to be cropped to the same rectangle
			 * before cvGuidedFilter can pair them.  A self-guided filter
			 * needs nothing: the hook already pointed the guide at the
			 * crop. */
			guide_roi = calloc(1, sizeof(fits));
			roi_fitting_needed = args->for_roi && guide != fit;
			if (roi_fitting_needed
			    && crop_fits_region(guide, &args->roi_rect, guide_roi)) {
				siril_log_error(_("Guide image does not cover the region of "
				                  "interest.\n"));
				clearfits(guide_roi);
				free(guide_roi);
				if (orig_type == DATA_USHORT)
					fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
				return 1;
			}
			guidance = roi_fitting_needed ? guide_roi : guide;
			cvGuidedFilter(fit, guidance, d, eps);
			clearfits(guide_roi);
			free(guide_roi);
			break;
	}
	if (mod < (1.0 - DBL_EPSILON)) {
		for (size_t j = 0 ; j < fit->ry ; j++) {
			size_t offset = j * fit->rx;
			float modrem = 1.0 - mod;
			if (fit->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread) collapse(2) if (com.max_thread > 1)
#endif
				for (size_t k = 0 ; k < fit->naxes[2] ; k++) {
					for (size_t i = 0 ; i < fit->rx ; i++) {
						fit->fpdata[k][i + offset] = (float) mod * fit->fpdata[k][i + offset] + modrem * orig.fpdata[k][i + offset];
					}
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread) collapse(2) if (com.max_thread > 1)
#endif
				for (size_t k = 0 ; k < fit->naxes[2] ; k++) {
					for (size_t i = 0 ; i < fit->rx ; i++) {
						fit->pdata[k][i + offset] = roundf_to_WORD((float)(mod * fit->pdata[k][i + offset] + modrem * orig.pdata[k][i + offset]));
					}
				}
			}
		}
		clearfits(&orig);
	}
	if (orig_type == DATA_USHORT) {
		fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
	}

	/* No display refresh here: generic_image_worker notifies once the
	 * result has been installed. */
	return 0;
}

/* The actual EPF processing hook */
int epf_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct epfargs *params = (struct epfargs *)args->user;
	if (!params)
		return 1;
	/* Plumb the worker's `fit` through, as the hook contract requires — the
	 * worker may hand us a private copy or a region crop rather than gfit.
	 * A self-guided filter is identified by guidefit == fit and must follow
	 * it, or a region run would pair a crop with a full-size guide.
	 * Restored afterwards: params outlives the hook (serialize reads it). */
	fits *saved_fit = params->fit, *saved_guide = params->guidefit;
	const gboolean self_guided = (params->guidefit == params->fit);
	params->fit = fit;
	if (self_guided)
		params->guidefit = fit;
	params->for_roi = args->for_roi;
	params->roi_rect = args->roi_rect;
	int retval = edge_preserving_filter(params);
	params->fit = saved_fit;
	params->guidefit = saved_guide;
	return retval;
}

