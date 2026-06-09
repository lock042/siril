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

#include <math.h>

#include "core/siril.h"
#include "core/arithm.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/gui_iface.h"
#include "io/image_format_fits.h"
#include "rgradient.h"

/*****************************************************************************
 *      R G R A D I E N T   A L L O C A T O R   A N D   D E S T R U C T O R  *
 ****************************************************************************/

/* Allocator for rgradient_data */
struct rgradient_data *new_rgradient_data() {
	struct rgradient_data *args = calloc(1, sizeof(struct rgradient_data));
	if (args) {
		args->destroy_fn = free_rgradient_data;
	}
	return args;
}

/* Destructor for rgradient_data */
void free_rgradient_data(void *ptr) {
	struct rgradient_data *args = (struct rgradient_data *)ptr;
	if (!args)
		return;
	free(ptr);
}

/*****************************************************************************
 *      R G R A D I E N T   C O R E   P R O C E S S I N G                    *
 ****************************************************************************/

static void to_polar(int x, int y, point center, double *r, double *theta) {
	double dx = x - center.x;
	double dy = y - center.y;
	*r = sqrt(dx * dx + dy * dy);
	*theta = atan2(dy, dx);
}

static void to_cartesian(double r, double theta, point center, point *p) {
	p->x = center.x + r * cos(theta);
	p->y = center.y + r * sin(theta);
}

gchar *rgradient_log_hook(gpointer p, log_hook_detail detail) {
	struct generic_img_args *ga = (struct generic_img_args *) p;
	struct rgradient_data *args = (struct rgradient_data *) ga->user;
	gchar *message = g_strdup_printf(_("Rotational gradient radial shift: %.3f, rotational shift: %.3f, centre (%.1f, %.1f)"),
			args->dR, args->da, args->xc, args->yc);
	return message;
}

static int apply_rgradient_filter(struct rgradient_data *args) {
	fits *fit = args->fit;
	int retval = 0;
	const point center = {args->xc, args->yc};
	const double dAlpha = M_PI / 180.0 * args->da;
	gboolean was_ushort;
	fits imA = { 0 }, imB = { 0 };

	int cur_nb = 0;
	const double total = fit->ry * fit->naxes[2];

	was_ushort = fit->type == DATA_USHORT;

	/* convenient transformation to not inverse y sign */
	fits_flip_top_to_bottom(fit);

	if (was_ushort) {
		const long n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
		float *newbuf = ushort_buffer_to_float(fit->data, n);
		if (!newbuf) {
			retval = 1;
			goto end_rgradient;
		}
		fit_replace_buffer(fit, newbuf, DATA_FLOAT);
	}

	retval = copyfits(fit, &imA, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (retval) {
		retval = 1;
		goto end_rgradient;
	}

	retval = copyfits(&imA, &imB, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (retval) {
		retval = 1;
		goto end_rgradient;
	}

	int layer, y;
	float global_min = FLT_MAX;

	// Processing the image layers
	for (layer = 0; layer < fit->naxes[2]; layer++) {
		float *gbuf = fit->fpdata[layer];
		float *Abuf = imA.fpdata[layer];
		float *Bbuf = imB.fpdata[layer];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(min:global_min)
#endif
		for (y = 0; y < fit->ry; y++) {
			size_t i = y * fit->rx;
#ifdef _OPENMP
			#pragma omp critical
#endif
			{
				gui_iface.set_progress(cur_nb / total, NULL);
				cur_nb++;
			}

			for (int x = 0; x < fit->rx; x++) {
				float buf = gbuf[i] + gbuf[i];

				double r, theta;
				point delta1, delta2;

				to_polar(x, y, center, &r, &theta);

				// Positive differential
				to_cartesian(r - args->dR, theta + dAlpha, center, &delta1);
				buf -= bilinear(Abuf, fit->rx, fit->ry, delta1.x, delta1.y);

				// Negative differential
				to_cartesian(r - args->dR, theta - dAlpha, center, &delta2);
				buf -= bilinear(Bbuf, fit->rx, fit->ry, delta2.x, delta2.y);

				gbuf[i] = buf > 1.f ? 1.f : buf;

				// Update global minimum value
				if (gbuf[i] < global_min) {
					global_min = gbuf[i];
				}

				i++;
			}
		}
	}

	retval = soper(fit, global_min, OPER_SUB, TRUE);
	if (retval) {
		retval = 1;
		goto end_rgradient;
	}

end_rgradient:
	fits_flip_top_to_bottom(fit);
	clearfits(&imA);
	clearfits(&imB);

	if (!retval && was_ushort) {
		const long n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
		WORD *newbuf = float_buffer_to_ushort(fit->fdata, n);
		if (!newbuf)
			retval = 1;
		else
			fit_replace_buffer(fit, newbuf, DATA_USHORT);
	}

	return retval;
}

/* The actual rgradient processing hook */
int rgradient_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct rgradient_data *params = (struct rgradient_data *)args->user;
	if (!params)
		return 1;
	params->fit = fit;
	return apply_rgradient_filter(params);
}

/* GUI callbacks moved to src/gui/rgradient.c */
