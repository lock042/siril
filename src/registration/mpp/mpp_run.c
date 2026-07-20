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
 *
 * This implementation of multipoint registration & stacking is based on
 * PlanetarySystemStacker by Rolf Hempel:
 *     https://github.com/Rolf-Hempel/PlanetarySystemStacker
 */

/*
 * mpp_run_t alloc / free. Lives in its own TU (free of Siril runtime
 * dependencies like siril_log_message) so the sidecar I/O module can link
 * it without dragging mpp.cpp.o (which references com / siril_log_message
 * via the orchestrator stages) into test binaries.
 */

#include <stdlib.h>

#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_shift.h"

mpp_run_t *mpp_run_alloc(void) {
	return (mpp_run_t *) calloc(1, sizeof(mpp_run_t));
}

void mpp_run_free(mpp_run_t *run) {
	if (!run) return;
	free(run->cfg);
	free(run->quality);
	free(run->frame_brightness);
	free(run->included);
	free(run->global_shifts);
	free(run->mean_frame_data);
	free(run->best_frame_indices);
	mpp_ap_free(run->aps);
	mpp_shift_free(run->shifts);
	mpp_run_drop_cache(run);
	free(run);
}

/* mpp_cache_free lives in mpp.cpp (needs the C++ cv::Mat dtor). Tests
 * link with --unresolved-symbols=ignore-all and never populate the
 * cache, so the guard skips the call before the unresolved symbol
 * would matter. */
void mpp_run_drop_cache(mpp_run_t *run) {
	if (!run || !run->cache) return;
	mpp_cache_free(run->cache);
	run->cache = NULL;
}
