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
 * Micro-benchmark for imgops::dft_patch, the per-patch real<->complex DFT used
 * by DA3D's hot loop. The Phase-0 spike concluded a dedicated, plan-caching
 * helper was needed (img_t's general fft allocates per call); this confirms the
 * helper amortises plan creation and runs a forward+inverse transform in the
 * low-microsecond range. Exits non-zero on a gross regression.
 */

#include <chrono>
#include <cstdio>
#include "algos/img_t/image.hpp"
#include "algos/img_t/image_dft.hpp"

cominfo com;
int sequence_is_running = 0;

static double bench(int w, int h, int d, int iters) {
	imgops::dft_patch p(w, h, d);
	for (int i = 0; i < d * h * w; ++i)
		p.space().data[i] = 0.001f * (i % 997);
	// warmup
	p.to_freq();
	p.to_space();
	auto t0 = std::chrono::steady_clock::now();
	for (int it = 0; it < iters; ++it) {
		p.to_freq();
		// trivial frequency-domain op so the loop can't be elided
		p.freq().data[0] *= 0.999f;
		p.to_space();
	}
	auto t1 = std::chrono::steady_clock::now();
	double us = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / 1000.0;
	return us / iters;  // microseconds per forward+inverse
}

int main(void) {
	com.max_thread = 1;
	const int iters = 20000;
	double mono = bench(64, 64, 1, iters);   // DA3D default patch size (r=31 -> s=64), mono
	double colour = bench(64, 64, 3, iters); // 3-channel
	printf("dft_patch 64x64x1: %.2f us / (fwd+inv)\n", mono);
	printf("dft_patch 64x64x3: %.2f us / (fwd+inv)\n", colour);
	// Generous ceiling: a 64x64 r2c+c2r is a few us on any modern CPU. Guard
	// only against a gross regression (e.g. per-call plan creation would be
	// orders of magnitude slower).
	if (mono > 500.0 || colour > 1500.0) {
		fprintf(stderr, "dft_patch perf regression: mono=%.1fus colour=%.1fus\n", mono, colour);
		return 1;
	}
	return 0;
}
