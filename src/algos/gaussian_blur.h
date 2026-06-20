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
#ifndef SRC_ALGOS_GAUSSIAN_BLUR_H_
#define SRC_ALGOS_GAUSSIAN_BLUR_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Separable Young-van Vliet recursive Gaussian blur, monochrome float.
 * src and dst are arrays of H row pointers each of width W. src == dst is
 * supported (in-place). Portable OpenMP implementation (no x86 intrinsics). */
void gaussian_blur_mono(float **src, float **dst, const int W, const int H, const double sigma, int threads);

#ifdef __cplusplus
}
#endif

#endif /* SRC_ALGOS_GAUSSIAN_BLUR_H_ */
