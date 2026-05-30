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
 * Boundary reflection helpers shared by the image tiling code paths
 * (deconvolution's img_t::process_in_slices and the imgops tiling in
 * image_ops.hpp). Kept dependency-free (plain int math, no img_t) so it can be
 * included by image.hpp itself without a circular include.
 *
 * There are two distinct conventions in use; both live here as the single
 * source of truth:
 *
 *  - half-sample: the edge pixel is NOT duplicated. -1 -> 0, size -> size-1.
 *    Used by da3d (SymmetricCoordinate) and nlbayes (symetrizeImage), and by
 *    the imgops symmetric padding / tiling.
 *
 *  - whole-sample: reflect about the edge pixel. -1 -> 1, size -> size-2.
 *    Used by deconvolution's process_in_slices overlap padding. Single
 *    reflection only (valid because the overlap is always smaller than the
 *    image), matching the original inline behaviour exactly.
 */

#pragma once

namespace imgops {

//! Half-sample symmetric reflection into [0, size). Edge pixel not repeated.
inline int symmetric_coordinate(int pos, int size) {
    if (pos < 0) pos = -pos - 1;
    if (pos >= 2 * size) pos %= 2 * size;
    if (pos >= size) pos = 2 * size - 1 - pos;
    return pos;
}

//! Whole-sample symmetric reflection (reflect about the edge pixel), single
//! bounce. Reproduces the inline reflection in img_t::process_in_slices.
inline int reflect_whole_sample(int pos, int size) {
    if (pos < 0) return -pos;
    if (pos >= size) return 2 * size - pos - 2;
    return pos;
}

}  // namespace imgops
