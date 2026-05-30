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
 * Shared image operations on img_t<T>, factored out so the da3d and nlbayes
 * denoisers can drop their bespoke image types (branch harmonize_img_t).
 *
 * img_t is planar (channel-major): operator()(x, y, c) -> data[x + w*(y + h*c)].
 * Coordinates follow img_t convention: x = column, y = row, c = channel.
 *
 * The boundary convention here is half-sample symmetric reflection (the edge
 * pixel is NOT duplicated): index -1 maps to 0, index w maps to w-1. This
 * matches both da3d's SymmetricCoordinate (Utils.cpp) and nlbayes'
 * symetrizeImage (LibImages.cpp), so a single primitive serves both.
 */

#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include "image.hpp"
#include "image_boundary.hpp"  // symmetric_coordinate / reflect_whole_sample
#include "image_expr.hpp"  // reduce_axis / broadcast_axis / AXIS_*

namespace imgops {

//! True if every pixel has identical values across all channels (so a
//! "colour" image actually carries monochrome content). Mono is trivially true.
template <typename T>
bool is_monochrome(const img_t<T>& u) {
    if (u.d <= 1)
        return true;
    for (int y = 0; y < u.h; ++y) {
        for (int x = 0; x < u.w; ++x) {
            const T v = u(x, y, 0);
            for (int c = 1; c < u.d; ++c) {
                if (u(x, y, c) != v)
                    return false;
            }
        }
    }
    return true;
}

//! Collapse all channels to a single-channel mean. Accumulates in double to
//! match da3d::makeMonochrome exactly.
template <typename T>
img_t<T> make_monochrome(const img_t<T>& u) {
    img_t<T> out(u.w, u.h, 1);
    for (int y = 0; y < u.h; ++y) {
        for (int x = 0; x < u.w; ++x) {
            double v = 0.0;
            for (int c = 0; c < u.d; ++c)
                v += u(x, y, c);
            out(x, y) = static_cast<T>(v / u.d);
        }
    }
    return out;
}

//! Pad every side by `border` pixels using half-sample symmetric reflection.
//! Equivalent to nlbayes symetrizeImage(forward) / addBoundary.
template <typename T>
img_t<T> pad_symmetric(const img_t<T>& in, int border) {
    img_t<T> out(in.w + 2 * border, in.h + 2 * border, in.d);
    for (int c = 0; c < in.d; ++c) {
        for (int y = 0; y < out.h; ++y) {
            const int sy = symmetric_coordinate(y - border, in.h);
            for (int x = 0; x < out.w; ++x) {
                const int sx = symmetric_coordinate(x - border, in.w);
                out(x, y, c) = in(sx, sy, c);
            }
        }
    }
    return out;
}

//! Inverse of pad_symmetric: crop `border` pixels off every side.
template <typename T>
img_t<T> unpad(const img_t<T>& in, int border) {
    img_t<T> out(in.w - 2 * border, in.h - 2 * border, in.d);
    for (int c = 0; c < in.d; ++c)
        for (int y = 0; y < out.h; ++y)
            for (int x = 0; x < out.w; ++x)
                out(x, y, c) = in(x + border, y + border, c);
    return out;
}

//! Mean of an image along one axis (AXIS_X/AXIS_Y/AXIS_D), materialised into a
//! new img_t whose extent along that axis is 1. Built on reduce_axis; this is
//! the NL-Bayes baricenter (mean over patches = mean along the patch axis).
template <typename T>
img_t<T> mean_axis(const img_t<T>& in, int axis) {
    const int n = (axis == AXIS_X) ? in.w : (axis == AXIS_Y) ? in.h : in.d;
    img_t<T> out(axis == AXIS_X ? 1 : in.w,
                 axis == AXIS_Y ? 1 : in.h,
                 axis == AXIS_D ? 1 : in.d);
    out.map(reduce_axis(in, axis,
                        std::function<T(T, T)>([](T a, T b) { return a + b; }),
                        T(0)));
    const T invn = T(1) / static_cast<T>(n);
    for (T& v : out.data)
        v *= invn;
    return out;
}

//! Tiling descriptor: number of tiles along y (rows) and x (columns).
struct tiling_t {
    int ny;  // tiles vertically (rows)
    int nx;  // tiles horizontally (columns)
};

//! Choose a square-ish tiling of an h x w image into at most `ntiles` tiles.
//! Faithful port of da3d::ComputeTiling (which works in rows/columns).
inline tiling_t compute_tiling(int h, int w, int ntiles) {
    const float best_r = std::sqrt(static_cast<float>(ntiles) * h / w);
    int r_low = static_cast<int>(best_r);
    int r_up = r_low + 1;
    if (r_low < 1) return {1, ntiles};       // single row of tiles
    if (r_up > ntiles) return {ntiles, 1};   // single column of tiles
    while (ntiles % r_low != 0) --r_low;
    while (ntiles % r_up != 0) ++r_up;
    if (r_up * r_low * w > ntiles * h)
        return {r_low, ntiles / r_low};
    return {r_up, ntiles / r_up};
}

//! Split into tiling.ny x tiling.nx tiles in lexicographic (row-major) order,
//! each grown by pad_before/pad_after pixels and filled with symmetric
//! reflection outside the source. Faithful port of da3d::SplitTiles.
template <typename T>
std::vector<img_t<T>> split_tiles(const img_t<T>& src, int pad_before,
                                  int pad_after, tiling_t tiling) {
    std::vector<img_t<T>> result;
    result.reserve(static_cast<size_t>(tiling.ny) * tiling.nx);
    for (int ty = 0; ty < tiling.ny; ++ty) {
        const int ystart = src.h * ty / tiling.ny - pad_before;
        const int yend = src.h * (ty + 1) / tiling.ny + pad_after;
        for (int tx = 0; tx < tiling.nx; ++tx) {
            const int xstart = src.w * tx / tiling.nx - pad_before;
            const int xend = src.w * (tx + 1) / tiling.nx + pad_after;
            img_t<T> tile(xend - xstart, yend - ystart, src.d);
            for (int c = 0; c < src.d; ++c)
                for (int y = ystart; y < yend; ++y)
                    for (int x = xstart; x < xend; ++x)
                        tile(x - xstart, y - ystart, c) =
                            src(symmetric_coordinate(x, src.w),
                                symmetric_coordinate(y, src.h), c);
            result.push_back(std::move(tile));
        }
    }
    return result;
}

//! Recompose tiles produced by split_tiles. Each entry is {value, weight}:
//! `value` holds the (already weight-multiplied) tile contribution and
//! `weight` is a single-channel per-pixel weight. Overlapping regions are
//! accumulated and divided by the summed weight. Faithful port of
//! da3d::MergeTiles. `out_w`/`out_h` are the full (unpadded) image size.
template <typename T>
img_t<T> merge_tiles(const std::vector<std::pair<img_t<T>, img_t<T>>>& tiles,
                     int out_w, int out_h, int pad_before, int pad_after,
                     tiling_t tiling) {
    const int channels = tiles[0].first.d;
    img_t<T> result(out_w, out_h, channels);
    img_t<T> weights(out_w, out_h, 1);
    result.set_value(T(0));
    weights.set_value(T(0));

    size_t idx = 0;
    for (int ty = 0; ty < tiling.ny; ++ty) {
        const int ystart = out_h * ty / tiling.ny - pad_before;
        const int yend = out_h * (ty + 1) / tiling.ny + pad_after;
        for (int tx = 0; tx < tiling.nx; ++tx, ++idx) {
            const int xstart = out_w * tx / tiling.nx - pad_before;
            const int xend = out_w * (tx + 1) / tiling.nx + pad_after;
            const img_t<T>& tval = tiles[idx].first;
            const img_t<T>& twgt = tiles[idx].second;
            const int y0 = std::max(0, ystart), y1 = std::min(out_h, yend);
            const int x0 = std::max(0, xstart), x1 = std::min(out_w, xend);
            for (int y = y0; y < y1; ++y) {
                for (int x = x0; x < x1; ++x) {
                    for (int c = 0; c < channels; ++c)
                        result(x, y, c) += tval(x - xstart, y - ystart, c);
                    weights(x, y) += twgt(x - xstart, y - ystart);
                }
            }
        }
    }
    for (int y = 0; y < out_h; ++y)
        for (int x = 0; x < out_w; ++x)
            for (int c = 0; c < channels; ++c)
                result(x, y, c) /= weights(x, y);
    return result;
}

}  // namespace imgops
