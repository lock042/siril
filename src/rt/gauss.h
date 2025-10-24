/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is siril_free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

enum eGaussType {GAUSS_STANDARD, GAUSS_MULT, GAUSS_DIV};

#ifdef __cplusplus
void gaussianBlur(float** src, float** dst, const int W, const int H, const double sigma, int max_threads, bool useBoxBlur = false, eGaussType gausstype = GAUSS_STANDARD, float** buffer2 = nullptr);

extern "C" {
void gaussianBlurC(float** src, float** dst, const int W, const int H, const double sigma, int threads);
}

#else

void gaussianBlurC(float** src, float** dst, const int W, const int H, const double sigma, int threads);

#endif
