/*
MIT License

Copyright (c) 2018 Jérémy Anger, Gabriele Facciolo, Mauricio Delbracio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef SHEAR_IMAGE_H
#define SHEAR_IMAGE_H

#include <vector>
#include "image.hpp"

#include "angleSet.hpp"

template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& imgX, const img_t<T>& imgY,
                const std::vector<angle_t>& angleSet);

template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& img,
                const std::vector<angle_t>& angleSet);

#include "projectImage.cpp"

#endif

