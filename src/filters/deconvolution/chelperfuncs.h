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

#ifndef CPLUSPLUS_HELPER_FUNCTIONS_H
#define CPLUSPLUS_HELPER_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

void magnify(float *y, const float *x, int W, int H, int pd, int w, int h, float n);
void shrink(float *y, float *x, int outw, int outh, int inw, int inh, float scale, float sigma);
void gaussblur(float*, float*, int, int, float);
void updateprogress(const char *text, double percent);
void sirillog(const char* text);
int is_thread_stopped();
int updatenoise(float *array, int nx, int ny, int nchans, double *noise);

#ifdef __cplusplus
}
#endif

#ifndef dontneedcppmaxthreads
static int cppmaxthreads;
#endif
#ifndef dontneedcppfftwflags
static unsigned cppfftwflags;
#endif
#ifndef dontneedcppfftwtimelimit
static double cppfftwtimelimit;
#endif
#ifndef dontneedcppfftwmultithreaded
static int cppfftwmultithreaded;
#endif

#endif // CPLUSPLUS_HELPER_FUNCTIONS_H
