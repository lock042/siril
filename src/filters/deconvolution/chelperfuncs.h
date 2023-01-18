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

// Helpful wrapper functions to call Siril functions without causing clashes with some of the img_expr functionality

#ifndef CPLUSPLUS_HELPER_FUNCTIONS_H
#define CPLUSPLUS_HELPER_FUNCTIONS_H

#ifdef __cplusplus
#define EXTERNC1 extern "C"
#else
#define EXTERNC1
#endif


EXTERNC1 void updateprogress(const char *text, double percent);
EXTERNC1 void sirillog(const char* text);
EXTERNC1 int is_thread_stopped();
EXTERNC1 int updatenoise(float *array, int nx, int ny, int nchans, double *noise);

#ifdef __cplusplus
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN int cppmaxthreads;
EXTERN unsigned cppfftwflags;
EXTERN double cppfftwtimelimit;
EXTERN int cppfftwmultithreaded;

#endif // CPLUSPLUS_HELPER_FUNCTIONS_H
