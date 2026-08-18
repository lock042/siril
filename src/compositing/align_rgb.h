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

#ifndef SRC_COMPOSITING_ALIGN_RGB_H_
#define SRC_COMPOSITING_ALIGN_RGB_H_

/* Values are indices into the reg_methods[] table in align_rgb.c */
typedef enum {
	RGBALIGN_PSF = 0,
	RGBALIGN_DFT = 1,
	RGBALIGN_GLOBAL = 2,
	RGBALIGN_KOMBAT = 3,
} rgb_align_method;

const char *rgb_align_method_name(rgb_align_method m);
gboolean rgb_align_prerequisites_met(rgb_align_method m);
int rgb_align(rgb_align_method m);

#endif /* SRC_COMPOSITING_ALIGN_RGB_H_ */
