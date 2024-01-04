/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_CORE_ARITHM_H_
#define SRC_CORE_ARITHM_H_

/* operations on image data */
typedef enum {
	OPER_ADD,
	OPER_SUB,
	OPER_MUL,
	OPER_DIV
} image_operator;

typedef struct blend_data {
	float sf[3]; // Luminance stretched values
	float tf[3]; // Independently stretched values
	gboolean do_channel[3]; // Whether or not to do each channel
} blend_data;

int soper(fits *a, float scalar, image_operator oper, gboolean conv_to_float);
int imoper(fits *a, fits *b, image_operator oper, gboolean allow_32bits);
int addmax(fits *a, fits *b);
int siril_fdiv(fits *a, fits *b, float scalar, gboolean allow_32bits);
int siril_ndiv(fits *a, fits *b);

int soper_unscaled_div_ushort_to_float(fits *a, int scalar);

void rgbblend(blend_data *data, float* r, float* g, float* b, float m_CB);
#endif /* SRC_CORE_ARITHM_H_ */
