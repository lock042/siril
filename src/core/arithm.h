/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"

struct arithm_data {
	image_operator oper;
	fits *operand;
	sequence *seq;
	const gchar *seqEntry;
};


int soper(fits *a, float scalar, image_operator oper, gboolean conv_to_float);
int imoper(fits *a, fits *b, image_operator oper, gboolean allow_32bits);
int addmax(fits *a, fits *b);
int siril_fdiv(fits *a, fits *b, float scalar, gboolean allow_32bits);
int siril_ndiv(fits *a, fits *b);

void apply_arithm_extraction_to_sequence(struct arithm_data *arithm_args);

#endif /* SRC_CORE_ARITHM_H_ */
