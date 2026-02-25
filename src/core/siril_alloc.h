// SPDX-License-Identifier: GPL-3.0-or-later
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


#ifndef SRC_CORE_SIRIL_ALLOC_H_
#define SRC_CORE_SIRIL_ALLOC_H_

void* siril_malloc(size_t size);
void* siril_calloc(size_t num, size_t size);
void siril_free(void* ptr);
// void* siril_realloc(void* ptr, size_t new_size);

#endif /* SRC_CORE_SIRIL_ALLOC_H_ */

