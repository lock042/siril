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

#include "HTMesh.h"
#include "MeshBuffer.h"
#include "htmesh_wrapper.h"
#include <stdlib.h>

int get_htm_index_for_coords(double ra, double dec, int levels) {
	HTMesh mesh = HTMesh(levels, 0, 1);
	return mesh.index(ra, dec);
}

int get_htm_indices_around_target(double ra, double dec, double radius, int levels, int **trixels, int *nb_trixels) {
	*trixels = NULL;
	HTMesh mesh = HTMesh(levels, 0, 1);
	// intersect also has an overload that takes the 4 corners of the image too
	mesh.intersect(ra, dec, radius, (BufNum)0);
	auto buffer = mesh.meshBuffer(0);
	*nb_trixels = buffer->size();
	if (*nb_trixels > 0) {
		*trixels = (int *)malloc(*nb_trixels * sizeof(int));
		if (!*trixels) {
			*nb_trixels = 0;
		} else {
			for (int i = 0; i < *nb_trixels; i++) {
				(*trixels)[i] = buffer->buffer()[i];
			}
		}
	}

	// an iterator could also have been used, see the example in
	// https://invent.kde.org/education/kstars/-/blob/master/kstars/skycomponents/skymesh.cpp
	return *nb_trixels <= 0;
}
