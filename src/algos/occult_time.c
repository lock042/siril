/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "algos/PSF.h"
#include "algos/photometry.h"
#include "algos/astrometry_solver.h"
#include "algos/statistics_float.h"
#include "algos/siril_wcs.h"
#include "algos/search_objects.h"
#include "algos/comparison_stars.h"
#include "gui/PSF_list.h"
#include "gui/plot.h"
#include "gui/image_display.h"
#include "gui/siril_plot.h"
#include "io/sequence.h"
#include "io/siril_plot.h"
#include "opencv/opencv.h"