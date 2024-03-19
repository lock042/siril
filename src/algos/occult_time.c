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
#include <gtk/gtk.h>
#include "gui/utils.h"

#include "comparison_stars.h"
#include "occult_time.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "algos/search_objects.h"
#include "io/siril_catalogues.h"
#include "io/remote_catalogues.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"


#define BORDER_RATIO 0.10 // the amount of image that is considered at border
#define MAX_VAR_CAT 3 // max number of variable stars catalogues
#define MAX_SEPARATION  2 * 0.000277778 // the max sepration to consider stars are the same
const siril_cat_index var_cat[MAX_VAR_CAT] = { CAT_VSX, CAT_GCVS, CAT_VARISUM }; // the variable stars catalogues, order to be consistent with discard_var_catalogues



