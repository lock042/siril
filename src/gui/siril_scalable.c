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

#include "core/siril.h"
#include "siril_scalable.h"

static double dpi = 0.;
static int scale = 0;

void set_DPInScale(const double newDPI, const int newScale) {
	if (!com.pref.pseudo_HiDPISupport) {
		scale = 1;
		dpi = baseDPI;
		return;
	}

	if (scale != newScale || (scale == 1 && dpi != newDPI)) {
		// reload all images
		scale = newScale;
		// HOMBRE: On windows, if scale = 2, the dpi is non significant, i.e. should be considered = 192 ; don't know for linux/macos
		dpi = newDPI;
		if (scale == 1) {
			if (dpi >= baseHiDPI) {
				scale = 2;
			}
		} else if (scale == 2) {
			if (dpi < baseHiDPI) {
				dpi *= 2.;
			}
		}
	}
}

double getDPI() {
	return dpi;
}

int getScale() {
	return scale;
}

void siril_scalable_init() {
	dpi = 0.;
	scale = 0;

	GdkMonitor *gdk_monitor = gdk_display_get_primary_monitor(gdk_display_get_default());
	set_DPInScale(gdk_screen_get_resolution(gdk_screen_get_default()),
			max((int )com.initial_GdkScale, gdk_monitor_get_scale_factor(gdk_monitor)));
	siril_debug_print("dpi=%lf, scale=%d\n", dpi, scale);
}
