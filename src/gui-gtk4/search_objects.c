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

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/siril_log.h"
#include "algos/search_objects.h"
#include "algos/astrometry_solver.h"
#include "algos/siril_wcs.h"
#include "io/siril_catalogues.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/search_objects.h"

void search_object(GtkEditable *entry) {
	if (!has_wcs(gfit))
		return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	sky_object_query_args *query_args = init_sky_object_query();
	query_args->name = g_strdup(gtk_editable_get_text(entry));
	query_args->fit = gfit;
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free_sky_object_query(query_args);
		PRINT_ALLOC_ERR;
		return;
	}
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->image_hook = catsearch_image_hook;
	args->description = _("Catalog search");
	args->verbose = TRUE;
	args->command_updates_gfit = FALSE;
	args->command = FALSE;
	args->idle_function = end_process_catsearch;
	args->user = query_args;
	args->log_hook = catsearch_log_hook;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		siril_log_error(_("Error running image worker for catsearch\n"));
		return;
	}
}

/* GTK4: GtkSearchEntry no longer derives from GtkEntry — both implement
 * GtkEditable.  The .ui binds this handler to a GtkSearchEntry::activate. */
void on_search_objects_entry_activate(GtkSearchEntry *entry, gpointer user_data) {
	(void) user_data;
	search_object(GTK_EDITABLE(entry));
}
