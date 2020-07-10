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

#ifdef _WIN32
#include <windows.h>
#endif
#include <stdio.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "gui/progress_and_log.h"

#include "FITS_symlink.h"

static gboolean end_symlink_idle(gpointer p) {
	struct _symlink_data *args = (struct _symlink_data *) p;
	struct timeval t_end;

	if (!args->retval && get_thread_run() && args->nb_renamed_files > 1) {
		// load the sequence
		char *renamed_seqname = NULL;
		renamed_seqname = malloc(strlen(args->destroot) + 5);
		sprintf(renamed_seqname, "%s.seq", args->destroot);
		check_seq(0);
		if (renamed_seqname) {
			update_sequences_list(renamed_seqname);
			free(renamed_seqname);
		}
	}

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	stop_processing_thread();
	g_free(args->destroot);
	free(args);
	return FALSE;
}

gpointer symlink_thread_worker(gpointer p) {
	double progress = 0.0;
	struct _symlink_data *args = (struct _symlink_data *) p;
	unsigned int frame_index = 0;

	args->nb_renamed_files = 0;
	args->retval = 0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) \
	if(!args->input_has_a_seq && fits_is_reentrant())
	// we should run in parallel only when images are converted, not sequences
#endif
	for (int i = 0; i < args->total; i++) {
		if (args->retval || !get_thread_run()) {
			continue;
		}

		gchar *src_filename = args->list[i];
		const char *src_ext = get_filename_ext(src_filename);
		int index = args->input_has_a_seq ? frame_index : args->start + i;

		gchar *name = g_utf8_strrchr(src_filename, strlen(src_filename), G_DIR_SEPARATOR);
		gchar *msg_bar;
		if (name)
			msg_bar = g_strdup_printf(_("Converting %s..."), name + 1);
		else msg_bar = g_strdup_printf(_("Converting %s..."), src_filename);

		image_type imagetype = get_type_for_extension(src_ext);
		com.filter = (int) imagetype;
		if (imagetype != TYPEFITS) {
			args->retval = 1;
			g_free(msg_bar);
			continue;
		} else {
			gchar *dest_filename = g_strdup_printf("%s%05d%s", args->destroot,
					index, com.pref.ext);
#ifdef _WIN32
			wchar_t *wsrc, *wdst;

			wsrc = g_utf8_to_utf16(src_filename, -1, NULL, NULL, NULL);
			wdst = g_utf8_to_utf16(dest_filename, -1, NULL, NULL, NULL);

			gboolean ret = CreateSymbolicLinkW(wsrc, wdst, SYMBOLIC_LINK_FLAG_ALLOW_UNPRIVILEGED_CREATE);
			if (!ret) {
				siril_log_color_message(_("You should enable the Developer Mode in order to make symbolic link "
						"instead of simply copying files."), "red");
				copy_fits_from_file(src_filename, dest_filename);
			}

			g_free(wsrc);
			g_free(wdst);
#else
			symlink(src_filename, dest_filename);
#endif

			g_free(dest_filename);
			frame_index++;
		}

#ifdef _OPENMP
#pragma omp atomic
#endif
		progress += 1.0;
		set_progress_bar_data(msg_bar, progress / ((double) args->total));
		g_free(msg_bar);
#ifdef _OPENMP
#pragma omp atomic
#endif
		args->nb_renamed_files++;
	}

	g_dir_close(args->dir);
	for (int i = 0; i < args->total; i++)
		g_free(args->list[i]);
	if (args->retval)
		siril_log_message(_("Conversion ended with error, %d/%d input files converted\n"), args->nb_renamed_files, args->total);
	else {
		if (args->nb_renamed_files == args->total)
			siril_log_message(_("Conversion succeeded, %d/%d input files converted\n"), args->nb_renamed_files, args->total);
		else siril_log_message(_("Conversion aborted, %d/%d input files converted\n"), args->nb_renamed_files, args->total);
	}
	siril_add_idle(end_symlink_idle, args);
	return NULL;
}
