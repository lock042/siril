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

#include <sys/time.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/utils.h"
#include "io/FITS_symlink.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "gui/image_display.h"

#define EXIT_TOKEN ":EXIT:"

static GFileMonitor *dirmon = NULL;
static gboolean do_links = FALSE;
static GThread *live_stacker_thread = NULL;
static GAsyncQueue *new_files_queue = NULL;

static gpointer live_stacker(gpointer arg);

void display(const char *str) {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("livest_label1"));
	gtk_label_set_text(label, str);
}

void on_livestacking_stop_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide((GtkWidget *)user_data);
}

void on_livestacking_player_hide(GtkWidget *widget, gpointer user_data) {
	GtkWindow *window;
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	gboolean is_fullscreen;

	window = GTK_WINDOW(lookup_widget("control_window"));

	GdkWindow *gdk_window = gtk_widget_get_window(GTK_WIDGET(window));
	is_fullscreen = gdk_window_get_state(gdk_window) & GDK_WINDOW_STATE_FULLSCREEN;

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);

	if (dirmon) {
		g_object_unref(dirmon);
		dirmon = NULL;
	}
	if (new_files_queue) {
		g_async_queue_push(new_files_queue, EXIT_TOKEN);
		g_async_queue_unref(new_files_queue);
		new_files_queue = NULL;
	}
	if (live_stacker_thread) {
		g_thread_join(live_stacker_thread);
		live_stacker_thread = NULL;
	}
}

static void file_changed(GFileMonitor *monitor, GFile *file, GFile *other,
		GFileMonitorEvent evtype, gpointer user_data) {
	if (evtype == G_FILE_MONITOR_EVENT_CREATED) {
		gchar *filename = g_file_get_basename(file);
		siril_debug_print("File %s created\n", filename);

		image_type type;
		stat_file(filename, &type, NULL);
		if (type != TYPEFITS) {
			gchar *str = g_strdup_printf(_("File not supported for live stacking: %s"), filename);
			display(str);
			g_free(str);
		} else {
			if (strncmp(filename, "live_stack", 10) &&
					strncmp(filename, "r_live_stack", 12) &&
					strncmp(filename, "result_live_stack", 17)) {
				if (!single_image_is_loaded())
					open_single_image(filename);
				g_async_queue_push(new_files_queue, filename);
			}
		}
		//g_free(filename);
	}
}

void on_livestacking_start() {
	display(_("Starting live stacking"));
	/* start monitoring CWD */
	GFile *cwd = g_file_new_for_path(com.wd);
	GError *err = NULL;
	dirmon = g_file_monitor_directory(cwd, G_FILE_MONITOR_WATCH_MOVES, NULL, &err);
	g_object_unref(cwd);
	if (err) {
		fprintf(stderr, "Unable to monitor CWD (%s): %s\n", com.wd, err->message);
		g_error_free(err);
		return;
	}

	g_signal_connect(G_OBJECT(dirmon), "changed", G_CALLBACK(file_changed), NULL);
	display(_("Live stacking waiting for files"));

	do_links = test_if_symlink_is_ok();
	reserve_thread();	// generic function fails otherwise

	new_files_queue = g_async_queue_new();
	live_stacker_thread = g_thread_new("live stacker", live_stacker, NULL);
}


static gpointer live_stacker(gpointer arg) {
	g_async_queue_ref(new_files_queue);
	int index = 1;
	gboolean first_loop = TRUE;	// only original images in the sequence
	do {
		gchar *filename = g_async_queue_pop(new_files_queue); // blocking
		if (!strcmp(filename, EXIT_TOKEN)) {
			siril_debug_print("Exiting thread\n");
			break;
		}
		siril_debug_print("Adding file to input sequence\n");
		gchar *target = g_strdup_printf("live_stack_%05d%s", index, com.pref.ext);
		if (symlink_uniq_file(filename, target, do_links)) {
			break;
		}
		g_free(filename);

		/* Create the sequence */
		siril_debug_print("Creating sequence %d\n", index);
		sequence seq;
		initialize_sequence(&seq, FALSE);

		seq.seqname = strdup("live_stack_");
		seq.beg = first_loop ? 1 : 0;
		seq.end = index;
		seq.fixed = 5;
		seq.reference_image = 0;
		if (first_loop) {
			if (buildseqfile(&seq, 1) || seq.number == 1) {
				index++;
				continue;
			}
			first_loop = FALSE;
		} else {
			seq.number = 2;
			seq.selnum = 2;
			seq.imgparam = malloc(2 * sizeof(imgdata));
			seq.imgparam[0].filenum = 0;
			seq.imgparam[0].incl = TRUE;
			seq.imgparam[0].date_obs = NULL;
			seq.imgparam[1].filenum = index;
			seq.imgparam[1].incl = TRUE;
			seq.imgparam[1].date_obs = NULL;
			writeseqfile(&seq);
		}
		/* TODO: cache this */
		if (seq_check_basic_data(&seq, FALSE) < 0) {
			break;
		}

		/* Register the sequence */
		struct registration_args reg_args = { 0 };
		reg_args.func = &register_star_alignment; // TODO: ability to choose a method
		reg_args.seq = &seq;
		reg_args.reference_image = 0;
		reg_args.process_all_frames = TRUE;
		reg_args.layer = (seq.nb_layers == 3) ? 1 : 0;
		reg_args.run_in_thread = FALSE;
		reg_args.follow_star = FALSE;
		reg_args.matchSelection = FALSE;
		//memcpy(&reg_args.selection, &com.selection, sizeof(rectangle));
		reg_args.x2upscale = FALSE;
		reg_args.cumul = FALSE;
		reg_args.min_pairs = 10;
		reg_args.translation_only = FALSE;
		reg_args.prefix = "r_";
		reg_args.load_new_sequence = FALSE;
		reg_args.interpolation = OPENCV_CUBIC;
		reg_args.func(&reg_args);
		if (reg_args.retval)
			break;

		gchar *result_filename = g_strdup_printf("live_stack_00000%s", com.pref.ext);

		free_sequence(&seq, FALSE);

		/* Stack the sequence */
		siril_debug_print("Stacking image %d\n", index);

		sequence *r_seq = readseqfile("r_live_stack_.seq");
		if (seq_check_basic_data(r_seq, FALSE) < 0) {
			free(r_seq);
			return FALSE;
		}

		struct stacking_args stackparam = { 0 };
		stackparam.method = stack_summing_generic; //stack_mean_with_rejection;
		stackparam.seq = r_seq;
		stackparam.ref_image = 0;
		// TODO Apply some filter on quality?
		stackparam.filtering_criterion = seq_filter_all;
		stackparam.filtering_parameter = 0.0;
		stackparam.nb_images_to_stack = 2;
		stack_fill_list_of_unfiltered_images(&stackparam);
		stackparam.description = describe_filter(r_seq, stackparam.filtering_criterion, stackparam.filtering_parameter);
		stackparam.output_filename = result_filename;
		stackparam.output_overwrite = TRUE;
		gettimeofday(&stackparam.t_start, NULL);

		stackparam.force_norm = FALSE;
		stackparam.output_norm = FALSE;
		stackparam.use_32bit_output = TRUE;
		stackparam.reglayer = (r_seq->nb_layers == 3) ? 1 : 0;

		main_stack(&stackparam);

		int retval = stackparam.retval;
		clean_end_stacking(&stackparam);
		free_sequence(r_seq, TRUE);
		free(stackparam.image_indices);
		free(stackparam.description);

		if (retval)
			break;

		//struct noise_data noise_args = { .fit = &gfit, .verbose = FALSE, .use_idle = FALSE };
		//noise(&noise_args);
		if (savefits(result_filename, &gfit)) {
			siril_log_color_message(_("Could not save the stacking result %s\n"),
					"red", result_filename);
			break;
		}

		/* Update display */
		queue_redraw(REMAP_ALL);

		index++;
	} while (1);

	siril_debug_print("===== exiting live stacking thread =====\n");

	// TODO: clean exit
	g_async_queue_unref(new_files_queue);
	return NULL;
}

