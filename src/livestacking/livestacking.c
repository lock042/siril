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
#include "core/preprocess.h"
#include "gui/utils.h"
#include "io/FITS_symlink.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
/* global registration */
#include "algos/star_finder.h"
#include "registration/registration.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/match.h"
#include "opencv/opencv.h"
/* ******************* */
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "gui/image_display.h"

#define EXIT_TOKEN ":EXIT:"

static GFileMonitor *dirmon = NULL;
static gboolean do_links = FALSE;
static GThread *live_stacker_thread = NULL;
static GAsyncQueue *new_files_queue = NULL;

//static int registration_layer = 0;
//static fitted_PSF **ref_stars = NULL;
//static int nb_ref_stars = 0;
//static sequence *registered_seq = NULL;
static int seq_rx = -1, seq_ry = -1;
static struct star_align_data *sadata = NULL;
static struct preprocessing_data *prepro = NULL;

static void init_preprocessing();
static gpointer live_stacker(gpointer arg);
int star_align_prepare_hook(struct generic_seq_args *args);

void display(const char *str) {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("livest_label1"));
	gtk_label_set_text(label, str);
}

static void create_seq_of_2(sequence *seq, char *name, int index) {
	seq->seqname = strdup(name);
	seq->beg = 0;
	seq->end = index;
	seq->fixed = 5;
	seq->reference_image = 0;
	seq->number = 2;
	seq->selnum = 2;
	seq->imgparam = malloc(2 * sizeof(imgdata));
	seq->imgparam[0].filenum = 1;
	seq->imgparam[0].incl = TRUE;
	seq->imgparam[0].date_obs = NULL;
	seq->imgparam[1].filenum = index;
	seq->imgparam[1].incl = TRUE;
	seq->imgparam[1].date_obs = NULL;
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
	}
	if (live_stacker_thread) {
		g_thread_join(live_stacker_thread);
		live_stacker_thread = NULL;
	}
	if (prepro) {
		clear_preprocessing_data(prepro);
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
			g_free(filename);
		} else {
			if (strncmp(filename, "live_stack", 10) &&
					strncmp(filename, "r_live_stack", 12) &&
					strncmp(filename, "result_live_stack", 17)) {
				usleep(250000);	// temporary
				if (!single_image_is_loaded())
					open_single_image(filename);
				g_async_queue_push(new_files_queue, filename);
			}
			else g_free(filename);
		}
	}
}

void on_livestacking_start() {
	display(_("Starting live stacking"));
	com.pref.force_to_16bit = TRUE;	// otherwise we'll register and stack different types of images
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

	init_preprocessing();

	g_signal_connect(G_OBJECT(dirmon), "changed", G_CALLBACK(file_changed), NULL);
	display(_("Live stacking waiting for files"));

	do_links = test_if_symlink_is_ok();
	reserve_thread();	// generic function fails otherwise

	new_files_queue = g_async_queue_new();
	live_stacker_thread = g_thread_new("live stacker", live_stacker, NULL);
}

/*static int register_image(fits *fit, char *output_name) {
	fitted_PSF **stars;
	int nb_stars;

	stars = peaker(fit, registration_layer, &com.starfinder_conf, &nb_stars, NULL, FALSE, TRUE);
	siril_debug_print("Found %d stars in new image\n", nb_stars);

	if (!ref_stars) {
		if (nb_stars < AT_MATCH_MINPAIRS || !stars) {
			display(_("not enough stars found in image, adjust settings\n"));
			return 1;
		}
		ref_stars = stars;
		nb_ref_stars = nb_stars;
	} else {
		nb_stars = min(nb_stars, nb_ref_stars);
		double scale_min = 0.9;
		double scale_max = 1.1;
		int attempt = 1, retvalue = 1, nobj = 0;
		Homography H = { 0 };
		while (retvalue && attempt < NB_OF_MATCHING_TRY){
			retvalue = new_star_match(stars, sadata->refstars, nb_stars, nobj, scale_min, scale_max, &H, FALSE);
			if (attempt == 1) {
				scale_min = -1.0;
				scale_max = -1.0;
			} else {
				nobj += 50;
			}
			attempt++;
		}
		free_fitted_stars(stars);
		if (retvalue) {
			siril_log_color_message(_("Cannot perform star matching: try #%d. Image skipped\n"),
					"red", attempt);
			return 1;
		}
		if (H.pair_matched < AT_MATCH_MINPAIRS) {
			siril_log_color_message(_("Not enough star pairs (%d): Image skipped\n"),
					"red", H.pair_matched);
			return 1;
		}

		if (cvTransformImage(fit, H, FALSE, OPENCV_CUBIC)) {
			return 1;
		}

		if (!savefits(output_name, fit)) {
			return 1;
		}
	}

	return 0;
}*/

static int live_stacking_star_align_prepare(struct generic_seq_args *args) {
	if (!sadata || !sadata->refstars) {
		return star_align_prepare_hook(args);
	}
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	regargs->imgparam = calloc(args->nb_filtered_images, sizeof(imgdata));
	regargs->regparam = calloc(args->nb_filtered_images, sizeof(regdata));
	if (!regargs->imgparam  || !regargs->regparam) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	return 0;
}

static int start_global_registration(sequence *seq) {
	/* Register the sequence */
	struct registration_args reg_args = { 0 };
	reg_args.func = &register_star_alignment; // TODO: ability to choose a method
	reg_args.seq = seq;
	reg_args.reference_image = 0;
	reg_args.process_all_frames = TRUE;
	reg_args.layer = (seq->nb_layers == 3) ? 1 : 0;
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
	/*reg_args.func(&reg_args);
	if (reg_args.retval)
		return 1;*/

	struct generic_seq_args *args = create_default_seqargs(seq);
	if (!reg_args.process_all_frames) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = seq->selnum;
	}
	//args->compute_mem_limits_hook = star_align_compute_mem_limits;
	args->prepare_hook = live_stacking_star_align_prepare;
	args->image_hook = star_align_image_hook;
	//args->finalize_hook = star_align_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("Global star registration");
	args->has_output = TRUE;
	args->output_type = get_data_type(seq->bitpix);
	args->upscale_ratio = 1.0;
	args->new_seq_prefix = "r_";
	args->load_new_sequence = FALSE;
	args->already_in_a_thread = TRUE;
	if (!sadata) {
		sadata = calloc(1, sizeof(struct star_align_data));
		sadata->regargs = &reg_args;
	}
	args->user = sadata;

	return GPOINTER_TO_INT(generic_sequence_worker(args));
}

static void init_preprocessing() {
	/* copied from core/preprocess.c test_for_master_files() but modified to not check
	 * for some options and not check for image properties against gfit, which is not
	 * already loaded here. Error management is different too */
	prepro = calloc(1, sizeof(struct preprocessing_data));
	GtkToggleButton *tbutton = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		GtkEntry *entry = GTK_ENTRY(lookup_widget("darkname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			set_progress_bar_data(_("Opening dark image..."), PROGRESS_NONE);
			prepro->dark = calloc(1, sizeof(fits));
			if (!readfits(filename, prepro->dark, NULL, FALSE)) {
				prepro->use_dark = TRUE;
			} else {
				display(_("NOT USING DARK: cannot open the file"));
				free(prepro->dark);
				gtk_entry_set_text(entry, "");
				prepro->use_dark = FALSE;
			}
		}

		if (prepro->use_dark) {
			// cosmetic correction
			tbutton = GTK_TOGGLE_BUTTON(lookup_widget("cosmEnabledCheck"));
			prepro->use_cosmetic_correction = gtk_toggle_button_get_active(tbutton);

			if (prepro->use_cosmetic_correction) {
				tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkSigCold"));
				if (gtk_toggle_button_get_active(tbutton)) {
					GtkSpinButton *sigCold = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox"));
					prepro->sigma[0] = gtk_spin_button_get_value(sigCold);
				} else prepro->sigma[0] = -1.0;

				tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkSigHot"));
				if (gtk_toggle_button_get_active(tbutton)) {
					GtkSpinButton *sigHot = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox"));
					prepro->sigma[1] = gtk_spin_button_get_value(sigHot);
				} else prepro->sigma[1] = -1.0;
			}
		}
	}

	if (prepro->use_dark) {
		struct generic_seq_args generic = { .user = prepro };
		if (prepro_prepare_hook(&generic)) {
			clear_preprocessing_data(prepro);
			free(prepro);
			prepro = NULL;
		}
	} else {
		free(prepro);
		prepro = NULL;
	}
	if (prepro) {
		char *msg = siril_log_message(_("Preprocessing is ready\n"));
		display(msg);
	} else {
		char *msg = siril_log_message(_("Preprocessing not used\n"));
		display(msg);
	}
}

static int preprocess_image(char *filename, char *target) {
	if (!prepro || !prepro->use_dark) return 1;

	int ret = 0;
	fits fit = { 0 };
	if (readfits(filename, &fit, NULL, FALSE)) {
		return 1;
	}
	struct generic_seq_args generic = { .user = prepro };
	ret = prepro_image_hook(&generic, 0, 0, &fit, NULL);
	if (!ret)
		ret = savefits(target, &fit);
	clearfits(&fit);
	if (ret) {
		char *msg = siril_log_message(_("preprocessing failed\n"));
		display(msg);
	}
	return ret;
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
		/* Preprocess image */
		if (!preprocess_image(filename, target)) {
			siril_log_message(_("Preprocessed image to %s\n"), target);
			g_free(filename);
			filename = target;
			target = NULL;
		}

		/* TODO: Demosaic image */

		if (target && symlink_uniq_file(filename, target, do_links)) {
			display(_("Failed to rename or make a symbolic link to the input file"));
			break;
		}
		g_free(filename);
		g_free(target);

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
				display(_("Waiting for second image"));
				continue;
			}
			first_loop = FALSE;
		} else {
			seq.number = 2;
			seq.selnum = 2;
			seq.imgparam = malloc(2 * sizeof(imgdata));
			seq.imgparam[0].filenum = 1;
			seq.imgparam[0].incl = TRUE;
			seq.imgparam[0].date_obs = NULL;
			seq.imgparam[1].filenum = index;
			seq.imgparam[1].incl = TRUE;
			seq.imgparam[1].date_obs = NULL;
			writeseqfile(&seq);
		}
		
		if (seq_check_basic_data(&seq, FALSE) < 0) {
			display(_("Failed to read the sequence, aborting."));
			break;
		}
		if (seq_rx <= 0) {
			seq_rx = seq.rx;
			seq_ry = seq.ry;
			if (prepro && prepro->dark && (prepro->dark->rx != seq_rx || prepro->dark->ry != seq_ry)) {
				char *msg = siril_log_color_message(_("Dark image is not the same size, not using (%dx%d)\n"), "salmon", prepro->dark->rx, prepro->dark->ry);
				display(msg);
				clearfits(prepro->dark);
				free(prepro);
				prepro = NULL;
			}
		} else {
			if (seq_rx != seq.rx || seq_ry != seq.ry) {
				char *msg = siril_log_color_message(_("Images must have same dimensions.\n"), "red");
				display(msg);
				break;
			}
		}

		if (start_global_registration(&seq))
			continue;

		gchar *result_filename = g_strdup_printf("live_stack_00001%s", com.pref.ext);

		free_sequence(&seq, FALSE);

		/* Stack the sequence */
		siril_debug_print("Stacking image %d\n", index);

		/*sequence *r_seq = readseqfile("r_live_stack_.seq");
		if (!r_seq || seq_check_basic_data(r_seq, FALSE) < 0) {
			free(r_seq);
			return FALSE;
		}*/
		sequence r_seq;
		initialize_sequence(&r_seq, FALSE);
		create_seq_of_2(&r_seq, "r_live_stack_", index);
		if (seq_check_basic_data(&r_seq, FALSE) < 0) {
			//free(r_seq);
			continue;
		}

		struct stacking_args stackparam = { 0 };
		stackparam.method = stack_summing_generic; //stack_mean_with_rejection;
		stackparam.seq = &r_seq;
		stackparam.ref_image = 0;
		// TODO Apply some filter on quality?
		stackparam.filtering_criterion = seq_filter_all;
		stackparam.filtering_parameter = 0.0;
		stackparam.nb_images_to_stack = 2;
		stack_fill_list_of_unfiltered_images(&stackparam);
		stackparam.description = describe_filter(&r_seq, stackparam.filtering_criterion, stackparam.filtering_parameter);
		stackparam.output_filename = result_filename;
		stackparam.output_overwrite = TRUE;
		gettimeofday(&stackparam.t_start, NULL);

		stackparam.force_norm = FALSE;
		stackparam.output_norm = FALSE;
		stackparam.use_32bit_output = FALSE;
		stackparam.reglayer = (r_seq.nb_layers == 3) ? 1 : 0;

		main_stack(&stackparam);

		int retval = stackparam.retval;
		clean_end_stacking(&stackparam);
		free_sequence(&r_seq, FALSE);
		free(stackparam.image_indices);
		free(stackparam.description);

		if (retval) {
			gchar *str = g_strdup_printf(_("Stacking failed for image %d"), index);
			display(str);
			g_free(str);
			break;
		}

		//struct noise_data noise_args = { .fit = &gfit, .verbose = FALSE, .use_idle = FALSE };
		//noise(&noise_args);
		if (savefits(result_filename, &gfit)) {
			char *msg = siril_log_color_message(_("Could not save the stacking result %s, aborting\n"),
					"red", result_filename);
			display(msg);
			break;
		}

		/* Update display */
		queue_redraw(REMAP_ALL);
		gchar *str = g_strdup_printf(_("Stacked image %d"), index);
		display(str);
		g_free(str);

		index++;
	} while (1);

	siril_debug_print("===== exiting live stacking thread =====\n");

	// TODO: clean exit
	g_async_queue_unref(new_files_queue);
	return NULL;
}

