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

#include <sys/time.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/preprocess.h"
#include "gui/utils.h"
#include "gui/PSF_list.h"	// clear_stars_list
#include "io/conversion.h"
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
#include "algos/noise.h"
#include "algos/statistics.h"
#include "algos/demosaicing.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui.h"

/* hard-coded configuration */
#define REGISTRATION_TYPE SIMILARITY_TRANSFORMATION
#define REGISTRATION_INTERPOLATION OPENCV_AREA
/****************************/


#define WAIT_FILE_WRITTEN_US 60000	// check file size or readbility every 60 ms
#define WAIT_FILE_WRITTEN_ITERS 11	// 60*11 = 660ms, after that it's a failure
#define EXIT_TOKEN ":EXIT:"

static GFileMonitor *dirmon = NULL;
static gboolean do_links = FALSE;
static GThread *live_stacker_thread = NULL;
static GAsyncQueue *new_files_queue = NULL;

static gboolean paused = FALSE;

//static int registration_layer = 0;
//static fitted_PSF **ref_stars = NULL;
//static int nb_ref_stars = 0;
//static sequence *registered_seq = NULL;
static int seq_rx = -1, seq_ry = -1;
static struct star_align_data *sadata = NULL;
static regdata *regparam_bkp = NULL;
static struct preprocessing_data *prepro = NULL;
static gboolean first_stacking_result = TRUE;
static imstats *refimage_stats[3] = { NULL };

typedef enum {
	BOOL_NOT_SET,
	BOOL_TRUE,
	BOOL_FALSE
} super_bool;
static super_bool use_demosaicing = BOOL_NOT_SET;

static void init_preprocessing();
static gpointer live_stacker(gpointer arg);
int star_align_prepare_hook(struct generic_seq_args *args);

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

void pause_live_stacking_engine() {
	paused = !paused;
}

void stop_live_stacking_engine() {
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
		free(prepro);
		prepro = NULL;
	}
	if (sadata) {
		free_fitted_stars(sadata->refstars);
		free(sadata);
		sadata = NULL;
	}
	if (regparam_bkp) {
		free(regparam_bkp);
		regparam_bkp = NULL;
	}

	for (int i = 0; i < 3; i++) {
		if (refimage_stats[i]) {
			free_stats(refimage_stats[i]);
			refimage_stats[i] = NULL;
		}
	}
	seq_rx = -1; seq_ry = -1;
	use_demosaicing = BOOL_NOT_SET;
	paused = FALSE;
	first_stacking_result = TRUE;

	show_hide_toolbox();
	set_cursor_waiting(FALSE);
	unreserve_thread();
}

static int wait_for_file_to_be_written(const gchar *filename) {
	int iter;
	guint64 last_size = 0;
	GFile *fd = g_file_new_for_path(filename);
	for (iter = 1; iter < WAIT_FILE_WRITTEN_ITERS; iter++) {
		g_usleep(WAIT_FILE_WRITTEN_US);
#ifdef _WIN32
		GFileInputStream *stream;
		if ((stream = g_file_read(fd, NULL, NULL))) {
			g_object_unref(stream);
			break;
		}
#else
		guint64 size;
		if (!g_file_measure_disk_usage(fd, G_FILE_MEASURE_NONE, NULL, NULL, NULL, &size, NULL, NULL, NULL)) {
			g_object_unref(fd);
			return 1;
		}
		siril_debug_print("image size: %d MB\n", (int)(size/1000000));
		if (last_size == 0 || size != last_size)
			last_size = size;
		else break;
#endif
	}
	g_object_unref(fd);
	return iter >= WAIT_FILE_WRITTEN_ITERS;
}

static void file_changed(GFileMonitor *monitor, GFile *file, GFile *other,
		GFileMonitorEvent evtype, gpointer user_data) {
	if (evtype == G_FILE_MONITOR_EVENT_CREATED) {
		gchar *filename = g_file_get_basename(file);
		siril_debug_print("File %s created\n", filename);
		if (paused) return;	// manage in https://gitlab.com/free-astro/siril/-/issues/786

		image_type type;
		stat_file(filename, &type, NULL);
		if (type != TYPEFITS) {
			siril_log_message(_("File not supported for live stacking: %s\n"), filename);
			g_free(filename);
		} else {
			if (strncmp(filename, "live_stack", 10) &&
					strncmp(filename, "r_live_stack", 12) &&
					strncmp(filename, "result_live_stack", 17)) {
				if (wait_for_file_to_be_written(filename)) {
					gchar *str = g_strdup_printf(_("Could not open file: %s"), filename);
					livestacking_display(str, TRUE);
					return;
				}
				g_async_queue_push(new_files_queue, filename);
			}
			else g_free(filename);
		}
	}
}

/* for fullscreen, see livestacking_action_activate() in core/siril_actions.c */
void on_livestacking_start() {
	if (live_stacker_thread)
		return;
	livestacking_display(_("Starting live stacking"), FALSE);
	gui.rendering_mode = STF_DISPLAY;
	set_display_mode();
	force_unlinked_channels();
	show_hide_toolbox();

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
	livestacking_display_config(prepro && prepro->use_dark, prepro && prepro->use_flat, REGISTRATION_TYPE);

	g_signal_connect(G_OBJECT(dirmon), "changed", G_CALLBACK(file_changed), NULL);
	livestacking_display(_("Live stacking waiting for files"), FALSE);

	do_links = test_if_symlink_is_ok();
	reserve_thread();	// hack: generic function fails otherwise

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
			livestacking_display(_("not enough stars found in image, adjust settings\n"));
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

		if (cvTransformImage(fit, H, FALSE, REGISTRATION_INTERPOLATION)) {
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
	struct registration_args *regargs = sadata->regargs;
	regargs->imgparam = calloc(args->nb_filtered_images, sizeof(imgdata));
	regargs->regparam = calloc(args->nb_filtered_images, sizeof(regdata));
	if (!regargs->imgparam  || !regargs->regparam) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	args->seq->regparam[regargs->layer] = regparam_bkp;
	sadata->success[0] = 0;
	sadata->success[1] = 0;
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
	reg_args.interpolation = REGISTRATION_INTERPOLATION;
	reg_args.type = REGISTRATION_TYPE;
	reg_args.max_stars_candidates = 444;
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

	int retval = GPOINTER_TO_INT(generic_sequence_worker(args));

	// hack to not free the regparam, because it's referenced by
	// sadata->current_regdata because of the first call to the prepare,
	// and used for reference frame params in the registration
	regparam_bkp = seq->regparam[reg_args.layer];
	seq->regparam[reg_args.layer] = NULL;
	free_sequence(seq, FALSE);

	return retval || !sadata->success[1];
}

static void init_preprocessing() {
	/* copied from core/preprocess.c test_for_master_files() but modified to not check
	 * for some options and not check for image properties against gfit, which is not
	 * already loaded here. Error management is different too.
	 * Some options also come from on_prepro_button_clicked()
	 */
	prepro = calloc(1, sizeof(struct preprocessing_data));

	/* checking for dark master and associated options */
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
				livestacking_display(_("NOT USING DARK: cannot open the file"), FALSE);
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

			GtkToggleButton *fix_xtrans = GTK_TOGGLE_BUTTON(lookup_widget("fix_xtrans_af"));
			prepro->fix_xtrans = gtk_toggle_button_get_active(fix_xtrans);
		}
	}

	/* checking for flat master and associated options */
	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		GtkEntry *entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
			prepro->flat = calloc(1, sizeof(fits));
			if (!readfits(filename, prepro->flat, NULL, !com.pref.force_to_16bit)) {
				prepro->use_flat = TRUE;
			} else {
				livestacking_display(_("NOT USING FLAT: cannot open the file"), FALSE);
				free(prepro->flat);
				gtk_entry_set_text(entry, "");
				prepro->use_flat = FALSE;
			}

			if (prepro->use_flat) {
				GtkToggleButton *autobutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto_evaluate"));
				prepro->autolevel = gtk_toggle_button_get_active(autobutton);
				if (!prepro->autolevel) {
					GtkEntry *norm_entry = GTK_ENTRY(lookup_widget("entry_flat_norm"));
					prepro->normalisation = g_ascii_strtod(gtk_entry_get_text(norm_entry), NULL);
				}


				GtkToggleButton *CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
				prepro->is_cfa = gtk_toggle_button_get_active(CFA);

				GtkToggleButton *equalize_cfa = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
				prepro->equalize_cfa = gtk_toggle_button_get_active(equalize_cfa);
			}
		}
	}
	prepro->is_sequence = FALSE;
	prepro->allow_32bit_output = !com.pref.force_to_16bit;

	if (prepro->use_dark || prepro->use_flat) {
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
		livestacking_display(msg, FALSE);
	} else {
		char *msg = siril_log_message(_("Preprocessing not used\n"));
		livestacking_display(msg, FALSE);
	}
}

static int preprocess_image(char *filename, char *target) {
	if (!prepro || (!prepro->use_dark && !prepro->use_flat)) return 1;

	int ret = 0;
	fits fit = { 0 };
	if (readfits(filename, &fit, NULL, FALSE)) {
		return 1;
	}
	struct generic_seq_args generic = { .user = prepro };
	ret = prepro_image_hook(&generic, 0, 0, &fit, NULL, com.max_thread);
	if (!ret)
		ret = savefits(target, &fit);
	clearfits(&fit);
	if (ret) {
		char *msg = siril_log_message(_("preprocessing failed\n"));
		livestacking_display(msg, FALSE);
	}
	return ret;
}

static gpointer live_stacker(gpointer arg) {
	g_async_queue_ref(new_files_queue);
	int index = 1, number_of_images_stacked = 1;
	gboolean first_loop = TRUE;	// only original images in the sequence
	do {
		gchar *filename = g_async_queue_pop(new_files_queue); // blocking
		if (!strcmp(filename, EXIT_TOKEN)) {
			siril_debug_print("Exiting thread\n");
			break;
		}

		/* init demosaicing (check if the incoming file is CFA) */
		if (use_demosaicing == BOOL_NOT_SET) {	// another kind of first_loop
			fits fit = { 0 };
			if (read_fits_metadata_from_path(filename, &fit)) {
				livestacking_display(_("Failed to open the first image\n"), FALSE);
				clearfits(&fit);
				break;
			}
			gboolean is_CFA = fit.bayer_pattern[0] != '\0';
			use_demosaicing = is_CFA ? BOOL_TRUE : BOOL_FALSE;
			if (prepro)
				prepro->debayer = is_CFA;

			enable_debayer(is_CFA);
			clearfits(&fit);
			update_debayer_button_status(is_CFA);
			if (is_CFA)
				siril_log_message(_("live stacking will debayer images\n"));
			else siril_log_message(_("live stacking will not debayer images\n"));

			/* we want to load the image while avoiding GTK+ calls
			 * from this thread, setting com.script is a hack to
			 * avoid this, but this doesn't display the loaded
			 * image, so we still need the extra idle in
			 * complete_image_loading
			 */
			com.script = TRUE;
			open_single_image(filename);
			com.script = FALSE;
			complete_image_loading();
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
		else if (use_demosaicing == BOOL_TRUE) {
			fits fit = { 0 };
			int retval = readfits(filename, &fit, NULL, !com.pref.force_to_16bit);
			if (!retval)
				retval = debayer_if_needed(TYPEFITS, &fit, TRUE);
			if (!retval)
				retval = savefits(target, &fit);
			if (!retval) {
				g_free(filename);
				filename = target;
				target = NULL;
			}
			clearfits(&fit);
		}

		if (target && symlink_uniq_file(filename, target, do_links)) {
			livestacking_display(_("Failed to rename or make a symbolic link to the input file"), FALSE);
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
				livestacking_display(_("Waiting for second image"), FALSE);
				livestacking_update_number_of_images(1, gfit.exposure, -1.0);
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
			livestacking_display(_("Failed to read the sequence, aborting."), FALSE);
			break;
		}
		if (seq_rx <= 0) {
			seq_rx = seq.rx;
			seq_ry = seq.ry;
			if (prepro && prepro->dark && (prepro->dark->rx != seq_rx || prepro->dark->ry != seq_ry)) {
				char *msg = siril_log_color_message(_("Dark image is not the same size, not using (%dx%d)\n"), "salmon", prepro->dark->rx, prepro->dark->ry);
				livestacking_display(msg, FALSE);
				clearfits(prepro->dark);
				prepro->use_dark = FALSE;
			}
			if (prepro && prepro->flat && (prepro->flat->rx != seq_rx || prepro->flat->ry != seq_ry)) {
				char *msg = siril_log_color_message(_("Flat image is not the same size, not using (%dx%d)\n"), "salmon", prepro->flat->rx, prepro->flat->ry);
				livestacking_display(msg, FALSE);
				clearfits(prepro->flat);
				prepro->use_flat = FALSE;
			}
			if (prepro && !prepro->dark && !prepro->flat) {
				free(prepro);
				prepro = NULL;
			}
		} else {
			if (seq_rx != seq.rx || seq_ry != seq.ry) {
				char *msg = siril_log_color_message(_("Images must have same dimensions.\n"), "red");
				livestacking_display(msg, FALSE);
				break;
			}
		}

		if (start_global_registration(&seq))
			continue;

		gchar *result_filename = g_strdup_printf("live_stack_00001%s", com.pref.ext);

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
			continue;
		}
		if (refimage_stats[0])
			for (int i = 0; i < r_seq.nb_layers; i++)
				add_stats_to_seq(&r_seq, 0, i, refimage_stats[i]);

		struct stacking_args stackparam = { 0 };
		stackparam.method = stack_mean_with_rejection;
		stackparam.seq = &r_seq;
		stackparam.ref_image = 0;
		// TODO Apply some filter on quality?
		stackparam.filtering_criterion = seq_filter_all;
		stackparam.filtering_parameter = 0.0;
		stackparam.nb_images_to_stack = 2;
		stack_fill_list_of_unfiltered_images(&stackparam);
		stackparam.description = describe_filter(&r_seq, stackparam.filtering_criterion, stackparam.filtering_parameter);
		stackparam.output_filename = result_filename; // not used with main_stack()
		stackparam.output_overwrite = TRUE;
		gettimeofday(&stackparam.t_start, NULL);

		stackparam.normalize = ADDITIVE_SCALING;	// TODO: toggle switch
		stackparam.force_norm = FALSE;
		stackparam.output_norm = FALSE;
		stackparam.equalizeRGB = FALSE;		// not possible currently
		stackparam.lite_norm = TRUE;
		stackparam.use_32bit_output = !com.pref.force_to_16bit;
		stackparam.reglayer = (r_seq.nb_layers == 3) ? 1 : 0;
		stackparam.apply_nbstack_weights = TRUE;

		main_stack(&stackparam);

		/* doing things similar to end_stacking */
		int retval = stackparam.retval;
		/* and hacking the stats for good normalization: the reference is the first image stacked */
		if (!retval && !refimage_stats[0]) {
			if (copy_cached_stats_for_image(&r_seq, 0, refimage_stats)) {
				siril_log_color_message(_("Reference image statistics not found\n"), "red");
				stackparam.normalize = NO_NORM;
			}
			else siril_debug_print("saved statistics of reference image, using normalization\n");
		}
		clean_end_stacking(&stackparam);
		free(stackparam.image_indices);
		free(stackparam.description);
		free_sequence(&r_seq, FALSE);

		if (retval) {
			gchar *str = g_strdup_printf(_("Stacking failed for image %d"), index);
			livestacking_display(str, TRUE);
			break;
		}
		clear_stars_list(FALSE);
		bgnoise_async(&gfit, FALSE);

		if (savefits(result_filename, &gfit)) {
			char *msg = siril_log_color_message(_("Could not save the stacking result %s, aborting\n"),
					"red", result_filename);
			livestacking_display(msg, FALSE);
			break;
		}

		/* Update display */
		if (first_stacking_result) {
			/* number of channels may have changed */
			com.seq.current = RESULT_IMAGE;
			com.uniq->nb_layers = gfit.naxes[2];
			com.uniq->fit = &gfit;
			siril_add_idle(livestacking_first_result_idle, NULL);
			first_stacking_result = FALSE;
		} else {
			queue_redraw(REMAP_ALL);
		}
		g_free(result_filename);
		gchar *str = g_strdup_printf(_("Stacked image %d"), index);
		livestacking_display(str, TRUE);

		index++;
		number_of_images_stacked++;
		double noise = bgnoise_await();
		livestacking_update_number_of_images(number_of_images_stacked, gfit.livetime, noise);
	} while (1);

	siril_debug_print("===== exiting live stacking thread =====\n");

	// TODO: clean exit
	g_async_queue_unref(new_files_queue);
	return NULL;
}

int get_paused_status() {
	return paused ? 1 : 0;
}

