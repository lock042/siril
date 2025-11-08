/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#ifdef _WIN32
#include <windows.h>
#endif

#include <sys/time.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/preprocess.h"
#include "core/siril_log.h"
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
#include "opencv/opencv.h"
/* ******************* */
#include "stacking/stacking.h"
#include "algos/noise.h"
#include "algos/statistics.h"
#include "algos/demosaicing.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "livestacking.h"
#include "gui.h"

/* hard-coded configuration */
#define REGISTRATION_INTERPOLATION OPENCV_AREA
/****************************/


#define WAIT_FILE_WRITTEN_US 80000	// check file size or readbility every 80 ms
#define WAIT_FILE_WRITTEN_ITERS 15	// .080*15 = 1.2s, after that it's a failure
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

/* config */
static gboolean use_32bits = FALSE;
static super_bool use_demosaicing = BOOL_NOT_SET;
static transformation_type reg_type = HOMOGRAPHY_TRANSFORMATION;
static gboolean reg_rotates = TRUE;

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
	seq->fz = com.pref.comp.fits_enabled;
}

void pause_live_stacking_engine() {
	paused = !paused;
}

void stop_live_stacking_engine() {
	siril_log_message(_("Stopping live stacking engine...\n"));
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
		new_files_queue = NULL;
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
	unreserve_thread();

	if (!com.headless) {
		GtkWidget *toolbar = lookup_widget("GtkToolMainBar");
		if (!gtk_widget_is_visible(toolbar)) show_hide_toolbox();
		set_cursor_waiting(FALSE);
	}
}

static int wait_for_file_to_be_written(const gchar *filename) {
	int iter;
#ifndef _WIN32
	guint64 last_size = 0;
#endif
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
		siril_debug_print("image size: %d MB\n", (int )(size / 1000000));
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
	if (evtype != G_FILE_MONITOR_EVENT_CREATED && evtype != G_FILE_MONITOR_EVENT_MOVED_IN) {
		return;
	}
	gchar *filename = g_file_get_basename(file);
	siril_debug_print("File %s added\n", filename);
	if (filename[0] == '.' || // hidden files
			paused)	{ // manage in https://gitlab.com/free-astro/siril/-/issues/786
		g_free(filename);
		return;
	}

	image_type type;
	if (stat_file(filename, &type, NULL)) {
		siril_debug_print("Filename is not canonical\n");
	}
	if (type != TYPEFITS) {
		if (type == TYPERAW) {
			if (!wait_for_file_to_be_written(filename)) {
				fits dest = { 0 };
				gchar *new = replace_ext(filename, com.pref.ext);
				any_to_fits(TYPERAW, filename, &dest, FALSE, !com.pref.force_16bit);
				savefits(new, &dest);
				clearfits(&dest);
			}
		}  else {
			siril_log_message(_("File not supported for live stacking: %s\n"), filename);
		}
		g_free(filename);
	} else {
		if (strncmp(filename, "live_stack", 10) &&
				strncmp(filename, "r_live_stack", 12) &&
				strncmp(filename, "result_live_stack", 17)) {
			if (wait_for_file_to_be_written(filename)) {
				gchar *str = g_strdup_printf(_("Could not open file: %s"), filename);
				livestacking_display(str, TRUE);
				g_free(filename);
				return;
			}
			g_async_queue_push(new_files_queue, filename);
		}
		else g_free(filename);
	}
}

void livestacking_queue_file(char *file) {
	g_async_queue_push(new_files_queue, file);
}

/* for fullscreen, see livestacking_action_activate() in core/siril_actions.c */
int start_livestacking(gboolean with_filewatcher) {
	if (live_stacker_thread)
		return 1;
	livestacking_display(_("Starting live stacking"), FALSE);
	if (!com.headless) {
		gui.rendering_mode = STF_DISPLAY;
		set_display_mode();
		force_unlinked_channels();
		GtkWidget *toolbar = lookup_widget("GtkToolMainBar");
		if (gtk_widget_is_visible(toolbar)) show_hide_toolbox();
		livestacking_display_config(prepro && prepro->use_dark, prepro && prepro->use_flat, reg_type);
	}

	do_links = test_if_symlink_is_ok(TRUE);

	new_files_queue = g_async_queue_new();
	if (with_filewatcher) {
		/* start monitoring CWD */
		GFile *cwd = g_file_new_for_path(com.wd);
		GError *err = NULL;
		dirmon = g_file_monitor_directory(cwd, G_FILE_MONITOR_WATCH_MOVES, NULL, &err);
		g_object_unref(cwd);
		if (err) {
			siril_log_message(_("Unable to monitor CWD (%s): %s\n"), com.wd, err->message);
			g_async_queue_unref(new_files_queue);
			new_files_queue = NULL;
			dirmon = NULL;
			g_error_free(err);
			return 1;
		}

		if (g_signal_connect(G_OBJECT(dirmon), "changed", G_CALLBACK(file_changed), NULL) <= 0) {
			siril_log_message(_("Unable to monitor CWD (%s): %s\n"), com.wd, "signal did not connect");
			new_files_queue = NULL;
			g_async_queue_unref(new_files_queue);
			g_object_unref(dirmon);
			dirmon = NULL;
			return 1;
		}

		siril_debug_print("file watcher active for CWD (%s)\n", com.wd);
	}

	live_stacker_thread = g_thread_new("live stacker", live_stacker, NULL);
	return 0;
}

static void init_preprocessing_from_command(char *dark, char *flat, gboolean use_32bits) {
	prepro = calloc(1, sizeof(struct preprocessing_data));
	if (dark) {
		prepro->dark = calloc(1, sizeof(fits));
		if (readfits(dark, prepro->dark, NULL, FALSE)) {
			siril_log_message(_("NOT USING DARK: cannot open file '%s'\n"), dark);
			free(prepro->dark);
			prepro->use_dark = FALSE;
			prepro->use_cosmetic_correction = FALSE;
		} else {
			prepro->use_dark = TRUE;
			if (prepro->dark->naxes[2] > 1) {
				clearfits(prepro->dark);
				free(prepro->dark);
				prepro->use_dark = FALSE;
				prepro->use_cosmetic_correction = FALSE;
				siril_log_message(_("Calibration with color images is not yet supported\n"));
			}
			else if (strlen(prepro->dark->keywords.bayer_pattern) > 4) {
				prepro->use_cosmetic_correction = FALSE;
				prepro->fix_xtrans = TRUE;
			} else {
				prepro->use_cosmetic_correction = TRUE;
				prepro->sigma[0] = -1.0;
				prepro->sigma[1] = 3.5;
				if (strlen(prepro->dark->keywords.bayer_pattern) >= 4)
					prepro->is_cfa = TRUE;
			}
			siril_log_message(_("Master dark %d x %d configured for live stacking (%s cosmetic correction)\n"), prepro->dark->rx, prepro->dark->ry, prepro->use_cosmetic_correction ? _("with") : _("without") );
		}
	}
	if (flat) {
		prepro->flat = calloc(1, sizeof(fits));
		if (readfits(flat, prepro->flat, NULL, FALSE)) {
			siril_log_message(_("NOT USING FLAT: cannot open file '%s'\n"), flat);
			free(prepro->flat);
			prepro->use_flat = FALSE;
		}
		else {
			if (prepro->flat->naxes[2] > 1) {
				clearfits(prepro->flat);
				free(prepro->flat);
				prepro->use_flat = FALSE;
				siril_log_message(_("Calibration with color images is not yet supported\n"));
			} else {
				prepro->use_flat = TRUE;
				prepro->autolevel = TRUE;
				if (strlen(prepro->flat->keywords.bayer_pattern) >= 4) {
					prepro->is_cfa = TRUE;
					prepro->equalize_cfa = TRUE;
				}
				siril_log_message(_("Master flat %d x %d configured for live stacking\n"), prepro->flat->rx, prepro->flat->ry);
			}
		}
	}

	init_preprocessing_finalize(prepro, use_32bits);
}

/* code common to GUI and CLI for prepro init */
void init_preprocessing_finalize(struct preprocessing_data *prepro_data, gboolean use_32b) {
	prepro = prepro_data;
	prepro->is_sequence = FALSE;
	use_32bits = use_32b;
	prepro->allow_32bit_output = use_32b;

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

	char *msg = prepro ?
		siril_log_message(_("Preprocessing is ready\n")) :
		siril_log_message(_("Preprocessing not used\n"));
	if (!com.headless)
		livestacking_display(msg, FALSE);
}

void init_registration_finalize(gboolean shift_only) {
	reg_type = shift_only ? SHIFT_TRANSFORMATION : SIMILARITY_TRANSFORMATION;
	reg_rotates = !shift_only;
}

int start_livestack_from_command(gchar *dark, gchar *flat, gboolean use_file_watcher/*, gboolean remove_gradient*/, gboolean shift_only, gboolean use_32b) {
	if (live_stacker_thread) {
		siril_log_message(_("live stacking is already running, stop it first\n"));
		return 1;
	}

	prepro = calloc(1, sizeof(struct preprocessing_data));
	init_preprocessing_from_command(dark, flat, use_32b);

	init_registration_finalize(shift_only);

	int retval = start_livestacking(use_file_watcher);
	if (retval && prepro) {
		clear_preprocessing_data(prepro);
		free(prepro);
		prepro = NULL;
	}
	return retval;
}

/*static int register_image(fits *fit, char *output_name) {
	fitted_PSF **stars;
	int nb_stars;

	stars = peaker(fit, registration_layer, &com.starfinder_conf, &nb_stars, NULL, FALSE, TRUE, MAXSTARS, com.pref.starfinder_conf.profile, com.max_thread);
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

		if (cvTransformImage(fit, H, FALSE, REGISTRATION_INTERPOLATION, TRUE, NULL)) {
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
	struct registration_args regargs = { 0 };
	regargs.func = &register_star_alignment;
	regargs.seq = seq;
	regargs.reference_image = 0;
	regargs.layer = (seq->nb_layers == 3) ? 1 : 0;
	regargs.run_in_thread = FALSE;
	regargs.follow_star = FALSE;
	regargs.matchSelection = FALSE;
	//memcpy(&regargs.selection, &com.selection, sizeof(rectangle));
	regargs.output_scale = 1.f;
	regargs.min_pairs = 10;
	regargs.no_output = !reg_rotates;
	regargs.prefix = "r_";
	regargs.load_new_sequence = FALSE;
	regargs.interpolation = REGISTRATION_INTERPOLATION;
	regargs.type = reg_type;
	regargs.max_stars_candidates = 200;
	cvGetEye(&regargs.framingd.Htransf);
	cvGetEye(&regargs.framingd.Hshift);
	regargs.framingd.roi_out = (framing_roi){ 0, 0, seq->rx, seq->ry};

	// preparing detection params
	regargs.sfargs = calloc(1, sizeof(struct starfinder_data));
	regargs.sfargs->im.from_seq = regargs.seq;
	regargs.sfargs->layer = regargs.layer;
	regargs.sfargs->keep_stars = TRUE;
	regargs.sfargs->save_to_file = FALSE;
	regargs.sfargs->max_stars_fitted = regargs.max_stars_candidates;

	struct generic_seq_args *args = create_default_seqargs(seq);
	if (regargs.filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = seq->selnum;
	}
	//args->compute_mem_limits_hook = star_align_compute_mem_limits;
	args->prepare_hook = live_stacking_star_align_prepare;
	args->image_hook = star_align_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Global star registration");
	args->has_output = reg_rotates;
	args->output_type = get_data_type(seq->bitpix);
	args->upscale_ratio = 1.0;
	args->new_seq_prefix = strdup(regargs.prefix);
	args->load_new_sequence = FALSE;
	args->already_in_a_thread = TRUE;
	if (!sadata) {
		sadata = calloc(1, sizeof(struct star_align_data));
		sadata->regargs = &regargs;
	}
	args->user = sadata;
	/* registration will work in 16 bits if input data is 16 bits */

	reserve_thread(); // hack: generic function fails otherwise
	int retval = GPOINTER_TO_INT(generic_sequence_worker(args));

	// hack to not free the regparam, because it's referenced by
	// sadata->current_regdata because of the first call to the prepare,
	// and used for reference frame params in the registration
	regparam_bkp = seq->regparam[regargs.layer];
	seq->regparam[regargs.layer] = NULL;
	free_sequence(seq, FALSE);
	free(regargs.sfargs);
	return retval || !sadata->success[1];
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
		msg[strlen(msg) - 1] = '\0';
		livestacking_display(msg, FALSE);
	}
	return ret;
}

static gpointer live_stacker(gpointer arg) {
	g_async_queue_ref(new_files_queue);
	int index = 1, number_of_images_stacked = 1;
	livestacking_display(_("Live stacking waiting for files"), FALSE);
	gboolean first_loop = TRUE;	// only original images in the sequence
	do {
		gchar *filename = g_async_queue_pop(new_files_queue); // blocking
		if (!strcmp(filename, EXIT_TOKEN)) {
			siril_debug_print("Exiting thread\n");
			break;
		}
		struct timeval tv_start, tv_tmp, tv_end;
		gettimeofday(&tv_start, NULL);

		/* init demosaicing (check if the incoming file is CFA) */
		if (use_demosaicing == BOOL_NOT_SET) {	// another kind of first_loop
			fits fit = { 0 };
			if (read_fits_metadata_from_path(filename, &fit)) {
				livestacking_display(_("Failed to open the first image"), FALSE);
				clearfits(&fit);
				break;
			}
			gboolean is_CFA = fit.keywords.bayer_pattern[0] != '\0';
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
			 * image, so we still need the extra idle called by
			 * execute_idle_and_wait_for_it() (this immediately
			 * returns if headless)
			 */
			gboolean script_bkp = com.script;
			com.script = TRUE;
			open_single_image(filename);
			com.script = script_bkp;
			if (!com.headless)
				execute_idle_and_wait_for_it(end_image_loading, NULL);
		}

		siril_debug_print("Adding file to input sequence\n");
		gchar *target = g_strdup_printf("live_stack_%05d%s", index, get_com_ext(com.pref.comp.fits_enabled));
		/* Preprocess image */
		if (!preprocess_image(filename, target)) {
			siril_log_message(_("Preprocessed image to %s\n"), target);
			g_free(filename);
			filename = target;
			target = NULL;
		}
		else if (use_demosaicing == BOOL_TRUE) {
			fits fit = { 0 };
			int retval = readfits(filename, &fit, NULL, !com.pref.force_16bit);
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
		gettimeofday(&tv_end, NULL);
		show_time_msg(tv_start, tv_end, "calibration and demosaicing");
		tv_tmp = tv_end;

		if (target && symlink_uniq_file(filename, target, do_links)) {
			g_free(target);
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
		seq.fz = com.pref.comp.fits_enabled;
		if (first_loop) {
			if (buildseqfile(&seq, 1) || seq.number == 1) {
				index++;
				livestacking_display(_("Waiting for second image"), FALSE);
				livestacking_update_number_of_images(1, gfit->keywords.exposure, -1.0, NULL);
				free(seq.seqname);
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
				msg[strlen(msg) - 1] = '\0';
				livestacking_display(msg, FALSE);
				clearfits(prepro->dark);
				prepro->use_dark = FALSE;
			}
			if (prepro && prepro->flat && (prepro->flat->rx != seq_rx || prepro->flat->ry != seq_ry)) {
				char *msg = siril_log_color_message(_("Flat image is not the same size, not using (%dx%d)\n"), "salmon", prepro->flat->rx, prepro->flat->ry);
				msg[strlen(msg) - 1] = '\0';
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
				msg[strlen(msg) - 1] = '\0';
				livestacking_display(msg, FALSE);
				break;
			}
		}

		if (start_global_registration(&seq))
			continue;

		gettimeofday(&tv_end, NULL);
		show_time_msg(tv_tmp, tv_end, "registration");
		tv_tmp = tv_end;

		gchar *result_filename = g_strdup_printf("live_stack_00001%s", get_com_ext(com.pref.comp.fits_enabled));

		/* Stack the sequence */
		siril_debug_print("Stacking image %d\n", index);

		/*sequence *r_seq = readseqfile("r_live_stack_.seq");
		if (!r_seq || seq_check_basic_data(r_seq, FALSE) < 0) {
			free(r_seq);
			return FALSE;
		}*/
		sequence r_seq;	// registered sequence
		initialize_sequence(&r_seq, FALSE);
		create_seq_of_2(&r_seq, reg_rotates ? "r_live_stack_" : "live_stack_", index);
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
		/* we should not use 32 bits for stack results if input files are 16 bits,
		 * otherwise we end up with a mixed sequence to process;
		 * inputs to stack are 16 bits when no preprocessing or debayer occur */
		stackparam.use_32bit_output = get_data_type(r_seq.bitpix) == DATA_FLOAT ||
			(use_32bits && (prepro || use_demosaicing == BOOL_TRUE));
		stackparam.reglayer = (r_seq.nb_layers == 3) ? 1 : 0;
		stackparam.weighting_type = NBSTACK_WEIGHT;

		reserve_thread(); // hack: generic function fails otherwise
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
		free_sequence(&r_seq, FALSE);
		free(stackparam.image_indices);
		g_free(stackparam.description);

		if (retval) {
			gchar *str = g_strdup_printf(_("Stacking failed for image %d"), index);
			livestacking_display(str, TRUE);
			break;
		}
		clear_stars_list(FALSE);
		bgnoise_async(&stackparam.result, TRUE);

		if (savefits(result_filename, &stackparam.result)) {
			char *msg = siril_log_color_message(_("Could not save the stacking result %s, aborting\n"),
					"red", result_filename);
			msg[strlen(msg) - 1] = '\0';
			livestacking_display(msg, FALSE);
			bgnoise_await();
			break;
		}

		if (!com.headless) {
			/* Update display */
			clearfits(gfit);
			memcpy(gfit, &stackparam.result, sizeof(fits));
			if (first_stacking_result) {
				/* number of channels may have changed */
				com.seq.current = RESULT_IMAGE;
				com.uniq->nb_layers = gfit->naxes[2];
				com.uniq->fit = gfit;
				gdk_threads_add_idle(livestacking_first_result_idle, NULL);
				first_stacking_result = FALSE;
			} else {
				queue_redraw(REMAP_ALL); // TODO: is this safe enough if the livestacking is running from a python command?
			}
		}
		g_free(result_filename);
		gchar *str = g_strdup_printf(_("Stacked image %d"), index);
		livestacking_display(str, TRUE);

		index++;
		number_of_images_stacked++;
		double noise = bgnoise_await();
		if (com.headless)
			clearfits(&stackparam.result);
		gettimeofday(&tv_end, NULL);
		show_time_msg(tv_tmp, tv_end, "stacking");
		const char *total_time = format_time_diff(tv_start, tv_end);
		siril_log_color_message(_("Time to process the last image for live stacking: %s\n"),
				"green", total_time);
		livestacking_update_number_of_images(number_of_images_stacked, gfit->keywords.livetime, noise, total_time);
	} while (1);

	siril_debug_print("===== exiting live stacking thread =====\n");

	// TODO: clean exit
	g_async_queue_unref(new_files_queue);
	return NULL;
}

int get_paused_status() {
	return paused ? 1 : 0;
}

gboolean livestacking_is_started() {
	return live_stacker_thread != NULL;
}

gboolean livestacking_uses_filewatcher() {
	return dirmon != NULL;
}
