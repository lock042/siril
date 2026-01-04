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
#else
#include <sys/socket.h>
#include <sys/wait.h>
#endif

#include <assert.h>
#include <string.h>
#include <glib.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/sequence_filtering.h"
#include "core/command_line_processor.h"
#include "core/masks.h"
#include "core/OS_utils.h"
#include "core/undo.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "io/single_image.h"
#include "gui/progress_and_log.h"
#include "gui/script_menu.h"
#include "gui/utils.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/seqwriter.h"
#include "io/fits_sequence.h"
#include "io/image_format_fits.h"
#include "algos/statistics.h"
#include "registration/registration.h"

static GMutex child_mutex = { 0 };

void child_mutex_lock() {
	g_mutex_lock(&child_mutex);
}

void child_mutex_unlock() {
	g_mutex_unlock(&child_mutex);
}

/* destroy_any_args() will destroy any operation-specific arguments struct that
 * contains a destructor as its first member.
 * The idea is that all structs passed to generic_seq_args or generic_img_args as
 * user_data fit this polymorphic pattern. If the destroy_fn is NULL we assume
 * there are no dynamically allocated members and the struct can just be freed.
 */
static void destroy_any_args(void *obj) {
    destructor *destroy_fn = obj; // first field assumed destroy
    if (*destroy_fn)
		(*destroy_fn)(obj);
	else
		free(obj);
}

// called in start_in_new_thread only
// works in parallel if the arg->parallel is TRUE for FITS or SER sequences
gpointer generic_sequence_worker(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *)p;
	struct timeval t_start, t_end;
	int frame;	// output frame index
	int input_idx;	// index of the frame being processed in the sequence
	int *index_mapping = NULL;
	int nb_frames, excluded_frames = 0, progress = 0;
	int abort = 0;	// variable for breaking out of loop
#ifdef _OPENMP
	int* threads_per_image = NULL;
#endif
	gboolean have_seqwriter = FALSE;

	assert(args);
	assert(args->seq);
	assert(args->image_hook);
	set_progress_bar_data(NULL, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	if (args->nb_filtered_images > 0)
		nb_frames = args->nb_filtered_images;
	else {
		nb_frames = compute_nb_filtered_images(args->seq, args->filtering_criterion, args->filtering_parameter);
		args->nb_filtered_images = nb_frames;
		if (nb_frames <= 0) {
			siril_log_message(_("No image selected for processing, aborting\n"));
			args->retval = 1;
			goto the_end;
		}
	}
	float nb_framesf = (float)nb_frames + 0.3f;	// leave margin for rounding errors and post processing
	args->retval = 0;

#ifdef HAVE_FFMS2
	// If the sequence is of type SEQ_AVI, lock out the sequence browser as it can cause a crash
	if (args->seq->type == SEQ_AVI && !com.headless) {
		gui_function(set_seq_browser_active, GINT_TO_POINTER(0));
	}
#endif
#ifdef _OPENMP
	// max_parallel_images can be computed by another function too, hence the check
	if (args->max_parallel_images < 1) {
		if (args->compute_mem_limits_hook)
			args->max_parallel_images = args->compute_mem_limits_hook(args, FALSE);
		else args->max_parallel_images = seq_compute_mem_limits(args, FALSE);
	}
	if (args->max_parallel_images < 1) {
		args->retval = 1;
		goto the_end;
	}
	/* if there are less images in the sequence than threads, we still
	 * distribute them */
	if (args->max_parallel_images > nb_frames)
		args->max_parallel_images = nb_frames;

	siril_log_message(_("%s: with the current memory and thread limits, up to %d thread(s) can be used\n"),
			args->description, args->max_parallel_images);

	// remaining threads distribution per image thread
	threads_per_image = compute_thread_distribution(args->max_parallel_images, com.max_thread);
#endif

	if (args->prepare_hook && args->prepare_hook(args)) {
		siril_log_message(_("Preparing sequence processing failed.\n"));
		args->retval = 1;
		goto the_end;
	}

	/* We have a sequence in which images can be filtered out. In order to
	 * distribute the workload fairly among all threads, the main iteration
	 * should not be on the list of images of the sequence, but on the
	 * actual list of selected images.
	 * Here we create this map of images to be processed, each cell of the
	 * array providing the image number in the input sequence. It cannot be
	 * done in parallel.
	 * This is mandatory for SER contiguous output. */
	if (args->filtering_criterion) {
		index_mapping = malloc(nb_frames * sizeof(int));
		if (!index_mapping) {
			PRINT_ALLOC_ERR;
			args->retval = 1;
			goto the_end;
		}
		for (input_idx = 0, frame = 0; input_idx < args->seq->number; input_idx++) {
			if (!args->filtering_criterion(args->seq, input_idx, args->filtering_parameter)) {
				continue;
			}
			index_mapping[frame++] = input_idx;
		}
		if (frame != nb_frames) {
			siril_log_message(_("Output index mapping failed (%d/%d).\n"), frame, nb_frames);
			args->retval = 1;
			goto the_end;
		}
	}

	if (args->has_output && !args->partial_image && !(args->seq->type == SEQ_INTERNAL)) {
		gint64 size;
		if (args->compute_size_hook)
			size = args->compute_size_hook(args, nb_frames);
		else size = seq_compute_size(args->seq, nb_frames, args->output_type);
		if (test_available_space(size)) {
			siril_log_color_message(_("Not enough space to save the output images, aborting\n"), "red");
			args->retval = 1;
			goto the_end;
		}
	}

	/* Output print of algorithm description */
	if (args->description) {
		gchar *desc = g_strdup_printf(_("%s: processing...\n"), args->description);
		siril_log_color_message(desc, "green");
		g_free(desc);
	}

	have_seqwriter = args->has_output &&
		((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) ||
		 (args->force_ser_output || args->seq->type == SEQ_SER));
#ifdef _OPENMP
	omp_init_lock(&args->lock);
	if (have_seqwriter)
		omp_set_schedule(omp_sched_dynamic, 1);
	else omp_set_schedule(omp_sched_guided, 0);
#ifdef HAVE_FFMS2
	// we don't want to enable parallel processing for films, as ffms2 is not thread-safe
#pragma omp parallel for num_threads(args->max_parallel_images) private(input_idx) schedule(runtime) \
	if(args->max_parallel_images > 1 && args->seq->type != SEQ_AVI && args->parallel && (args->seq->type == SEQ_SER || fits_is_reentrant()))
#else
#pragma omp parallel for num_threads(args->max_parallel_images) private(input_idx) schedule(runtime) \
	if(args->max_parallel_images > 1 && args->parallel && (args->seq->type == SEQ_SER || fits_is_reentrant()))
#endif // HAVE_FFMS2
#endif // _OPENMP
	for (frame = 0; frame < nb_frames; frame++) {
		if (abort) continue;

		char filename[256];
		rectangle area = { .x = args->area.x, .y = args->area.y,
			.w = args->area.w, .h = args->area.h };

		if (!get_thread_run()) {
			abort = 1;
			continue;
		}
		if (index_mapping)
			input_idx = index_mapping[frame];
		else input_idx = frame;

		if (!seq_get_image_filename(args->seq, input_idx, filename)) {
			abort = 1;
			continue;
		}

		int thread_id = -1;
		int nb_subthreads = 1;
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
		if (have_seqwriter) {
			seqwriter_wait_for_memory();
			if (abort) {
				seqwriter_release_memory();
				continue;
			}
		}
		nb_subthreads = threads_per_image[thread_id];
#endif

		gboolean read_image = args->image_read_hook ? args->image_read_hook(args, input_idx) : TRUE;

		fits *fit = calloc(1, sizeof(fits));
		if (!fit) {
			PRINT_ALLOC_ERR;
			abort = 1;
			continue;
		}

		if (args->partial_image) {
			regdata *regparam = NULL;
			if (args->regdata_for_partial)
				regparam = args->seq->regparam[args->layer_for_partial];
			if (regparam &&
					guess_transform_from_H(regparam[input_idx].H) > NULL_TRANSFORMATION &&
					guess_transform_from_H(regparam[args->seq->reference_image].H) > NULL_TRANSFORMATION) {
				// do not try to transform area if img matrix is null
				selection_H_transform(&area,
						regparam[args->seq->reference_image].H,
						regparam[input_idx].H);
			}
			// args->area may be modified in hooks

			/* We need to detect if the box has crossed the borders to invalidate
			 * the current frame in case the box position was computed from reg data.
			 */
			gboolean has_crossed = enforce_area_in_image(&area, args->seq, input_idx) && args->regdata_for_partial;

			if (has_crossed || (read_image && seq_read_frame_part(args->seq, args->layer_for_partial,
						input_idx, fit, &area,
						args->get_photometry_data_for_partial, thread_id)))
			{
				if (args->stop_on_error)
					abort = 1;
				else {
					g_atomic_int_inc(&excluded_frames);
				}
				clearfits(fit);
				g_free(fit);
				// TODO: for seqwriter, we need to notify the failed frame
				continue;
			}
			/*char tmpfn[100];	// this is for debug purposes
			  sprintf(tmpfn, "/tmp/partial_%d.fit", input_idx);
			  savefits(tmpfn, fit);*/
		} else {
			// image is read bottom-up here, while it's top-down for partial images
			if (read_image && seq_read_frame(args->seq, input_idx, fit, args->force_float, thread_id)) {
				abort = 1;
				clearfits(fit);
				free(fit);
				continue;
			}
			// TODO: for seqwriter, we need to notify the failed frame
		}
		// checking nb layers consistency, not for partial image
		if (read_image && !args->partial_image && (fit->naxes[2] != args->seq->nb_layers)) {
			siril_log_color_message(_("Image #%d: number of layers (%d) is not consistent with sequence (%d), aborting\n"),
						"red", input_idx + 1, fit->naxes[2], args->seq->nb_layers);
			abort = 1;
			clearfits(fit);
			free(fit);
			continue;
		}
		// If not reading image, we still load its metadata to fill imgparam
		if (!read_image && seq_read_frame_metadata(args->seq, input_idx, fit)) {
			abort = 1;
			clearfits(fit);
			free(fit);
			continue;
		}

		if (read_image && args->image_hook(args, frame, input_idx, fit, &area, nb_subthreads)) {
			if (args->stop_on_error)
				abort = 1;
			else {
				g_atomic_int_inc(&excluded_frames);
				g_atomic_int_inc(&progress);
				set_progress_bar_data(NULL, (float)progress / nb_framesf);
			}
			if (args->seq->type == SEQ_INTERNAL) {
				fit->data = NULL;
				fit->fdata = NULL;
			}
			clearfits(fit);
			free(fit);
			// for seqwriter, we need to notify the failed frame
			if (have_seqwriter) {
				int retval;
				if (args->save_hook)
					retval = args->save_hook(args, frame, input_idx, NULL);
				else retval = generic_save(args, frame, input_idx, NULL);
				if (retval)
					abort = 1;
			}
			continue;
		}

		if (args->has_output) {
			int retval;
			if (args->save_hook)
				retval = args->save_hook(args, frame, input_idx, fit);
			else retval = generic_save(args, frame, input_idx, fit);
			if (retval) {
				abort = 1;
				clearfits(fit);
				free(fit);
				continue;
			}
		} else {
			/* save stats that may have been computed for the first
			 * time, but if fit has been modified for the new
			 * sequence, we shouldn't save it for the old one.
			 */
			save_stats_from_fit(fit, args->seq, input_idx);
		}

		if (!have_seqwriter) {
			if (!(args->seq->type == SEQ_INTERNAL)) {
				clearfits(fit);
				free(fit);
			} else {
				clearfits_header(fit);
				free(fit);
			}
		}

		g_atomic_int_inc(&progress);
		gchar *msg = g_strdup_printf(_("%s. Processing image %d (%s)"), args->description, input_idx + 1, filename);
		set_progress_bar_data(msg, (float)progress / nb_framesf);
		g_free(msg);
	}

	/* the finalize hook contains the sequence writer synchronization, it
	 * should be called before outputting the logs */
	if (have_seqwriter && args->finalize_hook && args->finalize_hook(args)) {
		siril_log_message(_("Finalizing sequence processing failed.\n"));
		abort = 1;
	}
	if (abort || excluded_frames == nb_frames) {
		set_progress_bar_data(_("Sequence processing failed. Check the log."), PROGRESS_RESET);
		siril_log_color_message(_("Sequence processing failed.\n"), "red");
		if (!abort)
			abort = 1;
		args->retval = abort;
	}
	else {
		if (excluded_frames) {
			set_progress_bar_data(_("Sequence processing partially succeeded. Check the log."), PROGRESS_RESET);
			siril_log_color_message(_("Sequence processing partially succeeded, with %d images that failed.\n"), "salmon", excluded_frames);
		} else {
			set_progress_bar_data(_("Sequence processing succeeded."), PROGRESS_RESET);
			siril_log_color_message(_("Sequence processing succeeded.\n"), "green");
		}
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

#ifdef _OPENMP
	omp_destroy_lock(&args->lock);
#endif
the_end:
#ifdef _OPENMP
	free(threads_per_image);
#endif

	if (index_mapping) free(index_mapping);
	if (!have_seqwriter && args->finalize_hook && args->finalize_hook(args)) {
		siril_log_message(_("Finalizing sequence processing failed.\n"));
		args->retval = 1;
	}
	int retval = args->retval;	// so we can free args if needed in the idle
#ifdef HAVE_FFMS2
	if (args->seq->type == SEQ_AVI && !com.headless)
		gui_function(set_seq_browser_active, GINT_TO_POINTER(1)); // re-enable the sequence browser if necessary
#endif

	if (!args->already_in_a_thread) {
		gboolean run_idle;
		if (args->idle_function)
			run_idle = siril_add_idle(args->idle_function, args) > 0; // not python safe because of the call to
					// update_sequences_list, which may be running after the next python command has started
					// and clobbers things
		else run_idle = siril_add_idle(end_generic_sequence, args) > 0; // the generic idle

		if (!run_idle) {
			// some generic cleanup for scripts
			// here we make sure the new .seq is created if needed
			if (args->has_output && args->load_new_sequence &&
				args->new_seq_prefix && !args->retval) {
				const gchar *basename = g_path_get_basename(args->seq->seqname);
				const char *root = remove_ext_from_filename(basename);
				const gchar *seqname = g_strdup_printf("%s%s%s", args->new_seq_prefix, root, ".seq");
				gboolean seqfilecreated = create_one_seq(seqname, args->seq->type);
				if (!seqfilecreated) { // just a fallback
					check_seq();
				}
			}
			free_generic_seq_args(args, TRUE);
			// we then make sure we stop the thread and reset cursor to FALSE
			end_generic(NULL);
		}
	}
	return GINT_TO_POINTER(retval);
}

void free_generic_seq_args(struct generic_seq_args *args, gboolean free_seq) {
	free(args->new_seq_prefix);
	if (free_seq && !check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);
	free(args);
}

// default idle function (in GTK main thread) to run at the end of the generic sequence processing
gboolean end_generic_sequence(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;

	if (args->has_output && args->load_new_sequence &&
			args->new_seq_prefix && !args->retval) {
		gchar *basename = g_path_get_basename(args->seq->seqname);
		gchar *seqname = g_strdup_printf("%s%s.seq", args->new_seq_prefix, basename);
		check_seq();
		update_sequences_list(seqname);
		g_free(seqname);
		g_free(basename);
	}
	// frees new_seq_prefix. This means that new_seq_prefix (or the seqEntry
	// string in function-specific structs) *must* always be allocated using
	// a stdlib freeable function, i.e. char* seqEntry = strdup("prefix_");
	// rather than char* seqEntry = "prefix_";
	free_generic_seq_args(args, TRUE);
	return end_generic(NULL);
}

/* If for_writer is false, it computes how many images can be processed in
 * parallel, with regard to how many of them can fit in memory. It returns at
 * most com.max_thread.
 * If for_writer is true, it computes how many images can be stored in the
 * queue. It returns at most 3 times com.max_thread.
 */
int seq_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image, MB_avail;
	int thread_limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float, NULL, &MB_per_image, &MB_avail);
	int limit = thread_limit;
	if (limit == 0) {
		gchar *mem_per_image = g_format_size_full(MB_per_image * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_image, mem_available);

		g_free(mem_per_image);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		/* number of threads for the main computation */
		if (limit > com.max_thread)
			limit = com.max_thread;
		if (for_writer) {
			/* we take the remainder for the writer */
			int max_queue_size = com.max_thread * 3;
			if (thread_limit - limit > max_queue_size)
				limit = max_queue_size;
			else limit = thread_limit - limit;
			if (limit < 1)
				limit = 1;	// bug #777
		}
#else
		if (!for_writer)
			limit = 1;
		else {
			limit = max(1, limit-1);
			if (limit > 3)
				limit = 3;
		}
#endif
	}
	return limit;
}

int seq_prepare_hook(struct generic_seq_args *args) {
	int retval = 0;
	g_assert(args->has_output); // don't call this hook otherwise
	if (!(args->seq->type == SEQ_INTERNAL))
		g_assert(args->new_seq_prefix); // don't call this hook otherwise
	//removing existing files
	if (args->seq->type != SEQ_INTERNAL) {
		remove_prefixed_sequence_files(args->seq, args->new_seq_prefix);
		remove_prefixed_star_files(args->seq, args->new_seq_prefix);
		remove_prefixed_drizzle_files(args->seq, args->new_seq_prefix);
	}
	if (args->force_ser_output || (args->seq->type == SEQ_SER && !args->force_fitseq_output)) {
		gchar *dest;
		const char *ptr = strrchr(args->seq->seqname, G_DIR_SEPARATOR);
		if (ptr)
			dest = g_strdup_printf("%s%s.ser", args->new_seq_prefix, ptr + 1);
		else dest = g_strdup_printf("%s%s.ser", args->new_seq_prefix, args->seq->seqname);

		args->new_ser = calloc(1, sizeof(struct ser_struct));
		if (ser_create_file(dest, args->new_ser, TRUE, args->seq->ser_file)) {
			free(args->new_ser);
			args->new_ser = NULL;
			retval = 1;
		}
		g_free(dest);
	}
	else if (args->force_fitseq_output || (args->seq->type == SEQ_FITSEQ && !args->force_ser_output)) {
		gchar *dest;
		const char *ptr = strrchr(args->seq->seqname, G_DIR_SEPARATOR);
		if (ptr)
			dest = g_strdup_printf("%s%s%s", args->new_seq_prefix, ptr + 1, com.pref.ext);
		else dest = g_strdup_printf("%s%s%s", args->new_seq_prefix, args->seq->seqname, com.pref.ext);

		args->new_fitseq = calloc(1, sizeof(fitseq));
		if (fitseq_create_file(dest, args->new_fitseq, args->nb_filtered_images)) {
			free(args->new_fitseq);
			args->new_fitseq = NULL;
			retval = 1;
		}
		g_free(dest);
	}
	else return 0;

	if (!retval)
		retval = seq_prepare_writer(args);
	return retval;
}

int seq_prepare_writer(struct generic_seq_args *args) {
	int limit = 0;
	if (args->compute_mem_limits_hook)
		limit = args->compute_mem_limits_hook(args, TRUE);
	else limit = seq_compute_mem_limits(args, TRUE); // the default

	if (limit == 0) {
		if (args->force_ser_output || (args->seq->type == SEQ_SER && !args->force_fitseq_output)) {
			ser_close_file(args->new_ser);
			free(args->new_ser);
			args->new_ser = NULL;
		}
		else if (args->force_fitseq_output || (args->seq->type == SEQ_FITSEQ && !args->force_ser_output)) {
			fitseq_close_file(args->new_fitseq);
			free(args->new_fitseq);
			args->new_fitseq = NULL;
		}
		return 1;
	}
	seqwriter_set_max_active_blocks(limit);
	return 0;
}

int seq_finalize_hook(struct generic_seq_args *args) {
	int retval = 0;
	g_assert(args->has_output); // don't call this hook otherwise
	if ((args->force_ser_output || (args->seq->type == SEQ_SER && !args->force_fitseq_output)) && args->new_ser) {
		retval = ser_write_and_close(args->new_ser);
		free(args->new_ser);
	}
	else if ((args->force_fitseq_output || (args->seq->type == SEQ_FITSEQ && !args->force_ser_output)) && args->new_fitseq) {
		retval = fitseq_close_file(args->new_fitseq);
		free(args->new_fitseq);
	}
	return retval;
}

	/* Internal sequences work differently to all other sequence types.
	 * Where other sequence types save their results to a new file,
	 * internal sequences update their results in the internal_fits[]
	 * array. They therefore need different handling within the save hook.
	 *
	 * We copy the fits metadata back to the internal_fits[in_index]
	 * pointer.
	 * Only the pointer to data / fdata was copied across in seq_read_frame
	 * We re-copy the pointers back, because they are set to NULL in
	 * copyfits. We then set the pointers in fit to NULL so that the memory
	 * doesn't get freed by the generic sequence worker.
	*/

static int generic_save_internal(struct generic_seq_args *args, int in_index, fits *fit) {
	clearfits_header(args->seq->internal_fits[in_index]);
	copyfits(fit, args->seq->internal_fits[in_index], CP_FORMAT, -1);
	copy_fits_metadata(fit, args->seq->internal_fits[in_index]);
	if (fit->type == DATA_USHORT) {
		args->seq->internal_fits[in_index]->data = fit->data;
		args->seq->internal_fits[in_index]->pdata[0] = fit->pdata[0];
		args->seq->internal_fits[in_index]->pdata[1] = fit->pdata[1];
		args->seq->internal_fits[in_index]->pdata[2] = fit->pdata[2];
	} else if (fit->type == DATA_FLOAT) {
		args->seq->internal_fits[in_index]->fdata = fit->fdata;
		args->seq->internal_fits[in_index]->fpdata[0] = fit->fpdata[0];
		args->seq->internal_fits[in_index]->fpdata[1] = fit->fpdata[1];
		args->seq->internal_fits[in_index]->fpdata[2] = fit->fpdata[2];
	}
	fit->data = NULL;
	fit->fdata = NULL;
	clearfits(fit);
	return 0;
}

/* In SER, all images must be in a contiguous sequence, so we use the out_index.
 * In FITS sequences, to keep track of image across processings, we keep the
 * input file number all along (in_index is the index in the sequence, not the name).
 * The 2nd condition ensures that any force condition prevails over opposite
 * input-type condition
 */
int generic_save(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	if (args->seq->type == SEQ_INTERNAL) {
		return generic_save_internal(args, in_index, fit);
	} else {
		if (args->force_ser_output || (args->seq->type == SEQ_SER && !args->force_fitseq_output)) {
			return ser_write_frame_from_fit(args->new_ser, fit, out_index);
		} else if (args->force_fitseq_output || (args->seq->type == SEQ_FITSEQ && !args->force_ser_output)) {
			return fitseq_write_image(args->new_fitseq, fit, out_index);
		} else {
			char *dest = fit_sequence_get_image_filename_prefixed(args->seq,
					args->new_seq_prefix, in_index);
			fit->bitpix = fit->orig_bitpix;
			int retval = savefits(dest, fit);
			free(dest);
			return retval;
		}
	}
}

/* Functions for use in sequence operations that can generate multiple sequences
 * of output images
 */
int multi_prepare(struct generic_seq_args *args) {
	struct multi_output_data *multi_args = (struct multi_output_data *) args->user;
	int n = max(multi_args->n, 1);
	multi_args->n = n; // Deal with cases where multi_args->n has not been set after calloc
	multi_args->new_ser = calloc(n, sizeof(char*));
	multi_args->new_fitseq = calloc(n, sizeof(char*));
	// we call the generic prepare n times with different prefixes
	for (int i = 0 ; i < n ; i++) {
		args->new_seq_prefix = strdup(multi_args->prefixes[i]);
		if (seq_prepare_hook(args))
			return 1;
		// but we copy the result between each call
		multi_args->new_ser[i] = args->new_ser;
		multi_args->new_fitseq[i] = args->new_fitseq;
		free(args->new_seq_prefix);
		args->new_seq_prefix = NULL;
	}

	args->new_seq_prefix = strdup(multi_args->prefixes[multi_args->new_seq_index]);
	args->new_ser = NULL;
	args->new_fitseq = NULL;

	seqwriter_set_number_of_outputs(n);
	return 0;
}

int multi_save(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	struct multi_output_data *multi_args = (struct multi_output_data *) args->user;
	struct _multi_split *multi_data = NULL;
	// images are passed from the image_hook to the save in a list, because
	// there are n, which is unsupported by the generic arguments
#ifdef _OPENMP
	omp_set_lock(&args->lock);
#endif
	GList *list = multi_args->processed_images;
	while (list) {
		if (((struct _multi_split *)list->data)->index == out_index) {
			multi_data = list->data;
			break;
		}
		list = g_list_next(multi_args->processed_images);
	}
	if (multi_data)
		multi_args->processed_images = g_list_remove(multi_args->processed_images, multi_data);
#ifdef _OPENMP
	omp_unset_lock(&args->lock);
#endif
	if (!multi_data) {
		siril_log_color_message(_("Image %d not found for writing\n"), "red", in_index);
		return 1;
	}
	int n = multi_args->n;
	siril_debug_print("Multiple (%d) images to be saved (%d)\n", n, out_index);
	for (int i = 0 ; i < n ; i++) {
		if (multi_data->images[i]->naxes[0] == 0) {
			siril_debug_print("empty data\n");
			return 1;
		}
	}

	int retval = 0;

	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		for (int i = 0 ; i < n ; i++) {
			retval |= ser_write_frame_from_fit(multi_args->new_ser[i], multi_data->images[i], out_index);
		}
		// the two fits are freed by the writing thread
		if (!retval) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		for (int i = 0 ; i < n ; i++) {
			retval |= fitseq_write_image(multi_args->new_fitseq[i], multi_data->images[i], out_index);
		}
		// the two fits are freed by the writing thread
		if (!retval) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else {
		for (int i = 0 ; i < n ; i++) {
			char *dest = fit_sequence_get_image_filename_prefixed(args->seq, multi_args->prefixes[i], in_index);
			retval |= savefits(dest, multi_data->images[i]);
			free(dest);
		}
	}
	gboolean tally_fitseq = FALSE, tally_ser = FALSE;
	for (int i = 0 ; i < n ; i++) {
		if (multi_args->new_fitseq[i])
			tally_fitseq = TRUE;
		if (multi_args->new_ser[i])
			tally_ser = TRUE;
	}
	if (!tally_ser && !tally_fitseq) { // detect if there is a seqwriter
		for (int i = 0 ; i < n ; i++) {
			clearfits(multi_data->images[i]);
			free(multi_data->images[i]);
		}
	}
	free(multi_data->images);
	free(multi_data);
	return retval;
}
void free_multi_args(struct multi_output_data *multi_args) {
	for (int i = 0 ; i < multi_args->n ; i++) {
		g_free(multi_args->prefixes[i]);
	}
	free(multi_args->prefixes);
	free(multi_args);
}

int multi_finalize(struct generic_seq_args *args) {
	struct multi_output_data *multi_args = (struct multi_output_data *) args->user;
	int n = multi_args->n;
	int retval = 0;
	for (int i = 0 ; i < n ; i++) {
		args->new_ser = multi_args->new_ser[i];
		args->new_fitseq = multi_args->new_fitseq[i];
		retval |= seq_finalize_hook(args);
		multi_args->new_ser[i] = NULL;
		multi_args->new_fitseq[i] = NULL;
	}
	seqwriter_set_number_of_outputs(1);
	free(multi_args->prefixes);
	free(multi_args->new_ser);
	free(multi_args->new_fitseq);
	free(multi_args->user_data);
	free(multi_args);
	return retval;
}

/*****************************************************************************
 *      P R O C E S S I N G      T H R E A D      M A N A G E M E N T        *
 ****************************************************************************/

static void set_thread_run(gboolean b);

static gboolean thread_being_waited = FALSE;

// This function is reentrant. The pointer will be freed in the idle function,
// so it must be a proper pointer to an allocated memory chunk.
gboolean start_in_new_thread(gpointer (*f)(gpointer), gpointer p) {
	g_mutex_lock(&com.mutex);

	if (com.run_thread || com.python_claims_thread || com.thread) {
		fprintf(stderr, "The processing thread is busy, stop it first.\n");
		g_mutex_unlock(&com.mutex);
		// We can't free p here as it may have unknown members. We must
		// indicate failure and allow the caller to clean up
		return FALSE;
	}

	com.run_thread = TRUE;
	com.thread = g_thread_new("processing", f, p);

	// Add a fake "child" to the list of child processes. This doesn't represent an external
	// program but it allows selecting to stop the processing thread instead of
	// an external program if one happens to be running
	if (!add_child((GPid) -2, INT_PROC_THREAD, "Siril processing thread")) {
		siril_log_color_message(_("Warning: failed to add processing thread to child list\n"), "salmon");
	}

	g_mutex_unlock(&com.mutex);
	set_cursor_waiting(TRUE);
	return TRUE; // indicate success
}

gboolean start_in_reserved_thread(gpointer (*f)(gpointer), gpointer p) {
	g_mutex_lock(&com.mutex);
	if (com.thread || com.python_claims_thread) {
		fprintf(stderr, "The processing thread is busy, stop it first.\n");
		g_mutex_unlock(&com.mutex);
		return FALSE;
	}

	com.run_thread = TRUE;
	com.thread = g_thread_new("processing", f, p);

	// Prepend this "child" to the list of child processes
	if (!add_child((GPid) -2, INT_PROC_THREAD, "Siril processing thread")) {
		siril_log_color_message(_("Warning: failed to add processing thread to child list\n"), "salmon");
	}

	g_mutex_unlock(&com.mutex);
	set_cursor_waiting(TRUE);
	return TRUE;
}

gpointer waiting_for_thread() {
	gpointer retval = NULL;
	if (com.thread) {
		thread_being_waited = TRUE;
		retval = g_thread_join(com.thread);
	}
	com.thread = NULL;
	thread_being_waited = FALSE;
	set_thread_run(FALSE);	// do it anyway in case of wait without stop

	return GINT_TO_POINTER(GPOINTER_TO_INT(retval) & ~CMD_NOTIFY_GFIT_MODIFIED);
}

int claim_thread_for_python() {
	g_mutex_lock(&com.mutex);
	if (com.thread) {
		fprintf(stderr, "The processing thread is busy (com.thread is not NULL). "
					"It must stop before Python can claim the processing thread.\n");
		g_mutex_unlock(&com.mutex);
		return 1;
	}
	if (com.run_thread) {
		fprintf(stderr, "The processing thread is busy (com.run_thread is TRUE). "
					"It must stop before Python can claim the processing thread.\n");
		g_mutex_unlock(&com.mutex);
		return 1;
	}
	if (com.python_claims_thread) {
		fprintf(stderr, "The processing thread is already claimed by Python. It must be "
					"released before Python can claim the processing thread again.\n");
		g_mutex_unlock(&com.mutex);
		return 1;
	}
	if (is_an_image_processing_dialog_opened()) {
		fprintf(stderr, "A Siril image processing dialog is open. It must be closed before "
					"Python can claim the processing thread.\n");
		g_mutex_unlock(&com.mutex);
		return 2;
	}
	com.python_claims_thread = TRUE;
	g_mutex_unlock(&com.mutex);
	set_cursor_waiting(TRUE);
	return 0;
}

void python_releases_thread() {
	g_mutex_lock(&com.mutex);
	com.python_claims_thread = FALSE;
	g_mutex_unlock(&com.mutex);
	set_cursor_waiting(FALSE);
	return;
}


static gboolean stop_processing_requested = FALSE;

static gboolean stop_processing_thread_idle(gpointer user_data) {
    remove_child_from_children((GPid) -2); // magic number indicating the processing thread
    if (com.thread == NULL) {
        siril_debug_print("The processing thread is not running (may have already finished).\n");
        return FALSE;
    }
    set_thread_run(FALSE);
    if (!thread_being_waited)
        waiting_for_thread();
    set_cursor_waiting(FALSE);
    return FALSE;
}

static gboolean check_stop_processing_request(gpointer user_data) {
    if (stop_processing_requested) {
        stop_processing_requested = FALSE;
        stop_processing_thread_idle(NULL);
    }
    return FALSE; // Remove this idle callback
}

void stop_processing_thread() {
    // Check if we're in headless mode first
    if (com.headless || !g_main_context_is_owner(g_main_context_default())) {
        execute_idle_and_wait_for_it(stop_processing_thread_idle, NULL);
    } else {
        // We're in the main thread, but we might be inside execute_idle_and_wait_for_it
        // Set a flag and queue an idle to handle it asynchronously
        stop_processing_requested = TRUE;
        gdk_threads_add_idle(check_stop_processing_request, NULL);
    }
}

static void set_thread_run(gboolean b) {
	siril_debug_print("setting run %d to the processing thread\n", b);
	g_mutex_lock(&com.mutex);
	com.run_thread = b;
	g_mutex_unlock(&com.mutex);
}

gboolean get_thread_run() {
	gboolean retval;
	g_mutex_lock(&com.mutex);
	retval = com.run_thread;
	g_mutex_unlock(&com.mutex);
	return retval;
}

// equivalent to atomic get and set if not running
gboolean reserve_thread() {
	gboolean retval;
	g_mutex_lock(&com.mutex);
	retval = !com.run_thread;
	if (retval)
		com.run_thread = TRUE;
	g_mutex_unlock(&com.mutex);
	return retval;
}

void unreserve_thread() {
	set_thread_run(FALSE);
}

/* should be called in a threaded function if nothing special has to be done at the end.
 * siril_add_idle(end_generic, NULL);
 */
gboolean end_generic(gpointer arg) {
	stop_processing_thread();
	set_cursor_waiting(FALSE);
	return FALSE;
}

/* wrapper for gdk_threads_add_idle that deactivates idle functions when
 * running in script/console mode
 * return 0 if idle was not added
 */
guint siril_add_idle(GSourceFunc idle_function, gpointer data) {
	if (!com.script && !com.python_command && !com.headless)
		return gdk_threads_add_idle(idle_function, data);
	return 0;
}

/* Must only ever be used for GTK updates. Do not call any functions that mess about
 * with files using this idle, it can break python scripts */
guint siril_add_pythonsafe_idle(GSourceFunc idle_function, gpointer data) {
	if (!com.script && !com.headless)
		return gdk_threads_add_idle(idle_function, data);
	return 0;
}

gboolean get_script_thread_run() {
	/* we don't need mutex for this one because it's only checked by user
	 * input, which is much slower than the function execution */
	return com.script_thread != NULL;
}

void wait_for_script_thread() {
	if (com.script_thread) {
		g_thread_join(com.script_thread);
		com.script_thread = NULL;
	}
}

static void free_child(child_info *child) {
	g_free(child->name);
	g_date_time_unref(child->datetime);
	g_free(child);
}

static void child_watch_cb(GPid pid, gint status, gpointer user_data) {
#ifdef G_OS_WIN32
    DWORD exit_code;
    if (GetExitCodeProcess(pid, &exit_code)) {
        if (exit_code == 0) {
            siril_debug_print("Child process completed successfully\n");
        } else {
            siril_debug_print("Plate solver exited with code %lu\n", exit_code);
        }
    }
    CloseHandle(pid);
#else
    // POSIX: manually reap the child
    int wait_status;
    pid_t result;

    // Try to reap - if it fails with ECHILD, someone else already did
    result = waitpid(pid, &wait_status, WNOHANG);

    if (result == pid) {
        // We successfully reaped it
        if (WIFEXITED(wait_status)) {
            int exit_status = WEXITSTATUS(wait_status);
            if (exit_status == 0) {
                siril_debug_print("Child process completed successfully\n");
            } else {
                siril_debug_print("Child process exited with status %d\n", exit_status);
            }
        } else if (WIFSIGNALED(wait_status)) {
            siril_debug_print("Child process terminated by signal %d\n", WTERMSIG(wait_status));
        }
    } else if (result == 0) {
        // Process hasn't exited yet (shouldn't happen)
        siril_debug_print("Warning: child watch called but process still running\n");
    } else if (result == -1 && errno == ECHILD) {
        // Already reaped
        siril_debug_print("Warning: child process already reaped\n");
    } else {
        siril_debug_print("Error waiting for child: %s\n", strerror(errno));
    }
    g_spawn_close_pid(pid);
#endif
    remove_child_from_children(pid);
}

gboolean add_child(GPid child_pid, int program, const gchar *name) {
	// Add the callback to remove it on completion
	if (program != EXT_PYTHON && program != INT_PROC_THREAD)
		g_child_watch_add(child_pid, child_watch_cb, NULL);

	child_info *child = g_malloc(sizeof(child_info));
	if (!child) {
		return FALSE;
	}

	child->name = g_strdup(name);
	if (!child->name) {
		g_free(child);
		return FALSE;
	}

	child->datetime = g_date_time_new_now_local();
	if (!child->datetime) {
		g_free(child->name);
		g_free(child);
		return FALSE;
	}
	child->childpid = child_pid;
	child->program = program;

	child_mutex_lock();
	com.children = g_slist_prepend(com.children, child);
	child_mutex_unlock();

	return TRUE;
}

// This must be called in the g_spawn on_exit callback to remove the defunct child from the list
void remove_child_from_children(GPid pid) {
	child_mutex_lock();
	GSList *prev = NULL;
	GSList *iter = com.children;

	while (iter) {
		child_info *child = (child_info*) iter->data;

		if (!child) {
			// Save next pointer before freeing
			GSList *next = iter->next;

			if (prev)
				prev->next = next;
			else
				com.children = next;  // Don't forget to update head!

			g_slist_free_1(iter);
			iter = next;
			continue;
		}

		if (child->childpid == pid) {
			// Save next pointer before freeing
			GSList *next = iter->next;

			// Remove the node from the list
			if (prev == NULL) {
				com.children = next;
			} else {
				prev->next = next;
			}

			// Free the node and the child
			g_slist_free_1(iter);
			free_child(child);
			siril_debug_print("Removed GPid %d from com.children\n", pid);
			child_mutex_unlock();
			return;
		}

		// Move to next node
		prev = iter;
		iter = iter->next;
	}
	child_mutex_unlock();
}

// kills external calls
// if onexit is TRUE, also do some cleaning
void kill_child_process(GPid pid, gboolean onexit) {
	if (onexit)
		printf("Making sure no child is left behind...\n");
	// Find the correct child in the list
	// We can safely iterate over the list without worrying about simultaneously adding
	// to it from another thread as it would only ever be prepended, so the rest of the
	// list will still be fine. This means we can have the mutex control for child
	// removal in remove_child_from_children, which makes things considerably simpler.
	GSList *iter = com.children;
	gboolean success = FALSE;
	while (iter) {
		child_info *child = (child_info*) iter->data;
		if (!child)
			return;

		GSList *next = iter->next;

		if (child->childpid == pid || onexit) {
			if (child->program == INT_PROC_THREAD) {
				stop_processing_thread();
			} else {
				if (child->program == EXT_STARNET || child->program == EXT_PYTHON) {
#ifdef _WIN32
					TerminateProcess((void *) child->childpid, 1);
#else
					kill((pid_t) child->childpid, SIGKILL);
#endif
				} else if (child->program == EXT_ASNET) {
					FILE* fp = g_fopen("stop", "w");
					if (fp != NULL)
						fclose(fp);
					if (onexit) {
						g_usleep(1000);
						if (g_unlink("stop"))
							siril_debug_print("g_unlink() failed\n");
						siril_debug_print("asnet has been stopped on exit\n");
					}
					stop_processing_thread(); // we stop the sequence worker as well
				}
				// No need to manually remove the child as this is done by child_watch_cb
			}

			siril_log_color_message(_("Process aborted by user\n"), "red");
			success = TRUE;
			if (!onexit)
				break;
		}
		iter = next;
	}
	// If we get here without success, no matching PID was found
	if (!success && pid != (GPid) -1)
		siril_log_message(_("Failed to find GPid %d, it may already have exited...\n"), pid);
}

static void check_if_child_is_python(gpointer data, gpointer user_data) {
	// Cast the input pointers to their correct types
	child_info *child = (child_info *)data;
	if(!child) return;
	gboolean *is_ok = (gboolean *)user_data;

	// Convert name to lowercase for case-insensitive comparison
	gchar *lowercase_name = g_ascii_strdown(child->name, -1);

	// Check if the lowercased name contains "python"
	if (g_strstr_len(lowercase_name, -1, "python") != NULL) {
		// If "python" is found in the name, set is_ok to FALSE
		*is_ok = FALSE;
	}
	// Free the lowercase string
	g_free(lowercase_name);
}

void check_python_flag() {
	siril_debug_print("checking python flag\n");
	gboolean is_ok = TRUE;
	child_mutex_lock();
	g_slist_foreach(com.children, check_if_child_is_python, &is_ok);
	child_mutex_unlock();
	if (is_ok) {
		com.python_script = FALSE;
		siril_debug_print("com.python_script cleared\n");
	}
}

static void kill_if_child_is_python(gpointer data, gpointer user_data) {
	// Cast the input pointers to their correct types
	child_info *child = (child_info *)data;

	// Convert name to lowercase for case-insensitive comparison
	gchar *lowercase_name = g_ascii_strdown(child->name, -1);

	// Check if the lowercased name contains "python"
	if (g_strstr_len(lowercase_name, -1, "python") != NULL) {
		kill_child_process(child->childpid, FALSE);
	}
	// Free the lowercase string
	g_free(lowercase_name);
}

void kill_all_python_scripts() {
	g_slist_foreach(com.children, kill_if_child_is_python, NULL);
}

void on_processes_button_cancel_clicked(GtkButton *button, gpointer user_data) {
	child_mutex_lock();
	guint children = g_slist_length(com.children);
	if (children > 1) {
		child_mutex_unlock();
		GPid pid = show_child_process_selection_dialog(com.children);
		kill_child_process(pid, FALSE);
	} else if (children == 1) {
		child_info *child = (child_info*) com.children->data;
		GPid pid = child->childpid;
		child_mutex_unlock();
		kill_child_process(pid, FALSE);
	} else {
		child_mutex_unlock();
		com.stop_script = TRUE;
		stop_processing_thread();
		wait_for_script_thread();
	}

	if (!com.headless) {
		script_widgets_enable(TRUE);
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	}
}

struct generic_seq_args *create_default_seqargs(sequence *seq) {
	struct generic_seq_args *args = calloc(1, sizeof(struct generic_seq_args));
	args->seq = seq;
	args->filtering_criterion = seq_filter_all;
	args->nb_filtered_images = seq->number;
	args->stop_on_error = TRUE;
	args->upscale_ratio = 1.0;
	args->parallel = TRUE;
	return args;
}

gpointer generic_sequence_metadata_worker(gpointer arg) {
	struct generic_seq_metadata_args *args = (struct generic_seq_metadata_args *)arg;
	struct timeval t_start, t_end;
	set_progress_bar_data(NULL, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);
	int input_idx, frame, retval = 0;
	int *index_mapping = NULL;

	int nb_frames = compute_nb_filtered_images(args->seq, args->filtering_criterion, 1.0);
	if (nb_frames <= 0) {
		siril_log_message(_("No image selected for processing, aborting\n"));
		retval = 1;
		goto cleanup;
	}

	if (args->filtering_criterion) {
		index_mapping = malloc(nb_frames * sizeof(int));
		if (!index_mapping) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto cleanup;
		}
		for (input_idx = 0, frame = 0; input_idx < args->seq->number; input_idx++) {
			/* the third parameter, filtering_parameter, is not required here, so
			 * we set an arbitrary value of 1.0. The only relevant filtering is
			 * "by selection".
			 */
			if (!args->filtering_criterion(args->seq, input_idx, 1.0)) {
				continue;
			}
			index_mapping[frame++] = input_idx;
		}
		if (frame != nb_frames) {
			siril_log_message(_("Output index mapping failed (%d/%d).\n"), frame, nb_frames);
			retval = 1;
			goto cleanup;
		}
	}
	if (args->seq->type == SEQ_FITSEQ && (!args->seq->fitseq_file || !args->seq->fitseq_file->fptr)) {
		// List of extensions to check (both upper and lower case)
		const char *extensions[] = {".fit", ".fits", ".fts", ".FIT", ".FITS", ".FTS"};
		int num_extensions = 6;
		const gchar *com_ext = get_com_ext(args->seq->fz);

		retval = 1; // Default to fail status

		// First, try com_ext (the expected extension)
		gchar *seqfilename = g_strdup_printf("%s%s", args->seq->seqname, com_ext);
		if (g_file_test(seqfilename, G_FILE_TEST_EXISTS)) {
			retval = fitseq_open(seqfilename, args->seq->fitseq_file, 0);
			g_free(seqfilename);
			seqfilename = NULL;
			if (!retval) goto after_fitseq_check;
		}
		g_free(seqfilename);

		// If com_ext didn't match, try other extensions
		for (int i = 0; i < num_extensions; i++) {
			// Skip if this extension matches com_ext (already checked)
			if (args->seq->fz) {
				// For compressed files, check if extension + .fz matches com_ext
				char temp_ext[20];
				snprintf(temp_ext, 19, "%s.fz", extensions[i]);
				if (strcmp(temp_ext, com_ext) == 0) continue;
			} else {
				// For uncompressed files, check if extension matches com_ext
				if (strcmp(extensions[i], com_ext) == 0) continue;
			}

			if (args->seq->fz) {
				seqfilename = g_strdup_printf("%s%s.fz", args->seq->seqname, extensions[i]);
			} else {
				seqfilename = g_strdup_printf("%s%s", args->seq->seqname, extensions[i]);
			}

			if (g_file_test(seqfilename, G_FILE_TEST_EXISTS)) {
				retval = fitseq_open(seqfilename, args->seq->fitseq_file, 0);
				g_free(seqfilename);
				if (!retval) break;
			} else {
				g_free(seqfilename);
			}
		}

after_fitseq_check:
		if (retval) goto cleanup;
	}
	for (frame = 0; frame < nb_frames; frame++) {
		if (index_mapping)
			input_idx = index_mapping[frame];
		else input_idx = frame;

		fits fit = { 0 };
		if (args->seq->type == SEQ_REGULAR && seq_open_image(args->seq, input_idx)) {
			retval = 1;
			goto cleanup;
		} else if (args->seq->type == SEQ_FITSEQ) {
			// Need to move to the correct HDU
			fits_movabs_hdu(args->seq->fitseq_file->fptr, args->seq->fitseq_file->hdu_index[input_idx], NULL, &retval);
			if (retval) {
				siril_log_color_message(_("Error finding the HDU for frame %d\n"), "red", input_idx);
				goto cleanup;
			}
		}
		if (args->seq->type == SEQ_REGULAR)
			args->image_hook(args, args->seq->fptr[input_idx], input_idx);
		else args->image_hook(args, args->seq->fitseq_file->fptr, input_idx);
		if (args->seq->type == SEQ_REGULAR)
			seq_close_image(args->seq, input_idx);
		clearfits(&fit);
	}
cleanup:
	if (args->seq->type == SEQ_FITSEQ) {
		int status;
		fits_movabs_hdu(args->seq->fitseq_file->fptr, args->seq->fitseq_file->hdu_index[0], NULL, &status);
	}
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	if (!check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);
	g_slist_free_full(args->keys, g_free);
	g_free(args->header);
	if (args->output_stream)
		g_object_unref(args->output_stream);
	if (index_mapping)
		free(index_mapping);
	free(args);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

/********** per-image threading **********/
/* convention: call the variable threading if it can be MULTI_THREADED, call it
 * threads if it contains the other values, so the number of threads.
 * public functions should use threading, and use this function, private
 * functions can use threads if their caller has done the job.
 */
int check_threading(threading_type *threads) {
	if (*threads == MULTI_THREADED)
		*threads = com.max_thread;
	return *threads;
}

/* same as check_threading but also returns a number of threads adapted to a
 * set size */
int limit_threading(threading_type *threads, int min_iterations_per_thread, size_t total_iterations) {
	if (*threads == MULTI_THREADED)
		*threads = com.max_thread;
	int max_chunks = total_iterations / min_iterations_per_thread;
	if (max_chunks < 1)
		max_chunks = 1;
	if (max_chunks < *threads) {
		siril_debug_print("limiting operation to %d threads (%d allowed)\n", max_chunks, *threads);
		return max_chunks;
	}
	return *threads;
}

/* we could use several threads on a single image if there are less workers
 * than available threads on the system, due to memory constraints.
 * This computes the number of threads each worker can use, distributing
 * max threads to all workers.
 */
int *compute_thread_distribution(int nb_workers, int max) {
	int *threads = malloc(nb_workers * sizeof(int));
	if (!threads) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	int base = max / nb_workers;
	int rem = max % nb_workers;
	siril_debug_print("distributing %d threads to %d workers\n", max, nb_workers);
	for (int i = 0; i < nb_workers; i++) {
		threads[i] = base;
		if (rem > 0) {
			threads[i]++;
			rem--;
		}
		if (threads[i] == 0)
			threads[i] = 1;
		siril_debug_print("thread %d has %d subthreads\n", i, threads[i]);
	}
	return threads;
}

/* Image processing worker and hooks */

/** Default memory check hook - ensures enough memory for mem_ratio * image size */
static int default_img_mem_hook(struct generic_img_args *args) {
	if (!args || !args->fit)
		return 1;

	gint64 required_mem = (gint64)(args->mem_ratio *
		args->fit->rx * args->fit->ry * args->fit->naxes[2] *
		(args->fit->type == DATA_FLOAT ? sizeof(float) : sizeof(WORD)));

	gint64 available_mem = get_available_memory();

	if (required_mem > available_mem) {
		if (args->verbose)
			siril_log_color_message(_("Not enough memory for operation: need %.1f MB, have %.1f MB\n"),
				"red", (double)required_mem / BYTES_IN_A_MB, (double)available_mem / BYTES_IN_A_MB);
		return 1;
	}

	if (args->verbose)
		siril_debug_print("Memory check passed: need %.1f MB, have %.1f MB\n",
			(double)required_mem / BYTES_IN_A_MB, (double)available_mem / BYTES_IN_A_MB);
	return 0;
}

/** Default idle function to end generic image processing */
gboolean end_generic_image(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args*) p;
	stop_processing_thread();
	free_generic_img_args(args);
	return FALSE;
}

static gboolean end_generic_image_update_gfit(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args*) p;
	stop_processing_thread();
	notify_gfit_modified();
	free_generic_img_args(args);
	return FALSE;
}

static gboolean end_generic_image_update_gfit_mask(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args*) p;
	mask_t *mask = fits_to_mask(args->fit);
	clearfits(args->fit);
	free(args->fit);
	if (gfit->mask)
		free(gfit->mask);
	gfit->mask = mask;

	stop_processing_thread();

	free_generic_img_args(args);
	return FALSE;
}

/** Free generic_img_args structure */
void free_generic_img_args(struct generic_img_args *args) {
	if (!args)
		return;
	if (args->user)
		destroy_any_args(args->user);
	free(args);
}

FAST_MATH_PUSH
static inline void blend_fits_with_mask(fits* fit, fits* orig) {
	size_t npixels = fit->rx * fit->ry;
	mask_t *mask = fit->mask;

	if (!mask || !mask->data)
		return;

	uint8_t bitpix = mask->bitpix;

	if (fit->type == DATA_USHORT) {

		for (int chan = 0; chan < fit->naxes[2]; chan++) {
			WORD *ocdata = orig->pdata[chan];
			WORD *fcdata = fit->pdata[chan];

			switch (bitpix) {

			case 8: {
				uint8_t * restrict m = (uint8_t*)mask->data;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
				for (size_t i = 0; i < npixels; i++) {
					uint32_t temp =
						fcdata[i] * m[i] +
						ocdata[i] * (255 - m[i]);
					fcdata[i] = (temp * 257) >> 16;
				}
				break;
			}

			case 16: {
				uint16_t * restrict m = (uint16_t*)mask->data;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
				for (size_t i = 0; i < npixels; i++) {
					uint32_t alpha = m[i];          // 0..65535
					uint32_t temp =
						fcdata[i] * alpha +
						ocdata[i] * (65535 - alpha);
					fcdata[i] = temp >> 16;         // exact /65535
				}
				break;
			}

			case 32: {
				float * restrict m = (float*)mask->data;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
				for (size_t i = 0; i < npixels; i++) {
					float alpha = m[i];
					float origv = ocdata[i];
					float blend = fcdata[i] - origv;
					float out = origv + alpha * blend;
					if (out < 0.f) out = 0.f;
					if (out > 65535.f) out = 65535.f;
					fcdata[i] = (WORD)(out + 0.5f);
				}
				break;
			}

			default:
				return;
			}
		}

	} else {

		for (int chan = 0; chan < fit->naxes[2]; chan++) {
			float *ocdata = orig->fpdata[chan];
			float *fcdata = fit->fpdata[chan];

			switch (bitpix) {

			case 8: {
				uint8_t * restrict m = (uint8_t*)mask->data;
				const float scale = 1.0f / 255.0f;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
				for (size_t i = 0; i < npixels; i++) {
					float alpha = m[i] * scale;
					float origv = ocdata[i];
					fcdata[i] = origv + alpha * (fcdata[i] - origv);
				}
				break;
			}

			case 16: {
				uint16_t * restrict m = (uint16_t*)mask->data;
				const float scale = 1.0f / 65535.0f;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
				for (size_t i = 0; i < npixels; i++) {
					float alpha = m[i] * scale;
					float origv = ocdata[i];
					fcdata[i] = origv + alpha * (fcdata[i] - origv);
				}
				break;
			}

			case 32: {
				float * restrict m = (float*)mask->data;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread)
#endif
				for (size_t i = 0; i < npixels; i++) {
					float alpha = m[i];
					float origv = ocdata[i];
					fcdata[i] = origv + alpha * (fcdata[i] - origv);
				}
				break;
			}

			default:
				return;
			}
		}
	}
}
FAST_MATH_POP

/** Main generic image worker function
 * Works with a single image, optionally using multiple threads for processing
 */
gpointer generic_image_worker(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	struct timeval t_start, t_end;

	assert(args);
	assert(args->fit);
	assert(args->image_hook);

	gboolean verbose = args->verbose || !args->for_preview;
	gchar *history = NULL;
	gchar *summary = NULL;

	set_progress_bar_data(NULL, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);
	args->retval = 0;

	fits *orig = NULL;

	gboolean using_mask = args->mask_aware && args->fit->mask && args->fit->mask_active;
	// Create a copy so we still have the original fit for combining with the result
	// according to a mask
	if (using_mask) {
		orig = calloc(1, sizeof(fits));
		if (!orig) {
			PRINT_ALLOC_ERR;
			args->retval = 1;
			goto the_end_no_orig;
		}
		if (copyfits(args->fit, orig, CP_ALLOC | CP_FORMAT | CP_COPYA, -1)) {
			siril_log_color_message(_("Failed to copy original image.\n"), "red");
			args->retval = 1;
			goto the_end_no_orig;
		}
	}

	// Set default max_threads if not specified
	if (args->max_threads < 1)
		args->max_threads = com.max_thread;

	// Memory check
	if (args->mem_ratio > 0.0f) {
		if (default_img_mem_hook(args)) {
			args->retval = 1;
			goto the_end;
		}
	}

	// Output print of operation description
	if (args->description && verbose) {
		gchar *desc = g_strdup_printf(_("%s: processing...\n"), args->description);
		siril_log_color_message(desc, "green");
		g_free(desc);
	}

	set_progress_bar_data(_("Processing image..."), 0.1f);

	// Call the image processing hook - operates in-place on args->fit
	if (args->image_hook(args, args->fit, args->max_threads)) {
		siril_log_color_message(_("%s image processing failed.\n"), "red", args->description);
		args->retval = 1;
	} else {
		// Check for a mask hook and call it if one exists (we must do this even if not using the mask, to ensure
		// it is properly updated)
		if (args->fit->mask && args->mask_hook && args->mask_hook(args)) { // This hook updates the mask (eg to change geometry)
			// There is a mask_hook but it failed
			if (args->verbose)
				siril_log_color_message(_("%s mask processing failed.\n"), "red", args->description);
			args->retval = 1;
			// Put the image back how it was, as failing to update the mask could be catastrophic
			fits *tmp = gfit;
			gfit = orig;
			clearfits(tmp);
			free(tmp);
		} else { // Either no mask_hook or the mask_hook returned 0 (success)
			// Blend according to the mask
			if (using_mask) {
				blend_fits_with_mask(args->fit, orig);
			}
			// If there is a log_hook, set the HISTORY card and update the log as required
			// Generate the message used for undo label and HISTORY, ideally from the log hook but we use the simple description as a backup
			history = args->log_hook ? args->log_hook(args->user, DETAILED): g_strdup(args->description); // Dynamically allocates memory
			if (!args->for_preview)
				siril_log_message("%s\n", history); // Log the full detailed description

			// If we are being run from the GUI and not just updating a preview, set the undo state
			if (args->fit == gfit && !(args->custom_undo || args->for_preview || args->command)) {
				summary = args->log_hook ? args->log_hook(args->user, SUMMARY): g_strdup(args->description);
				undo_save_state(orig, summary); // We just use the short description here
				g_free(summary); // free the message
				g_free(history); // free the full description, we don't need it in this case
			} else {
				// Update the HISTORY card if it wasn't updated by undo_save_state'
				if ((args->custom_undo)) {
					// Do nothing, the custom undo will handle it
					g_free(history);
				} else if (args->command_updates_gfit) {
					args->fit->history = g_slist_append(args->fit->history, history);
					// args->fit->history now owns the allocated memory, we must not free it if this codepath is taken
					update_fits_header(args->fit); // update the header so the history is up to date and correctly ordered
				} else {
					g_free(history);
				}
			}
			if (verbose) {
				siril_log_color_message(_("%s succeeded.\n"), "green", args->description);
				gettimeofday(&t_end, NULL);
				show_time(t_start, t_end);
			}
		}
	}

the_end:
	// Always clean up orig if it was allocated
	if (orig) clearfits(orig);
	free(orig);

the_end_no_orig:
	if (args->retval) {
		set_progress_bar_data(_("Image processing failed. Check the log."), PROGRESS_RESET);
	} else {
		set_progress_bar_data(_("Image processing succeeded."), PROGRESS_DONE);
	}

	int retval = args->retval;

	/*
	 *
	 */
	if (args->updates_mask) {
		execute_idle_and_wait_for_it(end_generic_image_update_gfit_mask, args);
	} else if (args->command) {
		// commands do not use custom idles and the generic ones must run synchronously
		if (args->command_updates_gfit) {
			execute_idle_and_wait_for_it(end_generic_image_update_gfit, args);
		} else {
			execute_idle_and_wait_for_it(end_generic_image, args);
		}
	} else if (args->idle_function) {
		siril_add_idle(args->idle_function, args);
	} else {
		siril_add_idle(end_generic_image_update_gfit, args);
	}

	return GINT_TO_POINTER(retval);
}

