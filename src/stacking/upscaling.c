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

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/sequence_filtering.h"
#include "core/siril_log.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
#include "gui/image_interactions.h"
#include "gui/utils.h" // for delete_selected_area()
#include "opencv/opencv.h"

#include "stacking.h"

#define TMP_UPSCALED_PREFIX "tmp_upscaled_"

void remove_tmp_upscaled_files(struct stacking_args *args) {
	int i;
	if (!args->upscale_at_stacking)
		return;

	gchar *basename = g_path_get_basename(args->seq->seqname);
	if (!g_str_has_prefix(basename, TMP_UPSCALED_PREFIX)) {
		remove_prefixed_sequence_files(args->seq, TMP_UPSCALED_PREFIX);
		g_free(basename);
		return;
	}
	// else means that we are removing files after processing and that
	// we have access to files that were created for this processing

	gchar *seqname;
	int len = strlen(basename) + 5;
	seqname = malloc(len);
	g_snprintf(seqname, len, "%s.seq", basename);
	siril_debug_print("Removing %s\n", seqname);
	if (g_unlink(seqname))
		siril_debug_print("g_unlink() failed\n"); // removing the seqfile
	free(seqname);
	g_free(basename);

	switch (args->seq->type) {
	default:
	case SEQ_REGULAR:
		for (i = 0; i < args->seq->number; i++) {
			char filename[500];
			// FIXME: no preallocation of file name
			fit_sequence_get_image_filename(args->seq, args->image_indices[i], filename, TRUE);
			siril_debug_print("Removing %s\n", filename);
			if (g_unlink(filename))
				siril_debug_print("g_unlink() failed\n");
		}
		break;
	case SEQ_SER:
		siril_debug_print("Removing %s\n", args->seq->ser_file->filename);
		if (g_unlink(args->seq->ser_file->filename))
			siril_debug_print("g_unlink() failed\n");
		ser_close_file(args->seq->ser_file);
		break;
	case SEQ_FITSEQ:
		siril_debug_print("Removing %s\n", args->seq->fitseq_file->filename);
		if (g_unlink(args->seq->fitseq_file->filename))
			siril_debug_print("g_unlink() failed\n");
		fitseq_close_file(args->seq->fitseq_file);
		break;
	}
}

/*****************************************************************
 *      UP-SCALING A SEQUENCE: GENERIC FUNCTION IMPLEMENTATION   *
 *****************************************************************/

/* stacking an up-scaled sequence is a bit of a trick;
 * stacking a sequence is normally 3 steps (see stack_function_handler):
 * computing the normalization parameters, stacking the sequence, saving and
 * displaying the result. With the up-scale temporarily added in the middle, to
 * provide a cheap version of the drizzle algorithm, we have to create an
 * up-scaled sequence and pass it to the stacking operation seamlessly. The
 * problem with this is that at the end of the stacking, we have to close the
 * up-scaled sequence, maintain the original sequence as loaded, and display an
 * image, the result, that has a different size than the sequence's.
 */

struct upscale_args {
	double factor;
};

static int upscale_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	double factor = ((struct upscale_args *)args->user)->factor;
	/* updating pixel size if exist */
	fit->keywords.pixel_size_x /= factor;
	fit->keywords.pixel_size_y /= factor;

	return cvResizeGaussian(fit,
			round_to_int(fit->rx * factor),
			round_to_int(fit->ry * factor), OPENCV_NEAREST, FALSE);
}

int upscale_sequence(struct stacking_args *stackargs) {
	if (!stackargs->upscale_at_stacking)
		return 0;

	struct generic_seq_args *args = create_default_seqargs(stackargs->seq);
	struct upscale_args *upargs = calloc(1, sizeof(struct upscale_args));

	upargs->factor = 2.;

	args->filtering_criterion = stackargs->filtering_criterion;
	args->filtering_parameter = stackargs->filtering_parameter;
	args->nb_filtered_images = stackargs->nb_images_to_stack;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = seq_finalize_hook;
	args->image_hook = upscale_image_hook;
	args->description = _("Up-scaling sequence for stacking");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = upargs->factor;
	args->new_seq_prefix = strdup(TMP_UPSCALED_PREFIX);
	args->user = upargs;
	args->already_in_a_thread = TRUE;

	// check memory here, failure is not an error
	int nb_threads = seq_compute_mem_limits(args, FALSE);
	if (nb_threads == 0) {
		siril_log_color_message(_("Stacking will be done without up-scaling (disabling 'drizzle')\n"), "red");
		stackargs->upscale_at_stacking = FALSE;
		free(upargs);
		free_generic_seq_args(args, TRUE);
		return 0;
	}
	args->max_parallel_images = nb_threads;

	remove_tmp_upscaled_files(stackargs);

	generic_sequence_worker(args);

	stackargs->retval = args->retval;
	free(upargs);
	free(args);

	if (!stackargs->retval) {
		gchar *basename = g_path_get_basename(stackargs->seq->seqname);
		char *seqname = malloc(strlen(TMP_UPSCALED_PREFIX) + strlen(basename) + 5);
		sprintf(seqname, "%s%s.seq", TMP_UPSCALED_PREFIX, basename);
		if (g_unlink(seqname))
			siril_debug_print("g_unlink() failed\n");
		g_free(basename);

		// replace active sequence by upscaled
		if (check_seq()) {	// builds the new .seq
			free(seqname);
			return 1;
		}

		sequence *oldseq = stackargs->seq;
		sequence *newseq = readseqfile(seqname);
		if (!newseq) {
			free(seqname);
			return 1;
		}
		free(seqname);

		/* The original sequence and the up-scaled sequence differ by:
		 * - size, managed in the seq_check_basic_data() call below
		 * - images list, if there is a filter that excluded some images of the
		 *   sequence, excluded images are not up-scaled. resulting in the up-scaled
		 *   sequence having new contiguous image numbers.
		 *   Registration data is copied image per image below and the image_indices
		 *   array is rebuilt to identity in the stack_fill_list_of_unfiltered_images()
		 *   call after that.
		 * - registration data, since they are copied unmodified, when using the shifts
		 *   in stacking, they must be multiplied by the factor upscale_at_stacking.
		 */
		if (seq_check_basic_data(newseq, FALSE) == -1) {
			free(newseq);
			stackargs->retval = -1;
			return stackargs->retval;
		}
		stackargs->seq = newseq;
		stackargs->filtering_criterion = seq_filter_all;
		stackargs->filtering_parameter = 0.0;
		stackargs->nb_images_to_stack = newseq->number;

		newseq->reference_image = find_refimage_in_indices(stackargs->image_indices,
				stackargs->nb_images_to_stack, stackargs->ref_image);
		stackargs->ref_image = newseq->reference_image;
		newseq->regparam[stackargs->reglayer] = malloc(stackargs->nb_images_to_stack * sizeof(regdata));
		int i;
		for (i = 0; i < stackargs->nb_images_to_stack; i++) {
			regdata *data = &oldseq->regparam[stackargs->reglayer][stackargs->image_indices[i]];
			memcpy(&newseq->regparam[stackargs->reglayer][i], data, sizeof(regdata));
			// TODO: why don't we modify the shifts here already? indeed!
		}
		stackargs->retval = stack_fill_list_of_unfiltered_images(stackargs);

		// don't free oldseq, it's either still com.seq with GUI or freed in
		// stack_one_seq in scripts
		delete_selected_area();
	}
	return stackargs->retval;
}
