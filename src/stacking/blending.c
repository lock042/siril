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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "io/sequence.h"
#include "opencv/opencv.h"

#include "stacking/stacking.h"

#define MASK_SCALE 0.1
#define RAMP_PACE 1000

static float *ramp_array = NULL;

void init_ramp() {
	if (ramp_array)
		return;
	int nbpoints = RAMP_PACE + 1;
	ramp_array = malloc(nbpoints * sizeof(float));
	float norm = 1.f / (float)RAMP_PACE;
	for (int i = 0; i < nbpoints; i++) {
		float r = (float)i * norm;
		ramp_array[i] = r * r * r * (6.f * r * r - 15.f * r + 10.f);
	}
	siril_debug_print("ramp array initialized\n");
}

float get_ramped_value(float val) {
	int index = (int)(val * RAMP_PACE);
	return ramp_array[index];
}

void compute_downscaled_mask_size(int rx, int ry, int *rx_out, int *ry_out, double *fx, double *fy) {
	*rx_out = (int)(rx * MASK_SCALE);
	*ry_out = (int)(ry * MASK_SCALE);
	if (fx)
		*fx = (double)(*rx_out) / (double)rx;
	if (fy)
		*fy = (double)(*ry_out) / (double)ry;
}

// check if we need to create the mask or if it already exists
static gboolean compute_mask_read_hook(struct generic_seq_args *args, int i) {
	const gchar *mask_filename = get_cache_filename(args->seq, i, "msk", NULL);
	if (!mask_filename) {
		return TRUE;
	}
	if (!g_file_test(mask_filename, G_FILE_TEST_EXISTS)) { // mask file does not exist, we need to read and create the file
		return TRUE;
	}
	if (check_cachefile_date(args->seq, i, mask_filename)) { // the mask exists and is more recent than the img, we don't need to read again
		siril_log_message(_("Mask for image %d already exists, skipping\n"), i + 1);
		return FALSE;
	}
	return TRUE;
}

static int compute_mask_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_orig_image, MB_per_copies, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, 1., FALSE, &MB_per_orig_image, NULL, &MB_avail);

	/* The mask creation memory consumption is:
		* the original image of size O(w*h*l), with w/h the dimensions and l the number of layers
		* a 8b copy of one layer O(w*h as 8b) to feed opencv
		* a 32b copy of one layer at mask scale s, O(w*h*s^2 as 32b) out from opencv and saved
		* Internally in opencv, we also use (see opencv::cvDownscaleBlendMask):
		* a 8b copy of one layer at mask scale s, O((s*w+2)(s*h+2) as 8b)
		* a 32b copy of one layer at mask scale s, O((s*w+2)(s*h+2) as 32b)
	*/

	int seqrx = 0, seqry = 0, maskrx = 0, maskry = 0;
	size_t imgsize = get_max_seq_dimension(args->seq, &seqrx, &seqrx);
	compute_downscaled_mask_size(seqrx, seqry, &maskrx, &maskry, NULL, NULL);
	size_t memory_per_8b_orig_image = imgsize * sizeof(uint8_t);
	size_t memory_per_scaled_image = maskrx * maskry * sizeof(float);
	size_t memory_per_scaled_image_8b_opencv = (maskrx + 2) * (maskry + 2) * sizeof(uint8_t);
	size_t memory_per_scaled_image_32b_opencv = (maskrx + 2) * (maskry + 2) * sizeof(float);
	MB_per_copies = (memory_per_8b_orig_image + 
						memory_per_scaled_image + 
						memory_per_scaled_image_8b_opencv + 
						memory_per_scaled_image_32b_opencv) / 
						BYTES_IN_A_MB;

	unsigned int required = MB_per_orig_image + MB_per_copies;
	limit = MB_avail / required;

	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread) {
			thread_limit = com.max_thread;
		} else limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per thread, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);
		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, required, limit, for_writer ? "images" : "threads");
#else
		limit = 1;
#endif
	}
	return limit;
}

// this will create a downscaled 32b distance mask and save it directly (not using the save hook)
static int compute_mask_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	size_t nbpix = fit->naxes[0] * fit->naxes[1];
	int layer = (fit->naxes[2] == 3) ? GLAYER : RLAYER;
	uint8_t *buffer8in = calloc(nbpix, sizeof(uint8_t));
	float *buffer32out = NULL;

	if (!buffer8in) {
		PRINT_ALLOC_ERR;
		siril_debug_print("failed to allocate mask buffer\n");
		free(buffer8in);
		return 1;
	}
	// we convert the ref layer to 0 or 255 vals
	if (fit->type == DATA_FLOAT) {
		for (size_t i = 0; i < nbpix; i++) {
			if (fit->fpdata[layer][i] != 0.)
				buffer8in[i] = UCHAR_MAX;
		}
	} else if (fit->type == DATA_USHORT){
		for (size_t i = 0; i < nbpix; i++) {
			if (fit->pdata[layer][i] > 0)
				buffer8in[i] = UCHAR_MAX;
		}
	} else {
		siril_debug_print("wrong data type\n");
		free(buffer8in);
		return 1;
	}
	// we downscale the buffer
	int rx = fit->naxes[0];
	int ry = fit->naxes[1];
	int rx_out = 0, ry_out = 0;
	compute_downscaled_mask_size(rx, ry, &rx_out, &ry_out, NULL, NULL);
	buffer32out = calloc((size_t)(rx_out * ry_out), sizeof(float));
	if (!buffer32out) {
		PRINT_ALLOC_ERR;
		siril_debug_print("failed to allocate mask downscaled buffer\n");
		free(buffer8in);
		return 1;
	}
	cvDownscaleBlendMask(rx, ry, rx_out, ry_out, buffer8in, buffer32out);

	//we save the mask
	const gchar *mask_filename = get_cache_filename(args->seq, i, "msk", NULL);
	if (!mask_filename) {
		siril_debug_print("failed to create the mask filename");
		free(buffer8in);
		free(buffer32out);
		return 1;
	}
	if (save_mask_fits(rx_out, ry_out, buffer32out, mask_filename)) {
		siril_log_color_message(_("Failed to save mask for image %d\n"), "red", i + 1);
		free(buffer8in);
		free(buffer32out);
		return 1;
	}
	free(buffer8in);
	free(buffer32out);
	return 0;
}

gboolean end_compute_masks(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	free(args);
	return FALSE;
}

int compute_masks(struct stacking_args *args) {
	struct generic_seq_args *arg = create_default_seqargs(args->seq);
	arg->force_float = FALSE;
	arg->compute_mem_limits_hook = compute_mask_compute_mem_limits;
	arg->image_read_hook = compute_mask_read_hook;
	arg->filtering_criterion = args->filtering_criterion;
	arg->filtering_parameter = args->filtering_parameter;
	arg->nb_filtered_images = args->seq->selnum;
	arg->image_hook = compute_mask_image_hook;
	arg->description = _("Compute feathering masks");
	arg->has_output = FALSE;
	arg->already_in_a_thread = TRUE;
	arg->stop_on_error = TRUE;

	int retval = GPOINTER_TO_INT(generic_sequence_worker(arg));

	return retval;
}


