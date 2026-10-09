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

#include <string.h>

#include "core/gui_iface.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/proto.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "io/Astro-TIFF.h"
#include "io/avi_pipp/avi_writer.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/seqwriter.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#ifdef HAVE_FFMPEG
#include "io/mp4_output.h"
#endif
#include "algos/geometry.h"
#include "algos/siril_wcs.h"
#include "io/sequence_export.h"

/* film outputs expect as many channels as the sequence: mono frames must use
 * the gray sRGB profile, the RGB one would turn them into 3-channel images */
static void convert_to_srgb(fits *fit) {
	cmsHPROFILE srgb = (fit->naxes[2] == 1) ? gray_srgbtrc() : srgb_trc();
	siril_colorspace_transform(fit, srgb);
	cmsCloseProfile(srgb);
}

/* Used for avi exporter, creates buffer as BGRBGR from ushort FITS */
static uint8_t *fits_to_uint8(fits *fit) {
	size_t i, j;
	size_t n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	int channel = fit->naxes[2];
	int step = (channel == 3 ? 2 : 0);

	uint8_t *data = malloc(n);
	if (!data)
		return NULL;
	if (fit->orig_bitpix == BYTE_IMG) {
		for (i = 0, j = 0; i < n; i += channel, j++) {
			data[i + step] = (BYTE)fit->pdata[RLAYER][j];
			if (channel > 1) {
				data[i + 1] = (BYTE)fit->pdata[GLAYER][j];
				data[i + 2 - step] = (BYTE)fit->pdata[BLAYER][j];
			}
		}
	} else {
		siril_log_debug("converting from ushort\n");
		for (i = 0, j = 0; i < n; i += channel, j++) {
			data[i + step] = (BYTE)(fit->pdata[RLAYER][j] >> 8);
			if (channel > 1) {
				data[i + 1] = (BYTE)(fit->pdata[GLAYER][j] >> 8);
				data[i + 2 - step] = (BYTE)(fit->pdata[BLAYER][j] >> 8);
			}
		}
	}
	return data;
}

struct export_data {
	struct exportseq_args *ex;
	int output_bitpix;
	gboolean preserve_wcs;
	int reglayer;
	double dxref, dyref;
	unsigned int out_width, out_height;
	norm_coeff coeff;
	unsigned char *ref_icc;	// serialized: lcms profiles can't be shared between threads
	guint32 ref_icc_len;
	gint icc_msg_given;
	gchar *dest;	// single-file output, removed on abort
	fitseq *fitseq_file;
	struct ser_struct *ser_file;
	struct seqwriter_data *film_writer;	// keeps film frames ordered
	gboolean avi_opened;
#ifdef HAVE_FFMPEG
	struct mp4_struct *mp4_file;
#endif
};

static gboolean is_film(export_format output) {
	return output == EXPORT_AVI || output == EXPORT_MP4 ||
		output == EXPORT_MP4_H265 || output == EXPORT_WEBM_VP9;
}

static int get_output_bitpix(struct exportseq_args *ex) {
	int bitpix = ex->seq->bitpix;
	if (is_film(ex->output))
		return BYTE_IMG;
	// limit to 16 bits for TIFF and SER
	if ((ex->output == EXPORT_TIFF || ex->output == EXPORT_SER) && bitpix == FLOAT_IMG)
		return USHORT_IMG;
	if (bitpix == FLOAT_IMG && com.pref.force_16bit)
		return USHORT_IMG;
	return bitpix;
}

static int film_write_image(struct seqwriter_data *writer, fits *image, int index) {
	struct export_data *data = (struct export_data *)writer->sequence;
	if (data->ex->output == EXPORT_AVI) {
		uint8_t *buf = fits_to_uint8(image);
		if (!buf) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		int retval = avi_file_write_frame(0, buf);
		free(buf);
		return retval != 0;
	}
#ifdef HAVE_FFMPEG
	return mp4_add_frame(data->mp4_file, image) != 0;
#else
	return 1;
#endif
}

static int export_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	struct export_data *data = (struct export_data *)args->user;
	unsigned int MB_per_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, data->output_bitpix == FLOAT_IMG,
			NULL, &MB_per_image, &MB_avail);
	// the read frame and its shifted copy
	unsigned int required = 2 * MB_per_image;
	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;
		if (for_writer)
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_image;
		else limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		siril_log_error(_("%s: not enough memory to do this operation (%s required per thread, %s considered available)\n"),
				args->description, mem_per_thread, mem_available);
		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer && limit > com.max_thread * 3)
			limit = com.max_thread * 3;
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

static gint64 export_compute_size(struct generic_seq_args *args, int nb_frames) {
	struct export_data *data = (struct export_data *)args->user;
	export_format output = data->ex->output;
	if (output == EXPORT_MP4 || output == EXPORT_MP4_H265 || output == EXPORT_WEBM_VP9)
		return 1;	// compressed, size unknown
	gint64 frame_size = (gint64)data->out_width * data->out_height * args->seq->nb_layers;
	if (data->output_bitpix == FLOAT_IMG)
		frame_size *= sizeof(float);
	else if (data->output_bitpix != BYTE_IMG)
		frame_size *= sizeof(WORD);
	return frame_size * nb_frames;
}

static unsigned char *get_ref_icc(sequence *seq, int refindex, guint32 *len) {
	fits ref = { 0 };
	unsigned char *icc = NULL;
	if (!seq_read_frame(seq, refindex, &ref, FALSE, -1) && ref.icc_profile)
		icc = get_icc_profile_data(ref.icc_profile, len);
	clearfits(&ref);
	return icc;
}

static int export_prepare(struct generic_seq_args *args) {
	struct export_data *data = (struct export_data *)args->user;
	struct exportseq_args *ex = data->ex;
	unsigned int in_width, in_height;

	data->reglayer = get_registration_layer(args->seq);
	if (data->reglayer >= 0)
		siril_log_message(_("Using registration information from layer %d to export sequence\n"), data->reglayer);
	if (ex->crop) {
		in_width  = ex->crop_area.w;
		in_height = ex->crop_area.h;
	} else {
		in_width  = args->seq->rx;
		in_height = args->seq->ry;
	}

	if (ex->resample) {
		data->out_width = ex->dest_width;
		data->out_height = ex->dest_height;
		if (data->out_width == in_width && data->out_height == in_height)
			ex->resample = FALSE;
	} else {
		data->out_width = in_width;
		data->out_height = in_height;
	}
	data->preserve_wcs = (ex->output == EXPORT_FITS || ex->output == EXPORT_FITSEQ) && !ex->resample;

	gchar *filter_descr = describe_filter(args->seq, ex->filtering_criterion, ex->filtering_parameter);
	siril_log_message(filter_descr);
	g_free(filter_descr);

	int refindex = sequence_find_refimage(args->seq);
	if (data->reglayer != -1 && args->seq->regparam[data->reglayer]) {
		if (!test_regdata_is_valid_and_shift(args->seq, data->reglayer)) {
			siril_log_error(_("Export has detected registration data with more than simple shifts, this is not supported\n"));
			return 1;
		}
		translation_from_H(args->seq->regparam[data->reglayer][refindex].H, &data->dxref, &data->dyref);
	} else {
		data->reglayer = -1;
	}

	if (ex->normalize) {
		struct stacking_args stackargs = { 0 };
		stackargs.force_norm = FALSE;
		stackargs.seq = args->seq;
		stackargs.filtering_criterion = ex->filtering_criterion;
		stackargs.filtering_parameter = ex->filtering_parameter;
		stackargs.nb_images_to_stack = args->nb_filtered_images;
		stackargs.normalize = ADDITIVE_SCALING;
		stackargs.reglayer = data->reglayer;
		stackargs.use_32bit_output = (data->output_bitpix == FLOAT_IMG);
		stackargs.ref_image = args->seq->reference_image;
		stackargs.equalizeRGB = FALSE;
		stackargs.lite_norm = FALSE;

		// build image indices used by normalization
		if (stack_fill_list_of_unfiltered_images(&stackargs))
			return 1;
		int retval = do_normalization(&stackargs);
		free(stackargs.image_indices);
		free(stackargs.coeff.mul);
		data->coeff = stackargs.coeff;
		data->coeff.mul = NULL;
		if (retval)
			return 1;
	}

	if (args->seqwriter == SEQWRITER_ALWAYS) {
		int limit = export_compute_mem_limits(args, TRUE);
		if (limit == 0)
			return 1;
		seqwriter_set_max_active_blocks(limit);
	}

	/* possible output formats: FITS images, FITS cube, TIFF, SER, AVI, MP4, WEBM */
	// create the sequence file for single-file sequence formats
	switch (ex->output) {
		case EXPORT_FITS:
			data->ref_icc = get_ref_icc(args->seq, refindex, &data->ref_icc_len);
			if (data->ref_icc)
				siril_log_message(_("Reference frame has an ICC profile. Will assign / convert other frames to to match.\n"));
			break;
		case EXPORT_FITSEQ:
			data->dest = g_strdup_printf("%s%s", ex->basename, com.pref.ext);
			data->fitseq_file = calloc(1, sizeof(fitseq));
			if (!data->fitseq_file || fitseq_create_file(data->dest, data->fitseq_file, -1)) {
				free(data->fitseq_file);
				data->fitseq_file = NULL;
				return 1;
			}
			break;
		case EXPORT_TIFF:
#ifndef HAVE_LIBTIFF
			siril_log_error(_("TIFF output is not supported because siril was not compiled with libtiff support.\n"));
			return 1;
#endif
			break;
		case EXPORT_SER:
			data->dest = g_strdup_printf("%s.ser", ex->basename);
			data->ser_file = calloc(1, sizeof(struct ser_struct));
			if (!data->ser_file || ser_create_file(data->dest, data->ser_file, TRUE, args->seq->ser_file)) {
				free(data->ser_file);
				data->ser_file = NULL;
				return 1;
			}
			break;
		case EXPORT_AVI:
			// Check if the sequence has an ICC profile. If so, we should convert to sRGB
			// as that's really the only suitable option here
			data->ref_icc = get_ref_icc(args->seq, refindex, &data->ref_icc_len);
			if (data->ref_icc)
				siril_log_message(_("Reference frame has an ICC profile. Exporting as sRGB.\n"));

			data->dest = g_strdup_printf("%s.avi", ex->basename);
			if (avi_file_create(data->dest, data->out_width, data->out_height,
						args->seq->nb_layers == 1 ? AVI_WRITER_INPUT_FORMAT_MONOCHROME : AVI_WRITER_INPUT_FORMAT_COLOUR,
						AVI_WRITER_CODEC_DIB, ex->film_fps)) {
				siril_log_error(_("AVI file `%s' could not be created\n"), data->dest);
				return 1;
			}
			data->avi_opened = TRUE;
			break;
		case EXPORT_MP4:
		case EXPORT_MP4_H265:
		case EXPORT_WEBM_VP9:
#ifndef HAVE_FFMPEG
			siril_log_message(_("MP4 output is not supported because siril was not compiled with ffmpeg support.\n"));
			return 1;
#else
			data->ref_icc = get_ref_icc(args->seq, refindex, &data->ref_icc_len);
			if (data->ref_icc)
				siril_log_message(_("Reference frame has an ICC profile. Exporting as sRGB.\n"));

			/* resampling is managed by libswscale */
			data->dest = g_strdup_printf("%s.%s", ex->basename,
					ex->output == EXPORT_WEBM_VP9 ? "webm" : "mp4");

			if (in_width % 32 || data->out_height % 2 || data->out_width % 2) {
				siril_log_message(_("Film output needs to have a width that is a multiple of 32 and an even height, resizing selection.\n"));
				if (in_width % 32) in_width = (in_width / 32) * 32 + 32;
				if (in_height % 2) in_height++;
				if (ex->crop) {
					ex->crop_area.w = in_width;
					ex->crop_area.h = in_height;
				} else {
					ex->crop = TRUE;
					ex->crop_area.x = 0;
					ex->crop_area.y = 0;
					ex->crop_area.w = in_width;
					ex->crop_area.h = in_height;
				}
				compute_fitting_selection(&ex->crop_area, 32, 2, 0);
				memcpy(&com.selection, &ex->crop_area, sizeof(rectangle));
				fprintf(stdout, "final input area: %d,%d,\t%dx%d\n",
						ex->crop_area.x, ex->crop_area.y,
						ex->crop_area.w, ex->crop_area.h);
				in_width = ex->crop_area.w;
				in_height = ex->crop_area.h;
				if (!ex->resample) {
					data->out_width = in_width;
					data->out_height = in_height;
				} else {
					if (data->out_width % 2) data->out_width++;
					if (data->out_height % 2) data->out_height++;
				}
			}

			data->mp4_file = mp4_create(data->dest, data->out_width, data->out_height, ex->film_fps,
					args->seq->nb_layers, ex->film_quality, in_width, in_height, ex->output);
			if (!data->mp4_file)
				return 1;
			break;
#endif
	}

	if (is_film(ex->output)) {
		data->film_writer = calloc(1, sizeof(struct seqwriter_data));
		if (!data->film_writer) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		data->film_writer->write_image_hook = film_write_image;
		data->film_writer->sequence = data;
		data->film_writer->output_type = SEQ_AVI;
		start_writer(data->film_writer, -1);
	}
	return 0;
}

/* shifts the frame by the registration data and normalizes it, replacing its
 * buffer by one of the output data type */
static int shift_and_normalize(fits *fit, struct export_data *data, int out_index,
		int shiftx, int shifty, int threads) {
	data_type out_type = data->output_bitpix == FLOAT_IMG ? DATA_FLOAT : DATA_USHORT;
	gboolean normalize = data->ex->normalize;
	size_t nbpix = (size_t)fit->rx * fit->ry;
	void *out = calloc(nbpix * fit->naxes[2], out_type == DATA_FLOAT ? sizeof(float) : sizeof(WORD));
	if (!out) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	WORD *wout = (WORD *)out;
	float *fout = (float *)out;
	int x0 = max(0, -shiftx), x1 = min(fit->rx, fit->rx - shiftx);
	int y0 = max(0, -shifty), y1 = min(fit->ry, fit->ry - shifty);

	for (int layer = 0; layer < fit->naxes[2]; ++layer) {
		double scale = normalize ? data->coeff.pscale[layer][out_index] : 1.0;
		double offset = normalize ? data->coeff.poffset[layer][out_index] : 0.0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
#endif
		for (int y = y0; y < y1; ++y) {
			size_t src = (size_t)y * fit->rx;
			size_t dst = layer * nbpix + (size_t)(y + shifty) * fit->rx + shiftx;
			if (fit->type == DATA_USHORT) {
				for (int x = x0; x < x1; ++x) {
					WORD pixel = fit->pdata[layer][src + x];
					if (normalize && pixel > 0) // do not offset null pixels
						pixel = round_to_WORD(pixel * scale - offset);
					if (out_type == DATA_FLOAT)
						fout[dst + x] = pixel / USHRT_MAX_SINGLE;
					else wout[dst + x] = pixel;
				}
			} else {
				for (int x = x0; x < x1; ++x) {
					float pixel = fit->fpdata[layer][src + x];
					if (normalize && pixel != 0.f) { // do not offset null pixels
						pixel *= (float) scale;
						pixel -= (float) offset;
					}
					if (out_type == DATA_FLOAT)
						fout[dst + x] = pixel;
					else wout[dst + x] = roundf_to_WORD(pixel * USHRT_MAX_SINGLE);
				}
			}
		}
	}

	// for 8 bit output, the image will be transformed later with a linear scale
	int orig_bitpix = fit->orig_bitpix;
	if (out_type == DATA_FLOAT)
		free(fit->fdata);
	else free(fit->data);
	fit_replace_buffer(fit, out, out_type);
	fit->bitpix = data->output_bitpix;
	fit->orig_bitpix = orig_bitpix;
	return 0;
}

static int export_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct export_data *data = (struct export_data *)args->user;
	struct exportseq_args *ex = data->ex;

	if (fit->rx != args->seq->rx || fit->ry != args->seq->ry) {
		siril_log_error(_("An image of the sequence doesn't have the same dimensions\n"));
		return 1;
	}

	int shiftx = 0, shifty = 0;
	if (data->reglayer != -1) {
		double dx, dy;
		translation_from_H(args->seq->regparam[data->reglayer][i].H, &dx, &dy);
		shiftx = round_to_int(dx - data->dxref);
		shifty = round_to_int(dy - data->dyref);
		if (has_wcs(fit)) {
			if (data->preserve_wcs) {
				Homography H = { 0 };
				cvGetEye(&H);
				H.h02 = (double)shiftx;
				H.h12 = -(double)shifty;
				cvApplyFlips(&H, fit->ry, fit->ry);
				reframe_wcs(fit->keywords.wcslib, &H);
			} else
				free_wcs(fit);
		}
	}

	if (shift_and_normalize(fit, data, o, shiftx, shifty, threads))
		return 1;

	if (ex->crop) {
		rectangle area = ex->crop_area;	// crop() may round it
		if (crop(fit, &area))
			return 1;
	}

	cmsHPROFILE ref_icc = NULL;
	if (data->ref_icc)
		ref_icc = cmsOpenProfileFromMem(data->ref_icc, data->ref_icc_len);
	// frames without ICC profile get the one of the reference frame
	if (!fit->icc_profile && ref_icc)
		fit->icc_profile = copyICCProfile(ref_icc);
	color_manage(fit, fit->icc_profile != NULL);

	if (ex->output == EXPORT_FITS) {
		if (ref_icc)
			siril_colorspace_transform(fit, ref_icc);
		else if (fit->icc_profile && g_atomic_int_compare_and_exchange(&data->icc_msg_given, FALSE, TRUE))
			siril_log_message(_("Info: this frame has an ICC profile but the reference frame does not. Profile will be preserved...\n"));
	}
	else if (is_film(ex->output) && ref_icc)
		convert_to_srgb(fit);
	if (ref_icc)
		cmsCloseProfile(ref_icc);
	return 0;
}

/* frames are written in sequence order by a seqwriter for single-file outputs,
 * fit is NULL for a failed frame */
static int export_save_hook(struct generic_seq_args *args, int o, int i, fits *fit) {
	struct export_data *data = (struct export_data *)args->user;
	struct exportseq_args *ex = data->ex;
	int retval = 1;
	gchar *dest;

	switch (ex->output) {
		case EXPORT_FITS:
			dest = g_strdup_printf("%s%05d%s", ex->basename, i + 1, com.pref.ext);
			retval = savefits(dest, fit);
			g_free(dest);
			break;
		case EXPORT_TIFF:
#ifdef HAVE_LIBTIFF
			dest = g_strdup_printf("%s%05d", ex->basename, i + 1);
			gchar *astro_tiff = AstroTiff_build_header(fit);
			retval = savetif(dest, fit, 16, astro_tiff, com.pref.copyright, ex->tiff_compression, TRUE, TRUE);
			g_free(astro_tiff);
			g_free(dest);
#endif
			break;
		case EXPORT_FITSEQ:
			retval = fitseq_write_image(data->fitseq_file, fit, o);
			break;
		case EXPORT_SER:
			retval = ser_write_frame_from_fit(data->ser_file, fit, o);
			break;
		default:
			retval = seqwriter_append_write(data->film_writer, fit, o);
	}
	return retval;
}

static int export_finalize(struct generic_seq_args *args) {
	struct export_data *data = (struct export_data *)args->user;
	struct exportseq_args *ex = data->ex;
	gboolean aborted = !processing_should_continue();
	int retval = 0;

	if (data->fitseq_file) {
		retval = fitseq_close_file(data->fitseq_file);
		free(data->fitseq_file);
	}
	if (data->ser_file) {
		retval = ser_write_and_close(data->ser_file);
		free(data->ser_file);
	}
	if (data->film_writer) {
		retval = stop_writer(data->film_writer, aborted);
		free(data->film_writer);
	}
	if (data->avi_opened)
		avi_file_close(0);
#ifdef HAVE_FFMPEG
	if (data->mp4_file) {
		mp4_close(data->mp4_file, aborted);
		free(data->mp4_file);
	}
#endif
	// disposing of the file if Stop button was hit
	if (aborted && data->dest && g_unlink(data->dest))
		siril_log_debug("Failed to delete %s\n", data->dest);

	free(data->ref_icc);
	free(data->coeff.offset);
	free(data->coeff.scale);
	g_free(data->dest);
	free(ex->basename);
	free(ex);
	free(data);
	args->user = NULL;
	return retval;
}

gboolean sequence_export_start(struct exportseq_args *ex) {
	struct export_data *data = calloc(1, sizeof(struct export_data));
	struct generic_seq_args *args = create_default_seqargs(ex->seq);
	if (!data || !args) {
		PRINT_ALLOC_ERR;
		free(data);
		free(args);
		free(ex->basename);
		free(ex);
		return FALSE;
	}
	data->ex = ex;
	data->output_bitpix = get_output_bitpix(ex);

	args->filtering_criterion = ex->filtering_criterion;
	args->filtering_parameter = ex->filtering_parameter;
	args->nb_filtered_images = compute_nb_filtered_images(ex->seq,
			ex->filtering_criterion, ex->filtering_parameter);
	args->compute_mem_limits_hook = export_compute_mem_limits;
	args->compute_size_hook = export_compute_size;
	args->prepare_hook = export_prepare;
	args->image_hook = export_image_hook;
	args->save_hook = export_save_hook;
	args->finalize_hook = export_finalize;
	args->has_output = TRUE;
	args->output_type = data->output_bitpix == FLOAT_IMG ? DATA_FLOAT : DATA_USHORT;
	args->seqwriter = (ex->output == EXPORT_FITS || ex->output == EXPORT_TIFF) ?
		SEQWRITER_NONE : SEQWRITER_ALWAYS;
	// the generic check is on the input sequence, not on the FITS output
	args->parallel = ex->output != EXPORT_FITS || fits_is_reentrant();
	args->description = _("Sequence export");
	args->user = data;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(data);
		free_generic_seq_args(args, FALSE);
		free(ex->basename);
		free(ex);
		return FALSE;
	}
	return TRUE;
}
