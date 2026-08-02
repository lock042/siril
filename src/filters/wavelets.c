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
#include <stdlib.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/gui_iface.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "algos/Def_Wavelet.h"
#include "algos/wavelet_denoise.h"
#include "wavelets.h"
#include "core/op_descriptors.h"

/* Op descriptors — single source of truth for the wavelet ops. The wrecons
 * command / GUI-apply / preview sites use different progress labels, kept as
 * per-site description overrides. */
const op_descriptor op_desc_wrecons = {
	.id = "wavelets.wrecons", .version = 1,
	.image_hook = wrecons_image_hook,
	.log_hook = wrecons_log_hook,
	.description = N_("Wavelet reconstruction"),
	.mem_ratio = 0.0f,
	.flags = 0,
};

const op_descriptor op_desc_atrous = {
	.id = "wavelets.atrous", .version = 1,
	.image_hook = atrous_image_hook,
	.log_hook = atrous_log_hook,
	.description = N_("Wavelet transform"),
	.mem_ratio = 0.0f,
	.flags = 0,
};

/************* wavelet transform worker (wavelet command and GUI compute path) *************/

gpointer wavelet_transform_worker(gpointer p) {
	struct wavelet_transform_data *args = (struct wavelet_transform_data *)p;
	const char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" };
	const char *tmpdir = g_get_tmp_dir();
	int retval = 0;
	g_rw_lock_reader_lock(&gfit->rwlock);
	int nb_chan = gfit->naxes[2];

	if (gfit->type == DATA_USHORT) {
		float *Imag = f_vector_alloc((size_t) gfit->rx * (size_t) gfit->ry);
		if (!Imag) {
			PRINT_ALLOC_ERR;
			retval = 1;
		} else {
			for (int i = 0; i < nb_chan; i++) {
				gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
				wavelet_transform_file(Imag, gfit->ry, gfit->rx, dir,
						args->Type_Transform, args->Nbr_Plan, gfit->pdata[i], args->anscombe);
				g_free(dir);
			}
			free(Imag);
		}
	} else if (gfit->type == DATA_FLOAT) {
		for (int i = 0; i < nb_chan; i++) {
			gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
			wavelet_transform_file_float(gfit->fpdata[i], gfit->ry, gfit->rx, dir,
					args->Type_Transform, args->Nbr_Plan, args->anscombe);
			g_free(dir);
		}
	} else {
		retval = 1;
	}
	if (!retval)
		siril_log_message(_("Wavelet decomposition computed (%d plans)\n"), args->Nbr_Plan);
	g_rw_lock_reader_unlock(&gfit->rwlock);

	/* The .wave files are now written and gfit is unlocked. Hand off to the
	 * caller's completion idle (the GUI re-enables its widgets there); it owns
	 * args. Without one (command line, or on failure) free args and run the
	 * generic idle. */
	if (!retval && args->idle) {
		siril_add_idle(args->idle, args);
	} else {
		free(args);
		siril_add_idle(end_generic, NULL);
	}
	return GINT_TO_POINTER(retval);
}

/************* wrecons hook (shared with GUI OK path and process_wrecons) *************/

void free_wrecons_data(void *p) {
	free(p);
}

int wrecons_image_hook(struct generic_img_args *gargs, fits *fit, int threads) {
	struct wrecons_data *args = (struct wrecons_data *)gargs->user;
	const char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave", "b_rawdata.wave" };
	const char *tmpdir = g_get_tmp_dir();
	for (int i = 0; i < args->nb_chan; i++) {
		gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		int ret;
		if (args->for_roi) {
			/* Reconstruct only the selection: reads just the ROI rows of the
			 * transform and denoises/reconstructs that small window. */
			ret = wavelet_reconstruct_file_roi(dir, args->coef, &args->denoise,
					args->roi_x, args->roi_y, fit->rx, fit->ry, i, fit);
		} else if (fit->type == DATA_USHORT) {
			ret = wavelet_reconstruct_file(dir, args->coef, &args->denoise, fit->pdata[i]);
		} else if (fit->type == DATA_FLOAT) {
			ret = wavelet_reconstruct_file_float(dir, args->coef, &args->denoise, fit->fpdata[i]);
		} else {
			g_free(dir);
			return 1;
		}
		g_free(dir);
		if (ret) return 1;
	}
	return 0;
}

gchar *wrecons_log_hook(gpointer p, log_hook_detail detail) {
	return g_strdup(_("Wavelet reconstruction"));
}

/************* atrous: one-shot decompose + denoise + reconstruct *************/

void free_atrous_data(void *p) {
	struct atrous_data *args = (struct atrous_data *)p;
	if (!args)
		return;
	free(args->seqEntry);
	free(args);
}

/* The decomposition produced here is consumed by the reconstruction a few lines
 * below and then thrown away, so it is kept in memory rather than round-tripped
 * through a .wave file in the temporary directory: for a 24 Mpx frame at 6
 * layers that file is ~550 MB per channel, written and read back for nothing.
 * Keeping it in memory also makes this function re-entrant, which is what lets
 * seqatrous process several frames at once. */
int atrous_transform_image(fits *fit, const struct atrous_data *args) {
	const int Nl = fit->ry, Nc = fit->rx;
	const size_t n = (size_t) Nl * (size_t) Nc;
	const int nb_chan = fit->naxes[2];
	int retval = 0;
	/* per-layer weights are read (not written) by the reconstruction, but the
	 * API takes a non-const pointer, so work from a local copy */
	float coef[7];
	memcpy(coef, args->coef, sizeof coef);

	if (fit->type != DATA_USHORT && fit->type != DATA_FLOAT)
		return 1;

	/* USHORT works in a float scratch buffer; float decomposes straight from
	 * the image buffer, except under Anscombe where the variance-stabilised
	 * copy must not clobber the original pixels before the transform succeeds. */
	float *Imag = NULL;
	if (fit->type == DATA_USHORT || args->anscombe) {
		Imag = malloc(n * sizeof(float));
		if (!Imag) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	const double ans_scale = (fit->type == DATA_USHORT) ? ANSCOMBE_USHORT_SCALE
			: ANSCOMBE_FLOAT_SCALE;

	for (int i = 0; i < nb_chan; i++) {
		wave_transf_des wavelet = { 0 };
		float *src, *dst;

		if (fit->type == DATA_USHORT) {
			prepare_rawdata(Imag, Nl, Nc, fit->pdata[i]);
			if (args->anscombe)
				anscombe_forward(Imag, n, ans_scale);
			src = Imag;
			dst = Imag; /* reconstructed into the same scratch, then requantised */
		} else if (args->anscombe) {
			memcpy(Imag, fit->fpdata[i], n * sizeof(float));
			anscombe_forward(Imag, n, ans_scale);
			src = Imag;
			dst = fit->fpdata[i];
		} else {
			src = fit->fpdata[i];
			dst = fit->fpdata[i];
		}

		/* pave_2d_tfo copies its input before recursing, so src is intact here
		 * and dst may safely alias it */
		if (wavelet_transform_data(src, Nl, Nc, &wavelet, args->type,
				args->nbr_plan)) {
			retval = 1;
			break;
		}
		wavelet_denoise_planes(wavelet.Pave.Data, args->type, args->nbr_plan,
				Nl, Nc, &args->denoise);
		retval = wavelet_reconstruct_data(&wavelet, dst, coef);
		wave_io_free(&wavelet);
		if (retval)
			break;

		if (args->anscombe)
			anscombe_inverse(dst, n, ans_scale);
		if (fit->type == DATA_USHORT)
			reget_rawdata(dst, Nl, Nc, fit->pdata[i]);
	}
	free(Imag);
	return retval;
}

int atrous_image_hook(struct generic_img_args *gargs, fits *fit, int threads) {
	struct atrous_data *args = (struct atrous_data *)gargs->user;
	return atrous_transform_image(fit, args);
}

gchar *atrous_log_hook(gpointer p, log_hook_detail detail) {
	struct atrous_data *args = (struct atrous_data *)p;
	if (args->denoise.enabled)
		return g_strdup_printf(_("À trous wavelet transform (%d layers, denoise k=%.2f)"),
				args->nbr_plan, args->denoise.k);
	return g_strdup_printf(_("À trous wavelet transform (%d layers)"), args->nbr_plan);
}

/* Sequence hook: decompose + reconstruct one frame in place. */
static int atrous_seq_image_hook(struct generic_seq_args *args, int out_index,
		int in_index, fits *fit, rectangle *_, int threads) {
	struct atrous_data *a_args = (struct atrous_data *)args->user;
	return atrous_transform_image(fit, a_args);
}

/* Memory needed per frame, on top of the frame itself:
 *   the transform cube          -> nbr_plan planes of one channel
 *   pave_2d_tfo's two scratch planes
 *   the USHORT/Anscombe scratch -> one plane
 *   bivariate shrinkage         -> two more planes while it runs
 * all in float, all sized on a single channel. */
static int atrous_mem_limits_hook(struct generic_seq_args *args, gboolean for_writer) {
	struct atrous_data *a_args = (struct atrous_data *)args->user;
	unsigned int MB_per_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image,
			NULL, &MB_avail);
	unsigned int required = MB_per_image;

	if (limit > 0) {
		/* the transform runs on one channel at a time, whatever the frame's
		 * channel count, so the scratch is sized on a single channel */
		const guint64 chan_bytes = (guint64) args->seq->rx * args->seq->ry
				* sizeof(float);
		unsigned int MB_per_channel_float = max(1,
				(unsigned int) (chan_bytes / BYTES_IN_A_MB));
		int planes = a_args->nbr_plan + 2; /* cube + the two tfo scratch planes */
		if (a_args->anscombe || get_data_type(args->seq->bitpix) == DATA_USHORT)
			planes += 1;
		if (a_args->denoise.enabled && a_args->denoise.method == WD_BISHRINK)
			planes += 2;

		required = MB_per_image + planes * MB_per_channel_float;

		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* the already-accounted processing images, plus whatever else fits
			 * in what the main computation leaves unused */
			limit = thread_limit
					+ (MB_avail - required * thread_limit) / MB_per_image;
		} else
			limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB,
				G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB,
				G_FORMAT_SIZE_IEC_UNITS);

		siril_log_error(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_log_debug("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

static int atrous_finalize_hook(struct generic_seq_args *args) {
	struct atrous_data *a_args = (struct atrous_data *)args->user;
	int retval = seq_finalize_hook(args);
	free_atrous_data(a_args);
	return retval;
}

void apply_atrous_to_sequence(struct atrous_data *a_args) {
	struct generic_seq_args *args = create_default_seqargs(a_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = a_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = atrous_finalize_hook;
	args->image_hook = atrous_seq_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Wavelet transform");
	args->has_output = TRUE;
	args->output_type = get_data_type(a_args->seq->bitpix);
	args->new_seq_prefix = strdup(a_args->seqEntry);
	args->load_new_sequence = TRUE;
	args->user = a_args;
	/* The transform is held in memory, so frames can be processed in parallel
	 * up to whatever the transform cube plus scratch allows. */
	args->compute_mem_limits_hook = atrous_mem_limits_hook;

	a_args->fit = NULL; /* not used in sequence mode */

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free_atrous_data(a_args);
		free_generic_seq_args(args, TRUE); /* frees args->new_seq_prefix */
	}
}

/* This function computes wavelets with the number of Nbr_Plan and
 * extracts plan "Plan" in fit parameters */

int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer) {
	int chan, start, end, retval = 0;
	wave_transf_des wavelet[3] = { 0 };

	g_assert(fit->naxes[2] <= 3);

	float *Imag = NULL;
	size_t n = fit->naxes[0] * fit->naxes[1];
	if (fit->type == DATA_USHORT) {
		Imag = f_vector_alloc(n);
		if (!Imag) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	if (reqlayer < 0 || reqlayer >= 3) {
		start = 0;
		end = fit->naxes[2];
	}
	else {
		start = reqlayer;
		end = start + 1;
	}

	for (chan = start; chan < end; chan++) {
		int Nl, Nc;

		if (fit->type == DATA_USHORT) {
			/* float wavelet of data [0, 65535] */
			if (wavelet_transform(Imag, fit->ry, fit->rx, &wavelet[chan],
						Type, Nbr_Plan, fit->pdata[chan])) {
				retval = 1;
				break;
			}
		}
		else if (fit->type == DATA_FLOAT) {
			/* float wavelet of data [0, 1] */
			Imag = fit->fpdata[chan];
			if (wavelet_transform_float(Imag, fit->ry, fit->rx, &wavelet[chan],
						Type, Nbr_Plan)) {
				retval = 1;
				break;
			}
		} else { // Unknown fit->type
			retval = 1;
			break;
		}
		Nl = wavelet[chan].Nbr_Ligne;
		Nc = wavelet[chan].Nbr_Col;
		pave_2d_extract_plan(wavelet[chan].Pave.Data, Imag, Nl, Nc, Plan);
		if (fit->type == DATA_USHORT)
			reget_rawdata(Imag, Nl, Nc, fit->pdata[chan]);
		wave_io_free(&wavelet[chan]);
	}

	/* Free */
	if (fit->type == DATA_USHORT)
		free(Imag);
	return retval;
}

static gboolean end_wavelets_filter(gpointer p) {
	struct wavelets_filter_data *args = (struct wavelets_filter_data *) p;
	stop_processing_thread();
	gui_iface.set_progress(PROGRESS_DONE, PROGRESS_TEXT_RESET);
	gui_iface.set_busy(FALSE);
	free(args);
	return FALSE;
}

gpointer extract_plans(gpointer p) {
	int i;
	fits fit = { 0 };
	struct wavelets_filter_data *args = (struct wavelets_filter_data *) p;

	gui_iface.set_progress(PROGRESS_RESET, NULL);

	for (i = 0; i < args->Nbr_Plan; i++) {
		gchar *filename, *msg;
		/* Extracted layers are written to disk for later recombination
		 * (e.g. via PixelMath); keep full metadata so a recombined result
		 * can retain the source's WCS/header. Geometry is unchanged, so the
		 * solution stays valid. */
		if (copyfits(args->fit, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT | CP_METADATA_HEAP, -1)) {
			siril_log_error(_("Could not copy image, aborting\n"));
			siril_add_idle(end_wavelets_filter, args);
			return GINT_TO_POINTER(1);
		}
		filename = g_strdup_printf("layer%02d", i);
		msg = g_strdup_printf(_("Extracting %s..."), filename);
		gui_iface.set_progress((double)i / args->Nbr_Plan, msg);
		get_wavelet_layers(&fit, args->Nbr_Plan, i, args->Type, -1);
		savefits(filename, &fit);
		g_free(filename);
		g_free(msg);
	}
	clearfits(&fit);
	siril_add_idle(end_wavelets_filter, args);
	return GINT_TO_POINTER(0);
}

/* GUI callbacks moved to src/gui/wavelets.c */
