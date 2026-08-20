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
#include "core/nde_history.h"
#include "core/OS_utils.h"

/* Op descriptors — single source of truth for the wavelet ops. The wrecons
 * command / GUI-apply / preview sites use different progress labels, kept as
 * per-site description overrides. */

/* PERMANENT Tier B (maintainer verdict, NDE phase 4.5): wavelets.wrecons keeps
 * NO serializer BY DESIGN.  The hook reads a decomposition transform file
 * (r/g/b_rawdata.wave in the tmpdir) produced by a separate `wavelet`
 * decomposition step that is not itself a descriptor op — the pixel output is
 * not determined by wrecons_data alone, and that tmpdir multi-file lifecycle
 * cannot live in a params blob.  Checkpoints are the intended editability UX.
 * Do NOT add serialize/deserialize here — a future coverage sweep must not
 * "fix" this. */
const op_descriptor op_desc_wrecons = {
	.id = "wavelets.wrecons", .version = 1,
	.image_hook = wrecons_image_hook,
	.log_hook = wrecons_log_hook,
	.description = N_("Wavelet reconstruction"),
	.mem_ratio = 0.0f,
	/* ROI-capable in a different way from the pixel-local ops: the hook
	 * cannot compute from the crop it is handed, because its input is a
	 * decomposition of the whole image — but it CAN produce just that
	 * window, reading only the transform rows the window covers
	 * (wavelet_reconstruct_*_roi).  The rectangle comes from
	 * generic_img_args.roi_rect; full_rx/full_ry pin the geometry the
	 * transform was computed at. */
	.flags = OP_ROI_CAPABLE,
};

/* NDE serializers for wavelets.atrous.  atrous_transform_image is a
 * self-contained decompose+denoise+reconstruct that reads nbr_plan, type,
 * anscombe, the per-layer coef[7] and the denoise sub-struct; seq/fit/seqEntry
 * are runtime context and are skipped. */
static gchar *atrous_serialize(gconstpointer user) {
	const struct atrous_data *d = user;
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "nbr_plan", d->nbr_plan);
	nde_kv_add_int(kv, "type", d->type);
	nde_kv_add_bool(kv, "anscombe", d->anscombe);
	for (int i = 0; i < 7; i++) {
		gchar key[8];
		g_snprintf(key, sizeof key, "coef%d", i);
		nde_kv_add_float(kv, key, d->coef[i]);
	}
	/* denoise sub-struct (all POD) */
	nde_kv_add_bool(kv, "dn_enabled", d->denoise.enabled);
	/* on-disk value: enum order is frozen by the NDE format — do not reorder */
	nde_kv_add_int(kv, "dn_method", d->denoise.method);
	nde_kv_add_float(kv, "dn_k", d->denoise.k);
	for (int i = 0; i < WD_MAX_PLAN; i++) {
		gchar key[12];
		g_snprintf(key, sizeof key, "dn_f%d", i);
		nde_kv_add_float(kv, key, d->denoise.f[i]);
	}
	/* on-disk value: enum order is frozen by the NDE format — do not reorder */
	nde_kv_add_int(kv, "dn_sigma_source", d->denoise.sigma_source);
	nde_kv_add_bool(kv, "dn_soft", d->denoise.soft);
	nde_kv_add_bool(kv, "dn_anscombe", d->denoise.anscombe);
	return nde_kv_end(kv);
}

static gpointer atrous_deserialize(const gchar *blob, int version) {
	if (version > op_desc_atrous.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	gint64 nbr_plan, type, dn_method, dn_sigma_source;
	gboolean anscombe, dn_enabled, dn_soft, dn_anscombe;
	float coef[7], dn_k, dn_f[WD_MAX_PLAN];
	gboolean ok = nde_kv_get_int(kv, "nbr_plan", &nbr_plan) &&
	              nde_kv_get_int(kv, "type", &type) &&
	              nde_kv_get_bool(kv, "anscombe", &anscombe) &&
	              nde_kv_get_bool(kv, "dn_enabled", &dn_enabled) &&
	              nde_kv_get_int(kv, "dn_method", &dn_method) &&
	              nde_kv_get_float(kv, "dn_k", &dn_k) &&
	              nde_kv_get_int(kv, "dn_sigma_source", &dn_sigma_source) &&
	              nde_kv_get_bool(kv, "dn_soft", &dn_soft) &&
	              nde_kv_get_bool(kv, "dn_anscombe", &dn_anscombe);
	for (int i = 0; ok && i < 7; i++) {
		gchar key[8];
		g_snprintf(key, sizeof key, "coef%d", i);
		ok = nde_kv_get_float(kv, key, &coef[i]);
	}
	for (int i = 0; ok && i < WD_MAX_PLAN; i++) {
		gchar key[12];
		g_snprintf(key, sizeof key, "dn_f%d", i);
		ok = nde_kv_get_float(kv, key, &dn_f[i]);
	}
	NDE_KV_REQUIRE(kv, ok);
	struct atrous_data *d = calloc(1, sizeof(*d));
	if (d) {
		d->destroy_fn = free;
		d->nbr_plan = (int)nbr_plan;
		d->type = (int)type;
		d->anscombe = anscombe;
		memcpy(d->coef, coef, sizeof coef);
		d->denoise.enabled = dn_enabled;
		d->denoise.method = (int)dn_method;
		d->denoise.k = dn_k;
		memcpy(d->denoise.f, dn_f, sizeof dn_f);
		d->denoise.sigma_source = (int)dn_sigma_source;
		d->denoise.soft = dn_soft;
		d->denoise.anscombe = dn_anscombe;
	}
	g_hash_table_unref(kv);
	return d;
}

const op_descriptor op_desc_atrous = {
	.id = "wavelets.atrous", .version = 1,
	.image_hook = atrous_image_hook,
	.log_hook = atrous_log_hook,
	.description = N_("Wavelet transform"),
	.mem_ratio = 0.0f,
	.flags = 0,
	.serialize = atrous_serialize, .deserialize = atrous_deserialize,
};

/************* decomposition held for the tool window's lifetime *************/

/* The wavelets tool decomposes once and then reconstructs on every slider
 * move. Each of those reconstructions used to read the whole transform back
 * from the .wave files: 577 MB per channel for a 24 Mpx frame at 6 layers, so
 * 1.7 GB and ~1.3 s per preview refresh -- the ROI fast path only helps when
 * the user has actually made a selection, which is not the common case.
 *
 * So the tool can ask for the transform to be kept in memory for as long as
 * its window is open. The files are still written, so everything outside the
 * window (the wrecons command, a reopened dialog, a later session) behaves
 * exactly as it did. Peak memory during a reconstruction is unchanged -- the
 * file-backed path allocates the same cube to read into -- what is new is that
 * the cube stays allocated between reconstructions. */
static struct {
	wave_transf_des chan[3];
	int nb_chan;
	int rx, ry;
	gboolean held;
} wsession;
static GMutex wsession_mutex;

/* Caller must hold wsession_mutex. */
static void wsession_drop(void) {
	if (!wsession.held)
		return;
	for (int i = 0; i < wsession.nb_chan; i++)
		wave_io_free(&wsession.chan[i]);
	memset(&wsession, 0, sizeof(wsession));
}

void wavelet_session_release(void) {
	g_mutex_lock(&wsession_mutex);
	wsession_drop();
	g_mutex_unlock(&wsession_mutex);
}

/* Global noise sigma of a held channel, so the tool's readout does not have to
 * read the transform back from disk either. FALSE if nothing is held for it. */
gboolean wavelet_session_sigma(int chan, double *sigma_out, int threads) {
	gboolean ok = FALSE;
	g_mutex_lock(&wsession_mutex);
	if (wsession.held && chan >= 0 && chan < wsession.nb_chan)
		ok = !wavelet_sigma_from_data(&wsession.chan[chan], sigma_out, threads);
	g_mutex_unlock(&wsession_mutex);
	return ok;
}

/* Caller must hold wsession_mutex. The geometry check is the same guard the
 * file path gets from the dimensions in the .wave header. */
static gboolean wsession_matches(const fits *fit, int nb_chan) {
	return wsession.held && wsession.nb_chan == nb_chan
			&& wsession.rx == fit->rx && wsession.ry == fit->ry;
}

/* Holding the transform is only worth it if it comfortably fits the memory
 * budget: denoising a held transform works on a copy, so allow for two. */
static gboolean wsession_affordable(const fits *fit, int nbr_plan) {
	const guint64 cube = (guint64) fit->rx * (guint64) fit->ry
			* (guint64) fit->naxes[2] * (guint64) nbr_plan * sizeof(float);
	const guint64 budget = (guint64) get_max_memory_in_MB() * BYTES_IN_A_MB;
	return budget > 0 && cube * 2 <= budget;
}

/************* wavelet transform worker (wavelet command and GUI compute path) *************/

gpointer wavelet_transform_worker(gpointer p) {
	struct wavelet_transform_data *args = (struct wavelet_transform_data *)p;
	const char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" };
	const char *tmpdir = g_get_tmp_dir();
	int retval = 0;
	g_rw_lock_reader_lock(&gfit->rwlock);
	int nb_chan = gfit->naxes[2];
	const int Nl = gfit->ry, Nc = gfit->rx;
	const size_t n = (size_t) Nl * (size_t) Nc;

	/* whatever was held describes the previous decomposition, so it is stale
	 * from here on however this one turns out */
	g_mutex_lock(&wsession_mutex);
	wsession_drop();
	const gboolean keep = args->keep_in_memory
			&& wsession_affordable(gfit, args->Nbr_Plan);
	g_mutex_unlock(&wsession_mutex);

	float *Imag = NULL;
	if (gfit->type == DATA_USHORT || args->anscombe) {
		Imag = malloc(n * sizeof(float));
		if (!Imag) {
			PRINT_ALLOC_ERR;
			retval = 1;
		}
	}
	if (gfit->type != DATA_USHORT && gfit->type != DATA_FLOAT)
		retval = 1;

	for (int i = 0; i < nb_chan && !retval; i++) {
		wave_transf_des wavelet = { 0 };
		float *src;

		if (gfit->type == DATA_USHORT) {
			prepare_rawdata(Imag, Nl, Nc, gfit->pdata[i], com.max_thread);
			if (args->anscombe)
				anscombe_forward(Imag, n, ANSCOMBE_USHORT_SCALE, com.max_thread);
			src = Imag;
		} else if (args->anscombe) {
			memcpy(Imag, gfit->fpdata[i], n * sizeof(float));
			anscombe_forward(Imag, n, ANSCOMBE_FLOAT_SCALE, com.max_thread);
			src = Imag;
		} else {
			src = gfit->fpdata[i];
		}

		if (wavelet_transform_data(src, Nl, Nc, &wavelet, args->Type_Transform,
				args->Nbr_Plan, com.max_thread)) {
			retval = 1;
			break;
		}

		gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		snprintf(wavelet.Name_Imag, MAX_SIZE_NAME_IMAG - 1, "%s", dir);
		wavelet.Name_Imag[MAX_SIZE_NAME_IMAG - 1] = '\0';
		if (wave_io_write(dir, &wavelet))
			retval = 1;
		g_free(dir);

		if (!retval && keep) {
			g_mutex_lock(&wsession_mutex);
			wsession.chan[i] = wavelet; /* ownership moves to the session */
			wsession.nb_chan = i + 1;
			wsession.rx = Nc;
			wsession.ry = Nl;
			wsession.held = TRUE;
			g_mutex_unlock(&wsession_mutex);
		} else {
			wave_io_free(&wavelet);
		}
	}
	free(Imag);

	if (retval)
		wavelet_session_release();
	else
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
	int ret = 0;

	/* Held for the whole reconstruction so the tool window cannot free the
	 * transform underneath a preview that is still running. For an ROI
	 * preview the geometry to match is the full image, not the ROI fit. */
	const gboolean for_roi = gargs->for_roi;
	const int roi_x = gargs->roi_rect.x, roi_y = gargs->roi_rect.y;

	g_mutex_lock(&wsession_mutex);
	const gboolean cached = for_roi
			? (wsession.held && wsession.nb_chan == args->nb_chan
					&& wsession.rx == args->full_rx && wsession.ry == args->full_ry)
			: wsession_matches(fit, args->nb_chan);

	for (int i = 0; i < args->nb_chan && !ret; i++) {
		if (cached) {
			if (for_roi)
				ret = wavelet_reconstruct_data_roi(&wsession.chan[i], args->coef,
						&args->denoise, roi_x, roi_y, fit->rx,
						fit->ry, i, fit, threads);
			else if (fit->type == DATA_USHORT) {
				float *Imag = f_vector_alloc((size_t) fit->rx * (size_t) fit->ry);
				if (!Imag) {
					PRINT_ALLOC_ERR;
					ret = 1;
				} else {
					ret = wavelet_reconstruct_preserving(&wsession.chan[i], Imag,
							args->coef, &args->denoise, threads);
					if (!ret) {
						if (args->denoise.anscombe)
							anscombe_inverse(Imag,
									(size_t) fit->rx * (size_t) fit->ry,
									ANSCOMBE_USHORT_SCALE, threads);
						reget_rawdata(Imag, fit->ry, fit->rx, fit->pdata[i], threads);
					}
					free(Imag);
				}
			} else if (fit->type == DATA_FLOAT) {
				ret = wavelet_reconstruct_preserving(&wsession.chan[i],
						fit->fpdata[i], args->coef, &args->denoise, threads);
				if (!ret && args->denoise.anscombe)
					anscombe_inverse(fit->fpdata[i],
							(size_t) fit->rx * (size_t) fit->ry,
							ANSCOMBE_FLOAT_SCALE, threads);
			} else
				ret = 1;
			continue;
		}

		gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		/* The file reconstruction writers size their output from the .wave
		 * header, NOT from the destination: a transform of a larger image
		 * would write past the end of the destination buffer.  Stale files
		 * are easy to come by — they survive the image (and the FLIS active
		 * layer) that produced them — so refuse any geometry mismatch here.
		 * For an ROI preview the geometry to match is the full image. */
		const int exp_nl = for_roi ? args->full_ry : (int) fit->ry;
		const int exp_nc = for_roi ? args->full_rx : (int) fit->rx;
		int wnl = 0, wnc = 0;
		if (wave_io_read_header(dir, &wnl, &wnc, NULL)
				|| wnl != exp_nl || wnc != exp_nc) {
			siril_log_error(_("The stored wavelet transform (%d x %d) does not match "
					"the current image (%d x %d): run the wavelet transform again.\n"),
					wnc, wnl, exp_nc, exp_nl);
			ret = 1;
		} else if (for_roi) {
			/* Reconstruct only the selection: reads just the ROI rows of the
			 * transform and denoises/reconstructs that small window. */
			ret = wavelet_reconstruct_file_roi(dir, args->coef, &args->denoise,
					roi_x, roi_y, fit->rx, fit->ry, i, fit, threads);
		} else if (fit->type == DATA_USHORT) {
			ret = wavelet_reconstruct_file(dir, args->coef, &args->denoise, fit->pdata[i], threads);
		} else if (fit->type == DATA_FLOAT) {
			ret = wavelet_reconstruct_file_float(dir, args->coef, &args->denoise, fit->fpdata[i], threads);
		} else {
			ret = 1;
		}
		g_free(dir);
	}
	g_mutex_unlock(&wsession_mutex);
	return ret;
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
int atrous_transform_image(fits *fit, const struct atrous_data *args, int threads) {
	threads = wavelet_threads(threads);
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
			prepare_rawdata(Imag, Nl, Nc, fit->pdata[i], threads);
			if (args->anscombe)
				anscombe_forward(Imag, n, ans_scale, threads);
			src = Imag;
			dst = Imag; /* reconstructed into the same scratch, then requantised */
		} else if (args->anscombe) {
			memcpy(Imag, fit->fpdata[i], n * sizeof(float));
			anscombe_forward(Imag, n, ans_scale, threads);
			src = Imag;
			dst = fit->fpdata[i];
		} else {
			src = fit->fpdata[i];
			dst = fit->fpdata[i];
		}

		/* pave_2d_tfo copies its input before recursing, so src is intact here
		 * and dst may safely alias it */
		if (wavelet_transform_data(src, Nl, Nc, &wavelet, args->type,
				args->nbr_plan, threads)) {
			retval = 1;
			break;
		}
		wavelet_denoise_planes(wavelet.Pave.Data, args->type, args->nbr_plan,
				Nl, Nc, &args->denoise, threads);
		retval = wavelet_reconstruct_data(&wavelet, dst, coef, threads);
		wave_io_free(&wavelet);
		if (retval)
			break;

		if (args->anscombe)
			anscombe_inverse(dst, n, ans_scale, threads);
		if (fit->type == DATA_USHORT)
			reget_rawdata(dst, Nl, Nc, fit->pdata[i], threads);
	}
	free(Imag);
	return retval;
}

int atrous_image_hook(struct generic_img_args *gargs, fits *fit, int threads) {
	struct atrous_data *args = (struct atrous_data *)gargs->user;
	return atrous_transform_image(fit, args, threads);
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
	return atrous_transform_image(fit, a_args, threads);
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

int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer,
		int threads) {
	threads = wavelet_threads(threads);
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
						Type, Nbr_Plan, fit->pdata[chan], threads)) {
				retval = 1;
				break;
			}
		}
		else if (fit->type == DATA_FLOAT) {
			/* float wavelet of data [0, 1] */
			Imag = fit->fpdata[chan];
			if (wavelet_transform_float(Imag, fit->ry, fit->rx, &wavelet[chan],
						Type, Nbr_Plan, threads)) {
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
			reget_rawdata(Imag, Nl, Nc, fit->pdata[chan], threads);
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
		get_wavelet_layers(&fit, args->Nbr_Plan, i, args->Type, -1, com.max_thread);
		savefits(filename, &fit);
		g_free(filename);
		g_free(msg);
	}
	clearfits(&fit);
	siril_add_idle(end_wavelets_filter, args);
	return GINT_TO_POINTER(0);
}

/* GUI callbacks moved to src/gui/wavelets.c */
