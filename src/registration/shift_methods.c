
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

#include <complex.h>
#include <fftw3.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "drizzle/cdrizzleutil.h"
#include "algos/quality.h"
#include "algos/demosaicing.h"
#include "io/ser.h"
#include "opencv/opencv.h"
#include "opencv/kombat/kombat.h"


static void normalizeQualityData(struct registration_args *args, double q_min, double q_max) {
	int frame;
	double diff = q_max - q_min;

	/* this case occurs when all images but one are excluded */
	if (diff == 0) {
		q_min = 0;
		diff = (q_max == 0.) ? 1.: q_max;
	}

	for (frame = 0; frame < args->seq->number; ++frame) {
		args->seq->regparam[args->layer][frame].quality -= q_min;
		args->seq->regparam[args->layer][frame].quality /= diff;
		/* if thread has been manually stopped, some values will be < 0 */
		if ((args->seq->regparam[args->layer][frame].quality < 0)
				|| isnan(args->seq->regparam[args->layer][frame].quality))
			args->seq->regparam[args->layer][frame].quality = -1.0;
	}
}

/* Calculate shift in images to be aligned with the reference image, using
 * discrete Fourier transform on a square selected area and matching the
 * phases.
 */
int register_shift_dft(struct registration_args *args) {
	fits fit_ref = { 0 };
	struct timeval t_start, t_end;
	unsigned int frame, size, sqsize, j;
	fftwf_complex *ref, *in, *out, *convol;
	fftwf_plan p, q;
	int ret;
	int abort = 0;
	float nb_frames;
	int cur_nb = 0;
	int ref_image;
	regdata *current_regdata;
	double q_max = 0, q_min = DBL_MAX;

	/* the selection needs to be squared for the DFT */
	assert(args->selection.w == args->selection.h);
	size = args->selection.w;
	sqsize = size * size;

	if (args->filters.filter_included)
		nb_frames = (float) args->seq->selnum;
	else
		nb_frames = (float) args->seq->number;

	if (args->seq->regparam[args->layer]) {
		siril_log_message(_("Recomputing already existing registration for this layer\n"));
		current_regdata = args->seq->regparam[args->layer];
		/* we reset all values as we may register different images */
		memset(current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return -2;
		}
		args->seq->regparam[args->layer] = current_regdata;
	}

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);

	set_progress_bar_data(
			_("Register DFT: loading and processing reference frame"),
			PROGRESS_NONE);
	ret = seq_read_frame_metadata(args->seq, ref_image, &fit_ref);
	if (!ret)
		ret = seq_read_frame_part(args->seq, args->layer, ref_image, &fit_ref,
			&args->selection, FALSE, -1);

	if (ret || ((fit_ref.type == DATA_USHORT && !fit_ref.data) || (fit_ref.type == DATA_FLOAT && !fit_ref.fdata))) {
		siril_log_message(
				_("Register: could not load first image to register, aborting.\n"));
		args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		clearfits(&fit_ref);
		return ret;
	}
	if (args->seq->nb_layers == 1) {
		interpolate_nongreen(&fit_ref);
	}
#if BAYER_DEBUG
	const gchar *green_filename = get_sequence_cache_filename(args->seq, ref_image, "fit", "green_");
	savefits(green_filename, &fit_ref);
#endif
	gettimeofday(&t_start, NULL);

	ref = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	in = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	out = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	convol = fftwf_malloc(sizeof(fftwf_complex) * sqsize);

	set_wisdom_file();
	gchar* wisdomFile = com.pref.fftw_conf.wisdom_file;
#ifdef HAVE_FFTW3F_MULTITHREAD
	fftwf_plan_with_nthreads(com.max_thread);
#endif
	// test for available wisdom
	p = fftwf_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, FFTW_WISDOM_ONLY);
	if (!p) {
		// no wisdom available, load wisdom from file
		fftwf_import_wisdom_from_filename(wisdomFile);
		// test again for wisdom
		p = fftwf_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, FFTW_WISDOM_ONLY);
		if (!p) {
			// build plan with FFTW_MEASURE
			p = fftwf_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, FFTW_MEASURE);
			// save the wisdom
			fftwf_export_wisdom_to_filename(wisdomFile);
		}
	}

	q = fftwf_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
	if (!q) {
		// no wisdom available, load wisdom from file
		fftwf_import_wisdom_from_filename(wisdomFile);
		// test again for wisdom
		q = fftwf_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
		if (!q) {
			// build plan with FFTW_MEASURE
			q = fftwf_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, FFTW_MEASURE | FFTW_DESTROY_INPUT);
			// save the wisdom
			fftwf_export_wisdom_to_filename(wisdomFile);
		}
	}
	fftwf_free(out); // was needed to build the plan, can be freed now
	fftwf_free(convol); // was needed to build the plan, can be freed now

	// copying image selection into the fftw data
	if (fit_ref.type == DATA_USHORT) {
		for (j = 0; j < sqsize; j++)
			ref[j] = (float)fit_ref.data[j];
	}
	else if (fit_ref.type == DATA_FLOAT) {
		for (j = 0; j < sqsize; j++)
			ref[j] = fit_ref.fdata[j];
	}

	// We don't need fit_ref anymore, we can destroy it.
	current_regdata[ref_image].quality = QualityEstimate(&fit_ref, args->layer);
	clearfits(&fit_ref);
	fftwf_execute_dft(p, ref, in); /* repeat as needed */

	fftwf_free(ref); // not needed anymore

	set_shifts(args->seq, ref_image, args->layer, 0.0, 0.0, FALSE);

	q_min = q_max = current_regdata[ref_image].quality;
	int q_index = ref_image;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) \
	if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
#endif
	for (frame = 0; frame < args->seq->number; ++frame) {
		if (abort) continue;
		if (args->run_in_thread && !get_thread_run()) {
			abort = 1;
			continue;
		}
		if (frame == ref_image) continue;
		if (args->filters.filter_included && !args->seq->imgparam[frame].incl)
			continue;

		fits fit = { 0 };
		char tmpmsg[1024], tmpfilename[256];

		seq_get_image_filename(args->seq, frame, tmpfilename);
		g_snprintf(tmpmsg, 1024, _("Register: processing image %s"), tmpfilename);
		set_progress_bar_data(tmpmsg, PROGRESS_NONE);
		int thread_id = -1;
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif
		if (!seq_read_frame_metadata(args->seq, frame, &fit) && !seq_read_frame_part(args->seq, args->layer, frame, &fit,
						&args->selection, FALSE, thread_id)) {
			int x;
			if (args->seq->nb_layers == 1)
				interpolate_nongreen(&fit);
#if BAYER_DEBUG
			const gchar *green_filename = get_sequence_cache_filename(args->seq, frame, "fit", "green_");
			savefits(green_filename, &fit);
#endif

			fftwf_complex *img = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
			fftwf_complex *out2 = fftwf_malloc(sizeof(fftwf_complex) * sqsize);

			// copying image selection into the fftw data
			if (fit.type == DATA_USHORT) {
				for (x = 0; x < sqsize; x++)
					img[x] = (float)fit.data[x];
			}
			if (fit.type == DATA_FLOAT) {
				for (x = 0; x < sqsize; x++)
					img[x] = fit.fdata[x];
			}

			current_regdata[frame].quality = QualityEstimate(&fit, args->layer);
			// after this call, fit data is dead

#ifdef _OPENMP
#pragma omp critical
#endif
			{
				double qual = current_regdata[frame].quality;
				if (qual > q_max) {
					q_max = qual;
					q_index = frame;
				}
				q_min = min(q_min, qual);
			}

			fftwf_execute_dft(p, img, out2); /* repeat as needed */

			fftwf_complex *convol2 = img; // reuse buffer

			for (x = 0; x < sqsize; x++) {
				convol2[x] = in[x] * conjf(out2[x]);
			}

			fftwf_execute_dft(q, convol2, out2); /* repeat as needed */

			int shift = 0;
			for (x = 1; x < sqsize; ++x) {
				if (crealf(out2[x]) > crealf(out2[shift])) {
					shift = x;
					// break or get last value?
				}
			}
			int shifty = shift / size;
			int shiftx = shift % size;
			if (shifty > size / 2) {
				shifty -= size;
			}
			if (shiftx > size / 2) {
				shiftx -= size;
			}

			set_shifts(args->seq, frame, args->layer, (float)shiftx, (float)shifty,
					fit.top_down);

			// We don't need fit anymore, we can destroy it.
			clearfits(&fit);

			/* shiftx and shifty are the x and y values for translation that
			 * would make this image aligned with the reference image.
			 * WARNING: the y value is counted backwards, since the FITS is
			 * stored down from up.
			 */
#ifdef DEBUG
			fprintf(stderr,
					"reg: frame %d, shiftx=%f shifty=%f quality=%g\n",
					args->seq->imgparam[frame].filenum,
					current_regdata[frame].shiftx, current_regdata[frame].shifty,
					current_regdata[frame].quality);
#endif
			g_atomic_int_inc(&cur_nb);
			set_progress_bar_data(NULL, (float)cur_nb / nb_frames);
			fftwf_free(img);
			fftwf_free(out2);
		} else {
			//report_fits_error(ret, error_buffer);
			args->seq->regparam[args->layer] = NULL;
			free(current_regdata);
			abort = ret = 1;
			continue;
		}
	}

	fftwf_destroy_plan(p);
	fftwf_destroy_plan(q);
	fftwf_free(in);
	if (!ret) {
		normalizeQualityData(args, q_min, q_max);

		siril_log_message(_("Registration finished.\n"));
		siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index + 1);
	} else {
		free(args->seq->regparam[args->layer]);
		args->seq->regparam[args->layer] = NULL;
	}
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	return ret;
}

/* simple tool function for KOMBAT registration below; might be
  * adopted by other functions too? */
void _register_kombat_disable_frame(struct registration_args *args, regdata *current_regdata, int frame) {
	args->seq->imgparam[frame].incl = FALSE;
	current_regdata[frame].quality = -1;
	g_atomic_int_dec_and_test(&args->seq->selnum);
}

/* register images using KOMBAT algorithm: calculate shift in images to be
  * aligned with the reference image; images are not modified, only shift
  * parameters are saved in regparam in the sequence.
  *
  * TODO: if multiple selections would be possible, this algorithm can be
  * extended to build mosaics (by searching for corresponding patterns
  * iteratively).
  *
  * TODO: missing optimisations, as many components are computed
  * again and again for each images (as for others registration
  * algorithms, it seems).
  */
int register_kombat(struct registration_args *args) {
	fits fit_templ = { 0 };
	fits fit_ref = { 0 };

	regdata *current_regdata;
	struct timeval t_start, t_end;
	float nb_frames;
	int cur_nb = 0;
	double q_max = 0, q_min = DBL_MAX;
	double qual;
	int frame;
	int ref_idx;
	int ret_ref, ret_templ;
	int abort = 0;
	rectangle full = { 0 };

	reg_kombat ref_align;

	if (args->seq->regparam[args->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		current_regdata = args->seq->regparam[args->layer];
		/* we reset all values as we may register different images */
		memset(current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		current_regdata = (regdata*) calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		args->seq->regparam[args->layer] = current_regdata;
	}


	if (args->filters.filter_included)
		nb_frames = (float) args->seq->selnum;
	else
		nb_frames = (float) args->seq->number;

	/* loading reference frame */
	ref_idx = sequence_find_refimage(args->seq);
	int q_index = ref_idx;
	full.w = args->seq->rx;
	full.h = args->seq->ry;

	set_progress_bar_data(
			_("Register using KOMBAT: loading and processing reference frame"),
			PROGRESS_NONE);

	/* we want pattern (the selection) to be located on each image */
	ret_templ = seq_read_frame_part(args->seq, args->layer, ref_idx, &fit_templ,
			&args->selection, FALSE, -1);
	ret_ref = seq_read_frame_part(args->seq, args->layer, ref_idx, &fit_ref,
			&full, FALSE, -1);

	if (ret_templ || ret_ref) {
		siril_log_message(
				_("Register: could not load first image to register, aborting.\n"));
		args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		clearfits(&fit_ref);
		clearfits(&fit_templ);
		return 1;
	}
	if (args->seq->nb_layers == 1) {
		interpolate_nongreen(&fit_templ);
		interpolate_nongreen(&fit_ref);
#if BAYER_DEBUG
		const gchar *green_filename = get_sequence_cache_filename(args->seq, ref_idx, "fit", "green_");
		savefits(green_filename, &fit_ref);
#endif
	}

	gettimeofday(&t_start, NULL);
	set_shifts(args->seq, ref_idx, args->layer, 0.0, 0.0, FALSE);

	/* we want pattern position on the reference image */
	kombat_find_template(ref_idx, args, &fit_templ, &fit_ref, &ref_align, NULL, NULL);
	/* main part starts here */
	int max_threads;
#ifdef _OPENMP
	max_threads = omp_get_max_threads() + 1;
#else
		max_threads = 1;
#endif

	void *caches[max_threads];
	for (int i = 0; i < max_threads; i++)
		caches[i] = NULL;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ || args->seq->type == SEQ_INTERNAL) && fits_is_reentrant()))
#endif
		for (frame = 0; frame < args->seq->number; frame++) {
			if (abort)
				continue;
			if (args->run_in_thread && !get_thread_run()) {
				abort = 1;
				continue;
			}

			if (args->filters.filter_included && !args->seq->imgparam[frame].incl)
				continue;

			fits cur_fit = { 0 };
			char tmpmsg[1024], tmpfilename[256];

			seq_get_image_filename(args->seq, frame, tmpfilename);
			g_snprintf(tmpmsg, 1024, _("Register: processing image %s"),
					tmpfilename);
			set_progress_bar_data(tmpmsg, PROGRESS_NONE);
			int thread_id = -1;
#ifdef _OPENMP
			thread_id = omp_get_thread_num();
#endif

			if (seq_read_frame_part(args->seq, args->layer, frame, &cur_fit, &full, FALSE, thread_id)) {
				siril_log_message(_("Cannot perform KOMBAT alignment for frame %d\n"), frame + 1);
				/* we exclude this frame */
				_register_kombat_disable_frame(args, current_regdata, frame);
				continue;
			} else {
				if (args->seq->nb_layers == 1) {
					interpolate_nongreen(&cur_fit);
				}
				qual = QualityEstimate(&cur_fit, args->layer);
				current_regdata[frame].quality = qual;
				if (frame == ref_idx) {
					clearfits(&cur_fit);
					continue;
				}
			}

			/* we want pattern position on any single image */
			reg_kombat cur_align;
			if (kombat_find_template(frame, args, &fit_templ, &cur_fit,
					&cur_align, &ref_align, &(caches[thread_id + 1]))) {
				siril_log_color_message(_("Register: KOMBAT could not find alignment pattern on image #%d.\n"),	"red", frame);
				/* we exclude this frame too */
				_register_kombat_disable_frame(args, current_regdata, frame);
			} else {
				set_shifts(args->seq, frame, args->layer,
						ref_align.dx - cur_align.dx,
						(ref_align.dy - cur_align.dy), cur_fit.top_down);
#if BAYER_DEBUG
				const gchar *green_filename = get_sequence_cache_filename(args->seq, frame, "fit", "green_");
				savefits(green_filename, &cur_fit);
#endif
			}

			g_atomic_int_inc(&cur_nb);
			set_progress_bar_data(NULL, (float)cur_nb / nb_frames);

			// We don't need fit anymore, we can destroy it.
			clearfits(&cur_fit);
		}

	siril_log_message(_("Registration finished.\n"));

	for (frame = 0; frame < args->seq->number; frame++) {
		if (args->seq->imgparam[frame].incl) {
			qual = current_regdata[frame].quality;
			q_min = min(q_min, qual);
			if (qual>q_max) {
				q_max = qual;
				q_index = frame;
			}
		}
	}
	normalizeQualityData(args, q_min, q_max);

	clearfits(&fit_templ);
	clearfits(&fit_ref);

	siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index + 1);

	for (int i = 0; i < max_threads; i++)
		kombat_done(&caches[i]);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	return 0;
}

/* register images: calculate shift in images to be aligned with the reference image;
 * images are not modified, only shift parameters are saved in regparam in the sequence.
 * layer is the layer on which the registration will be done, green by default (set in siril_init())
 */
int register_shift_fwhm(struct registration_args *args) {
	int frame, ref_image;
	float nb_frames, cur_nb = 0.f;
	double reference_xpos, reference_ypos;
	double fwhm_min = DBL_MAX;
	int fwhm_index = -1;
	regdata *current_regdata;

	framing_mode framing = ORIGINAL_FRAME;
	if (args->follow_star)
		framing = FOLLOW_STAR_FRAME;

	/* First and longest step: get the minimization data on one star for all
	 * images to register, which provides FWHM but also star coordinates */
	// TODO: detect that it was already computed, and don't do it again
	// -> should be done at a higher level and passed in the args
	if (seqpsf(args->seq, args->layer, TRUE, FALSE, !args->filters.filter_included, framing, FALSE, TRUE))
		return 1;

	// regparam is managed in seqpsf idle function already
	current_regdata = args->seq->regparam[args->layer];

	if (args->filters.filter_included)
		nb_frames = (float)args->seq->selnum;
	else nb_frames = (float)args->seq->number;

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);
	if (!current_regdata[ref_image].fwhm_data) {
		siril_log_message(
				_("Registration PSF: failed to compute PSF for reference frame at least\n"));
		return -1;
	}
	reference_xpos = current_regdata[ref_image].fwhm_data->xpos;
	reference_ypos = current_regdata[ref_image].fwhm_data->ypos;

	fwhm_min = current_regdata[ref_image].fwhm_data->fwhmx;

	fwhm_index = ref_image;

	/* Second step: align image by aligning star coordinates together */
	for (frame = 0; frame < args->seq->number; frame++) {
		if (args->run_in_thread && !get_thread_run())
			break;
		if (args->filters.filter_included && !args->seq->imgparam[frame].incl) {
			// current_regdata was set with identity matrices
			// need to return null matcices for unselected frames
			SetNullH(&current_regdata[frame].H);
			continue;
		}
		if (frame == ref_image || !current_regdata[frame].fwhm_data) {
			set_shifts(args->seq, frame, args->layer, 0.0, 0.0, FALSE);
			continue;
		}
		if (current_regdata[frame].fwhm < fwhm_min
				&& current_regdata[frame].fwhm > 0.0) {
			fwhm_min = current_regdata[frame].fwhm;
			fwhm_index = frame;
		}
		double shiftx = reference_xpos - current_regdata[frame].fwhm_data->xpos;
		double shifty = current_regdata[frame].fwhm_data->ypos - reference_ypos;
		current_regdata[frame].H = H_from_translation(shiftx, shifty);

		fprintf(stderr, "reg: file %d, shiftx=%f shifty=%f\n",
				args->seq->imgparam[frame].filenum, shiftx, shifty);
		cur_nb += 1.f;
		set_progress_bar_data(NULL, cur_nb / nb_frames);
	}

	siril_log_message(_("Registration finished.\n"));
	siril_log_color_message(_("Best frame: #%d with fwhm=%.3g.\n"), "bold",
			fwhm_index + 1, fwhm_min);
	return 0;
}
