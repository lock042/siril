
/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "gui/sequence_list.h"
#include "core/initfile.h"
#include "core/OS_utils.h"
#include "core/siril_app_dirs.h"
#include "registration/registration.h"
#include "registration/matching/misc.h"
#include "registration/matching/match.h"
#include "registration/matching/atpmatch.h"
#include "stacking/stacking.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"
#include "algos/quality.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "opencv/kombat/kombat.h"

static gboolean keep_noout_state = FALSE;

#undef DEBUG

static char *tooltip_text[] = {
	N_("<b>One Star Registration</b>: This is the simplest method to register deep-sky images. "
		"Because only one star is concerned for register, images are aligned using shifting "
		"(at a fraction of pixel). No rotation or scaling are performed. "
		"Shifts at pixel precision are saved in seq file."),
	N_("<b>Two or Three Stars Registration</b>: This method looks like the one star registration "
		"except one need to select two or three stars. This is very useful for field with a "
		"few stars."),
	N_("<b>Global Star Alignment</b>: This is a more powerful and accurate algorithm (but also "
		"slower) to perform deep-sky images. The global matching is based on triangle "
		"similarity method for automatically identify common stars in each image. A new "
		"sequence is created with the prefix of your choice (r_ by default)."),
	N_("<b>Two-Pass Global Star Alignment</b>: The global star alignment is done in two passes, "
		"allowing the reference frame to be chosen from detected star information instead of "
		"automatically choosing the first frame of the sequence."),
	N_("<b>Image Pattern Alignment</b>: This is a simple registration by translation method "
		"using cross correlation in the spatial domain. This method is fast and is used to "
		"register planetary movies. It can also be used for some deep-sky images registration. "
		"Shifts at pixel precision are saved in seq file."),
	N_("<b>KOMBAT</b>: This simple algorithm tries to locate a single pattern on images and to "
		"align them accordingly. Only translation is taken into account yet."),
	N_("<b>Comet/Asteroid Registration</b>: This algorithm is dedicated to the comet and asteroid "
		"registration. It is necessary to have timestamps stored in FITS header and to load a "
		"sequence of star aligned images. This methods makes a translation of a certain number "
		"of pixels depending on the timestamp of each images and the global shift of the "
		"object between the first and the last image."),
	N_("<b>Apply existing registration</b>: This is not an algorithm but rather a commodity to "
		"apply previously computed registration data stored in the sequence file. The "
		"interpolation method and simplified drizzle can be selected in the Output "
		"Registration section and it can be applied on selected images only, to avoid saving "
		"unnecessary images.")
};

static char *reg_frame_registration[] = {
	"framing-default.svg",
	"framing-max.svg",
	"framing-min.svg",
	"framing-cog.svg"
};

/*Possible values for max stars combo box
Needs to be consistent with list in comboreg_maxstars*/
static int maxstars_values[] = { 100, 200, 500, 1000, 2000 };

int register_kombat(struct registration_args *args);

/* callback for the selected area event */
void _reg_selected_area_callback() {
	if (!com.headless)
		update_reg_interface(TRUE);
}

static struct registration_method *reg_methods[NUMBER_OF_METHODS + 1];
static gboolean end_register_idle(gpointer p);

struct registration_method *new_reg_method(const char *name, registration_function f,
		selection_type s, registration_type t) {
	struct registration_method *reg = malloc(sizeof(struct registration_method));
	reg->name = strdup(name);
	reg->method_ptr = f;
	reg->sel = s;
	reg->type = t;
	return reg;
}

void initialize_registration_methods() {
	GtkComboBoxText *regcombo;
	int i = 0, j = 0;
	GString *tip;
	gchar *ctip;

	reg_methods[i++] = new_reg_method(_("One Star Registration (deep-sky)"),
			&register_shift_fwhm, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Two or Three Stars Registration (deep-sky)"),
			&register_3stars, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Global Star Alignment (deep-sky)"),
			&register_star_alignment, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Two-Pass Global Star Alignment (deep-sky)"),
			&register_multi_step_global, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Image Pattern Alignment (planetary - full disk)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
	reg_methods[i++] = new_reg_method(_("KOMBAT (planetary surfaces or full disk)"),
			&register_kombat, REQUIRES_ANY_SELECTION, REGTYPE_PLANETARY);
	reg_methods[i++] = new_reg_method(_("Comet/Asteroid Registration"),
			&register_comet, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Apply Existing Registration"),
			&register_apply_reg, REQUIRES_NO_SELECTION, REGTYPE_APPLY);
	reg_methods[i] = NULL;

	tip = g_string_new ("");
	for (j = 0; j < i; j ++) {
		tip = g_string_append(tip, _(tooltip_text[j]));
		if (j < i - 1)
			tip = g_string_append(tip, "\n\n");
	}
	ctip = g_string_free (tip, FALSE);
	gtk_widget_set_tooltip_markup(lookup_widget("comboboxregmethod"), ctip);
	g_free(ctip);

	/* fill comboboxregmethod */
	regcombo = GTK_COMBO_BOX_TEXT(lookup_widget("comboboxregmethod"));
	gtk_combo_box_text_remove_all(regcombo);
	i = 0;
	while (reg_methods[i] != NULL) {
		gtk_combo_box_text_append_text(regcombo, reg_methods[i]->name);
		siril_log_message(_("Loading registration method: %s\n"),
				reg_methods[i]->name);
		i++;
	}
	if (i > 0) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(regcombo), com.pref.gui.reg_settings);
	}

	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("ComboBoxRegInter")), com.pref.gui.reg_interpolation);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_reg_clamp")), com.pref.gui.reg_clamping);
	gtk_widget_set_sensitive(lookup_widget("toggle_reg_clamp"),
			com.pref.gui.reg_interpolation == OPENCV_LANCZOS4 || com.pref.gui.reg_interpolation == OPENCV_CUBIC);

	/* register to the new area selected event */
	register_selection_update_callback(_reg_selected_area_callback);
}

struct registration_method *get_selected_registration_method() {
	GtkComboBoxText *regcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(gui.builder, "comboboxregmethod"));
	int index = 0;

	gchar *text = gtk_combo_box_text_get_active_text (regcombo);
	while (reg_methods[index] && text != NULL) {
		if (!strcmp(reg_methods[index]->name, text))
			break;
		index++;
	}
	g_free(text);
	return reg_methods[index];
}

static void normalizeQualityData(struct registration_args *args, double q_min, double q_max) {
	int frame;
	double diff = q_max - q_min;

	/* this case occurs when all images but one are excluded */
	if (diff == 0) {
		q_min = 0;
		diff = q_max;
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
	float nb_frames, cur_nb;
	int ref_image;
	regdata *current_regdata;
	double q_max = 0, q_min = DBL_MAX;
	int q_index = -1;

	/* the selection needs to be squared for the DFT */
	assert(args->selection.w == args->selection.h);
	size = args->selection.w;
	sqsize = size * size;

	if (args->filters.filter_included)
		nb_frames = (float) args->seq->selnum;
	else
		nb_frames = (float) args->seq->number;

	if (args->seq->regparam[args->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
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
	ret = seq_read_frame_part(args->seq, args->layer, ref_image, &fit_ref,
			&args->selection, FALSE, -1);

	if (ret) {
		siril_log_message(
				_("Register: could not load first image to register, aborting.\n"));
		args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		clearfits(&fit_ref);
		return ret;
	}

	gettimeofday(&t_start, NULL);

	ref = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	in = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	out = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	convol = fftwf_malloc(sizeof(fftwf_complex) * sqsize);

	gchar* wisdomFile = g_build_filename(g_get_user_cache_dir(), "siril_fftw.wisdom", NULL);

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
	g_free(wisdomFile);
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
	q_index = ref_image;

	cur_nb = 0.f;

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
		if (!seq_read_frame_part(args->seq, args->layer, frame, &fit,
						&args->selection, FALSE, thread_id)) {
			int x;
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
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / nb_frames);
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
		if (args->x2upscale)
			args->seq->upscale_at_stacking = 2.0;
		else
			args->seq->upscale_at_stacking = 1.0;
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
	#ifdef _OPENMP
	#pragma omp atomic
	#endif
	args->seq->selnum--;
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
	float nb_frames, cur_nb = 0;
	double q_max = 0, q_min = DBL_MAX;
	double qual;
	int frame;
	int ref_idx;
	int ret;
	int abort = 0;
	rectangle full = { 0 };
	int q_index = -1;

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
	q_index = ref_idx;

	set_progress_bar_data(
			_("Register using KOMBAT: loading and processing reference frame"),
			PROGRESS_NONE);

	/* we want pattern (the selection) to be located on each image */
	ret = seq_read_frame_part(args->seq, args->layer, ref_idx, &fit_templ,
			&args->selection, FALSE, -1);

	/* we load reference image just to get dimensions of images,
	 in order to call seq_read_frame_part() and use only the desired layer, over each full image */
	seq_read_frame(args->seq, ref_idx, &fit_ref, FALSE, -1);
	full.x = full.y = 0;
	full.w = fit_ref.rx;
	full.h = fit_ref.ry;

	if (ret || seq_read_frame_part(args->seq, args->layer, ref_idx, &fit_ref, &full, FALSE, -1)) {
		siril_log_message(
				_("Register: could not load first image to register, aborting.\n"));
		args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		clearfits(&fit_templ);
		return 1;
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
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
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
			}

#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / nb_frames);

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
	if (seqpsf(args->seq, args->layer, TRUE, !args->filters.filter_included, framing, FALSE, FALSE))
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

	if (args->x2upscale)
		args->seq->upscale_at_stacking = 2.0;
	else
		args->seq->upscale_at_stacking = 1.0;

	siril_log_message(_("Registration finished.\n"));
	siril_log_color_message(_("Best frame: #%d with fwhm=%.3g.\n"), "bold",
			fwhm_index + 1, fwhm_min);
	return 0;
}

void on_comboboxregmethod_changed(GtkComboBox *box, gpointer user_data) {
	int index = 0;
	gchar *text = gtk_combo_box_text_get_active_text (GTK_COMBO_BOX_TEXT(box));

	while (reg_methods[index] && text != NULL) {
		if (!strcmp(reg_methods[index]->name, text))
			break;
		index++;
	}
	g_free(text);

	com.pref.gui.reg_settings = index;
	reset_3stars();
	update_reg_interface(TRUE);
}

void on_toggle_reg_clamp_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean active = gtk_toggle_button_get_active(button);
	gtk_widget_set_sensitive(lookup_widget("spin_reg_clamp"), active);

	com.pref.gui.reg_clamping = active;
}

void on_ComboBoxRegInter_changed(GtkComboBox *box, gpointer user_data) {
	com.pref.gui.reg_interpolation = gtk_combo_box_get_active(box);
	gtk_widget_set_sensitive(lookup_widget("toggle_reg_clamp"),
			com.pref.gui.reg_interpolation == OPENCV_LANCZOS4
					|| com.pref.gui.reg_interpolation == OPENCV_CUBIC);
}

void on_comboreg_transfo_changed(GtkComboBox *box, gpointer user_data) {
	GtkAdjustment *register_minpairs = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "register_minpairs"));
	double val = gtk_adjustment_get_value(register_minpairs);

	switch(gtk_combo_box_get_active(box)) {
	case SHIFT_TRANSFORMATION:
	case SIMILARITY_TRANSFORMATION:
	case AFFINE_TRANSFORMATION:
		gtk_adjustment_set_lower (register_minpairs, 3);
		break;
	case HOMOGRAPHY_TRANSFORMATION:
		gtk_adjustment_set_lower (register_minpairs, 4);
		if (val < 4)
			gtk_adjustment_set_value(register_minpairs, 4);
		break;
	default:
		printf("on_comboreg_transfo_changed: Value not handled.\n");
	}
}

void on_comboreg_sel_all_combobox_changed(GtkComboBox *box, gpointer user_data) {
	update_reg_interface(TRUE);
}

int get_registration_layer(sequence *seq) {
	if (!com.script && seq == &com.seq) {
		GtkComboBox *registbox = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
		int reglayer = gtk_combo_box_get_active(registbox);
		if (!seq || !seq->regparam || seq->nb_layers < 0 || seq->nb_layers <= reglayer)
			return -1;
		return reglayer;
	} else {
		// find first available regdata
		if (!seq || !seq->regparam || seq->nb_layers < 0)
			return -1;
		int i;
		for (i = 0; i < seq->nb_layers; i++)
			if (seq->regparam[i])
				return i;
		return -1;
	}
}

int get_first_selected(sequence *seq) {
	if (!seq || !seq->imgparam) return -1;
	fix_selnum(seq, TRUE);
	for (int i = 0; i < seq->number; i++)
		if (seq->imgparam[i].incl)
			return i;
	return -1;
}

gboolean layer_has_registration(sequence *seq, int layer) {
	if (!seq || layer < 0 || !seq->regparam || seq->nb_layers < 0 || layer >= seq->nb_layers || !seq->regparam[layer] ) return FALSE;
	return TRUE;
}
gboolean layer_has_usable_registration(sequence *seq, int layer) {
	int min, max;
	guess_transform_from_seq(seq, layer, &min, &max, FALSE); // will check first that layer_has_registration
	if (max <= -1) return FALSE; // max <= -1 means all H matrices are identity or null
	return TRUE;
}

int seq_has_any_regdata(sequence *seq) {
	if (!seq || !seq->regparam || seq->nb_layers < 0)
		return -1;
	int i;
	for (i = 0; i < seq->nb_layers; i++)
		if (seq->regparam[i])
			return i;
	return -1;
}

/****************************************************************/

#define MAX_FILTERS 5
static struct filtering_tuple regfilters[MAX_FILTERS] = { 0 };

static void update_filters_registration();

void on_regsel_changed(GtkComboBox *widget, gpointer user_data) {
	update_filters_registration();
}

void on_regspinbut_percent_change(GtkSpinButton *spinbutton, gpointer user_data) {
	update_filters_registration();
}

void on_filter_add4_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter5"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin5"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_add5"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem5"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter5"), TRUE);
	update_filters_registration();
}

void on_filter_add5_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter6"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin6"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem6"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter6"), TRUE);
	update_filters_registration();
}

void on_filter_rem5_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter5"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin5"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_add5"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem5"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter5"), FALSE);
	update_filters_registration();
}

void on_filter_rem6_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter6"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin6"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem6"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter6"), FALSE);
	update_filters_registration();
}

static void get_reg_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter) {
	int filter, guifilter, channel = 0, type;
	double percent = 0.0;
	static GtkComboBox *filter_combo[3] = { NULL };
	static GtkAdjustment *stackadj[3] = { NULL };
	static GtkWidget *spin[3] = { NULL };
	if (!spin[0]) {
		spin[0] = lookup_widget("stackspin4");
		spin[1] = lookup_widget("stackspin5");
		spin[2] = lookup_widget("stackspin6");
		stackadj[0] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[0]));
		stackadj[1] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[1]));
		stackadj[2] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[2]));
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter4"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter5"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter6"));
	}
	for (filter = 0, guifilter = 0; guifilter < 3; guifilter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[guifilter]))) {
			continue;
		}

		type = gtk_combo_box_get_active(filter_combo[guifilter]);
		if (type != ALL_IMAGES && type != SELECTED_IMAGES) {
			channel = get_registration_layer(&com.seq);
			percent = gtk_adjustment_get_value(stackadj[guifilter]);
		}

		switch (type) {
			default:
			case ALL_IMAGES:
				regfilters[filter].filter = seq_filter_all;
				regfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				break;
			case SELECTED_IMAGES:
				regfilters[filter].filter = seq_filter_included;
				regfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				break;
			case BEST_PSF_IMAGES:
				regfilters[filter].filter = seq_filter_fwhm;
				regfilters[filter].param = compute_highest_accepted_fwhm(
						&com.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_WPSF_IMAGES:
				regfilters[filter].filter = seq_filter_weighted_fwhm;
				regfilters[filter].param = compute_highest_accepted_weighted_fwhm(
						&com.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_ROUND_IMAGES:
				regfilters[filter].filter = seq_filter_roundness;
				regfilters[filter].param = compute_lowest_accepted_roundness(
						&com.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_BKG_IMAGES:
				regfilters[filter].filter = seq_filter_background;
				regfilters[filter].param = compute_highest_accepted_background(
						&com.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_NBSTARS_IMAGES:
				regfilters[filter].filter = seq_filter_nbstars;
				regfilters[filter].param = compute_lowest_accepted_nbstars(
						&com.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_QUALITY_IMAGES:
				regfilters[filter].filter = seq_filter_quality;
				regfilters[filter].param = compute_lowest_accepted_quality(
						&com.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
		}
		filter++;
	}
	regfilters[filter].filter = NULL;

	if (filter == 1) {
		*filtering_criterion = regfilters[0].filter;
		*filtering_parameter = regfilters[0].param;
	} else {
		*filtering_criterion = create_multiple_filter_from_list(regfilters);
		*filtering_parameter = 0.0;
	}
}

static void update_filter_label(seq_image_filter filtering_criterion, double filtering_parameter) {
	static GtkComboBox *filter_combo[3] = { NULL };
	static GtkLabel *filter_label[3] = { NULL };
	static GtkLabel *result_label = NULL;
	if (!filter_combo[0]) {
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter4"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter5"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter6"));
		filter_label[0] = GTK_LABEL(lookup_widget("labelfilter4"));
		filter_label[1] = GTK_LABEL(lookup_widget("labelfilter5"));
		filter_label[2] = GTK_LABEL(lookup_widget("labelfilter6"));
		result_label = GTK_LABEL(lookup_widget("regfilter_label"));
	}

	for (int filter = 0; filter < 3; filter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[filter]))) {
			break;
		}

		int type = gtk_combo_box_get_active(filter_combo[filter]);
		double param = regfilters[filter].param;
		gchar *filter_str;
		if (param == DBL_MIN || param == DBL_MAX || param == 0.0) {
			if (type == ALL_IMAGES || type == SELECTED_IMAGES)
				filter_str = g_strdup("");
			else filter_str = g_strdup("N/A");
		} else {
			switch (type) {
			default:
			case ALL_IMAGES:
			case SELECTED_IMAGES:
				filter_str = g_strdup("");
				break;
			case BEST_PSF_IMAGES:
			case BEST_WPSF_IMAGES:
				filter_str = g_strdup_printf("< %.2lf", param);
				break;
			case BEST_BKG_IMAGES :
				filter_str = (param < 1.) ? g_strdup_printf("< %.5lf", param) : g_strdup_printf("< %d", (int)param);
				break;
			case BEST_ROUND_IMAGES:
			case BEST_QUALITY_IMAGES:
				filter_str = g_strdup_printf("> %.3lf", param);
				break;
			case BEST_NBSTARS_IMAGES:
				filter_str = g_strdup_printf("> %d", (int)param);
				break;
			}
		}
		gtk_label_set_text(filter_label[filter], filter_str);
		g_free(filter_str);
	}

	int nb_images_to_stack = compute_nb_filtered_images(&com.seq,
			filtering_criterion, filtering_parameter);
	gchar *labelbuffer = g_strdup_printf(_("Registering %d images of the %d of the sequence"),
			nb_images_to_stack, com.seq.number);
	gtk_label_set_text(result_label, labelbuffer);
	g_free(labelbuffer);
}

static void update_filters_registration() {
	if (!sequence_is_loaded())
		return;
	siril_debug_print("updating registration filters GUI\n");
	seq_image_filter criterion;
	double param;
	get_reg_sequence_filtering_from_gui(&criterion, &param);
	update_filter_label(criterion, param);
}

/* Selects the "register all" or "register selected" according to the number of
 * selected images, if argument is false.
 * Verifies that enough images are selected and an area is selected.
 */
void update_reg_interface(gboolean dont_change_reg_radio) {
	static GtkWidget *go_register = NULL, *follow = NULL, *cumul_data = NULL, *noout = NULL;
	static GtkLabel *labelreginfo = NULL;
	static GtkComboBox *reg_all_sel_box = NULL, *reglayer = NULL, *filter_combo = NULL;
	static GtkNotebook *notebook_reg = NULL;
	int nb_images_reg; /* the number of images to register */
	struct registration_method *method;
	gboolean selection_is_done;
	gboolean has_reg, ready;

	if (!go_register) {
		go_register = lookup_widget("goregister_button");
		follow = lookup_widget("followStarCheckButton");
		reg_all_sel_box = GTK_COMBO_BOX(lookup_widget("reg_sel_all_combobox"));
		labelreginfo = GTK_LABEL(lookup_widget("labelregisterinfo"));
		notebook_reg = GTK_NOTEBOOK(lookup_widget("notebook_registration"));
		cumul_data = lookup_widget("check_button_comet");
		noout = lookup_widget("regNoOutput");
		reglayer = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
		filter_combo = GTK_COMBO_BOX(lookup_widget("combofilter4"));
	}

	if (!dont_change_reg_radio) {
		if (com.seq.selnum < com.seq.number)
			gtk_combo_box_set_active(reg_all_sel_box, 1);
		else
			gtk_combo_box_set_active(reg_all_sel_box, 0);
	}

	selection_is_done = (com.selection.w > 0 && com.selection.h > 0);

	/* initialize default */
	gtk_notebook_set_current_page(notebook_reg, REG_PAGE_MISC);
	gtk_widget_set_visible(cumul_data, FALSE);

	/* getting the selected registration method */
	method = get_selected_registration_method();

	/* show the appropriate frame selection widgets */
	gboolean isapplyreg = method->method_ptr == &register_apply_reg;
	gtk_widget_set_visible(GTK_WIDGET(reg_all_sel_box), !isapplyreg);
	gtk_widget_set_visible(lookup_widget("seq_filters_box_reg"), isapplyreg);
	if (isapplyreg) {
		if (!dont_change_reg_radio && com.seq.selnum < com.seq.number) {
			gtk_combo_box_set_active(filter_combo, SELECTED_IMAGES);
		}
		update_filters_registration();
	}

	/* number of registered image */
	nb_images_reg = gtk_combo_box_get_active(reg_all_sel_box) == 0 ? com.seq.number : com.seq.selnum;

	/* registration data exists for the selected layer */
	has_reg = layer_has_registration(&com.seq, gtk_combo_box_get_active(reglayer));

	if (method && nb_images_reg > 1 && (selection_is_done || method->sel == REQUIRES_NO_SELECTION) && (has_reg || method->type != REGTYPE_APPLY) ) {
		if (method->method_ptr == &register_star_alignment || method->method_ptr == &register_multi_step_global) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_GLOBAL);
		} else if (method->method_ptr == &register_comet) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_COMET);
		} else if (method->method_ptr == &register_3stars) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_3_STARS);
		} else if (method->method_ptr == &register_kombat) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_KOMBAT);
		} else if (method->method_ptr == &register_apply_reg) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_APPLYREG);
		}
		gtk_widget_set_visible(follow, (method->method_ptr == &register_shift_fwhm) || (method->method_ptr == &register_3stars));
		gtk_widget_set_visible(cumul_data, method->method_ptr == &register_comet);
		ready = TRUE;
		if (method->method_ptr == &register_3stars || method->method_ptr == &register_shift_fwhm) {
			ready = _3stars_check_selection(); // checks that the right image is loaded based on doall and dofollow
		}
		else if (gfit.naxes[2] == 1 && gfit.bayer_pattern[0] != '\0') {
			gtk_label_set_text(labelreginfo, _("Debayer the sequence for registration"));
			ready = FALSE;
		}
		else gtk_label_set_text(labelreginfo, "");
		// the 3 stars method has special GUI requirements
		if (method->method_ptr == &register_3stars) {
			if (!ready) gtk_widget_set_sensitive(go_register,FALSE);
			else _3stars_check_registration_ready();
		} else gtk_widget_set_sensitive(go_register, ready);
	} else {
		gtk_widget_set_sensitive(go_register, FALSE);
		if (nb_images_reg <= 1 && !selection_is_done) {
			if (sequence_is_loaded()) {
				if (method && method->sel == REQUIRES_NO_SELECTION) {
					gtk_label_set_text(labelreginfo, _("Select images in the sequence"));
				} else {
					gtk_label_set_text(labelreginfo, _("Select an area in image first, and select images in the sequence"));
				}
			}
			else {
				gtk_label_set_text(labelreginfo, _("Load a sequence first"));
			}
		} else if (nb_images_reg <= 1) {
			gtk_label_set_text(labelreginfo, _("Select images in the sequence"));
		} else if (method->type == REGTYPE_APPLY && !has_reg){
			gtk_label_set_text(labelreginfo, _("Select a layer with existing registration"));
		} else {
			gtk_label_set_text(labelreginfo, _("Select an area in image first"));
		}
	}
	/* we temporary save value as keep_noout_state will be changed in the callback */
	gboolean save_state = keep_noout_state;
	// for now, methods which do not save images but only shift in seq files are constrained to this option (no_output is true and unsensitive)

	if (method && ((method->method_ptr == &register_comet) ||
			(method->method_ptr == &register_kombat) ||
			(method->method_ptr == &register_shift_fwhm) ||
			(method->method_ptr == &register_shift_dft) ||
			(method->method_ptr == &register_multi_step_global))) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(noout), TRUE);
		gtk_widget_set_sensitive(noout, FALSE);
	} else if (method && method->method_ptr == &register_apply_reg) { // cannot have no output with apply registration method
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(noout), FALSE);
		gtk_widget_set_sensitive(noout, FALSE);
	} else {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(noout), save_state);
		gtk_widget_set_sensitive(noout, TRUE);
	}
	keep_noout_state  = save_state;
}

/* try to maximize the area within the image size (based on gfit)
 * hsteps and vsteps are used to resize the selection zone when it is larger than the image
 * they must be at least 2 */
void compute_fitting_selection(rectangle *area, int hsteps, int vsteps, int preserve_square) {
	//fprintf(stdout, "function entry: %d,%d,\t%dx%d\n", area->x, area->y, area->w, area->h);
	if (area->x >= 0 && area->x + area->w <= gfit.rx && area->y >= 0
			&& area->y + area->h <= gfit.ry)
		return;

	if (area->x < 0) {
		area->x++;
		if (area->x + area->w > gfit.rx) {
			/* reduce area */
			area->w -= hsteps;
			if (preserve_square) {
				area->h -= vsteps;
				area->y++;
			}
		}
	} else if (area->x + area->w > gfit.rx) {
		area->x--;
		if (area->x < 0) {
			/* reduce area */
			area->x++;
			area->w -= hsteps;
			if (preserve_square) {
				area->h -= vsteps;
				area->y++;
			}
		}
	}

	if (area->y < 0) {
		area->y++;
		if (area->y + area->h > gfit.ry) {
			/* reduce area */
			area->h -= hsteps;
			if (preserve_square) {
				area->w -= vsteps;
				area->x++;
			}
		}
	} else if (area->y + area->h > gfit.ry) {
		area->y--;
		if (area->y < 0) {
			/* reduce area */
			area->y++;
			area->h -= vsteps;
			if (preserve_square) {
				area->w -= hsteps;
				area->x++;
			}
		}
	}

	return compute_fitting_selection(area, hsteps, vsteps, preserve_square);
}

void get_the_registration_area(struct registration_args *reg_args,
		struct registration_method *method) {
	int max;
	switch (method->sel) {
		/* even in the case of REQUIRES_NO_SELECTION selection is needed for MatchSelection of starAlignment */
		case REQUIRES_NO_SELECTION:
		case REQUIRES_ANY_SELECTION:
			memcpy(&reg_args->selection, &com.selection, sizeof(rectangle));
			break;
		case REQUIRES_SQUARED_SELECTION:
			/* Passed arguments are X,Y of the center of the square and the size of
			 * the square. */
			if (com.selection.w > com.selection.h)
				max = com.selection.w;
			else
				max = com.selection.h;

			reg_args->selection.x = com.selection.x + com.selection.w / 2 - max / 2;
			reg_args->selection.w = max;
			reg_args->selection.y = com.selection.y + com.selection.h / 2 - max / 2;
			reg_args->selection.h = max;
			compute_fitting_selection(&reg_args->selection, 2, 2, 1);

			/* save it back to com.selection do display it properly */
			memcpy(&com.selection, &reg_args->selection, sizeof(rectangle));
			fprintf(stdout, "final area: %d,%d,\t%dx%d\n", reg_args->selection.x,
					reg_args->selection.y, reg_args->selection.w,
					reg_args->selection.h);
			redraw(REDRAW_OVERLAY);
			break;
	}
}

/* callback for no output button */
void on_regNoOutput_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *Algo = lookup_widget("ComboBoxRegInter");
	GtkWidget *Prefix = lookup_widget("regseqname_entry");

	gboolean toggled = gtk_toggle_button_get_active(togglebutton);

	gtk_widget_set_sensitive(Algo, !toggled);
	gtk_widget_set_sensitive(Prefix, !toggled);

	keep_noout_state = toggled;
}

void on_regfollowStar_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	update_reg_interface(TRUE);
}

/* callback for 'Go register' button, GTK thread */
void on_seqregister_button_clicked(GtkButton *button, gpointer user_data) {
	struct registration_args *reg_args;
	struct registration_method *method;
	char *msg;
	GtkToggleButton *follow, *matchSel, *x2upscale, *cumul;
	GtkComboBox *cbbt_layers, *reg_all_sel_box;
	GtkComboBoxText *ComboBoxRegInter, *ComboBoxTransfo, *ComboBoxMaxStars, *ComboBoxFraming;
	GtkSpinButton *minpairs, *percent_moved;

	if (!reserve_thread()) {	// reentrant from here
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (!com.seq.regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		// means that a call to seq_check_basic_data() or
		// check_or_allocate_regparam() is missing somewhere else
		unreserve_thread();
		return;
	}

	method = get_selected_registration_method();

	if (com.selection.w <= 0 && com.selection.h <= 0
			&& method->sel != REQUIRES_NO_SELECTION) {
		msg = siril_log_message(
				_("All prerequisites are not filled for registration. Select a rectangle first.\n"));
		siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), msg);
		unreserve_thread();
		return;
	}

	reg_args = calloc(1, sizeof(struct registration_args));

	control_window_switch_to_tab(OUTPUT_LOGS);

	/* filling the arguments for registration */
	follow = GTK_TOGGLE_BUTTON(lookup_widget("followStarCheckButton"));
	matchSel = GTK_TOGGLE_BUTTON(lookup_widget("checkStarSelect"));
	x2upscale = GTK_TOGGLE_BUTTON(lookup_widget("upscaleCheckButton"));
	cbbt_layers = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
	ComboBoxRegInter = GTK_COMBO_BOX_TEXT(lookup_widget("ComboBoxRegInter"));
	cumul = GTK_TOGGLE_BUTTON(lookup_widget("check_button_comet"));
	minpairs = GTK_SPIN_BUTTON(lookup_widget("spinbut_minpairs"));
	percent_moved = GTK_SPIN_BUTTON(lookup_widget("spin_kombat_percent"));
	ComboBoxMaxStars = GTK_COMBO_BOX_TEXT(lookup_widget("comboreg_maxstars"));
	ComboBoxTransfo = GTK_COMBO_BOX_TEXT(lookup_widget("comboreg_transfo"));
	ComboBoxFraming = GTK_COMBO_BOX_TEXT(lookup_widget("comboreg_framing"));
	reg_all_sel_box = GTK_COMBO_BOX(GTK_COMBO_BOX_TEXT(lookup_widget("reg_sel_all_combobox")));
	reg_args->clamp = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_reg_clamp")));

	reg_args->func = method->method_ptr;
	reg_args->seq = &com.seq;
	reg_args->reference_image = sequence_find_refimage(&com.seq);
	reg_args->follow_star = gtk_toggle_button_get_active(follow);
	reg_args->matchSelection = gtk_toggle_button_get_active(matchSel);
	reg_args->no_output = keep_noout_state;
	reg_args->x2upscale = gtk_toggle_button_get_active(x2upscale);
	reg_args->cumul = gtk_toggle_button_get_active(cumul);
	reg_args->prefix = gtk_entry_get_text(GTK_ENTRY(lookup_widget("regseqname_entry")));
	reg_args->min_pairs = gtk_spin_button_get_value_as_int(minpairs);
	reg_args->percent_moved = (float) gtk_spin_button_get_value(percent_moved) / 100.f;
	int starmaxactive = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxMaxStars));
	reg_args->max_stars_candidates = (starmaxactive == -1) ? MAX_STARS_FITTED : maxstars_values[starmaxactive];
	reg_args->type = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxTransfo));
	reg_args->framing = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxFraming));
#ifndef HAVE_CV44
	if (reg_args->type == SHIFT_TRANSFORMATION) {
		siril_log_color_message(_("Shift-only registration is only possible with OpenCV 4.4\n"), "red");
		free(reg_args);
		unreserve_thread();
		return;
	}
#endif
	if (method->method_ptr == register_apply_reg) {
		get_reg_sequence_filtering_from_gui(
				&reg_args->filtering_criterion, &reg_args->filtering_parameter);
	} else {
		reg_args->filters.filter_included = gtk_combo_box_get_active(reg_all_sel_box);
	}

	/* We check that available disk space is enough when
	the registration method produces a new sequence
	*/
	if (!reg_args->no_output && method->method_ptr == register_star_alignment) {
		// first, remove the files that we are about to create
		remove_prefixed_sequence_files(reg_args->seq, reg_args->prefix);

		int nb_frames = reg_args->filters.filter_included ? reg_args->seq->selnum : reg_args->seq->number;
		gint64 size = seq_compute_size(reg_args->seq, nb_frames, get_data_type(reg_args->seq->bitpix));
		if (reg_args->x2upscale)
			size *= 4;
		if (test_available_space(size)) {
			siril_log_color_message(_("Not enough space to save the output images, aborting\n"), "red");
			free(reg_args);
			unreserve_thread();
			return;
		}
	} else if (method->method_ptr == register_comet) {
		pointf velocity = get_velocity();
		if ((velocity.x == 0.0 && velocity.y == 0.0)
				|| isinf(velocity.x) || isinf(velocity.y)) {
			msg = siril_log_color_message(_("The object is not moving, please check your registration data.\n"), "red");
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), msg);
			free(reg_args);
			unreserve_thread();
			return;
		}
	}
	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	reg_args->interpolation = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxRegInter));
	get_the_registration_area(reg_args, method);	// sets selection
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE; // only TRUE for global registration. Will be updated in this case

	if (method->method_ptr == register_star_alignment) { // seqpplyreg case is dealt with in the sanity checks of the method
		if (reg_args->interpolation == OPENCV_NONE && (reg_args->x2upscale || com.seq.is_variable)) {
			siril_log_color_message(_("When interpolation is set to None, the images must be of same size and no upscaling can be applied. Aborting\n"), "red");
			free(reg_args);
			unreserve_thread();
			return;
		}
	}
	if (((method->method_ptr == register_star_alignment || method->method_ptr == register_3stars || method->method_ptr == register_apply_reg) &&
		(reg_args->interpolation == OPENCV_AREA || reg_args->interpolation == OPENCV_LINEAR || reg_args->interpolation == OPENCV_NEAREST || reg_args->interpolation == OPENCV_NONE)) ||
		reg_args->no_output)
		reg_args->clamp = FALSE;

	if (method->method_ptr != register_3stars) clear_stars_list(TRUE); //to avoid problems with com.stars later on in the process

	msg = siril_log_color_message(_("Registration: processing using method: %s\n"),
			"green", method->name);
	msg[strlen(msg) - 1] = '\0';
	if (reg_args->clamp)
		siril_log_message(_("Interpolation clamping active\n"));
	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_reserved_thread(register_thread_func, reg_args);
}

// worker thread function for the registration
gpointer register_thread_func(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	int retval;
	args->seq->reg_invalidated = TRUE;

	args->retval = args->func(args);

	if (args->seq->reference_image == -1) {
		// set new reference image: should we do it all the time?
		// also done in generated sequence in global.c
		args->seq->reference_image = sequence_find_refimage(args->seq);
	}
	if (!args->retval) writeseqfile(args->seq);
	retval = args->retval;
	if (!siril_add_idle(end_register_idle, args)) {
		free_sequence(args->seq, TRUE);
		free(args);
	}
	return GINT_TO_POINTER(retval);
}

// end of registration, GTK thread. Executed when started from the GUI and in
// the graphical command line but not from a script (headless mode)
static gboolean end_register_idle(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	stop_processing_thread();

	if (!args->retval) {
		if (!args->load_new_sequence && sequence_is_loaded()) {
			int chan = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboboxreglayer")));
			update_seqlist(chan);
			fill_sequence_list(args->seq, chan, FALSE);
			set_layers_for_registration();	// update display of available reg data
		}
		else {
			check_seq();
			update_sequences_list(args->new_seq_name);
		}
	}
	set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);
	args->seq->reg_invalidated = FALSE;
	drawPlot();
	update_stack_interface(TRUE);
	adjust_sellabel();

	set_cursor_waiting(FALSE);
	if (args->func == &register_3stars) reset_3stars();

	free(args->new_seq_name);
	free(args);
	return FALSE;
}

/* Moves the selection x, and y after transformation by Href^-1*Him */
void selection_H_transform(rectangle *selection, Homography Href, Homography Himg) {
	double xc, yc;
	xc = (double)selection->x + (double)selection->w * 0.5;
	yc = (double)selection->y + (double)selection->h * 0.5;
	cvTransfPoint(&xc, &yc, Href, Himg);
	selection->x = (int)(xc - (double)selection->w * 0.5);
	selection->y = (int)(yc - (double)selection->h * 0.5);
	// siril_log_message(_("boxselect %d %d %d %d\n"), selection->x, selection->y, selection->w, selection->h);
}

void translation_from_H(Homography H, double *dx, double *dy) {
	*dx = H.h02;
	*dy = -H.h12;
}

Homography H_from_translation(double dx, double dy) {
	Homography H;
	cvGetEye(&H);
	H.h02 = dx;
	H.h12 = -dy;
	return H;
}

void SetNullH(Homography *H) {
	cvGetEye(H);
	H->h00 = 0.0;
	H->h11 = 0.0;
	H->h22 = 0.0;
}
// this transform only works if the source and dest fits have the same size
int shift_fit_from_reg(fits *fit, Homography H) {
	fits *destfit = NULL;
	if (new_fit_image(&destfit, fit->rx, fit->ry, fit->naxes[2], fit->type)) {
		return 1;
	}
	destfit->bitpix = fit->bitpix;
	destfit->orig_bitpix = fit->orig_bitpix;
	int nbpix = fit->naxes[0] * fit->naxes[1];
	if (destfit->type == DATA_FLOAT) {
		memset(destfit->fdata, 0, nbpix * fit->naxes[2] * sizeof(float));
		if (fit->naxes[2] == 3) {
			destfit->fpdata[1] = destfit->fdata + nbpix;
			destfit->fpdata[2] = destfit->fdata + nbpix * 2;
		}
	} else {
		memset(destfit->data, 0, nbpix * fit->naxes[2] * sizeof(WORD));
		if (fit->naxes[2] == 3) {
			destfit->pdata[1] = destfit->data + nbpix;
			destfit->pdata[2] = destfit->data + nbpix * 2;
		}
	}
	copy_fits_metadata(fit, destfit);
	destfit->rx = destfit->naxes[0] = fit->rx;
	destfit->ry = destfit->naxes[1] = fit->ry;
	int shiftx, shifty;
	/* load registration data for current image */
	double dx, dy;
	translation_from_H(H, &dx, &dy);
	shiftx = round_to_int(dx);
	shifty = round_to_int(dy);
	for (int layer = 0; layer < fit->naxes[2]; ++layer) {
		for (int y = 0; y < destfit->ry; ++y) {
			for (int x = 0; x < destfit->rx; ++x) {
				int nx = x + shiftx;
				int ny = y + shifty;
				if (nx >= 0 && nx < destfit->rx && ny >= 0 && ny < destfit->ry) {
					if (destfit->type == DATA_USHORT) {
						destfit->pdata[layer][nx + ny * destfit->rx] = fit->pdata[layer][x + y * fit->rx];
					} else if (destfit->type == DATA_FLOAT) {
						destfit->fpdata[layer][nx + ny * destfit->rx] = fit->fpdata[layer][x + y * fit->rx];
					}
				}
			}
		}
	}
	copyfits(destfit, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	clearfits(destfit);
	return 0;
}

void on_comboreg_framing_changed(GtkComboBox *box, gpointer user_data) {
	gchar *name;
	GtkImage *image = GTK_IMAGE(lookup_widget("framing-image"));
	int i = gtk_combo_box_get_active(box);

	if (i >= 0 && i < G_N_ELEMENTS(reg_frame_registration)) {
		name = g_build_filename(siril_get_system_data_dir(), "pixmaps", reg_frame_registration[i], NULL);
		gtk_image_set_from_file(image, name);

		g_free(name);
	}
}
