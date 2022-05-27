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

#include <stdlib.h>
#include <string.h>
#include "core/siril.h"
#include "core/proto.h"

#include "registration.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "core/processing.h"
#include "opencv/opencv.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/utils.h"

static int awaiting_star = 0;
static int selected_stars = 0;

static GtkWidget *three_buttons[3] = { 0 };
static GtkWidget *go_register = NULL;
static GtkLabel *labelreginfo = NULL;
static GtkImage *image_3stars[3] = { NULL };

struct _3psf {
	psf_star *stars[3];
};

static struct _3psf *results;
static int results_size;

// local functions
static int _3stars_alignment(struct registration_args *regargs, regdata *current_regdata);

static void set_registration_ready(gboolean ready) {
	if (!go_register)
		go_register = lookup_widget("goregister_button");
	gtk_widget_set_sensitive(go_register, ready);
}

static void update_label(gchar* str) {
	if (!labelreginfo)
		labelreginfo = GTK_LABEL(lookup_widget("labelregisterinfo"));
	gtk_label_set_text(labelreginfo, str);
}

static void update_icons(int idx, gboolean OK) {
	if (!image_3stars[0]) {
		image_3stars[0] = GTK_IMAGE(lookup_widget("3stars-image1"));
		image_3stars[1] = GTK_IMAGE(lookup_widget("3stars-image2"));
		image_3stars[2] = GTK_IMAGE(lookup_widget("3stars-image3"));
	}
	gtk_image_set_from_icon_name(image_3stars[idx],
			OK ? "gtk-yes" : "gtk-no", GTK_ICON_SIZE_LARGE_TOOLBAR);

}

static void reset_icons() {
	for (int i = 0; i < 3; i++) {
		update_icons(i, FALSE);
	}
}

static int _3stars_seqpsf_finalize_hook(struct generic_seq_args *args) {
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;

	if (args->retval) {
		if (args->seq->current != 0)
			update_label(_("Make sure you load the first image"));
		else update_label(_("Star analysis failed"));
		goto psf_end;
	}

	GSList *iterator;
	for (iterator = spsfargs->list; iterator; iterator = iterator->next) {
		struct seqpsf_data *data = iterator->data;
		results[data->image_index].stars[awaiting_star - 1] = data->psf;
	}
	g_slist_free(spsfargs->list);

	//should not happen as the stars were confirmed through psf on the ref image
	int refimage = sequence_find_refimage(&com.seq);
	if (!results[refimage].stars[awaiting_star - 1]) {
		siril_log_color_message(_("The star was not found in the reference image. Change the selection or the reference image\n"), "red");
		for (int i = 0 ; i < com.seq.number; i++)
			results[i].stars[awaiting_star - 1] = NULL;
		args->retval = 1;
		goto psf_end;
	}

	com.stars = realloc(com.stars, 4 * sizeof(psf_star *)); // to be sure...
	com.stars[awaiting_star - 1] = duplicate_psf(results[args->seq->current].stars[awaiting_star - 1]);


psf_end:
	free(spsfargs);
	return args->retval;
}

static int start_seqpsf(struct registration_args *regargs) {
	struct seqpsf_args *spsfargs = malloc(sizeof(struct seqpsf_args));
	spsfargs->for_registration = TRUE; // if false, photometry is computed
	spsfargs->framing = FOLLOW_STAR_FRAME;
	spsfargs->list = NULL;	// GSList init is NULL
	struct generic_seq_args *args = calloc(1, sizeof(struct generic_seq_args));
	if (!regargs->process_all_frames) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->seq = &com.seq;
	args->partial_image = TRUE;
	args->layer_for_partial = get_registration_layer(&com.seq);
	args->regdata_for_partial = FALSE;
	args->get_photometry_data_for_partial = FALSE;
	args->image_hook = seqpsf_image_hook;
	args->finalize_hook = _3stars_seqpsf_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("PSF on area for 2 or 3 stars");
	args->upscale_ratio = 1.0;
	args->user = spsfargs;
	args->already_in_a_thread = TRUE;
	args->parallel = FALSE;	// follow star implies not parallel
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	if (!results) {
		results = calloc(com.seq.number, sizeof(struct _3psf));
		if (!results) {
			PRINT_ALLOC_ERR;
			free(spsfargs);
			free(args);
			return 1;
		}
		results_size = com.seq.number;
	}

	generic_sequence_worker(args);
	regargs->retval = args->retval;
	free(args);
	return regargs->retval;
}

void on_select_star_button_clicked(GtkButton *button, gpointer user_data) {
	if (!three_buttons[0]) {
		three_buttons[0] = lookup_widget("pickstar1");
		three_buttons[1] = lookup_widget("pickstar2");
		three_buttons[2] = lookup_widget("pickstar3");
	}
	if (!com.selection.w || !com.selection.h) {
		update_label(_("Draw a selection around the star"));
		return;
	}
	/* for now we force the first image to be loaded - will need to update
	this if we allow to use exsting registration data
	We may want to force the first selected image instead, in case we need
	to discard a few images at the beginning of the series
	*/
	if (com.seq.current != 0) { 
		update_label(_("Make sure you load the first image"));
		return;
	}

	GtkWidget *widget = GTK_WIDGET(button);
	if (three_buttons[0] == widget)
		awaiting_star = 1;
	else if (three_buttons[1] == widget)
		awaiting_star = 2;
	else if (three_buttons[2] == widget)
		awaiting_star = 3;
	else {
		fprintf(stderr, "unknown button clicked\n");
		return;
	}

	if (!com.stars)
		com.stars = calloc(4, sizeof(psf_star *)); // don't use new_psf_star. It is a bit different

	int index;
	int layer = get_registration_layer(&com.seq);
	add_star(&gfit, layer, &index);
	if (index == -1) {
		update_label(_("No star found, make another selection"));
	} else {
		selected_stars ++;
		unset_suggested(three_buttons[awaiting_star - 1]);
		if (awaiting_star < 3) set_suggested(three_buttons[awaiting_star]);
		update_icons(awaiting_star - 1, TRUE);
		delete_selected_area();
		set_registration_ready((awaiting_star >= 2) ? TRUE : FALSE);
	}

}

int register_3stars(struct registration_args *regargs) {

	// we need first to run seqpsf for the 2 or 3 stars selected in the first images
	// and then to proceed with the registration (matching + saving images if !no_ouput)
	// for the selection, we use the com.stars x/y pos to redraw a box
	for (int i = 0; i < selected_stars; i++) {
		// TODO - save initial selection in case the boxes were made larger to avoid follow_star
		delete_selected_area();
		com.selection.w = (int)com.stars[i]->fwhmx * 4;
		com.selection.h = (int)com.stars[i]->fwhmx * 4;
		com.selection.x = com.stars[i]->xpos - com.selection.w * 0.5;
		com.selection.y = com.stars[i]->ypos - com.selection.w * 0.5;
		new_selection_zone();
		awaiting_star = i + 1;
		siril_log_color_message(_("Processing star #%d\n"), "salmon", awaiting_star);
		if (start_seqpsf(regargs)) return 1;
	}
	delete_selected_area();

	int refimage = regargs->reference_image;
	if (!results[refimage].stars[0] || !results[refimage].stars[1]) {
		siril_log_color_message(_("Less than two stars were found in the reference image, try setting another as reference?\n"), "red");
		return 1;
	}

	regdata *current_regdata = star_align_get_current_regdata(regargs);
	if (!current_regdata) return -2;

	/* set regparams for current sequence before closing it */
	for (int i = 0; i < regargs->seq->number; i++) {
		double sumx = 0.0, sumy = 0.0, sumb = 0.0;
		int nb_stars = 0;
		if (results[i].stars[0]) {
			sumx += results[i].stars[0]->fwhmx;
			sumy += results[i].stars[0]->fwhmy;
			sumb += results[i].stars[0]->B;
			nb_stars++;
		}
		if (results[i].stars[1]) {
			sumx += results[i].stars[1]->fwhmx;
			sumy += results[i].stars[1]->fwhmy;
			sumb += results[i].stars[1]->B;
			nb_stars++;
		}
		if (results[i].stars[2]) {
			sumx += results[i].stars[2]->fwhmx;
			sumy += results[i].stars[2]->fwhmy;
			sumb += results[i].stars[2]->B;
			nb_stars++;
		}
		if (nb_stars >= 2) {
			double fwhm = sumx / nb_stars;
			current_regdata[i].roundness = sumy / sumx;
			current_regdata[i].fwhm = fwhm;
			current_regdata[i].weighted_fwhm = 2. * fwhm * (double)(selected_stars - nb_stars) / (double)nb_stars + fwhm; 
			current_regdata[i].background_lvl = sumb / nb_stars;
			current_regdata[i].number_of_stars = nb_stars;
		}
	}

	return _3stars_alignment(regargs, current_regdata);
}

/* image sequence processing */
static int _3stars_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int refimage = regargs->reference_image;
	Homography H = { 0 };
	if (regargs->no_output) {
		/* if "save transformation only", we choose to initialize all frames
		 * to exclude status. If registration is ok, the status is
		 * set to include */
		args->seq->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
	}
	int nb_stars = sadata->current_regdata[in_index].number_of_stars;
	if (in_index != refimage) {
		if (nb_stars < 2) return 1;
		struct s_star *arrayref, *arraycur, *starsin, *starsout;
		arrayref = (s_star *) shMalloc(nb_stars * sizeof(s_star));
		arraycur = (s_star *) shMalloc(nb_stars * sizeof(s_star));
		int j = 0;
		unsigned char *mask;
		for (int i = 0; i < 3; i++) {
			starsin = &(arrayref[j]);
			starsout = &(arraycur[j]);
			g_assert(starsin != NULL);
			g_assert(starsout != NULL);
			if (results[in_index].stars[i] != NULL && results[refimage].stars[i] != NULL) {
				starsin->x = results[refimage].stars[i]->xpos;
				starsin->y = results[refimage].stars[i]->ypos;
				starsout->x = results[in_index].stars[i]->xpos;
				starsout->y = results[in_index].stars[i]->ypos;
				j++;
			}
		}
		mask = cvCalculH(arraycur, arrayref, nb_stars, &H, SIMILARITY_TRANSFORMATION);
		if (!mask || H.Inliers < 2) {
			siril_log_color_message(_("Cannot perform star matching: Image %d skipped\n"), "red",  args->seq->imgparam[in_index].filenum);
			return 1;
		}
		sadata->current_regdata[in_index].H = H;

		if (!regargs->no_output) {
			if (regargs->interpolation <= OPENCV_LANCZOS4) {
				if (cvTransformImage(fit, sadata->ref.x, sadata->ref.y, H, regargs->x2upscale, regargs->interpolation)) {
					return 1;
				}
			} else {
				fits *destfit = NULL;
				if (new_fit_image(&destfit, fit->rx, fit->ry, fit->naxes[2], fit->type)) {
					return 1;
				}
				destfit->bitpix = fit->type;
				destfit->orig_bitpix = fit->orig_bitpix;
				int nbpix = fit->naxes[0] * fit->naxes[1] * (regargs->x2upscale ? 4 : 1);
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
				double scale = regargs->x2upscale ? 2. : 1.;
				destfit->rx = destfit->naxes[0] = fit->rx * scale;
				destfit->ry = destfit->naxes[1] = fit->ry * scale;
				int shiftx, shifty;
				/* load registration data for current image */
				double dx, dy;
				translation_from_H(H, &dx, &dy);
				shiftx = round_to_int(dx * scale);
				shifty = round_to_int(dy * scale);
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
			}
		}
	} else {
		// reference image
		cvGetEye(&H);
		sadata->current_regdata[in_index].H = H;
		if (regargs->x2upscale && !regargs->no_output) {
			if (cvResizeGaussian(fit, fit->rx * 2, fit->ry * 2, OPENCV_NEAREST))
				return 1;
		}
	}

	if (!regargs->no_output) {
		regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
		regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
		regargs->imgparam[out_index].rx = args->seq->imgparam[in_index].rx;
		regargs->imgparam[out_index].ry = args->seq->imgparam[in_index].ry;
		regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;
		regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
		regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;
		regargs->regparam[out_index].background_lvl = sadata->current_regdata[in_index].background_lvl;
		regargs->regparam[out_index].number_of_stars = sadata->current_regdata[in_index].number_of_stars;
		cvGetEye(&regargs->regparam[out_index].H);

		if (regargs->x2upscale) {
			fit->pixel_size_x /= 2;
			fit->pixel_size_y /= 2;
			regargs->regparam[out_index].fwhm *= 2.0;
			regargs->regparam[out_index].weighted_fwhm *= 2.0;
		}
	} else {
		// TODO: check if H matrix needs to include a flip or not based on fit->top_down
		// seems like not but this could backfire at some point
		args->seq->imgparam[in_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	}
	sadata->success[out_index] = 1;
	return 0;
}

static int _3stars_align_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, FALSE,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	unsigned int required = MB_per_scaled_image;
	if (limit > 0) {
		/* The registration memory consumption, n is original image size:
		 * Monochrome: O(n) for loaded image, O(nscaled) for output image,
		 *             so O(2n) for unscaled, O(n+nscaled) for scaled
		 * Color:
		 * 	allocations				sum
		 *	O(n) for loaded image			O(n) as input
		 *	O(n) for bgr image			O(2n)
		 *	-O(n) for input				O(n)
		 *	O(nscaled) for output			O(n+nscaled)
		 *	-O(n) for bgr				O(nscaled)
		 *	O(nscaled) pour alloc de data		O(2nscaled)
		 *	-O(nscaled) for output			O(nscaled) as output
		 * so maximum is O(2n) for unscaled and O(2nscaled) for scaled
		 */
		if (args->upscale_ratio == 1.0)
			required = MB_per_orig_image * 2;
		else if (args->seq->nb_layers == 3)
			required = MB_per_scaled_image * 2;
		else required = MB_per_orig_image + MB_per_scaled_image;

		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_scaled_image;
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
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_scaled_image, limit, for_writer ? "images" : "threads");
#else
		/* we still want the check of limit = 0 above */
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

static int _3stars_alignment(struct registration_args *regargs, regdata *current_regdata) {
	struct generic_seq_args *args = create_default_seqargs(&com.seq);
	args->stop_on_error = FALSE;
	if (!regargs->process_all_frames) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->compute_mem_limits_hook = _3stars_align_compute_mem_limits;
	args->prepare_hook = star_align_prepare_results;
	args->image_hook = _3stars_align_image_hook;
	args->finalize_hook = star_align_finalize_hook;	// from global registration
	args->stop_on_error = FALSE;
	args->description = (!regargs->no_output) ? _("Creating the rotated image sequence") : _("Saving the transformation matrices");
	args->has_output = !regargs->no_output;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = regargs->x2upscale ? 2.0 : 1.0;
	args->new_seq_prefix = regargs->prefix;
	args->load_new_sequence = !regargs->no_output;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		free(args);
		return -1;
	}
	sadata->regargs = regargs;
	// we pass the regdata just to avoid recomputing it for the new sequence
	sadata->current_regdata = current_regdata; 
	args->user = sadata;

	// some prep work done in star_align_prepare_hook for global
	// need to duplicate it here
	int refimage = args->seq->reference_image;
	sadata->ref.x = args->seq->imgparam[refimage].rx;
	sadata->ref.y = args->seq->imgparam[refimage].ry;

	if (regargs->x2upscale) {
		if (regargs->no_output) {
			args->seq->upscale_at_stacking = 2.0;
		} else {
			sadata->ref.x *= 2.0;
			sadata->ref.y *= 2.0;
		}
	}
	else {
		if (regargs->no_output) {
			args->seq->upscale_at_stacking = 1.0;
		}
	}

	generic_sequence_worker(args);
	
	for (int i = 0; i < results_size; i++) {
		for (int s = 0; s < 3; s++)
			if (results[i].stars[s])
				free(results[i].stars[s]);
	}
	free(results);
	results = NULL;
	reset_3stars();
	return args->retval;
}

void reset_3stars(){
    if (!GTK_IS_WIDGET(three_buttons[0])) return;
    reset_icons();
    for (int i = 0; i < 3; i++)
        unset_suggested(three_buttons[i]);
    set_suggested(three_buttons[0]);
    set_registration_ready(FALSE);
	awaiting_star = 0;
	selected_stars = 0;
}

