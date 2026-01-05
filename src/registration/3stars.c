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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "core/siril.h"
#include "core/proto.h"

#include "registration.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "opencv/opencv.h"
#include "gui/image_interactions.h"
#include "gui/utils.h"
#include "gui/PSF_list.h"	// clear_stars_list

static int awaiting_star = 0;
static int selected_stars = 0;

static GtkWidget *three_buttons[3] = { 0 };
static GtkWidget *go_register = NULL;
static GtkLabel *labelreginfo = NULL;
static GtkImage *image_3stars[3] = { NULL };
static GtkWidget *follow = NULL, *onlyshift = NULL;
static GtkComboBox *reg_all_sel_box = NULL, *comboboxreglayer = NULL;

struct _3psf {
	psf_star *stars[3];
};

static struct _3psf *results;
static int results_size;

static rectangle _3boxes[3];

/* UI functions */
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

void reset_3stars(){
	if (!GTK_IS_WIDGET(three_buttons[0])) return;
	reset_icons();
	for (int i = 1; i < 3; i++) {
		unset_suggested(three_buttons[i]);
		gtk_widget_set_sensitive(three_buttons[i], FALSE);
	}
	set_suggested(three_buttons[0]);
	gtk_widget_set_sensitive(three_buttons[0], TRUE);
	set_registration_ready(FALSE);
	clear_stars_list(TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(onlyshift), FALSE);
	gtk_widget_set_sensitive(onlyshift, FALSE);
	awaiting_star = 0;
	selected_stars = 0;
}

int _3stars_get_number_selected_stars() {
	return selected_stars;
}

gboolean _3stars_check_selection() {
	if (com.seq.current < 0)
		return FALSE;

	if (!follow) {
		follow = lookup_widget("followStarCheckButton");
		reg_all_sel_box = GTK_COMBO_BOX(GTK_COMBO_BOX_TEXT(lookup_widget("reg_sel_all_combobox")));
		labelreginfo = GTK_LABEL(lookup_widget("labelregisterinfo"));
		onlyshift = lookup_widget("onlyshift_checkbutton");
	}
	gboolean dofollow = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(follow));
	gboolean doall = !gtk_combo_box_get_active(reg_all_sel_box);

	if (dofollow) {
		if (doall && com.seq.current != 0) {
			gtk_label_set_text(labelreginfo, _("Make sure you load the first image"));
			return FALSE;
		} else if (!doall && com.seq.current != get_first_selected(&com.seq)) {
			gtk_label_set_text(labelreginfo, _("Make sure you load the first selected image"));
			return FALSE;
		}
	}
	if (!doall && !com.seq.imgparam[com.seq.current].incl) {
		update_label(_("Make sure you load an image which is included"));
		return FALSE;
	}
	return TRUE;
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

	if (!_3stars_check_selection()) return;

	GtkWidget *widget = GTK_WIDGET(button);
	if (three_buttons[0] == widget) {
		reset_icons();
		clear_stars_list(TRUE);
		selected_stars = 0;
		awaiting_star = 1;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(onlyshift), TRUE);
		gtk_widget_set_sensitive(onlyshift, FALSE);
	} else if (three_buttons[1] == widget) {
		awaiting_star = 2;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(onlyshift), FALSE);
		gtk_widget_set_sensitive(onlyshift, TRUE);
	} else if (three_buttons[2] == widget) {
		awaiting_star = 3;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(onlyshift), FALSE);
		gtk_widget_set_sensitive(onlyshift, TRUE);
	} else {
		fprintf(stderr, "unknown button clicked\n");
		return;
	}

	if (!com.stars)
		com.stars = calloc(4, sizeof(psf_star *)); // don't use new_psf_star. It is a bit different

	int index;
	if (!comboboxreglayer)
		comboboxreglayer = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "comboboxreglayer"));
	int layer = gtk_combo_box_get_active(comboboxreglayer);
	if (layer < 0) {
		fprintf(stderr, "invalid registration layer\n");
		return;
	}
	add_star(gfit, layer, &index);
	if (index == -1) {
		update_label(_("No star found, make another selection"));
	} else {
		memcpy(&_3boxes[selected_stars], &com.selection, sizeof(rectangle));
		selected_stars ++;
		unset_suggested(three_buttons[awaiting_star - 1]);
		gtk_widget_set_sensitive(three_buttons[awaiting_star - 1], FALSE);
		if (awaiting_star < 3) {
			set_suggested(three_buttons[awaiting_star]);
			gtk_widget_set_sensitive(three_buttons[awaiting_star], TRUE);
		}
		update_icons(awaiting_star - 1, TRUE);
		delete_selected_area();
	}
}

/* seqpsf hooks and main process */
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

	int refimage = sequence_find_refimage(&com.seq);
	if (!results[refimage].stars[awaiting_star - 1]) {
		siril_log_color_message(_("The star was not found in the reference image. Change the selection or the reference image\n"), "red");
		for (int i = 0 ; i < com.seq.number; i++)
			results[i].stars[awaiting_star - 1] = NULL;
		args->retval = 1;
		goto psf_end;
	}

	com.stars = realloc(com.stars, 4 * sizeof(psf_star *)); // to be sure...
	com.stars[3] = NULL;
	com.stars[awaiting_star - 1] = duplicate_psf(results[args->seq->current].stars[awaiting_star - 1]);

psf_end:
	g_slist_free(spsfargs->list);
	free(spsfargs);
	args->user = NULL;
	return args->retval;
}

static void _3stars_free_results() {
	if (!results) return;
	for (int i = 0; i < results_size; i++) {
		for (int s = 0; s < 3; s++)
			if (results[i].stars[s])
				free(results[i].stars[s]);
	}
	free(results);
	results = NULL;
}

static int _3stars_seqpsf(struct registration_args *regargs) {
	struct seqpsf_args *spsfargs = calloc(1, sizeof(struct seqpsf_args));
	struct generic_seq_args *args = calloc(1, sizeof(struct generic_seq_args));
	spsfargs->for_photometry = FALSE;
	spsfargs->allow_use_as_regdata = BOOL_FALSE;
	spsfargs->list = NULL;	// GSList init is NULL
	spsfargs->framing = (regargs->follow_star) ? FOLLOW_STAR_FRAME : REGISTERED_FRAME;
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	// making sure we can use registration data - maybe we could have done that beforehand...
	if (spsfargs->framing == REGISTERED_FRAME && !layer_has_usable_registration(regargs->seq, regargs->layer)) {
		spsfargs->framing = ORIGINAL_FRAME;
	}
	if (spsfargs->framing == REGISTERED_FRAME) {
		if (regargs->seq->reference_image < 0) regargs->seq->reference_image = sequence_find_refimage(regargs->seq);
		if (guess_transform_from_H(regargs->seq->regparam[regargs->layer][regargs->seq->reference_image].H) == NULL_TRANSFORMATION) {
			siril_log_color_message(_("The reference image has a null matrix and was not previously registered. Please select another one.\n"), "red");
			free(args);
			free(spsfargs);
			return 1;
		}
		if (regargs->seq->current != regargs->seq->reference_image) {
			// transform selection back from current to ref frame coordinates
			if (guess_transform_from_H(regargs->seq->regparam[regargs->layer][regargs->seq->current].H) == NULL_TRANSFORMATION) {
				siril_log_color_message(_("The current image has a null matrix and was not previously registered. Please load another one to select the stars.\n"), "red");
				free(args);
				free(spsfargs);
				return 1;
			}
			selection_H_transform(&args->area, regargs->seq->regparam[regargs->layer][regargs->seq->current].H, regargs->seq->regparam[regargs->layer][regargs->seq->reference_image].H);
			if (args->area.x < 0 || args-> area.x > regargs->seq->rx - args->area.w ||
					args->area.y < 0 || args->area.y > regargs->seq->ry - args->area.h) {
				siril_log_color_message(_("This area is outside of the reference image. Please select the reference image to select another star.\n"), "red");
				free(args);
				free(spsfargs);
				return 1;
			}
		}
	}

	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	} else {
		args->filtering_criterion = seq_filter_all;
		args->nb_filtered_images = regargs->seq->number;
	}
	args->seq = regargs->seq;
	args->partial_image = TRUE;
	args->layer_for_partial = (spsfargs->framing == REGISTERED_FRAME) ? get_registration_layer(&com.seq) : gtk_combo_box_get_active(comboboxreglayer);
	args->regdata_for_partial = spsfargs->framing == REGISTERED_FRAME;
	args->get_photometry_data_for_partial = FALSE;
	args->image_hook = seqpsf_image_hook;
	args->finalize_hook = _3stars_seqpsf_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("PSF on area for 2 or 3 stars");
	args->upscale_ratio = 1.0;
	args->user = spsfargs;
	args->already_in_a_thread = TRUE;
	args->parallel = !regargs->follow_star;	// follow star implies not parallel
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


/*
This function runs seqpsfs on 1/2/3 stars as selected by the user
then computes transformation matrix to the ref image
Registration data is saved to the input sequence
*/
int register_3stars(struct registration_args *regargs) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	int nb_stars_ref = 0;
	int refimage = regargs->reference_image;
	Homography H = { 0 };
	delete_selected_area();
	gboolean onestar = selected_stars == 1;
	// for the selection, we use the com.stars x/y pos to redraw a box
	for (int i = 0; i < selected_stars; i++) {
		// delete_selected_area();
		memcpy(&com.selection, &_3boxes[awaiting_star - 1], sizeof(rectangle));
		// gui_function(new_selection_zone, NULL);
		awaiting_star = i + 1;
		siril_log_color_message(_("Processing star #%d\n"), "salmon", awaiting_star);
		if (_3stars_seqpsf(regargs)) return 1;
		if (results[refimage].stars[i] != NULL) nb_stars_ref++;
		// Determine if it's worth going on, i.e. if enough stars were found in ref image
		// before we proceed with next star
		if (!onestar && ((selected_stars == 2 && nb_stars_ref <= i) || (selected_stars == 3 && i == 1 && nb_stars_ref == 0))) {
			siril_log_color_message(_("Less than two stars were found in the reference image, try setting another as reference?\n"), "red");
			_3stars_free_results();
			return 1;
		}
		if (onestar && nb_stars_ref <= i) {
			siril_log_color_message(_("No star was found in the reference image, try setting another as reference?\n"), "red");
			_3stars_free_results();
			return 1;
		}
	}
	delete_selected_area();

	regdata *current_regdata = registration_get_current_regdata(regargs);
	if (!current_regdata) return -2;

	char *msg;
	msg = siril_log_message(_("Saving the transformation matrices\n"));
	msg[strlen(msg)-1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	int processed = 0, failed = 0;

	// local flag accounting both for process_all_frames flag and collecting failures along the process
	gboolean *included = calloc(regargs->seq->number, sizeof(gboolean));
	if (!included) {
		PRINT_ALLOC_ERR;
		_3stars_free_results();
		return 1;
	}
	float *scores = calloc(regargs->seq->number, sizeof(float));
	if (!scores) {
		PRINT_ALLOC_ERR;
		free(included);
		_3stars_free_results();
		return 1;
	}

	/* set regparams for current sequence before closing it */
	for (int i = 0; i < regargs->seq->number; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included) continue;
		processed++;
		double sumx = 0.0, sumy = 0.0, sumb = 0.0;
		int nb_stars = 0;
		if (!(i % 32)) {
			set_progress_bar_data(NULL, (double)i / regargs->seq->number);
		}

		/* we choose to initialize all frames
		* to exclude status. If registration is ok, the status is
		* set to include */
		regargs->seq->imgparam[i].incl = !SEQUENCE_DEFAULT_INCLUDE;

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
		if ((!onestar && nb_stars >= 2) || (onestar && nb_stars == 1)) {
			double fwhm = sumx / nb_stars;
			current_regdata[i].roundness = sumy / sumx;
			current_regdata[i].fwhm = fwhm;
			current_regdata[i].weighted_fwhm = 2. * fwhm * (double)(nb_stars_ref - nb_stars) / (double)nb_stars + fwhm;
			current_regdata[i].background_lvl = sumb / nb_stars;
			current_regdata[i].number_of_stars = nb_stars;
			included[i] = TRUE;
			scores[i] = current_regdata[i].weighted_fwhm;
		} else {
			siril_log_color_message(_("Cannot perform star matching: Image %d skipped\n"), "red",  regargs->seq->imgparam[i].filenum);
			failed++;
			continue;
		}
	}

	// setting the new reference
	int best_index = minidx(scores, included, regargs->seq->number, NULL);
	regargs->seq->reference_image = best_index;
	int reffilenum = regargs->seq->imgparam[best_index].filenum;	// for display purposes
	siril_log_message(_("Trial #%d: After sequence analysis, we are choosing image %d as new reference for registration\n"), 1, reffilenum);

	// computing the transformation matrices
	for (int i = 0; i < regargs->seq->number; i++) {
		if (!included[i]) continue;
		// Determine number of stars present in both in image and ref
		int nb_stars = 0;
		for (int j = 0; j < selected_stars; j++) {
			if (results[i].stars[j] != NULL && results[refimage].stars[j] != NULL) nb_stars++;
		}
		if (i != refimage) {
			if (regargs->type == SHIFT_TRANSFORMATION) { // shift only 2-3 stars or onestar
				if (nb_stars == 0) {
					siril_log_color_message(_("Cannot perform star matching: Image %d skipped\n"), "red",  regargs->seq->imgparam[i].filenum);
					failed++;
					continue;
				}
				double shiftx = 0., shifty = 0.;
				int k = 0;
				for (int j = 0; j < selected_stars; j++) {
					if (results[i].stars[j] != NULL && results[refimage].stars[j] != NULL) {
						shiftx += results[refimage].stars[j]->xpos - results[i].stars[j]->xpos;
						shifty += results[i].stars[j]->ypos - results[refimage].stars[j]->ypos;
						k++;
					}
				}
				shiftx /= (double)k;
				shifty /= (double)k;
				if (selected_stars > 1 && nb_stars > 1) { // error checking can only be computed if more than one star
					double err = 0., tmp_err = 0.;;
					for (int j = 0; j < selected_stars; j++) {
						if (results[i].stars[j] != NULL && results[refimage].stars[j] != NULL) {
							tmp_err += SQR(results[refimage].stars[j]->xpos - results[i].stars[j]->xpos - shiftx);
							tmp_err += SQR(results[i].stars[j]->ypos - results[refimage].stars[j]->ypos - shifty);
							if (tmp_err > err) err = tmp_err;
						}
					}
					err = pow(err, 0.5);
					if (err > current_regdata[i].fwhm) {
						siril_log_color_message(_("Cannot perform star matching: Image %d skipped\n"), "red",  regargs->seq->imgparam[i].filenum);
						printf("Image %d max_error : %3.2f > fwhm: %3.2f\n", regargs->seq->imgparam[i].filenum, err, current_regdata[i].fwhm);
						failed++;
						continue;
					}
				}
				current_regdata[i].H = H_from_translation(shiftx, shifty);
				fprintf(stderr, "reg: file %d, shiftx=%f shifty=%f\n",
				regargs->seq->imgparam[i].filenum, shiftx, shifty);
			} else { // 2-3 stars reg with rotation
				if (nb_stars < 2) {
					siril_log_color_message(_("Cannot perform star matching: Image %d skipped\n"), "red",  regargs->seq->imgparam[i].filenum);
					failed++;
					continue;
				}
				struct s_star *arrayref, *arraycur, *starsin, *starsout;
				arrayref = (s_star *) shMalloc(nb_stars * sizeof(s_star));
				arraycur = (s_star *) shMalloc(nb_stars * sizeof(s_star));
				int k = 0;
				for (int j = 0; j < selected_stars; j++) {
					starsin = &(arrayref[k]);
					starsout = &(arraycur[k]);
					g_assert(starsin != NULL);
					g_assert(starsout != NULL);
					if (results[i].stars[j] != NULL && results[refimage].stars[j] != NULL) {
						starsin->x = results[refimage].stars[j]->xpos;
						starsin->y = results[refimage].stars[j]->ypos;
						starsout->x = results[i].stars[j]->xpos;
						starsout->y = results[i].stars[j]->ypos;
						k++;
					}
				}
				double err = cvCalculRigidTransform(arrayref, arraycur, nb_stars, &H);
				free(arrayref);
				free(arraycur);
				if (err > current_regdata[i].fwhm) {
					siril_log_color_message(_("Cannot perform star matching: Image %d skipped\n"), "red",  regargs->seq->imgparam[i].filenum);
					printf("Image %d max_error : %3.2f > fwhm: %3.2f\n", regargs->seq->imgparam[i].filenum, err, current_regdata[i].fwhm);
					failed++;
					continue;
				}
				current_regdata[i].H = H;
			}
		} else {
			cvGetEye(&current_regdata[i].H);
		}
		// H computation was sucessful, include the image
		regargs->seq->imgparam[i].incl = SEQUENCE_DEFAULT_INCLUDE;
	}
	// cleaning
	regargs->new_total = processed - failed;
	free(included);
	free(scores);
	_3stars_free_results();

	fix_selnum(regargs->seq, FALSE);
	siril_log_message(_("Registration finished.\n"));
	siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, regargs->new_total);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	return 0;

}
