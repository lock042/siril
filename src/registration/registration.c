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



#include "core/proto.h"
#include "core/siril_log.h"
#include "gui/image_display.h"
#include "gui/registration.h"
#include "opencv/opencv.h"
#include "drizzle/cdrizzleutil.h"
#include "algos/siril_wcs.h"

int get_registration_layer(const sequence *seq) {
	if (!com.script && seq == &com.seq) {
		return get_registration_layer_from_GUI(seq);
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

gboolean layer_has_registration(const sequence *seq, int layer) {
	if (!seq || layer < 0 || !seq->regparam || seq->nb_layers < 0 || layer >= seq->nb_layers || !seq->regparam[layer] ) return FALSE;
	return TRUE;
}

gboolean layer_has_usable_registration(sequence *seq, int layer) {
	transformation_type min, max;
	guess_transform_from_seq(seq, layer, &min, &max, FALSE); // will check first that layer_has_registration
	if (max <= IDENTITY_TRANSFORMATION)
		return FALSE; // max <= -1 means all H matrices are identity or null
	return TRUE;
}

int seq_has_any_regdata(const sequence *seq) {
	if (!seq || !seq->regparam || seq->nb_layers < 0)
		return -1;
	int i;
	for (i = 0; i < seq->nb_layers; i++)
		if (seq->regparam[i])
			return i;
	return -1;
}

gboolean layer_has_distortion(const sequence *seq, int layer) {
	if (!seq || layer < 0 || !seq->distoparam || seq->nb_layers < 0 || layer >= seq->nb_layers || seq->distoparam[layer].index == DISTO_UNDEF) return FALSE;
	return TRUE;
}

gboolean seq_has_any_distortion(const sequence *seq) {
	for (int i = 0; i < seq->nb_layers; i++) {
		if (layer_has_distortion(seq, i))
			return TRUE;
	}
	return FALSE;
}

regdata *registration_get_current_regdata(struct registration_args *regargs) {
	regdata *current_regdata;
	if (regargs->seq->regparam[regargs->layer]) {
		siril_log_message(_("Recomputing already existing registration for this layer\n"));
		current_regdata = regargs->seq->regparam[regargs->layer];
		/* we reset all values as we may register different images */
		memset(current_regdata, 0, regargs->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(regargs->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		regargs->seq->regparam[regargs->layer] = current_regdata;
	}
	return current_regdata;
}

void create_output_sequence_for_registration(struct registration_args *args, int refindex) {
	sequence seq = { 0 };
	initialize_sequence(&seq, TRUE);

	/* we are not interested in the whole path */
	gchar *seqname = g_path_get_basename(args->seq->seqname);
	char *rseqname = malloc(strlen(args->prefix) + strlen(seqname) + 5);
	sprintf(rseqname, "%s%s.seq", args->prefix, seqname);
	g_unlink(rseqname);	// remove previous to overwrite
	args->new_seq_name = remove_ext_from_filename(rseqname);
	free(rseqname);
	seq.seqname = strdup(args->new_seq_name);
	seq.number = args->new_total;
	seq.selnum = args->new_total;
	seq.fixed = args->seq->fixed;
	seq.nb_layers = (args->driz && args->driz->is_bayer) ? 3 : args->seq->nb_layers;
	seq.imgparam = args->imgparam;
	seq.regparam = calloc(seq.nb_layers, sizeof(regdata*));
	seq.regparam[args->layer] = args->regparam;
	seq.beg = seq.imgparam[0].filenum;
	seq.end = seq.imgparam[seq.number-1].filenum;
	seq.type = args->seq->type;
	seq.current = -1;
	seq.is_variable = check_seq_is_variable(&seq);
	if (!seq.is_variable) {
		seq.rx = args->seq->rx;
		seq.ry = args->seq->ry;
	}
	seq.fz = com.pref.comp.fits_enabled;
	// don't copy from old sequence, it may not be the same image
	if (refindex == -1)
		seq.reference_image = sequence_find_refimage(&seq); //global
	else
		seq.reference_image = refindex; //applyreg
	seq.needs_saving = TRUE;
	writeseqfile(&seq);
	g_free(seqname);
	free_sequence(&seq, FALSE);
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

void get_the_registration_area(struct registration_args *regargs, const struct registration_method *method) {
	int max;
	switch (method->sel) {
		/* even in the case of REQUIRES_NO_SELECTION selection is needed for MatchSelection of starAlignment */
		case REQUIRES_NO_SELECTION:
		case REQUIRES_ANY_SELECTION:
			memcpy(&regargs->selection, &com.selection, sizeof(rectangle));
			break;
		case REQUIRES_SQUARED_SELECTION:
			/* Passed arguments are X,Y of the center of the square and the size of
			 * the square. */
			if (com.selection.w > com.selection.h)
				max = com.selection.w;
			else
				max = com.selection.h;

			regargs->selection.x = com.selection.x + com.selection.w / 2 - max / 2;
			regargs->selection.w = max;
			regargs->selection.y = com.selection.y + com.selection.h / 2 - max / 2;
			regargs->selection.h = max;
			compute_fitting_selection(&regargs->selection, 2, 2, 1);

			/* save it back to com.selection do display it properly */
			memcpy(&com.selection, &regargs->selection, sizeof(rectangle));
			fprintf(stdout, "final area: %d,%d,\t%dx%d\n", regargs->selection.x,
					regargs->selection.y, regargs->selection.w,
					regargs->selection.h);
			redraw(REDRAW_OVERLAY);
			break;
	}
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
	if (!args->retval)
		writeseqfile(args->seq);
	retval = args->retval;
	if (args->disto) {
		free_disto_args(args->disto);
		free(args->disto);
	}
	if (args->driz) {
		free(args->driz);
	}
	if (args->reference_date)
		g_date_time_unref(args->reference_date);
	if (args->wcsref)
		wcsfree(args->wcsref);
	if (!siril_add_idle(end_register_idle, args)) {
		free_sequence(args->seq, TRUE);
		free(args);
	}
	return GINT_TO_POINTER(retval);
}

/* Moves the selection x, and y after transformation by Href^-1*Him */
void selection_H_transform(rectangle *selection, Homography Href, Homography Himg) {
	double xc = selection->x + selection->w * 0.5;
	double yc = selection->y + selection->h * 0.5;
	cvTransfPoint(&xc, &yc, Href, Himg, 1.);
	selection->x = round_to_int(xc - selection->w * 0.5);
	selection->y = round_to_int(yc - selection->h * 0.5);
	siril_debug_print("boxselect %d %d %d %d\n",
			selection->x, selection->y, selection->w, selection->h);
}

void translation_from_H(Homography H, double *dx, double *dy) {
	*dx = H.h02;
	*dy = -H.h12;
}

Homography H_from_translation(double dx, double dy) {
	Homography H = { 0 }; // cvGetEye() cannot fail, but doesn't initialize H.pair_matched,
			      // hence it is initialized here before the call to cvGetEye()
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
	copy_fits_metadata(destfit, fit);
	clearfits(destfit);
	return 0;
}

struct registration_method *new_reg_method(const char *name, registration_function f,
		selection_type s, registration_type t) {
	struct registration_method *reg = malloc(sizeof(struct registration_method));
	reg->name = strdup(name);
	reg->method_ptr = f;
	reg->sel = s;
	reg->type = t;
	return reg;
}
