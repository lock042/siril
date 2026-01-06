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

#include "io/FITS_symlink.h" // needs to be included before siril.h to avoid type redefinition
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "algos/siril_wcs.h"
#include "gui/progress_and_log.h"
#include "registration.h"
#include "opencv/opencv.h"
#include "drizzle/cdrizzleutil.h"

static int new_ref_index = -1;

pointf compute_velocity(GDateTime *t1, GDateTime *t2, point d1, point d2) {
	pointf delta_d, px_per_hour = { 0.f, 0.f };

	if (t1 && t2) {
		float delta_t = g_date_time_difference(t2, t1);
		delta_d.x = (float)(d2.x - d1.x);
		delta_d.y = (float)(d2.y - d1.y);

		px_per_hour.x = delta_d.x / delta_t * 3600000000.f;
		px_per_hour.y = delta_d.y / delta_t * 3600000000.f;
	}

	return px_per_hour;
}

int get_comet_shift(GDateTime *ref, GDateTime *img, pointf px_per_hour, pointf *reg) {
	if (img && ref) {
		float delta_t = (float) g_date_time_difference(img, ref);
		delta_t /= 3600000000.f;
		reg->x = delta_t * px_per_hour.x;
		reg->y = delta_t * px_per_hour.y;
	}
	return 0;
}

/***** generic moving object registration *****/

struct comet_align_data {
	struct registration_args *regargs;
	regdata *current_regdata;
	disto_params *distoparam;
	gboolean flipped;
};

static void create_output_sequence_for_comet(struct registration_args *args, int refindex, disto_params *distoparam) {
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
	seq.number = args->seq->number;
	seq.selnum = args->new_total;
	seq.fixed = args->seq->fixed;
	seq.nb_layers = args->seq->nb_layers;
	seq.imgparam = args->imgparam;
	seq.regparam = calloc(seq.nb_layers, sizeof(regdata*));
	seq.regparam[args->layer] = args->regparam;
	seq.beg = seq.imgparam[0].filenum;
	seq.end = seq.imgparam[seq.number - 1].filenum;
	seq.type = args->seq->type;
	seq.current = -1;
	seq.is_variable = FALSE;
	seq.fz = args->seq->fz;
	seq.reference_image = refindex;
	seq.needs_saving = TRUE;
	if (distoparam)
		seq.distoparam = distoparam;
	fix_selnum(&seq, FALSE);
	writeseqfile(&seq);
	g_free(seqname);
	free_sequence(&seq, FALSE);
}


static int comet_align_prepare_hook(struct generic_seq_args *args) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;

	// allocate destination sequence data
	regargs->imgparam = malloc(args->seq->number * sizeof(imgdata));
	regargs->regparam = calloc(args->seq->number, sizeof(regdata));
	if (!regargs->imgparam  || !regargs->regparam) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	// we copy imgparam from the originating sequence and reset incl flags
	memcpy(regargs->imgparam, args->seq->imgparam, args->seq->number * sizeof(imgdata));
	for (int i = 0; i < args->seq->number; i++)
		regargs->imgparam[i].incl = !SEQUENCE_DEFAULT_INCLUDE; // this will be set by the image hook

	// we copy regparam from the originating sequence if it exists
	if (args->seq->regparam[regargs->layer]) {
		memcpy(regargs->regparam, args->seq->regparam[regargs->layer], args->seq->number * sizeof(regdata));
		cadata->current_regdata = args->seq->regparam[regargs->layer];
	} else {
		cadata->current_regdata = NULL;
	}

	// if FITS, we will symlink to originating images
	// if FITSEQ or SER, we will symlink directly to the orginal file in the finalize hook
	if (args->seq->type == SEQ_REGULAR && seq_prepare_hook(args))
		return 1;

	// we fetch the reference date
	fits ref = { 0 };
	if (seq_read_frame_metadata(args->seq, regargs->reference_image, &ref)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(cadata->current_regdata);
		return 1;
	}
	regargs->reference_date = g_date_time_ref(ref.keywords.date_obs);

	// we must copy the disto data from the originating sequence
	if (layer_has_distortion(args->seq, regargs->layer)) {
		disto_params *disto_orig = args->seq->distoparam;
		cadata->distoparam = calloc(args->seq->nb_layers, sizeof(disto_params));
		if (disto_orig[regargs->layer].index != DISTO_FILES) {
			cadata->distoparam[regargs->layer].index = disto_orig[regargs->layer].index;
			if (cadata->distoparam[regargs->layer].filename)
				cadata->distoparam[regargs->layer].filename = g_strdup(disto_orig[regargs->layer].filename);
		} else {
		// if registration originates from astrometry, we can't keep it as is
		// otherwise, registration will be recomputed from wcs at applyreg step (so without the shifts)
		// instead, we change to DISTO_FILE_COMET and use the refimage as distortion master (if it has distortion)
			if (ref.keywords.wcslib && ref.keywords.wcslib->lin.dispre != NULL) {
				char buffer[256];
				fit_sequence_get_image_filename_checkext(args->seq, args->seq->reference_image, buffer);
				cadata->distoparam[regargs->layer].filename = g_strdup(buffer);
			}
			if (image_is_flipped_from_wcs(ref.keywords.wcslib)) // and if astrometry is flipped
				cadata->flipped = TRUE;
		}
	}
	if (cadata->flipped)
		regargs->velocity.y *= -1; // we correct the velocity due to the filp from astrometry which will not be applied during applyreg

	if (layer_has_registration(args->seq, regargs->layer)) { // we must keep track to correctly recompute astrometry during applyreg
		if (!cadata->distoparam) // there was no distorsion in any layer
			cadata->distoparam = calloc(args->seq->nb_layers, sizeof(disto_params));
		cadata->distoparam[regargs->layer].index = DISTO_FILE_COMET;
		cadata->distoparam[regargs->layer].velocity = regargs->velocity;
		if (cadata->flipped)
			cadata->distoparam[regargs->layer].velocity.y *= -1;
	}

	clearfits(&ref);
	return 0;
}

static int comet_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;

	pointf reg = { 0.f, 0.f };
	get_comet_shift(regargs->reference_date, fit->keywords.date_obs, regargs->velocity, &reg);
	if (cadata->current_regdata) {
		if (guess_transform_from_H(cadata->current_regdata[in_index].H) == NULL_TRANSFORMATION) {
			siril_log_color_message(_("Image %d has no registration data, removing\n"), "red", in_index + 1);
			regargs->imgparam[in_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
			return 1;
		}
		regargs->regparam[in_index].H = cadata->current_regdata[in_index].H; // previous reg exists, we set H with the initial H matrix
	} else {
		cvGetEye(&regargs->regparam[in_index].H); // previous reg does not exist, we set H with identity
	}
	cum_shifts(regargs->regparam, in_index, -reg.x, -reg.y); // we left-compose with the additional shift

	regargs->imgparam[in_index].incl = SEQUENCE_DEFAULT_INCLUDE;

	if (in_index == regargs->reference_image)
		new_ref_index = in_index; // keeping track of the new ref index in output sequence

	return 0;
}

static int comet_save_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;
	if (args->seq->type == SEQ_REGULAR) {
		char src[PATH_MAX];
		seq_get_image_filename(args->seq, in_index, src);
		gchar *dest = fit_sequence_get_image_filename_prefixed(args->seq, regargs->prefix, in_index);
		int res = symlink_uniq_file(src, dest, TRUE);
		g_free(dest);
		return res;
	}
	return 0;
}

static int comet_align_finalize_hook(struct generic_seq_args *args) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;

	regargs->new_total = args->nb_filtered_images;

 	// for fitseq or ser, we copy a symlink to the whole file (or hard-copy if symlink not possible)
	if (!args->retval && (args->seq->type == SEQ_FITSEQ || args->seq->type == SEQ_SER)) {
		gchar *basename = (args->seq->type == SEQ_FITSEQ) ? args->seq->fitseq_file->filename : args->seq->ser_file->filename;
		GString *prefixedname = g_string_new(basename);
		prefixedname = g_string_prepend(prefixedname, regargs->prefix);
		gchar *outname = g_string_free(prefixedname, FALSE);
		args->retval = symlink_uniq_file(basename, outname, TRUE);
		g_free(outname);
	}

	if (!args->retval) {
		siril_log_message(_("Applying registration completed.\n"));
		// explicit sequence creation to copy imgparam and regparam
		create_output_sequence_for_comet(regargs, new_ref_index, cadata->distoparam);
		// will be loaded in the idle function if (load_new_sequence)
		regargs->load_new_sequence = TRUE;
	}

	free(cadata);
	args->user = NULL;
	new_ref_index = -1;
	return 0;
}

int register_comet(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	/* we don't need to read image data, for simplicity we just read one
	 * pixel from it, making sure the header is read */
	args->partial_image = TRUE;
	args->area.x = 0; args->area.y = 0;
	args->area.w = 1; args->area.h = 1;
	args->layer_for_partial = 0;
	args->get_photometry_data_for_partial = TRUE;

	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->prepare_hook = comet_align_prepare_hook;
	args->image_hook = comet_align_image_hook;
	args->finalize_hook = comet_align_finalize_hook;
	if (args->seq->type == SEQ_REGULAR) {
		args->save_hook = comet_save_hook; // saves a symlink to original images
		args->has_output = TRUE;
		args->new_seq_prefix = strdup(regargs->prefix);
	}
	args->description = _("Moving object registration");
	args->already_in_a_thread = TRUE;
	args->stop_on_error = FALSE;

	struct comet_align_data *cadata = calloc(1, sizeof(struct comet_align_data));
	if (!cadata) {
		free_generic_seq_args(args, FALSE);
		return -1;
	}
	cadata->regargs = regargs;
	args->user = cadata;

	generic_sequence_worker(args);

	regargs->retval = args->retval;
	free(args);
	return regargs->retval;
}
