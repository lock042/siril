/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_world_cs.h"
#include "core/siril_log.h"
#include "algos/astrometry_solver.h"
#include "algos/siril_wcs.h"
#include "algos/sorting.h"
#include "io/sequence.h"
#include "io/siril_catalogues.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/degtorad.h"

static int astrometric_alignment(struct registration_args *regargs, struct astrometric_args *astargs, regdata *current_regdata);

// computing center cog dealing with range jumps
// https://stackoverflow.com/questions/6671183/calculate-the-center-point-of-multiple-latitude-longitude-coordinate-pairs
// ra is lon, dec is lat
static void compute_center_cog(double *ra, double *dec, int n, gboolean *incl, double *RA, double *DEC) {
	double x = 0.;
	double y = 0.;
	double z = 0.;
	int nb = 0;
	for (int i = 0; i < n; i++) {
		if (!incl[i])
			continue;
		double ra_rad = ra[i] * DEGTORAD;
		double dec_rad = dec[i] * DEGTORAD;
		x += cos(dec_rad) * cos(ra_rad);
		y += cos(dec_rad) * sin(ra_rad);
		z += sin(dec_rad);
		nb++;
	}
	x /= (double)nb;
	y /= (double)nb;
	z /= (double)nb;
	double norm = sqrt( x * x + y * y + z * z);
	x /= norm;
	y /= norm;
	z /= norm;
	*RA = atan2(y, x) * RADTODEG;
	*DEC = atan2(z, sqrt(x * x + y * y)) * RADTODEG;
}

static gboolean get_scales_and_framing(struct wcsprm *WCSDATA, Homography *K, double *framing) {
	K->h00 = -1;
	K->h11 = -1;
	*framing = -1;
	double pc[2][2];
	double cdelt[2];
	double *pcij = WCSDATA->pc;
	for (int i = 0; i < 2; i++) {
		cdelt[i] = WCSDATA->cdelt[i];
		for (int j = 0; j <2; j++) {
			pc[i][j] = *(pcij++);
		}
	}
	int flipped = (image_is_flipped_from_wcs(WCSDATA)) ? -1 : 1;

	K->h00 = -RADTODEG / cdelt[0];
	K->h11 = flipped * RADTODEG / cdelt[1];
	double rotation1 = atan2(flipped * pc[0][1], pc[0][0]);
	double rotation2 = atan2(-pc[1][0], flipped * pc[1][1]);
	*framing = -0.5 * (rotation1 + rotation2) * RADTODEG;
	return TRUE;
}

int register_astrometric(struct registration_args *regargs) {

	int retval = 0;
	struct timeval t_start, t_end;

	regdata *current_regdata = apply_reg_get_current_regdata(regargs); // clean the structure if it exists, allocates otherwise
	if (!current_regdata)
		return -2;
	regargs->type = HOMOGRAPHY_TRANSFORMATION; // forcing homography calc for the input sequence

	gettimeofday(&t_start, NULL);
	fits fit = { 0 };
	double ra0 = 0., dec0 = 0.;
	int n = regargs->seq->number;
	double *RA = malloc(n * sizeof(double));
	double *DEC = malloc(n * sizeof(double));
	double *dist = malloc(n * sizeof(double));
	struct wcsprm *WCSDATA = calloc(n, sizeof(struct wcsprm));
	astrometric_roi *rois = malloc(n * sizeof(astrometric_roi));
	int failed = 0, nb_aligned = 0;
	Homography *Rs = NULL, *Ks = NULL;
	Rs = calloc(n, sizeof(Homography)); // camera rotation matrices
	Ks = calloc(n, sizeof(Homography)); // camera intrinsic matrices
	gboolean *incl = calloc(n, sizeof(gboolean));

	if (!RA || !DEC || !dist || !WCSDATA || !rois || !Rs || !Ks || !incl) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}

	for (int i = 0; i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		if (seq_read_frame_metadata(regargs->seq, i, &fit)) {
			siril_log_message(_("Could not load image %d from sequence %s\n"),
			i + 1, regargs->seq->seqname);
			retval = 1;
			goto free_all;
		}
		if (!has_wcs(&fit)) {
			siril_log_message(_("Image %d has not been plate-solved, unselecting\n"), i + 1);
			regargs->seq->imgparam[i].incl = FALSE;
			regargs->filters.filter_included = TRUE;
			failed++;
			continue;
		}
		center2wcs(&fit, RA + i, DEC + i);
		WCSDATA[i].flag = -1;
		wcssub(1, fit.keywords.wcslib, NULL, NULL, WCSDATA + i); // copying wcsprm structure for each fit to avoid reopening
		clearfits(&fit);
		siril_log_message("Image #%2d - RA:%.3f - DEC:%.3f\n", i + 1, RA[i], DEC[i]);
		incl[i] = TRUE;
	}

	// find sequence cog and closest image to use as reference
	compute_center_cog(RA, DEC, n, incl, &ra0, &dec0);
	ra0 = (ra0 < 0) ? ra0 + 360. : ra0;
	siril_log_message("Sequence COG - RA:%.3f - DEC:%.3f\n", ra0, dec0);
	int refindex = -1;
	double mindist = DBL_MAX;
	for (int i = 0; i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		dist[i] = compute_coords_distance(RA[i], DEC[i], ra0, dec0);
		if (dist[i] < mindist) {
			mindist = dist[i];
			refindex = i;
		}
	}
	siril_debug_print("Closest image  is #%d - RA:%.3f - DEC:%.3f\n", refindex + 1, RA[refindex], DEC[refindex]);

	// Obtaining Camera extrinsic and instrinsic matrices (resp. R and K)
	// ##################################################################
	// We will define the rotations made by the camera wrt World:
	// - We start with the camera pointing towards North Pole, X (sensor width) aligned with RA = 0deg
	// - We first make a rotation around cam Z axis of 90 + RA degrees
	// - We then make a rotation around cam X axis of 90 - DEC degrees (0 to 180 from N to S pole)
	// - We finally rotate the camera around its Z axis to get the correct framing (called "framing" angle below)
	// We repeat this calc for each image
	// We can then compute the relative rotation matrix between the images and refimage
	// This is made by computing Rref.t()*Rimg, this is the extrinsic matrix for a pure rotation wrt. mosaic ref
	// The intrisic camera matrix is obtained thanks to the CDij or PCij+CDELTi WCS values:
	// The focal length is 180 / pi / average of abs(CDELTi) - Note: CDELTi * 3600 = image sampling in "/px
	// The "framing" angle is obtained by making the PC matrix a true rotation matrix (averaging of CROTAi)

	rotation_type rottypes[3] = { ROTZ, ROTX, ROTZ};
	double framingref = 0.;
	int framingref_index = (regargs->framing == FRAMING_CURRENT) ? regargs->seq->reference_image : refindex;
	for (int i = 0; i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		double angles[3];
		angles[0] = 90. - RA[i]; // TODO: need to understand why 90- instead of 90 + .... opencv vs wcs convention?
		angles[1] = 90. - DEC[i];
		// initializing K with center point
		cvGetEye(Ks + i); // initializing to unity
		Ks[i].h02 = 0.5 * ((regargs->seq->is_variable) ? regargs->seq->imgparam[i].rx : regargs->seq->rx);
		Ks[i].h12 = 0.5 * ((regargs->seq->is_variable) ? regargs->seq->imgparam[i].ry : regargs->seq->ry);
		if (!get_scales_and_framing(WCSDATA + i, Ks + i, angles + 2)) {
			siril_log_message(_("Could not compute camera parameters of image %d from sequence %s\n"),
			i + 1, regargs->seq->seqname);
			retval = 1;
			goto free_all;
		}
		cvRotMat3(angles, rottypes, TRUE, Rs + i);
		siril_debug_print("Image #%d - rot:%.3f\n", i + 1, angles[2]);
		if (i == framingref_index)
			framingref = angles[2];
	}

	// computing relative rotations wrt to proj center + central image framing
	// this is important for warping by the spherical projector
	Homography Rref = { 0 };
	double angles[3];
	switch (regargs->framing) {
		case FRAMING_CURRENT:
		default:
			angles[0] = 90. - RA[regargs->seq->reference_image];
			angles[1] = 90. - DEC[regargs->seq->reference_image];
			break;
		case FRAMING_COG:
		case FRAMING_MAX:
		case FRAMING_MIN:
			angles[0] = 90. - ra0;
			angles[1] = 90. - dec0;
			break;
	}
	angles[2] = framingref;
	cvRotMat3(angles, rottypes, TRUE, &Rref);
	for (int i = 0; i < n; i++) {
		cvRelRot(&Rref, Rs + i);
	}

	// We compute the H matrices wrt to ref as Kref * Rrel * Kimg^-1
	// The K matrices have the focals on the diag and (rx/2, ry/2) for the translation terms
	// We add to the seq in order to display some kind of alignment, sufficient if small FOV
	for (int i = 0;  i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		Homography H = { 0 };
		cvcalcH_fromKKR(Ks[refindex], Ks[i], Rs[i], &H);
		if (!regargs->no_output)
			current_regdata[i].H = H;
		nb_aligned++;
	}
	regargs->seq->reference_image = refindex;

	// Give the full output size
	float scale = 0.5 * (fabs(Ks[refindex].h00) + fabs(Ks[refindex].h11)) * regargs->astrometric_scale;
	pointi tl, br; // top left and bottom-right
	if (regargs->framing == FRAMING_MIN) {
		tl = (pointi){ INT_MIN, INT_MIN };
		br = (pointi){ INT_MAX, INT_MAX };
	} else {
		tl = (pointi){ INT_MAX, INT_MAX };
		br = (pointi){ INT_MIN, INT_MIN };
	}
	gboolean savewarped = !regargs->no_output;
	struct astrometric_args *astargs = NULL;
	// if we have some output sequence, prepare the astrometric_args
	if (savewarped) {
		astargs = malloc(sizeof(struct astrometric_args));
		astargs->nb = n;
		astargs->refindex = refindex;
		astargs->Ks = Ks;
		astargs->Rs = Rs;
		astargs->scale = scale;
		gboolean found = FALSE;
		if (regargs->undistort) {
			astargs->disto = calloc(n, sizeof(disto_data));
			for (int i = 0;  i < n; i++) {
				if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
					continue;
				if (WCSDATA[i].lin.dispre) {
					double A[MAX_SIP_SIZE][MAX_SIP_SIZE], B[MAX_SIP_SIZE][MAX_SIP_SIZE]; // we won't need them
					astargs->disto[i].order = extract_SIP_order_and_matrices(WCSDATA[i].lin.dispre, A, B, astargs->disto[i].AP, astargs->disto[i].BP);
					found = TRUE;
					astargs->disto[i].xref = WCSDATA[i].crpix[0] - 1.; // -1 comes from the difference of convention between opencv and wcs
					astargs->disto[i].yref = WCSDATA[i].crpix[1] - 1.;
				}
			}
			if (!found) {
				siril_debug_print("no distorsion terms found in any of the images, disabling undistort\n");
				free(astargs->disto);
				astargs->disto = NULL;
				regargs->undistort = FALSE;
			}

		}
	}
	for (int i = 0;  i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
		int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
		astrometric_roi roi_in = {.x = 0, .y = 0, .w = rx, .h = ry};
		cvWarp_fromKR(NULL, &roi_in, Ks[i], Rs[i], scale, regargs->projector, OPENCV_NONE, FALSE, NULL, rois + i);
		// first determine the corners
		printf("%d,%d,%d,%d,%d\n", i + 1, rois[i].x, rois[i].y, rois[i].w, rois[i].h);
		switch (regargs->framing) {
			case FRAMING_MAX:
				if (rois[i].x < tl.x) tl.x = rois[i].x;
				if (rois[i].y < tl.y) tl.y = rois[i].y;
				if (rois[i].x + rois[i].w > br.x) br.x = rois[i].x + rois[i].w;
				if (rois[i].y + rois[i].h > br.y) br.y = rois[i].y + rois[i].h;
				break;
			case FRAMING_CURRENT:
			default:
				if (i == regargs->reference_image) {
					tl.x = rois[i].x;
					tl.y = rois[i].y;
					br.x = rois[i].x + rois[i].w;
					br.y = rois[i].y + rois[i].h;
				}
				break;
			case FRAMING_MIN:
				if (rois[i].x > tl.x) tl.x = rois[i].x;
				if (rois[i].y > tl.y) tl.y = rois[i].y;
				if (rois[i].x + rois[i].w < br.x) br.x = rois[i].x + rois[i].w;
				if (rois[i].y + rois[i].h < br.y) br.y = rois[i].y + rois[i].h;
				break;
			case FRAMING_COG:
				break;
		}
	}
	if (regargs->framing == FRAMING_COG) {
		int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
		int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
		astrometric_roi roi_in = {.x = 0, .y = 0, .w = rx, .h = ry};
		Homography R = { 0 };
		cvGetEye(&R); // initializing to unity
		cvWarp_fromKR(NULL, &roi_in, Ks[refindex], R, scale, regargs->projector, OPENCV_NONE, FALSE, NULL, rois);
		tl.x = rois[0].x;
		tl.y = rois[0].y;
		br.x = rois[0].x + rois[0].w;
		br.y = rois[0].y + rois[0].h;
	}
	// then compute the roi size and full image space requirement
	int imagew = br.x - tl.x;
	int imageh = br.y - tl.y;
	const gchar *proj = (regargs->projector == OPENCV_PLANE) ? _("plane") : _("spherical");
	gchar *downscale = (regargs->astrometric_scale != 1.f) ? g_strdup_printf(_(" and a scaling factor of %.2f"), regargs->astrometric_scale) : g_strdup("");
	siril_log_color_message(_("Output image: %d x %d pixels (assuming %s projection%s)\n"), "salmon", imagew, imageh, proj, downscale);
	g_free(downscale);
	int64_t frame_size = (int64_t)imagew * imageh * regargs->seq->nb_layers;
	frame_size *= (get_data_type(regargs->seq->bitpix) == DATA_USHORT) ? sizeof(WORD) : sizeof(float);
	gchar *mem = g_format_size_full(frame_size, G_FORMAT_SIZE_IEC_UNITS);
	siril_log_color_message(_("Space required for storage: %s\n"), "salmon", mem);
	g_free(mem);
	if (savewarped) {
		astargs->roi = (astrometric_roi) {.x = tl.x, .y = tl.y, .w = imagew, .h = imageh};
	}

free_all:
	free(RA);
	free(DEC);
	free(dist);
	free(rois);
	free(incl);
	for (int i = 0; i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		wcsfree(WCSDATA + i);
	}
	free(WCSDATA);
	if (!retval && !regargs->no_output) {
		return astrometric_alignment(regargs, astargs, current_regdata);
	}
	free(Ks);
	free(Rs);
	fix_selnum(regargs->seq, FALSE);
	siril_log_message(_("Registration finished.\n"));
	siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, nb_aligned);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	return retval;
}

/* image alignment hooks and main process */
static int astrometric_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct astrometric_args *astargs = sadata->astargs;
	Homography *Rs = astargs->Rs;
	Homography *Ks = astargs->Ks;
	astrometric_roi roi = { 0 };
	Homography H = { 0 };
	disto_data *disto = NULL;
	if (regargs->undistort && astargs->disto)
		disto = &astargs->disto[in_index];
	cvGetEye(&H);

	sadata->success[out_index] = 0;
	// TODO: find in opencv codebase if smthg smart can be done with K/R to avoid the double-flip
	fits_flip_top_to_bottom(fit);
	int status = cvWarp_fromKR(fit,  &astargs->roi, Ks[in_index], Rs[in_index], astargs->scale ,regargs->projector, regargs->interpolation, regargs->clamp, disto, &roi);
	if (!status) {
		fits_flip_top_to_bottom(fit);
		// TODO: keeping this for later when max framing won't save large black borders...
		// H.h02 = roi.x - astargs->roi.x;
		// H.h12 = roi.y - astargs->roi.y;
		free_wcs(fit); // we remove astrometric solution
	}
	regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	regargs->imgparam[out_index].incl = (int)(!status);
	regargs->imgparam[out_index].rx = roi.w;
	regargs->imgparam[out_index].ry = roi.h;
	regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;
	regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
	regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;
	regargs->regparam[out_index].background_lvl = sadata->current_regdata[in_index].background_lvl;
	regargs->regparam[out_index].number_of_stars = sadata->current_regdata[in_index].number_of_stars;
	regargs->regparam[out_index].H = H;
	sadata->success[out_index] = (int)(!status);
	if (astargs->scale != 1.f) {
		fit->keywords.pixel_size_x /= regargs->astrometric_scale;
		fit->keywords.pixel_size_y /= regargs->astrometric_scale;
		regargs->regparam[out_index].fwhm *= regargs->astrometric_scale;
		regargs->regparam[out_index].weighted_fwhm *= regargs->astrometric_scale;
	}
	return status;
}


static int astrometric_alignment(struct registration_args *regargs, struct astrometric_args *astargs, regdata *current_regdata) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	// args->compute_mem_limits_hook =  astrometric_mem_hook; // TODO
	args->prepare_hook = star_align_prepare_results; // from global registration
	args->image_hook = astrometric_image_hook;
	args->finalize_hook = star_align_finalize_hook;	// from global registration
	args->stop_on_error = FALSE;
	args->description = _("Creating the warped images sequence");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = g_strdup(regargs->prefix);
	args->load_new_sequence = TRUE;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		free(args);
		return -1;
	}
	sadata->regargs = regargs;
	sadata->astargs = astargs;
	sadata->current_regdata = current_regdata;
	sadata->fitted_stars = 0;
	sadata->refstars = NULL;
	args->user = sadata;

	generic_sequence_worker(args);
	regargs->retval = args->retval;
	free(args);
	return regargs->retval;
}

void free_astrometric_args(struct astrometric_args *astargs) {
	if (!astargs)
		return;
	free(astargs->Ks);
	free(astargs->Rs);
	free(astargs->disto);
	free(astargs);
}



// FOR LATER

	//Composing the mosaic

	// seq_read_frame_metadata(regargs->seq, refindex, &fit);
	// cvmosaiccompose(regargs->seq, Ks, Rs, n, scale, 0.2f, 1.f, &fit);
	// savefits("mosaic.fit", &fit);
	// clearfits(&fit);


// just in case

	// // now we store the relative Rs and cameras K
	// // writing outputs to *.smf file, a.k.a Siril Mosaic File
	// char *filename;
	// FILE *mscfile;
	// if (!regargs->seq->seqname || regargs->seq->seqname[0] == '\0') {
	// 	retval = 1;
	// 	goto free_all;
	// }
	// filename = malloc(strlen(regargs->seq->seqname)+5);
	// sprintf(filename, "%s.smf", regargs->seq->seqname);
	// mscfile = g_fopen(filename, "w+t"); // TODO deal with errors
	// for (int i = 0; i < n; i++) {
	// 	fprintf(mscfile, "%2d K %12.1f %12.1f %8.1f %8.1f R %g %g %g %g %g %g %g %g %g\n",
	// 	i + 1,
	// 	Ks[i].h00,
	// 	Ks[i].h11,
	// 	Ks[i].h02,
	// 	Ks[i].h12,
	// 	Rs[i].h00,
	// 	Rs[i].h01,
	// 	Rs[i].h02,
	// 	Rs[i].h10,
	// 	Rs[i].h11,
	// 	Rs[i].h12,
	// 	Rs[i].h20,
	// 	Rs[i].h21,
	// 	Rs[i].h22
	// 	);
	// }
	// fclose(mscfile);
