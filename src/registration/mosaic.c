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

static int mosaic_alignment(struct registration_args *regargs, struct mosaic_args *mosargs);

// computing center cog dealing with range jumps
// https://stackoverflow.com/questions/6671183/calculate-the-center-point-of-multiple-latitude-longitude-coordinate-pairs
// ra is lon, dec is lat
static void compute_center_cog(double *ra, double *dec, int n, double *RA, double *DEC) {
	double x = 0.;
	double y = 0.;
	double z = 0.;
	for (int i = 0; i < n; i++) {
		double ra_rad = ra[i] * DEGTORAD;
		double dec_rad = dec[i] * DEGTORAD;
		x += cos(dec_rad) * cos(ra_rad);
		y += cos(dec_rad) * sin(ra_rad);
		z += sin(dec_rad);
	}
	x /= (double)n;
	y /= (double)n;
	z /= (double)n;
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

	K->h00 = -RADTODEG / cdelt[0];
	K->h11 = RADTODEG / cdelt[1];
	double rotation1 = atan2(pc[0][1], pc[0][0]);
	double rotation2 = atan2(-pc[1][0], pc[1][1]);
	*framing = -0.5 * (rotation1 + rotation2) * RADTODEG;
	return TRUE;
}

int register_mosaic(struct registration_args *regargs) {

	int retval = 0;
	struct timeval t_start, t_end;

	regdata *current_regdata = star_align_get_current_regdata(regargs); // clean the structure if it exists, allocates otherwise
	if (!current_regdata)
		return -2;
	regargs->type = HOMOGRAPHY_TRANSFORMATION; // forcing homography calc

	gettimeofday(&t_start, NULL);
	fits fit = { 0 };
	double ra0 = 0., dec0 = 0.;
	int n = regargs->seq->number;
	double *RA = malloc(n * sizeof(double));
	double *DEC = malloc(n * sizeof(double));
	double *dist = malloc(n * sizeof(double));
	struct wcsprm *WCSDATA = malloc(n * sizeof(struct wcsprm));
	mosaic_roi *rois = malloc(n * sizeof(mosaic_roi));
	int failed = 0, nb_aligned = 0;

	if (!RA || !DEC || !dist || !WCSDATA || !rois) {
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
			siril_log_message(_("Image %d has not been plate-solved, aborting\n"), i + 1);
			retval = 1;
			goto free_all;
		}
		center2wcs(&fit, RA + i, DEC + i);
		WCSDATA[i].flag = -1;
		wcssub(1, fit.wcslib, NULL, NULL, WCSDATA + i); // copying wcsprm structure for each fit to avoid reopening
		clearfits(&fit);
		siril_log_message("Image #%2d - RA:%.3f - DEC:%.3f\n", i + 1, RA[i], DEC[i]);
	}

	// find sequence cog and closest image to use as reference
	compute_center_cog(RA, DEC, n, &ra0, &dec0);
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

	// TODO: should check allocation success
	Homography *Rs = calloc(n, sizeof(Homography)); // camera rotation matrices
	Homography *Ks = calloc(n, sizeof(Homography)); // camera intrinsic matrices

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
	for (int i = 0; i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		double angles[3];
		angles[0] = 90. - RA[i]; 
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
	}

	// computing relative rotations wrt to ref image
	Homography Rref = { 0 };
	memcpy(&Rref, Rs + refindex, sizeof(Homography));
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
		if (i != refindex) {
			cvcalcH_fromKKR(Ks[refindex], Ks[i], Rs[i], &H);
		} else {
			cvGetEye(&H);
		}
		current_regdata[i].roundness = 1.;
		current_regdata[i].fwhm = 1.;
		current_regdata[i].weighted_fwhm = 1.;
		current_regdata[i].background_lvl = 0.;
		current_regdata[i].number_of_stars = 1;
		current_regdata[i].H = H;
		nb_aligned++;
	}
	regargs->seq->reference_image = refindex;

	// Give the full mosaic output size
	float scale = 0.5 * (fabs(Ks[refindex].h00) + fabs(Ks[refindex].h11));
	pointi tl = { INT_MAX, INT_MAX }, br = { INT_MIN, INT_MIN }; // top left and bottom-right
	gboolean savewarped = !regargs->no_output;
	struct mosaic_args *margs = NULL;
	// if we have some output sequence, prepare the mosaic_args
	if (savewarped) {
		margs = malloc(sizeof(struct mosaic_args));
		margs->nb = n;
		margs->refindex = refindex;
		margs->Ks = Ks;
		margs->Rs = Rs;
		margs->scale = scale;
	}
	for (int i = 0;  i < n; i++) {
		if (!regargs->seq->imgparam[i].incl && regargs->filters.filter_included)
			continue;
		seq_read_frame_metadata(regargs->seq, i, &fit);
		cvWarp_fromKR(&fit, Ks[i], Rs[i], scale, rois + i, OPENCV_NONE, FALSE);
		clearfits(&fit);
		// first determine the corners
		if (rois[i].x < tl.x) tl.x = rois[i].x;
		if (rois[i].y < tl.y) tl.y = rois[i].y;
		if (rois[i].x + rois[i].w > br.x) br.x = rois[i].x + rois[i].w;
		if (rois[i].y + rois[i].h > br.y) br.y = rois[i].y + rois[i].h;
	}
	// then compute the roi size and full mosaic space requirement
	int mosaicw = br.x - tl.x;
	int mosaich = br.y - tl.y;
	siril_log_message(_("Full size output mosaic: %d x %d pixels (assuming spherical projection)\n"), mosaicw, mosaich);
	int64_t frame_size = mosaicw * mosaich * regargs->seq->nb_layers;
	frame_size *= (get_data_type(regargs->seq->bitpix) == DATA_USHORT) ? sizeof(WORD) : sizeof(float);
	gchar *mem = g_format_size_full(frame_size, G_FORMAT_SIZE_IEC_UNITS);
	siril_log_message(_("Space required for full size storage: %s\n"), mem);
	g_free(mem);
	if (savewarped) {
		margs->tl = tl;
	}

free_all:
	free(RA);
	free(DEC);
	free(dist);
	free(rois);
	for (int i = 0; i < n; i++) {
		wcsfree(WCSDATA + i);
	}
	free(WCSDATA);
	if (!retval && !regargs->no_output) {
		return mosaic_alignment(regargs, margs);
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
static int mosaic_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct mosaic_args *mosargs = sadata->mosargs;
	Homography *Rs = mosargs->Rs;
	Homography *Ks = mosargs->Ks;
	mosaic_roi roi = { 0 };
	Homography H = { 0 };
	cvGetEye(&H);

	sadata->success[out_index] = 0;
	// TODO: find in opencv codebase if smthg smart can be done with K/R to avoid the double-flip
	fits_flip_top_to_bottom(fit);
	int status = cvWarp_fromKR(fit, Ks[out_index], Rs[out_index], mosargs->scale, &roi, regargs->interpolation, regargs->clamp);
	if (!status) {
		fits_flip_top_to_bottom(fit);
		H.h02 = roi.x - mosargs->tl.x;
		H.h12 = roi.y - mosargs->tl.y;
	}
	regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	regargs->imgparam[out_index].incl = (int)(!status);
	regargs->imgparam[out_index].rx = roi.w;
	regargs->imgparam[out_index].ry = roi.h;
	regargs->regparam[out_index].fwhm = (in_index == mosargs->refindex) ? 0.5 : 1;
	regargs->regparam[out_index].weighted_fwhm = 1.;
	regargs->regparam[out_index].roundness = 1.;
	regargs->regparam[out_index].background_lvl = 0.;
	regargs->regparam[out_index].number_of_stars = 1;
	regargs->regparam[out_index].H = H;
	sadata->success[out_index] = (int)(!status);
	return status;
}


static int mosaic_alignment(struct registration_args *regargs, struct mosaic_args *mosargs) {
	struct generic_seq_args *args = create_default_seqargs(&com.seq);
	if (regargs->filters.filter_included) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->prepare_hook = star_align_prepare_results; // from global registration
	args->image_hook = mosaic_image_hook;
	args->finalize_hook = star_align_finalize_hook;	// from global registration
	args->stop_on_error = FALSE;
	args->description = _("Creating the warped images sequence");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = g_strdup("mosaic_");
	args->load_new_sequence = TRUE;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		free(args);
		return -1;
	}
	g_free(regargs->prefix);
	regargs->prefix = g_strdup("mosaic_");
	sadata->regargs = regargs;
	sadata->mosargs = mosargs;
	sadata->current_regdata = NULL;
	sadata->fitted_stars = 0;
	sadata->refstars = NULL;
	args->user = sadata;

	generic_sequence_worker(args);
	regargs->retval = args->retval;
	free(args);
	return regargs->retval;
}

void free_mosaic_args(struct mosaic_args *mosargs) {
	if (!mosargs)
		return;
	free(mosargs->Ks);
	free(mosargs->Rs);
	free(mosargs);
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