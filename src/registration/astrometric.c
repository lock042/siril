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
#include "algos/statistics.h"
#include "drizzle/cdrizzlebox.h"
#include "drizzle/cdrizzlemap.h"
#include "drizzle/cdrizzleutil.h"
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
	if (nb == 0)
		return;
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
	if (!regargs->filtering_criterion &&
			convert_parsed_filter_to_filter(&regargs->filters,
				regargs->seq, &regargs->filtering_criterion,
				&regargs->filtering_parameter))
		return 1;
	float scale = (regargs->driz) ? regargs->driz->scale : regargs->astrometric_scale;

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
		if (regargs->filtering_criterion && !regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
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
			convert_parsed_filter_to_filter(&regargs->filters,
				regargs->seq, &regargs->filtering_criterion,
				&regargs->filtering_parameter);
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
		if (!incl[i])
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
		if (!incl[i])
			continue;
		double angles[3];
		angles[0] = 90. - RA[i]; // TODO: need to understand why 90- instead of 90 + .... opencv vs wcs convention?
		angles[1] = 90. - DEC[i];
		// initializing K with center point (-0.5 is for opencv convention)
		cvGetEye(Ks + i); // initializing to unity
		Ks[i].h02 = 0.5 * ((regargs->seq->is_variable) ? regargs->seq->imgparam[i].rx : regargs->seq->rx) - 0.5;
		Ks[i].h12 = 0.5 * ((regargs->seq->is_variable) ? regargs->seq->imgparam[i].ry : regargs->seq->ry) - 0.5;
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
		if (!incl[i])
			continue;
		cvRelRot(&Rref, Rs + i);
	}

	// We compute the H matrices wrt to ref as Kref * Rrel * Kimg^-1
	// The K matrices have the focals on the diag and (rx/2, ry/2) for the translation terms
	// We add to the seq in order to display some kind of alignment, sufficient if small FOV
	for (int i = 0;  i < n; i++) {
		if (!incl[i])
			continue;
		Homography H = { 0 };
		cvcalcH_fromKKR(Ks[refindex], Ks[i], Rs[i], &H);
		if (!regargs->no_output)
			current_regdata[i].H = H;
		nb_aligned++;
	}
	regargs->seq->reference_image = refindex;

	// Give the full output size
	float fscale = 0.5 * (fabs(Ks[refindex].h00) + fabs(Ks[refindex].h11)) * scale;
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
		astargs->scale = fscale;
		astargs->disto = NULL;
		gboolean found = FALSE;
		if (regargs->undistort) {
			if (regargs->seq->is_variable) {
				// sequence has variable size, we need to use each image undistorsion coeffs
				astargs->disto = calloc(n, sizeof(disto_data));
				for (int i = 0;  i < n; i++) {
					if (!incl[i])
						continue;
					if (WCSDATA[i].lin.dispre) {
						found = TRUE;
						astargs->disto[i].order = extract_SIP_order_and_matrices(WCSDATA[i].lin.dispre, astargs->disto[i].A, astargs->disto[i].B, astargs->disto[i].AP, astargs->disto[i].BP);
						astargs->disto[i].xref = WCSDATA[i].crpix[0] - 1.; // -1 comes from the difference of convention between opencv and wcs
						astargs->disto[i].yref = WCSDATA[i].crpix[1] - 1.;
						astargs->disto[i].dtype = DISTO_D2S;
					} else {
						astargs->disto[i].dtype = DISTO_NONE;
					}
				}
			} else {
				// sequence has same size images, we search for the first image that has wcs data
				astargs->disto = calloc(1, sizeof(disto_data));
				for (int i = 0;  i < n; i++) {
					if (!incl[i])
						continue;
					if (WCSDATA[i].lin.dispre) { // we find the first image that has disto data and break afterwards
						astargs->disto[0].order = extract_SIP_order_and_matrices(WCSDATA[i].lin.dispre, astargs->disto[i].A, astargs->disto[i].B, astargs->disto[i].AP, astargs->disto[i].BP);
						astargs->disto[0].xref = WCSDATA[i].crpix[0] - 1.; // -1 comes from the difference of convention between opencv and wcs
						astargs->disto[0].yref = WCSDATA[i].crpix[1] - 1.;
						astargs->disto[0].dtype = (regargs->driz) ? DISTO_MAP_S2D: DISTO_MAP_D2S;
						found = TRUE;
						break;
					}
				}
				if (found) { // and we prepare the mapping that will be used for all images
					siril_log_message(_("Computing distortion mapping\n"));
					init_disto_map(regargs->seq->rx, regargs->seq->ry, astargs->disto);
					siril_log_message(_("Done\n"));
				}
			}
			if (!found) {
				siril_debug_print("no distorsion terms found in any of the images, disabling undistort\n");
				free(astargs->disto); // we don't need to call free_disto_args as maps have not been allocated
				astargs->disto = NULL;
				regargs->undistort = FALSE;
			}
		}
	}

	if (regargs->framing == FRAMING_COG || regargs->framing == FRAMING_CURRENT) {
		int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
		int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
		astrometric_roi roi_in = {.x = 0, .y = 0, .w = rx, .h = ry};
		Homography R = { 0 };
		cvGetEye(&R); // initializing to unity
		cvWarp_fromKR(NULL, &roi_in, Ks[refindex], R, fscale, OPENCV_NONE, FALSE, NULL, rois + refindex);
		tl.x = rois[refindex].x;
		tl.y = rois[refindex].y;
		br.x = rois[refindex].x + rois[refindex].w;
		br.y = rois[refindex].y + rois[refindex].h;
		siril_debug_print("ref:%d,%d,%d,%d,%d\n", refindex + 1, rois[refindex].x, rois[refindex].y, rois[refindex].w, rois[refindex].h);
	}

	gboolean error = FALSE;

	for (int i = 0;  i < n; i++) {
		if (!incl[i])
			continue;
		int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].rx : regargs->seq->rx;
		int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].ry : regargs->seq->ry;
		astrometric_roi roi_in = {.x = 0, .y = 0, .w = rx, .h = ry};
		cvWarp_fromKR(NULL, &roi_in, Ks[i], Rs[i], fscale, OPENCV_NONE, FALSE, NULL, rois + i);
		// first determine the corners
		siril_debug_print("%d,%d,%d,%d,%d\n", i + 1, rois[i].x, rois[i].y, rois[i].w, rois[i].h);
		switch (regargs->framing) {
			case FRAMING_MAX:
				if (rois[i].x < tl.x) tl.x = rois[i].x;
				if (rois[i].y < tl.y) tl.y = rois[i].y;
				if (rois[i].x + rois[i].w > br.x) br.x = rois[i].x + rois[i].w;
				if (rois[i].y + rois[i].h > br.y) br.y = rois[i].y + rois[i].h;
				break;
			case FRAMING_MIN:
				if (rois[i].x > tl.x) tl.x = rois[i].x;
				if (rois[i].y > tl.y) tl.y = rois[i].y;
				if (rois[i].x + rois[i].w < br.x) br.x = rois[i].x + rois[i].w;
				if (rois[i].y + rois[i].h < br.y) br.y = rois[i].y + rois[i].h;
				break;
			case FRAMING_COG:
			case FRAMING_CURRENT:
			default:
				// we just perform checks to make sure there is some overlap
				if (i != regargs->reference_image) {
					if (rois[i].x > br.x || rois[i].x + rois[i].w < tl.x ||
						rois[i].y > br.y || rois[i].y + rois[i].h < tl.y) {
							error = TRUE;
							siril_log_color_message(_("Image %d has no overlap with the reference\n"), "red", i + 1);
						}
				}
				break;
		}
	}
	// then compute the roi size and full image space requirement
	int imagew = br.x - tl.x;
	int imageh = br.y - tl.y;
	if  (regargs->framing == FRAMING_MIN && (imageh <= 0 || imagew <= 0)) {
		siril_log_color_message(_("The intersection of all images is null or negative\n"), "red");
		error = TRUE;
	}
	if (error) {
		siril_log_color_message(_("Unselect the images generating the error or change framing method to max\n"), "red");
		retval = 1;
		goto free_all;
	}

	gchar *downscale = (scale != 1.f) ? g_strdup_printf(_(" (assuming a scaling factor of %.2f)"), scale) : g_strdup("");
	siril_log_color_message(_("Output image: %d x %d pixels%s\n"), "salmon", imagew, imageh, downscale);
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
	for (int i = 0; i < n; i++) {
		if (!incl[i])
			continue;
		wcsfree(WCSDATA + i);
	}
	free(WCSDATA);
	free(incl);
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
	struct driz_args_t *driz = regargs->driz;
	float scale = 1.f;
	struct driz_param_t *p = NULL;
	int status = 0;
	Homography *Rs = astargs->Rs;
	Homography *Ks = astargs->Ks;
	astrometric_roi roi = { 0 };
	Homography H = { 0 };
	disto_data *disto = NULL;
	if (regargs->undistort && astargs->disto) {
		if (astargs->disto[0].dtype == DISTO_MAP_D2S || astargs->disto[0].dtype == DISTO_MAP_S2D) {
			disto = &astargs->disto[0];
		} else {
			disto = &astargs->disto[in_index];
		}
	}
	cvGetEye(&H);

	sadata->success[out_index] = 0;
	if (!regargs->driz) {
		// TODO: find in opencv codebase if smthg smart can be done with K/R to avoid the double-flip
		fits_flip_top_to_bottom(fit);
		astrometric_roi *roi_in = (regargs->framing != FRAMING_MAX) ?  &astargs->roi : NULL;
		status = cvWarp_fromKR(fit, roi_in, Ks[in_index], Rs[in_index], astargs->scale, regargs->interpolation, regargs->clamp, disto, &roi);
		if (!status) {
			fits_flip_top_to_bottom(fit);
			if (regargs->framing == FRAMING_MAX) {
				H.h02 = roi.x - astargs->roi.x;
				H.h12 = roi.y - astargs->roi.y;
			}
			free_wcs(fit); // we remove astrometric solution
		}
	} else {
		scale = driz->scale;
		p = calloc(1, sizeof(struct driz_param_t));
		driz_param_init(p);
		p->kernel = driz->kernel;
		p->driz = driz;
		p->error = malloc(sizeof(struct driz_error_t));
		p->scale = driz->scale;
		p->pixel_fraction = driz->pixel_fraction;
		p->cfa = driz->cfa;
		// Set bounds equal to whole image
		p->xmin = p->ymin = 0;
		p->xmax = fit->rx - 1;
		p->ymax = fit->ry - 1;
		p->pixmap = calloc(1, sizeof(imgmap_t));
		p->pixmap->rx = fit->rx;
		p->pixmap->ry = fit->ry;
		p->threads = threads;
		Homography Hs = { 0 }, Himg = { 0 }, Htransf = { 0 };
		cvGetEye(&Htransf);
		cvGetEye(&Hs);
		cvcalcH_fromKKR(Ks[regargs->seq->reference_image], Ks[in_index], Rs[in_index], &Himg);
		int dst_rx, dst_ry;
		if (regargs->framing != FRAMING_MAX) {
			cvTransfH(Himg, Htransf, &H);
			dst_rx = astargs->roi.w;
			dst_ry = astargs->roi.h;
		} else {
			compute_Hmax(&Himg, &Htransf, fit->rx, fit->ry, scale, &H, &Hs, &dst_rx, &dst_ry);
		}
		status = map_image_coordinates_h(fit, H, p->pixmap, dst_ry, driz->scale, disto, threads);
		H = Hs;

		if (status) {
			siril_log_color_message(_("Error generating mapping array.\n"), "red");
			// TODO: need to handle this, cannot go on
			free(p->error);
			free(p->pixmap);
			free(p);
		}
		p->data = fit;
		// Convert fit to 32-bit float if required
		float *newbuf = NULL;
		if (fit->type == DATA_USHORT) {
			siril_debug_print("Replacing ushort buffer for drizzling\n");
			size_t ndata = fit->rx * fit->ry * fit->naxes[2];
			newbuf = malloc(ndata * sizeof(float));
			float invnorm = 1.f / USHRT_MAX_SINGLE;
			for (size_t i = 0 ; i < ndata ; i++) {
				newbuf[i] = fit->data[i] * invnorm;
			}
			fit_replace_buffer(fit, newbuf, DATA_FLOAT);
		}
		/* Set up output fits */
		fits out = { 0 };
		copyfits(fit, &out, CP_FORMAT, -1);
		// copy the DATE_OBS
		out.keywords.date_obs = g_date_time_ref(fit->keywords.date_obs);
		out.rx = out.naxes[0] = dst_rx;
		out.ry = out.naxes[1] = dst_ry;
		out.naxes[2] = driz->is_bayer ? 3 : 1;
		size_t chansize = out.rx * out.ry * sizeof(float);
		out.fdata = calloc(out.naxes[2] * chansize, 1);
		out.fpdata[0] = out.fdata;
		out.fpdata[1] = out.naxes[2] == 1 ? out.fdata : out.fdata + chansize;
		out.fpdata[2] = out.naxes[2] == 1 ? out.fdata : out.fdata + 2 * chansize;
		p->output_data = &out;

		// Set up the output_counts fits to store pixel hit counts
		fits *output_counts = calloc(1, sizeof(fits));
		copyfits(&out, output_counts, CP_FORMAT, -1);
		output_counts->fdata = calloc(output_counts->rx * output_counts->ry * output_counts->naxes[2], sizeof(float));
		p->output_counts = output_counts;

		p->weights = driz->flat;

		if ((status = dobox(p))) { // Do the drizzle
			siril_log_color_message("s\n", p->error->last_message);
			// TODO: need to handle this, cannot go on
		}
		clearfits(fit);
		// copy the DATE_OBS
		copyfits(&out, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		fit->keywords.date_obs = g_date_time_ref(out.keywords.date_obs);
		clearfits(&out);
		if (args->seq->type == SEQ_SER || com.pref.force_16bit) {
			fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, fit->rx * fit->ry * fit->naxes[2]), DATA_USHORT);
		}
		if (driz->is_bayer) {
			/* we need to do something special here because it's a 1-channel sequence and
			* image that will become 3-channel after this call, so stats caching will
			* not be correct. Preprocessing does not compute new stats so we don't need
			* to save them into cache now.
			* Destroying the stats from the fit will also prevent saving them in the
			* sequence, which may be done automatically by the caller.
			*/
			full_stats_invalidation_from_fit(fit);
			fit->keywords.lo = 0;
		}

		free(p->pixmap->xmap);
		free(p->pixmap);

		// Get rid of the output_counts image, no longer required
		clearfits(output_counts);
		free(output_counts);
		output_counts = NULL;
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
	args->filtering_criterion = regargs->filtering_criterion;
	args->filtering_parameter = regargs->filtering_parameter;
	args->nb_filtered_images = compute_nb_filtered_images(regargs->seq,
			regargs->filtering_criterion, regargs->filtering_parameter);

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

static void free_disto_args(disto_data *disto) {
	if (!disto)
		return;
	// we only need to free the maps for the 2 types which store them (disto has only one element in that case)
	if (disto->dtype == DISTO_MAP_D2S || disto->dtype == DISTO_MAP_S2D) {
		free(disto->xmap);
		free(disto->ymap);
	}
}

void free_astrometric_args(struct astrometric_args *astargs) {
	if (!astargs)
		return;
	free(astargs->Ks);
	free(astargs->Rs);
	free_disto_args(astargs->disto);
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
