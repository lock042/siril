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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/siril_wcs.h"
#include "algos/PSF.h"
#include "io/fits_keywords.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "distorsion.h"

// applies a distortion correction to the xpos/ypos members of psf_star list
int disto_correct_stars(psf_star **stars, disto_data *disto) {
	if(!stars || !disto || disto->dtype == DISTO_NONE)
		return 1;
	double U, V, x, y;
	double U2, V2, U3, V3, U4, V4, U5, V5;
	int i = 0;
	while (stars && stars[i]) {
		U = stars[i]->xpos - disto->xref;
		V = disto->yref - stars[i]->ypos; // opencv convention is y down while wcs is y up
		// V = stars[i]->ypos - disto->yref; // opencv convention is y down while wcs is y up
		x = U + disto->A[0][0] + disto->A[1][0] * U + disto->A[0][1] * V;
		y = V + disto->B[0][0] + disto->B[1][0] * U + disto->B[0][1] * V;
		if (disto->order >= 2) {
			U2 = U * U;
			V2 = V * V;
			double UV = U * V;
			x += disto->A[2][0] * U2 + disto->A[1][1] * UV + disto->A[0][2] * V2;
			y += disto->B[2][0] * U2 + disto->B[1][1] * UV + disto->B[0][2] * V2;
			if (disto->order >= 3) {
				U3 = U2 * U;
				V3 = V2 * V;
				double U2V = U2 * V;
				double UV2 = U * V2;
				x += disto->A[3][0] * U3 + disto->A[2][1] * U2V + disto->A[1][2] * UV2 + disto->A[0][3] * V3;
				y += disto->B[3][0] * U3 + disto->B[2][1] * U2V + disto->B[1][2] * UV2 + disto->B[0][3] * V3;
				if (disto->order >= 4) {
					U4 = U3 * U;
					V4 = V3 * V;
					double U3V = U3 * V;
					double U2V2 = U2 * V2;
					double UV3 = U * V3;
					x += disto->A[4][0] * U4 + disto->A[3][1] * U3V + disto->A[2][2] * U2V2 + disto->A[1][3] * UV3 + disto->A[0][4] * V4;
					y += disto->B[4][0] * U4 + disto->B[3][1] * U3V + disto->B[2][2] * U2V2 + disto->B[1][3] * UV3 + disto->B[0][4] * V4;
					if (disto->order >= 5) {
						U5 = U4 * U;
						V5 = V4 * V;
						double U4V = U4 * V;
						double U3V2 = U3 * V2;
						double U2V3 = U2 * V3;
						double UV4 = U * V4;
						x += disto->A[5][0] * U5 + disto->A[4][1] * U4V + disto->A[3][2] * U3V2 + disto->A[2][3] * U2V3 + disto->A[1][4] * UV4 + disto->A[0][5] * V5;
						y += disto->B[5][0] * U5 + disto->B[4][1] * U4V + disto->B[3][2] * U3V2 + disto->B[2][3] * U2V3 + disto->B[1][4] * UV4 + disto->B[0][5] * V5;
					}
				}
			}

		}
		stars[i]->xpos = (x + disto->xref);
		stars[i]->ypos = (disto->yref - y);
		i++;
	}
	return 0;
}

// maps undistortion dst to src (for interpolation)
void map_undistortion_D2S(disto_data *disto, int rx, int ry, float *xmap, float *ymap) {
	g_assert(disto != NULL);
	g_assert(xmap != NULL);
	g_assert(ymap != NULL);
	double U, V, x, y;
	double U2, V2, U3, V3, U4, V4, U5, V5;
	int r = 0;
	for (int v = 0; v < ry; ++v) {
		float *rxptr = xmap + r;
		float *ryptr = ymap + r;
		for (int u = 0; u < rx; ++u) {
			U = (double)rxptr[u] - disto->xref;
			V = (double)ryptr[u] - disto->yref;
			x = U + disto->AP[0][0] + disto->AP[1][0] * U + disto->AP[0][1] * V;
			y = V + disto->BP[0][0] + disto->BP[1][0] * U + disto->BP[0][1] * V;
			if (disto->order >= 2) {
				U2 = U * U;
				V2 = V * V;
				double UV = U * V;
				x += disto->AP[2][0] * U2 + disto->AP[1][1] * UV + disto->AP[0][2] * V2;
				y += disto->BP[2][0] * U2 + disto->BP[1][1] * UV + disto->BP[0][2] * V2;
				if (disto->order >= 3) {
					U3 = U2 * U;
					V3 = V2 * V;
					double U2V = U2 * V;
					double UV2 = U * V2;
					x += disto->AP[3][0] * U3 + disto->AP[2][1] * U2V + disto->AP[1][2] * UV2 + disto->AP[0][3] * V3;
					y += disto->BP[3][0] * U3 + disto->BP[2][1] * U2V + disto->BP[1][2] * UV2 + disto->BP[0][3] * V3;
					if (disto->order >= 4) {
						U4 = U3 * U;
						V4 = V3 * V;
						double U3V = U3 * V;
						double U2V2 = U2 * V2;
						double UV3 = U * V3;
						x += disto->AP[4][0] * U4 + disto->AP[3][1] * U3V + disto->AP[2][2] * U2V2 + disto->AP[1][3] * UV3 + disto->AP[0][4] * V4;
						y += disto->BP[4][0] * U4 + disto->BP[3][1] * U3V + disto->BP[2][2] * U2V2 + disto->BP[1][3] * UV3 + disto->BP[0][4] * V4;
						if (disto->order >= 5) {
							U5 = U4 * U;
							V5 = V4 * V;
							double U4V = U4 * V;
							double U3V2 = U3 * V2;
							double U2V3 = U2 * V3;
							double UV4 = U * V4;
							x += disto->AP[5][0] * U5 + disto->AP[4][1] * U4V + disto->AP[3][2] * U3V2 + disto->AP[2][3] * U2V3 + disto->AP[1][4] * UV4 + disto->AP[0][5] * V5;
							y += disto->BP[5][0] * U5 + disto->BP[4][1] * U4V + disto->BP[3][2] * U3V2 + disto->BP[2][3] * U2V3 + disto->BP[1][4] * UV4 + disto->BP[0][5] * V5;
						}
					}
				}
			}
			rxptr[u] = (float)(x + disto->xref);
			ryptr[u] = (float)(y + disto->yref);
 		}
		r += rx;
 	}
}

// maps undistortion src to dst (for drizzle)
void map_undistortion_S2D(disto_data *disto, int rx, int ry, float *xmap, float *ymap) {
	g_assert(disto != NULL);
	g_assert(xmap != NULL);
	g_assert(ymap != NULL);
	double x, y;
	double *U = malloc(rx * sizeof(double));
	double *V = malloc(ry * sizeof(double));
	double *U2 = NULL, *V2 = NULL, *U3 = NULL, *V3 = NULL, *U4 = NULL, *V4 = NULL, *U5 = NULL, *V5 = NULL;
	if (disto->order >= 2) {
		U2 = malloc(rx * sizeof(double));
		V2 = malloc(ry * sizeof(double));
	}
	if (disto->order >= 3) {
		U3 = malloc(rx * sizeof(double));
		V3 = malloc(ry * sizeof(double));
	}
	if (disto->order >= 4) {
		U4 = malloc(rx * sizeof(double));
		V4 = malloc(ry * sizeof(double));
	}
	if (disto->order >= 5) {
		U5 = malloc(rx * sizeof(double));
		V5 = malloc(ry * sizeof(double));
	}
	for (int u = 0; u < rx; ++u) {
		U[u] = (double)u - disto->xref;
		if (disto->order >= 2) {
			U2[u] = U[u] * U[u];
			if (disto->order >= 3) {
				U3[u] = U2[u] * U[u];
				if (disto->order >= 4) {
					U4[u] = U3[u] * U[u];
					if (disto->order >= 5) {
						U5[u] = U4[u] * U[u];
					}
				}
			}
		}
	}

	for (int v = 0; v < ry; ++v) {
		V[v] = (double)v - disto->yref;
		if (disto->order >= 2) {
			V2[v] = V[v] * V[v];
			if (disto->order >= 3) {
				V3[v] = V2[v] * V[v];
				if (disto->order >= 4) {
					V4[v] = V3[v] * V[v];
					if (disto->order >= 5) {
						V5[v] = V4[v] * V[v];
					}
				}
			}
		}
	}

	int r = 0;
	for (int v = 0; v < ry; ++v) {
		for (int u = 0; u < rx; ++u) {
			x = U[u] + disto->A[0][0] + disto->A[1][0] * U[u] + disto->A[0][1] * V[v];
			y = V[v] + disto->B[0][0] + disto->B[1][0] * U[u] + disto->B[0][1] * V[v];
			if (disto->order >= 2) {
				double UV = U[u] * V[v];
				x += disto->A[2][0] * U2[u] + disto->A[1][1] * UV + disto->A[0][2] * V2[v];
				y += disto->B[2][0] * U2[u] + disto->B[1][1] * UV + disto->B[0][2] * V2[v];
				if (disto->order >= 3) {
					double U2V = U2[u] * V[v];
					double UV2 = U[u] * V2[v];
					x += disto->A[3][0] * U3[u] + disto->A[2][1] * U2V + disto->A[1][2] * UV2 + disto->A[0][3] * V3[v];
					y += disto->B[3][0] * U3[u] + disto->B[2][1] * U2V + disto->B[1][2] * UV2 + disto->B[0][3] * V3[v];
					if (disto->order >= 4) {
						double U3V = U3[u] * V[v];
						double U2V2 = U2[u] * V2[v];
						double UV3 = U[u] * V3[v];
						x += disto->A[4][0] * U4[u] + disto->A[3][1] * U3V + disto->A[2][2] * U2V2 + disto->A[1][3] * UV3 + disto->A[0][4] * V4[v];
						y += disto->B[4][0] * U4[u] + disto->B[3][1] * U3V + disto->B[2][2] * U2V2 + disto->B[1][3] * UV3 + disto->B[0][4] * V4[v];
						if (disto->order >= 5) {
							double U4V = U4[u] * V[v];
							double U3V2 = U3[u] * V2[v];
							double U2V3 = U2[u] * V3[v];
							double UV4 = U[u] * V4[v];
							x += disto->A[5][0] * U5[u] + disto->A[4][1] * U4V + disto->A[3][2] * U3V2 + disto->A[2][3] * U2V3 + disto->A[1][4] * UV4 + disto->A[0][5] * V5[v];
							y += disto->B[5][0] * U5[u] + disto->B[4][1] * U4V + disto->B[3][2] * U3V2 + disto->B[2][3] * U2V3 + disto->B[1][4] * UV4 + disto->B[0][5] * V5[v];
						}
					}
				}
			}
			xmap[r + u] = (float)(x + disto->xref);
			ymap[r + u] = (float)(y + disto->yref);
 		}
		r += rx;
 	}
	free(U);
	free(U2);
	free(U3);
	free(U4);
	free(U5);
	free(V);
	free(V2);
	free(V3);
	free(V4);
	free(V5);
}


// Computes the distortion map and stores it in the disto structure
int init_disto_map(int rx, int ry, disto_data *disto) {
	if (disto == NULL ||(disto->dtype != DISTO_MAP_D2S && disto->dtype != DISTO_MAP_S2D && disto->dtype != DISTO_S2D)) //nothing to do
		return 0;

	if (!disto->xmap) {
		disto->xmap = malloc(rx * ry *sizeof(float));
		disto->ymap = malloc(rx * ry *sizeof(float));
	}

	if (!disto->xmap || !disto->ymap) {
		free(disto->xmap);
		free(disto->ymap);
		disto->xmap = NULL;
		disto->ymap = NULL;
		return 2;
	}

	size_t s = 0;
	for (int j = 0; j < ry; j++) {
		for (int i = 0; i < rx; i++) {
			disto->xmap[s] = (float)i;
			disto->ymap[s] = (float)j;
			s++;
		}
	}

	if (disto->dtype == DISTO_MAP_D2S) {
		map_undistortion_D2S(disto, rx, ry, disto->xmap, disto->ymap);
	} else
		map_undistortion_S2D(disto, rx, ry, disto->xmap, disto->ymap);
	return 0;
}

// interpolates a disto map from the one saved in disto structure
static void map_undistortion_interp(disto_data *disto, int rx_in, int ry_in, int rx_out, int ry, float *xmap, float *ymap) {
	g_assert(disto != NULL);
	g_assert(xmap != NULL);
	g_assert(ymap != NULL);
	int r = 0;
	for (int v = 0; v < ry; ++v) {
		float *rxptr = xmap + r;
		float *ryptr = ymap + r;
		for (int u = 0; u < rx_out; ++u) {
			int i = floor(rxptr[u]);
			int j = floor(ryptr[u]);
			if (i < 0 || i > rx_in - 2 || j < 0 || j > ry_in - 2) {
				rxptr[u] = -1.f;
				ryptr[u] = -1.f;
			} else {
				int s = j * rx_in + i;
				float c1 = rxptr[u] - (float)i;
				float c2 = ryptr[u] - (float)j;
				float w11 = (1.f - c1) * (1.f - c2);
				float w12 = c1 * (1.f - c2);
				float w21 = (1.f - c1) * c2;
				float w22 = c1 * c2;
				rxptr[u] = w11 * disto->xmap[s] + w12 * disto->xmap[s + 1] + w21 * disto->xmap[s + rx_in] + w22 * disto->xmap[s + rx_in + 1];
				ryptr[u] = w11 * disto->ymap[s] + w12 * disto->ymap[s + 1] + w21 * disto->ymap[s + rx_in] + w22 * disto->ymap[s + rx_in + 1];
			}
		}
		r += rx_out;
	}
}

// computes the dst->src map from homography and corrects the mapping for distortions
void prepare_H_with_disto_4remap(double *H, int rx_in, int ry_in, int rx_out, int ry_out, disto_data *disto, float *xmap, float *ymap) {
	size_t index = 0;
	for (int y = 0; y < ry_out; y++) {
		for (int x = 0; x < rx_out; x++) {
			float x0 = (float)x;
			float y0 = (float)y;
			float z = 1. / (x0 * H[6] + y0 * H[7] + H[8]);
			xmap[index] = (x0 * H[0] + y0 * H[1] + H[2]) * z;
			ymap[index++] = (x0 * H[3] + y0 * H[4] + H[5]) * z;
		}
	}
	if (disto->dtype == DISTO_D2S) {
		map_undistortion_D2S(disto, rx_out, ry_out, xmap, ymap);
	} else if (disto->dtype == DISTO_MAP_D2S){
		map_undistortion_interp(disto, rx_in, ry_in, rx_out, ry_out, xmap, ymap);
	}
}

// get the master disto name if set
// otherwise, returns seqname.wcs
static gchar *get_wcs_filename(pathparse_mode mode, sequence *seq) {
	gchar *wcsname = NULL;
	gboolean found = FALSE;
	if (com.pref.prepro.disto_lib && com.pref.prepro.use_disto_lib) { //we have a distortion master and we should use it
		int status = 0; 
		wcsname = path_parse(&gfit, com.pref.prepro.disto_lib, mode, &status);
		if (status) {
			siril_log_color_message(_("Could not parse the distortion master, aborting\n"), "red");
			g_free(wcsname);
			wcsname = NULL;
		} else {
			found = TRUE;
		}
	}
	if (!found && seq) {
		char *namewoext = remove_ext_from_filename(seq->seqname);
		wcsname = g_strdup_printf("%s%s", namewoext, ".wcs");
		free(namewoext);
	}
	return wcsname;
}

disto_data *init_disto_data(disto_params *distoparam, sequence *seq, struct wcsprm *WCSDATA, gboolean drizzle, int *status) {
	*status = 1;
	if (!distoparam)
		return NULL;
	struct wcsprm *wcs = NULL;
	disto_data *disto = NULL;
	switch (distoparam->index) {
		case DISTO_IMAGE:
			wcs = wcs_deepcopy(gfit.keywords.wcslib, NULL);
			gchar *wcsname = get_wcs_filename(PATHPARSE_MODE_WRITE, seq);
			if (!wcsname || save_wcs_fits(&gfit, wcsname)) {
				siril_log_color_message(_("Could not save WCS file for distortion\n"), "red");
				wcsfree(wcs);
				return NULL;
			}
			distoparam->index = DISTO_FILE; // this will be written to the seq file
			distoparam->filename = wcsname;
			break;
		case DISTO_FILE_COMET:
			if (!distoparam->filename) // if no filename is passed, we can go on, otherwise, we do as DISTO_FILE below
				break;
		case DISTO_FILE:;
			fits fit = { 0 };
			int statusread = read_fits_metadata_from_path_first_HDU(distoparam->filename, &fit);
			if (statusread) {
				siril_log_color_message(_("Could not load FITS file for distortion\n"), "red");
				clearfits(&fit);
				return NULL;
			}
			wcs = wcs_deepcopy(fit.keywords.wcslib, &statusread);
			if (statusread) {
				siril_log_color_message(_("Could not copy WCS information for distortion\n"), "red");
				clearfits(&fit);
				return NULL;
			}
			if (!g_str_has_suffix(distoparam->filename, ".wcs")) { // we will refer to this wcs when saving
				g_free(distoparam->filename);
				gchar *wcsname = get_wcs_filename(PATHPARSE_MODE_WRITE, seq);
				distoparam->filename = wcsname;
				if (!wcsname || save_wcs_fits(&fit, wcsname)) {
					siril_log_color_message(_("Could not save WCS file for distortion\n"), "red");
					wcsfree(wcs);
					return NULL;
				}
			}
			clearfits(&fit);
			break;
		case DISTO_MASTER:
			if (!distoparam->filename) {
				if (!com.pref.prepro.disto_lib) {
					siril_log_color_message(_("Sequence file points to master distorsion file but its specification is empty in the preferences\n"), "red");
					return NULL;
				} else {
					distoparam->filename = g_strdup(com.pref.prepro.disto_lib);
				}
			} // otherwise, it was already set, we reuse the pattern saved in the seqfile
			break;
		case DISTO_FILES:
			break;
		default:
			return NULL;
	}

	if (distoparam->index == DISTO_MASTER) {
		fits fit  = { 0 };
		disto = calloc(seq->number, sizeof(disto_data));
		gboolean found = FALSE;
		for (int i = 0;  i < seq->number; i++) {
			if (!seq->imgparam[i].incl)
				continue;
			if (seq_read_frame_metadata(seq, i, &fit)) {
				siril_log_color_message(_("Could not load image# %d, aborting\n"), "red", i + 1);
				seq->imgparam[i].incl = FALSE;
				clearfits(&fit);
				free_disto_args(disto);
				return NULL;
			}
			int statusread = 0;
			gchar *wcsname = path_parse(&fit, distoparam->filename, PATHPARSE_MODE_READ, &statusread);
			clearfits(&fit);
			if (statusread) {
				siril_log_color_message(_("Could not parse master file name for distortion, aborting\n"), "red");
				free_disto_args(disto);
				return NULL;
			}
			statusread = read_fits_metadata_from_path_first_HDU(wcsname, &fit);
			if (statusread) {
				siril_log_color_message(_("Could not load master file for distortion, aborting\n"), "red");
				clearfits(&fit);
				free_disto_args(disto);
				return NULL;
			}
			wcs = wcs_deepcopy(fit.keywords.wcslib, &statusread);
			clearfits(&fit);
			if (statusread) {
				siril_log_color_message(_("Could not copy WCS information for distortion, aborting\n"), "red");
				free_disto_args(disto);
				return NULL;
			}
			if (wcs->lin.dispre) {
				found = TRUE;
				disto[i].order = extract_SIP_order_and_matrices(wcs->lin.dispre, disto[i].A, disto[i].B, disto[i].AP, disto[i].BP);
				disto[i].xref = wcs->crpix[0] - 1.; // -1 comes from the difference of convention between opencv and wcs
				disto[i].yref = wcs->crpix[1] - 1.;
				disto[i].dtype = (drizzle) ? DISTO_S2D: DISTO_D2S;
			} else {
				disto[i].dtype = DISTO_NONE;
			}
			wcsfree(wcs);
			wcs = NULL;
		}
		if (!found) {
			free_disto_args(disto);
			distoparam->index = DISTO_UNDEF;
			return NULL;
		}
	}

	// we only have one disto spec, we can init the disto structure
	if (distoparam->index == DISTO_IMAGE || distoparam->index == DISTO_FILE || (distoparam->index == DISTO_FILE_COMET && distoparam->filename != NULL)) {
		if (!wcs)
			return disto;
		if (!wcs->lin.dispre) {
			siril_log_color_message(_("Selected file has no distortion information\n"), "red");
			wcsfree(wcs);
			g_free(distoparam->filename);
			distoparam->index = DISTO_UNDEF;
			return NULL;
		}
		disto = calloc(1, sizeof(disto_data));
		disto[0].order = extract_SIP_order_and_matrices(wcs->lin.dispre, disto[0].A, disto[0].B, disto[0].AP, disto[0].BP);
		disto[0].xref = wcs->crpix[0] - 1.; // -1 comes from the difference of convention between opencv and wcs
		disto[0].yref = wcs->crpix[1] - 1.;
		disto[0].dtype = (drizzle) ? DISTO_MAP_S2D: DISTO_MAP_D2S;
		wcsfree(wcs);
	}
	if (distoparam->index == DISTO_FILES) {
		disto = calloc(seq->number, sizeof(disto_data));
		gboolean found = FALSE;
		for (int i = 0;  i < seq->number; i++) {
			if (!seq->imgparam[i].incl)
				continue;
			if (WCSDATA[i].lin.dispre) {
				found = TRUE;
				disto[i].order = extract_SIP_order_and_matrices(WCSDATA[i].lin.dispre, disto[i].A, disto[i].B, disto[i].AP, disto[i].BP);
				disto[i].xref = WCSDATA[i].crpix[0] - 1.; // -1 comes from the difference of convention between opencv and wcs
				disto[i].yref = WCSDATA[i].crpix[1] - 1.;
				disto[i].dtype = (drizzle) ? DISTO_S2D: DISTO_D2S;
			} else {
				disto[i].dtype = DISTO_NONE;
			}
			wcsfree(WCSDATA + i);
		}
		// we haven't found any disto
		// this case happens when we have used astrometric registration with no SIP
		// it's ok, we'll just set regargs->undistort to DISTO_UNDEF
		// so as not to use maps and use optimized image transform instead
		if (!found) {
			free_disto_args(disto);
			distoparam->index = DISTO_UNDEF;
			disto = NULL;
		}
	}
	*status = 0;
	if (distoparam->index != DISTO_FILE_COMET || distoparam->filename != NULL)
		siril_log_color_message(_("Distortion data is valid and will be used\n"), "green");
	return disto;
}

gboolean validate_disto_params(fits *reffit, const gchar *text, disto_source index, gchar **msg1, gchar **msg2) {
	if (index == DISTO_UNDEF)
		return TRUE;
	if (index == DISTO_IMAGE) {
		if (!has_wcs(reffit)) {
			*msg1 = g_strdup(_("You have selected undistortion from current image but it is not platesolved, perform astrometry first or disable distortion"));
			if (msg2)
				*msg2 = g_strdup(_("Platesolve current image"));
			return FALSE;
		}
		if (!reffit->keywords.wcslib->lin.dispre) {
			*msg1 = g_strdup(_("You have selected undistortion from current image but it is not platesolved, perform astrometry first or disable distortion"));
			if (msg2)
				*msg2 = g_strdup(_("Platesolve current image with distortions"));
			return FALSE;
		}
	}
	if (index == DISTO_FILE) {
		if (!text || *text == '\0') {
			*msg1 = g_strdup(_("You have selected undistortion from an existing file, you need to specify the filename"));
			if (msg2)
				*msg2 = g_strdup(_("Load a FITS/WCS file for distortion"));
			return FALSE;
		} else {
			fits fit = { 0 };
			if (read_fits_metadata_from_path_first_HDU(text, &fit)) {
				*msg1 = g_strdup(_("Could not load FITS image for distortion"));
				if (msg2)
					*msg2 = g_strdup(_("Could not load FITS image for distortion"));
				clearfits(&fit);
				return FALSE;
			}
			if (!has_wcs(&fit)) {
				*msg1 = g_strdup(_("You have selected undistortion from file but it is not platesolved, perform astrometry first or disable distortion"));
				if (msg2)
					*msg2 = g_strdup(_("Selected file has no WCS information"));
				clearfits(&fit);
				return FALSE;
			}
			if (!fit.keywords.wcslib->lin.dispre) {
				*msg1 = g_strdup(_("You have selected undistortion from file but it is has no distortion terms, perform astrometry with SIP enabled or disable distortion"));
				if (msg2)
					*msg2 = g_strdup(_("Selected file has no distortion information"));
				clearfits(&fit);
				return FALSE;
			}
			int rx = fit.rx;
			int ry = fit.ry;
			if ((rx == 0 || ry == 0) && parse_wcs_image_dimensions(&fit, &rx, &ry)) {
				*msg1 = g_strdup(_("Selected file has no dimensions information"));
				if (msg2)
					*msg2 = g_strdup(_("Selected file has no dimensions information"));
				clearfits(&fit);
				return FALSE;
			}
			if (rx != reffit->rx || ry != reffit->ry) {
				*msg1 = g_strdup(_("Selected file and current sequence do not have the same size"));
				if (msg2)
					*msg2 = g_strdup(_("Selected file and current sequence do not have the same size"));
				clearfits(&fit);
				return FALSE;
			}
			clearfits(&fit);
		}
	}
	if (index == DISTO_MASTER) {
		if (!com.pref.prepro.disto_lib || com.pref.prepro.disto_lib[0] == '\0') {
			*msg1 = g_strdup(_("You need to set a distortion master template in Preferences/Preprocessing"));
			if (msg2)
				*msg2 = g_strdup(_("You need to set a distortion master template in Preferences/Preprocessing"));
			return FALSE;
		}
		gchar *wcsfilename = get_wcs_filename(PATHPARSE_MODE_READ, NULL);
		if (!wcsfilename) {
			*msg1 = g_strdup(_("Distortion master could not be parsed"));
			if (msg2)
				*msg2 = g_strdup(_("Distortion master could not be parsed"));
			return FALSE;
		}
		g_free(wcsfilename);
	}
	return TRUE;
}

void free_disto_args(disto_data *disto) {
	if (!disto)
		return;
	// we only need to free the maps for the 2 types which store them (disto has only one element in that case)
	if (disto->dtype == DISTO_MAP_D2S || disto->dtype == DISTO_MAP_S2D) {
		free(disto->xmap);
		free(disto->ymap);
	}
	free(disto);
}

void copy_disto(disto_data *disto_in, disto_data *disto_out) {
	if (!disto_in || !disto_out)
		return;
	disto_out->dtype = disto_in->dtype;
	memcpy(disto_out->A,  disto_in->A,  sizeof(double) * MAX_DISTO_SIZE * MAX_DISTO_SIZE);
	memcpy(disto_out->B,  disto_in->B,  sizeof(double) * MAX_DISTO_SIZE * MAX_DISTO_SIZE);
	memcpy(disto_out->AP, disto_in->AP, sizeof(double) * MAX_DISTO_SIZE * MAX_DISTO_SIZE);
	memcpy(disto_out->BP, disto_in->BP, sizeof(double) * MAX_DISTO_SIZE * MAX_DISTO_SIZE);
	disto_out->order = disto_in->order;
	disto_out->xref = disto_in->xref;
	disto_out->yref = disto_in->yref;
	disto_out->xmap = NULL;
	disto_out->ymap = NULL;
}
