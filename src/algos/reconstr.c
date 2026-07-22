/* @(#)reconstr.c	19.1 (ES0-DMD) 02/25/03 13:34:40 */
/*===========================================================================
 Copyright (C) 1995 European Southern Observatory (ESO)

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 2 of
 the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public
 License along with this program;
 If not, see <http://www.gnu.org/licenses/>.

 Corresponding concerning ESO-MIDAS should be addressed as follows:
 Internet e-mail: midas@eso.org
 Postal address: European Southern Observatory
 Data Management Division
 Karl-Schwarzschild-Strasse 2
 D 85748 Garching bei Muenchen
 GERMANY
 ===========================================================================*/

/******************************************************************************
 **              Copyright (C) 1993 by European Southern Observatory
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 19.1
 **
 **    Author: Jean-Luc Starck
 **
 **    Date:  03/02/25
 **
 **    File:  reconstr.c
 **
 *******************************************************************************
 **
 **    DESCRIPTION  reconstruction routines
 **    -----------
 **
 ******************************************************************************
 **
 ** wavelet_reconstruct_file (File_Name_Imag,
 **                           File_Name_Transform, Build_Direct_Ok)
 ** char *File_Name_Imag, *File_Name_Transform;
 ** Build_Direct_Ok;
 **
 ** Reconstructs an image from its the wavelet transform of file name
 ** File_Name_Transform and write the result in the file File_Name_Imag
 **
 ** File_Name_Imag = File name of the output image
 ** File_Name_Transform = File name of the input wavelet transform
 ** Build_Direct_Ok = input parameter (TRUE=1 or FALSE=0)
 **    if the wavelet transform algorithm used the FFT and is pyramidal
 **        (Type_Wave_Transform =  TO_PYR_FFT_DIFF_RESOL
 **                             or TO_PYR_FFT_DIFF_SQUARE_RESOL)
 **        then if Build_Direct_Ok = 1 (TRUE) then
 **                       the reconstruction is done by addition
 **                       of the wavelet coefficient in the Fourier space
 **             else (Build_Direct_Ok = 0 (FALSE))
 **                       the reconstruction is done from a least mean square
 **                       estimation
 **
 ******************************************************************************
 **
 ** wavelet_reconstruct_data (Wavelet, Imag, Build_Direct_Ok)
 ** float *Imag;
 ** wave_transf_des *Wavelet;
 ** int Build_Direct_Ok;
 **
 ** Reconstructs an image from its the wavelet transform
 **
 ** Imag = OUTPUT:image
 ** Wavelet = INPUT:wavelet
 ** Build_Direct_Ok = input parameter (TRUE=1 or FALSE=0)
 **    if Wavelet->Type_Wave_Transform =  TO_PYR_FFT_DIFF_RESOL
 **                                    or TO_PYR_FFT_DIFF_SQUARE_RESOL
 **        then if Build_Direct_Ok = 1 (TRUE) then
 **                       the reconstruction is done by addition
 **                       of the wavelet coefficient in the Fourier space
 **             else (Build_Direct_Ok = 0 (FALSE))
 **                       the reconstruction is done from a least mean square
 **                       estimation

 ******************************************************************************
 **
 ** int W_Pyr_Rec_Iter_Number = 1;
 **
 ** Global variable used to define the number of iterations in
 ** Van Ciitert's iterative reconstruction. This parameter is used if
 **  Wavelet->Type_Wave_Transform = TO_PYR_BSPLINE
 **                              or TO_PYR_LINEAR
 **
 ** This parameter can be modified by an external program by:
 **
 ** extern int W_Pyr_Rec_Iter_Number;
 ** W_Pyr_Rec_Iter_Number = Value;
 **
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <glib/gstdio.h>

#include "core/siril_log.h"
#include "core/siril.h"
#include "core/proto.h"
#include "algos/Def_Math.h"
#include "algos/Def_Mem.h"
#include "algos/Def_Wavelet.h"
#include "algos/wavelet_denoise.h"

int reget_rawdata(float *Imag, int Nl, int Nc, WORD *buf) {
	float *im = Imag;
	float maximum = 0.f;
	double ratio;
	int i;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(max:maximum)
#endif
	for (i = 0; i < Nl * Nc; ++i) {
		if (im[i] > maximum)
			maximum = im[i];
	}
	ratio = (maximum > USHRT_MAX) ? USHRT_MAX_DOUBLE / maximum : 1.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (i = 0; i < Nl * Nc; ++i) {
		buf[i] = round_to_WORD(im[i] * ratio);
	}
	return 0;
}

/*****************************************************************************/

int wavelet_reconstruct_file(char *File_Name_Transform, float *coef, const struct denoise_params *dp, WORD *data) {
	float *Imag;
	wave_transf_des Wavelet;
	int Nl, Nc;

	/* read the wavelet file */
	if (wave_io_read(File_Name_Transform, &Wavelet))
		return 1;

	Nl = Wavelet.Nbr_Ligne;
	Nc = Wavelet.Nbr_Col;
	if (Nl < 1 || Nc < 1 || Nl > MAX_IMAGE_DIM || Nc > MAX_IMAGE_DIM) {
		siril_log_error(_("Error: dimensions reported by wavelets file are negative, zero or excessive.\n"));
		return 1;
	}
	if (Wavelet.Nbr_Plan < 1 || Wavelet.Nbr_Plan > 6) {
		siril_log_error(_("Error: number of plans reported by wavelets file is out of bounds.\n"));
		return 1;
	}
	Imag = f_vector_alloc(Nl * Nc);
	if (Imag == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	/* Optional per-scale denoising on the coefficient planes before synthesis;
	 * the per-scale amplitude coef[] is still applied afterwards. */
	wavelet_denoise_planes(Wavelet.Pave.Data, Wavelet.Type_Wave_Transform,
			Wavelet.Nbr_Plan, Nl, Nc, dp);
	wavelet_reconstruct_data(&Wavelet, Imag, coef);

	/* invert the Anscombe VST applied at decomposition, back to linear ADU */
	if (dp && dp->anscombe)
		anscombe_inverse(Imag, (size_t) Nl * Nc, ANSCOMBE_USHORT_SCALE);

	/* get and view result */
	reget_rawdata(Imag, Nl, Nc, data);

	wave_io_free(&Wavelet);
	free(Imag);
	return 0;
}

int wavelet_reconstruct_file_float(char *File_Name_Transform, float *coef, const struct denoise_params *dp, float *data) {
	wave_transf_des Wavelet;

	/* read the wavelet file */
	if (wave_io_read(File_Name_Transform, &Wavelet))
		return 1;

	// Sanity check before using values read from file
	if (Wavelet.Nbr_Plan < 1 || Wavelet.Nbr_Plan > 6) {
		siril_log_error(_("Error: number of plans reported by wavelets file is out of bounds.\n"));
		return 1;
	}
	if (Wavelet.Nbr_Ligne < 1 || Wavelet.Nbr_Ligne > MAX_IMAGE_DIM || Wavelet.Nbr_Col < 1 || Wavelet.Nbr_Col > MAX_IMAGE_DIM) {
		siril_log_error(_("Error: dimensions reported by wavelets file are negative, zero or excessive.\n"));
		return 1;
	}

	wavelet_denoise_planes(Wavelet.Pave.Data, Wavelet.Type_Wave_Transform,
			Wavelet.Nbr_Plan, Wavelet.Nbr_Ligne, Wavelet.Nbr_Col, dp);
	wavelet_reconstruct_data(&Wavelet, data, coef);

	/* invert the Anscombe VST applied at decomposition, back to linear [0,1] */
	if (dp && dp->anscombe)
		anscombe_inverse(data, (size_t) Wavelet.Nbr_Ligne * Wavelet.Nbr_Col,
				ANSCOMBE_FLOAT_SCALE);

	wave_io_free(&Wavelet);
	return 0;
}

/* Reconstruct only a region of interest, bounding both I/O and compute to the
 * selection. Only the ROI rows (plus a small margin so the bivariate-shrinkage
 * local window and the cropped reconstruction have real neighbours) are read
 * from the transform file; denoising and reconstruction then run on that small
 * sub-transform. roi_x/roi_y are the selection's top-left in top-down display
 * coordinates; the result is written into roifit's channel buffer in roifit's
 * (top-down) layout, matching populate_roi(). */
int wavelet_reconstruct_file_roi(char *File_Name_Transform, float *coef,
		const struct denoise_params *dp, int roi_x, int roi_y, int roi_w,
		int roi_h, int chan, fits *roifit) {
	const int margin = WD_BISHRINK_MARGIN;
	wave_transf_des hdr;
	FILE *f = g_fopen(File_Name_Transform, "rb");
	if (!f)
		return 1;
	if (fread(&hdr, sizeof(wave_transf_des), 1, f) != 1) {
		fclose(f);
		return 1;
	}
	const int Nl = hdr.Nbr_Ligne, Nc = hdr.Nbr_Col, Nbr_Plan = hdr.Nbr_Plan;
	const int type = hdr.Type_Wave_Transform;
	if (Nbr_Plan < 1 || Nbr_Plan > 6 || Nl < 1 || Nc < 1
			|| Nl > MAX_IMAGE_DIM || Nc > MAX_IMAGE_DIM
			|| (type != TO_PAVE_LINEAR && type != TO_PAVE_BSPLINE)) {
		fclose(f);
		return 1;
	}

	/* ROI bounding box in FITS (bottom-up) rows, padded by the margin */
	int cy0 = (Nl - roi_y - roi_h) - margin;
	int cy1 = (Nl - 1 - roi_y) + margin;
	int cx0 = roi_x - margin;
	int cx1 = roi_x + roi_w - 1 + margin;
	if (cy0 < 0) cy0 = 0;
	if (cy1 > Nl - 1) cy1 = Nl - 1;
	if (cx0 < 0) cx0 = 0;
	if (cx1 > Nc - 1) cx1 = Nc - 1;
	const int ch = cy1 - cy0 + 1, cw = cx1 - cx0 + 1;
	if (ch < 1 || cw < 1) {
		fclose(f);
		return 1;
	}

	const off_t data_off = (off_t) sizeof(wave_transf_des);
	float *rows = malloc((size_t) ch * Nc * sizeof(float)); /* one plane's ROI rows */
	float *sub = f_vector_alloc((size_t) cw * ch * Nbr_Plan);
	float *out = f_vector_alloc((size_t) cw * ch);
	if (!rows || !sub || !out) {
		free(rows); if (sub) free(sub); if (out) free(out);
		fclose(f);
		return 1;
	}

	for (int p = 0; p < Nbr_Plan; p++) {
		off_t off = data_off + ((off_t) p * Nl * Nc + (off_t) cy0 * Nc) * (off_t) sizeof(float);
		if (fseeko(f, off, SEEK_SET) != 0
				|| fread(rows, sizeof(float), (size_t) ch * Nc, f) != (size_t) ch * Nc) {
			free(rows); free(sub); free(out);
			fclose(f);
			return 1;
		}
		float *dstplane = sub + (size_t) cw * ch * p;
		for (int r = 0; r < ch; r++)
			memcpy(dstplane + (size_t) r * cw, rows + (size_t) r * Nc + cx0,
					(size_t) cw * sizeof(float));
	}
	fclose(f);
	free(rows);

	/* denoise + reconstruct only the cropped sub-transform */
	wavelet_denoise_planes(sub, type, Nbr_Plan, ch, cw, dp);
	pave_2d_build(sub, out, ch, cw, Nbr_Plan, coef);
	free(sub);
	if (dp && dp->anscombe)
		anscombe_inverse(out, (size_t) cw * ch,
				(roifit->type == DATA_USHORT) ? ANSCOMBE_USHORT_SCALE : ANSCOMBE_FLOAT_SCALE);

	/* copy the inner ROI into roifit (top-down), flipping FITS rows */
	for (int y = 0; y < roi_h; y++) {
		const int sr = (Nl - 1 - roi_y - y) - cy0; /* sub-image row (FITS order) */
		const float *srow = out + (size_t) sr * cw + (roi_x - cx0);
		if (roifit->type == DATA_USHORT) {
			WORD *drow = roifit->pdata[chan] + (size_t) y * roi_w;
			for (int x = 0; x < roi_w; x++)
				drow[x] = round_to_WORD(srow[x]);
		} else {
			memcpy(roifit->fpdata[chan] + (size_t) y * roi_w, srow,
					(size_t) roi_w * sizeof(float));
		}
	}
	free(out);
	return 0;
}

/*****************************************************************************/

int wavelet_reconstruct_data(wave_transf_des *Wavelet, float *Imag, float *coef) {
	float *Pave;
	int Nl, Nc, Nbr_Plan;

	Nl = Wavelet->Nbr_Ligne;
	Nc = Wavelet->Nbr_Col;
	Nbr_Plan = Wavelet->Nbr_Plan;
	switch (Wavelet->Type_Wave_Transform) {
	case TO_PAVE_LINEAR:
	case TO_PAVE_BSPLINE:
		Pave = Wavelet->Pave.Data;
		pave_2d_build(Pave, Imag, Nl, Nc, Nbr_Plan, coef);
		break;
	default:
		siril_log_message(_("Unknown transform\n"));
		return 1;
		break;
	}
	return 0;
}

/*****************************************************************************/
