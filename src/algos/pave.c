/* @(#)pave.c	19.1 (ES0-DMD) 02/25/03 13:34:39 */
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
**                   Copyright (C) 1993 by European Southern Observatory
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
**    File:  pave.c
**
*******************************************************************************
**
**    DESCRIPTION  routines for the a trous algorithme
**    ----------
*******************************************************************************
**
** pave_2d_tfo (Imag, Pave, Nl, Nc, Nbr_Plan, Type_To)
** float *Imag, *Pave;
** int Nl, Nc, Nbr_Plan;
** int Type_To;
**
** computes the wavelet transform without reduction of the sampling
**
** Type_To = TO_PAVE_LINEAR for a linear scaling function
** Type_To = TO_PAVE_BSPLINE for a b3-spline scaling function
**
*******************************************************************************
**
** pave_2d_build (Pave, Imag, Nl, Nc, Nbr_Plan)
** float *Imag, *Pave;
** int Nl, Nc, Nbr_Plan;
**
** reconstruction of the image from its wavelet transform
**
*******************************************************************************
**
** pave_2d_extract_plan (Pave, Imag, Nl, Nc, Num_Plan)
** float *Imag, *Pave;
** int Nl, Nc, Num_Plan;
**
** extracts a plan from the wavelet transform
**
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "core/siril.h"
#include "algos/Def_Math.h"
#include "algos/Def_Mem.h"
#include "algos/Def_Wavelet.h"

/****************************************************************************/

int static test_ind(int ind, int N) {
	int Val;

	if (ind < 0)
		Val = -0;
	else {
		if (ind >= N)
			Val = N - 1;
		else
			Val = ind;
	}
	return (Val);
}

/****************************************************************************/

int pave_2d_linear_smooth(const float *Imag, float *Smooth, int Nl, int Nc,
		int Num_Plan) {

    int Step = pow(2.f, (float) Num_Plan) + 0.5f;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static) collapse(2)
#endif
	for (int i = 0; i < Nl; i++) {
		for (int j = 0; j < Nc; j++) {
            int indi1 = test_ind (i - Step, Nl);
            int indj1 = test_ind (j - Step, Nc);
            int indi2 = test_ind (i + Step, Nl);
            int indj2 = test_ind (j + Step,Nc);
            Smooth [i * Nc + j] = 1.f/16.f * (   Imag [indi1 * Nc + indj1]
                                           + Imag [indi1 * Nc + indj2]
                                           + Imag [indi2 * Nc + indj1]
                                           + Imag [indi2 * Nc + indj2])
                               + 1.f/8.f * (   Imag [indi1 * Nc + j]
                                           + Imag [i * Nc + indj1]
                                           + Imag [i * Nc + indj2]
                                           + Imag [indi2 * Nc + j])
                               + 1.f/4.f * Imag [i * Nc + j];

        }
    }
	return 0;
}

/***************************************************************************/

int pave_2d_tfo(float *Pict, float *Pave, int Nl, int Nc, int Nbr_Plan,
		int Type_To) {
	int Num_Plan, i;
	float *Imag, *Plan;

	Imag = f_vector_alloc(Nl * Nc);
	if (Imag == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	memcpy(Imag, Pict, Nl * Nc * sizeof(float));

	for (Num_Plan = 0; Num_Plan < Nbr_Plan - 1; Num_Plan++) {
		Plan = Pave + (Nl * Nc * Num_Plan);

		/* Copy */
		memcpy(Plan, Imag, Nl * Nc * sizeof(float));

		/* we smooth the image */
		switch (Type_To) {
		case TO_PAVE_LINEAR:
			pave_2d_linear_smooth(Plan, Imag, Nl, Nc, Num_Plan);
			break;
		case TO_PAVE_BSPLINE:
			pave_2d_bspline_smooth(Plan, Imag, Nl, Nc, Num_Plan);
			break;
		default:
			fprintf(stderr, "pave_2d.c: unknown transform\n");
			exit(-1);
			break;
		}

		/* computes the wavelet transform */
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (i = 0; i < Nl * Nc; i++)
			Plan[i] -= Imag[i];
	}

	/* copy the low resolution image in the transform */
	Plan = Pave + (Nl * Nc * (Nbr_Plan - 1));
	memcpy(Plan, Imag, Nl * Nc * sizeof(float));

	free((char *) Imag);

	return 0;
}

/***************************************************************************/

int pave_2d_build(float *Pave, float *Imag, int Nl, int Nc, int Nbr_Plan,
		const float *coef) {
	int Num_Plan;

	memset(Imag, 0.0, Nl * Nc * sizeof(float));

	for (Num_Plan = Nbr_Plan - 1; Num_Plan >= 0; Num_Plan--) {
		int Pos = Nl * Nc * Num_Plan;
		float *Plan = Pave + Pos;

		int i;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (i = 0; i < Nl * Nc; i++)
			Imag[i] += coef[Num_Plan] * Plan[i];
	}
	return 0;
}

/***************************************************************************/

int pave_2d_extract_plan(float *Pave, float *Imag, int Nl, int Nc, int Num_Plan) {
	int Pos;
	float *Plan;

	Pos = Nl * Nc * Num_Plan;
	Plan = Pave + Pos;

	memcpy(Imag, Plan, Nl * Nc * sizeof(float));

	return 0;
}

/***************************************************************************/

int pave_2d_bspline_smooth(const float *Imag, float *Smooth, int Nl, int Nc,
		int Num_Plan) {

	int Step = pow(2., (float) Num_Plan) + 0.5;

#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static) collapse(2)
#endif
	for (int i = 0; i < Nl; i++) {
		for (int j = 0; j < Nc; j++) {
			int indi1 = test_ind(i - Step, Nl);
			int indj1 = test_ind (j - Step, Nc);
			int indi2 = test_ind (i + Step, Nl);
			int indj2 = test_ind (j + Step,Nc);
			int indi3 = test_ind (i - 2 * Step, Nl);
			int indj3 = test_ind (j - 2 * Step,Nc);
			int indi4 = test_ind (i + 2 * Step, Nl);
			int indj4 = test_ind (j + 2 * Step,Nc);


			Smooth [i * Nc + j] = 0.00390625 * ( Imag [indi3 * Nc + indj3]
										+ Imag [indi3 * Nc + indj4]
										+ Imag [indi4 * Nc + indj3]
										+ Imag [indi4 * Nc + indj4])
							+ 0.015625 * ( Imag [indi4 * Nc + indj2]
										+ Imag [indi3 * Nc + indj2]
										+ Imag [indi4 * Nc + indj1]
										+ Imag [indi3 * Nc + indj1]

										+ Imag [indi2 * Nc + indj3]
										+ Imag [indi2 * Nc + indj4]
										+ Imag [indi1 * Nc + indj3]
										+ Imag [indi1 * Nc + indj4])

							+ 0.0234375 * ( Imag [indi3 * Nc + j]
										+ Imag [indi4 * Nc + j]
										+ Imag [i * Nc + indj3]
										+ Imag [i * Nc + indj4])

							+ 0.06250 * ( Imag [indi1 * Nc + indj1]
										+ Imag [indi1 * Nc + indj2]
										+ Imag [indi2 * Nc + indj1]
										+ Imag [indi2* Nc + indj2])

							+ 0.09375 * ( Imag [indi1 * Nc + j]
										+ Imag [indi2 * Nc + j]
										+ Imag [i * Nc + indj1]
										+ Imag [i * Nc + indj2])

							+ 0.140625 * Imag [i * Nc + j];

		}
	}
	return 0;
}

/****************************************************************************/
