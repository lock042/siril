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
 * along with Siril. If not, see <http://www.gnu.org/licenses/ >.
 */

/* Regression tests for the homography point transform used to reframe
 * astrometric solutions (reframe_wcs() and every geometry operation that
 * follows it: crop, rotate, mirror, binning).
 *
 * cvTransformImageRefPoint() used to build the homogeneous input point with
 * the implicit Mat(Point3d) constructor, whose shape is 3x1 with OpenCV 4 but
 * 1x3 with OpenCV 5. Against OpenCV 5 the matrix product therefore threw an
 * uncaught cv::Exception from gemm(), aborting Siril at the end of a plate
 * solve. These tests pin both the shape and the operand order.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "core/siril.h"
#include "opencv/opencv.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image (now a pointer)

#define TOL 1e-9

/* identity must be a no-op */
Test(homography, refpoint_identity) {
	Homography H = { 0 };
	cvGetEye(&H);
	point in = { 1503.5, 1481.2 }, out = { 0., 0. };

	cvTransformImageRefPoint(H, in, &out);

	cr_assert_float_eq(out.x, in.x, TOL, "x changed under identity");
	cr_assert_float_eq(out.y, in.y, TOL, "y changed under identity");
}

/* the vertical flip built by flip_bottom_up_astrometry_data(): y -> ry - y.
 * This is the transform that crashed in issue "Astrometry causes Siril to
 * freeze or crash", so it is the direct regression case. */
Test(homography, refpoint_bottom_up_flip) {
	const double ry = 3008.;
	Homography H = { 0 };
	cvGetEye(&H);
	H.h11 = -1.;
	H.h12 = ry;
	point in = { 1503.5, 1481.2 }, out = { 0., 0. };

	cvTransformImageRefPoint(H, in, &out);

	cr_assert_float_eq(out.x, 1503.5, TOL, "x must be unaffected by a y flip");
	cr_assert_float_eq(out.y, ry - 1481.2, TOL, "y not mirrored about ry");

	/* flipping twice returns to the starting point */
	point back = { 0., 0. };
	cvTransformImageRefPoint(H, out, &back);
	cr_assert_float_eq(back.x, in.x, TOL, "x not restored by double flip");
	cr_assert_float_eq(back.y, in.y, TOL, "y not restored by double flip");
}

/* a fully asymmetric projective transform: catches a transposed H (which an
 * identity or pure-flip case would not) and exercises the w division */
Test(homography, refpoint_projective) {
	Homography H = { 0 };
	H.h00 =  2.0;   H.h01 = 0.5;    H.h02 = 10.0;
	H.h10 = -0.25;  H.h11 = 1.5;    H.h12 = -4.0;
	H.h20 =  1e-4;  H.h21 = 2e-4;   H.h22 = 1.0;
	point in = { 100., 200. }, out = { 0., 0. };

	/* u = 310, v = 271, w = 1.05 */
	const double expected_x = 310. / 1.05;
	const double expected_y = 271. / 1.05;

	cvTransformImageRefPoint(H, in, &out);

	cr_assert_float_eq(out.x, expected_x, TOL, "x: got %.12f, expected %.12f",
			out.x, expected_x);
	cr_assert_float_eq(out.y, expected_y, TOL, "y: got %.12f, expected %.12f",
			out.y, expected_y);
}

/* transforming by H then by H^-1 must be the identity */
Test(homography, refpoint_inverse_round_trip) {
	Homography H = { 0 };
	H.h00 =  2.0;   H.h01 = 0.5;    H.h02 = 10.0;
	H.h10 = -0.25;  H.h11 = 1.5;    H.h12 = -4.0;
	H.h20 =  1e-4;  H.h21 = 2e-4;   H.h22 = 1.0;
	Homography Hinv = H;
	cvInvertH(&Hinv);

	point in = { 100., 200. }, mid = { 0., 0. }, back = { 0., 0. };
	cvTransformImageRefPoint(H, in, &mid);
	cvTransformImageRefPoint(Hinv, mid, &back);

	cr_assert_float_eq(back.x, in.x, 1e-6, "x not restored by H^-1");
	cr_assert_float_eq(back.y, in.y, 1e-6, "y not restored by H^-1");
}
