/*
 *  match: a package to match lists of stars (or other items)
 *  Copyright (C) 2000  Michael William Richmond
 *
 *  Contact: Michael William Richmond
 *           Physics Department
 *           Rochester Institute of Technology
 *           85 Lomb Memorial Drive
 *           Rochester, NY  14623-5603
 *           E-mail: mwrsps@rit.edu
 *
 *
 *  This program is siril_free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/*
 * <AUTO>
 * FILE: apply_match.c
 *
 * <HTML>
 * Given
 *    - an ASCII file consisting of a list of stars, with one line per
 *        star and multiple columns of information separated by white space
 *    - the numbers of the columns containing "X" and "Y" coords of stars
 *    - a central RA and Dec, each in decimal degrees
 *    - coefficients of a TRANS structure which convert "X" and "Y"
 *        coordinates from the ASCII file's system to plate coordinates
 *        (xi, eta), which are projections onto the tangent plane
 *        centered on the central RA and Dec.
 *
 * run through the data file.  For each entry, calculate the (RA, Dec) of
 * the star, then replace the "X" and "Y" values with (RA, Dec).  Leave all
 * other information in the ASCII file as-is.
 *
 * The TRANS structure converts "X" and "Y" to (xi, eta) in one of three
 * ways.  If the user specifies a 'linear' transformation (which is the
 * default), then
 *
 *     xi = A + Bx + Cy
 *    eta = D + Ex + Fy
 *
 * In the case of 'quadratic',
 *
 *     xi =  A + Bx + Cy + Dxx + Exy + Fyy
 *    eta =  G + Hx + Iy + Jxx + Kxy + Lyy
 *
 * In the case of 'cubic',
 *
 *     xi =  A + Bx + Cy + Dxx + Exy + Fyy + Gx(xx+yy) + Hy(xx+yy)
 *    eta =  I + Jx + Ky + Lxx + Mxy + Nyy + Ox(xx+yy) + Py(xx+yy)
 *
 * where "xi" and "eta" are in radians, and measure the distance
 * of the star from the center of field from which the TRANS was calculated.
 * We assume that the given "ra" and "dec" values are the same as this
 * central position.
 *
 * We force all values of RA to lie between 0 < RA < 360
 *
 * Print the results to stdout, or place them into the file given
 * by the optional "outfile" command-line argument.
 *
 * Usage: apply_match starfile1 xcol ycol ra dec linear|quadratic|cubic
 *                    a b c d e f [g h i j k [l m n o ]] [outfile=]
 *
 * </HTML>
 * </AUTO>
 *
 * modified to assume that the TRANS structure takes (x, y) coords and
 *   turns them into (xi, eta), in radians, rather than in arcseconds.
 *   MWR 5/24/2000
 *
 * modified to handle the three cases of linear, quadratic, or cubic
 *   TRANSformations.
 *   MWR 6/11/2000
 *
 * fixed equations in proc_star_file() so that they handle properly
 *   the coordinate transformations near the celestial poles.
 *   MWR 5/19/2003
 *
 * added 10 more %s in the "sscanf" statement in proc_star_file()
 *   MWR 6/11/2008
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "apply_match.h"
#include "misc.h"
#include "degtorad.h"

#undef DEBUG           /* get some of diagnostic output */

void apply_match(double ra, double dec, double xval, double yval, TRANS *trans, double *a, double *d) {
	double r_dec;
	double z, zz, alpha, delta;
	double delta_ra, delta_dec;
	r_dec = dec * DEGTORAD;

	/*
	 * let's transform from (x,y) to (delta_ra, delta_dec),
	 * using either a linear, quadratic, or cubic transformation
	 * (signalled by the 'order' field of the TRANS)
	 */
	switch (trans->order) {
	default:
	case AT_TRANS_LINEAR:
		delta_ra  = trans->x00 + trans->x10 * xval + trans->x01 * yval;
		delta_dec = trans->y00 + trans->y10 * xval + trans->y01 * yval;
		break;
	case AT_TRANS_QUADRATIC:
		delta_ra  = trans->x00 + trans->x10 * xval + trans->x01 * yval
				  + trans->x20 * xval * xval + trans->x11 * xval * yval + trans->x02 * yval * yval;
		delta_dec = trans->y00 + trans->y10 * xval + trans->y01 * yval
				  + trans->y20 * xval * xval + trans->y11 * xval * yval + trans->y02 * yval * yval;
		break;
	case AT_TRANS_CUBIC:
		delta_ra  = trans->x00 + trans->x10 * xval + trans->x01 * yval
				  + trans->x20 * xval * xval + trans->x11 * xval * yval + trans->x02 * yval * yval
				  + trans->x30 * xval * xval * xval + trans->x21 * xval * xval * yval
				  + trans->x12 * xval * yval * yval + trans->x03 * yval * yval * yval;
		delta_dec = trans->y00 + trans->y10 * xval + trans->y01 * yval
				  + trans->y20 * xval * xval + trans->y11 * xval * yval + trans->y02 * yval * yval
				  + trans->y30 * xval * xval * xval + trans->y21 * xval * xval * yval
				  + trans->y12 * xval * yval * yval + trans->y03 * yval * yval * yval;
		break;
	case AT_TRANS_QUARTIC:
		delta_ra  = trans->x00 + trans->x10 * xval + trans->x01 * yval
				  + trans->x20 * xval * xval + trans->x11 * xval * yval + trans->x02 * yval * yval
				  + trans->x30 * xval * xval * xval + trans->x21 * xval * xval * yval
				  + trans->x12 * xval * yval * yval + trans->x03 * yval * yval * yval
				  + trans->x40 * xval * xval * xval * xval + trans->x31 * xval * xval * xval * yval
				  + trans->x22 * xval * xval * yval * yval + trans->x13 * xval * yval * yval * yval
				  + trans->x04 * yval * yval * yval * yval;
		delta_dec = trans->y00 + trans->y10 * xval + trans->y01 * yval
				  + trans->y20 * xval * xval + trans->y11 * xval * yval + trans->y02 * yval * yval
				  + trans->y30 * xval * xval * xval + trans->y21 * xval * xval * yval
				  + trans->y12 * xval * yval * yval + trans->y03 * yval * yval * yval
				  + trans->y40 * xval * xval * xval * xval + trans->y31 * xval * xval * xval * yval
				  + trans->y22 * xval * xval * yval * yval + trans->y13 * xval * yval * yval * yval
				  + trans->y04 * yval * yval * yval * yval;
		break;
	}

	/*
	 * and now convert from arcseconds to radians
	 * (convenient for calculations)
	 */
	delta_ra = (delta_ra / 3600.0) * DEGTORAD;
	delta_dec = (delta_dec / 3600.0) * DEGTORAD;

	/*
	 * we have (delta_ra, delta_dec), in radians; these give the distance
	 * of this star from the central (RA, Dec).  Now we can de-project from
	 * the tangent plane (centered on RA,Dec) and calculate the actual
	 * RA, Dec of the star (in degrees)
	 */
	z = cos(r_dec) - delta_dec * sin(r_dec);
	zz = atan2(delta_ra, z) * RADTODEG;
	alpha = zz + ra;
	delta = asin((sin(r_dec) + delta_dec * cos(r_dec)) / sqrt( 1. + delta_ra * delta_ra + delta_dec * delta_dec)) * RADTODEG;


	if (alpha < 0) {
		alpha += 360.0;
	}
	if (alpha >= 360.0) {
		alpha -= 360.0;
	}
	// this avoids converging exactly to 90. which creates a mess in wcslib
	// (change of convention when native pole is at celestial pole)
	if (delta == 90.) {
		delta -= 1.e-8;
	}
#ifdef DEBUG
	fprintf(stdout, "new RA = %10.5f, new dec = %10.5f\n", alpha, delta);
#endif

	*a = alpha;
	*d = delta;
}
