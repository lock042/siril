/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#define _USE_MATH_DEFINES
#include <math.h>
#ifndef SELF_CONTAINED
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#ifdef __cplusplus
extern "C" {
#endif
#include <gsl/gsl_histogram.h>
#include "algos/sorting.h"
#ifdef __cplusplus
}
#endif
#include "tracks.h"
#endif // SELF_CONTAINED

#include "core/siril_log.h"

#ifdef SELF_CONTAINED
#include <vector>
#include <cfloat>
typedef struct {
	int x, y;
} pointi;

// from siril.h
#include <stdio.h>
#undef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#undef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define DEBUG_TEST 1
#define siril_debug_print(fmt, ...) \
	do { if (DEBUG_TEST) fprintf(stdout, fmt, ##__VA_ARGS__); } while (0)

using namespace std;
#else
using namespace cv;
#endif

class Segment {
	private:
		double theta;	// in radians
		double theta_deg;
		double len;

	public:
		int x1, y1, x2, y2; // don't asign new values, theta is computed from this

		Segment(int startx, int starty, int endx, int endy) {
			x1 = startx;
			y1 = starty;
			x2 = endx;
			y2 = endy;

			double dx = x1 - x2;
			double dy = y1 - y2;
			theta = atan2(dy, dx);
			theta_deg = theta / M_PI * 180.0;
			len = sqrt(dx * dx + dy * dy);
		}
		Segment() { };	// for vector allocation

		double toDegrees() { return theta_deg; }
		double angle() { return theta; }
		double length() { return len; }

		pointi start() { pointi p = { x1, y1 }; return p; }
		pointi end() { pointi p = { x2, y2 }; return p; }
};

class Rectangle {
	public:
		int x1, y1, x2, y2;

		Rectangle(Segment &s, int inflation);

		pointi start() { pointi p = { x1, y1 }; return p; }
		pointi end() { pointi p = { x2, y2 }; return p; }

		bool point_is_inside(pointi p);
};


// compute the signed triangle area
static int orientation(pointi a, pointi b, pointi c) {
	int area = (b.y - a.y) * (c.x - b.x) -
		(b.x - a.x) * (c.y - b.y);
	if (area == 0)
		return 0;
	return area > 0 ? 1 : -1;
}

static bool onSegment(pointi a, pointi b, pointi c) {
	return b.x <= max(a.x, c.x) && b.x >= min(a.x, c.x) &&
		b.y <= max(a.y, c.y) && b.y >= min(a.y, c.y);
}

#define ANGLES_EPSILON	10.0	// degrees
// TODO: ANGLES_EPSILON should be function of minimal segment length
#define POS_EPSILON 20	// pixels
// TODO: make POS_EPSILON variable, based on sampling
#define MAX_MEDIAN_LENGTH_FACTOR 1.35

// https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// prerequisite: segments[i] and segments[j] have similar angles
static bool segments_intersect(Segment &s1, Segment &s2) {
	pointi a = s1.start();
	pointi b = s1.end();
	pointi c = s2.start();
	pointi d = s2.end();
	int orientABC = orientation(a, b, c);
	int orientABD = orientation(a, b, d);
	int orientCDA = orientation(c, d, a);
	int orientCDB = orientation(c, d, b);
	if (orientABC != orientABD && orientCDA != orientCDB)
		return true;
	// handle colinearity cases too
	if (orientABC == 0 && onSegment(a, c, b)) return true;
	if (orientABD == 0 && onSegment(a, d, b)) return true;
	if (orientCDA == 0 && onSegment(c, a, d)) return true;
	if (orientCDB == 0 && onSegment(c, b, d)) return true;
	return false;
}

static bool segments_have_similar_angles(Segment &s1, Segment &s2) {
	double angle1 = s1.toDegrees(), angle2 = s2.toDegrees();
	if (angle1 < 0.0)
		return fabs(angle1 - angle2) <= ANGLES_EPSILON || fabs(angle1 + 180.0 - angle2) <= ANGLES_EPSILON;
	return fabs(angle1 - angle2) <= ANGLES_EPSILON || fabs(angle1 - 180.0 - angle2) <= ANGLES_EPSILON;
}

static bool segments_share_an_end(Segment &s1, Segment &s2) {
	int x1 = s1.x1, x2 = s1.x2, y1 = s1.y1, y2 = s1.y2;
	int x3 = s2.x1, x4 = s2.x2, y3 = s2.y1, y4 = s2.y2;
	// check if segments have almost the same start or end
	return ((abs(x1 - x3) <= POS_EPSILON && abs(y1 - y3) <= POS_EPSILON) ||
			(abs(x2 - x4) <= POS_EPSILON && abs(y2 - y4) <= POS_EPSILON) ||
			// check if they almost start at the end or end at the start of the other
			(abs(x1 - x4) <= POS_EPSILON && abs(y1 - y4) <= POS_EPSILON) ||
			(abs(x2 - x3) <= POS_EPSILON && abs(y2 - y3) <= POS_EPSILON));
}

// we project point A of s2 on s1, check if it falls inside s1 and if the projected length
// is <= POS_EPSILON
// from https://stackoverflow.com/a/73437072
static bool is_projection_close(Segment &s1, pointi s2a, pointi s2b) {
	int s1x = s1.x2 - s1.x1;
	int s1y = s1.y2 - s1.y1;
	int s2x = s2b.x - s2a.x;
	int s2y = s2b.y - s2a.y;
	int denom = s1x * s2x + s1y * s2y;
	//siril_debug_print("s1x: %d, s1y: %d, s2x: %d, s2y: %d\n", s1x, s1y, s2x, s2y);
	if (!denom)
		return false; // segments are orthogonal
	double t = (s2x * (s2a.x - s1.x1) + s2y * (s2a.y - s1.y1)) / (double)denom;
	//siril_debug_print("t: %f\n", t);
	if (t < 0. || t > 1.)
		return false; // projected point is outside s1
	int qx = (int)(s1.x1 + s1x * t + 0.5);
	int qy = (int)(s1.y1 + s1y * t + 0.5);

	double dist = sqrt((s2a.x - qx) * (s2a.x - qx) + (s2a.y - qy) * (s2a.y - qy));
	//siril_debug_print("distance: %f\n", dist);
	return dist <= (double)POS_EPSILON;
}

static bool check_distance_by_projection(Segment &s1, Segment &s2) {
	// we project the two points of s1 on s2 and the two points of s2 on s1 and if they
	// fall inside the segments, we compare the length of the projection with POS_EPSILON
	return is_projection_close(s1, s2.start(), s2.end()) ||
		is_projection_close(s1, s2.end(), s2.start()) ||
		is_projection_close(s2, s1.end(), s1.start()) ||
		is_projection_close(s2, s1.start(), s1.end());
}

// this only changes coordinates, doesn't recompute theta
// warning, they may be out of image bounds at this stage
Rectangle::Rectangle(Segment &s, int amount) {
	if (s.x1 <= s.x2) {
		if (s.y1 <= s.y2) {
			x1 = s.x1 - amount; y1 = s.y1 - amount;
			x2 = s.x2 + amount; y2 = s.y2 + amount;
		}
		else {
			x1 = s.x1 - amount; y1 = s.y1 + amount;
			x2 = s.x2 + amount; y2 = s.y2 - amount;
		}
	}
	else if (s.x2 < s.x1) {
		if (s.y1 <= s.y2) {
			x1 = s.x1 + amount; y1 = s.y1 - amount;
			x2 = s.x2 - amount; y2 = s.y2 + amount;
		}
		else {
			x1 = s.x1 + amount; y1 = s.y1 + amount;
			x2 = s.x2 - amount; y2 = s.y2 - amount;
		}
	}
}

// check if a point is inside the square containing a segment
bool Rectangle::point_is_inside(pointi p) {
	if (x1 <= x2) {
		if (y1 <= y2)
			return p.x >= x1 && p.x <= x2 && p.y >= y1 && p.y <= y2;
		else 	return p.x >= x1 && p.x <= x2 && p.y <= y1 && p.y >= y2;
	} else {
		if (y1 <= y2)
			return p.x <= x1 && p.x >= x2 && p.y >= y1 && p.y <= y2;
		else 	return p.x <= x1 && p.x >= x2 && p.y <= y1 && p.y >= y2;
	}
}

// check if two segments are close enough and parallel or almost to be considered one
// prerequisite: segments have similar angles
static bool segments_are_closely_colinear(Segment &s1, Segment &s2) {
	Rectangle gs1 = Rectangle(s1, POS_EPSILON);
	Rectangle gs2 = Rectangle(s2, POS_EPSILON);
	//siril_debug_print("grown s1 (%d,%d)->(%d,%d)\n", gs1.x1, gs1.y1, gs1.x2, gs1.y2);
	//siril_debug_print("grown s2 (%d,%d)->(%d,%d)\n", gs2.x1, gs2.y1, gs2.x2, gs2.y2);

	if (gs1.point_is_inside(s2.start()) || gs1.point_is_inside(s2.end()) ||
			gs2.point_is_inside(s1.start()) || gs2.point_is_inside(s1.end())) {
		bool test = check_distance_by_projection(s1, s2);
		if (test)
			siril_debug_print("  one point of a segment is inside another and segments project closely\n");
		else siril_debug_print("  one point of a segment is inside another but segments are not close\n");
		return test;
	}
	return false;
}

// assuming similar lengths of segment
static bool overlapping_too_much(Segment &s1, Segment &s2) {
	int x1 = s1.x1, x2 = s1.x2, y1 = s1.y1, y2 = s1.y2;
	int x3 = s2.x1, x4 = s2.x2, y3 = s2.y1, y4 = s2.y2;
	/*      *---------*                 *-------------*
	 *        *---------*    *------------*
	 *        ^ yes ^                ^ no ^
	 */
	int threshold = s1.length() / 2;
	if (abs(x1 - x2) > abs(y1 - y2)) {
		// check X
		if (x1 < x2 && x3 < x4 && (x1 + threshold > x3)) return true;
		if (x1 < x2 && x4 < x3 && (x1 + threshold > x4)) return true;
		if (x2 < x1 && x3 < x4 && (x2 + threshold > x3)) return true;
		if (x2 < x1 && x4 < x3 && (x2 + threshold > x4)) return true;
	} else {
		// check Y
		if (y1 < y2 && y3 < y4 && (y1 + threshold > y3)) return true;
		if (y1 < y2 && y4 < y3 && (y1 + threshold > y4)) return true;
		if (y2 < y1 && y3 < y4 && (y2 + threshold > y3)) return true;
		if (y2 < y1 && y4 < y3 && (y2 + threshold > y4)) return true;
	}
	return false;
}

// returns: if the merged segment was too long or not
static bool merge_two_segments(std::vector<Segment> &segments, size_t i, size_t j, double median_length) {
	int x1 = segments[i].x1, x2 = segments[i].x2, y1 = segments[i].y1, y2 = segments[i].y2;
	int x3 = segments[j].x1, x4 = segments[j].x2, y3 = segments[j].y1, y4 = segments[j].y2;
	int newx1, newx2, newy1, newy2;
	// FIXME this approach is not ideal because it can change X without changing Y, resulting
	// in a different orientation if the test that requests the merge is too loose
	if (x1 < x2) {
		newx1 = min(x1, min(x3, x4));
		newx2 = max(x2, max(x3, x4));
	} else {
		newx1 = max(x1, max(x3, x4));
		newx2 = min(x2, min(x3, x4));
	}
	if (y1 < y2) {
		newy1 = min(y1, min(y3, y4));
		newy2 = max(y2, max(y3, y4));
	} else {
		newy1 = max(y1, max(y3, y4));
		newy2 = min(y2, min(y3, y4));
	}
	Segment s = Segment(newx1, newy1, newx2, newy2);
	if (s.length() > median_length * MAX_MEDIAN_LENGTH_FACTOR) {
		siril_debug_print("   not updating segment %zd, too long: %f\n", i, s.length());
		return true;	// keep or remove both segments
	}
	segments[i] = s;
	siril_debug_print("   updated segment %zd to (%d,%d) -> (%d,%d)\n", i, newx1, newy1, newx2, newy2);
	return false;
}

#ifndef SELF_CONTAINED

static size_t remove_duplicate_segments(std::vector<Segment> &segments, bool assume_main_population) {
	size_t nb_lines = segments.size();

	// compute median segment length
	double *lengths = (double *)malloc(nb_lines * sizeof(double));
	for (size_t i = 0; i < nb_lines; i++)
		lengths[i] = segments[i].length();
	double median_length = quickmedian_double(lengths, nb_lines);
	siril_log_message(_("Median length of %zd detected segments: %f\n"), nb_lines, median_length);
	if (!assume_main_population)
		median_length = DBL_MAX; // deactivates the check
	free(lengths);

	// O(n^3), this can be very slow
	size_t new_size = nb_lines;
	do {
		nb_lines = new_size;
		std::vector<bool> kept(nb_lines, true);
		std::vector<Segment> tmpsegments(segments); // copy that will contain segments for next iteration

		for (size_t i = 0; i < nb_lines; i++) {
			//if (!kept[i]) continue;	// can still be useful to check those
			if (segments[i].length() > median_length * MAX_MEDIAN_LENGTH_FACTOR) {
				kept[i] = false;
				siril_debug_print(" removing %zd, too long\n", i);
				continue;
			}
			/*int x1 = segments[i].x1, x2 = segments[i].x2, y1 = segments[i].y1, y2 = segments[i].y2;
			siril_debug_print("considering %ssegment %zd (%d,%d) -> (%d,%d)\n",
					kept[i]?"":"removed ", i, x1, y1, x2, y2);*/
			for (size_t j = i+1; j < nb_lines; j++) {
				if (!kept[j]) continue;
				if (segments_have_similar_angles(segments[i], segments[j])) {
					/*int x3 = segments[j].x1, x4 = segments[j].x2, y3 = segments[j].y1, y4 = segments[j].y2;
					siril_debug_print("  looking at %ssegment %zd (%d,%d) -> (%d,%d)\n",
							kept[j]?"":"removed ", j, x3, y3, x4, y4);*/
					// main test for two segments with similar angle: share a point, intersect or close
					if (segments_share_an_end(segments[i], segments[j]) ||
							segments_intersect(segments[i], segments[j]) ||
							segments_are_closely_colinear(segments[i], segments[j])) {
						bool newtoolong = merge_two_segments(tmpsegments, i, j, median_length);
						if (!kept[i] && !newtoolong) {
							kept[i] = true;
							kept[j] = false;
							siril_debug_print("  segment %zd was merged in %zd, reincluding it\n", j, i);
						}
						else if (newtoolong && overlapping_too_much(segments[i], segments[j])) {
							// if they are overlapping each other more than 50% we ignore one
							kept[j] = false;
							siril_debug_print("  removing segments %zd (overlapping %zd)\n", j, i);
						}
						else {
							kept[j] = newtoolong; // if not too long, it was merged, it's useless
							if (!kept[j])
								siril_debug_print("  merged and removed %zd\n", j);
						}
					}
				}
			}
		}

		// remove unkept segments and truncate the list
		size_t j = 0;
		for (size_t i = 0; i < nb_lines; i++) {
			if (kept[i])
				segments[j++] = tmpsegments[i];
		}
		new_size = j;
		for (; j < nb_lines; j++) {
			segments.pop_back();
		}
		siril_debug_print("end of iteration: %zd lines kept\n", new_size);
	} while (nb_lines != new_size);
	siril_debug_print("kept %zd unique lines\n", new_size);
	return new_size;
}

static size_t remove_segments_on_the_sides(std::vector<Vec4i> &lines, int rx, int ry) {
	size_t nb_lines = lines.size();
	size_t j = 0;
	for (size_t i = 0; i < nb_lines; i++) {
		int x1 = lines[i][0], y1 = lines[i][1], x2 = lines[i][2], y2 = lines[i][3];
		if (x1 > 2 && x2 > 2 && x1 < rx-3 && x2 < rx-3 && y1 > 2 && y2 > 2 && y1 < ry-3 && y2 < ry-3) {
			if (i != j) {
				lines[j] = lines[i];
			}
			j++;
		}
	}
	size_t new_size = j;
	for (; j < nb_lines; j++)
		lines.pop_back();
	siril_debug_print("kept %zd lines not on the sides\n", new_size);
	return new_size;
}

int cvHoughLines(fits *image, int layer, float threshvalue, int minlen, gboolean assume_main_population, struct track **tracks) {
	if (layer < 0 || layer >= image->naxes[2]) {
		layer = (image->naxes[2] == 3) ? 1 : 0;
		siril_log_message(_("Using layer %d\n"), layer);
	}

	// making the binary image using the passed threshold
	size_t nbpixels = image->naxes[0] * image->naxes[1];
	BYTE *buffer = (BYTE *)malloc(nbpixels);
	if (image->type == DATA_USHORT) {
		for (size_t i = 0; i < nbpixels; i++)
			buffer[i] = image->data[i] > threshvalue ? 255 : 0;
	} else {
		for (size_t i = 0; i < nbpixels; i++)
			buffer[i] = image->fdata[i] > threshvalue ? 255 : 0;
	}
	Mat binary = Mat(image->ry, image->rx, CV_8UC1, buffer);

	std::vector<Vec4i> lines; // will hold the results of the detection
	HoughLinesP(binary, lines, 1.0, CV_PI / 360.0, minlen, (double)minlen, (double)POS_EPSILON);
	size_t nb_lines = lines.size();
	if (nb_lines) {
		nb_lines = remove_segments_on_the_sides(lines, image->rx, image->ry);
		siril_debug_print("we have %zd candidate segments\n", nb_lines);

		// convert to C++ segment objects
		std::vector<Segment> segments = std::vector<Segment>(nb_lines);
		const int ry = image->ry;
		std::transform(lines.begin(), lines.end(), segments.begin(), [=](Vec4i v) {
				return Segment(v[0], ry - v[1] - 1, v[2], ry - v[3] - 1);
				});
		nb_lines = remove_duplicate_segments(segments, assume_main_population); // in-place edit of segments

		// print lines and the distribution by angle
		const int nb_bins = 72;		// 5 degree steps
		gsl_histogram *hist = gsl_histogram_alloc(nb_bins);
		gsl_histogram_set_ranges_uniform(hist, -180.0, 180.0);
		for (size_t i = 0; i < nb_lines; i++) {
			double degs = segments[i].toDegrees();
			gsl_histogram_increment(hist, degs);
			double len = segments[i].length();
			int x1 = segments[i].x1, x2 = segments[i].x2, y1 = segments[i].y1, y2 = segments[i].y2;
			siril_debug_print("Line detected (%d,%d)->(%d,%d), angle %f, length %f\n",
					x1, y1, x2, y2, degs, len);
		}
		for (size_t i = 0; i < gsl_histogram_bins(hist); i++) {
			int inbin = (int)gsl_histogram_get(hist, i);
			if (inbin > 0) {
				double bin_size = 360.0 / nb_bins;
				double start = -180.0 + i * bin_size;
				siril_log_message(_("%d objects with an angle in [%f, %f]\n"), inbin, start, start + bin_size);
			}
		}
		gsl_histogram_free(hist);

		if (tracks) {	// return the segments in C form
			*tracks = (struct track *)malloc(nb_lines * sizeof(struct track));
			if (!*tracks) {
				PRINT_ALLOC_ERR;
				nb_lines = 0;
			} else {
				for (size_t i = 0; i < nb_lines; i++) {
					(*tracks)[i].start = segments[i].start();
					(*tracks)[i].end = segments[i].end();
					(*tracks)[i].angle = segments[i].toDegrees();
				}
			}
		}
	}

	binary.release();
	free(buffer);
	return nb_lines;
}

#else
#define RESULT(text) \
	if (!res) { \
		printf(text " failed\n"); \
		return 1; \
	} \
	printf(text " PASS\n");

int main() {
	// two intersecting segments, orthogonal
	Segment s1 = Segment(1000, 1000, 1100, 1000);
	Segment s2 = Segment(1050, 1050, 1050, 950);
	bool res = segments_intersect(s1, s2);
	RESULT("intersection");

	// two touching segments
	s1 = Segment(1000, 1000, 1100, 1000);
	s2 = Segment(1000, 1000, 1050, 950);
	res = segments_intersect(s1, s2);
	RESULT("touching1");

	// very close inverted
	s1 = Segment(3530, 1157, 3532, 964);
	s2 = Segment(3530, 1158, 3533, 962);
	res = segments_share_an_end(s1, s2);
	RESULT("share ends");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(1309,591,1343,622);
	s2 = Segment(1322,600,1365,643);
	res = segments_have_similar_angles(s1, s2);
	RESULT("similar angle");

	// dissimilar angle orthogonal
	s1 = Segment(1000, 1000, 1100, 1000);
	s2 = Segment(1050, 1050, 1050, 950);
	res = !segments_have_similar_angles(s1, s2);
	RESULT("dissimilar angles ortho");

	// similar angle reverted
	s1 = Segment(3530, 1157, 3532, 964);
	s2 = Segment(3530, 1158, 3533, 962);
	res = segments_have_similar_angles(s1, s2);
	RESULT("similar angle inv");

	// two touching segments, orthogonal
	s1 = Segment(1000, 1000, 1100, 1000);
	s2 = Segment(1100, 1000, 1100, 950);
	res = segments_intersect(s1, s2);
	RESULT("touching2");

	// projection, horizontal segment
	s1 = Segment(1000, 1000, 1100, 1000);
	pointi s2a = { 1030, 1003 };
	pointi s2b = { 1080, 1005 };
	res = is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection h");

	// projection, horizontal segment, backwards
	s1 = Segment(1100, 1000, 1000, 1000);
	s2a = { 1030, 997 };
	s2b = { 1080, 1005 };
	res = is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection -h");

	// projection, horizontal segment, far away
	s1 = Segment(1000, 1000, 1100, 1000);
	s2a = { 1080, 1003 };
	s2b = { 1030, 1100 };
	res = is_projection_close(s1, s2a, s2b) && !is_projection_close(s1, s2b, s2a);
	RESULT("projection h far");

	// projection, vertical segment
	s1 = Segment(1000, 1000, 1000, 1100);
	s2a = { 1003, 1030 };
	s2b = { 1005, 1080 };
	res = is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection v");

	// projection, vertical segment, far away
	s1 = Segment(1000, 1000, 1000, 1100);
	s2a = { 1070, 1080 };
	s2b = { 995, 1030 };
	res = !is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection v far");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(1309,591,1343,622);
	s2 = Segment(1322,600,1365,643);
	res = segments_are_closely_colinear(s1, s2);
	RESULT("proximity 1");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(2203,3198,2369,3136);
	s2 = Segment(2274,3172,2341,3148);
	res = segments_are_closely_colinear(s1, s2);
	RESULT("proximity 2");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(2608,2023,2792,1942);
	s2 = Segment(2644,2005,2669,1994);
	res = segments_are_closely_colinear(s1, s2);
	RESULT("proximity 3");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(1984,2241,2069,2200);
	s2 = Segment(2056,2204,2106,2182);
	res = segments_are_closely_colinear(s1, s2);
	RESULT("proximity 4");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(2397,2393,2556,2387);
	s2 = Segment(2467,2393,2641,2385);
	res = segments_are_closely_colinear(s1, s2);
	RESULT("proximity 5");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(2203,3198,2369,3136);
	s2 = Segment(2274,3172,2341,3148);
	res = overlapping_too_much(s1, s2);
	RESULT("overlap 1");

	return 0;
}
#endif
