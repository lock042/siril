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
#ifdef __cplusplus
}
#endif
#include "tracks.h"
#endif // SELF_CONTAINED
#include "core/siril_log.h"

#ifdef SELF_CONTAINED
#include <vector>
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
		}
		Segment() { };	// for vector allocation

		double toDegrees() { return theta_deg; }
		double angle() { return theta; }

		pointi start() { pointi p = { x1, y1 }; return p; }
		pointi end() { pointi p = { x2, y2 }; return p; }
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

#define ANGLES_EPSILON	15.0	// degrees
#define POS_EPSILON 11	// pixels
#define RHO_EPSILON 7.0	// something like pixels perpendicular to lines
// TODO: make POS_EPSILON variable, based on sampling

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

// we project point A of s2 on s1, check if it falls inside s1 and if the projected length
// is <= POS_EPSILON
// from https://stackoverflow.com/a/73437072
static bool is_projection_close(Segment &s1, pointi s2a, pointi s2b) {
	int s1x = s1.x2 - s1.x1;
	int s1y = s1.y2 - s1.y1;
	int s2x = s2a.x - s2b.x;
	int s2y = s2a.y - s2b.y;
	int denom = s1x * s2x + s1y * s2y;
	if (!denom)
		return false; // segments are orthogonal
	double t = (s2x * (s2a.x - s1.x1) + s2y * (s2a.y - s1.y1)) / (double)denom;
	siril_debug_print("t: %f\n", t);
	if (t < 0. || t > 1.)
		return false; // projected point is outside s1
	int qx = (int)(s1.x1 + s1x * t + 0.5);
	int qy = (int)(s1.y1 + s1y * t + 0.5);

	double dist = sqrt((s2a.x - qx) * (s2a.x - qx) + (s2a.y - qy) * (s2a.y - qy));
	siril_debug_print("distance: %f\n", dist);
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

// prerequisite: segments[i] and segments[j] have similar angles
static bool segments_are_closely_colinear(Segment &s1, Segment &s2) {
	pointi s1a = s1.start(), s1b = s1.end(), s2a = s2.start(), s2b = s2.end();
	if (s1a.x < s1b.x) {
		s1a.x -= POS_EPSILON; s1b.x += POS_EPSILON;
	}
	else if (s1a.x > s1b.x) {
		s1a.x += POS_EPSILON; s1b.x -= POS_EPSILON;
	}
	if (s1a.y < s1b.y) {
		s1a.y -= POS_EPSILON; s1b.y += POS_EPSILON;
	}
	else if (s1a.y > s1b.y) {
		s1a.y += POS_EPSILON; s1b.y -= POS_EPSILON;
	}
	if (s2a.x < s2b.x) {
		s2a.x -= POS_EPSILON; s2b.x += POS_EPSILON;
	}
	else if (s2a.x > s2b.x) {
		s2a.x += POS_EPSILON; s2b.x -= POS_EPSILON;
	}
	if (s2a.y < s2b.y) {
		s2a.y -= POS_EPSILON; s2b.y += POS_EPSILON;
	}
	else if (s2a.y > s2b.y) {
		s2a.y += POS_EPSILON; s2b.y -= POS_EPSILON;
	}
	// warning, they may be out of image bounds at this stage

	if ((s2a.x >= s1a.x && s2a.x <= s1b.x && s2a.y >= s1a.y && s2a.y <= s1b.y) ||
			(s2b.x >= s1a.x && s2b.x <= s1b.x && s2b.y >= s1a.y && s2b.y <= s1b.y) ||
			(s1a.x >= s2a.x && s1a.x <= s2b.x && s1a.y >= s2a.y && s1a.y <= s2b.y) ||
			(s1b.x >= s2a.x && s1b.x <= s2b.x && s1b.y >= s2a.y && s1b.y <= s2b.y)) {
		bool test = check_distance_by_projection(s1, s2);
		if (test)
			siril_debug_print("  one point of a segment is inside another and segments project closely\n");
		else siril_debug_print("  one point of a segment is inside another but segments are not close\n");
		return test;
	}
	return false;
}

/*static size_t remove_segments_on_the_sides(std::vector<Segment> &segments, size_t nb_lines, int rx, int ry) {
	size_t j = 0;
	for (size_t i = 0; i < nb_lines; i++) {
		int x1 = segments[i].x1, x2 = segments[i].x2, y1 = segments[i].y1, y2 = segments[i].y2;
		if (x1 > 2 && x2 > 2 && x1 < rx-3 && x2 < rx-3 && y1 > 2 && y2 > 2 && y1 < ry-3 && y2 < ry-3) {
			if (i != j) {
				segments[j] = segments[i];
			}
			j++;
		}
	}
	size_t new_size = j;
	for (; j < nb_lines; j++)
		segments.pop_back();
	siril_debug_print("kept %zd lines not on the sides\n", new_size);
	return new_size;
}*/

// doing this in a separate function allows the main loop coordinates to not be
// changed and the global value to not be overwritten without check
static void merge_two_segments(std::vector<Segment> &segments, size_t i, size_t j) {
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
	//siril_debug_print("   rho for merged was %f and %f\n", segments[i].rho, segments[j].rho);
	segments[i] = Segment(newx1, newy1, newx2, newy2);
	siril_debug_print("   updated segment %zd to (%d,%d) -> (%d,%d)\n", i, newx1, newy1, newx2, newy2);
}

static size_t remove_duplicate_segments(std::vector<Segment> &segments) {
	size_t nb_lines = segments.size();
	// compute segment lengths
	/*std::vector<double> lengths(nb_lines);
	for (size_t i = 0; i < nb_lines; i++) {
		int x1 = lines[i][0], y1 = lines[i][1], x2 = lines[i][2], y2 = lines[i][3];
		lengths[i] = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	}*/

	// O(n^3), this can be very slow
	size_t new_size = nb_lines;
	do {
		nb_lines = new_size;
		std::vector<bool> kept(nb_lines, true);
		std::vector<Segment> tmpsegments(segments); // copy that will contain segments for next iteration
		for (size_t i = 0; i < nb_lines; i++) {
			//if (!kept[i]) continue;
			int x1 = segments[i].x1, x2 = segments[i].x2, y1 = segments[i].y1, y2 = segments[i].y2;
			siril_debug_print("considering %ssegment %zd (%d,%d) -> (%d,%d)\n",
					kept[i]?"":"removed ", i, x1, y1, x2, y2);
			for (size_t j = i+1; j < nb_lines; j++) {
				if (fabs(segments[i].toDegrees() - segments[j].toDegrees()) <= ANGLES_EPSILON) {
					int x3 = segments[j].x1, x4 = segments[j].x2, y3 = segments[j].y1, y4 = segments[j].y2;
					siril_debug_print("  looking at %ssegment %zd (%d,%d) -> (%d,%d)\n",
							kept[j]?"":"removed ", j, x3, y3, x4, y4);
					// check if segments have almost the same start or end
					if ((abs(x1 - x3) <= POS_EPSILON && abs(y1 - y3) <= POS_EPSILON) ||
							(abs(x2 - x4) <= POS_EPSILON && abs(y2 - y4) <= POS_EPSILON) ||
							// check if they almost start at the end or end at the start of the other
							(abs(x1 - x4) <= POS_EPSILON && abs(y1 - y4) <= POS_EPSILON) ||
							(abs(x2 - x3) <= POS_EPSILON && abs(y2 - y3) <= POS_EPSILON) ||
							// check if they intersect
							segments_intersect(segments[i], segments[j]) ||
							// check if they are close enough
							segments_are_closely_colinear(segments[i], segments[j])) {
						if (kept[i])
							merge_two_segments(tmpsegments, i, j);
						kept[j] = false;
						siril_debug_print("  removing %zd: (%d,%d) -> (%d,%d)\n", j, x3, y3, x4, y4);
					}
				}
			}
		}

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


#ifndef SELF_CONTAINED

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

int cvHoughLines(fits *image, int layer, float threshvalue, int minlen, struct track **tracks) {
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
	HoughLinesP(binary, lines, 1.0, CV_PI / 180.0, minlen, (double)minlen, 5.0);
	size_t nb_lines = lines.size();
	if (nb_lines) {
		nb_lines = remove_segments_on_the_sides(lines, image->rx, image->ry);
		siril_debug_print("we have %zd candidate segments\n", nb_lines);

		std::vector<Segment> segments = std::vector<Segment>(nb_lines);
		const int ry = image->ry;
		std::transform(lines.begin(), lines.end(), segments.begin(), [=](Vec4i v) {
				return Segment(v[0], ry - v[1] - 1, v[2], ry - v[3] - 1);
				});
		nb_lines = remove_duplicate_segments(segments); // in-place edit of segments

		// display lines and the distribution by angle
		const int nb_bins = 72;		// 5 degree steps
		gsl_histogram *hist = gsl_histogram_alloc(nb_bins);
		gsl_histogram_set_ranges_uniform(hist, -180.0, 180.0);
		for (size_t i = 0; i < nb_lines; i++) {
			double degs = segments[i].toDegrees();
			gsl_histogram_increment(hist, degs);
			int x1 = segments[i].x1, x2 = segments[i].x2, y1 = segments[i].y1, y2 = segments[i].y2;
			siril_debug_print("Line detected (%d,%d)->(%d,%d), angle %f\n", x1, y1, x2, y2, degs);
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
	pointi s2a = { 1030, -1003 };
	pointi s2b = { 1080, 1005 };
	res = is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection -h");

	// projection, horizontal segment, far away
	s1 = Segment(1000, 1000, 1100, 1000);
	pointi s2a = { 1080, 1003 };
	pointi s2b = { 1030, 1500 };
	res = !is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection h far");

	// projection, vertical segment
	s1 = Segment(1000, 1000, 1000, 1100);
	pointi s2a = { 1003, 1030 };
	pointi s2b = { 1005, 1080 };
	res = is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection v");

	// projection, vertical segment, far away
	s1 = Segment(1000, 1000, 1000, 1100);
	pointi s2a = { 1003, 1080 };
	pointi s2b = { 1005, 1030 };
	res = !is_projection_close(s1, s2a, s2b) && is_projection_close(s1, s2b, s2a);
	RESULT("projection v far");

	// two segments slightly separated, not intersecting, close to parallel
	s1 = Segment(1309,591,1343,622);
	s2 = Segment(1322,600,1365,643);
	res = segments_are_closely_colinear(s1, s2);
	RESULT("proximity");

	return 0;
}
#endif
