#define _USE_MATH_DEFINES
#include <math.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "tracks.h"
#include "core/siril_log.h"

#ifdef __cplusplus
extern "C" {
#endif
//#include "algos/statistics.h"
#include <gsl/gsl_histogram.h>
#ifdef __cplusplus
}
#endif

using namespace cv;

class Segment {
	private:
		double theta;	// in radians
		double sin_theta;
		double theta_deg;
		double m;

	public:
		int x1, y1, x2, y2;
		double rho;

		Segment(int startx, int starty, int endx, int endy) {
			x1 = startx;
			y1 = starty;
			x2 = endx;
			y2 = endy;

			double dx = x1 - x2;
			double dy = y1 - y2;
			m = dy / dx;
			theta = atan2(dy, dx);
			theta_deg = theta / M_PI * 180.0;
			sin_theta = ::sin(theta);
			rho = fabs(y1 - m * x1) / sqrt(m * m + 1);
		}
		Segment() { };	// for vector allocation

		double toDegrees() { return theta_deg; }
		double angle() { return theta; }
		double sin() { return sin_theta; }

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
#define RHO_EPSILON 2.0	// something like pixels perpendicular to lines
// TODO: make POS_EPSILON variable, based on sampling

// https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// prerequisite: segments[i] and segments[j] have similar angles
static bool segments_intersect_or_almost(Segment s1, Segment s2) {
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
	// and close colinear cases (same theta and rho, one point inside the segment of the other)
	return false;
	// FIXME	v  this formula is wrong  v
	if (fabs((s1.rho - s2.rho) * s1.sin()) > RHO_EPSILON) return false;
	if (s2.x1 >= s1.x1 && s2.x1 <= s1.x2 && s2.y1 >= s1.y1 && s2.y1 <= s1.y2) return true;
	if (s2.x2 >= s1.x1 && s2.x2 <= s1.x2 && s2.y2 >= s1.y1 && s2.y2 <= s1.y2) return true;
	if (s1.x1 >= s2.x1 && s1.x1 <= s2.x2 && s1.y1 >= s2.y1 && s1.y1 <= s2.y2) return true;
	if (s1.x2 >= s2.x1 && s1.x2 <= s2.x2 && s1.y2 >= s2.y1 && s1.y2 <= s2.y2) return true;
	return false;
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
	siril_debug_print("   rho for merged was %f and %f, sin delta rho: %f\n", segments[i].rho,
			segments[j].rho, fabs((segments[i].rho - segments[j].rho) * segments[i].sin()));
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
		for (size_t i = 0; i < nb_lines; i++) {
			if (!kept[i]) continue;
			int x1 = segments[i].x1, x2 = segments[i].x2, y1 = segments[i].y1, y2 = segments[i].y2;
			siril_debug_print("considering %ssegment %zd (%d,%d) -> (%d,%d)\n",
					kept[i]?"":"removed ", i, x1, y1, x2, y2);
			for (size_t j = i+1; j < nb_lines; j++) {
				//if (!kept[j]) continue; // this may cause segments similar to removed to be kept
				//if (i == 0 && j == 12)
				//	siril_debug_print("rho: %f and %f, sin delta rho: %f\n", segments[i].rho,
				//			segments[j].rho, fabs((segments[i].rho - segments[j].rho) * segments[i].sin()));
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
							// WARNING: this uses the new merged coordinates for i
							segments_intersect_or_almost(segments[i], segments[j])) {
						if (kept[i])
							merge_two_segments(segments, i, j);
						kept[j] = false;
						siril_debug_print("  removing %zd: (%d,%d) -> (%d,%d)\n", j, x3, y3, x4, y4);
					}
				}
			}
		}

		size_t j = 0;
		for (size_t i = 0; i < nb_lines; i++) {
			if (kept[i])
				segments[j++] = segments[i];
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

#define PROBABILITSIC_HOUGH

int cvHoughLines(fits *image, int layer, float threshvalue, int minlen, struct track **tracks) {
	if (layer < 0 || layer >= image->naxes[2]) {
		layer = (image->naxes[2] == 3) ? 1 : 0;
		siril_log_message(_("Using layer %d\n"), layer);
	}

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

#ifdef PROBABILITSIC_HOUGH
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

#else
	std::vector<Vec2f> lines; // will hold the results of the detection
	HoughLines(binary, lines, 1.0, CV_PI / 180.0, minlen);

#ifdef SIRIL_OUTPUT_DEBUG
	for (size_t i = 0; i < lines.size(); i++) {
		float rho = lines[i][0], theta = lines[i][1];
		// rho is the distance to origin in pixels, theta the angle to X axis in radian
		Point pt1, pt2;
		double a = cos(theta), b = sin(theta);
		double x0 = a*rho, y0 = b*rho;
		pt1.x = cvRound(x0 + 1000*(-b));
		pt1.y = cvRound(y0 + 1000*(a));
		pt2.x = cvRound(x0 - 1000*(-b));
		pt2.y = cvRound(y0 - 1000*(a));
		siril_debug_print("Line detected (%d,%d)->(%d,%d)\n", pt1.x, pt1.y, pt2.x, pt2.y);
	}
#endif
#endif

	binary.release();
	free(buffer);
	return nb_lines;
}

