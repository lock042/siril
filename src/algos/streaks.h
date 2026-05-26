#ifndef _TRACKING_H
#define _TRACKING_H

#include <glib.h>
#include "algos/PSF.h"
#include "opencv/tracks.h"

struct streak_result {
	int filenum;
	int target_idx;
	float angle;
	double start_ra, start_dec;	// or middle
	GDateTime *start_date;
	double middle_ra, middle_dec;
	GDateTime *middle_date;		// when computed from GPS, it's the timestamp of the middle
	gboolean middle_date_is_gps;	// of the exposure at the row where the target was found
	double end_ra, end_dec;
	GDateTime *end_date;
	double image_center_ra, image_center_dec, image_rotation;
	pointi streak_start, streak_end;
	float magnitude;
	float magnitude_error;
	gboolean magnitude_reliable;	// if false, target was out of the linear response
	gboolean magnitude_absolute;	// if true, magnitude is calibrated with stars
	float snr;
};

/* the final results, with equatorial coordinates for each frame and target */
struct result_set {
	int size;	// the allocation size of data
	GSList **data;
	gboolean has_data;	// some was added in the set
};

typedef enum { MODE_STARSTRACKED, MODE_SATTRACKED } trackmode;
typedef enum { SELECTION_NO_CHANGE, SELECTION_INCLUDE, SELECTION_EXCLUDE } selectmode;

struct streak_detection_conf {
	fits *fit;
	gboolean free_fit;
	int im_idx;	// the index in the sequence
	int filenum;	// the index in the file name, for the CSV of the sequence
	gchar *filename;
	int layer;
	int initial_segment_length; // <= 0 for unset, defaults to 60 (pixels)
	int minimum_segment_length; // <= 0 for unset, chosen from image size (pixels)
	gboolean bright_target;     // use simple linear thresholds instead of ksigma

	gboolean display_streaks;

	gboolean use_idle;	// not a sequence operation, GUI operation
	int nb_threads;

	struct result_set *results;
	float fwhm;
};

gpointer streak_detection_worker(gpointer ptr);

struct result_set *alloc_results(int size);
void dump_results(struct result_set *set, const char *filename);
void free_results(struct result_set *set);

void display_streaks(struct track *tracks, int nblines);

gboolean has_streaks(fits *fit, int layer, int nb_threads);

void ssr_internal(fits *fit, int layer, double median, double bgnoise, double pixvalue, psf_star **stars, int nb_stars);

#endif

