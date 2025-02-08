#ifndef _SIRIL_STATS_H
#define _SIRIL_STATS_H

struct stat_data {
	fits *fit;
	int option;
	rectangle selection;
	gchar **list;
	sequence *seq;
	gchar *csv_name;
	const gchar *seqEntry;	// not used for stats
	gboolean cfa;
};

#define NULL_STATS -999999.0

#define STATS_MINMAX	(1 << 0)    // min, max
#define STATS_SIGMEAN	(1 << 1)    // noise, mean, sigma
#define STATS_BASIC	(1 << 2)    // median, mean, sigma, noise, min, max
#define STATS_MAD	(1 << 3)    // median absolute deviation
#define STATS_AVGDEV	(1 << 4)    // average absolute deviation
#define STATS_BWMV	(1 << 5)    // bidweight midvariance
#define STATS_IKSS	(1 << 6)    // IKSS, used for normalisation

#define STATS_FOR_CFA	(1 << 7)    // work on CFA channels

#define STATS_MAIN	(STATS_BASIC | STATS_AVGDEV | STATS_MAD | STATS_BWMV)
#define STATS_EXTRA	(STATS_MAIN | STATS_IKSS)
#define STATS_NORM	(STATS_BASIC | STATS_MAD | STATS_IKSS)   // IKSSlite needs ngoodpix, median, mad and IKSS
#define STATS_LITENORM	(STATS_BASIC | STATS_MAD) // for faster normalization

#include "core/siril.h"

imstats* statistics(sequence *seq, int image_index, fits *fit, int layer,
		rectangle *selection, int option, threading_type threads);

int compute_means_from_flat_cfa(fits *fit, double mean[36]);

void allocate_stats(imstats **stat);
imstats* free_stats(imstats *stat);
void clear_stats(sequence *seq, int layer);
void clear_stats_bkp(sequence *seq, int layer);

void add_stats_to_fit(fits *fit, int layer, imstats *stat);
void add_stats_to_seq(sequence *seq, int image_index, int layer, imstats *stat);
void add_stats_to_seq_backup(sequence *seq, int image_index, int layer, imstats *stat);

void copy_seq_stats_to_fit(sequence *seq, int index, fits *fit);
void save_stats_from_fit(fits *fit, sequence *seq, int index);
void invalidate_stats_from_fit(fits *fit);
void full_stats_invalidation_from_fit(fits *fit);

void apply_stats_to_sequence(struct stat_data *stat_args);

float siril_stats_ushort_sd_64(const WORD data[], const int N);
float siril_stats_ushort_sd_32(const WORD data[], const int N);
float siril_stats_ushort_mad(const WORD* data, const size_t n, const double m, threading_type threads);
float siril_stats_float_sd(const float data[], const int N, float *mean);
double siril_stats_float_mad(const float *data, const size_t n, const double m, threading_type threads, float *buffer);
float siril_stats_trmean_from_sorted_data(const float trim, const float sorted_data[], const size_t stride, const size_t size);

int compute_all_channels_statistics_seqimage(sequence *seq, int image_index, fits *fit, int option,
		threading_type threading, int image_thread_id, imstats **stats);
int compute_all_channels_statistics_single_image(fits *fit, int option,
		threading_type threading, imstats **stats);

int copy_cached_stats_for_image(sequence *seq, int image, imstats **channels);

int sos_update_noise_float(float *array, long nx, long ny, long nchans, double *noise);

double robust_median_w(fits *fit, rectangle *area, int chan, float lower, float upper);
double robust_median_f(fits *fit, rectangle *area, int chan, float lower, float upper);
int quick_minmax(fits *fit, double *minval, double *maxval);
#endif
