#ifndef SRC_ALGOS_BACKGROUND_EXTRACTION_H_
#define SRC_ALGOS_BACKGROUND_EXTRACTION_H_

#include "core/siril.h"

/* ORDER OF POLYNOMES */
typedef enum {
	BACKGROUND_POLY_1,
	BACKGROUND_POLY_2,
	BACKGROUND_POLY_3,
	BACKGROUND_POLY_4,
} poly_order;

typedef enum {
	BACKGROUND_CORRECTION_SUBTRACT,
	BACKGROUND_CORRECTION_DIVIDE
} background_correction;

typedef enum {
	BACKGROUND_INTER_RBF = 0,
	BACKGROUND_INTER_POLY = 1,
} background_interpolation;

struct background_data {
	int nb_of_samples;
	double tolerance;
	background_correction correction;
	background_interpolation interpolation_method;
	poly_order degree;
	double smoothing;
	int threads;
	gboolean dither;
	fits *fit;
	gboolean from_ui;
	sequence *seq;
	char *seqEntry;
	gboolean is_cfa;
};

typedef struct sample {
	double median[3]; // median of each channel of the sample (if color)
	double mean; // mean of the 3 channel of the sample (if color)
	double min, max;
	size_t size;
	point position;
	gboolean valid;
} background_sample;

int get_background_sample_radius();
void free_background_sample_list(GSList *list);
GSList *generate_samples(fits *fit, int nb_per_line, double tolerance, int size, const char **error, threading_type threads);
GSList* add_background_sample(GSList *list, fits *fit, point pt);
GSList *add_background_samples(GSList *orig, fits *fit, GSList *pts);
GSList* remove_background_sample(GSList *orig, fits *fit, point pt);
int generate_background_samples(int nb_of_samples, double tolerance);
gpointer remove_gradient_from_image(gpointer p);
gpointer remove_gradient_from_cfa_image(gpointer p);
void apply_background_extraction_to_sequence(struct background_data *background_args);

gboolean background_sample_is_valid(background_sample *sample);
gdouble background_sample_get_size(background_sample *sample);
point background_sample_get_position(background_sample *sample);

void apply_background_cancel();

#endif /* SRC_ALGOS_BACKGROUND_EXTRACTION_H_ */

