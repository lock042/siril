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

enum {
	INTER_POLY = 0,
	INTER_RBF = 1
};

struct background_data {
	int nb_of_samples;
	double tolerance;
	int correction;
    int interpolation_method;
	poly_order degree;
	double smoothing;
	int threads;
	gboolean dither;
	fits *fit;
	gboolean from_ui;
	sequence *seq;
	const gchar *seqEntry;
};

typedef struct sample background_sample;

int get_sample_radius();
void free_background_sample_list(GSList *list);
GSList* add_background_sample(GSList *list, fits *fit, point pt);
GSList* remove_background_sample(GSList *orig, fits *fit, point pt);
void generate_background_samples(int nb_of_samples, double tolerance);
gpointer remove_gradient_from_image(gpointer p);
void apply_background_extraction_to_sequence(struct background_data *background_args);

gboolean background_sample_is_valid(background_sample *sample);
gdouble background_sample_get_size(background_sample *sample);
point background_sample_get_position(background_sample *sample);

void apply_background_cancel();

#endif /* SRC_ALGOS_BACKGROUND_EXTRACTION_H_ */
