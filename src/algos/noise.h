#ifndef SRC_ALGOS_NOISE_H_
#define SRC_ALGOS_NOISE_H_

#include "core/siril.h"

/* Noise data from GUI */
struct noise_data {
	destructor destroy_fn;
	gboolean display_start_end;
	gboolean display_results;
	gboolean use_idle; // will free this struct, display things and call stop_processing_thread()
	fits *fit;

	double bgnoise[3];
	double mean_noise; // mean of bgnoise across channels, normalized in [0, 65535]
	struct timeval t_start;
	int retval;
};

gpointer noise_worker(gpointer p);
void evaluate_noise_in_image();
void bgnoise_async(fits *fit, gboolean display_values);
double bgnoise_await();

#endif /* SRC_ALGOS_NOISE_H_ */
