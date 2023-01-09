#include "gui/progress_and_log.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "core/siril_log.h"
// Helpful wrapper functions to call Siril functions without causing clashes with some of the img_expr functionality

void updateprogress(const char *text, double percent) {
	set_progress_bar_data(text, percent);
}

int is_thread_stopped() {
	return (get_thread_run() ? 0 : 1);
}

void sirillog(const char* text) {
	siril_log_message(_(text));
}

int updatenoise(float *array, int nx, int ny, int nchans, double *noise) {
	return sos_update_noise_float(array, (long) nx, (long) ny, (long) nchans, noise);
}
