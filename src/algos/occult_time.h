#ifndef SRC_ALGOS_OCCULT_TIME_H_
#define SRC_ALGOS_OCCULT_TIME_H_

//#include <glib.h>
//#include <gsl/gsl_matrix.h>
//#include "core/siril.h"
//#include "core/settings.h"
//#include "algos/astrometry_solver.h"
//#include "algos/PSF.h"
//#include "algos/photometry.h"
//#include "io/siril_plot.h"

struct occultation_args {
	int start_ind_inseq;	// Index of the first image of the pulse in the original sequence
	int start_ind;	// Index of the first image of the pulse in the sorted list
	int pls_nbr;	// Number of usefull images in the 100ms pulse
	double sum_flux;	// Total flux during the pulse
	double delay_comp;	// delay time between the real PPS and the forseen PPS
};

// temporary structure. Data will be used later 
struct occ_res {
	double median_seq;	// median over the sequence
	double std_seq;		// sigma over the sequence
	double exposure;	// Computed exposure time
	int th_pls_nbr;	// Theorical pulse number in the sequence
	int det_pulses; 	// Number of detected pulses
	double hi_val;	// Maximum value in the sequence
	double lo_val;	// Minimum value in the sequence
	int valid_images;	// Number of valid images in the sequence
};

gpointer occultation_worker(gpointer arg);

int occult_curve(struct light_curve_args *lcargs);

#endif /* SRC_ALGOS_SRC_ALGOS_OCCULT_TIME_H__H_ */