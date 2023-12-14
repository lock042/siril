#ifndef SRC_GUI_PHOTOMETRIC_CC_H_
#define SRC_GUI_PHOTOMETRIC_CC_H_

#include <stdio.h>
#include <glib.h>

void initialize_photometric_cc_dialog();
int get_photometry_catalog_from_GUI();
gpointer plot_pcc_results(double* x_pre, double* y_pre, double* x_post, double* y_post, double* x_ref, double* y_ref, uint32_t nbr_points);

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */
