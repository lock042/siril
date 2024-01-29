#ifndef SRC_GUI_PHOTOMETRIC_CC_H_
#define SRC_GUI_PHOTOMETRIC_CC_H_

#include <stdio.h>
#include <glib.h>

void get_spectrum_from_ui(xpsampled *spectrum, int chan);
void initialize_photometric_cc_dialog();
void initialize_spectrophotometric_cc_dialog();
int get_photometry_catalog_from_GUI();
void populate_spcc_combos();
void reset_spcc_filters();

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */
