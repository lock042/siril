#ifndef SRC_GUI_PHOTOMETRIC_CC_H_
#define SRC_GUI_PHOTOMETRIC_CC_H_

#include <stdio.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/PSF.h"
#include "algos/photometry.h"

typedef struct struct_coeff {
	float value;
	int channel;
} coeff;

typedef enum {
	CHANNEL_HIGHEST,
	CHANNEL_MIDDLE,
	CHANNEL_LOWEST
} normalization_channel;

struct photometric_cc_data {
	fits *fit;			// the image to process
	gboolean bg_auto;		// automatically select an area for bkg neutralization
	rectangle bg_area;		// the area for background if not bg_auto
	normalization_channel n_channel;// the reference channel for the white balance

	pcc_star *stars;		// the list of stars with BV index in the image
	int nb_stars;			// the number of stars in the array
	float fwhm;			// representative FWHM for stars
};

void initialize_photometric_cc_dialog();
int photometric_cc(struct photometric_cc_data *args);
gpointer photometric_cc_standalone(gpointer p);
int get_photometry_catalog();

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */
