#ifndef SRC_GUI_PHOTOMETRIC_CC_H_
#define SRC_GUI_PHOTOMETRIC_CC_H_

#include <stdio.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/PSF.h"

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
	GInputStream *bv_stream;	// the stream containing the star list with B-V values
	gboolean bg_auto;		// automatically select an area for bkg neutralization
	rectangle bg_area;		// the area for background if not bg_auto
	normalization_channel n_channel;// the reference channel for the white balance
};

void initialize_photometric_cc_dialog();
int apply_photometric_cc(struct photometric_cc_data *args);
int get_photometry_catalog();

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */
