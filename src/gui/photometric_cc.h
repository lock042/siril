#ifndef SRC_GUI_PHOTOMETRIC_CC_H_
#define SRC_GUI_PHOTOMETRIC_CC_H_

#include <stdio.h>
#include <glib.h>

typedef enum {
	IMX571M
} mono_sensor_t;

typedef enum {
	CANONT3I
} rgb_sensor_t;

typedef enum {
	FILTER_NONE,
	FILTER_L_ENHANCE,
	FILTER_DUAL,
	FILTER_QUAD,
	ANTLIA,
	ASTRODON,
	ASTRONOMIK,
	BAADER,
	CHROMA,
	OPTOLONG,
	ZWO
} filter_t;

#define MAX_OSC_FILTER FILTER_QUAD

void get_spectrum_from_ui(xpsampled *spectrum, int chan);
void initialize_photometric_cc_dialog();
void initialize_spectrophotometric_cc_dialog();
int get_photometry_catalog_from_GUI();

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */
