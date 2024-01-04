#ifndef SRC_ALGOS_SPCC_H
#define SRC_ALGOS_SPCC_H

// Define the SPCC filters
spectral_intensity Johnson_B, Johnson_V, Optolong_Blue, Optolong_Green, Optolong_Red;

// Define the SPCC sensors
spectral_intensity Sony_IMX571M;

// SPCC functions
void init_spcc_filters();
void init_xpsampled(xpsampled *xps);
void init_xpsampled_from_library(xpsampled *out, spectral_intensity *in);
void multiply_xpsampled(xpsampled *result, const xpsampled *a, const xpsampled *b);
double integrate_xpsampled(const xpsampled *xps);
void si_free(spectral_intensity *si, gboolean free_struct);
gpointer spectrophotometric_cc_standalone(gpointer p);
cmsCIExyY xpsampled_to_xyY(xpsampled* xps, const int cmf);
void get_spectrum_from_args(struct photometric_cc_data *args, xpsampled* spectrum, int chan);
int spcc_colorspace_transform(struct photometric_cc_data *args);

enum {
	CMF_1931, CMF_1964
};

enum {
	RED, GREEN, BLUE
};

#endif
