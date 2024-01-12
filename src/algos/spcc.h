#ifndef SRC_ALGOS_SPCC_H
#define SRC_ALGOS_SPCC_H

// SPCC functions
void init_spcc_filters();
xpsampled init_xpsampled();
void init_xpsampled_from_library(xpsampled *out, spcc_object *in);
void multiply_xpsampled(xpsampled *result, const xpsampled *a, const xpsampled *b);
double integrate_xpsampled(const xpsampled *xps);
gpointer spectrophotometric_cc_standalone(gpointer p);
cmsCIExyY xpsampled_to_xyY(xpsampled* xps, const int cmf);
void get_spectrum_from_args(struct photometric_cc_data *args, xpsampled* spectrum, int chan);
int spcc_set_source_profile(struct photometric_cc_data *args);

// In io/spcc_json.c
void spcc_object_free(spcc_object *data, gboolean free_struct);
void osc_sensor_free(osc_sensor *data, gboolean free_struct);
void spcc_object_free_arrays(spcc_object *data);
void load_all_spcc_metadata();
gboolean load_spcc_object_arrays(spcc_object *data);

enum {
	CMF_1931, CMF_1964
};

#endif
