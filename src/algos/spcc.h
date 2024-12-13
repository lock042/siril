#ifndef SRC_ALGOS_SPCC_H
#define SRC_ALGOS_SPCC_H

#define XPSAMPLED_MIN_WL 337.0
#define XPSAMPLED_MAX_WL 1019.0

#include "algos/photometric_cc.h"

// SPCC functions
void init_spcc_filters();
xpsampled init_xpsampled();
double compute_airmass(double z);
void fill_xpsampled_from_atmos_model(xpsampled *out, struct photometric_cc_data *args);
void init_xpsampled_from_library(xpsampled *out, spcc_object *in);
void multiply_xpsampled(xpsampled *result, const xpsampled *a, const xpsampled *b);
void multiply_xpsampled_scalar(xpsampled *a, const float b);
double integrate_xpsampled(const xpsampled *xps, const double minimum, const double maximum);
double xpsampled_wl_weighted_sum(xpsampled *a);
void flux_to_relcount(xpsampled *xps);
gpointer spectrophotometric_cc_standalone(gpointer p);
cmsCIExyY xpsampled_to_xyY(xpsampled* xps, const cmf_pref cmf, const double minwl, const double maxwl);
void get_spectrum_from_args(struct photometric_cc_data *args, xpsampled* spectrum, int chan);
int spcc_set_source_profile(struct photometric_cc_data *args);

// In io/spcc_json.c
void spcc_object_free(spcc_object *data, gboolean free_struct);
void osc_sensor_free(osc_sensor *data, gboolean free_struct);
void spcc_object_free_arrays(spcc_object *data);
void load_all_spcc_metadata();
void load_spcc_metadata_if_needed();
gboolean load_spcc_object_arrays(spcc_object *data);

enum {
	CMF_1931, CMF_1964
};

enum {
	MONO_SENSORS = 1,
	OSC_SENSORS,
	MONO_FILTERS,
	OSC_FILTERS,
	OSC_LPFS,
	WB_REFS
};

#endif
