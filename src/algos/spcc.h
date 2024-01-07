#ifndef SRC_ALGOS_SPCC_H
#define SRC_ALGOS_SPCC_H

// Define the SPCC filters
spcc_object	Johnson_B, Johnson_V,
			Optolong_Blue, Optolong_Green, Optolong_Red,
			Chroma_Red, Chroma_Green, Chroma_Blue,
			Astrodon_RE, Astrodon_RI, Astrodon_GE, Astrodon_GI, Astrodon_B;

// Define the SPCC sensors
spcc_object Sony_IMX571M, ZWO_1600M, KAF1603ME, KAF3200, KAF8300, Sony_ICX694;

// SPCC functions
void init_spcc_filters();
void init_xpsampled(xpsampled *xps);
void init_xpsampled_from_library(xpsampled *out, spcc_object *in);
void multiply_xpsampled(xpsampled *result, const xpsampled *a, const xpsampled *b);
double integrate_xpsampled(const xpsampled *xps);
gpointer spectrophotometric_cc_standalone(gpointer p);
cmsCIExyY xpsampled_to_xyY(xpsampled* xps, const int cmf);
void get_spectrum_from_args(struct photometric_cc_data *args, xpsampled* spectrum, int chan);
int spcc_colorspace_transform(struct photometric_cc_data *args);
int check_prior_spcc(fits *fit);

// In io/spcc_json.c
void spcc_object_free(spcc_object *data, gboolean free_struct);
void spcc_object_free_arrays(spcc_object *data);
void load_all_spcc_metadata(gchar *path);
gboolean load_spcc_object_arrays(spcc_object *data);

enum {
	CMF_1931, CMF_1964
};

enum {
	RED, GREEN, BLUE
};

#endif
