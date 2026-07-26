#ifndef SRC_SYNTHSTAR_H_
#define SRC_SYNTHSTAR_H_

#include "core/processing.h"

/* Destructor-first params struct for star.synthstar / star.unclip (formerly
 * paramless): carries the NDE replay provenance of the consumed star set.
 *
 * Two provenance modes, decided at capture (NDE Convention 2 refinement):
 *  - EXPLICIT (star_auto == FALSE): the op consumed a pre-existing com.stars
 *    the user detected earlier.  That choice is pinned as stars_blob and
 *    reinstalled verbatim at replay — bit-exact.
 *  - DELEGATED (star_auto == TRUE): com.stars was empty, so the op auto-
 *    detected.  The user's intent was "find stars automatically", so replay
 *    RE-RUNS detection against the (possibly amended) upstream pixels using
 *    the detection parameters recorded here — not a frozen list.  Not bit-
 *    exact, but faithful to intent (an upstream amend re-flows into the star
 *    set, as it would had the user done that step differently originally).
 * star.unclip is always DELEGATED (reprofile never consumes com.stars). */
struct synthstar_data {
	destructor destroy_fn;   /* must be first member */
	gchar *stars_blob;       /* EXPLICIT: "x,y,A,fwhmx,fwhmy,beta,sat,profile" tuples, ':'-joined */
	gboolean star_auto;      /* DELEGATED: op auto-detected → re-detect at replay */
	star_finder_params star_conf; /* DELEGATED: detection parameters to reproduce */
};

struct synthstar_data *new_synthstar_data(void);

/* Star-list <-> blob codec (Convention 2), shared with unpurple. */
gchar *synthstar_stars_to_blob(psf_star **stars, int nb_stars);
psf_star **synthstar_stars_from_blob(const char *blob, int *nb_out);

/* star_finder_params <-> kv codec (DELEGATED provenance), shared with unpurple.
 * synthstar_conf_from_kv returns FALSE if any required key is absent. */
void synthstar_conf_to_kv(GString *kv, const star_finder_params *sf);
gboolean synthstar_conf_from_kv(GHashTable *kv, star_finder_params *out);

/* Run @body with @conf installed as the global starfinder configuration,
 * restoring the user's afterwards.  The DELEGATED replay path uses it so a
 * delegated star op re-detects with the recorded parameters, not current
 * prefs.  Replay owns the exclusive processing slot, so the transient global
 * swap races no concurrent detection. */
#define WITH_STARFINDER_CONF(conf, body) do {                       \
		star_finder_params _saved = com.pref.starfinder_conf;       \
		com.pref.starfinder_conf = *(conf);                         \
		body;                                                       \
		com.pref.starfinder_conf = _saved;                          \
	} while (0)

/* generic_image_worker hooks and log_hooks */
int synthstar_image_hook(struct generic_img_args *args, fits *fit, int threads);
gchar *synthstar_log_hook(gpointer p, log_hook_detail detail);
int unclip_image_hook(struct generic_img_args *args, fits *fit, int threads);
gchar *unclip_log_hook(gpointer p, log_hook_detail detail);

void makeairy(float *psf, const int size, const float lum, const float xoff, const float yoff, const float wavelength, const float aperture, const float focal_length, const float pixel_scale, const float obstruction);

void makegaussian(float *psf, int size, float fwhm, float lum, float xoffset, float yoffset, float ratio, float angle);
void makemoffat(float *psf, const int size, const float fwhm, const float lum, const float xoff, const float yoff, const float beta, float ratio, float angle);
void makedisc(float *psf, int size, float radius, float lum, float xoffset, float yoffset);
int starcount(psf_star **stars);

#define SYNTHESIZE_GAUSSIAN 0
#define SYNTHESIZE_MOFFAT 1

#endif /* SRC_SYNTHSTAR_H_ */
