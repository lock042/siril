#ifndef COSMETIC_CORRECTION_H_
#define COSMETIC_CORRECTION_H_

#include "core/siril.h"

typedef struct deviant_struct deviant_pixel;

/* Cosmetic data from GUI */
struct cosmetic_data {
	destructor destroy_fn;
	fits *fit;
	sequence *seq;
	double sigma[2];
	double amount;
	gboolean is_cfa;
	threading_type threading;
	char *seqEntry;
};

/* structure for cosme command */
struct cosme_data {
	destructor destroy_fn;  // Must be first member
	fits *fit;
	sequence *seq;
	int is_cfa;
	GFile *file;
	char *prefix;
};

typedef enum {
	COLD_PIXEL, HOT_PIXEL
} typeOfDeviant;

struct deviant_struct {
	point p;
	typeOfDeviant type;
};

deviant_pixel *find_deviant_pixels(fits *fit, const double sig[2], long *icold, long *ihot, gboolean eval_only);
void apply_cosmetic_to_sequence(struct cosmetic_data *cosme_args);
int apply_cosme_to_image(fits *fit, GFile *file, int is_cfa);
void apply_cosme_to_sequence(struct cosme_data *cosme_args);
gpointer autoDetectThreaded(gpointer p);
int cosmeticCorrection(fits *fit, deviant_pixel *dev, int size, gboolean is_CFA);
int cosmeticCorrOneLine(fits *fit, deviant_pixel dev, gboolean is_cfa);
int cosmeticCorrOnePoint(fits *fit, deviant_pixel dev, gboolean is_cfa);

int denoise_hook_cosmetic(fits *fit);

struct cosme_data *new_cosme_data();
void free_cosme_data(void *args);
struct cosmetic_data *new_cosmetic_data();
void free_cosmetic_data(void *args);
gchar* cosmetic_log_hook(gpointer p, log_hook_detail detail);

/* Image processing hook for generic_image_worker */
int cosme_image_hook_generic(struct generic_img_args *args, fits *fit, int nb_threads);

/* Idle function for generic_image_worker */
gboolean cosme_idle(gpointer p);

#endif /* COSMETIC_CORRECTION_H_ */
