#ifndef SRC_ALGOS_GEOMETRY_H_
#define SRC_ALGOS_GEOMETRY_H_

#include "core/siril.h"

/* Parameter structs for single image geometry operations */

struct binning_args {
	destructor destroy_fn;  // Must be first member
	int factor;
	gboolean mean;
};

struct crop_args {
	destructor destroy_fn;  // Must be first member
	rectangle area;
};

struct mirror_args {
	destructor destroy_fn;  // Must be first member
	gboolean x_axis;  // TRUE for mirrorx, FALSE for mirrory
};

struct resample_args {
	destructor destroy_fn;  // Must be first member
	int toX;
	int toY;
	opencv_interpolation interpolation;
	gboolean clamp;
};

struct rotation_args {
	destructor destroy_fn;  // Must be first member
	rectangle area;
	double angle;
	int interpolation;
	int cropped;
	gboolean clamp;
};

/* Allocator and destructor functions */
struct binning_args *new_binning_args();
void free_binning_args(void *args);

struct crop_args *new_crop_args();
void free_crop_args(void *args);

struct mirror_args *new_mirror_args();
void free_mirror_args(void *args);

struct resample_args *new_resample_args();
void free_resample_args(void *args);

struct rotation_args *new_rotation_args();
void free_rotation_args(void *args);

/* Mask hooks for generic_image_worker */
int binning_mask_hook(struct generic_img_args *args);
int crop_mask_hook(struct generic_img_args *args);
int resample_mask_hook(struct generic_img_args *args);
int rotation_mask_hook(struct generic_img_args *args);
int mirrorx_mask_hook(struct generic_img_args *args);
int mirrory_mask_hook(struct generic_img_args *args);

/* Log hooks for generic_image_worker */
gchar *binning_log_hook(gpointer p, log_hook_detail detail);
gchar *crop_log_hook(gpointer p, log_hook_detail detail);
gchar *resample_log_hook(gpointer p, log_hook_detail detail);
gchar *rotation_log_hook(gpointer p, log_hook_detail detail);
// No mirror_log_hook as the description provides all that's needed

/* Image processing hooks for generic_image_worker */
int binning_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
int crop_image_hook_single(struct generic_img_args *args, fits *fit, int nb_threads);
int mirrorx_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
int mirrory_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
int resample_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
int rotation_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

/* Legacy functions - still used for sequences */
int fits_binning(fits *fit, int factor, gboolean mean);
int crop(fits *fit, rectangle *bounds);
void mirrorx(fits *fit, gboolean verbose);
void mirrory(fits *fit, gboolean verbose);
int verbose_resize_gaussian(fits *image, int toX, int toY, opencv_interpolation interpolation, gboolean clamp);
int verbose_rotate_fast(fits *image, int angle);
int verbose_rotate_image(fits *image, rectangle area, double angle, int interpolation, int cropped, gboolean clamp);

gboolean crop_gui_updates(gpointer user); // used by the command

/* crop sequence data from GUI */
struct crop_sequence_data {
	sequence *seq;
	rectangle area;
	char *prefix;
	int retvalue;
};

struct scale_sequence_data {
	sequence *seq;
	char *prefix;
	double scale;
	opencv_interpolation interpolation;
	gboolean clamp;
	int retvalue;
};

int fits_binning(fits *fit, int factor, gboolean mean);

int verbose_resize_gaussian(fits *image, int toX, int toY, opencv_interpolation interpolation, gboolean clamp);

int verbose_rotate_image(fits *, rectangle, double, int, int, gboolean);
int verbose_rotate_fast(fits *image, int angle);

void mirrorx(fits *fit, gboolean verbose);
void mirrory(fits *fit, gboolean verbose);

int crop(fits *fit, rectangle *bounds);

gpointer crop_sequence(struct crop_sequence_data *crop_sequence_data);

gpointer scale_sequence(struct scale_sequence_data *scale_sequence_data);

int eqcrop(double ra1, double dec1, double ra2, double dec2, int margin_px, double margin_asec, int minsize, fits *fit);

const char *interp_to_str(int interpolation);

#endif /* SRC_ALGOS_GEOMETRY_H_ */
