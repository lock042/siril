#ifndef SRC_ALGOS_GEOMETRY_H_
#define SRC_ALGOS_GEOMETRY_H_

#include "core/siril.h"

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

int eqcrop(double ra1, double dec1, double ra2, double dec2, int margin_px, double margin_asec, int minsize, fits *fit, int *newx, int *newy);

const char *interp_to_str(opencv_interpolation interpolation);

#endif /* SRC_ALGOS_GEOMETRY_H_ */
