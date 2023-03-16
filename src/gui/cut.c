#include "core/siril.h"
#include "gui/image_interactions.h"
#include "core/processing.h"

float interp(fits* fit, float x, float y) {
}

gpointer cut_profile() {
	int retval = 0;

	pointi start, finish, delta;
	start.x = gui.start.x;
	start.y = gui.start.y;
	finish.x = com.cut_point.x;
	finish.y = com.cut_point.y;
	delta.x = finish.x - start.x;
	delta.y = finish.y - start.y;
	float length = sqrtf(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.f) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	float point_spacing = length / nbr_points;
	float point_spacing_x = delta.x / nbr_points;
	float point_spacing_y = delta_y / nbr_points;

	float *y = malloc(nbr_points * sizeof(float);
	float *x = malloc(nbr_points * sizeof(float);
	for (int i = 0 ; i < nbr_points ; i++) {
		x[i] = start.x + i * point_spacing;
		y[i] = interp(&gfit, (float) (start.x + ((delta.x / nbr_points) * point_spacing)), (float) (start.y + ((delta.y / nbr_points) * point_spacing)));
	}
	// TODO: Plot y = f(x) using kplot

END:
	// Clean up
	com.cut_point.x = -1;
	com.cut_point.y = -1;
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}
