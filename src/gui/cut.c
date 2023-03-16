#include <math.h>
#include "core/siril.h"
#include "gui/image_interactions.h"
#include "core/processing.h"

float interpf(fits* fit, float x, float y, int chan) {
	if (chan >= fit->naxes[2])
		return -9999.f;
	int width = fit->rx;
	int x0 = floor(x);
	int x1 = ceil(x);
	int y0 = floor(y);
	int y1 = ceil(y);
	float val00 = fit->fpdata[chan][x0 + y0 * width];
	float val01 = fit->fpdata[chan][x1 + y0 * width];
	float val10 = fit->fpdata[chan][x0 + y1 * width];
	float val11 = fit->fpdata[chan][x1 + y1 * width];
	float interp1 = (x - x0) * val00 + (x1 - x) * val01;
	float interp2 = (x - x0) * val10 + (x1 - x) * val11;
	float interp = (y - y0) * interp1 + (y1 - y) * interp2;
	return interp;
}

float interpw(fits* fit, float x, float y, int chan) {
	if (chan >= fit->naxes[2])
		return -9999.f;
	int width = fit->rx;
	int x0 = floor(x);
	int x1 = ceil(x);
	int y0 = floor(y);
	int y1 = ceil(y);
	float val00 = (float) fit->pdata[chan][x0 + y0 * width];
	float val01 = (float) fit->pdata[chan][x1 + y0 * width];
	float val10 = (float) fit->pdata[chan][x0 + y1 * width];
	float val11 = (float) fit->pdata[chan][x1 + y1 * width];
	float interp1 = (x - x0) * val00 + (x1 - x) * val01;
	float interp2 = (x - x0) * val10 + (x1 - x) * val11;
	float interp = (y - y0) * interp1 + (y1 - y) * interp2;
	return interp;
}

float interp(fits *fit, float x, float y, int chan) {
	switch (fit->type) {
		case DATA_FLOAT:
			return interpf(fit, x, y, chan);
			break;
		case DATA_USHORT:
			return interpw(fit, x, y, chan);
			break;
		default:
			return -9999.f;
			break;
	}
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
	float *r = NULL, *g = NULL, *b = NULL;
	float length = sqrtf(delta.x * delta.x + delta.y * delta.y);
	if (length < 1.f) {
		retval = 1;
		goto END;
	}
	int nbr_points = (int) length;
	float point_spacing = length / nbr_points;
	float point_spacing_x = delta.x / nbr_points;
	float point_spacing_y = delta.y / nbr_points;

	r = malloc(nbr_points * sizeof(float));
	if (gfit.naxes[2] > 1) {
		g = malloc(nbr_points * sizeof(float));
		b = malloc(nbr_points * sizeof(float));
	}
	float *x = malloc(nbr_points * sizeof(float));
	for (int i = 0 ; i < nbr_points ; i++) {
		x[i] = start.x + i * point_spacing;
		r[i] = interp(&gfit, (float) (start.x + ((delta.x / nbr_points) * point_spacing)), (float) (start.y + ((delta.y / nbr_points) * point_spacing)), 0);
		if (gfit.naxes[2] > 1) {
			g[i] = interp(&gfit, (float) (start.x + ((delta.x / nbr_points) * point_spacing)), (float) (start.y + ((delta.y / nbr_points) * point_spacing)), 1);
			b[i] = interp(&gfit, (float) (start.x + ((delta.x / nbr_points) * point_spacing)), (float) (start.y + ((delta.y / nbr_points) * point_spacing)), 2);
		}
	}
	// TODO: Plot r{,g,b} against x using kplot

END:
	// Clean up
	com.cut_point.x = -1;
	com.cut_point.y = -1;
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}
