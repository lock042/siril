#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../core/siril.h"
#include "../core/proto.h"
#include "../io/sequence.h"
#include "../io/single_image.h"
#include "../planetary/planetary.h"
#include "../gui/zones.h"
#include "../opencv/opencv.h"
#include "../algos/statistics.h"
#include "../registration/registration.h"

#define NB_DISPLAYS 3
#define MAX_RADIUS 50 // in pixels, maximum local shift
#define DEV_STRIDE 1 // use every DEV_STRIDE pixels to compute the deviation

/* global variables for siril, not used here, but required for linking against libsiril.a */
cominfo com;	// the main data struct
fits gfit;	// currently loaded image
//GtkBuilder *builder;	// get widget references anywhere
char **supported_extensions;
/*************************************************************/

GtkBuilder *builder = NULL;	// get widget references anywhere
fits fit[NB_DISPLAYS];
float *ref, *im; // gaussian-filtered, normalized data of fit[0] and fit[1]
sequence *seq = NULL;
cairo_surface_t *surface[NB_DISPLAYS] = { 0 };
guchar *graybuf[NB_DISPLAYS];	// R=G=B 8bit version
int surface_stride[NB_DISPLAYS];
int changed = 1;
double zoom[NB_DISPLAYS] = { 1.0 };
stacking_zone zone = { .centre = {.x = -1.0 } };
enum { MODE_DIRECT, MODE_SQUARED, MODE_MEAN } mode = MODE_DIRECT;
int adaptive_rendering = 1;
double local_shift_x = 0.0, local_shift_y = 0.0;
int kernel_size = 13;

#define LAYER 0

void filter_and_normalize_fit(fits *fit, float **buf, gboolean flip);
void remap(int image);
void update_comparison();
void update_shift_display();

int main(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stderr, "Usage: %s reference_image tested_sequence\n", *argv);
		exit(1);
	}

	gtk_init(&argc, &argv);
	builder = gtk_builder_new();
	GError *err = NULL;
	if (!gtk_builder_add_from_file (builder, "zone_compare.glade", &err)) {
		fprintf(stderr, "could not open glade: %s\n", err->message);
		g_error_free(err);
		exit(1);
	}
	gtk_builder_connect_signals (builder, NULL);


	/* open files */
	if (readfits(argv[1], &fit[0], NULL)) {
		fprintf(stderr, "%s not found\n", argv[1]);
		exit(1);
	}
	image_find_minmax(&fit[0]);
	remap(0);
	changed = 1;

	if (!(seq = readseqfile(argv[2])))
		exit(1);
	if (seq_read_frame(seq, 0, &fit[1]))
		exit(1);
	if (fit[0].rx != fit[1].rx || fit[0].ry != fit[1].ry) {
		fprintf(stderr, "reference frame and images of the sequence must be of the same size\n");
		exit(1);
	}
	filter_and_normalize_fit(&fit[0], &ref, FALSE);
	filter_and_normalize_fit(&fit[1], &im, FALSE);
	seq->selnum = 0;
	GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2"));
	gtk_adjustment_set_upper(adj, (double)(seq->number-1));
	image_find_minmax(&fit[1]);
	update_shift_display();
	remap(1);

	gtk_main();
	return 0;
}

void filter_and_normalize_fit(fits *fit, float **buf, gboolean flip) {
	unsigned int nbdata = fit->rx * fit->ry;
	WORD min, max;
	WORD *gauss = malloc(nbdata * sizeof(WORD));
	*buf = realloc(*buf, nbdata * sizeof(float));
	// TODO: RGB to monochrome
	cvGaussian(fit->data, fit->rx, fit->ry, kernel_size, gauss);
	if (flip)
		cvFlip_siril(gauss, fit->rx, fit->ry);
	siril_stats_ushort_minmax(&min, &max, gauss, nbdata);
	normalize_data(gauss, nbdata, min, max, *buf);
	free(gauss);
}

/* image display functions */
void remap(int image) {
	if (changed) {
		fprintf(stdout, "remap %d\n", image);
		if (image > 2) return;
		surface_stride[image] = cairo_format_stride_for_width(
				CAIRO_FORMAT_RGB24, fit[image].rx);
		graybuf[image] = realloc(graybuf[image],
				surface_stride[image] * fit[image].ry * sizeof(guchar));
		if (surface[image])
			cairo_surface_destroy(surface[image]);
		surface[image] = cairo_image_surface_create_for_data(graybuf[image],
				CAIRO_FORMAT_RGB24, fit[image].rx, fit[image].ry,
				surface_stride[image]);
		if (cairo_surface_status(surface[image])
				!= CAIRO_STATUS_SUCCESS) {
			fprintf(stderr, "Error creating the Cairo image surface\n");
			cairo_surface_destroy(surface[image]);
			surface[image] = NULL;
			return;
		}

		// screen transfer function, linear
		int y;
		double pente = UCHAR_MAX_DOUBLE / (fit[image].maxi - fit[image].mini);
		//fprintf(stdout, "remap pente: %g, mini: %g, maxi: %g\n", pente, fit[image].mini, fit[image].maxi);
		for (y = 0; y < fit[image].ry; y++) {
			int x;
			for (x = 0; x < fit[image].rx; x++) {
				guint dst_index = ((fit[image].ry - 1 - y) * fit[image].rx + x) * 4;
				BYTE pix = round_to_BYTE((fit[image].data[y*fit[image].rx + x] - fit[image].mini) * pente);
				graybuf[image][dst_index++] = pix;
				graybuf[image][dst_index++] = pix;
				graybuf[image][dst_index++] = pix;
			}
		}

		changed = 0;
	}
}

/* callbacks */
int image_for_widget(GtkWidget *widget) {
	if (widget == GTK_WIDGET(gtk_builder_get_object(builder, "draw1")))
		return 0;
	if (widget == GTK_WIDGET(gtk_builder_get_object(builder, "draw2")))
		return 1;
	if (widget == GTK_WIDGET(gtk_builder_get_object(builder, "draw3")))
		return 2;
	fprintf(stderr, "unknown widget\n");
	return -1;
}

GtkWidget *widget_for_image(int image) {
	char wname[10];
	sprintf(wname, "draw%d", image+1);
	return GTK_WIDGET(gtk_builder_get_object(builder, wname));
}

gboolean redraw_drawingarea(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int window_width = gtk_widget_get_allocated_width(widget);
	int window_height = gtk_widget_get_allocated_height(widget);

	int image = image_for_widget(widget);
	if (image < 0 || image > 2) return TRUE;
	if (zone.half_side <= 0 && image == 2) return TRUE;
	remap(image); // does nothing if !changed

	// autofit zoom
	double wtmp = (double)window_width / (double)fit[image].rx;
	double htmp = (double)window_height / (double)fit[image].ry;
	zoom[image] = min(wtmp, htmp);
	cairo_scale(cr, zoom[image], zoom[image]);

	cairo_set_source_surface(cr, surface[image], 0, 0);
	cairo_paint(cr);

	if (image == 0 && zone.centre.x >= 0.0) {
		double side = (double)get_side(&zone);
		// draw the zone
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_line_width(cr, 1.5);
		cairo_set_source_rgba(cr, 0.0, 1.0, 1.0, 1.0);	// cyan
		cairo_rectangle(cr, zone.centre.x - zone.half_side,
				zone.centre.y - zone.half_side,
				side, side);
		cairo_stroke(cr);

		// draw points at the centre
		cairo_set_line_width(cr, 4.0);
		cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, 1.0);	// red 
		cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
		cairo_move_to(cr, zone.centre.x, zone.centre.y);
		cairo_close_path(cr);
		cairo_stroke(cr);
		cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
	}

	return FALSE;
}

void update_shift_display() {
	char shifts[50];
	GtkLabel *label = GTK_LABEL(gtk_builder_get_object(builder, "labelgshifts"));
	regdata *regparam = seq->regparam[0];
	if (!regparam || seq->selnum < 0)
		gtk_label_set_text(label, "not available");
	else {
		sprintf(shifts, "%.3f, %.3f",
				regparam[seq->selnum].shiftx, regparam[seq->selnum].shifty);
		gtk_label_set_text(label, shifts);
	}
}

void update_error_display(double err) {
	char error[50];
	GtkLabel *label = GTK_LABEL(gtk_builder_get_object(builder, "labelerr"));
	sprintf(error, "%.4f", err);
	gtk_label_set_text(label, error);
}

void on_scale_seqnumber_value_changed(GtkRange *range, gpointer user_data) {
	int image = round_to_int(gtk_range_get_value(range));		
	clearfits(&fit[1]);
	if (seq_read_frame(seq, image, &fit[1])) {
		fprintf(stderr, "could not open image %d from sequence\n", image);
		return;
	}
	seq->selnum = image;
	filter_and_normalize_fit(&fit[1], &im, FALSE);
	image_find_minmax(&fit[1]);
	changed = 1;
	remap(1);
	gtk_widget_queue_draw(widget_for_image(1));
	update_shift_display();
	update_comparison();
}

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	int image = image_for_widget(widget);
	zone.centre.x = event->x / zoom[image];
	zone.centre.y = event->y / zoom[image];
	GtkRange *range = GTK_RANGE(gtk_builder_get_object(builder, "scale_zonesize"));
	zone.half_side = round_to_int(gtk_range_get_value(range)) / 2;
	zone.mpregparam = NULL;
	gtk_widget_queue_draw(widget);	// display the selection
	fprintf(stdout, "zone at %d, %d, half_side: %d\n", zone.centre.x, zone.centre.y, zone.half_side);

	update_comparison();
	return FALSE;
}

void update_comparison() {
	if (zone.half_side <= 0)
		return;
	int side = get_side(&zone), nb_pix = side * side;
	float *ref_zone = malloc(nb_pix * sizeof(float));
	if (copy_image_buffer_zone_to_buffer_float(ref, fit[0].rx, fit[0].ry, &zone, ref_zone)) {
		fprintf(stderr, "zone is outside image\n");
		free(ref_zone);
		return;
	}

	fprintf(stdout, "ref_zone: %g %g %g %g %g\n", ref_zone[0],
			ref_zone[1], ref_zone[2], ref_zone[3], ref_zone[4]); 
	
	// get the area from the sequence image, with regdata
	float *im_zone = malloc(nb_pix * sizeof(float));
	regdata *regparam = seq->regparam[0];
	if (regparam) {
		stacking_zone shifted_zone = { .centre =
			{ .x = round_to_int(zone.centre.x - regparam[seq->selnum].shiftx + local_shift_x),
				.y = round_to_int(zone.centre.y + regparam[seq->selnum].shifty + local_shift_y) },
			.half_side = zone.half_side };
		if (copy_image_buffer_zone_to_buffer_float(im, fit[1].rx, fit[1].ry, &shifted_zone, im_zone)) {
			fprintf(stderr, "zone is outside image\n");
			free(ref_zone); free(im_zone);
			return;
		}
		fprintf(stdout, "comparing zones centered on %d, %d of reference to %d, %d of image\n",
				zone.centre.x, zone.centre.y,
				shifted_zone.centre.x, shifted_zone.centre.y);
	}
	else {
		if (copy_image_buffer_zone_to_buffer_float(im, fit[1].rx, fit[1].ry, &zone, im_zone)) {
			fprintf(stderr, "zone is outside image\n");
			free(ref_zone); free(im_zone);
			return;
		}
		fprintf(stdout, "comparing zones centered on %d, %d of reference to the same of image\n",
				zone.centre.x, zone.centre.y);
	}

	// compute the displayed patch
	double mini = DBL_MAX, maxi = DBL_MIN, sum = 0.0;
	double *diff = malloc(nb_pix * sizeof(double));
	int i;
	for (i = 0; i < nb_pix; i++) {
		switch (mode) {
		case MODE_DIRECT:
		case MODE_MEAN:
			diff[i] = ref_zone[i] - im_zone[i];
			sum += fabs(diff[i]);
			break;
		case MODE_SQUARED:
			diff[i] = (ref_zone[i] - im_zone[i]) * (ref_zone[i] - im_zone[i]);
			sum += diff[i];
			break;
		}
		if (diff[i] < mini)
			mini = diff[i];
		if (diff[i] > maxi)
			maxi = diff[i];
	}
	if (mode == MODE_MEAN) {
		double add = 0.0;
		for (i = 0; i < nb_pix; i++)
			add += fabs(diff[i]);
		double mean = add / nb_pix;
		double max = 0.0;
		for (i = 0; i < nb_pix; i++) {
			double dev = (diff[i] - mean) * (diff[i] - mean);
			if (dev > max)
				max = dev;
		}
		sum = max;	// for display
	}

	free(im_zone);
	free(ref_zone);
	fprintf(stdout, "mini: %g, maxi: %g\n", mini, maxi);
	update_error_display(sum);

	// transfer to ushort range for display
	if (!adaptive_rendering) {
		// force the range
		if (mode == MODE_DIRECT) {
			mini = -0.1; maxi = 0.1;
		} else {
			mini = 0.0; maxi = 0.006;
		}
	}
	double pente = USHRT_MAX_DOUBLE / (maxi - mini);
	WORD *result = malloc(nb_pix * sizeof(WORD));
	for (i = 0; i < nb_pix; i++) {
		result[i] = round_to_WORD((diff[i] - mini) * pente);
	}
	free(diff);

	if (fit[2].rx > 0)
		clearfits(&fit[2]);
	fits *newfit = &fit[2];
	new_fit_image(&newfit, side, side, 1, result);
	newfit->maxi = 0;
	newfit->mini = 65535;

	changed = 1;
	remap(2);
	gtk_widget_queue_draw(widget_for_image(2));
}

void on_radiodirect_toggled(GtkRadioButton *button, gpointer user_data) {
	mode = MODE_DIRECT;
	update_comparison();
}

void on_radiosquared_toggled(GtkRadioButton *button, gpointer user_data) {
	mode = MODE_SQUARED;
	update_comparison();
}

void on_radiovariance_toggled(GtkRadioButton *button, gpointer user_data) {
	mode = MODE_MEAN;
	update_comparison();
}

void on_localXspin_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	local_shift_x = gtk_spin_button_get_value(spinbutton);
	update_comparison();
}

void on_localYspin_value_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	local_shift_y = gtk_spin_button_get_value(spinbutton);
	update_comparison();
}

static void zone_to_rectangle(stacking_zone *zone, rectangle *rectangle) {
	int side = get_side(zone);
	rectangle->x = round_to_int(zone->centre.x - zone->half_side);
	rectangle->y = round_to_int(zone->centre.y - zone->half_side);
	rectangle->w = side;
	rectangle->h = side;
}

void on_computebutton_clicked(GtkButton *button, gpointer user_data) {
	if (zone.half_side <= 0)
		return;

	float *ref_flip = NULL, *im_flip = NULL;
	filter_and_normalize_fit(&fit[0], &ref_flip, TRUE);
	filter_and_normalize_fit(&fit[1], &im_flip, TRUE);

	rectangle ref_area, im_area;
	zone_to_rectangle(&zone, &ref_area);
	zone_to_rectangle(&zone, &im_area);
	regdata *regparam = seq->regparam[0];
	if (regparam) {
		stacking_zone shifted_zone = { .centre =
			{ .x = round_to_int(zone.centre.x - regparam[seq->selnum].shiftx),
				.y = round_to_int(zone.centre.y + regparam[seq->selnum].shifty) },
			.half_side = zone.half_side };
		zone_to_rectangle(&shifted_zone, &im_area);
	}

	int ref_i = ref_area.y * fit[0].rx + ref_area.x;
	fprintf(stdout, "ref_flip: %g %g %g %g %g\n", ref_flip[ref_i], ref_flip[ref_i+1],
			ref_flip[ref_i+2], ref_flip[ref_i+3], ref_flip[ref_i+4]);
	int shiftx = 0, shifty = 0;
	int error = search_local_match_gradient_float(ref_flip, im_flip,
			fit[0].rx, fit[0].ry, &ref_area, &im_area,
			MAX_RADIUS, DEV_STRIDE,
			&shiftx, &shifty);
	fprintf(stdout, "Found: %d, %d\n", shiftx, shifty);
	if (regparam) {
		// result is relative to the provided area, but our local shift is added
		// to the area
/*		local_shift_x = -regparam[seq->selnum].shiftx - (double)shiftx;
		local_shift_y = regparam[seq->selnum].shifty + (double)shifty;
	} else {*/
		local_shift_x = (double)shiftx;
		local_shift_y = -(double)shifty;
	}

	if (error) {
		fprintf(stderr, "FINDING THE LOCAL SHIFT FAILED\n");
	} else {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "localXspin")), local_shift_x);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "localYspin")), local_shift_y);
	}
}

void on_spingauss_value_changed(GtkSpinButton *spin, gpointer user_data) {
	kernel_size = gtk_spin_button_get_value_as_int(spin);
	fprintf(stdout, "new kernel size: %d\n", kernel_size);
	update_comparison();
}

void on_adaptive_render_checkbox_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	adaptive_rendering = gtk_toggle_button_get_active(togglebutton);
	update_comparison();
}

