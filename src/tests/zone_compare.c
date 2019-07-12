#include <stdio.h>
#include <stdlib.h>
#include "../core/siril.h"
#include "../core/proto.h"
#include "../io/sequence.h"
#include "../io/single_image.h"
#include "../planetary/planetary.h"
#include "../gui/zones.h"

#define NB_DISPLAYS 3

/* global variables for siril, not used here, but required for linking against libsiril.a */
cominfo com;	// the main data struct
fits gfit;	// currently loaded image
//GtkBuilder *builder;	// get widget references anywhere
char **supported_extensions;
/*************************************************************/

GtkBuilder *builder = NULL;	// get widget references anywhere
fits fit[NB_DISPLAYS];
sequence *seq = NULL;
cairo_surface_t *surface[NB_DISPLAYS];
guchar *graybuf[NB_DISPLAYS];	// R=G=B 8bit version
int surface_stride[NB_DISPLAYS];
int changed = 1;
double zoom[NB_DISPLAYS] = { 1.0 };
stacking_zone zone = { .centre = {.x = -1.0 } };

#define LAYER 0

void remap(int image);
void update_comparison();

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
	GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2"));
	gtk_adjustment_set_upper(adj, (double)(seq->number-1));
	image_find_minmax(&fit[1]);
	remap(1);

	gtk_main();
	return 0;
}

/* image modification functions */


/* image matching functions */


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

		int y;
		float pente = UCHAR_MAX_SINGLE / (float)(fit[image].maxi - fit[image].mini);
		//float pente = 1.0f;
		for (y = 0; y < fit[image].ry; y++) {
			int x;
			for (x = 0; x < fit[image].rx; x++) {
				guint dst_index = ((fit[image].ry - 1 - y) * fit[image].rx + x) * 4;
				BYTE pix = round_to_BYTE((float)fit[image].data[y*fit[image].rx + x] * pente);
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

void on_scale_seqnumber_value_changed(GtkRange *range, gpointer user_data) {
	int image = round_to_int(gtk_range_get_value(range));		
	clearfits(&fit[1]);
	if (seq_read_frame(seq, image, &fit[1])) {
		fprintf(stderr, "could not open image %d from sequence\n", image);
		return;
	}
	seq->number = image;
	image_find_minmax(&fit[1]);
	changed = 1;
	remap(1);
	gtk_widget_queue_draw(widget_for_image(1));
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

	int i, side = get_side(&zone), nb_pix = side * side;

	WORD *ref, *im, *diff;
	// get the area from the ref image
	ref = malloc(nb_pix * sizeof(WORD));
	copy_image_zone_to_buffer(&fit[0], &zone, ref, LAYER);
	
	// get the area from the sequence image, with regdata
	im = malloc(nb_pix * sizeof(WORD));
	regdata *regparam = seq->regparam[0];
	if (regparam) {
		stacking_zone shifted_zone = { .centre =
			{ .x = round_to_int(zone.centre.x - regparam[seq->number].shiftx),
				.y = round_to_int(zone.centre.y + regparam[seq->number].shifty) },
			.half_side = zone.half_side };
		copy_image_zone_to_buffer(&fit[1], &shifted_zone, im, LAYER);
	}
	else copy_image_zone_to_buffer(&fit[1], &zone, im, LAYER);
	// FIXME ^ apparently it's not working ^

	// compute the displayed patch
	diff = malloc(nb_pix * sizeof(WORD));
	WORD mini = 65535, maxi = 0;
	for (i = 0; i < nb_pix; i++) {
		// use the abs(ref_i - im_i)
		if (ref[i] < im[i])
			diff[i] = im[i] - ref[i];
		else diff[i] = ref[i] - im[i];
		if (diff[i] < mini)
			mini = diff[i];
		if (diff[i] > maxi)
			maxi = diff[i];
	}
	free(im);
	free(ref);

	if (fit[2].rx > 0)
		clearfits(&fit[2]);
	fits *newfit = &fit[2];
	new_fit_image(&newfit, side, side, 1, diff);
	newfit->maxi = maxi;
	newfit->mini = mini;

	changed = 1;
	remap(2);
	gtk_widget_queue_draw(widget_for_image(2));
}
