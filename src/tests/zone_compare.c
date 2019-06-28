#include <stdio.h>
#include <stdlib.h>
#include "../core/siril.h"
#include "../core/proto.h"
#include "../io/sequence.h"
#include "../io/single_image.h"

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
double zoom = 1.0;
stacking_zone zone = { .centre = {.x = -1.0 } };


void remap(int image);

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
		surface_stride[image] = cairo_format_stride_for_width(
				CAIRO_FORMAT_RGB24, fit[image].rx);
		graybuf[image] = realloc(graybuf[image],
				surface_stride[image] * fit[image].ry * sizeof(guchar));
		if (surface[image])
			cairo_surface_destroy(surface[image]);
		surface[image] = cairo_image_surface_create_for_data(graybuf[image],
				CAIRO_FORMAT_RGB24, fit->rx, fit->ry,
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
	if (image < 0 || image >= 2) return TRUE;
	remap(image);

	// autofit zoom
	double wtmp = (double)window_width / (double)fit[image].rx;
	double htmp = (double)window_height / (double)fit[image].ry;
	zoom = min(wtmp, htmp);
	cairo_scale(cr, zoom, zoom);

	cairo_set_source_surface(cr, surface[image], 0, 0);
	cairo_paint(cr);

	if (image == 0 && zone.centre.x >= 0.0) {
		// draw the zone
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_line_width(cr, 1.5);
		cairo_set_source_rgba(cr, 0.0, 1.0, 1.0, 1.0);	// cyan
		cairo_rectangle(cr, zone.centre.x - zone.half_side,
				zone.centre.y - zone.half_side,
				2.0 * zone.half_side, 2.0 * zone.half_side);
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
	if (seq_read_frame(seq, image, &fit[1]))
		fprintf(stderr, "could not open image %d from sequence\n", image);
	image_find_minmax(&fit[1]);
	changed = 1;
	remap(1);
	gtk_widget_queue_draw(widget_for_image(1));
}

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	zone.centre.x = event->x / zoom;
	zone.centre.y = event->y / zoom;
	GtkRange *range = GTK_RANGE(gtk_builder_get_object(builder, "scale_zonesize"));
	zone.half_side = round_to_int(gtk_range_get_value(range));
	zone.mpregparam = NULL;
	gtk_widget_queue_draw(widget);
	return FALSE;
}

