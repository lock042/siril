#include <stdio.h>
#include <stdlib.h>
#include "../core/siril.h"
#include "../core/proto.h"

#define NB_DISPLAYS 3

GtkBuilder *builder = NULL;	// get widget references anywhere
fits fit[NB_DISPLAYS];
cairo_surface_t *surface[NB_DISPLAYS];
guchar *graybuf[NB_DISPLAYS];	// R=G=B 8bit version
int surface_stride[NB_DISPLAYS];
int changed = 1;

int main(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stderr, "Usage: %s reference_image tested_image\n", *argv);
		exit(1);
	}

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
	if (readfits(argv[2], &fit[1], NULL)) {
		fprintf(stderr, "%s not found\n", argv[2]);
		exit(1);
	}

	gtk_main();
	return 0;
}

/* image modification functions */


/* image matching functions */


/* image display functions */
void remap(int image) {
	if (changed) {
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

		int i;
		float pente = UCHAR_MAX_SINGLE / (float)(fit[image].maxi - fit[image].mini);
		for (i = 0; i < fit[image].rx * fit[image].ry; i++) {
			 graybuf[image][i] = round_to_BYTE((float)fit[image].data[i] * pente);
		}

		changed = 0;
	}
}

/* callbacks */
int image_for_widget(GtkWidget *widget) {
	if (widget == GTK_WIDGET(gtk_builder_get_object(builder, "drawing1")))
		return 0;
	if (widget == GTK_WIDGET(gtk_builder_get_object(builder, "drawing2")))
		return 1;
	if (widget == GTK_WIDGET(gtk_builder_get_object(builder, "drawing3")))
		return 2;
	fprintf(stderr, "unknown widget\n");
	return -1;
}

gboolean redraw_drawingarea1(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int window_width = gtk_widget_get_allocated_width(widget);
	int window_height = gtk_widget_get_allocated_height(widget);

	int image = image_for_widget(widget);
	remap(image);
	//cairo_scale(cr, zoom, zoom);
	cairo_set_source_surface(cr, com.surface[image], 0, 0);

	return FALSE;
}
