#ifndef SRC_GUI_SAVE_DIALOG_H_
#include <gtk/gtk.h>
#define SRC_GUI_SAVE_DIALOG_H_

/* Savedialog data from GUI */
struct savedial_data {
	GtkEntry *entry;
	gint bitspersamples;
	gboolean tiff_compression;
	gchar *description;
	gchar *copyright;
	gint quality;
	gint jxl_effort;
	gdouble jxl_quality;
	gboolean jxl_force_8bit;
	gboolean lossless;
	const gchar *filename;
	int bitpix;
	gboolean update_hilo;
	gboolean checksum;
	int retval;
};

enum {
	PAGE_TIFF, PAGE_JPG, PAGE_FITS, PAGE_JXL, PAGE_MISC
};

void on_header_save_as_button_clicked();
void on_header_snapshot_button_clicked(gboolean clipboard);
void on_header_save_button_clicked();

/* TIFF UI helpers (moved from io/image_formats_libraries.c) */
gboolean get_tiff_compression(void);
void get_tif_data_from_ui(fits *fit, gchar **description, gchar **copyright);

#endif /* SRC_GUI_SAVE_DIALOG_H_ */
