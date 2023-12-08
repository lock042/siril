#ifndef SRC_GUI_SAVE_DIALOG_H_
#define SRC_GUI_SAVE_DIALOG_H_

/* Savedialog data from GUI */
struct savedial_data {
	GtkEntry *entry;
	gint bitspersamples;
	gboolean tiff_compression;
	char *description;
	char *copyright;
	gint quality;
	gint jxl_effort;
	gdouble heif_quality;
	gdouble jxl_quality;
	gboolean jxl_force_8bit;
	gboolean lossless;
	int avif_bitdepth;
	const gchar *filename;
	int bitpix;
	gboolean update_hilo;
	int retval;
};

enum {
	PAGE_TIFF, PAGE_JPG, PAGE_FITS, PAGE_JXL, PAGE_AVIF, PAGE_MISC
};

void on_header_save_as_button_clicked();
void on_header_snapshot_button_clicked(gboolean clipboard);
void on_header_save_button_clicked();

#endif /* SRC_GUI_SAVE_DIALOG_H_ */
