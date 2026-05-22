#ifndef SRC_GUI_SAVE_DIALOG_H_
#include <gtk/gtk.h>
#define SRC_GUI_SAVE_DIALOG_H_

/* User's choice from the multi-layer FLIS save confirmation dialog.
 * The dialog is shown on the main thread (GTK can only be called from
 * the GUI thread); the result is stashed in savedial_data and read by
 * the worker thread when it dispatches the actual save.  Default
 * FLIS_SAVE_AUTOMATIC means "no choice needed" — the worker decides
 * based on data state (single-layer FLIS → save_flis; plain FITS →
 * savefits; multi-layer FLIS reaching the worker with AUTOMATIC is
 * a bug — the dialog should have been shown by the caller). */
typedef enum {
	FLIS_SAVE_AUTOMATIC = 0,  /* not multi-layer or no dialog needed */
	FLIS_SAVE_AS_FLIS,        /* preserve layers via save_flis */
	FLIS_SAVE_FLATTEN,        /* flis_flatten_all then savefits */
} flis_save_choice_t;

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
	flis_save_choice_t flis_save_choice;  /* set on main thread before worker starts */
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
