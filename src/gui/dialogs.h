#ifndef SRC_GUI_DIALOGS_H_
#define SRC_GUI_DIALOGS_H_

#include "core/siril.h"
#include "core/proto.h"

typedef enum {
	NO_DIALOG = -1,
	INFORMATION_DIALOG,
	IMAGE_PROCESSING_DIALOG,
	SEARCH_ENTRY_DIALOG,
	OTHER_DIALOG
} DialogType;

struct _SirilDialogEntry {
	gchar *identifier;	// the identifier from the glade file if defined in it, or any if not
	GtkWidget *(*get_window)(void);	// if not defined in glade, this is a getter for the window
	DialogType type;
	gboolean has_preview;
	void (*apply_function)(void);	// if has_preview
};

void siril_open_dialog(gchar *id);
void siril_close_dialog(gchar *id);
void siril_close_preview_dialogs();
gboolean siril_widget_hide_on_delete(GtkWidget *widget);
gboolean is_a_dialog_opened();
gboolean is_an_image_processing_dialog_opened();
void mark_imgproc_dialog_closed();

void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title, const gchar *pattern, gboolean set_default);
SirilWidget* siril_file_chooser_open(GtkWindow *parent, GtkFileChooserAction action);
SirilWidget* siril_file_chooser_add(GtkWindow *parent, GtkFileChooserAction action);
SirilWidget* siril_file_chooser_save(GtkWindow *parent, GtkFileChooserAction action);
gint siril_dialog_run(SirilWidget *widgetdialog);
void siril_widget_destroy(SirilWidget *widgetdialog);

#endif /* SRC_GUI_DIALOGS_H_ */
