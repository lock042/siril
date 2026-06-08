#ifndef SRC_GUI_DIALOGS_H_
#define SRC_GUI_DIALOGS_H_
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui-gtk4/gui_state.h"

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

/* SirilFileChooser: thin synchronous wrapper around GtkFileDialog so the
 * existing call sites can keep their create → configure → run → read →
 * destroy flow.  GtkFileDialog itself is async-only; siril_dialog_run()
 * spins a local GMainLoop until the underlying _open/_save/_select_folder
 * finish callback fires.  Distinct `siril_fc_*` names avoid colliding
 * with the inline-button shims in utils.c (siril_file_chooser_*). */
typedef struct _SirilFileChooser SirilFileChooser;

SirilFileChooser* siril_fc_open(GtkWindow *parent, GtkFileChooserAction action);
SirilFileChooser* siril_fc_add (GtkWindow *parent, GtkFileChooserAction action);
SirilFileChooser* siril_fc_save(GtkWindow *parent, GtkFileChooserAction action);

/* Configuration (apply before siril_fc_run). */
void siril_fc_set_current_name        (SirilFileChooser *fc, const gchar *name);
void siril_fc_set_current_folder_path (SirilFileChooser *fc, const gchar *path);
void siril_fc_set_filename            (SirilFileChooser *fc, const gchar *path);
void siril_fc_set_select_multiple     (SirilFileChooser *fc, gboolean multi);
void siril_fc_add_filter              (SirilFileChooser *fc, GtkFileFilter *filter, gboolean set_default);
void siril_fc_add_filter_pattern      (SirilFileChooser *fc, const gchar *title, const gchar *pattern, gboolean set_default);
GtkFileFilter *siril_fc_get_filter    (SirilFileChooser *fc);

/* Run synchronously.  Returns GTK_RESPONSE_ACCEPT, GTK_RESPONSE_CANCEL or
 * GTK_RESPONSE_NONE (on error other than dismissal). */
gint siril_fc_run(SirilFileChooser *fc);

/* Read the result after a successful run. */
gchar  *siril_fc_get_filename (SirilFileChooser *fc);
GSList *siril_fc_get_filenames(SirilFileChooser *fc);

void siril_fc_destroy(SirilFileChooser *fc);

int number_of_dialogs();

#endif /* SRC_GUI_DIALOGS_H_ */
