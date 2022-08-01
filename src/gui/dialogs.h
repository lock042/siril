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
void siril_close_all_standard_dialogs();

#endif /* SRC_GUI_DIALOGS_H_ */
