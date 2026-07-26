#ifndef SRC_GUI_CLAHE_H_
#define SRC_GUI_CLAHE_H_

#include <glib.h>

void apply_clahe_cancel(void);

/* Open the dialog in amend mode for history record @record_id
 * (nde_editors registry, convergence C5b). */
gboolean clahe_open_amend(gint64 record_id);

#endif /* SRC_GUI_CLAHE_H_ */
