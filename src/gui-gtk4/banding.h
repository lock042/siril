#ifndef SRC_GUI_BANDING_H_
#define SRC_GUI_BANDING_H_

#include <glib.h>

/* Open the dialog in amend mode for history record @record_id
 * (nde_editors registry, convergence C5b). */
gboolean banding_open_amend(gint64 record_id);

#endif /* SRC_GUI_BANDING_H_ */
