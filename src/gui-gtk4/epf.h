#ifndef SRC_GUI_BILAT_H_
#define SRC_GUI_BILAT_H_

#include <glib.h>

void epf_change_between_roi_and_image();
void apply_epf_cancel();

/* Open the dialog in amend mode for history record @record_id
 * (nde_editors registry, convergence C5b). */
gboolean epf_open_amend(gint64 record_id);

#endif /* SRC_GUI_ASINH_H_ */
