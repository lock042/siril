#ifndef SRC_GUI_SATURATION_H_
#define SRC_GUI_SATURATION_H_

#include <glib.h>

void apply_satu_cancel(void);
void satu_change_between_roi_and_image(void);

/* Open the dialog in amend mode for history record @record_id
 * (nde_editors registry, convergence C5). */
void satu_open_amend(gint64 record_id);

#endif /* SRC_GUI_SATURATION_H_ */
