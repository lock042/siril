#ifndef SRC_GUI_ASINH_H_
#define SRC_GUI_ASINH_H_

#include <glib.h>

void apply_asinh_cancel(void);
void asinh_change_between_roi_and_image(void);

/* Open the dialog in amend mode for history record @record_id
 * (nde_editors registry, convergence C4). */
gboolean asinh_open_amend(gint64 record_id);

#endif /* SRC_GUI_ASINH_H_ */
