#ifndef GUI_COLORS_H_
#define GUI_COLORS_H_

#include <glib.h>

void initialize_calibration_interface();
void negative_processing();

/* Open the CCM dialog in amend mode for history record @record_id
 * (nde_editors registry, convergence C5b). */
gboolean ccm_open_amend(gint64 record_id);

#endif /* !GUI_COLORS_H */
