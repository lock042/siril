#ifndef SRC_GUI_MEDIAN_H_
#define SRC_GUI_MEDIAN_H_

#include <glib.h>

void median_close(void);

/* Open the dialog in amend mode for history record @record_id
 * (nde_editors registry, convergence C5). */
void median_open_amend(gint64 record_id);

#endif /* SRC_GUI_MEDIAN_H_ */
