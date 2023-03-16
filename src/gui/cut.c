#include "core/siril.h"
#include "gui/image_interactions.h"
#include "core/processing.h"

gpointer cut_profile(gpointer p) {
	cut_args *args = (cut_args *) p;
	int retval = 0;
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}
