#ifndef SRC_GUI_GEOMETRY_H_
#define SRC_GUI_GEOMETRY_H_

#include "core/siril.h"

void siril_rotate90();
void siril_rotate270();

void mirrorx_gui(fits *fit);
void mirrory_gui(fits *fit);

void siril_crop();

#endif

