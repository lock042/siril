#define MAIN
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "core/siril.h"

/* the global variables of the whole project */
cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
fits gfit;	// currently loaded image

