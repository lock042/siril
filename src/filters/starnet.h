#ifndef SRC_FILTERS_STARNET_H_
#define SRC_FILTERS_STARNET_H_

#include "core/siril.h"

#ifdef _WIN32
#define STARNET_BIN "starnet++.exe"
#define STARNET_RGB "rgb_starnet++.exe"
#define STARNET_MONO "mono_starnet++.exe"
#else
#define STARNET_BIN "starnet++"
#define STARNET_RGB "rgb_starnet++"
#define STARNET_MONO "mono_starnet++"
#endif

gboolean starnet_executablecheck();
gpointer do_starnet(gpointer p);

typedef struct starnet_data {
	struct timeval t_start;
	gchar stride[6];
	gboolean linear;
	gboolean customstride;
	gboolean upscale;
	gboolean starmask;
} starnet_data;

#endif /* SRC_FILTERS_STARNET_H_ */
