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


#ifndef _WIN32
long get_pathmax(void);
#endif

typedef struct starnet_data {
	struct ser_struct *new_ser_starless;
	fitseq *new_fitseq_starless;
	struct ser_struct *new_ser_starmask;
	fitseq *new_fitseq_starmask;
	gboolean force_ser;
	sequence *seq;
	fits *starnet_fit;
	fits *starmask_fit;
	struct timeval t_start;
	gchar stride[6];
	gboolean linear;
	gboolean customstride;
	gboolean upscale;
	gboolean starmask;
	gboolean follow_on;
	const gchar *seqname;
	const gchar *seqEntry;
	int imgnumber;
} starnet_data;

typedef struct remixargs {
	fits *fit1;
	fits *fit2;
} remixargs;

gboolean starnet_executablecheck();
gpointer do_starnet(gpointer p);
void apply_starnet_to_sequence(struct starnet_data *seqdata);

#endif /* SRC_FILTERS_STARNET_H_ */
