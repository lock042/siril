#ifndef SRC_FILTERS_SCNR_H_
#define SRC_FILTERS_SCNR_H_

#include "core/siril.h"

typedef enum {
	SCNR_AVERAGE_NEUTRAL,
	SCNR_MAXIMUM_NEUTRAL,
	SCNR_MAXIMUM_MASK,
	SCNR_ADDITIVE_MASK
} scnr_type;

/* scnr data from GUI */
struct scnr_data {
	fits *fit;
	scnr_type type;
	double amount;
	gboolean preserve;
};

gpointer scnr(gpointer p);

const char *scnr_type_to_string(scnr_type t);

#endif /* SRC_FILTERS_SCNR_H_ */
