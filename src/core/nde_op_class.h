#ifndef _NDE_OP_CLASS_H_
#define _NDE_OP_CLASS_H_

#include <glib.h>
#include "core/op_descriptor.h"

#ifdef __cplusplus
extern "C" {
#endif

/* What KIND of record an op_id produces — asked by the replay, the graph, the
 * history and the GUI, and until this module answered by each of them running
 * its own list of string comparisons.
 *
 * Three families have no op_descriptor at all (composite, compositing-state,
 * structural): they are captured by nde_capture_structural() and replayed by
 * nde_composite.c / nde_compositing.c, which is why classification cannot live
 * in op_descriptor_flags and needs a registry of its own.
 *
 * If you are adding an ordinary image operation you do not need this file.
 * Your op gets NDE_OPC_PIXEL from its descriptor with nothing to declare.
 */
typedef enum {
	NDE_OPC_PIXEL = 0,   /* ordinary image op, applied by the image worker */
	NDE_OPC_MASK,        /* mask op, applied by the mask worker */
	NDE_OPC_COMPOSITE,   /* layer.merge_down / document.flatten */
	NDE_OPC_COMPOSITING, /* layer.set_* — an input to the compositor, not a
	                      * pixel operation, so never a chain member */
	NDE_OPC_STRUCTURAL,  /* layer add/remove/reorder/duplicate, canvas.*,
	                      * image.origin */
	NDE_OPC_JOINT,       /* the multi-layer analyses (nde_joint.h) */
	NDE_OPC_DOCUMENT,    /* document-wide pixel op, applied by the layer
	                      * worker (icc.convert) */
	NDE_OPC_ANALYSIS,    /* measures or exports; never modifies the item, and
	                      * so must never reach the NDE capture */
	NDE_OPC_UNKNOWN,     /* an id this build does not know: from a newer Siril,
	                      * or a corrupt log.  Carries NO traits on purpose */
} nde_op_family;

/* Properties that cut across families.  Each one exists because some site
 * needs exactly it; there is no trait here without a reader. */
typedef enum {
	NDE_OPT_DESTRUCTIVE       = 1 << 0, /* consumes items and replaces pixels */
	NDE_OPT_GEOMETRIC         = 1 << 1, /* moves the layer on the canvas */
	NDE_OPT_COMPOSITING_RESET = 1 << 2, /* after it, compositing sits at defaults */
	NDE_OPT_DELETES_ITEM      = 1 << 3, /* layer.remove: a dead end, not a consumption */
	NDE_OPT_INSERT_DISTURBS   = 1 << 4, /* an armed insertion cannot survive it */
	NDE_OPT_CHAIN_IGNORE      = 1 << 5, /* neither a chain member nor a blocker */
} nde_op_trait;

typedef struct {
	const char          *id;
	nde_op_family        family;
	guint32              traits;
	const op_descriptor *desc;   /* NULL for the descriptor-less families */
	const char          *label;  /* N_()-marked; translate at the point of use */
} nde_op_class;

/**
 * Classify @op_id.  Never NULL: an id this build does not know resolves to a
 * shared NDE_OPC_UNKNOWN class with no traits, which is what makes unknown ops
 * fail CLOSED — lacking NDE_OPT_CHAIN_IGNORE, a document-scope one still
 * blocks a chain exactly as it did when that was a string comparison.
 *
 * Thread-safe; the table is built once and is read-only thereafter.
 */
const nde_op_class *nde_op_class_for(const char *op_id);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_OP_CLASS_H_ */
