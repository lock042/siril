#ifndef _OP_DESCRIPTOR_H_
#define _OP_DESCRIPTOR_H_

/**
 * \file op_descriptor.h
 * \brief Per-operation descriptor structs — the single source of truth for the
 * generic image-processing framework.
 *
 * A single \c op_descriptor captures the per-op invariants (image or mask hook,
 * log hook, progress description, memory ratio, capability flags) that were
 * previously duplicated across every GUI + command construction site of a
 * \c generic_img_args / \c generic_mask_args.  Descriptors are \c const data
 * defined next to the hooks they reference; call sites point \c args->op at the
 * relevant descriptor and the worker fills the args fields from it (see
 * op_descriptor_fill_img_args() / op_descriptor_fill_mask_args()).
 *
 * The \c serialize / \c deserialize members are reserved for the
 * nondestructive-editing branch (flis-nde-sketch.md §11) and are always NULL in
 * this MR.
 */

#include "core/siril.h"   /* fits, log_hook_detail */

#ifdef __cplusplus
extern "C" {
#endif

struct generic_img_args;
struct generic_mask_args;

typedef enum {
	OP_MASK_CAPABLE      = 1 << 0,  /* op supports mask-aware application */
	OP_GEOMETRY_CHANGING = 1 << 1,  /* changes image dimensions (consumed by FLIS branch) */
	OP_EXPENSIVE         = 1 << 2,  /* reserved: NDE checkpoint policy */
	OP_REQ_RGB           = 1 << 3,  /* reserved: replay-time validation */
	OP_REQ_MONO          = 1 << 4,  /* reserved */
} op_descriptor_flags;

typedef struct op_descriptor {
	const char *id;           /* stable identity, "area.op", never reused */
	int         version;      /* param-format version; 1 for all ops in this MR */

	/* per-op invariants; the worker fills args from these */
	int       (*image_hook)(struct generic_img_args *, fits *, int);
	int       (*mask_hook)(struct generic_mask_args *);  /* mask ops only; else NULL */
	gchar    *(*log_hook)(gpointer, log_hook_detail);
	const char *description;  /* N_() msgid — translated by the worker at fill time */
	float       mem_ratio;    /* default; args->mem_ratio != 0 overrides */
	guint32     flags;

	/* Reserved for nondestructive-editing work (flis-nde-sketch.md §11): always
	 * NULL in this MR.  Present so the NDE branch adds implementations without
	 * changing this struct or any call site. */
	gchar    *(*serialize)(gconstpointer user);
	gpointer  (*deserialize)(const gchar *blob, int version);
} op_descriptor;

/* Fill args fields from args->op, if set.  No-op when args->op == NULL, so the
 * legacy (un-migrated) path is byte-for-byte unchanged.  Factored out of the
 * worker so the unit test can exercise the fill semantics directly. */
void op_descriptor_fill_img_args(struct generic_img_args *args);
void op_descriptor_fill_mask_args(struct generic_mask_args *args);

/* Enumerate every descriptor (for tests and, later, the NDE by-id registry). */
const op_descriptor *const *op_descriptor_all(size_t *count);

#ifdef __cplusplus
}
#endif

#endif /* _OP_DESCRIPTOR_H_ */
