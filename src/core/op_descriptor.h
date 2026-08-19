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
struct generic_layer_args;

typedef enum {
	OP_MASK_CAPABLE      = 1 << 0,  /* op supports mask-aware application */
	OP_GEOMETRY_CHANGING = 1 << 1,  /* changes image dimensions (consumed by FLIS branch) */
	OP_EXPENSIVE         = 1 << 2,  /* reserved: NDE checkpoint policy */
	OP_REQ_RGB           = 1 << 3,  /* reserved: replay-time validation */
	OP_REQ_MONO          = 1 << 4,  /* reserved */
	OP_MASK_FROM_IMAGE   = 1 << 5,  /* mask op that DERIVES from the image
	                                 * (image -> MASK), as opposed to editing
	                                 * an existing mask.  Its NDE record pins
	                                 * the image it read, so the replay can
	                                 * re-derive the mask (nde_replay.h). */
	OP_ROI_CAPABLE       = 1 << 6,  /* the op can be computed on a
	                                 * sub-rectangle of the image.  Declares
	                                 * capability only, NOT locality: how much
	                                 * surrounding context the hook needs is a
	                                 * separate question (roi-nde-plan.md
	                                 * §1.1), and a pixel-local op needing none
	                                 * is still ROI-capable.  This is what
	                                 * generic_image_worker consults to decide
	                                 * whether a preview run may be computed on
	                                 * the ROI rectangle, and what the dialogs
	                                 * project into the overlay colour. */
} op_descriptor_flags;

typedef struct op_descriptor {
	const char *id;           /* stable identity, "area.op", never reused */
	int         version;      /* param-format version; 1 for all ops in this MR */

	/* per-op invariants; the worker fills args from these */
	int       (*image_hook)(struct generic_img_args *, fits *, int);
	int       (*mask_hook)(struct generic_mask_args *);  /* mask ops only; else NULL */
	int       (*layer_hook)(struct generic_layer_args *); /* FLIS-document ops; else NULL */
	gchar    *(*log_hook)(gpointer, log_hook_detail);
	const char *description;  /* N_() msgid — translated by the worker at fill time */
	float       mem_ratio;    /* default; args->mem_ratio != 0 overrides */
	guint32     flags;

	/* Nondestructive-editing serialization (flis-nde-sketch.md §11): both
	 * set for Tier-A (replayable) ops, both NULL otherwise. */
	gchar    *(*serialize)(gconstpointer user);
	gpointer  (*deserialize)(const gchar *blob, int version);

	/* Optional replay preparation (nde-phase2-3-plan.md decision 6): called
	 * by the replay driver after deserialize and before applying the record,
	 * with the parsed params table and the fits the record will run on.
	 * For ops whose inputs live outside the params struct — background
	 * extraction reinstalls its recorded sample positions here.  Return
	 * non-zero to abort the replay.  NULL for everything else. */
	int (*replay_pre)(gpointer user, GHashTable *kv, fits *target);

	/* Region previews (the ROI model, above roi_t in core/siril.h).
	 *
	 * Pixels of INPUT CONTEXT the hook needs beyond the requested region for
	 * its result INSIDE that region to match a full-image run.  The worker
	 * grows the crop by this much and writes back only the requested
	 * rectangle, so the halo never shrinks what the user sees.
	 *
	 * NULL means zero — a pixel-local op, which is fully region-capable.  It
	 * is NOT "cannot do regions"; that is the absence of OP_ROI_CAPABLE.
	 * Conflating the two would disable region previews for every stretch.
	 *
	 * @user is the params struct — the same argument serialize() takes — so
	 * one implementation serves both the interactive path (args->user) and a
	 * replay (the struct from deserialize).
	 *
	 * roi_halo_exact says the value is a true bound on the op's SPATIAL
	 * SUPPORT: no pixel further away can influence the result.  It does not
	 * claim bit-identical arithmetic, and for OpenCV-backed ops it is not —
	 * the same pixels summed in a different order (different image width,
	 * different SIMD blocking) differ by a few units in the last place, about
	 * 0.01 ADU in 65535 measured.  Iterative and non-local ops (deconvolution,
	 * NL-means) have no such bound at all: the agreed policy for them is a
	 * flat one radius, documented as an approximation, and FALSE here.
	 *
	 * The property test (src/tests/roi_halo_test.c) checks the claim rather
	 * than trusting it: negligible deviation at the declared halo AND a large
	 * one below it, so a halo that is too small fails and one that is merely
	 * generous cannot pass by accident. */
	int      (*roi_halo)(gconstpointer user);
	gboolean roi_halo_exact;
} op_descriptor;

/* Fill args fields from args->op, if set.  No-op when args->op == NULL, so the
 * legacy (un-migrated) path is byte-for-byte unchanged.  Factored out of the
 * worker so the unit test can exercise the fill semantics directly. */
/* TRUE when @op declares OP_ROI_CAPABLE.  NULL-safe: an un-migrated call site
 * with no descriptor is not ROI-capable. */
static inline gboolean op_descriptor_is_roi_capable(const op_descriptor *op) {
	return op && (op->flags & OP_ROI_CAPABLE) != 0;
}

/* Halo in pixels for @op with @user params; 0 when it declares none.  Negative
 * returns are clamped away, so a buggy implementation cannot shrink the crop
 * below the requested region. */
static inline int op_descriptor_roi_halo(const op_descriptor *op, gconstpointer user) {
	if (!op || !op->roi_halo)
		return 0;
	const int h = op->roi_halo(user);
	return h > 0 ? h : 0;
}

void op_descriptor_fill_img_args(struct generic_img_args *args);
void op_descriptor_fill_mask_args(struct generic_mask_args *args);
void op_descriptor_fill_layer_args(struct generic_layer_args *args);

/* Enumerate every descriptor (for tests and the NDE by-id registry). */
const op_descriptor *const *op_descriptor_all(size_t *count);

/* By-id lookup over op_descriptor_all(), NULL for unknown ids.  Needed only
 * where an op arrives as a string (FLIS_HIST deserialization) — live capture
 * reads args->op directly.  Thread-safe; the table is built on first use. */
const op_descriptor *op_descriptor_by_id(const char *id);

#ifdef __cplusplus
}
#endif

#endif /* _OP_DESCRIPTOR_H_ */
