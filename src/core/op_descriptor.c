/*
 * This file is part of Siril, an astronomy image processor.
 *
 * Per-operation descriptors: the single source of truth for the generic
 * image-processing framework.  See op_descriptor.h for the design rationale.
 */

#include "core/siril.h"
#include "core/processing.h"
#include "core/op_descriptor.h"

/* Fill args fields from args->op, if set.  Called at the top of
 * generic_image_worker before any of the affected fields are read.  When
 * args->op is NULL (an un-migrated site) this is a no-op and the worker behaves
 * exactly as before. */
void op_descriptor_fill_img_args(struct generic_img_args *args) {
	if (!args->op)
		return;
	const op_descriptor *op = args->op;
	g_assert(op->image_hook != NULL);
	/* Setting both a descriptor and a per-site hook/log_hook is a programming
	 * error: migrated sites hand these to the descriptor exclusively. */
	g_assert(args->image_hook == NULL);
	g_assert(args->log_hook == NULL);
	args->image_hook = op->image_hook;
	args->log_hook = op->log_hook;             /* may be NULL */
	if (!args->description)
		args->description = _(op->description); /* site may pre-set to override */
	if (args->mem_ratio == 0.0f)
		args->mem_ratio = op->mem_ratio;       /* 0 in args means "use default" */
	if (args->mask_aware)
		g_warn_if_fail(op->flags & OP_MASK_CAPABLE);
}

void op_descriptor_fill_mask_args(struct generic_mask_args *args) {
	if (!args->op)
		return;
	const op_descriptor *op = args->op;
	g_assert(op->mask_hook != NULL);
	g_assert(args->mask_hook == NULL);
	g_assert(args->log_hook == NULL);
	args->mask_hook = op->mask_hook;
	args->log_hook = op->log_hook;             /* may be NULL */
	if (!args->description)
		args->description = _(op->description);
	if (args->mem_ratio == 0.0f)
		args->mem_ratio = op->mem_ratio;
}

/* ---------------------------------------------------------------------------
 * Registry — every descriptor in the codebase, alphabetised by id.  Descriptors
 * are defined next to their hooks and declared extern in the owning module's
 * header; list them here so op_descriptor_all() can enumerate them (used by the
 * unit test today and by the NDE by-id registry later).
 * ------------------------------------------------------------------------- */

/* Descriptors are defined next to their hooks in the owning modules; declared
 * here so the registry can reference them without pulling in every module
 * header.  Keep both this block and the array alphabetised by id. */
extern const op_descriptor op_desc_binning;
extern const op_descriptor op_desc_crop;
extern const op_descriptor op_desc_mirrorx;
extern const op_descriptor op_desc_mirrory;
extern const op_descriptor op_desc_resample;
extern const op_descriptor op_desc_rotation;

static const op_descriptor *const descriptors[] = {
	&op_desc_binning,       /* geometry.binning */
	&op_desc_crop,          /* geometry.crop */
	&op_desc_mirrorx,       /* geometry.mirrorx */
	&op_desc_mirrory,       /* geometry.mirrory */
	&op_desc_resample,      /* geometry.resample */
	&op_desc_rotation,      /* geometry.rotation */
	NULL
};

const op_descriptor *const *op_descriptor_all(size_t *count) {
	/* -1 for the trailing NULL sentinel */
	if (count)
		*count = (sizeof(descriptors) / sizeof(descriptors[0])) - 1;
	return descriptors;
}
