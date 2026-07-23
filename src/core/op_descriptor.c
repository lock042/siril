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
extern const op_descriptor op_desc_addmax;
extern const op_descriptor op_desc_asinh;
extern const op_descriptor op_desc_atrous;
extern const op_descriptor op_desc_autoghs;
extern const op_descriptor op_desc_autoghs_unlinked;
extern const op_descriptor op_desc_banding;
extern const op_descriptor op_desc_bg;
extern const op_descriptor op_desc_bgnoise;
extern const op_descriptor op_desc_binning;
extern const op_descriptor op_desc_catsearch;
extern const op_descriptor op_desc_ccm;
extern const op_descriptor op_desc_cdg;
extern const op_descriptor op_desc_cfa_extract_green;
extern const op_descriptor op_desc_cfa_extract_ha;
extern const op_descriptor op_desc_cfa_extract_haoiii;
extern const op_descriptor op_desc_cfa_split;
extern const op_descriptor op_desc_clahe;
extern const op_descriptor op_desc_cosme;
extern const op_descriptor op_desc_cosmetic;
extern const op_descriptor op_desc_crop;
extern const op_descriptor op_desc_curves;
extern const op_descriptor op_desc_ddp;
extern const op_descriptor op_desc_deconvolve;
extern const op_descriptor op_desc_denoise;
extern const op_descriptor op_desc_entropy;
extern const op_descriptor op_desc_epf;
extern const op_descriptor op_desc_fdiv;
extern const op_descriptor op_desc_ffill;
extern const op_descriptor op_desc_fft;
extern const op_descriptor op_desc_fill;
extern const op_descriptor op_desc_findhot;
extern const op_descriptor op_desc_fix_xtrans;
extern const op_descriptor op_desc_fmul;
extern const op_descriptor op_desc_gauss;
extern const op_descriptor op_desc_ghs;
extern const op_descriptor op_desc_grey_flat;
extern const op_descriptor op_desc_icc_assign;
extern const op_descriptor op_desc_icc_convert;
extern const op_descriptor op_desc_icc_remove;
extern const op_descriptor op_desc_imoper;
extern const op_descriptor op_desc_limit;
extern const op_descriptor op_desc_linear_match;
extern const op_descriptor op_desc_logstretch;
extern const op_descriptor op_desc_mask_autostretch;
extern const op_descriptor op_desc_mask_bitpix;
extern const op_descriptor op_desc_mask_blur;
extern const op_descriptor op_desc_mask_clear;
extern const op_descriptor op_desc_mask_feather;
extern const op_descriptor op_desc_mask_from_channel;
extern const op_descriptor op_desc_mask_from_color;
extern const op_descriptor op_desc_mask_from_gradient;
extern const op_descriptor op_desc_mask_from_luminance;
extern const op_descriptor op_desc_mask_from_stars;
extern const op_descriptor op_desc_mask_invert;
extern const op_descriptor op_desc_mask_mtf_autostretch;
extern const op_descriptor op_desc_mask_multiply;
extern const op_descriptor op_desc_mask_threshold;
extern const op_descriptor op_desc_median;
extern const op_descriptor op_desc_mirrorx;
extern const op_descriptor op_desc_mirrory;
extern const op_descriptor op_desc_mtf;
extern const op_descriptor op_desc_mtf_inverse;
extern const op_descriptor op_desc_neg;
extern const op_descriptor op_desc_nozero;
extern const op_descriptor op_desc_offset;
extern const op_descriptor op_desc_photometric_cc;
extern const op_descriptor op_desc_psf_estimate;
extern const op_descriptor op_desc_remove_gradient;
extern const op_descriptor op_desc_resample;
extern const op_descriptor op_desc_rgradient;
extern const op_descriptor op_desc_rotation;
extern const op_descriptor op_desc_saturation;
extern const op_descriptor op_desc_scnr;
extern const op_descriptor op_desc_stat;
extern const op_descriptor op_desc_synthstar;
extern const op_descriptor op_desc_thresh;
extern const op_descriptor op_desc_unclip;
extern const op_descriptor op_desc_unpurple;
extern const op_descriptor op_desc_unsharp;
extern const op_descriptor op_desc_wrecons;

static const op_descriptor *const descriptors[] = {
	&op_desc_addmax,              /* arith.addmax */
	&op_desc_fdiv,                /* arith.fdiv */
	&op_desc_ffill,               /* arith.ffill */
	&op_desc_fill,                /* arith.fill */
	&op_desc_fmul,                /* arith.fmul */
	&op_desc_imoper,              /* arith.imoper */
	&op_desc_limit,               /* arith.limit */
	&op_desc_neg,                 /* arith.neg */
	&op_desc_nozero,              /* arith.nozero */
	&op_desc_offset,              /* arith.offset */
	&op_desc_thresh,              /* arith.thresh */
	&op_desc_remove_gradient,     /* bkg.remove_gradient */
	&op_desc_catsearch,           /* catalog.search */
	&op_desc_cfa_extract_green,   /* cfa.extract_green */
	&op_desc_cfa_extract_ha,      /* cfa.extract_ha */
	&op_desc_cfa_extract_haoiii,  /* cfa.extract_haoiii */
	&op_desc_findhot,             /* cfa.findhot */
	&op_desc_fix_xtrans,          /* cfa.fix_xtrans */
	&op_desc_cfa_split,           /* cfa.split */
	&op_desc_ccm,                 /* color.ccm */
	&op_desc_grey_flat,           /* color.grey_flat */
	&op_desc_linear_match,        /* color.linear_match */
	&op_desc_photometric_cc,      /* color.photometric_cc */
	&op_desc_saturation,          /* color.saturation */
	&op_desc_banding,             /* filters.banding */
	&op_desc_clahe,               /* filters.clahe */
	&op_desc_cosme,               /* filters.cosme */
	&op_desc_cosmetic,            /* filters.cosmetic */
	&op_desc_ddp,                 /* filters.ddp */
	&op_desc_deconvolve,          /* filters.deconvolve */
	&op_desc_denoise,             /* filters.denoise */
	&op_desc_epf,                 /* filters.epf */
	&op_desc_fft,                 /* filters.fft */
	&op_desc_gauss,               /* filters.gauss */
	&op_desc_median,              /* filters.median */
	&op_desc_rgradient,           /* filters.rgradient */
	&op_desc_scnr,                /* filters.scnr */
	&op_desc_unpurple,            /* filters.unpurple */
	&op_desc_unsharp,             /* filters.unsharp */
	&op_desc_binning,             /* geometry.binning */
	&op_desc_crop,                /* geometry.crop */
	&op_desc_mirrorx,             /* geometry.mirrorx */
	&op_desc_mirrory,             /* geometry.mirrory */
	&op_desc_resample,            /* geometry.resample */
	&op_desc_rotation,            /* geometry.rotation */
	&op_desc_icc_assign,          /* icc.assign */
	&op_desc_icc_convert,         /* icc.convert */
	&op_desc_icc_remove,          /* icc.remove */
	&op_desc_mask_autostretch,    /* mask.autostretch */
	&op_desc_mask_bitpix,         /* mask.bitpix */
	&op_desc_mask_blur,           /* mask.blur */
	&op_desc_mask_clear,          /* mask.clear */
	&op_desc_mask_feather,        /* mask.feather */
	&op_desc_mask_from_channel,   /* mask.from_channel */
	&op_desc_mask_from_color,     /* mask.from_color */
	&op_desc_mask_from_gradient,  /* mask.from_gradient */
	&op_desc_mask_from_luminance, /* mask.from_luminance */
	&op_desc_mask_from_stars,     /* mask.from_stars */
	&op_desc_mask_invert,         /* mask.invert */
	&op_desc_mask_mtf_autostretch, /* mask.mtf_autostretch */
	&op_desc_mask_multiply,       /* mask.multiply */
	&op_desc_mask_threshold,      /* mask.threshold */
	&op_desc_psf_estimate,        /* psf.estimate */
	&op_desc_synthstar,           /* star.synthstar */
	&op_desc_unclip,              /* star.unclip */
	&op_desc_bg,                  /* stats.bg */
	&op_desc_bgnoise,             /* stats.bgnoise */
	&op_desc_cdg,                 /* stats.cdg */
	&op_desc_entropy,             /* stats.entropy */
	&op_desc_stat,                /* stats.stat */
	&op_desc_asinh,               /* stretch.asinh */
	&op_desc_autoghs,             /* stretch.autoghs */
	&op_desc_autoghs_unlinked,    /* stretch.autoghs_unlinked */
	&op_desc_curves,              /* stretch.curves */
	&op_desc_ghs,                 /* stretch.ghs */
	&op_desc_logstretch,          /* stretch.log */
	&op_desc_mtf,                 /* stretch.mtf */
	&op_desc_mtf_inverse,         /* stretch.mtf_inverse */
	&op_desc_atrous,              /* wavelets.atrous */
	&op_desc_wrecons,             /* wavelets.wrecons */
	NULL
};

const op_descriptor *const *op_descriptor_all(size_t *count) {
	/* -1 for the trailing NULL sentinel */
	if (count)
		*count = (sizeof(descriptors) / sizeof(descriptors[0])) - 1;
	return descriptors;
}
