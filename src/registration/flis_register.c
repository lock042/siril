/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * flis_register — layer-to-layer registration for FLIS images.
 *
 * The Siril registration pipeline (register_star_alignment,
 * register_shift_dft, register_kombat, register_apply_reg, …) operates
 * on a `sequence`.  FLIS layer registration reuses that machinery by
 * building an internal sequence whose entries are the layer fits
 * pointers, computing the transforms, applying them (FRAMING_MAX so
 * each layer is resampled into its own bounding box), then writing
 * the per-layer offsets back to each layer's position_x/y on the
 * canvas.
 *
 * Method choice is exposed to the caller (dialog + command).  The
 * default — single-pass global star alignment — works on a freshly-
 * loaded FLIS without any user input.  DFT shift and KOMBAT require
 * the user to first make a selection on the image; the validation is
 * centralised here so both code paths (panel dialog, CLI command)
 * share the same refusal logic.
 *
 * Adapted from the flis branch's gui/flis_gui.c flis_register_worker,
 * stripped of all GUI state so the same primitive serves both the
 * layers-panel "Register Layers…" dialog and the `flis_register_layers`
 * command (§C.2).
 */

#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing_thread.h"
#include "io/sequence.h"
#include "io/image_format_flis.h"
#include "registration/registration.h"
#include "opencv/opencv.h"

#include "flis_register.h"

registration_function flis_register_resolve_method(flis_reg_method_id id,
                                                    selection_type *out_sel,
                                                    transformation_type *out_tx) {
	switch (id) {
		case FLIS_REG_GLOBAL:
			/* Two-pass global alignment, same as FLIS_REG_2PASS.  The
			 * single-pass register_star_alignment is NOT usable here: it
			 * resamples the layer pixels inline (apply_reg_image_hook runs
			 * inside the method), so the register_apply_reg() pass that
			 * flis_register_layers performs afterwards would warp the
			 * already-warped pixels a second time and every layer would
			 * end up misplaced.  register_multi_step_global only computes
			 * the transforms; the single apply happens in
			 * flis_register_layers. */
			if (out_sel) *out_sel = REQUIRES_NO_SELECTION;
			if (out_tx)  *out_tx  = HOMOGRAPHY_TRANSFORMATION;
			return register_multi_step_global;
		case FLIS_REG_2PASS:
			if (out_sel) *out_sel = REQUIRES_NO_SELECTION;
			if (out_tx)  *out_tx  = HOMOGRAPHY_TRANSFORMATION;
			return register_multi_step_global;
		case FLIS_REG_DFT:
			if (out_sel) *out_sel = REQUIRES_SQUARED_SELECTION;
			if (out_tx)  *out_tx  = SHIFT_TRANSFORMATION;
			return register_shift_dft;
		case FLIS_REG_KOMBAT:
			if (out_sel) *out_sel = REQUIRES_ANY_SELECTION;
			if (out_tx)  *out_tx  = SHIFT_TRANSFORMATION;
			return register_kombat;
		default:
			return NULL;
	}
}

/* Centralised selection-requirement check shared between dialog and
 * command paths.  Logs a user-facing message and returns FALSE on
 * failure; TRUE when the selection (if any) satisfies @sel_req. */
static gboolean check_selection_requirement(selection_type sel_req) {
	const gboolean has_sel = (com.selection.w > 0 && com.selection.h > 0);
	switch (sel_req) {
		case REQUIRES_NO_SELECTION:
			return TRUE;
		case REQUIRES_ANY_SELECTION:
			if (!has_sel) {
				siril_log_error(_("flis_register_layers: this method requires "
				                  "an image selection — drag a rectangle on "
				                  "the image and try again\n"));
				return FALSE;
			}
			return TRUE;
		case REQUIRES_SQUARED_SELECTION:
			if (!has_sel) {
				siril_log_error(_("flis_register_layers: this method requires "
				                  "a square image selection — drag a square "
				                  "rectangle on the image and try again\n"));
				return FALSE;
			}
			if (com.selection.w != com.selection.h) {
				siril_log_error(_("flis_register_layers: this method requires "
				                  "a SQUARE image selection — current "
				                  "selection is %dx%d\n"),
				                com.selection.w, com.selection.h);
				return FALSE;
			}
			return TRUE;
	}
	return TRUE;
}

int flis_register_layers(flis_layer_t *ref_lay,
                         GSList       *target_layers,
                         registration_function method,
                         selection_type        sel_req,
                         transformation_type   tx_type,
                         opencv_interpolation interpolation,
                         gboolean       clamp) {
	if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) {
		siril_log_error(_("flis_register_layers: no FLIS image loaded\n"));
		return 1;
	}
	GSList *targets = target_layers ? target_layers : com.uniq->layers;
	gint n_layers = (gint)g_slist_length(targets);
	if (n_layers < 2) {
		siril_log_error(_("flis_register_layers: need at least two layers\n"));
		return 1;
	}

	/* Method defaults: single-pass global star alignment — the safest
	 * choice for a freshly-loaded FLIS (no selection required, works on
	 * any deep-sky content). */
	if (!method) {
		selection_type s;
		transformation_type t;
		method  = flis_register_resolve_method(FLIS_REG_GLOBAL, &s, &t);
		sel_req = s;
		tx_type = t;
	}

	if (!check_selection_requirement(sel_req))
		return 1;

	/* The registration method and register_apply_reg() below drive
	 * generic_sequence_worker synchronously (already_in_a_thread), whose
	 * frame loop polls processing_should_continue().  When we are NOT
	 * already inside the processing worker (the command / script path),
	 * we must claim the job slot via reserve_thread(): it clears any
	 * stale cancel_flag left by a previous job's stop_processing_thread()
	 * — without this the very first frame aborts and registration fails
	 * with a bare "Sequence processing failed".  The GUI panel path runs
	 * inside a queued generic_layer_worker job which already cleared the
	 * flag and owns the slot, so reserving again would (rightly) fail. */
	gboolean reserved = FALSE;
	if (!processing_in_worker_thread()) {
		if (!reserve_thread()) {
			siril_log_error(_("flis_register_layers: another processing "
			                  "task is running, try again later\n"));
			return 1;
		}
		reserved = TRUE;
	}

	flis_layer_t *canvas_lay = (flis_layer_t *)com.uniq->layers->data;
	flis_layer_t *ref = ref_lay ? ref_lay : (flis_layer_t *)targets->data;
	if (!ref || !ref->fit) {
		siril_log_error(_("flis_register_layers: reference layer has no pixel data\n"));
		if (reserved) unreserve_thread();
		return 1;
	}

	const gint ref_orig_x = ref->position_x;
	const gint ref_orig_y = ref->position_y;

	/* Build an internal sequence: slot 0 = reference, remaining slots in
	 * target_layers order (skipping the ref).  Track the canvas (base)
	 * layer's slot so we can normalise offsets relative to it after the
	 * registration completes. */
	sequence *seq = create_internal_sequence(n_layers);
	if (!seq) {
		siril_log_error(_("flis_register_layers: could not build internal sequence\n"));
		if (reserved) unreserve_thread();
		return 1;
	}
	seq->bitpix = ref->fit->bitpix;
	seq->rx     = ref->fit->rx;
	seq->ry     = ref->fit->ry;

	const int ref_seq_idx = 0;
	int canvas_seq_idx = (ref == canvas_lay) ? 0 : -1;

	internal_sequence_set(seq, 0, ref->fit);
	seq->imgparam[0].rx = ref->fit->rx;
	seq->imgparam[0].ry = ref->fit->ry;

	int i = 1;
	gboolean is_variable = FALSE;
	for (GSList *l = targets; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay || !lay->fit) continue;
		if (lay == ref) continue;
		if (lay == canvas_lay) canvas_seq_idx = i;
		if (lay->fit->rx != ref->fit->rx || lay->fit->ry != ref->fit->ry)
			is_variable = TRUE;
		seq->imgparam[i].rx = lay->fit->rx;
		seq->imgparam[i].ry = lay->fit->ry;
		internal_sequence_set(seq, i++, lay->fit);
	}
	seq->reference_image = ref_seq_idx;
	seq->is_variable     = is_variable;

	/* Registration arguments. */
	struct registration_args regargs = { 0 };
	regargs.seq                  = seq;
	regargs.layer                = 0;
	regargs.reference_image      = ref_seq_idx;
	regargs.run_in_thread        = FALSE;
	regargs.interpolation        = interpolation;
	regargs.output_scale         = 1.f;
	regargs.clamp                = clamp;
	regargs.framing              = FRAMING_MAX;
	regargs.type                 = tx_type;
	regargs.max_stars_candidates = MAX_STARS_FITTED;
	regargs.two_pass             = (method == register_multi_step_global);
	regargs.percent_moved        = 0.50f;
	/* register_star_alignment / register_multi_step_global do an
	 * unconditional strdup(regargs.prefix) for the seq-output naming
	 * scaffolding.  Internal sequences never produce an on-disk file so
	 * the prefix is academic — but the strdup() must see a non-NULL
	 * pointer.  Without this the registration thread segfaults on the
	 * first call, easy to miss because the user hasn't named anything. */
	regargs.prefix               = "flis_";
	/* Selection is forwarded for methods that need it; the requirement
	 * check above already gated this. */
	if (sel_req != REQUIRES_NO_SELECTION) {
		regargs.selection = com.selection;
	} else {
		regargs.selection.x = regargs.selection.y = 0;
		regargs.selection.w = regargs.selection.h = 0;
	}

	int ret = method(&regargs);
	if (ret) {
		siril_log_error(_("flis_register_layers: alignment step failed\n"));
		free(regargs.imgparam); regargs.imgparam = NULL;
		free(regargs.regparam); regargs.regparam = NULL;
		free_sequence(seq, TRUE);
		if (reserved) unreserve_thread();
		return 1;
	}
	free(regargs.imgparam); regargs.imgparam = NULL;

	/* Apply transforms — resamples layer pixel data + records bounding-box
	 * offsets in regparam[i].H.h02 / h12. */
	regargs.framing = FRAMING_MAX;
	ret = register_apply_reg(&regargs);
	if (!ret && regargs.regparam) {
		double cx, cy;
		if (canvas_seq_idx >= 0) {
			cx = regargs.regparam[canvas_seq_idx].H.h02;
			cy = regargs.regparam[canvas_seq_idx].H.h12;
			canvas_lay->position_x = 0;
			canvas_lay->position_y = 0;
		} else {
			cx = regargs.regparam[0].H.h02 - ref_orig_x;
			cy = regargs.regparam[0].H.h12 - ref_orig_y;
		}

		int k = 1;
		for (GSList *l = targets; l; l = l->next) {
			flis_layer_t *lay = (flis_layer_t *)l->data;
			if (!lay || !lay->fit) continue;
			if (lay == canvas_lay) {
				if (lay != ref) k++;
				continue;
			}
			int seq_idx = (lay == ref) ? 0 : k++;
			double dx = regargs.regparam[seq_idx].H.h02 - cx;
			double dy = regargs.regparam[seq_idx].H.h12 - cy;
			lay->position_x = (gint)round(dx);
			lay->position_y = (gint)round(dy);
		}

		/* FRAMING_MAX resamples the canvas (base) layer into its own
		 * bounding box, which can be larger than the document canvas
		 * (rotation, or a reference whose frame extends past it).  On
		 * the original flis branch the canvas implicitly followed the
		 * base layer's dimensions, so registration never clipped it;
		 * with the §7 canvas-independent model we reproduce that by
		 * growing the canvas to the registered base layer's new size
		 * (still pinned at the origin, no layer shifts). */
		if (canvas_seq_idx >= 0 && canvas_lay->fit &&
		    (canvas_lay->fit->rx != com.uniq->canvas_w ||
		     canvas_lay->fit->ry != com.uniq->canvas_h)) {
			siril_log_message(_("flis_register_layers: canvas resized to "
			                    "%ux%u to match the registered base layer\n"),
			                  canvas_lay->fit->rx, canvas_lay->fit->ry);
			flis_canvas_resize(canvas_lay->fit->rx, canvas_lay->fit->ry, 0, 0);
		}
	}

	free(regargs.imgparam); regargs.imgparam = NULL;
	free(regargs.regparam); regargs.regparam = NULL;
	free_sequence(seq, TRUE);
	if (reserved) unreserve_thread();
	return ret;
}
