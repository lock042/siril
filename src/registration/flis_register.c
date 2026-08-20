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
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/undo.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_joint.h"
#include "core/nde_script_scope.h"
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

gboolean flis_register_selection_ok(selection_type sel_req,
                                    const rectangle *sel,
                                    int img_rx, int img_ry, gchar **err) {
	if (err) *err = NULL;
	if (sel_req == REQUIRES_NO_SELECTION)
		return TRUE;

	const rectangle empty = { 0 };
	if (!sel) sel = &empty;
	if (sel->w <= 0 || sel->h <= 0) {
		if (err)
			*err = g_strdup(sel_req == REQUIRES_SQUARED_SELECTION
					? _("this method needs a square image selection — drag a "
					    "square on the image and try again")
					: _("this method needs an image selection — drag a "
					    "rectangle on the image and try again"));
		return FALSE;
	}
	if (sel_req == REQUIRES_SQUARED_SELECTION && sel->w != sel->h) {
		if (err)
			*err = g_strdup_printf(_("this method needs a SQUARE image "
			                         "selection — this one is %dx%d"),
			                       sel->w, sel->h);
		return FALSE;
	}
	/* The template is cut from the reference frame, so one that runs off the
	 * edge is as fatal as one that is empty. */
	if (img_rx > 0 && img_ry > 0 &&
	    (sel->x < 0 || sel->y < 0 ||
	     sel->x + sel->w > img_rx || sel->y + sel->h > img_ry)) {
		if (err)
			*err = g_strdup_printf(_("the selection (%d,%d %dx%d) does not fit "
			                         "inside the %dx%d reference image"),
			                       sel->x, sel->y, sel->w, sel->h,
			                       img_rx, img_ry);
		return FALSE;
	}
	return TRUE;
}

/* Centralised selection-requirement check shared between dialog and
 * command paths.  Logs a user-facing message and returns FALSE on
 * failure; TRUE when com.selection satisfies @sel_req.  Bounds are not
 * checked here: the GUI already confines com.selection to the loaded image,
 * which is the canvas the layers are cut from. */
static gboolean check_selection_requirement(selection_type sel_req) {
	gchar *err = NULL;
	if (flis_register_selection_ok(sel_req, &com.selection, 0, 0, &err))
		return TRUE;
	siril_log_error(_("flis_register_layers: %s\n"), err);
	g_free(err);
	return FALSE;
}

/* The transform register_apply_reg ACTUALLY hands to cvTransformImage under
 * FRAMING_MAX, reconstructed for the NDE record.
 *
 * apply_reg_image_hook composes the per-frame transform @Himg with the
 * framing reference @Href and then post-multiplies the bounding-box origin
 * back out, so the warp lands at (0,0) of its own box.  It keeps NEITHER
 * result: regparam[i].H is overwritten with the residual SHIFT (which is all
 * flis_register_layers needs for the canvas offsets) and the composed matrix
 * dies with the frame.  A replay needs the composed matrix and the box size,
 * so we recompute both from data that DOES survive the call —
 * seq->regparam (the input transforms, untouched by apply) and
 * regargs->framingd.Htransf (the framing reference the apply computed).
 *
 * This mirrors compute_Hmax() in applyreg.c, which is static; @scale is
 * fixed at 1 here because flis_register_layers never sets output_scale.
 * If that function's framing arithmetic ever changes, the L1 oracle test in
 * test_flis_register_nde.c is what will catch the drift. */
static void flis_reg_compute_framed_transform(const Homography *Himg,
                                              const Homography *Href,
                                              int src_rx, int src_ry,
                                              Homography *H_out,
                                              int *dst_rx, int *dst_ry) {
	double xs[4] = { 0., (double)src_rx - 1., (double)src_rx - 1., 0. };
	double ys[4] = { 0., 0., (double)src_ry - 1., (double)src_ry - 1. };
	double xmin = DBL_MAX, xmax = -DBL_MAX, ymin = DBL_MAX, ymax = -DBL_MAX;
	for (int j = 0; j < 4; j++) {
		cvTransfPoint(&xs[j], &ys[j], *Himg, *Href, 1.);
		if (xmin > xs[j]) xmin = xs[j];
		if (xmax < xs[j]) xmax = xs[j];
		if (ymin > ys[j]) ymin = ys[j];
		if (ymax < ys[j]) ymax = ys[j];
	}
	*dst_rx = (int)xmax - (int)xmin + 1;
	*dst_ry = (int)ymax - (int)ymin + 1;
	Homography H = { 0 }, Hcorr = { 0 };
	cvTransfH((Homography *)Himg, (Homography *)Href, &H);
	cvGetEye(&Hcorr);
	Hcorr.h02 = -round(xmin);
	Hcorr.h12 = -round(ymin);
	cvMultH(Hcorr, H, H_out);
}

/* Flatten a Homography into the row-major nine the record stores. */
static void flis_reg_H_to_array(const Homography *H, double out[9]) {
	out[0] = H->h00; out[1] = H->h01; out[2] = H->h02;
	out[3] = H->h10; out[4] = H->h11; out[5] = H->h12;
	out[6] = H->h20; out[7] = H->h21; out[8] = H->h22;
}

int flis_register_layers(flis_layer_t *ref_lay,
                         GSList       *target_layers,
                         registration_function method,
                         selection_type        sel_req,
                         transformation_type   tx_type,
                         opencv_interpolation interpolation,
                         gboolean       clamp,
                         flis_reg_method_id method_id) {
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

	/* The NDE record needs, per participant, the sequence slot its transform
	 * lands in and the dimensions it went IN with — register_apply_reg
	 * replaces the layer's pixels, so its input size is gone by the time we
	 * would want it (see flis_reg_compute_framed_transform). */
	struct flis_reg_part {
		flis_layer_t *lay;
		int   seq_idx;
		int   src_rx, src_ry;
	} *parts = g_new0(struct flis_reg_part, n_layers);
	int n_parts = 0;

	internal_sequence_set(seq, 0, ref->fit);
	seq->imgparam[0].rx = ref->fit->rx;
	seq->imgparam[0].ry = ref->fit->ry;
	parts[n_parts++] = (struct flis_reg_part){ ref, 0,
	                                           (int)ref->fit->rx, (int)ref->fit->ry };

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
		parts[n_parts++] = (struct flis_reg_part){ lay, i,
		                                           (int)lay->fit->rx, (int)lay->fit->ry };
		internal_sequence_set(seq, i++, lay->fit);
	}
	seq->reference_image = ref_seq_idx;
	seq->is_variable     = is_variable;

	/* Every participant needs a replay starting point BEFORE its pixels are
	 * resampled — pixels and position both, because registration moves the
	 * layer as well as warping it and a replay has nothing to anchor the
	 * move to without the starting offset (nde_checkpoint.h).  Recorded even
	 * when we go on to skip capture: a replay run may itself be the thing
	 * that first needs the baseline, and ensure() is a no-op if it exists. */
	const gboolean nde_replaying = processing_is_reserved_for_replay();
	if (!nde_replaying) {
		for (int k = 0; k < n_parts; k++)
			nde_checkpoint_baseline_ensure(parts[k].lay->fit, parts[k].lay->item_id);
	}

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
		g_free(parts);
		if (reserved) unreserve_thread();
		return 1;
	}
	free(regargs.imgparam); regargs.imgparam = NULL;

	/* The alignment method owns the reference, and it does not tell us in the
	 * field we set: register_multi_step_global scores the frames and writes
	 * its own pick to seq->reference_image (global.c), leaving
	 * regargs.reference_image on the slot WE asked for.  That field is the one
	 * compute_framing() reads to build the framing homography (applyreg.c), so
	 * when the two diverge the apply pass can be handed the H of a slot the
	 * method EXCLUDED — an all-zero matrix.  cvTransfH() inverts it to zeros,
	 * every layer is then warped by a null homography, and since a null
	 * homography maps every destination pixel onto the same source pixel each
	 * layer comes out a flat constant (the reference layer, skipped by the
	 * apply, is the only one left intact).  Re-sync before applying. */
	flis_layer_t *ref_used = ref;
	if (seq->reference_image >= 0 && seq->reference_image < n_layers &&
	    seq->reference_image != regargs.reference_image) {
		for (int k = 0; k < n_parts; k++)
			if (parts[k].seq_idx == seq->reference_image)
				ref_used = parts[k].lay;
		siril_log_message(_("flis_register_layers: the alignment method chose "
		                    "'%s' as its reference instead of '%s'\n"),
		                  ref_used->layer_name ? ref_used->layer_name : "?",
		                  ref->layer_name ? ref->layer_name : "?");
		regargs.reference_image = seq->reference_image;
	}

	/* Layers the alignment could not match are dropped from the apply pass —
	 * their pixels are left alone, which is fine, but they also get NO entry
	 * in regargs.regparam.  That array is indexed by position in the FILTERED
	 * output sequence (applyreg.c, apply_reg_finalize_hook), not by input
	 * slot, so every layer after a dropped one shifts down by one and reading
	 * it at its slot index returns another layer's offset — or runs off the
	 * end of the allocation.  Build the slot → output-index map here, while
	 * the sequence still describes what the apply is about to do. */
	int *out_of_slot = g_new(int, n_layers);
	int n_registered = 0;
	GString *unaligned = NULL;
	for (int s = 0; s < n_layers; s++) {
		gboolean ok = seq->imgparam[s].incl &&
		              guess_transform_from_H(seq->regparam[0][s].H) != NULL_TRANSFORMATION;
		out_of_slot[s] = ok ? n_registered++ : -1;
		if (ok) continue;
		for (int k = 0; k < n_parts; k++) {
			if (parts[k].seq_idx != s) continue;
			const gchar *nm = parts[k].lay->layer_name ? parts[k].lay->layer_name : "?";
			if (!unaligned) unaligned = g_string_new(nm);
			else g_string_append_printf(unaligned, ", %s", nm);
		}
	}
	if (n_registered < 2) {
		siril_log_error(_("flis_register_layers: only %d layer(s) could be "
		                  "aligned — nothing to do.  Try another reference "
		                  "layer or alignment method\n"), n_registered);
		if (unaligned) g_string_free(unaligned, TRUE);
		g_free(out_of_slot);
		free(regargs.regparam); regargs.regparam = NULL;
		free_sequence(seq, TRUE);
		g_free(parts);
		if (reserved) unreserve_thread();
		return 1;
	}
	if (unaligned) {
		siril_log_warning(_("flis_register_layers: could not align %s — "
		                    "left untouched, in place\n"), unaligned->str);
		g_string_free(unaligned, TRUE);
	}

	/* Apply transforms — resamples layer pixel data + records bounding-box
	 * offsets in regparam[i].H.h02 / h12. */
	regargs.framing = FRAMING_MAX;
	ret = register_apply_reg(&regargs);
	if (!ret && regargs.regparam) {
		/* Anchor: one layer stays put and every other layer's new offset is
		 * expressed against it.  Normally that is the canvas (base) layer,
		 * pinned back to the origin.  When the canvas layer is one of the
		 * unaligned ones its pixels were never warped, so it keeps its place
		 * and the alignment's own reference layer anchors instead. */
		const int canvas_out = (canvas_seq_idx >= 0) ? out_of_slot[canvas_seq_idx] : -1;
		double cx, cy;
		if (canvas_out >= 0) {
			cx = regargs.regparam[canvas_out].H.h02;
			cy = regargs.regparam[canvas_out].H.h12;
			canvas_lay->position_x = 0;
			canvas_lay->position_y = 0;
		} else {
			/* The anchor is the layer the alignment referenced; it is aligned
			 * by construction, and holding it at the offset it already had is
			 * what keeps the rest of the document from jumping. */
			int anchor_slot = (out_of_slot[regargs.reference_image] >= 0)
			                ? regargs.reference_image : 0;
			while (anchor_slot < n_layers && out_of_slot[anchor_slot] < 0)
				anchor_slot++;
			const flis_layer_t *anchor = ref;
			for (int k = 0; k < n_parts; k++)
				if (parts[k].seq_idx == anchor_slot)
					anchor = parts[k].lay;
			cx = regargs.regparam[out_of_slot[anchor_slot]].H.h02 -
			     (anchor == ref ? ref_orig_x : anchor->position_x);
			cy = regargs.regparam[out_of_slot[anchor_slot]].H.h12 -
			     (anchor == ref ? ref_orig_y : anchor->position_y);
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
			int out_idx = out_of_slot[seq_idx];
			if (out_idx < 0)
				continue;   /* unaligned: pixels untouched, so is the offset */
			double dx = regargs.regparam[out_idx].H.h02 - cx;
			double dy = regargs.regparam[out_idx].H.h12 - cy;
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

		/* Provenance: ONE joint record for the whole registration
		 * (nde_joint.h).  Registration's meaning is multi-layer by
		 * construction — every participant's transform is defined against
		 * the same reference — so recording it as N independent steps would
		 * both lie about the operation and make undo N steps.  Same guard as
		 * the other joint captures (flis_group_apply_channel_calibration_full):
		 * a replay must not re-record what it is reproducing, and under a
		 * python provenance scope the scope's own record subsumes this one.
		 *
		 * Only the layers that were actually warped are participants: an
		 * unaligned one has a null H, so reconstructing its framed transform
		 * would store nonsense, and a replay would then move a layer the live
		 * run left exactly where it was. */
		int n_kept = 0;
		for (int k = 0; k < n_parts; k++) {
			if (out_of_slot[parts[k].seq_idx] < 0) continue;
			if (n_kept != k) parts[n_kept] = parts[k];
			n_kept++;
		}
		n_parts = n_kept;

		if (nde_replaying) {
			/* no capture */
		} else if (nde_script_scope_active()) {
			nde_script_scope_mark_pixels_dirty();
		} else {
			struct nde_joint_register_data *p =
					nde_joint_register_data_new((guint)n_parts);
			nde_pin_spec *pins = g_new0(nde_pin_spec, n_parts);
			gchar **roles = g_new0(gchar *, n_parts);
			p->method        = (gint)method_id;
			p->tx_type       = (gint)tx_type;
			p->interpolation = (gint)interpolation;
			p->clamp         = clamp;
			p->ref_item      = ref_used->item_id;
			p->selection     = regargs.selection;
			p->canvas_w      = (gint)com.uniq->canvas_w;
			p->canvas_h      = (gint)com.uniq->canvas_h;
			for (int k = 0; k < n_parts; k++) {
				flis_layer_t *lay = parts[k].lay;
				Homography H = { 0 };
				int drx = 0, dry = 0;
				flis_reg_compute_framed_transform(
						&seq->regparam[0][parts[k].seq_idx].H,
						&regargs.framingd.Htransf,
						parts[k].src_rx, parts[k].src_ry, &H, &drx, &dry);
				flis_reg_H_to_array(&H, p->parts[k].H);
				p->parts[k].item_id = lay->item_id;
				p->parts[k].name = g_strdup(lay->layer_name ? lay->layer_name : "");
				p->parts[k].pos_x = lay->position_x;
				p->parts[k].pos_y = lay->position_y;
				/* Trust the pixels over the arithmetic: the warp has already
				 * run, so the layer's own dimensions ARE the framed output
				 * size.  The reconstruction above should agree; when it does
				 * not, storing what the warp produced keeps L1 exact. */
				p->parts[k].out_rx = lay->fit ? (gint)lay->fit->rx : drx;
				p->parts[k].out_ry = lay->fit ? (gint)lay->fit->ry : dry;
				g_free(p->parts[k].geom_sig);
				/* 0 = "everything in the log so far": the record is not
				 * appended yet, so capture time IS the end of the log. */
				p->parts[k].geom_sig =
						nde_joint_geometry_signature(lay->item_id, 0);
				roles[k] = g_strdup_printf("in%d", k);
				pins[k].role = roles[k];
				pins[k].item_id = lay->item_id;
				pins[k].record_id =
						nde_history_last_record_for_item(lay->item_id);
			}
			gchar *summary = g_strdup_printf("%s (%d %s)", _("Register layers"),
			                                 n_parts, _("layers"));
			gint64 rid = nde_capture_from_descriptor_pinned(&op_desc_flis_register,
					p, summary, NULL, ref_used->item_id, pins, (guint)n_parts);
			for (int k = 0; k < n_parts; k++)
				g_free(roles[k]);
			g_free(roles);
			g_free(pins);
			g_free(summary);
			nde_joint_register_data_free(p);
			/* Couple the undo entry the caller pushed ("Register layers") to
			 * the record, so undo and history stay one step. */
			if (rid > 0 && !com.script)
				undo_tag_top_nde_record(rid);
		}
	}

	free(regargs.imgparam); regargs.imgparam = NULL;
	free(regargs.regparam); regargs.regparam = NULL;
	free_sequence(seq, TRUE);
	g_free(out_of_slot);
	g_free(parts);
	if (reserved) unreserve_thread();
	return ret;
}

/* ======================================================================= */
/* Offline re-solve (the L2 half of the flis.register NDE replay)          */
/* ======================================================================= */

int flis_register_solve(fits **fits_in, int n, int ref_index,
                        registration_function method,
                        selection_type sel_req,
                        transformation_type tx_type,
                        rectangle selection,
                        int canvas_index,
                        int ref_pos_x, int ref_pos_y,
                        flis_reg_solution *out) {
	if (!fits_in || !out || n < 2 || ref_index < 0 || ref_index >= n || !method)
		return 1;
	for (int k = 0; k < n; k++)
		if (!fits_in[k])
			return 1;

	/* A record's selection and its method can disagree, and by the time the
	 * disagreement reaches OpenCV it is an abort, not an error return: the
	 * amend dialog lets the method be changed, and a record whose original
	 * method needed no selection stored an empty one, so switching it to
	 * KOMBAT or DFT hands cv::matchTemplate a 0x0 template.  Refusing here
	 * turns that into the L3 fallback — the stored transforms are replayed
	 * and the user is told why. */
	{
		gchar *why = NULL;
		if (!flis_register_selection_ok(sel_req, &selection,
		                                (int)fits_in[ref_index]->rx,
		                                (int)fits_in[ref_index]->ry, &why)) {
			siril_log_error(_("Register layers: cannot solve the alignment — "
			                  "%s\n"), why);
			g_free(why);
			return 1;
		}
	}

	/* register_apply_reg resamples IN PLACE, and the pixels it would be
	 * resampling are the caller's replay scratch — which the caller still
	 * has to warp itself, with the transform this solve is here to produce.
	 * So the sequence runs over copies and the originals are never touched. */
	fits **copies = g_new0(fits *, n);
	sequence *seq = create_internal_sequence(n);
	int retval = 1;
	if (!seq)
		goto out_free;

	/* Slot 0 is the reference, exactly as the live path builds it: the
	 * registration methods take the reference from seq->reference_image but
	 * the framing arithmetic is written around slot ordering, and matching
	 * the live layout is the only way the two paths can agree. */
	int *slot_of = g_new0(int, n);
	int slot = 1;
	for (int k = 0; k < n; k++)
		slot_of[k] = (k == ref_index) ? 0 : slot++;

	gboolean is_variable = FALSE;
	for (int k = 0; k < n; k++) {
		copies[k] = calloc(1, sizeof(fits));
		if (!copies[k] || copyfits(fits_in[k], copies[k], CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
			g_free(slot_of);
			goto out_free;
		}
		seq->imgparam[slot_of[k]].rx = copies[k]->rx;
		seq->imgparam[slot_of[k]].ry = copies[k]->ry;
		if (copies[k]->rx != fits_in[ref_index]->rx ||
		    copies[k]->ry != fits_in[ref_index]->ry)
			is_variable = TRUE;
		internal_sequence_set(seq, slot_of[k], copies[k]);
	}
	seq->bitpix          = fits_in[ref_index]->bitpix;
	seq->rx              = fits_in[ref_index]->rx;
	seq->ry              = fits_in[ref_index]->ry;
	seq->reference_image = 0;
	seq->is_variable     = is_variable;

	struct registration_args regargs = { 0 };
	regargs.seq                  = seq;
	regargs.layer                = 0;
	regargs.reference_image      = 0;
	regargs.run_in_thread        = FALSE;
	regargs.interpolation        = OPENCV_LANCZOS4;
	regargs.output_scale         = 1.f;
	regargs.clamp                = TRUE;
	regargs.framing              = FRAMING_MAX;
	regargs.type                 = tx_type;
	regargs.max_stars_candidates = MAX_STARS_FITTED;
	regargs.two_pass             = (method == register_multi_step_global);
	regargs.percent_moved        = 0.50f;
	regargs.prefix               = "flis_";
	if (sel_req != REQUIRES_NO_SELECTION)
		regargs.selection = selection;

	/* The input dimensions each transform was solved against; the apply pass
	 * below replaces the copies' pixels and with them their sizes. */
	int *src_rx = g_new0(int, n), *src_ry = g_new0(int, n);
	for (int k = 0; k < n; k++) {
		src_rx[k] = (int)copies[k]->rx;
		src_ry[k] = (int)copies[k]->ry;
	}

	if (method(&regargs)) {
		free(regargs.imgparam); regargs.imgparam = NULL;
		free(regargs.regparam); regargs.regparam = NULL;
		g_free(slot_of); g_free(src_rx); g_free(src_ry);
		goto out_free;
	}
	free(regargs.imgparam); regargs.imgparam = NULL;

	/* Same two traps as the live path (see flis_register_layers): the method
	 * may have picked its own reference without updating the field
	 * compute_framing() reads, and any frame it could not match is dropped
	 * from the apply — which both shifts regargs.regparam's indexing off the
	 * input slots and leaves a null matrix behind.  A solve that cannot cover
	 * every input is no solve at all here: returning failure hands the replay
	 * back to its L3 fallback, which reproduces the recorded transforms. */
	if (seq->reference_image >= 0 && seq->reference_image < n)
		regargs.reference_image = seq->reference_image;
	for (int k = 0; k < n; k++) {
		if (seq->imgparam[slot_of[k]].incl &&
		    guess_transform_from_H(seq->regparam[0][slot_of[k]].H) != NULL_TRANSFORMATION)
			continue;
		siril_log_error(_("Register layers: cannot solve the alignment — "
		                  "layer %d could not be matched to the reference\n"), k + 1);
		free(regargs.regparam); regargs.regparam = NULL;
		g_free(slot_of); g_free(src_rx); g_free(src_ry);
		goto out_free;
	}

	regargs.framing = FRAMING_MAX;
	int ret = register_apply_reg(&regargs);
	if (!ret && regargs.regparam) {
		double cx, cy;
		if (canvas_index >= 0) {
			cx = regargs.regparam[slot_of[canvas_index]].H.h02;
			cy = regargs.regparam[slot_of[canvas_index]].H.h12;
		} else {
			cx = regargs.regparam[0].H.h02 - ref_pos_x;
			cy = regargs.regparam[0].H.h12 - ref_pos_y;
		}
		for (int k = 0; k < n; k++) {
			Homography H = { 0 };
			int drx = 0, dry = 0;
			flis_reg_compute_framed_transform(&seq->regparam[0][slot_of[k]].H,
			                                  &regargs.framingd.Htransf,
			                                  src_rx[k], src_ry[k],
			                                  &H, &drx, &dry);
			flis_reg_H_to_array(&H, out[k].H);
			out[k].out_rx = copies[k] ? (int)copies[k]->rx : drx;
			out[k].out_ry = copies[k] ? (int)copies[k]->ry : dry;
			if (k == canvas_index) {
				out[k].pos_x = 0;
				out[k].pos_y = 0;
			} else {
				out[k].pos_x = (int)round(regargs.regparam[slot_of[k]].H.h02 - cx);
				out[k].pos_y = (int)round(regargs.regparam[slot_of[k]].H.h12 - cy);
			}
		}
		retval = 0;
	}
	free(regargs.imgparam); regargs.imgparam = NULL;
	free(regargs.regparam); regargs.regparam = NULL;
	g_free(slot_of); g_free(src_rx); g_free(src_ry);

out_free:
	/* free_sequence() deliberately leaves internal_fits CONTENTS alone (the
	 * live path hands it live LAYER pixels), so the copies are ours. */
	if (seq)
		free_sequence(seq, TRUE);
	for (int k = 0; k < n; k++)
		if (copies[k]) { clearfits(copies[k]); free(copies[k]); }
	g_free(copies);
	return retval;
}
