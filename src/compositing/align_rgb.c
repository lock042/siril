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

/* This file is currently not used by compositing, only by the RGB align menu
 * entry in the RGB image popup and by the rgbalign command. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/gui_iface.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "registration/registration.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "compositing/align_rgb.h"
#include "core/op_descriptors.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/undo.h"

#define REGLAYER 0

/* the largest selection accepted by the one-star method: it must contain a
 * single star for the PSF fit to be meaningful */
#define MAX_PSF_SELECTION 300

static sequence *seq = NULL;		// the sequence of channels
static struct registration_method *reg_methods[5];


static void initialize_methods() {
	reg_methods[0] = new_reg_method(_("One star registration (deep-sky)"),
			&register_shift_fwhm, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[1] = new_reg_method(_("Image pattern alignment (planetary/deep-sky)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
	reg_methods[2] = new_reg_method(_("Global star registration (deep-sky)"),
			&register_multi_step_global, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[3] = new_reg_method(_("KOMBAT registration (planetary / deep-sky)"),
			&register_kombat, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[4] = NULL;
}

const char *rgb_align_method_name(rgb_align_method m) {
	switch (m) {
		case RGBALIGN_PSF:
			return _("one star");
		case RGBALIGN_DFT:
			return _("image pattern");
		case RGBALIGN_KOMBAT:
			return _("KOMBAT");
		default:
			return _("global star");
	}
}

/* Checks that the loaded image and the current selection are suitable for the
 * requested alignment method, logging the reason if they are not. Shared by the
 * RGB align menu entries and the rgbalign command. */
gboolean rgb_align_prerequisites_met(rgb_align_method m) {
	if (!single_image_is_loaded() || !isrgb(gfit)) {
		siril_log_message(_("RGB alignment requires a loaded colour image.\n"));
		return FALSE;
	}
	/* the global star method registers on the whole image, the others all
	 * need an area to work on */
	if (m != RGBALIGN_GLOBAL && (com.selection.w <= 0 || com.selection.h <= 0)) {
		siril_log_message(_("The %s alignment method requires a selection. Make one in the "
					"image or use the boxselect command.\n"), rgb_align_method_name(m));
		return FALSE;
	}
	if (m == RGBALIGN_PSF && (com.selection.w > MAX_PSF_SELECTION || com.selection.h > MAX_PSF_SELECTION)) {
		siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
		return FALSE;
	}
	return TRUE;
}

// We cannot currently do this in free_sequence() because compositing still
// uses the references, so we have to do it here as a special case
static void free_internal_sequence(sequence *seq) {
	if (seq) {
		for (int i = 0; i < seq->number; i++)
			clearfits(internal_sequence_get(seq, i));
		free_sequence(seq, TRUE);
	}
	clear_stars_list(TRUE);
}

static int initialize_internal_rgb_sequence(fits *source) {
	if (seq) free_internal_sequence(seq);

	seq = create_internal_sequence(3);
	for (int i = 0; i < 3; i++) {
		fits *fit = calloc(1, sizeof(fits));
		if (extract_fits(source, fit, i, FALSE)) {
			free(fit);
			free_sequence(seq, TRUE);
			return -1;
		}
		internal_sequence_set(seq, i, fit);
	}
	seq->rx = source->rx;
	seq->ry = source->ry;
	seq->bitpix = source->bitpix;

	return 0;
}

static void compose(fits *target) {
	size_t npixels = target->rx * target->ry;
	fits *fit[3];
	/* Channels the registration failed on are excluded from the sequence and
	 * so keep their original data: the result is only partially aligned, which
	 * is easy to miss in a script if we don't say so. */
	for (int i = 0; i < 3; i++) {
		if (!seq->imgparam[i].incl)
			siril_log_warning(_("The %s channel could not be aligned and is left unchanged.\n"),
					channel_number_to_name(i));
	}
	for (int i = 0 ; i < 3 ; i++) {
		fit[i] = internal_sequence_get(seq, i);
	}
	if (target->type == DATA_FLOAT) {
		for (int i = 0 ; i < 3 ; i++) {
			memcpy(target->fpdata[i], fit[i]->fdata, sizeof(float) * npixels);
		}
	} else {
		for (int i = 0 ; i < 3 ; i++) {
			memcpy(target->pdata[i], fit[i]->data, sizeof(WORD) * npixels);
		}
	}
}

/* The alignment proper, against whichever image it is handed.
 *
 * @area, when non-NULL, is used verbatim instead of asking
 * get_the_registration_area() — which reads com.selection, a live GUI state
 * with no meaning at replay time.  @area_used, when non-NULL, receives
 * whichever area was in force, so a capture can record it.
 *
 * @claim_thread says whether to take the processing slot.  The registration
 * methods call generic_sequence_worker synchronously and poll
 * processing_should_continue(), so the slot must be held and cancel_flag
 * cleared for the duration — but reserve_thread() is a plain compare-and-swap,
 * not a recursive lock, so only the caller that is NOT already inside a worker
 * may take it.  The menu and the command are; replay, arriving through
 * generic_image_worker, is not. */
static int rgb_align_core(fits *target, rgb_align_method m, const rectangle *area,
                          rectangle *area_used, gboolean claim_thread) {
	struct registration_args regargs = { 0 };
	struct registration_method *method;
	framing_type framing = FRAMING_COG;
	int retval1 = 0, retval2 = 0;

	initialize_methods();
	if (initialize_internal_rgb_sequence(target)) {
		siril_log_message(_("Could not extract the channels of the loaded image.\n"));
		return 1;
	}
	gui_iface.set_progress(PROGRESS_RESET, NULL);

	/* align it */
	method = reg_methods[m];
	regargs.seq = seq;
	regargs.no_output = FALSE;
	if (area)
		regargs.selection = *area;
	else
		get_the_registration_area(&regargs, method);
	if (area_used)
		*area_used = regargs.selection;
	regargs.layer = REGLAYER;
	seq->reference_image = 0;
	regargs.seq->nb_layers = 1;
	regargs.max_stars_candidates = MAX_STARS_FITTED;
	regargs.run_in_thread = FALSE;
	regargs.interpolation = OPENCV_LANCZOS4;
	regargs.clamp = TRUE;
	regargs.framing = framing;
	regargs.output_scale = 1.f;
	regargs.percent_moved = 0.50f; // Only needed for KOMBAT
	regargs.two_pass = TRUE;
	if (method->method_ptr == register_shift_fwhm || method->method_ptr == register_shift_dft)
		regargs.type = SHIFT_TRANSFORMATION;
	else
		regargs.type = HOMOGRAPHY_TRANSFORMATION;

	if (claim_thread && !reserve_thread()) {
		siril_log_message(_("A processing operation is already running.\n"));
		free_internal_sequence(seq);
		seq = NULL;
		return 1;
	}
	retval1 = method->method_ptr(&regargs);
	free(regargs.imgparam);
	regargs.imgparam = NULL;
	free(regargs.regparam);
	regargs.regparam = NULL;
	if (retval1) {
		gui_iface.set_progress(PROGRESS_DONE, _("Error in channels alignment."));
		if (claim_thread)
			unreserve_thread();
		free_internal_sequence(seq);
		seq = NULL;
		return retval1;
	}
	retval2 = register_apply_reg(&regargs);
	compose(target);
	if (claim_thread)
		unreserve_thread();

	if (retval2)
		gui_iface.set_progress(PROGRESS_DONE, _("Error in layers alignment."));
	else {
		gui_iface.set_progress(PROGRESS_DONE, _("Registration complete."));
		siril_log_message(_("Aligned RGB channels using the %s method\n"), rgb_align_method_name(m));
	}
	free_internal_sequence(seq);
	seq =  NULL;
	return retval2;
}

int rgb_align(rgb_align_method m, rectangle *area_used) {
	gui_iface.set_busy(TRUE);
	int retval = rgb_align_core(gfit, m, NULL, area_used, TRUE);
	if (!retval) {
		notify_gfit_data_modified();
		gfit_modified_update_gui();
	}
	gui_iface.set_busy(FALSE);
	return retval;
}

/* RGB channel alignment shares a direct-apply path (no generic_image_worker),
 * so the caller's undo_save_state is the sole commit point.  The save is
 * taken before rgb_align() (which has failure/early-return paths), so
 * provenance is recorded ONLY once rgb_align returns success (0), tagging the
 * entry saved just above (@undo_err gates the tag).
 *
 * The record is Tier A: the method is the user's choice and the registration
 * area is whatever com.selection resolved to, which rgb_align reports back
 * because by replay time the selection is long gone. */
int rgb_align_and_capture(rgb_align_method m, int undo_err, const char *summary) {
	gint target = nde_checkpoint_active_item_id();
	/* NDE baseline (phase 2): snapshot gfit's pre-op pixels BEFORE the
	 * direct-apply mutation.  A stray baseline for an item with no live
	 * record (rgb_align failure) is never persisted — the saver/loader gate
	 * on live records — so ensuring first is safe. */
	nde_checkpoint_baseline_ensure(gfit, target);
	rectangle area_used = { 0 };
	int retval = rgb_align(m, &area_used);
	if (retval)
		return retval;   /* failure: pixels unchanged, no provenance */
	struct rgb_align_data *p = new_rgb_align_data();
	if (!p)
		return 0;
	p->method = m;
	p->area = area_used;
	p->have_area = TRUE;
	gint64 rid = nde_capture_from_descriptor(&op_desc_rgb_align, p, summary,
			gfit, FALSE);
	free_rgb_align_data(p);
	if (!undo_err)
		undo_tag_top_nde_record(rid);
	return 0;
}

/* Replay path: the worker hands us a private fits and already holds the
 * processing slot, and the recorded area stands in for the live selection. */
int rgb_align_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	(void)nb_threads;
	const struct rgb_align_data *p = args->user;
	if (!p)
		return 1;
	return rgb_align_core(fit, p->method, p->have_area ? &p->area : NULL,
	                      NULL, FALSE);
}

/* ---- params struct + NDE serializers ------------------------------------ */

void free_rgb_align_data(void *p) {
	free(p);
}

struct rgb_align_data *new_rgb_align_data(void) {
	struct rgb_align_data *p = calloc(1, sizeof(struct rgb_align_data));
	if (p)
		p->destroy_fn = free_rgb_align_data;
	return p;
}

static gchar *rgb_align_log_hook(gpointer p, log_hook_detail detail) {
	const struct rgb_align_data *d = p;
	(void)detail;
	if (!d)
		return g_strdup(_("RGB alignment"));
	return g_strdup_printf(_("RGB alignment (%s)"),
			rgb_align_method_name((rgb_align_method)d->method));
}

static gchar *rgb_align_serialize(gconstpointer user) {
	const struct rgb_align_data *p = user;
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "method", p->method);
	if (p->have_area) {
		nde_kv_add_int(kv, "sel_x", p->area.x);
		nde_kv_add_int(kv, "sel_y", p->area.y);
		nde_kv_add_int(kv, "sel_w", p->area.w);
		nde_kv_add_int(kv, "sel_h", p->area.h);
	}
	return nde_kv_end(kv);
}

static gpointer rgb_align_deserialize(const gchar *blob, int version) {
	if (version > op_desc_rgb_align.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	gint64 method = 0;
	if (!nde_kv_get_int(kv, "method", &method) || method < 0 || method > 3) {
		g_hash_table_unref(kv);
		return NULL;
	}
	gint64 x, y, w, h;
	gboolean have_area = nde_kv_get_int(kv, "sel_x", &x) &&
	                     nde_kv_get_int(kv, "sel_y", &y) &&
	                     nde_kv_get_int(kv, "sel_w", &w) &&
	                     nde_kv_get_int(kv, "sel_h", &h);
	struct rgb_align_data *p = new_rgb_align_data();
	if (p) {
		p->method = (int)method;
		p->have_area = have_area;
		if (have_area) {
			p->area.x = (int)x; p->area.y = (int)y;
			p->area.w = (int)w; p->area.h = (int)h;
		}
	}
	g_hash_table_unref(kv);
	return p;
}

/* Op descriptor — single source of truth for this operation (op_descriptor.h) */
const op_descriptor op_desc_rgb_align = {
	.id = "color.rgb_align", .version = 1,
	.image_hook = rgb_align_image_hook,
	.log_hook = rgb_align_log_hook,
	.description = N_("RGB Alignment"),
	.mem_ratio = 3.0f,
	.flags = 0,
	.serialize = rgb_align_serialize, .deserialize = rgb_align_deserialize,
};
