/*
 * STACK_MPP method handler — Stage C of the PSS multipoint pipeline.
 *
 * Reads the <seqname>.mpp sidecar written by REG_MPP (or the CLI
 * register_mpp command), applies any stack-side widget overrides from
 * args->mpp_cfg, runs mpp_stack_apply, and hands the result back to the
 * stacking framework via args->result. The rest of the framework
 * (gfit copy, savefits, GUI refresh) is unchanged.
 */

#include <glib.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "core/gui_iface.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"

#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"        /* mpp_aps_t — full type for aps->count */
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_sidecar.h"

#include "stacking/stacking.h"

int stack_mpp_handler(struct stacking_args *args) {
	if (!args || !args->seq) {
		if (args) args->retval = ST_GENERIC_ERROR;
		return ST_GENERIC_ERROR;
	}

	siril_log_status(_("Stack (mpp): reading sidecar\n"));
	gui_iface.set_progress(PROGRESS_RESET, _("Stack (mpp): reading sidecar"));

	gchar *sidecar_path = g_strdup_printf("%s.mpp", args->seq->seqname);
	mpp_run_t *run = NULL;
	int rc = mpp_sidecar_read(sidecar_path, &run);
	if (rc != MPP_OK || !run) {
		siril_log_error(_("Stack (mpp): cannot read sidecar %s (code %d) — "
		                          "run \"Multipoint Registration\" from the registration "
		                          "tab first.\n"), sidecar_path, rc);
		g_free(sidecar_path);
		args->retval = ST_GENERIC_ERROR;
		gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
		return ST_GENERIC_ERROR;
	}
	g_free(sidecar_path);

	/* Apply stack-side widget overrides on top of the cfg the sidecar
	 * persisted at register time. Only stack-side fields are touched;
	 * register-side decisions (AP grid, per-AP shifts) are already baked
	 * into the sidecar and not adjustable here. */
	if (args->mpp_cfg) {
		const mpp_config_t *gui = (const mpp_config_t *) args->mpp_cfg;
		run->cfg->drizzle_scale                           = gui->drizzle_scale;
		run->cfg->drizzle_mode                            = gui->drizzle_mode;
		run->cfg->drizzle_pixfrac                         = gui->drizzle_pixfrac;
		run->cfg->drizzle_kernel                          = gui->drizzle_kernel;
		run->cfg->alignment_points_frame_percent          = gui->alignment_points_frame_percent;
		run->cfg->alignment_points_frame_number           = gui->alignment_points_frame_number;
		run->cfg->stack_frames_background_fraction        = gui->stack_frames_background_fraction;
		run->cfg->stack_frames_background_blend_threshold = gui->stack_frames_background_blend_threshold;
	}

	/* Apply the generic Stack-tab filter (all / selected / criteria) on
	 * top of run->included. The framework hands us image_indices listing
	 * the frames that passed args->filtering_criterion; everything else
	 * gets excluded. AND with whatever run->included already had (Stage A
	 * fills it with 1s, sidecar can persist user exclusions). When the
	 * generic filter narrows the set AND the MPP per-AP %/num control is
	 * also narrower than all-frames, emit a one-line note so the user
	 * sees how the two stack. */
	int generic_filter_active = 0;
	if (args->image_indices && args->nb_images_to_stack > 0
	 && args->nb_images_to_stack < run->num_frames) {
		int *new_inc = calloc((size_t) run->num_frames, sizeof(int));
		if (new_inc) {
			for (int i = 0; i < args->nb_images_to_stack; ++i) {
				const int idx = args->image_indices[i];
				if (idx >= 0 && idx < run->num_frames)
					new_inc[idx] = 1;
			}
			for (int i = 0; i < run->num_frames; ++i)
				run->included[i] = new_inc[i] && run->included[i];
			free(new_inc);
			generic_filter_active = 1;
		}
	}

	const int mpp_filter_narrows =
	    (run->cfg->alignment_points_frame_number > 0
	     && run->cfg->alignment_points_frame_number < run->num_frames)
	 || (run->cfg->alignment_points_frame_number <= 0
	     && run->cfg->alignment_points_frame_percent < 100);
	if (generic_filter_active && mpp_filter_narrows) {
		siril_log_warning(_("Stack (mpp): both filters active — Siril's "
		                    "selection narrows the eligible set to %d/%d "
		                    "frames; the MPP per-AP control then picks the "
		                    "top quality slice within that set. If you only "
		                    "want one filter, set the other to all-frames.\n"),
		                  args->nb_images_to_stack, run->num_frames);
	}

	siril_log_message(_("Stack (mpp): %d frames, %dx%d, %d APs, bitdepth=%d%s\n"),
	                  run->num_frames, run->frame_cols, run->frame_rows,
	                  run->aps->count, run->bitdepth,
	                  run->cfg->drizzle_scale > 1.001 ? " (drizzle)" : "");

	fits stacked = { 0 };
	rc = mpp_stack_apply(args->seq, run->cfg, run, &stacked);
	if (rc == MPP_EINTR) {
		siril_log_message(_("Stack (mpp): cancelled by user.\n"));
		clearfits(&stacked);
		mpp_run_free(run);
		args->retval = ST_GENERIC_ERROR;
		gui_iface.set_progress(PROGRESS_DONE, _("Cancelled"));
		return ST_GENERIC_ERROR;
	}
	if (rc != MPP_OK) {
		siril_log_error(_("Stack (mpp): Stage C failed (code %d)\n"), rc);
		clearfits(&stacked);
		mpp_run_free(run);
		args->retval = ST_GENERIC_ERROR;
		gui_iface.set_progress(PROGRESS_DONE, _("Failed"));
		return ST_GENERIC_ERROR;
	}
	mpp_run_free(run);

	/* Shallow-copy the stacked fits into args->result; ownership of
	 * `data` transfers to args->result. stack_function_handler will
	 * then memcpy it into gfit (which clears the source). */
	memcpy(&args->result, &stacked, sizeof(fits));
	args->retval = ST_OK;
	gui_iface.set_progress(PROGRESS_DONE, _("Stack (mpp): done"));
	return ST_OK;
}
