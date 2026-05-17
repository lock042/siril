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
#include "registration/mpp_ap.h"        /* mpp_aps_t — full type for aps->count */
#include "registration/mpp_config.h"
#include "registration/mpp_sidecar.h"

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

	siril_log_message(_("Stack (mpp): %d frames, %dx%d, %d APs, bitdepth=%d%s\n"),
	                  run->num_frames, run->frame_cols, run->frame_rows,
	                  run->aps->count, run->bitdepth,
	                  run->cfg->drizzle_scale > 1.001 ? " (drizzle)" : "");

	fits stacked = { 0 };
	rc = mpp_stack_apply(args->seq, run->cfg, run, &stacked);
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
