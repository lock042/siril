/*
 * AP editor dialog handlers — Phase 9 slice 9.
 *
 * Edits the cached MPP run (com.mpp_run) in place. Dialog is non-modal so
 * the user can interact with the image. Slice 9 ships Auto-place + Clear
 * (which mutate run->aps via mpp_ap_replace / mpp_ap_clear_all);
 * slice 10 adds mouse interactivity (left=add, right=remove, drag=move)
 * and the snapshot/revert plumbing for Cancel.
 *
 * Spinners share GtkAdjustments with the MPP register sub-panel — values
 * stay in sync automatically. The editor reads them at Auto-place time.
 */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/mpp_ap_editor.h"
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"

static GtkWidget *dialog = NULL;
static GtkLabel  *count_label = NULL;
static GtkSpinButton *spin_half_box = NULL;
static GtkSpinButton *spin_brightness = NULL;
static GtkSpinButton *spin_contrast = NULL;
static GtkSpinButton *spin_structure = NULL;

/* AP currently being dragged by left-mouse-down; -1 when not dragging.
 * Set by main_action_click on a hit; read by motion_notify; cleared by
 * the release callback. */
static int g_ap_drag_idx = -1;

/* AP currently under the cursor for hover highlighting; -1 if none.
 * Updated by motion_notify while in MOUSE_ACTION_EDIT_APS, cleared on
 * dialog hide. */
static int g_ap_hover_idx = -1;

/* Snapshot of the AP grid taken when the dialog opens. Used by Cancel to
 * revert edits made between show and Cancel. Commit drops it. */
static mpp_aps_t *g_aps_snapshot = NULL;

static void editor_init_statics(void) {
	if (dialog) return;
	dialog          = GTK_WIDGET     (gtk_builder_get_object(gui.builder, "mpp_ap_editor_dialog"));
	count_label     = GTK_LABEL      (gtk_builder_get_object(gui.builder, "label_mpp_ap_editor_count"));
	spin_half_box   = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_half_box"));
	spin_brightness = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_brightness"));
	spin_contrast   = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_contrast"));
	spin_structure  = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_structure"));
}

/* Refresh the "Current APs: N" label from the cached run. */
static void editor_refresh_count_label(void) {
	if (!count_label) return;
	const mpp_run_t *run = mpp_get_cached_run();
	const int n = (run && run->aps) ? run->aps->count : 0;
	gchar *txt = g_strdup_printf(_("Current APs: %d"), n);
	gtk_label_set_text(count_label, txt);
	g_free(txt);
}

/* Public hook for mouse handlers: refresh the count label after add/remove
 * from mouse_action_functions.c. No-op if dialog hasn't been shown yet. */
void mpp_ap_editor_refresh_count_label(void) {
	editor_refresh_count_label();
}

int  mpp_ap_editor_get_drag_idx(void)      { return g_ap_drag_idx; }
void mpp_ap_editor_set_drag_idx(int idx)   { g_ap_drag_idx = idx; }

int  mpp_ap_editor_get_hover_idx(void)     { return g_ap_hover_idx; }
void mpp_ap_editor_set_hover_idx(int idx)  { g_ap_hover_idx = idx; }

gboolean mpp_ap_editor_is_open(void) {
	return dialog && gtk_widget_get_visible(dialog);
}

/* Open the editor dialog from the "Edit APs…" button on the MPP sub-panel. */
void on_seqmpp_edit_aps_button_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (!mpp_get_cached_run()) {
		siril_log_color_message(_("AP editor: no analysis result is cached. "
		                          "Run Analyze first.\n"), "salmon");
		return;
	}
	editor_init_statics();
	if (!dialog) return;
	gtk_widget_show(dialog);
}

/* GtkDialog "delete-event" -> revert + hide. Closing the window with the
 * X (or Escape if/when bound) is semantically a Cancel — discard edits
 * made since the dialog was opened. */
gboolean on_mpp_ap_editor_dialog_delete_event(GtkWidget *widget, GdkEvent *event,
                                              gpointer user_data) {
	(void) event; (void) user_data;
	mpp_run_t *run = mpp_get_cached_run();
	if (run && g_aps_snapshot) {
		mpp_aps_restore(run, g_aps_snapshot);
		g_aps_snapshot = NULL;
		siril_log_message(_("AP editor: edits reverted (window closed).\n"));
		redraw(REDRAW_OVERLAY);
	}
	gtk_widget_hide(widget);
	return TRUE;   /* swallow the delete so the toplevel isn't destroyed */
}

void on_mpp_ap_editor_dialog_show(GtkWidget *widget, gpointer user_data) {
	(void) widget; (void) user_data;
	editor_init_statics();
	editor_refresh_count_label();
	/* Snapshot APs so Cancel can revert. */
	const mpp_run_t *run = mpp_get_cached_run();
	if (run && run->aps) {
		if (g_aps_snapshot) mpp_ap_free(g_aps_snapshot);
		g_aps_snapshot = mpp_aps_snapshot(run->aps);
	}
	/* Flip mouse mode so left=add, right=remove, drag=move (handlers in
	 * mouse_action_functions.c). */
	mouse_status = MOUSE_ACTION_EDIT_APS;
	siril_log_message(_("AP editor: left-click to add, right-click to remove, "
	                    "drag to move an AP.\n"));
}

void on_mpp_ap_editor_dialog_hide(GtkWidget *widget, gpointer user_data) {
	(void) widget; (void) user_data;
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	g_ap_drag_idx = -1;
	g_ap_hover_idx = -1;
	redraw(REDRAW_OVERLAY);   /* clear hover highlight */
	/* Don't auto-free the snapshot here — Cancel and Commit handle it
	 * explicitly; this path can run via window-manager close, in which
	 * case we discard (same as Cancel). */
	if (g_aps_snapshot) {
		mpp_ap_free(g_aps_snapshot);
		g_aps_snapshot = NULL;
	}
}

/* Read AP-placement params from the dialog's spinners into a config snapshot.
 * Other cfg fields (search_width, etc.) are taken from the cached run. */
static void cfg_from_spinners(const mpp_run_t *run, mpp_config_t *out) {
	*out = *run->cfg;
	if (spin_half_box)
		out->alignment_points_half_box_width   = gtk_spin_button_get_value_as_int(spin_half_box);
	if (spin_brightness)
		out->alignment_points_brightness_threshold = gtk_spin_button_get_value_as_int(spin_brightness);
	if (spin_contrast)
		out->alignment_points_contrast_threshold = gtk_spin_button_get_value_as_int(spin_contrast);
	if (spin_structure)
		out->alignment_points_structure_threshold = gtk_spin_button_get_value(spin_structure);
}

void on_mpp_ap_editor_auto_place_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	mpp_run_t *run = mpp_get_cached_run();
	if (!run) {
		siril_log_color_message(_("AP editor: no cached run.\n"), "red");
		return;
	}
	mpp_config_t cfg;
	cfg_from_spinners(run, &cfg);
	const int rc = mpp_ap_replace(run, &cfg);
	switch (rc) {
		case MPP_OK:
			siril_log_color_message(_("AP editor: auto-placed %d APs. "
			                          "Re-Analyze before Register to refresh per-AP qualities.\n"),
			                        "green", run->aps->count);
			break;
		case MPP_ENODATA:
			siril_log_color_message(_("AP editor: zero APs survived the threshold filters. "
			                          "Loosen min-brightness/contrast/structure and try again.\n"),
			                        "red");
			break;
		case MPP_ENOMEM:
			siril_log_color_message(_("AP editor: out of memory.\n"), "red");
			break;
		default:
			siril_log_color_message(_("AP editor: auto-place failed (code %d).\n"), "red", rc);
			break;
	}
	editor_refresh_count_label();
	redraw(REDRAW_OVERLAY);
}

void on_mpp_ap_editor_clear_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	mpp_run_t *run = mpp_get_cached_run();
	if (!run) return;
	mpp_ap_clear_all(run);
	siril_log_color_message(_("AP editor: cleared all APs.\n"), "salmon");
	editor_refresh_count_label();
	redraw(REDRAW_OVERLAY);
}

void on_mpp_ap_editor_commit_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	/* Drop the snapshot (no revert) before hiding so the show/hide
	 * snapshot logic doesn't restore it on a subsequent open. */
	if (g_aps_snapshot) {
		mpp_ap_free(g_aps_snapshot);
		g_aps_snapshot = NULL;
	}
	if (dialog) gtk_widget_hide(dialog);
	const mpp_run_t *run = mpp_get_cached_run();
	if (run && run->aps && !run->best_frame_indices) {
		siril_log_color_message(_("AP editor: edits committed. Click Analyze again to refresh "
		                          "per-AP frame qualities, then Register.\n"), "blue");
	}
}

void on_mpp_ap_editor_cancel_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	mpp_run_t *run = mpp_get_cached_run();
	if (run && g_aps_snapshot) {
		mpp_aps_restore(run, g_aps_snapshot);
		g_aps_snapshot = NULL;   /* ownership transferred */
		siril_log_message(_("AP editor: edits reverted.\n"));
		redraw(REDRAW_OVERLAY);
	}
	if (dialog) gtk_widget_hide(dialog);
}
