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
#include <gdk/gdkkeysyms.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/mpp_ap_editor.h"
#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_config.h"

static GtkWidget *dialog = NULL;
static GtkLabel  *count_label = NULL;
static GtkSpinButton *spin_half_box = NULL;
static GtkSpinButton *spin_brightness = NULL;
static GtkSpinButton *spin_contrast = NULL;
static GtkSpinButton *spin_structure = NULL;
/* "Show stacking patches" toggle. State is read by draw_mpp_aps as a
 * persistent flag: outlines stay drawn whenever the AP overlay is
 * visible (Registration tab, matching dims) even after the dialog
 * closes, until the user reopens the editor and unchecks. */
static GtkToggleButton *check_show_patches = NULL;

/* AP currently being dragged by left-mouse-down; -1 when not dragging.
 * Set by main_action_click on a hit; read by motion_notify; cleared by
 * the release callback. */
static int g_ap_drag_idx = -1;

/* AP currently under the cursor for hover highlighting; -1 if none.
 * Updated by motion_notify while in MOUSE_ACTION_EDIT_APS, cleared on
 * dialog hide. */
static int g_ap_hover_idx = -1;

/* AP currently selected (by clicking it); -1 if none. The resize controls
 * (mouse scroll, +/- keys) act on this AP. Cleared when the AP set changes
 * (add/remove/clear/auto-place) so the index can't go stale, and on close. */
static int g_ap_selected_idx = -1;

/* Half-box of the selected AP at the moment it was selected. Escape reverts
 * the AP to this size. Captured by mpp_ap_editor_set_selected_idx. */
static int g_ap_selected_orig_hb = -1;

/* Within-session AP undo/redo. While the editor is open the global Undo/Redo
 * (Ctrl-Z/Y, the toolbar buttons, the menu) are diverted here. Each entry is
 * a deep-copied AP grid state; the stacks are cleared when the dialog opens
 * and closes. */
static GPtrArray *g_ap_undo = NULL;   /* prior states, most recent last */
static GPtrArray *g_ap_redo = NULL;
static int g_ap_undo_coalesce = -1;   /* key of the last coalescing edit, else -1 */

/* Snapshot of the AP grid taken when the dialog opens. Used by Cancel to
 * revert edits made between show and Cancel. Commit drops it. */
static mpp_aps_t *g_aps_snapshot = NULL;

/* +/- (and =, KP +/-) resize the selected AP while the editor is open.
 * Connected to both the editor dialog and the control window so it fires
 * whichever has keyboard focus (the user typically mouses the image, i.e.
 * the control window). Inert unless the editor is open, and never steals
 * the keys from a focused entry/spin-button so numeric fields stay usable. */
static gboolean editor_key_press(GtkWidget *widget, GdkEventKey *event,
                                 gpointer user_data) {
	(void) user_data;
	if (!mpp_ap_editor_is_open()) return FALSE;
	GtkWidget *top = gtk_widget_get_toplevel(widget);
	GtkWidget *focus = GTK_IS_WINDOW(top) ? gtk_window_get_focus(GTK_WINDOW(top)) : NULL;
	if (focus && (GTK_IS_ENTRY(focus) || GTK_IS_SPIN_BUTTON(focus)))
		return FALSE;
	mpp_run_t *run = mpp_get_cached_run();
	const int sel = mpp_ap_editor_get_selected_idx();
	const gboolean have_sel = (run && run->aps && sel >= 0 && sel < run->aps->count);

	/* Escape reverts the selected AP to its size at selection (stays
	 * selected); Enter fixes the current size and deselects. With no
	 * selection both fall through so Escape can still close the dialog. */
	switch (event->keyval) {
		case GDK_KEY_Escape:
			if (have_sel) {
				if (g_ap_selected_orig_hb > 0) {
					mpp_ap_editor_record_undo(-1);
					mpp_ap_set_half_box(run, sel, g_ap_selected_orig_hb);
				}
				redraw(REDRAW_OVERLAY);
				return TRUE;
			}
			return FALSE;
		case GDK_KEY_Return:
		case GDK_KEY_KP_Enter:
			if (have_sel) {
				mpp_ap_editor_set_selected_idx(-1);
				redraw(REDRAW_OVERLAY);
				return TRUE;
			}
			return FALSE;
		default:
			break;
	}

	/* +/- resize — leave Ctrl/Alt/Super combos alone so Ctrl+plus/minus
	 * still zoom. */
	if (event->state & (GDK_CONTROL_MASK | GDK_MOD1_MASK | GDK_SUPER_MASK))
		return FALSE;
	int step = 0;
	switch (event->keyval) {
		case GDK_KEY_plus:  case GDK_KEY_equal: case GDK_KEY_KP_Add:
			step = 2; break;
		case GDK_KEY_minus: case GDK_KEY_underscore: case GDK_KEY_KP_Subtract:
			step = -2; break;
		default:
			return FALSE;
	}
	if (have_sel) {
		mpp_ap_editor_record_undo(sel);   /* coalesce the resize burst */
		mpp_ap_resize(run, sel, step);
		redraw(REDRAW_OVERLAY);
		return TRUE;   /* consumed */
	}
	return FALSE;   /* no selection — let the key propagate */
}

static void editor_init_statics(void) {
	if (dialog) return;
	dialog             = GTK_WIDGET     (gtk_builder_get_object(gui.builder, "mpp_ap_editor_dialog"));
	count_label        = GTK_LABEL      (gtk_builder_get_object(gui.builder, "label_mpp_ap_editor_count"));
	spin_half_box      = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_half_box"));
	spin_brightness    = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_brightness"));
	spin_contrast      = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_contrast"));
	spin_structure     = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_ap_editor_structure"));
	check_show_patches = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "check_mpp_ap_editor_show_patches"));
	/* Resize-by-key: connect once. The handler self-gates on
	 * mpp_ap_editor_is_open() so it's inert when the dialog is closed. */
	if (dialog) {
		/* Ensure no default widget, so Enter never activates a button — the
		 * key handler uses Enter to deselect a selected AP, and with no
		 * selection Enter should simply do nothing. */
		gtk_window_set_default(GTK_WINDOW(dialog), NULL);
		g_signal_connect(dialog, "key-press-event", G_CALLBACK(editor_key_press), NULL);
	}
	GtkWidget *ctrl = GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window"));
	if (ctrl)
		g_signal_connect(ctrl, "key-press-event", G_CALLBACK(editor_key_press), NULL);
}

/* Persistent flag — true iff the user has ticked "Show stacking
 * patches" in the AP editor (state survives dialog close). Read by
 * draw_mpp_aps to decide whether to paint the larger dashed-cyan patch
 * outlines alongside the boxes. Must lazily init the widget pointer
 * because the overlay can repaint before the user has ever opened the
 * editor (count_label etc. are filled by editor_init_statics on first
 * dialog-open, which may not have happened yet). */
gboolean mpp_ap_editor_show_patches(void) {
	if (!check_show_patches) {
		if (!gui.builder) return FALSE;
		check_show_patches = GTK_TOGGLE_BUTTON(
		    gtk_builder_get_object(gui.builder, "check_mpp_ap_editor_show_patches"));
		if (!check_show_patches) return FALSE;
	}
	return gtk_toggle_button_get_active(check_show_patches);
}

void on_mpp_ap_editor_show_patches_toggled(GtkToggleButton *toggle, gpointer user_data) {
	(void) toggle; (void) user_data;
	redraw(REDRAW_OVERLAY);
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

int  mpp_ap_editor_get_selected_idx(void)    { return g_ap_selected_idx; }
void mpp_ap_editor_set_selected_idx(int idx) {
	g_ap_selected_idx = idx;
	/* Capture the AP's current size as the Escape-revert baseline. */
	if (idx >= 0)
		g_ap_selected_orig_hb = mpp_ap_get_half_box(mpp_get_cached_run(), idx);
}

static void ap_undo_free_all(GPtrArray *a) {
	if (!a) return;
	for (guint i = 0; i < a->len; i++)
		mpp_ap_free((mpp_aps_t *) g_ptr_array_index(a, i));
	g_ptr_array_set_size(a, 0);
}

static void ap_undo_clear(void) {
	ap_undo_free_all(g_ap_undo);
	ap_undo_free_all(g_ap_redo);
	g_ap_undo_coalesce = -1;
}

gboolean mpp_ap_editor_can_undo(void) { return g_ap_undo && g_ap_undo->len > 0; }
gboolean mpp_ap_editor_can_redo(void) { return g_ap_redo && g_ap_redo->len > 0; }

/* Keep the win.undo/redo actions enabled to reflect the AP undo stacks while
 * the editor is open (otherwise they'd be disabled whenever the image has no
 * processing history, and Ctrl-Z / the toolbar buttons would never reach the
 * diversion). Restores the image-history state when the editor is closed. */
static void ap_sync_undo_actions(void) {
	if (!gui.builder) return;
	/* Editor closed: hand the Undo/Redo actions and tooltips back to the
	 * image-history state via the central updater. */
	if (!mpp_ap_editor_is_open()) {
		update_MenuItem(NULL);
		return;
	}
	const gboolean can_u = mpp_ap_editor_can_undo();
	const gboolean can_r = mpp_ap_editor_can_redo();
	GObject *w = gtk_builder_get_object(gui.builder, "control_window");
	if (w) {
		GActionMap *map = G_ACTION_MAP(w);
		GAction *u = g_action_map_lookup_action(map, "undo");
		GAction *r = g_action_map_lookup_action(map, "redo");
		if (u) g_simple_action_set_enabled(G_SIMPLE_ACTION(u), can_u);
		if (r) g_simple_action_set_enabled(G_SIMPLE_ACTION(r), can_r);
	}
	GObject *ub = gtk_builder_get_object(gui.builder, "header_undo_button");
	GObject *rb = gtk_builder_get_object(gui.builder, "header_redo_button");
	if (ub) gtk_widget_set_tooltip_text(GTK_WIDGET(ub),
	        can_u ? _("Undo alignment-point edit") : _("Nothing to undo"));
	if (rb) gtk_widget_set_tooltip_text(GTK_WIDGET(rb),
	        can_r ? _("Redo alignment-point edit") : _("Nothing to redo"));
}

void mpp_ap_editor_record_undo(int coalesce_key) {
	if (!mpp_ap_editor_is_open()) return;
	const mpp_run_t *run = mpp_get_cached_run();
	if (!run || !run->aps) return;
	if (!g_ap_undo) g_ap_undo = g_ptr_array_new();
	if (!g_ap_redo) g_ap_redo = g_ptr_array_new();
	/* Coalesce consecutive edits sharing a non-negative key (e.g. a burst of
	 * scroll/+/- resizes, or a single drag) into one undo step. */
	if (coalesce_key >= 0 && coalesce_key == g_ap_undo_coalesce && g_ap_undo->len > 0)
		return;
	mpp_aps_t *snap = mpp_aps_snapshot(run->aps);
	if (!snap) return;
	g_ptr_array_add(g_ap_undo, snap);
	ap_undo_free_all(g_ap_redo);   /* a fresh edit invalidates the redo stack */
	g_ap_undo_coalesce = coalesce_key;
	ap_sync_undo_actions();
}

gboolean mpp_ap_editor_undo(void) {
	mpp_run_t *run = mpp_get_cached_run();
	if (!run || !run->aps || !g_ap_undo || g_ap_undo->len == 0) return FALSE;
	g_ptr_array_add(g_ap_redo, mpp_aps_snapshot(run->aps));   /* current → redo */
	mpp_aps_t *prev = (mpp_aps_t *) g_ptr_array_remove_index(g_ap_undo, g_ap_undo->len - 1);
	mpp_aps_restore(run, prev);   /* run takes ownership */
	g_ap_undo_coalesce = -1;
	g_ap_selected_idx = -1;       /* AP indices may have shifted */
	editor_refresh_count_label();
	ap_sync_undo_actions();
	redraw(REDRAW_OVERLAY);
	return TRUE;
}

gboolean mpp_ap_editor_redo(void) {
	mpp_run_t *run = mpp_get_cached_run();
	if (!run || !run->aps || !g_ap_redo || g_ap_redo->len == 0) return FALSE;
	g_ptr_array_add(g_ap_undo, mpp_aps_snapshot(run->aps));   /* current → undo */
	mpp_aps_t *next = (mpp_aps_t *) g_ptr_array_remove_index(g_ap_redo, g_ap_redo->len - 1);
	mpp_aps_restore(run, next);
	g_ap_undo_coalesce = -1;
	g_ap_selected_idx = -1;
	editor_refresh_count_label();
	ap_sync_undo_actions();
	redraw(REDRAW_OVERLAY);
	return TRUE;
}

gboolean mpp_ap_editor_is_open(void) {
	return dialog && gtk_widget_get_visible(dialog);
}

/* Open the editor dialog from the "Edit APs…" button on the MPP sub-panel.
 * Set transient_for the control window so the editor stays above the
 * main window when the user clicks on the image to edit APs. */
void on_seqmpp_edit_aps_button_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (!mpp_get_cached_run()) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No analysis result"),
		                     _("There is no multipoint analysis to edit yet. "
		                       "Run Analyze first."));
		return;
	}
	editor_init_statics();
	if (!dialog) return;
	GtkWindow *parent = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	if (parent)
		gtk_window_set_transient_for(GTK_WINDOW(dialog), parent);
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
	ap_undo_clear();   /* fresh undo history for this editing session */
	ap_sync_undo_actions();   /* take over the Undo/Redo actions (empty -> disabled) */
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
	g_ap_selected_idx = -1;
	ap_undo_clear();          /* drop the session's undo history */
	ap_sync_undo_actions();   /* hand the Undo/Redo actions back to image history */
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
		siril_log_error(_("AP editor: no cached run.\n"));
		return;
	}
	mpp_config_t cfg;
	cfg_from_spinners(run, &cfg);
	mpp_ap_editor_record_undo(-1);
	const int rc = mpp_ap_replace(run, &cfg);
	g_ap_selected_idx = -1;   /* AP set replaced — old index is meaningless */
	switch (rc) {
		case MPP_OK:
			siril_log_status(_("AP editor: auto-placed %d APs.\n"), run->aps->count);
			break;
		case MPP_ENODATA:
			siril_log_error(_("AP editor: zero APs survived the threshold filters. "
			                          "Loosen min-brightness/contrast/structure and try again.\n"));
			break;
		case MPP_ENOMEM:
			siril_log_error(_("AP editor: out of memory.\n"));
			break;
		default:
			siril_log_error(_("AP editor: auto-place failed (code %d).\n"), rc);
			break;
	}
	editor_refresh_count_label();
	redraw(REDRAW_OVERLAY);
}

void on_mpp_ap_editor_clear_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	mpp_run_t *run = mpp_get_cached_run();
	if (!run) return;
	mpp_ap_editor_record_undo(-1);
	mpp_ap_clear_all(run);
	g_ap_selected_idx = -1;
	siril_log_warning(_("AP editor: cleared all APs.\n"));
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
		siril_log_info(_("AP editor: edits committed. Register will refresh "
		                          "per-AP qualities for the edited grid automatically.\n"));
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
