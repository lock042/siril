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
 */

/*
 * flis_gui.c — FLIS layers panel.  Stage 4 §4.1 + §4.2.
 *
 * Organisation:
 *   • Row-model GObject (FlisRowItem) — wraps either a flis_layer_t* or
 *     flis_group_t* so a single GListStore can drive the GtkListView.
 *   • struct flis_panel — singleton holding every widget pointer the
 *     refresh path needs to touch.  Allocated lazily on first show.
 *   • Build helpers — one per major section (toolbar, list, property
 *     panel, mask sub-frame, context menu).  Each appends to the
 *     panel struct and registers its own signal handlers.
 *   • Refresh — rebuild the GListStore from com.uniq state and sync
 *     property widgets to the currently-selected layer.
 *   • Action handlers — every widget signal routes through a
 *     dispatch_op() helper that builds a generic_layer_args, sets
 *     invalidate_flags appropriately (per §3.4), and submits to
 *     generic_layer_worker.
 *
 * Why procedural: the layer list and property panel morph with
 * selection (no layer / single layer / group), the mask sub-frame
 * appears/disappears based on whether the layer has an lmask, and the
 * tint frame's sensitivity flips with mono/RGB.  Driving all that from
 * a static .ui file would require either layering signals onto
 * GtkBuilder templates with substantial code-side rewiring, or
 * embedding so many runtime widget swaps that the static description
 * adds little.  Construction and refresh stay co-located here.
 */

#include "flis_gui.h"

#include <math.h>
#include <float.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"          /* gui_function */
#include "core/processing.h"
#include "core/undo.h"
#include "core/nde_history.h"
#include "core/nde_replay.h"    /* nde_record_amendable/deletable, amend/delete_start */
#include "gui-gtk4/nde_editors.h"
#include "core/op_descriptor.h" /* op_descriptor_by_id for the edit round-trip check */
#include "core/gui_iface.h"
#include "io/image_format_flis.h"
#include "io/single_image.h"
#include "gui-gtk4/image_interactions.h"   /* mouse_status, MOUSE_ACTION_FLIS_DRAG_LAYER */
#include "gui-gtk4/gui_state.h"            /* gui.flis_layer_dragging */
#include "gui-gtk4/callbacks.h"            /* set_GUI_CWD (header bar refresh) */
#include "gui-gtk4/message_dialog.h"       /* siril_confirm_dialog */
#include "core/initfile.h"                 /* writeinitfile — persist "don't ask again" */
#include "gui-gtk4/utils.h"                /* siril_toggle_*, siril_drop_down_* */
#include "registration/flis_register.h"    /* flis_register_layers primitive */
#include "algos/geometry.h"                /* verbose_rotate_image for canvas dialog */

extern GtkWidget *lookup_widget(const gchar *widget_name);
extern gboolean is_current_image_flis(void);

/* =========================================================================
 * Row model: one GObject per layer or group, suitable for GListStore.
 * The model is reset on every refresh — layers / groups themselves are
 * borrowed pointers (lifetime owned by com.uniq).
 * ========================================================================= */

#define FLIS_ROW_KIND_LAYER  0
#define FLIS_ROW_KIND_GROUP  1

#define FLIS_TYPE_ROW_ITEM (flis_row_item_get_type())
G_DECLARE_FINAL_TYPE(FlisRowItem, flis_row_item, FLIS, ROW_ITEM, GObject)

struct _FlisRowItem {
	GObject parent_instance;
	int   kind;          /* FLIS_ROW_KIND_LAYER / FLIS_ROW_KIND_GROUP */
	gint  item_id;       /* layer or group item_id (look up live each access) */
	int   indent_level;  /* 0 = flat row, 1 = grouped-layer (indent for hierarchy) */
};

G_DEFINE_TYPE(FlisRowItem, flis_row_item, G_TYPE_OBJECT)

static void flis_row_item_init       (FlisRowItem *self)       { (void)self; }
static void flis_row_item_class_init (FlisRowItemClass *klass) { (void)klass; }

static FlisRowItem *flis_row_item_new(int kind, gint item_id, int indent) {
	FlisRowItem *r = g_object_new(FLIS_TYPE_ROW_ITEM, NULL);
	r->kind         = kind;
	r->item_id      = item_id;
	r->indent_level = indent;
	return r;
}

/* Resolve the live layer/group from a row item.  Safe to call after
 * mutations — returns NULL if the item was removed in the interim. */
static flis_layer_t *row_layer(FlisRowItem *row) {
	if (!row || row->kind != FLIS_ROW_KIND_LAYER) return NULL;
	return flis_layer_get_by_id(row->item_id);
}
static flis_group_t *row_group(FlisRowItem *row) {
	if (!row || row->kind != FLIS_ROW_KIND_GROUP) return NULL;
	return flis_group_get_by_id(row->item_id);
}

/* ---- NdeHistRowItem: main-thread display mirror of one nde_record ------
 * Holds copied strings only (sketch §6a/§16): the bind path never touches
 * the canonical log, and the mirror survives records whose target layer
 * has since been deleted. */
#define NDE_TYPE_HIST_ROW_ITEM (nde_hist_row_item_get_type())
G_DECLARE_FINAL_TYPE(NdeHistRowItem, nde_hist_row_item, NDE, HIST_ROW_ITEM, GObject)

struct _NdeHistRowItem {
	GObject parent_instance;
	gchar   *badge;      /* "#7" (Tier B marked "#7·B") */
	gchar   *summary;
	gchar   *target;     /* layer name / "canvas" / "document" / "deleted layer" */
	gchar   *detail;     /* popover + tooltip text */
	gboolean dead;       /* beyond live_count (undone) → dimmed */
	/* Edit/Delete UI state (P3.D) — resolved at refresh from the snapshot. */
	gint64   record_id;
	gboolean amendable;  /* nde_record_amendable() — Tier A + deserializer */
	gboolean deletable;  /* nde_record_deletable() — not structural/compositing */
	gchar   *params;     /* raw kv blob copy (NULL when the record has none) */
	gchar   *op_id;      /* descriptor id — native-editor registry key (C4) */
	/* Phase-4 barrier freeze state (nde-phase4-plan.md P4.4), from the
	 * per-item chain verdicts computed once per distinct item in
	 * refresh_history — no per-row replay. */
	gboolean frozen;     /* at-or-before the item's last barrier (chain member) */
	gboolean barrier;    /* an opaque restart-point member */
	gboolean last_barrier;/* the most recent barrier: Edit off, Delete unlocks */
	/* Reorder state (C5): the adjacent CHAIN members of the same item (not the
	 * adjacent panel rows, which interleave other layers / structural records).
	 * record_id 0 means no such neighbour.  *_in_tail says the neighbour lies
	 * in the editable tail; in_tail says this row does. */
	gboolean in_tail;
	gint64   prev_member_id, next_member_id;
	gboolean prev_in_tail,  next_in_tail;
};

G_DEFINE_TYPE(NdeHistRowItem, nde_hist_row_item, G_TYPE_OBJECT)

static void nde_hist_row_item_finalize(GObject *obj) {
	NdeHistRowItem *self = (NdeHistRowItem *)obj;
	g_free(self->badge);
	g_free(self->summary);
	g_free(self->target);
	g_free(self->detail);
	g_free(self->params);
	g_free(self->op_id);
	G_OBJECT_CLASS(nde_hist_row_item_parent_class)->finalize(obj);
}
static void nde_hist_row_item_class_init(NdeHistRowItemClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = nde_hist_row_item_finalize;
}
static void nde_hist_row_item_init(NdeHistRowItem *self) { (void)self; }

/* =========================================================================
 * Panel singleton
 * ========================================================================= */

struct flis_panel {
	GtkWidget *window;
	GtkWidget *mode_label;          /* FITS / FLIS / FLIS·active group */

	/* Toolbar */
	GtkWidget *btn_add, *btn_remove, *btn_duplicate, *btn_group;
	GtkWidget *btn_drag, *btn_move_up, *btn_move_down;
	GtkWidget *btn_canvas;   /* opens the canvas properties dialog */
	GtkWidget *btn_menu;

	/* List */
	GListStore *list_store;         /* of FlisRowItem* */
	GtkSingleSelection *list_sel;
	GtkWidget *list_view;
	GtkWidget *edge_drop_top;       /* drop zone above the list — "to top of stack" */
	GtkWidget *edge_drop_bottom;    /* drop zone below the list — "to bottom of stack" */

	/* Property panel */
	GtkWidget *prop_frame;
	GtkWidget *name_entry;
	GtkWidget *blend_dropdown;
	GtkWidget *opacity_scale;
	GtkWidget *opacity_spin;
	GtkWidget *opacity_apply;
	GtkAdjustment *opacity_adj;

	/* Tint sub-frame */
	GtkWidget *tint_frame;
	GtkWidget *tint_check;
	GtkWidget *tint_color_btn;
	GtkWidget *tint_hint_label;

	/* Mask sub-frame */
	GtkWidget *mask_frame;
	GtkWidget *mask_status_btn;
	GtkWidget *mask_toggle_btn;
	GtkWidget *mask_move_btn;
	GtkWidget *mask_view_row;
	GtkWidget *mask_view_proc_radio;
	GtkWidget *mask_view_layer_radio;

	/* History section (NDE provenance, sketch §16) — a large popover under
	 * the header-bar clock button since nde-phase5. */
	GtkWidget  *hist_pop;           /* the popover itself */
	GtkWidget  *hist_header;        /* "Nondestructive History (n)" heading */
	GtkWidget  *hist_stale_revealer;
	GtkWidget  *hist_fits_hint;     /* plain-FITS: history is in-memory only */
	GListStore *hist_store;          /* of NdeHistRowItem* */
	GtkWidget  *hist_list;
	GtkWidget  *hist_popover;        /* lazy; parented to hist_list */
	GtkWidget  *hist_popover_label;
	GtkWidget  *hist_edit_btn;       /* "Edit parameters…" in the popover */
	GtkWidget  *hist_delete_btn;     /* "Delete step" in the popover */
	GtkWidget  *hist_move_up_btn;    /* "Move up" in the popover (C5) */
	GtkWidget  *hist_move_down_btn;  /* "Move down" in the popover (C5) */
	NdeHistRowItem *hist_popover_row;/* row the popover currently targets (ref held) */

	/* Refresh suppression: TRUE while we're programmatically syncing
	 * widget state to the model, so signal handlers don't loop back. */
	gboolean   refreshing;
	/* Pending idle refresh ID (0 = none).  Coalesces multiple
	 * flis_gui_update_from_idle calls during a worker burst. */
	/* What the property panel's widgets currently target.  Set by
	 * on_selection_changed; read by name/blend/opacity handlers to
	 * dispatch to either flis_layer_set_* or flis_group_set_*. */
	int        selected_kind;       /* -1, FLIS_ROW_KIND_LAYER, FLIS_ROW_KIND_GROUP */
	gint       selected_item_id;
};

static struct flis_panel *g_panel;  /* NULL until first show */

/* Forward decls of the section builders / handlers. */
static void   build_panel(void);
static void   build_toolbar  (GtkWidget *box);
static void   build_list     (GtkWidget *box);
static void   build_property (GtkWidget *box);
static void   build_mask     (GtkWidget *box);
static void   build_history_popover(void);
static void   refresh_history(void);
static GMenu *build_context_menu(void);
static void   register_panel_actions(void);
static void   refresh_panel(void);
static void   sync_property_widgets(flis_layer_t *lay);
static void   update_toolbar_sensitivity(void);
static void   ensure_flis_css(void);
static flis_layer_t *current_selected_layer(void);
static void dispatch_ungroup(gint layer_id);
static flis_group_t *current_selected_group(void);

/* =========================================================================
 * Public entry points
 * ========================================================================= */

void flis_gui_toggle_visible(void) {
	if (!g_panel) build_panel();
	if (!g_panel || !g_panel->window) return;
	if (gtk_widget_get_visible(g_panel->window)) {
		gtk_widget_set_visible(g_panel->window, FALSE);
	} else {
		refresh_panel();
		gtk_widget_set_visible(g_panel->window, TRUE);
		gtk_window_present(GTK_WINDOW(g_panel->window));
	}
}

void flis_gui_present_if_flis(void) {
	if (!is_current_image_flis()) return;
	if (!g_panel) build_panel();
	if (!g_panel || !g_panel->window) return;
	if (gtk_widget_get_visible(g_panel->window)) {
		/* Already up — just refresh to pick up the freshly-loaded stack. */
		refresh_panel();
		return;
	}
	refresh_panel();
	gtk_widget_set_visible(g_panel->window, TRUE);
	gtk_window_present(GTK_WINDOW(g_panel->window));
}

static void canvas_dialog_refresh_external(void);

/* Startup hook (initialize_all_GUI): build the history popover and attach it
 * to the header-bar clock button (nde_history_button), which then opens and
 * closes it natively.  Builds the panel structures too; the layers window
 * itself stays hidden. */
void flis_gui_history_popover_init(void) {
	if (!g_panel) build_panel();
	if (!g_panel) return;
	build_history_popover();
}

/* Programmatic toggle (win.show-nde-history). */
void flis_gui_history_toggle_visible(void) {
	flis_gui_history_popover_init();
	if (!g_panel || !g_panel->hist_pop) return;
	if (gtk_widget_get_visible(g_panel->hist_pop))
		gtk_popover_popdown(GTK_POPOVER(g_panel->hist_pop));
	else
		gtk_popover_popup(GTK_POPOVER(g_panel->hist_pop));
}

/* Coalescing flag for refresh_idle_cb.  Atomic because
 * flis_gui_update_from_idle() is documented safe from worker threads: with
 * the old per-panel source-id scheme, the idle could run and zero the id
 * on the main loop BEFORE the worker's g_idle_add return value was stored,
 * leaving a stale nonzero id that silenced every future refresh.  The flag
 * is claimed with a CAS before queueing and released as the idle's first
 * action, so a state change landing mid-refresh queues a fresh idle. */
static gint refresh_pending = 0;

static gboolean refresh_idle_cb(gpointer p) {
	(void)p;
	g_atomic_int_set(&refresh_pending, 0);
	/* Title bar reflects FLIS active-layer state independently of the panel.
	 * Refresh it on every coalesced tick so name / count / active changes
	 * land in the header even when the layers panel is hidden. */
	gui_function(set_GUI_CWD, NULL);
	if (g_panel) {
		if (gtk_widget_get_visible(g_panel->window))
			refresh_panel();
	}
	/* Same idle path repaints the canvas-properties dialog so it
	 * reflects any FLIS state change that happened underneath it
	 * (e.g. an undo / redo of a canvas op). */
	canvas_dialog_refresh_external();
	return G_SOURCE_REMOVE;
}

void flis_gui_update_from_idle(void) {
	if (!g_atomic_int_compare_and_exchange(&refresh_pending, 0, 1))
		return;  /* already pending */
	g_idle_add(refresh_idle_cb, NULL);
}

/* =========================================================================
 * Build helpers
 * ========================================================================= */

static void build_panel(void) {
	if (g_panel) return;
	g_panel = g_new0(struct flis_panel, 1);
	ensure_flis_css();

	GtkWidget *w = gtk_window_new();
	g_panel->window = w;
	gtk_window_set_title(GTK_WINDOW(w), _("FLIS Layers"));
	/* Compact default size; the layer list inside the scrolled window
	 * absorbs vertical growth, so we don't need a tall initial window. */
	gtk_window_set_default_size(GTK_WINDOW(w), 320, 540);
	gtk_window_set_hide_on_close(GTK_WINDOW(w), TRUE);
	GtkWidget *main_w = lookup_widget("control_window");
	if (main_w)
		gtk_window_set_transient_for(GTK_WINDOW(w), GTK_WINDOW(main_w));

	GtkWidget *outer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
	gtk_widget_set_margin_start (outer, 4);
	gtk_widget_set_margin_end   (outer, 4);
	gtk_widget_set_margin_top   (outer, 4);
	gtk_widget_set_margin_bottom(outer, 4);
	gtk_window_set_child(GTK_WINDOW(w), outer);

	/* Mode indicator — folded into the toolbar to save a row.  Kept as a
	 * member for the refresh path which writes "FITS" / "FLIS" / "FLIS ·
	 * group X selected". */
	g_panel->mode_label = gtk_label_new("FITS");
	gtk_label_set_xalign(GTK_LABEL(g_panel->mode_label), 0.0f);
	gtk_widget_add_css_class(g_panel->mode_label, "dim-label");

	build_toolbar (outer);
	build_list    (outer);
	build_property(outer);
	build_mask    (outer);

	register_panel_actions();
}

/* ---- Toolbar -------------------------------------------------------- */

static GtkWidget *icon_button(const char *icon, const char *tooltip) {
	GtkWidget *b = gtk_button_new_from_icon_name(icon);
	gtk_widget_set_tooltip_text(b, tooltip);
	return b;
}

static GtkWidget *icon_toggle(const char *icon, const char *tooltip) {
	GtkWidget *b = gtk_toggle_button_new();
	gtk_button_set_icon_name(GTK_BUTTON(b), icon);
	gtk_widget_set_tooltip_text(b, tooltip);
	return b;
}

/* Forward decls of toolbar handlers — defined below after dispatch_op. */
static void on_add_clicked       (GtkButton *b, gpointer u);
static void on_remove_clicked    (GtkButton *b, gpointer u);
static void on_duplicate_clicked (GtkButton *b, gpointer u);
static void on_group_clicked     (GtkButton *b, gpointer u);
static void on_canvas_clicked    (GtkButton *b, gpointer u);
static void on_drag_toggled      (GtkToggleButton *b, gpointer u);
static void on_move_up_clicked   (GtkButton *b, gpointer u);
static void on_move_down_clicked (GtkButton *b, gpointer u);

static void build_toolbar(GtkWidget *box) {
	GtkWidget *bar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_widget_add_css_class(bar, "toolbar");

	g_panel->btn_add       = icon_button("list-add-symbolic",        _("Add layer"));
	g_panel->btn_remove    = icon_button("list-remove-symbolic",     _("Remove selected layer"));
	g_panel->btn_duplicate = icon_button("edit-copy-symbolic",       _("Duplicate selected layer"));
	g_panel->btn_group     = icon_button("folder-new-symbolic",      _("Create group"));
	g_panel->btn_drag      = icon_toggle("input-mouse-symbolic",     _("Drag layer in canvas"));
	g_panel->btn_move_up   = icon_button("go-up-symbolic",           _("Move layer up"));
	g_panel->btn_move_down = icon_button("go-down-symbolic",         _("Move layer down"));
	g_panel->btn_canvas    = icon_button("document-page-setup-symbolic",
	                                      _("Canvas properties — resize, fit, rotate, mirror…"));

	g_signal_connect(g_panel->btn_add,       "clicked", G_CALLBACK(on_add_clicked),       NULL);
	g_signal_connect(g_panel->btn_remove,    "clicked", G_CALLBACK(on_remove_clicked),    NULL);
	g_signal_connect(g_panel->btn_duplicate, "clicked", G_CALLBACK(on_duplicate_clicked), NULL);
	g_signal_connect(g_panel->btn_group,     "clicked", G_CALLBACK(on_group_clicked),     NULL);
	g_signal_connect(g_panel->btn_drag,      "toggled", G_CALLBACK(on_drag_toggled),      NULL);
	g_signal_connect(g_panel->btn_move_up,   "clicked", G_CALLBACK(on_move_up_clicked),   NULL);
	g_signal_connect(g_panel->btn_move_down, "clicked", G_CALLBACK(on_move_down_clicked), NULL);
	g_signal_connect(g_panel->btn_canvas,    "clicked", G_CALLBACK(on_canvas_clicked),    NULL);

	gtk_box_append(GTK_BOX(bar), g_panel->btn_add);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_remove);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_duplicate);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_group);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_drag);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_move_up);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_move_down);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_canvas);

	/* Mode label between the toolbar buttons and the right-side menu —
	 * keeps the panel header to a single row. */
	gtk_widget_set_margin_start(g_panel->mode_label, 8);
	gtk_widget_set_hexpand(g_panel->mode_label, TRUE);
	gtk_box_append(GTK_BOX(bar), g_panel->mode_label);

	g_panel->btn_menu = gtk_menu_button_new();
	gtk_menu_button_set_icon_name(GTK_MENU_BUTTON(g_panel->btn_menu), "open-menu-symbolic");
	gtk_widget_set_tooltip_text(g_panel->btn_menu, _("Layer actions"));
	GMenu *ctx = build_context_menu();
	gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(g_panel->btn_menu), G_MENU_MODEL(ctx));
	g_object_unref(ctx);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_menu);

	gtk_box_append(GTK_BOX(box), bar);
}

/* ---- Layer list ---------------------------------------------------- */

static void on_row_visible_toggled  (GtkToggleButton *btn, gpointer u);
static void on_row_lock_toggled     (GtkToggleButton *btn, gpointer u);
static void on_group_expander_clicked(GtkButton *btn, gpointer u);
static void on_selection_changed    (GtkSelectionModel *sel, guint pos, guint nitems, gpointer u);

/* The row template — built procedurally rather than from a .ui file
 * since it's small.  Each row is a horizontal box of: expander (groups
 * only, invisible for layers), visibility toggle, lock toggle,
 * thumbnail, and name label. */
typedef struct {
	GtkWidget *row_box;
	GtkWidget *expander;          /* chevron button on group rows; hidden on layer rows */
	GtkWidget *visible_toggle;
	GtkWidget *lock_toggle;
	GtkWidget *thumb;
	GtkWidget *name_label;
	GtkWidget *kind_badge;        /* "group" badge for group rows */
	gulong     vis_handler_id;
	gulong     lock_handler_id;
	gulong     exp_handler_id;
	/* Drag-to-reorder state — cached so bind can update item_id without
	 * tearing down the controllers each time. */
	GtkDragSource *drag_source;
	GtkDropTarget *drop_target;
	gint           current_item_id;
	int            current_kind;
} row_widgets_t;

/* Load global CSS for the panel — runs once per process.  Currently
 * styles the active-layer-row marker (subtle highlight tint) and the
 * "drop into stack edges" zones. */
static GtkCssProvider *g_flis_css_provider = NULL;
static void ensure_flis_css(void) {
	if (g_flis_css_provider) return;
	g_flis_css_provider = gtk_css_provider_new();
	gtk_css_provider_load_from_string(g_flis_css_provider,
		".flis-active-layer-row {"
		"  background-color: alpha(@accent_bg_color, 0.18);"
		"  border-left: 3px solid @accent_color;"
		"}"
		".flis-edge-drop {"
		"  min-height: 18px;"
		"  background: transparent;"
		"  border: 1px dashed alpha(@borders, 0.6);"
		"  border-radius: 4px;"
		"  margin: 2px 4px;"
		"  color: alpha(@view_fg_color, 0.55);"
		"  font-size: 90%;"
		"}"
		".flis-edge-drop.drop-active {"
		"  background-color: alpha(@accent_bg_color, 0.25);"
		"  border-style: solid;"
		"  color: @accent_fg_color;"
		"}"
	);
	gtk_style_context_add_provider_for_display(gdk_display_get_default(),
		GTK_STYLE_PROVIDER(g_flis_css_provider),
		GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
}

/* Drag-source prepare callback: builds a content provider carrying the
 * source row's item_id as a G_TYPE_INT. */
static GdkContentProvider *row_drag_prepare(GtkDragSource *src,
                                             double x, double y, gpointer ud) {
	(void)src; (void)x; (void)y;
	row_widgets_t *rw = ud;
	if (!rw || rw->current_kind != FLIS_ROW_KIND_LAYER || rw->current_item_id == 0)
		return NULL;     /* groups don't drag-reorder yet */
	GValue v = G_VALUE_INIT;
	g_value_init(&v, G_TYPE_INT);
	g_value_set_int(&v, rw->current_item_id);
	GdkContentProvider *p = gdk_content_provider_new_for_value(&v);
	g_value_unset(&v);
	return p;
}

/* Drop target callback: invoked when a layer is dropped on this row.
 * Source item_id arrives in @value (G_TYPE_INT); target is rw's. */
static gboolean row_drop_received(GtkDropTarget *tgt, const GValue *value,
                                   double x, double y, gpointer ud) {
	(void)tgt; (void)x;
	row_widgets_t *rw = ud;
	if (!rw || rw->current_kind != FLIS_ROW_KIND_LAYER || rw->current_item_id == 0)
		return FALSE;     /* can't drop on group rows */
	if (!G_VALUE_HOLDS_INT(value)) return FALSE;
	const gint src_id = g_value_get_int(value);
	if (src_id == rw->current_item_id) return FALSE;  /* dropped on self */

	/* Drop-zone heuristic: top half of the row → place ABOVE target
	 * (higher z); bottom half → place BELOW.  Matches Photoshop's
	 * intuitive feel where the dropped layer ends up where the cursor
	 * landed in the list. */
	const int h = gtk_widget_get_height(rw->row_box);
	const gboolean place_above = (h > 0) ? (y < (double)h * 0.5) : TRUE;

	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_reorder_args_free;
	payload->target_id   = rw->current_item_id;
	payload->place_above = place_above;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->description        = g_strdup(_("Reorder layer"));
	args->invalidate_flags   = FLIS_INV_STACK;
	args->invalidate_item_id = src_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
	return TRUE;
}

static void on_row_setup(GtkListItemFactory *f, GtkListItem *item, gpointer u) {
	(void)f; (void)u;
	row_widgets_t *rw = g_new0(row_widgets_t, 1);

	rw->row_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
	gtk_widget_set_margin_start (rw->row_box, 2);
	gtk_widget_set_margin_end   (rw->row_box, 2);
	gtk_widget_set_margin_top   (rw->row_box, 1);
	gtk_widget_set_margin_bottom(rw->row_box, 1);

	/* Plain GtkButton (not toggle) — we don't track state on the widget,
	 * grp->collapsed is the source of truth; on click we just flip it
	 * and refresh.  Using a toggle button raced with the bind path's
	 * set_active call and produced a "two clicks needed" behaviour. */
	rw->expander = gtk_button_new_from_icon_name("pan-down-symbolic");
	gtk_widget_set_tooltip_text(rw->expander, _("Toggle group collapse"));
	gtk_widget_add_css_class(rw->expander, "flat");

	rw->visible_toggle = gtk_toggle_button_new();
	gtk_button_set_icon_name(GTK_BUTTON(rw->visible_toggle), "view-reveal-symbolic");
	gtk_widget_set_tooltip_text(rw->visible_toggle, _("Toggle visibility"));
	gtk_widget_add_css_class(rw->visible_toggle, "flat");

	rw->lock_toggle = gtk_toggle_button_new();
	/* Initial icon = "allow edits" (unlocked).  sync_property_widgets
	 * swaps to "prevent" when the layer is locked. */
	gtk_button_set_icon_name(GTK_BUTTON(rw->lock_toggle), "changes-allow-symbolic");
	gtk_widget_set_tooltip_text(rw->lock_toggle, _("Unlocked — click to lock"));
	gtk_widget_add_css_class(rw->lock_toggle, "flat");

	rw->thumb = gtk_picture_new();
	gtk_widget_set_size_request(rw->thumb, 32, 32);

	rw->name_label = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(rw->name_label), 0.0f);
	gtk_widget_set_hexpand(rw->name_label, TRUE);

	rw->kind_badge = gtk_label_new(NULL);
	gtk_widget_add_css_class(rw->kind_badge, "dim-label");

	gtk_box_append(GTK_BOX(rw->row_box), rw->expander);
	gtk_box_append(GTK_BOX(rw->row_box), rw->visible_toggle);
	gtk_box_append(GTK_BOX(rw->row_box), rw->lock_toggle);
	gtk_box_append(GTK_BOX(rw->row_box), rw->thumb);
	gtk_box_append(GTK_BOX(rw->row_box), rw->name_label);
	gtk_box_append(GTK_BOX(rw->row_box), rw->kind_badge);

	/* Drag source + drop target for drag-to-reorder (§4.3 slice 4).
	 * Attached once per row widget; bind updates rw->current_item_id
	 * so the callbacks see the right layer.  No teardown needed —
	 * controllers are owned by the row_box widget. */
	rw->drag_source = gtk_drag_source_new();
	gtk_drag_source_set_actions(rw->drag_source, GDK_ACTION_MOVE);
	g_signal_connect(rw->drag_source, "prepare",
	                 G_CALLBACK(row_drag_prepare), rw);
	gtk_widget_add_controller(rw->row_box, GTK_EVENT_CONTROLLER(rw->drag_source));

	rw->drop_target = gtk_drop_target_new(G_TYPE_INT, GDK_ACTION_MOVE);
	g_signal_connect(rw->drop_target, "drop",
	                 G_CALLBACK(row_drop_received), rw);
	gtk_widget_add_controller(rw->row_box, GTK_EVENT_CONTROLLER(rw->drop_target));

	g_object_set_data_full(G_OBJECT(item), "row-widgets", rw, g_free);
	gtk_list_item_set_child(item, rw->row_box);
}

/* Build a 32×32 BGRA thumbnail of a layer's pixels.  Cheap nearest-
 * neighbour downsample; called once per layer per panel refresh.  The
 * texture is freshly allocated on every call — GdkPicture takes its
 * own ref so it's safe to unref locally. */
static GdkTexture *build_layer_thumb(flis_layer_t *lay) {
	if (!lay || !lay->fit) return NULL;
	const fits *f = lay->fit;
	const int sw = (int)f->rx;
	const int sh = (int)f->ry;
	if (sw <= 0 || sh <= 0) return NULL;
	const int tw = 32, th = 32;
	uint8_t *bgra = g_malloc((gsize)tw * th * 4);

	const gboolean is_rgb = (f->naxes[2] >= 3);
	for (int y = 0; y < th; y++) {
		const int sy_disp = y * sh / th;
		const int sy_fits = sh - 1 - sy_disp;
		for (int x = 0; x < tw; x++) {
			const int sx = x * sw / tw;
			const size_t si = (size_t)sy_fits * sw + (size_t)sx;
			uint8_t r, g, b;
			if (f->type == DATA_FLOAT) {
				float fr = f->fpdata[0][si];
				float fg = is_rgb ? f->fpdata[1][si] : fr;
				float fb = is_rgb ? f->fpdata[2][si] : fr;
				r = (uint8_t)CLAMP((int)(fr * 255.f), 0, 255);
				g = (uint8_t)CLAMP((int)(fg * 255.f), 0, 255);
				b = (uint8_t)CLAMP((int)(fb * 255.f), 0, 255);
			} else {
				WORD wr = f->pdata[0][si];
				WORD wg = is_rgb ? f->pdata[1][si] : wr;
				WORD wb = is_rgb ? f->pdata[2][si] : wr;
				r = (uint8_t)(wr >> 8);
				g = (uint8_t)(wg >> 8);
				b = (uint8_t)(wb >> 8);
			}
			const size_t di = ((size_t)y * tw + (size_t)x) * 4;
			bgra[di + 0] = b;
			bgra[di + 1] = g;
			bgra[di + 2] = r;
			bgra[di + 3] = 255;
		}
	}
	GBytes *bytes = g_bytes_new_take(bgra, (gsize)tw * th * 4);
	GdkTexture *tex = gdk_memory_texture_new(tw, th, GDK_MEMORY_B8G8R8A8,
	                                          bytes, (gsize)tw * 4);
	g_bytes_unref(bytes);
	return tex;
}

static void on_row_bind(GtkListItemFactory *f, GtkListItem *item, gpointer u) {
	(void)f; (void)u;
	row_widgets_t *rw = g_object_get_data(G_OBJECT(item), "row-widgets");
	FlisRowItem   *ri = gtk_list_item_get_item(item);
	if (!rw || !ri) return;

	/* Disconnect previous handlers (rebind on scroll reuses rows). */
	if (rw->vis_handler_id) {
		g_signal_handler_disconnect(rw->visible_toggle, rw->vis_handler_id);
		rw->vis_handler_id = 0;
	}
	if (rw->lock_handler_id) {
		g_signal_handler_disconnect(rw->lock_toggle, rw->lock_handler_id);
		rw->lock_handler_id = 0;
	}
	if (rw->exp_handler_id) {
		g_signal_handler_disconnect(rw->expander, rw->exp_handler_id);
		rw->exp_handler_id = 0;
	}

	/* Reset row CSS classes (rebind may receive a row that previously
	 * displayed a group header or a different indent level / active
	 * state). */
	gtk_widget_remove_css_class(rw->row_box, "flis-group-header");
	gtk_widget_remove_css_class(rw->row_box, "flis-active-layer-row");

	/* Drag-reorder identity: callbacks read these to know which layer
	 * the row currently represents. */
	rw->current_item_id = ri->item_id;
	rw->current_kind    = ri->kind;

	/* Left padding gives the visual indent for grouped layers.  Group
	 * headers themselves sit flush left at indent_level=0. */
	const int indent_px = (ri->indent_level > 0) ? (16 * ri->indent_level + 4) : 4;
	gtk_widget_set_margin_start(rw->row_box, indent_px);

	if (ri->kind == FLIS_ROW_KIND_LAYER) {
		flis_layer_t *lay = row_layer(ri);
		if (!lay) {
			gtk_label_set_text(GTK_LABEL(rw->name_label), "?");
			return;
		}
		/* Highlight the active layer's row (the one bound to gfit).
		 * This may diverge from selection if the user clicks a group
		 * header, so it's a separate visual signal from selection. */
		flis_layer_t *active = flis_active_layer();
		if (active == lay)
			gtk_widget_add_css_class(rw->row_box, "flis-active-layer-row");
		gtk_widget_set_visible(rw->expander,       FALSE);
		gtk_widget_set_visible(rw->visible_toggle, TRUE);
		gtk_widget_set_visible(rw->lock_toggle,    TRUE);
		gtk_widget_set_visible(rw->thumb,          TRUE);
		gtk_widget_set_sensitive(rw->lock_toggle, TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rw->visible_toggle), lay->visible);
		/* Open / crossed-out eye to make the state unambiguous at a
		 * glance (the toggle button's pressed/unpressed visual alone
		 * is too subtle, especially in flat-button style). */
		gtk_button_set_icon_name(GTK_BUTTON(rw->visible_toggle),
			lay->visible ? "view-reveal-symbolic" : "view-conceal-symbolic");
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rw->lock_toggle),    lay->locked);
		/* Same treatment as visibility: swap closed-padlock / open-padlock
		 * icons so the state is obvious at a glance.  changes-prevent =
		 * "edits prevented" (locked); changes-allow = "edits allowed". */
		gtk_button_set_icon_name(GTK_BUTTON(rw->lock_toggle),
			lay->locked ? "changes-prevent-symbolic" : "changes-allow-symbolic");
		gtk_widget_set_tooltip_text(rw->lock_toggle,
			lay->locked ? _("Locked — click to unlock")
			            : _("Unlocked — click to lock"));
		gtk_label_set_text(GTK_LABEL(rw->name_label),
		                   lay->layer_name ? lay->layer_name : _("(unnamed)"));

		GdkTexture *thumb = build_layer_thumb(lay);
		gtk_picture_set_paintable(GTK_PICTURE(rw->thumb), GDK_PAINTABLE(thumb));
		if (thumb) g_object_unref(thumb);

		gtk_label_set_text(GTK_LABEL(rw->kind_badge), "");

		rw->vis_handler_id  = g_signal_connect(rw->visible_toggle, "toggled",
		                          G_CALLBACK(on_row_visible_toggled),
		                          GINT_TO_POINTER(ri->item_id));
		rw->lock_handler_id = g_signal_connect(rw->lock_toggle, "toggled",
		                          G_CALLBACK(on_row_lock_toggled),
		                          GINT_TO_POINTER(ri->item_id));
	} else {
		/* Group header row: distinguished visually so the hierarchy is
		 * obvious.  Bold label with a folder symbol, no lock toggle,
		 * no thumbnail, header CSS class for any future theme styling. */
		flis_group_t *grp = row_group(ri);
		if (!grp) {
			gtk_label_set_text(GTK_LABEL(rw->name_label), "?");
			return;
		}
		gtk_widget_add_css_class(rw->row_box, "flis-group-header");
		gtk_widget_set_visible(rw->expander,    TRUE);
		gtk_widget_set_visible(rw->lock_toggle, FALSE);
		gtk_widget_set_visible(rw->thumb,       FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rw->visible_toggle), grp->visible);
		gtk_button_set_icon_name(GTK_BUTTON(rw->visible_toggle),
			grp->visible ? "view-reveal-symbolic" : "view-conceal-symbolic");
		/* Chevron orientation: pan-down when expanded (the contents
		 * are below), pan-end when collapsed (click to expand). */
		gtk_button_set_icon_name(GTK_BUTTON(rw->expander),
			grp->collapsed ? "pan-end-symbolic" : "pan-down-symbolic");
		gchar *markup = g_markup_printf_escaped(
			"\xF0\x9F\x97\x80 <b>%s</b>",   /* U+1F5C0 file folder */
			grp->name ? grp->name : _("(group)"));
		gtk_label_set_markup(GTK_LABEL(rw->name_label), markup);
		g_free(markup);
		gtk_label_set_text(GTK_LABEL(rw->kind_badge), "");

		rw->vis_handler_id = g_signal_connect(rw->visible_toggle, "toggled",
		                         G_CALLBACK(on_row_visible_toggled),
		                         GINT_TO_POINTER(-ri->item_id));  /* negative encodes group */
		rw->exp_handler_id = g_signal_connect(rw->expander, "clicked",
		                         G_CALLBACK(on_group_expander_clicked),
		                         GINT_TO_POINTER(ri->item_id));
	}
}

/* Edge drop zones (top / bottom of list).  Always visible so users
 * can drag a layer out of a group with no external siblings (the row-
 * level drop targets all sit on layer rows, so there's no way to
 * express "drop outside any layer" without these).  Drops here also
 * force-clear the source's group_id, regardless of what the nearest
 * layer's group is. */
static gboolean edge_drop_received(GtkDropTarget *tgt, const GValue *value,
                                    double x, double y, gpointer ud) {
	(void)tgt; (void)x; (void)y;
	if (!G_VALUE_HOLDS_INT(value) || !com.uniq || !com.uniq->layers)
		return FALSE;
	const gint src_id = g_value_get_int(value);
	const gboolean to_top = (ud != NULL);  /* user_data: NULL = bottom, !NULL = top */

	/* Find the topmost / bottommost layer.  com.uniq->layers is sorted
	 * ascending by layer_order, so head = bottom, tail = top. */
	GSList *first = com.uniq->layers;
	GSList *last  = g_slist_last(com.uniq->layers);
	flis_layer_t *anchor = (flis_layer_t *)(to_top ? last->data : first->data);
	if (!anchor) return FALSE;
	if (anchor->item_id == src_id) {
		/* Source is already at this edge — anchor would be itself.
		 * We still need to handle the "drag out of group" case: even
		 * if the layer is at the right end of the stack, its group_id
		 * might be non-zero and the user dragged it here precisely to
		 * unset that.  Use the anchor as a self-target with
		 * force_ungroup; the hook handles src == tgt as a no-op so we
		 * have to special-case the ungroup here. */
		flis_layer_t *src = flis_layer_get_by_id(src_id);
		if (src && src->group_id != 0) {
			/* Dispatch through the worker: the lock check in
			 * flis_layer_set_group, the stack writer lock, and the
			 * invalidation all come from the shared path (the old
			 * inline group_id write bypassed all three). */
			dispatch_ungroup(src_id);
		}
		return TRUE;
	}

	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn    = flis_reorder_args_free;
	payload->target_id     = anchor->item_id;
	payload->place_above   = to_top;
	payload->force_ungroup = TRUE;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->description        = g_strdup(to_top ? _("Move layer to top")
	                                            : _("Move layer to bottom"));
	args->invalidate_flags   = FLIS_INV_STACK;
	args->invalidate_item_id = src_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
	return TRUE;
}

/* GtkDropTarget signals "enter"/"leave" let us add a hover CSS class so
 * the dashed zone visibly fills with the accent colour when a drag
 * is overhead.  Plain UX nicety, no functional effect. */
static GdkDragAction edge_drop_motion(GtkDropTarget *tgt, double x, double y, gpointer ud) {
	(void)x; (void)y; (void)ud;
	GtkWidget *w = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(tgt));
	if (w) gtk_widget_add_css_class(w, "drop-active");
	return GDK_ACTION_MOVE;
}
static void edge_drop_leave(GtkDropTarget *tgt, gpointer ud) {
	(void)ud;
	GtkWidget *w = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(tgt));
	if (w) gtk_widget_remove_css_class(w, "drop-active");
}

static GtkWidget *make_edge_drop_zone(const gchar *label_text, gboolean to_top) {
	GtkWidget *lbl = gtk_label_new(label_text);
	gtk_widget_add_css_class(lbl, "flis-edge-drop");
	gtk_widget_set_halign(lbl, GTK_ALIGN_FILL);
	GtkDropTarget *dt = gtk_drop_target_new(G_TYPE_INT, GDK_ACTION_MOVE);
	g_signal_connect(dt, "drop",  G_CALLBACK(edge_drop_received), to_top ? GINT_TO_POINTER(1) : NULL);
	g_signal_connect(dt, "motion", G_CALLBACK(edge_drop_motion), NULL);
	g_signal_connect(dt, "leave",  G_CALLBACK(edge_drop_leave),  NULL);
	gtk_widget_add_controller(lbl, GTK_EVENT_CONTROLLER(dt));
	return lbl;
}

static void build_list(GtkWidget *box) {
	g_panel->list_store = g_list_store_new(FLIS_TYPE_ROW_ITEM);
	g_panel->list_sel   = gtk_single_selection_new(G_LIST_MODEL(g_panel->list_store));
	gtk_single_selection_set_autoselect(g_panel->list_sel, FALSE);
	gtk_single_selection_set_can_unselect(g_panel->list_sel, TRUE);

	GtkListItemFactory *factory = gtk_signal_list_item_factory_new();
	g_signal_connect(factory, "setup", G_CALLBACK(on_row_setup), NULL);
	g_signal_connect(factory, "bind",  G_CALLBACK(on_row_bind),  NULL);

	g_panel->list_view = gtk_list_view_new(GTK_SELECTION_MODEL(g_panel->list_sel), factory);
	gtk_list_view_set_show_separators(GTK_LIST_VIEW(g_panel->list_view), TRUE);
	gtk_widget_set_vexpand(g_panel->list_view, TRUE);

	GtkWidget *sw = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
	                                GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(sw), g_panel->list_view);
	gtk_widget_set_vexpand(sw, TRUE);

	/* Edge drop zones flanking the list — sit always-visible so the
	 * user can drag a layer to the top / bottom of the stack and
	 * (most importantly) drag a layer OUT of a group even when no
	 * other ungrouped layer exists to drop on. */
	g_panel->edge_drop_top    = make_edge_drop_zone(
		_("⤴ drop here to move to top of stack (outside any group)"),  TRUE);
	g_panel->edge_drop_bottom = make_edge_drop_zone(
		_("⤵ drop here to move to bottom of stack (outside any group)"), FALSE);
	gtk_box_append(GTK_BOX(box), g_panel->edge_drop_top);
	gtk_box_append(GTK_BOX(box), sw);
	gtk_box_append(GTK_BOX(box), g_panel->edge_drop_bottom);

	g_signal_connect(g_panel->list_sel, "selection-changed",
	                 G_CALLBACK(on_selection_changed), NULL);
}

/* ---- Property panel ------------------------------------------------ */

static void on_name_activate     (GtkEntry *e,        gpointer u);
static void on_name_focus_leave  (GtkEventControllerFocus *c, gpointer u);
static void on_drag_toggled      (GtkToggleButton *b,  gpointer u);
static void on_blend_changed     (GtkDropDown *dd, GParamSpec *p, gpointer u);
static void on_opacity_changed   (GtkAdjustment *adj, gpointer u);
static void on_opacity_apply_clicked(GtkButton *btn, gpointer u);
static void opacity_edit_revert  (void);
static void on_tint_check_toggled(GtkCheckButton *btn, gpointer u);
static void on_tint_color_chosen (GtkColorDialogButton *btn, GParamSpec *p, gpointer u);

/* 18 GSK-equivalent FLIS blend modes — the dropdown's index order maps
 * to the FLIS_BLEND_* values via blend_mode_for_index(). */
static const char *blend_names[] = {
	N_("Normal"),       N_("Multiply"),   N_("Screen"),      N_("Overlay"),
	N_("Soft Light"),   N_("Hard Light"), N_("Color Dodge"), N_("Color Burn"),
	N_("Darken"),       N_("Lighten"),    N_("Difference"),  N_("Exclusion"),
	N_("Hue"),          N_("Saturation"), N_("Color"),       N_("Luminosity"),
	N_("Chroma (LRGB)"), N_("Pass-through"),
};
static const flis_blend_mode_t blend_mode_for_index_table[] = {
	FLIS_BLEND_NORMAL,      FLIS_BLEND_MULTIPLY,    FLIS_BLEND_SCREEN,      FLIS_BLEND_OVERLAY,
	FLIS_BLEND_SOFT_LIGHT,  FLIS_BLEND_HARD_LIGHT,  FLIS_BLEND_COLOR_DODGE, FLIS_BLEND_COLOR_BURN,
	FLIS_BLEND_DARKEN,      FLIS_BLEND_LIGHTEN,     FLIS_BLEND_DIFFERENCE,  FLIS_BLEND_EXCLUSION,
	FLIS_BLEND_HUE,         FLIS_BLEND_SATURATION,  FLIS_BLEND_COLOR,       FLIS_BLEND_LUMINOSITY,
	FLIS_BLEND_CHROMA,      FLIS_BLEND_PASS_THROUGH,
};
#define N_BLENDS (int)(sizeof(blend_names) / sizeof(blend_names[0]))

static int index_for_blend_mode(flis_blend_mode_t m) {
	for (int i = 0; i < N_BLENDS; i++)
		if (blend_mode_for_index_table[i] == m) return i;
	return 0;
}

static void build_property(GtkWidget *box) {
	g_panel->prop_frame = gtk_frame_new(_("Layer properties"));

	GtkWidget *grid = gtk_grid_new();
	gtk_grid_set_row_spacing   (GTK_GRID(grid), 3);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 6);
	gtk_widget_set_margin_start (grid, 6);
	gtk_widget_set_margin_end   (grid, 6);
	gtk_widget_set_margin_top   (grid, 4);
	gtk_widget_set_margin_bottom(grid, 4);

	int row = 0;

	/* Name */
	gtk_grid_attach(GTK_GRID(grid), gtk_label_new(_("Name")), 0, row, 1, 1);
	g_panel->name_entry = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(g_panel->name_entry), 32);
	{
		/* Commit the edit when focus leaves the entry, not only on
		 * Enter — clicking elsewhere silently discarded the new name. */
		GtkEventController *fc = gtk_event_controller_focus_new();
		g_signal_connect(fc, "leave", G_CALLBACK(on_name_focus_leave), NULL);
		gtk_widget_add_controller(g_panel->name_entry, fc);
	}
	gtk_widget_set_hexpand(g_panel->name_entry, TRUE);
	gtk_grid_attach(GTK_GRID(grid), g_panel->name_entry, 1, row, 2, 1);
	g_signal_connect(g_panel->name_entry, "activate",
	                 G_CALLBACK(on_name_activate), NULL);
	row++;

	/* Blend mode */
	gtk_grid_attach(GTK_GRID(grid), gtk_label_new(_("Blend")), 0, row, 1, 1);
	GtkStringList *bl = gtk_string_list_new(NULL);
	for (int i = 0; i < N_BLENDS; i++)
		gtk_string_list_append(bl, _(blend_names[i]));
	g_panel->blend_dropdown = gtk_drop_down_new(G_LIST_MODEL(bl), NULL);
	gtk_widget_set_hexpand(g_panel->blend_dropdown, TRUE);
	gtk_grid_attach(GTK_GRID(grid), g_panel->blend_dropdown, 1, row, 2, 1);
	g_signal_connect(g_panel->blend_dropdown, "notify::selected",
	                 G_CALLBACK(on_blend_changed), NULL);
	row++;

	/* Opacity */
	gtk_grid_attach(GTK_GRID(grid), gtk_label_new(_("Opacity")), 0, row, 1, 1);
	g_panel->opacity_adj = gtk_adjustment_new(100, 0, 100, 1, 10, 0);
	g_panel->opacity_scale = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, g_panel->opacity_adj);
	gtk_widget_set_hexpand(g_panel->opacity_scale, TRUE);
	gtk_scale_set_draw_value(GTK_SCALE(g_panel->opacity_scale), FALSE);
	g_panel->opacity_spin  = gtk_spin_button_new(g_panel->opacity_adj, 1.0, 0);
	g_panel->opacity_apply = gtk_button_new_from_icon_name("object-select-symbolic");
	gtk_widget_set_tooltip_text(g_panel->opacity_apply,
		_("Apply the opacity change (adds an undo step). "
		  "Selecting another layer without applying reverts it."));
	gtk_widget_set_sensitive(g_panel->opacity_apply, FALSE);
	gtk_grid_attach(GTK_GRID(grid), g_panel->opacity_scale, 1, row, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), g_panel->opacity_spin,  2, row, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), g_panel->opacity_apply, 3, row, 1, 1);
	g_signal_connect(g_panel->opacity_adj, "value-changed",
	                 G_CALLBACK(on_opacity_changed), NULL);
	g_signal_connect(g_panel->opacity_apply, "clicked",
	                 G_CALLBACK(on_opacity_apply_clicked), NULL);
	row++;

	/* Tint frame */
	g_panel->tint_frame = gtk_frame_new(_("Tint"));
	GtkWidget *tbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_margin_start (tbox, 4);
	gtk_widget_set_margin_end   (tbox, 4);
	gtk_widget_set_margin_top   (tbox, 3);
	gtk_widget_set_margin_bottom(tbox, 3);
	gtk_frame_set_child(GTK_FRAME(g_panel->tint_frame), tbox);
	g_panel->tint_check    = gtk_check_button_new_with_label(_("Apply"));
	GtkColorDialog *cd     = gtk_color_dialog_new();
	g_panel->tint_color_btn = gtk_color_dialog_button_new(cd);
	g_panel->tint_hint_label = gtk_label_new(_("e.g. red for Hα"));
	gtk_widget_add_css_class(g_panel->tint_hint_label, "dim-label");
	gtk_box_append(GTK_BOX(tbox), g_panel->tint_check);
	gtk_box_append(GTK_BOX(tbox), g_panel->tint_color_btn);
	gtk_box_append(GTK_BOX(tbox), g_panel->tint_hint_label);
	g_signal_connect(g_panel->tint_check,     "toggled",        G_CALLBACK(on_tint_check_toggled), NULL);
	g_signal_connect(g_panel->tint_color_btn, "notify::rgba",   G_CALLBACK(on_tint_color_chosen),  NULL);
	gtk_grid_attach(GTK_GRID(grid), g_panel->tint_frame, 0, row, 3, 1);
	row++;

	gtk_frame_set_child(GTK_FRAME(g_panel->prop_frame), grid);
	gtk_box_append(GTK_BOX(box), g_panel->prop_frame);
}

/* ---- Mask sub-frame ------------------------------------------------ */

static void on_mask_status_clicked  (GtkButton *b, gpointer u);
static void on_mask_toggle_clicked  (GtkButton *b, gpointer u);
static void on_mask_move_clicked    (GtkButton *b, gpointer u);
static void on_mask_view_radio_toggled(GtkCheckButton *btn, gpointer u);

static void build_mask(GtkWidget *box) {
	g_panel->mask_frame = gtk_frame_new(_("Layer mask"));
	GtkWidget *mbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
	gtk_widget_set_margin_start (mbox, 6);
	gtk_widget_set_margin_end   (mbox, 6);
	gtk_widget_set_margin_top   (mbox, 4);
	gtk_widget_set_margin_bottom(mbox, 4);
	gtk_frame_set_child(GTK_FRAME(g_panel->mask_frame), mbox);

	GtkWidget *row1 = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	g_panel->mask_status_btn = gtk_button_new_with_label(_("(no mask)"));
	gtk_widget_set_hexpand(g_panel->mask_status_btn, TRUE);
	gtk_widget_set_tooltip_text(g_panel->mask_status_btn,
		_("Click to toggle mask active; right-click for more options"));
	g_panel->mask_toggle_btn = gtk_button_new_with_label(_("Add…"));
	g_panel->mask_move_btn   = gtk_button_new_with_label(_("Move…"));
	gtk_box_append(GTK_BOX(row1), g_panel->mask_status_btn);
	gtk_box_append(GTK_BOX(row1), g_panel->mask_toggle_btn);
	gtk_box_append(GTK_BOX(row1), g_panel->mask_move_btn);
	gtk_box_append(GTK_BOX(mbox), row1);

	g_panel->mask_view_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	GtkWidget *vl = gtk_label_new(_("Show mask:"));
	g_panel->mask_view_proc_radio  = gtk_check_button_new_with_label(_("Processing"));
	g_panel->mask_view_layer_radio = gtk_check_button_new_with_label(_("Layer"));
	gtk_check_button_set_group(GTK_CHECK_BUTTON(g_panel->mask_view_layer_radio),
	                            GTK_CHECK_BUTTON(g_panel->mask_view_proc_radio));
	gtk_box_append(GTK_BOX(g_panel->mask_view_row), vl);
	gtk_box_append(GTK_BOX(g_panel->mask_view_row), g_panel->mask_view_proc_radio);
	gtk_box_append(GTK_BOX(g_panel->mask_view_row), g_panel->mask_view_layer_radio);
	gtk_box_append(GTK_BOX(mbox), g_panel->mask_view_row);

	g_signal_connect(g_panel->mask_status_btn, "clicked", G_CALLBACK(on_mask_status_clicked), NULL);
	g_signal_connect(g_panel->mask_toggle_btn, "clicked", G_CALLBACK(on_mask_toggle_clicked), NULL);
	g_signal_connect(g_panel->mask_move_btn,   "clicked", G_CALLBACK(on_mask_move_clicked),   NULL);
	g_signal_connect(g_panel->mask_view_proc_radio,  "toggled",
	                 G_CALLBACK(on_mask_view_radio_toggled), GINT_TO_POINTER(0));
	g_signal_connect(g_panel->mask_view_layer_radio, "toggled",
	                 G_CALLBACK(on_mask_view_radio_toggled), GINT_TO_POINTER(1));

	gtk_box_append(GTK_BOX(box), g_panel->mask_frame);
}

/* ---- History section (NDE provenance, sketch §16) -------------------
 * Read-only mirror of the document's edit log.  All content is copied
 * strings prepared in refresh_history() from a snapshot — GTK callbacks
 * never read the canonical log. */

/* Popover/tooltip text for one record: identity, provenance metadata and
 * the parsed params, keys sorted for a stable display. */
static gchar *build_hist_detail(const nde_record *rec) {
	GString *s = g_string_new(NULL);
	const gchar *tier_label = rec->tier == NDE_TIER_A ? _("replayable") :
	                          rec->tier == NDE_TIER_C ? _("script, replayable") :
	                          _("opaque");
	g_string_append_printf(s, "record %" G_GINT64_FORMAT ": %s v%d (%s)\n",
	                       rec->record_id, rec->op_id ? rec->op_id : "?",
	                       rec->op_version, tier_label);
	if (rec->timestamp)
		g_string_append_printf(s, "%s\n", rec->timestamp);
	if (rec->impl)
		g_string_append_printf(s, "%s\n", rec->impl);
	if (rec->mask_active)
		g_string_append_printf(s, "%s\n", _("a mask was active"));
	if (rec->params && *rec->params) {
		g_string_append_c(s, '\n');
		GHashTable *kv = nde_kv_parse(rec->params);
		GList *keys = g_list_sort(g_hash_table_get_keys(kv), (GCompareFunc)g_strcmp0);
		for (GList *l = keys; l; l = l->next)
			g_string_append_printf(s, "%s: %s\n", (char *)l->data,
			                       (char *)g_hash_table_lookup(kv, l->data));
		g_list_free(keys);
		g_hash_table_unref(kv);
	}
	if (s->len && s->str[s->len - 1] == '\n')
		g_string_truncate(s, s->len - 1);
	return g_string_free(s, FALSE);
}

static gchar *hist_target_label(const nde_record *rec) {
	if (rec->target_item_id >= 0) {
		flis_layer_t *lay = flis_layer_get_by_id(rec->target_item_id);
		if (lay)
			return g_strdup(lay->layer_name ? lay->layer_name : "");
		return g_strdup(_("deleted layer"));
	}
	if (rec->scope == NDE_SCOPE_CANVAS)
		return g_strdup(_("canvas"));
	if (rec->scope == NDE_SCOPE_DOCUMENT)
		return g_strdup(_("document"));
	return g_strdup("");
}

static void on_hist_row_setup(GtkListItemFactory *f, GtkListItem *item, gpointer u) {
	(void)f; (void)u;
	GtkWidget *row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_margin_start(row, 2);
	gtk_widget_set_margin_end  (row, 2);

	GtkWidget *badge = gtk_label_new(NULL);
	gtk_widget_add_css_class(badge, "dim-label");

	/* Lock indicator for frozen/barrier rows — matches the layer-lock icon
	 * used elsewhere in the panel.  Hidden on ordinary rows. */
	GtkWidget *lock = gtk_image_new_from_icon_name("changes-prevent-symbolic");
	gtk_widget_add_css_class(lock, "dim-label");
	gtk_widget_set_visible(lock, FALSE);

	GtkWidget *summary = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(summary), 0.0f);
	gtk_label_set_ellipsize(GTK_LABEL(summary), PANGO_ELLIPSIZE_END);
	gtk_widget_set_hexpand(summary, TRUE);

	GtkWidget *target = gtk_label_new(NULL);
	gtk_widget_add_css_class(target, "dim-label");

	gtk_box_append(GTK_BOX(row), badge);
	gtk_box_append(GTK_BOX(row), lock);
	gtk_box_append(GTK_BOX(row), summary);
	gtk_box_append(GTK_BOX(row), target);
	g_object_set_data(G_OBJECT(item), "hist-badge",   badge);
	g_object_set_data(G_OBJECT(item), "hist-lock",    lock);
	g_object_set_data(G_OBJECT(item), "hist-summary", summary);
	g_object_set_data(G_OBJECT(item), "hist-target",  target);
	gtk_list_item_set_child(item, row);
}

static void on_hist_row_bind(GtkListItemFactory *f, GtkListItem *item, gpointer u) {
	(void)f; (void)u;
	NdeHistRowItem *r = gtk_list_item_get_item(item);
	GtkWidget *row     = gtk_list_item_get_child(item);
	GtkWidget *badge   = g_object_get_data(G_OBJECT(item), "hist-badge");
	GtkWidget *lock    = g_object_get_data(G_OBJECT(item), "hist-lock");
	GtkWidget *summary = g_object_get_data(G_OBJECT(item), "hist-summary");
	GtkWidget *target  = g_object_get_data(G_OBJECT(item), "hist-target");
	if (!r || !row) return;
	gtk_label_set_text(GTK_LABEL(badge),   r->badge   ? r->badge   : "");
	gtk_label_set_text(GTK_LABEL(summary), r->summary ? r->summary : "");
	gtk_label_set_text(GTK_LABEL(target),  r->target  ? r->target  : "");
	gtk_widget_set_tooltip_text(row, r->detail);
	/* Barrier and frozen rows carry the lock icon (subtle, next to the
	 * badge); a dead row keeps its dimming — dead wins over frozen. */
	gtk_widget_set_visible(lock, r->barrier || r->frozen);
	/* Undone records (§16) and frozen prefixes both read as dimmed. */
	if (r->dead || r->frozen)
		gtk_widget_add_css_class(row, "dim-label");
	else
		gtk_widget_remove_css_class(row, "dim-label");
}

/* Consequence line shown before any amend/delete runs: an edit recomputes
 * the image from the baseline and clears the undo stack (no meta-undo,
 * sketch §7 decision 1).  Shared by the delete confirmation and the edit
 * dialog so the wording stays identical. */
#define HIST_EDIT_CONSEQUENCE \
	_("This recomputes the image from the edit history and clears the undo history.")

/* ---- parameter editor (amend) --------------------------------------------
 * A modal window with one GtkEntry per key (keys sorted), an inline error
 * label, and Cancel / Apply.  Apply folds the confirmation in: the
 * consequence line is shown permanently in the dialog, so a click on Apply
 * validates, then amends directly — no second confirm popup.  The window is
 * destroyed on close; no state outlives it. */
struct hist_edit_ctx {
	GtkWidget *window;
	GtkWidget *error_label;
	gint64     record_id;
	GPtrArray *keys;      /* gchar*, sorted — parallel to entries */
	GPtrArray *entries;   /* GtkWidget* (GtkEntry), same order */
};

static void hist_edit_ctx_free(gpointer p) {
	struct hist_edit_ctx *ctx = p;
	g_ptr_array_unref(ctx->keys);
	g_ptr_array_unref(ctx->entries);
	g_free(ctx);
}

static void on_hist_edit_cancel(GtkButton *b, gpointer u) {
	(void)b;
	struct hist_edit_ctx *ctx = u;
	gtk_window_destroy(GTK_WINDOW(ctx->window));
}

/* Round-trip the assembled blob through the record's op deserializer as a
 * client-side check (mirrors destroy_user in nde_replay.c). */
static gboolean hist_params_valid(gint64 record_id, const gchar *blob) {
	const nde_record *rec = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (r->record_id == record_id) { rec = r; break; }
	}
	gboolean ok = FALSE;
	if (rec) {
		const op_descriptor *op = op_descriptor_by_id(rec->op_id);
		if (op && op->deserialize) {
			gpointer user = op->deserialize(blob, rec->op_version);
			if (user) {
				/* destructor-first convention (see nde_replay.c) */
				void (*d)(void *) = *(void (**)(void *))user;
				if (d) d(user); else free(user);
				ok = TRUE;
			}
		}
	}
	if (snap)
		g_ptr_array_unref(snap);
	return ok;
}

static void on_hist_edit_apply(GtkButton *b, gpointer u) {
	(void)b;
	struct hist_edit_ctx *ctx = u;
	/* Rebuild the blob over the entries in the SAME sorted key order; the
	 * builder re-escapes the values. */
	GString *kv = nde_kv_start();
	for (guint i = 0; i < ctx->keys->len; i++) {
		const char *key = g_ptr_array_index(ctx->keys, i);
		GtkWidget  *entry = g_ptr_array_index(ctx->entries, i);
		const char *val = gtk_editable_get_text(GTK_EDITABLE(entry));
		nde_kv_add_str(kv, key, val ? val : "");
	}
	gchar *blob = nde_kv_end(kv);

	if (!hist_params_valid(ctx->record_id, blob)) {
		gtk_label_set_text(GTK_LABEL(ctx->error_label),
		                   _("These parameters are not valid for this operation"));
		gtk_widget_set_visible(ctx->error_label, TRUE);
		g_free(blob);
		return;   /* keep the dialog open */
	}

	gint64 record_id = ctx->record_id;
	gtk_window_destroy(GTK_WINDOW(ctx->window));   /* frees ctx via close */
	nde_amend_start(record_id, blob);              /* refresh via nde_history_changed */
	g_free(blob);
}

static void open_hist_edit_dialog(NdeHistRowItem *r) {
	GtkWidget *window = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(window),
	                     r->summary ? r->summary : _("Edit parameters"));
	gtk_window_set_modal(GTK_WINDOW(window), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(window),
	                             GTK_WINDOW(lookup_widget("control_window")));
	gtk_window_set_default_size(GTK_WINDOW(window), 360, -1);

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_widget_set_margin_start (box, 12);
	gtk_widget_set_margin_end   (box, 12);
	gtk_widget_set_margin_top   (box, 12);
	gtk_widget_set_margin_bottom(box, 12);
	gtk_window_set_child(GTK_WINDOW(window), box);

	struct hist_edit_ctx *ctx = g_new0(struct hist_edit_ctx, 1);
	ctx->window    = window;
	ctx->record_id = r->record_id;
	ctx->keys      = g_ptr_array_new_with_free_func(g_free);
	ctx->entries   = g_ptr_array_new();

	GtkWidget *grid = gtk_grid_new();
	gtk_grid_set_row_spacing   (GTK_GRID(grid), 6);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 8);

	GHashTable *kv = nde_kv_parse(r->params);   /* values already unescaped */
	GList *keys = g_list_sort(g_hash_table_get_keys(kv), (GCompareFunc)g_strcmp0);
	int row = 0;
	for (GList *l = keys; l; l = l->next, row++) {
		const char *key = l->data;
		const char *val = g_hash_table_lookup(kv, key);
		GtkWidget *label = gtk_label_new(key);
		gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
		GtkWidget *entry = gtk_entry_new();
		gtk_editable_set_text(GTK_EDITABLE(entry), val ? val : "");
		gtk_widget_set_hexpand(entry, TRUE);
		gtk_grid_attach(GTK_GRID(grid), label, 0, row, 1, 1);
		gtk_grid_attach(GTK_GRID(grid), entry, 1, row, 1, 1);
		g_ptr_array_add(ctx->keys, g_strdup(key));
		g_ptr_array_add(ctx->entries, entry);
	}
	g_list_free(keys);
	g_hash_table_unref(kv);
	gtk_box_append(GTK_BOX(box), grid);

	ctx->error_label = gtk_label_new(NULL);
	gtk_label_set_wrap(GTK_LABEL(ctx->error_label), TRUE);
	gtk_label_set_xalign(GTK_LABEL(ctx->error_label), 0.0f);
	gtk_widget_add_css_class(ctx->error_label, "error");
	gtk_widget_set_visible(ctx->error_label, FALSE);
	gtk_box_append(GTK_BOX(box), ctx->error_label);

	GtkWidget *consequence = gtk_label_new(HIST_EDIT_CONSEQUENCE);
	gtk_label_set_wrap(GTK_LABEL(consequence), TRUE);
	gtk_label_set_xalign(GTK_LABEL(consequence), 0.0f);
	gtk_widget_add_css_class(consequence, "dim-label");
	gtk_box_append(GTK_BOX(box), consequence);

	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *cancel = gtk_button_new_with_label(_("Cancel"));
	GtkWidget *apply  = gtk_button_new_with_label(_("Apply"));
	gtk_widget_add_css_class(apply, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), cancel);
	gtk_box_append(GTK_BOX(bbox), apply);
	gtk_box_append(GTK_BOX(box), bbox);

	/* ctx lifetime follows the window. */
	g_object_set_data_full(G_OBJECT(window), "hist-edit-ctx", ctx,
	                       hist_edit_ctx_free);
	g_signal_connect(cancel, "clicked", G_CALLBACK(on_hist_edit_cancel), ctx);
	g_signal_connect(apply,  "clicked", G_CALLBACK(on_hist_edit_apply),  ctx);
	gtk_window_present(GTK_WINDOW(window));
}

static void on_hist_edit_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	if (!g_panel || !g_panel->hist_popover_row)
		return;
	NdeHistRowItem *r = g_panel->hist_popover_row;
	if (!r->amendable || !r->params)
		return;
	gtk_popover_popdown(GTK_POPOVER(g_panel->hist_popover));
	/* Native editor first (amend mode with live preview); the kv grid
	 * remains the fallback for ops without one. */
	if (nde_editor_open(r->op_id, r->record_id))
		return;
	open_hist_edit_dialog(r);
}

static const char *hist_move_disabled_reason(const NdeHistRowItem *r,
                                             gboolean up);

static void on_hist_delete_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	if (!g_panel || !g_panel->hist_popover_row)
		return;
	NdeHistRowItem *r = g_panel->hist_popover_row;
	if (!r->deletable)
		return;
	gtk_popover_popdown(GTK_POPOVER(g_panel->hist_popover));

	gchar *msg = g_strdup_printf("%s\n\n%s",
	                             r->summary ? r->summary : "",
	                             HIST_EDIT_CONSEQUENCE);
	gboolean confirmed = siril_confirm_dialog(_("Delete this step?"), msg,
	                                          _("_Delete"));
	g_free(msg);
	if (confirmed)
		nde_delete_start(r->record_id);   /* refresh via nde_history_changed */
}

/* Move the popover's row before its previous chain member (@up) or after its
 * next chain member.  The anchor is the ADJACENT CHAIN MEMBER of the same
 * item, not the adjacent panel row.  Confirms first (a move rewrites the
 * lineage, like Delete). */
static void hist_move_clicked(gboolean up) {
	if (!g_panel || !g_panel->hist_popover_row)
		return;
	NdeHistRowItem *r = g_panel->hist_popover_row;
	if (hist_move_disabled_reason(r, up) != NULL)
		return;
	gint64 anchor = up ? r->prev_member_id : r->next_member_id;
	if (!anchor)
		return;
	gint64 record_id = r->record_id;
	gtk_popover_popdown(GTK_POPOVER(g_panel->hist_popover));

	gboolean confirmed = TRUE;
	if (!com.pref.gui.silent_hist_move) {
		gchar *msg = g_strdup_printf("%s\n\n%s",
		                             r->summary ? r->summary : "",
		                             HIST_EDIT_CONSEQUENCE);
		confirmed = siril_confirm_dialog_and_remember(
				up ? _("Move this step up?") : _("Move this step down?"),
				msg, _("_Move"), &com.pref.gui.silent_hist_move);
		g_free(msg);
		if (com.pref.gui.silent_hist_move)
			writeinitfile();   /* persist "don't ask again" right away */
	}
	if (confirmed)
		/* Move up = before the previous member; Move down = after the next
		 * member.  Refresh arrives via nde_history_changed. */
		nde_reorder_start(record_id, anchor, /*after=*/ up ? FALSE : TRUE);
}

static void on_hist_move_up_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	hist_move_clicked(TRUE);
}

static void on_hist_move_down_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	hist_move_clicked(FALSE);
}

/* First applicable reason a row's Edit/Delete button is greyed out, in plain
 * language; NULL when the button is enabled.  @is_edit selects the extra
 * amendable/deletable requirement. */
static const char *hist_action_disabled_reason(const NdeHistRowItem *r,
                                                gboolean is_edit) {
	if (r->dead)
		return _("This step has been undone");
	if (nde_history_is_stale())
		return _("The image was modified outside the recorded history");
	if (processing_is_job_active())
		return _("Not available while processing");
	if (nde_amend_preview_active())
		return _("Another history step is being edited");
	/* A frozen prefix step (before the last opaque step) is fully locked;
	 * the last opaque step itself can still be deleted to unlock the rest. */
	if (r->frozen && !r->last_barrier)
		return _("Locked by a later opaque step");
	if (r->last_barrier && is_edit)
		return _("Opaque steps cannot be edited");
	if (is_edit && (!r->amendable || !r->params))
		return _("Opaque steps cannot be edited");
	if (!is_edit && !r->deletable)
		return _("Structural steps cannot be deleted");
	return NULL;
}

/* First applicable reason a row's Move up/down button is greyed out, in plain
 * language; NULL when the button is enabled.  A move re-orders the record
 * relative to its ADJACENT CHAIN MEMBER (@up: the previous member; else the
 * next member).  The engine re-validates anyway; this mirrors it so the
 * tooltip explains the first blocker.  Shares the generic disables with the
 * Edit/Delete path. */
static const char *hist_move_disabled_reason(const NdeHistRowItem *r,
                                             gboolean up) {
	if (r->dead)
		return _("This step has been undone");
	if (nde_history_is_stale())
		return _("The image was modified outside the recorded history");
	if (processing_is_job_active())
		return _("Not available while processing");
	if (nde_amend_preview_active())
		return _("Another history step is being edited");
	/* A move rewrites the pixel lineage, so it obeys the amend/delete policy:
	 * the record must be amendable and deletable, and both it and the target
	 * neighbour must lie in the editable tail (crossing an opaque barrier is
	 * refused). */
	if (!r->amendable || !r->deletable || !r->in_tail)
		return _("Only steps after the last opaque step can be moved");
	gint64 neighbour = up ? r->prev_member_id : r->next_member_id;
	gboolean neighbour_in_tail = up ? r->prev_in_tail : r->next_in_tail;
	if (!neighbour)
		return up ? _("This is already the first movable step")
		          : _("This is already the last movable step");
	if (!neighbour_in_tail)
		return _("Only steps after the last opaque step can be moved");
	return NULL;
}

static void on_hist_row_activate(GtkListView *lv, guint position, gpointer u) {
	(void)u;
	GListModel *model = G_LIST_MODEL(gtk_list_view_get_model(lv));
	NdeHistRowItem *r = g_list_model_get_item(model, position);
	if (!r) return;
	if (!g_panel->hist_popover) {
		g_panel->hist_popover = gtk_popover_new();
		gtk_widget_set_parent(g_panel->hist_popover, GTK_WIDGET(lv));

		GtkWidget *pbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
		gtk_widget_set_margin_start (pbox, 6);
		gtk_widget_set_margin_end   (pbox, 6);
		gtk_widget_set_margin_top   (pbox, 6);
		gtk_widget_set_margin_bottom(pbox, 6);

		g_panel->hist_popover_label = gtk_label_new(NULL);
		gtk_label_set_selectable(GTK_LABEL(g_panel->hist_popover_label), TRUE);
		gtk_label_set_xalign(GTK_LABEL(g_panel->hist_popover_label), 0.0f);
		gtk_widget_add_css_class(g_panel->hist_popover_label, "monospace");
		gtk_box_append(GTK_BOX(pbox), g_panel->hist_popover_label);

		GtkWidget *actions = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		gtk_widget_set_halign(actions, GTK_ALIGN_END);
		g_panel->hist_move_up_btn   = gtk_button_new_with_label(_("Move up"));
		g_panel->hist_move_down_btn = gtk_button_new_with_label(_("Move down"));
		g_panel->hist_edit_btn   = gtk_button_new_with_label(_("Edit parameters…"));
		g_panel->hist_delete_btn = gtk_button_new_with_label(_("Delete step"));
		gtk_widget_add_css_class(g_panel->hist_delete_btn, "destructive-action");
		gtk_box_append(GTK_BOX(actions), g_panel->hist_move_up_btn);
		gtk_box_append(GTK_BOX(actions), g_panel->hist_move_down_btn);
		gtk_box_append(GTK_BOX(actions), g_panel->hist_edit_btn);
		gtk_box_append(GTK_BOX(actions), g_panel->hist_delete_btn);
		gtk_box_append(GTK_BOX(pbox), actions);

		g_signal_connect(g_panel->hist_move_up_btn,   "clicked",
		                 G_CALLBACK(on_hist_move_up_clicked),   NULL);
		g_signal_connect(g_panel->hist_move_down_btn, "clicked",
		                 G_CALLBACK(on_hist_move_down_clicked), NULL);
		g_signal_connect(g_panel->hist_edit_btn,   "clicked",
		                 G_CALLBACK(on_hist_edit_clicked),   NULL);
		g_signal_connect(g_panel->hist_delete_btn, "clicked",
		                 G_CALLBACK(on_hist_delete_clicked), NULL);

		gtk_popover_set_child(GTK_POPOVER(g_panel->hist_popover), pbox);
	}
	gtk_label_set_text(GTK_LABEL(g_panel->hist_popover_label),
	                   r->detail ? r->detail : "");

	/* Remember the targeted row for the action callbacks. */
	g_clear_object(&g_panel->hist_popover_row);
	g_panel->hist_popover_row = g_object_ref(r);

	const char *edit_off = hist_action_disabled_reason(r, TRUE);
	const char *del_off  = hist_action_disabled_reason(r, FALSE);
	const char *up_off   = hist_move_disabled_reason(r, TRUE);
	const char *down_off = hist_move_disabled_reason(r, FALSE);
	gtk_widget_set_sensitive(g_panel->hist_edit_btn,   edit_off == NULL);
	gtk_widget_set_sensitive(g_panel->hist_delete_btn, del_off == NULL);
	gtk_widget_set_sensitive(g_panel->hist_move_up_btn,   up_off == NULL);
	gtk_widget_set_sensitive(g_panel->hist_move_down_btn, down_off == NULL);
	gtk_widget_set_tooltip_text(g_panel->hist_move_up_btn,   up_off);
	gtk_widget_set_tooltip_text(g_panel->hist_move_down_btn, down_off);
	gtk_widget_set_tooltip_text(g_panel->hist_edit_btn,   edit_off);
	/* Deleting the last opaque step is the moment the freeze feature is
	 * discoverable — spell out that it unlocks the earlier steps. */
	gtk_widget_set_tooltip_text(g_panel->hist_delete_btn,
	                            del_off ? del_off :
	                            (r->last_barrier
	                             ? _("Deleting this opaque step unlocks the steps before it")
	                             : NULL));

	gtk_popover_popup(GTK_POPOVER(g_panel->hist_popover));
	g_object_unref(r);
}

/* Per-record freeze marks derived from the P4.1 chain verdicts (P4.4).
 * frozen/barrier/last_barrier mirror NdeHistRowItem; keyed by record_id.
 * The reorder fields (C5) carry the adjacent same-chain members. */
typedef struct {
	gboolean frozen, barrier, last_barrier;
	gboolean in_tail;
	gint64   prev_member_id, next_member_id;
	gboolean prev_in_tail, next_in_tail;
} hist_freeze_mark;

/* Fold one item's chain into @map (record_id → hist_freeze_mark).  A chain
 * member at index < tail_start is frozen; the member at tail_start-1 that is
 * itself a barrier is "the last barrier" (frozen for Edit, still deletable).
 * A mark is stashed for EVERY member (C5 reorder needs the tail members too);
 * it records each member's adjacent chain members and their tail membership.
 * Non-member records (structural/document-wide) never appear here. */
static void hist_fold_chain_marks(GHashTable *map, const nde_chain *chain) {
	if (!chain || !chain->records)
		return;
	guint tail = chain->tail_start;
	guint len = chain->records->len;
	for (guint i = 0; i < len; i++) {
		const nde_record *rec = g_ptr_array_index(chain->records, i);
		guint8 flags = (chain->member_flags && i < chain->member_flags->len)
		             ? g_array_index(chain->member_flags, guint8, i) : 0;
		gboolean is_barrier = (flags & NDE_CHAIN_MEMBER_BARRIER) != 0;
		hist_freeze_mark *m = g_new0(hist_freeze_mark, 1);
		m->frozen  = (i < tail);
		m->barrier = is_barrier;
		/* The barrier immediately before the tail is the one whose deletion
		 * regains editability of the prefix. */
		m->last_barrier = is_barrier && (i == tail - 1);
		m->in_tail = (i >= tail);
		if (i > 0) {
			const nde_record *p = g_ptr_array_index(chain->records, i - 1);
			m->prev_member_id = p->record_id;
			m->prev_in_tail   = ((i - 1) >= tail);
		}
		if (i + 1 < len) {
			const nde_record *n = g_ptr_array_index(chain->records, i + 1);
			m->next_member_id = n->record_id;
			m->next_in_tail   = ((i + 1) >= tail);
		}
		g_hash_table_insert(map, g_memdup2(&rec->record_id, sizeof(gint64)), m);
	}
}

/* Rebuild the mirror store from a snapshot.  Main thread only. */
static void refresh_history(void) {
	if (!g_panel || !g_panel->hist_pop || !g_panel->hist_store)
		return;
	guint live = 0;
	GPtrArray *snap = nde_history_snapshot_all(&live);
	guint total = snap ? snap->len : 0;
	g_list_store_remove_all(g_panel->hist_store);

	/* Build the freeze-verdict map once per DISTINCT target item present in
	 * the snapshot (records are tens, items a handful): reuse the real P4.1
	 * chain verdicts rather than re-deriving the barrier rules here.  This is
	 * a snapshot+validate pass only (nde_chain_build never replays). */
	GHashTable *freeze = g_hash_table_new_full(g_int64_hash, g_int64_equal,
	                                           g_free, g_free);
	if (snap) {
		GHashTable *seen = g_hash_table_new(g_direct_hash, g_direct_equal);
		for (guint i = 0; i < snap->len; i++) {
			const nde_record *rec = g_ptr_array_index(snap, i);
			gpointer key = GINT_TO_POINTER(rec->target_item_id);
			if (g_hash_table_contains(seen, key))
				continue;
			g_hash_table_add(seen, key);
			nde_chain *chain = nde_chain_build(rec->target_item_id);
			hist_fold_chain_marks(freeze, chain);
			nde_chain_free(chain);
		}
		g_hash_table_destroy(seen);
	}

	if (snap) {
		for (guint i = 0; i < snap->len; i++) {
			const nde_record *rec = g_ptr_array_index(snap, i);
			NdeHistRowItem *item = g_object_new(NDE_TYPE_HIST_ROW_ITEM, NULL);
			item->badge   = g_strdup_printf(rec->tier == NDE_TIER_B ? "#%u·B" :
			                                rec->tier == NDE_TIER_C ? "#%u·C" : "#%u", i + 1);
			item->summary = g_strdup(rec->summary ? rec->summary : rec->op_id);
			item->target  = hist_target_label(rec);
			item->detail  = build_hist_detail(rec);
			item->dead    = (i >= live);
			item->record_id = rec->record_id;
			item->amendable = nde_record_amendable(rec);
			item->deletable = nde_record_deletable(rec);
			item->params    = (rec->params && *rec->params) ? g_strdup(rec->params) : NULL;
			item->op_id     = g_strdup(rec->op_id);
			const hist_freeze_mark *m = g_hash_table_lookup(freeze, &rec->record_id);
			if (m) {
				item->frozen       = m->frozen;
				item->barrier      = m->barrier;
				item->last_barrier = m->last_barrier;
				item->in_tail        = m->in_tail;
				item->prev_member_id = m->prev_member_id;
				item->next_member_id = m->next_member_id;
				item->prev_in_tail   = m->prev_in_tail;
				item->next_in_tail   = m->next_in_tail;
			}
			/* Frozen rows explain why they are locked in the detail text. */
			if (item->frozen && !item->last_barrier) {
				gchar *d = g_strconcat(item->detail, "\n",
				                       _("Locked by a later opaque step"), NULL);
				g_free(item->detail);
				item->detail = d;
			}
			g_list_store_append(g_panel->hist_store, item);
			g_object_unref(item);
		}
		g_ptr_array_unref(snap);
	}
	g_hash_table_destroy(freeze);
	gtk_revealer_set_reveal_child(GTK_REVEALER(g_panel->hist_stale_revealer),
	                              nde_history_is_stale());
	gtk_widget_set_visible(g_panel->hist_fits_hint,
	                       total > 0 && !is_current_image_flis());
	gchar *lbl = g_strdup_printf(_("Nondestructive History (%u)"), live);
	gtk_label_set_text(GTK_LABEL(g_panel->hist_header), lbl);
	g_free(lbl);
}

/* Coalesced cross-thread refresh — same CAS pattern as refresh_pending
 * above; gui_iface.nde_history_changed routes here from any thread. */
static gint history_refresh_pending = 0;

static gboolean history_refresh_idle_cb(gpointer p) {
	(void)p;
	g_atomic_int_set(&history_refresh_pending, 0);
	if (g_panel && g_panel->hist_pop &&
	    gtk_widget_get_visible(g_panel->hist_pop))
		refresh_history();
	return G_SOURCE_REMOVE;
}

void flis_gui_history_update_from_idle(void) {
	if (!g_atomic_int_compare_and_exchange(&history_refresh_pending, 0, 1))
		return;  /* already pending */
	g_idle_add(history_refresh_idle_cb, NULL);
}

/* Refresh when the popover opens — while it is closed the coalesced idle
 * skips refreshes, so the store may be stale. */
static void on_hist_pop_show(GtkWidget *pop, gpointer user_data) {
	(void)pop; (void)user_data;
	refresh_history();
}

/* The Nondestructive History panel: a large popover anchored below the
 * header-bar clock button (nde_history_button).  Separate from the layers
 * panel because the edit history is document-wide, not layer-specific — it
 * applies to plain FITS images too. */
static void build_history_popover(void) {
	if (g_panel->hist_pop) return;

	GtkWidget *pop = gtk_popover_new();
	g_panel->hist_pop = pop;
	gtk_popover_set_position(GTK_POPOVER(pop), GTK_POS_BOTTOM);

	GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
	gtk_widget_set_size_request(vbox, 340, 480);
	gtk_widget_set_margin_start (vbox, 4);
	gtk_widget_set_margin_end   (vbox, 4);
	gtk_widget_set_margin_top   (vbox, 4);
	gtk_widget_set_margin_bottom(vbox, 4);
	gtk_popover_set_child(GTK_POPOVER(pop), vbox);

	g_panel->hist_header = gtk_label_new(_("Nondestructive History"));
	gtk_label_set_xalign(GTK_LABEL(g_panel->hist_header), 0.0f);
	gtk_widget_add_css_class(g_panel->hist_header, "heading");
	gtk_box_append(GTK_BOX(vbox), g_panel->hist_header);

	g_panel->hist_stale_revealer = gtk_revealer_new();
	gtk_revealer_set_transition_type(GTK_REVEALER(g_panel->hist_stale_revealer),
	                                 GTK_REVEALER_TRANSITION_TYPE_SLIDE_DOWN);
	GtkWidget *stale = gtk_label_new(_("Pixels were modified outside the recorded history"));
	gtk_label_set_wrap(GTK_LABEL(stale), TRUE);
	gtk_label_set_xalign(GTK_LABEL(stale), 0.0f);
	gtk_widget_add_css_class(stale, "warning");
	gtk_revealer_set_child(GTK_REVEALER(g_panel->hist_stale_revealer), stale);
	gtk_box_append(GTK_BOX(vbox), g_panel->hist_stale_revealer);

	/* Plain FITS images keep their history in memory only (sketch §13.2);
	 * shown whenever records exist for a non-FLIS image. */
	g_panel->hist_fits_hint = gtk_label_new(
			_("Kept in memory only — save as FLIS to persist this history"));
	gtk_label_set_wrap(GTK_LABEL(g_panel->hist_fits_hint), TRUE);
	gtk_label_set_xalign(GTK_LABEL(g_panel->hist_fits_hint), 0.0f);
	gtk_widget_add_css_class(g_panel->hist_fits_hint, "dim-label");
	gtk_widget_set_visible(g_panel->hist_fits_hint, FALSE);
	gtk_box_append(GTK_BOX(vbox), g_panel->hist_fits_hint);

	g_panel->hist_store = g_list_store_new(NDE_TYPE_HIST_ROW_ITEM);
	GtkSingleSelection *sel = gtk_single_selection_new(
			G_LIST_MODEL(g_object_ref(g_panel->hist_store)));
	gtk_single_selection_set_autoselect(sel, FALSE);
	gtk_single_selection_set_can_unselect(sel, TRUE);

	GtkListItemFactory *factory = gtk_signal_list_item_factory_new();
	g_signal_connect(factory, "setup", G_CALLBACK(on_hist_row_setup), NULL);
	g_signal_connect(factory, "bind",  G_CALLBACK(on_hist_row_bind),  NULL);
	g_panel->hist_list = gtk_list_view_new(GTK_SELECTION_MODEL(sel), factory);
	g_signal_connect(g_panel->hist_list, "activate",
	                 G_CALLBACK(on_hist_row_activate), NULL);

	GtkWidget *sw = gtk_scrolled_window_new();
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
	                               GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(sw), 140);
	gtk_widget_set_vexpand(sw, TRUE);
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(sw), g_panel->hist_list);
	gtk_box_append(GTK_BOX(vbox), sw);

	g_signal_connect(pop, "show", G_CALLBACK(on_hist_pop_show), NULL);
	GtkWidget *btn = lookup_widget("nde_history_button");
	if (btn && GTK_IS_MENU_BUTTON(btn))
		gtk_menu_button_set_popover(GTK_MENU_BUTTON(btn), pop);
}

/* ---- Context menu -------------------------------------------------- */

/* Stub action callback shared by all context menu items.  The body is
 * replaced piecewise as the wiring lands; for now each prints a
 * targeted "TODO: <verb>" log so users see why the menu item didn't
 * do what they expected.  Merge Down and Flatten Image have working
 * primitives (flis_merge_down_layer / flis_flatten_all) so they
 * dispatch through the worker directly. */
static void on_ctx_merge_down(GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_flatten   (GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_move_to_group (GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_export_layer  (GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_layers_match  (GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_register_layers(GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_delete_group  (GSimpleAction *a, GVariant *v, gpointer u);

static GMenu *build_context_menu(void) {
	GMenu *m = g_menu_new();
	g_menu_append(m, _("Export current layer as FITS…"), "win.flis-export-layer");
	g_menu_append(m, _("Register layers…"),              "win.flis-register-layers");
	g_menu_append(m, _("Layers match…"),                 "win.flis-layers-match");
	g_menu_append(m, _("Move layer to group…"),          "win.flis-move-to-group");
	g_menu_append(m, _("Delete group…"),                 "win.flis-delete-group");
	g_menu_append(m, _("Merge Down"),                    "win.flis-merge-down");
	g_menu_append(m, _("Flatten Image"),                 "win.flis-flatten");
	return m;
}

static const GActionEntry flis_panel_actions[] = {
	{ "flis-export-layer",   on_ctx_export_layer, NULL, NULL, NULL },
	{ "flis-register-layers",on_ctx_register_layers, NULL, NULL, NULL },
	{ "flis-layers-match",   on_ctx_layers_match,    NULL, NULL, NULL },
	{ "flis-move-to-group",  on_ctx_move_to_group, NULL, NULL, NULL },
	{ "flis-delete-group",   on_ctx_delete_group,  NULL, NULL, NULL },
	{ "flis-merge-down",     on_ctx_merge_down, NULL, NULL, NULL },
	{ "flis-flatten",        on_ctx_flatten,    NULL, NULL, NULL },
};

static void register_panel_actions(void) {
	/* GTK4 resolves "win.*" action names by walking up from the
	 * widget invoking the action (the GtkMenuButton inside the
	 * panel) to its root window — which is g_panel->window.  Our
	 * panel is a plain GtkWindow (not GtkApplicationWindow) so it
	 * isn't a GActionMap by default; install a GSimpleActionGroup
	 * under the "win" prefix so "win.flis-*" resolves locally. */
	if (!g_panel || !g_panel->window) return;
	GSimpleActionGroup *group = g_simple_action_group_new();
	g_action_map_add_action_entries(G_ACTION_MAP(group),
	                                 flis_panel_actions,
	                                 G_N_ELEMENTS(flis_panel_actions),
	                                 NULL);
	gtk_widget_insert_action_group(g_panel->window, "win",
	                                G_ACTION_GROUP(group));
	g_object_unref(group);
}

/* =========================================================================
 * Refresh: rebuild the GListStore from com.uniq, sync property widgets
 * to the currently-selected layer.
 * ========================================================================= */

static flis_layer_t *current_selected_layer(void) {
	if (!g_panel || !g_panel->list_sel) return NULL;
	guint pos = gtk_single_selection_get_selected(g_panel->list_sel);
	if (pos == GTK_INVALID_LIST_POSITION) return NULL;
	GObject *obj = g_list_model_get_item(G_LIST_MODEL(g_panel->list_store), pos);
	if (!obj) return NULL;
	FlisRowItem *ri = (FlisRowItem *)obj;
	flis_layer_t *lay = row_layer(ri);
	g_object_unref(obj);
	return lay;
}

/* Compute toolbar button sensitivity from current selection + image
 * state.  Called from refresh_panel and on_selection_changed so the
 * buttons re-enable as soon as the user picks a different row, not
 * only after a panel rebuild. */
static void update_toolbar_sensitivity(void) {
	if (!g_panel) return;
	const gboolean is_flis  = is_current_image_flis();
	flis_layer_t  *sel_lay  = current_selected_layer();
	flis_group_t  *sel_grp  = current_selected_group();
	const gboolean have_sel = (sel_lay != NULL) || (sel_grp != NULL);

	gtk_widget_set_sensitive(g_panel->btn_remove,    have_sel);
	gtk_widget_set_sensitive(g_panel->btn_duplicate, sel_lay != NULL);
	gtk_widget_set_sensitive(g_panel->btn_move_up,   have_sel);
	gtk_widget_set_sensitive(g_panel->btn_move_down, have_sel);
	gtk_widget_set_sensitive(g_panel->btn_drag,      sel_lay != NULL);
	gtk_widget_set_sensitive(g_panel->btn_add,       com.uniq != NULL);
	gtk_widget_set_sensitive(g_panel->btn_group,     is_flis);

	/* Resync the drag toggle to the actual mouse mode — the mode can be
	 * switched away from FLIS drag via the mouse-actions UI, which left
	 * the toggle visually stuck on. */
	gboolean drag_on = (mouse_status == MOUSE_ACTION_FLIS_DRAG_LAYER);
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(g_panel->btn_drag)) != drag_on) {
		g_signal_handlers_block_by_func(g_panel->btn_drag, on_drag_toggled, NULL);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g_panel->btn_drag), drag_on);
		g_signal_handlers_unblock_by_func(g_panel->btn_drag, on_drag_toggled, NULL);
	}
}

static flis_group_t *current_selected_group(void) {
	if (!g_panel || !g_panel->list_sel) return NULL;
	guint pos = gtk_single_selection_get_selected(g_panel->list_sel);
	if (pos == GTK_INVALID_LIST_POSITION) return NULL;
	GObject *obj = g_list_model_get_item(G_LIST_MODEL(g_panel->list_store), pos);
	if (!obj) return NULL;
	FlisRowItem *ri = (FlisRowItem *)obj;
	flis_group_t *grp = row_group(ri);
	g_object_unref(obj);
	return grp;
}

gboolean flis_panel_group_is_selected(void) {
	return g_panel && g_panel->selected_kind == FLIS_ROW_KIND_GROUP;
}

static void refresh_panel(void) {
	if (!g_panel) return;
	/* M-F12: the rebuild walks com.uniq->layers/groups and reads layer
	 * fields; hold the stack reader lock so a worker hook can't free a
	 * layer mid-walk.  Main thread only; no fits/cairo locks and no
	 * blocking calls are made inside (GTK widget setters only). */
	flis_stack_reader_lock();
	g_panel->refreshing = TRUE;

	const gboolean is_flis = is_current_image_flis();

	/* Mode label */
	if (!com.uniq) {
		gtk_label_set_text(GTK_LABEL(g_panel->mode_label), _("(no image)"));
	} else if (!is_flis) {
		gtk_label_set_text(GTK_LABEL(g_panel->mode_label), "FITS");
	} else {
		/* Active-group mode (FLIS·group_name) is a §5 addition once
		 * the group-edit-context lives somewhere on com.uniq; for
		 * now the panel just signals plain "FLIS" mode. */
		gtk_label_set_text(GTK_LABEL(g_panel->mode_label), "FLIS");
	}

	/* Remember which row was selected so we can restore it after the
	 * rebuild — otherwise every refresh (triggered by any mutation)
	 * silently drops the user's selection. */
	gint prev_kind = -1;
	gint prev_id   = 0;
	{
		guint psel = gtk_single_selection_get_selected(g_panel->list_sel);
		if (psel != GTK_INVALID_LIST_POSITION) {
			GObject *po = g_list_model_get_item(G_LIST_MODEL(g_panel->list_store), psel);
			if (po) {
				FlisRowItem *pri = (FlisRowItem *)po;
				prev_kind = pri->kind;
				prev_id   = pri->item_id;
				g_object_unref(po);
			}
		}
	}

	/* Rebuild the row model.  Cheap (2-8 layers typical).
	 *
	 * Display convention (Photoshop / GIMP standard):
	 *   - Top of stack (highest layer_order) at the VISUAL top of list
	 *   - Base layer (lowest layer_order) at the VISUAL bottom
	 *
	 * com.uniq->layers is sorted ascending by layer_order — i.e.,
	 * com.uniq->layers->data is the BASE.  To render top-down we walk
	 * back-to-front: collect into an array, then iterate from index
	 * (n-1) down to 0.
	 *
	 * Group rendering: when we cross from one group_id to a different
	 * one (going top-down), we emit a group HEADER row right above
	 * the first child layer we see for that group.  Child layers get
	 * indent_level=1.  A layer with group_id==0 stays flat
	 * (indent_level=0, no header). */
	g_list_store_remove_all(g_panel->list_store);
	if (is_flis && com.uniq && com.uniq->layers) {
		/* Snapshot to a pointer array for index-based reverse walk. */
		const guint n = g_slist_length(com.uniq->layers);
		flis_layer_t **arr = g_new(flis_layer_t *, n);
		guint i = 0;
		for (GSList *l = com.uniq->layers; l; l = l->next, i++)
			arr[i] = (flis_layer_t *)l->data;

		GHashTable *emitted_groups = g_hash_table_new(g_direct_hash, g_direct_equal);

		/* Walk top-down (highest layer_order → lowest). */
		for (gint k = (gint)n - 1; k >= 0; k--) {
			flis_layer_t *lay = arr[k];
			if (!lay) continue;
			const gint gid = lay->group_id;
			flis_group_t *grp = (gid != 0) ? flis_group_get_by_id(gid) : NULL;
			if (gid != 0 &&
			    !g_hash_table_contains(emitted_groups, GINT_TO_POINTER(gid))) {
				FlisRowItem *hri = flis_row_item_new(FLIS_ROW_KIND_GROUP, gid, 0);
				g_list_store_append(g_panel->list_store, hri);
				g_object_unref(hri);
				g_hash_table_insert(emitted_groups, GINT_TO_POINTER(gid),
				                     GINT_TO_POINTER(1));
			}
			/* Skip member rows when the group is collapsed — the header
			 * row's chevron is the user's affordance to expand. */
			if (grp && grp->collapsed) continue;
			FlisRowItem *ri = flis_row_item_new(FLIS_ROW_KIND_LAYER, lay->item_id,
			                                     gid != 0 ? 1 : 0);
			g_list_store_append(g_panel->list_store, ri);
			g_object_unref(ri);
		}
		g_hash_table_destroy(emitted_groups);
		g_free(arr);
	}

	/* Restore the previously-selected row by item_id (matches across
	 * the rebuild — the same layer/group keeps its identity even if
	 * its index in the store moved due to other layers being added /
	 * removed / reordered). */
	if (prev_id != 0) {
		const guint nstore = g_list_model_get_n_items(G_LIST_MODEL(g_panel->list_store));
		for (guint i = 0; i < nstore; i++) {
			GObject *o = g_list_model_get_item(G_LIST_MODEL(g_panel->list_store), i);
			FlisRowItem *ri = (FlisRowItem *)o;
			gboolean match = (ri->kind == prev_kind && ri->item_id == prev_id);
			g_object_unref(o);
			if (match) {
				gtk_single_selection_set_selected(g_panel->list_sel, i);
				break;
			}
		}
	}

	/* Sync property widgets to selection. */
	sync_property_widgets(current_selected_layer());

	/* Sensitivity of toolbar buttons. */
	update_toolbar_sensitivity();

	g_panel->refreshing = FALSE;
	flis_stack_reader_unlock();

	/* History section mirrors the provenance log; refreshed with the stack
	 * lock released (its snapshot takes only the nde leaf mutex; target
	 * name resolution re-reads the live stack lock-free like the rest of
	 * the panel's main-thread accessors). */
	refresh_history();
}

/* Set the property panel to a "nothing selected" state. */
static void sync_property_widgets_empty(void) {
	gtk_widget_set_sensitive(g_panel->prop_frame, FALSE);
	gtk_widget_set_sensitive(g_panel->mask_frame, FALSE);
	gtk_editable_set_text(GTK_EDITABLE(g_panel->name_entry), "");
	gtk_button_set_label(GTK_BUTTON(g_panel->mask_status_btn), _("(no mask)"));
	gtk_widget_set_sensitive(g_panel->tint_frame, FALSE);
}

/* Fill the property panel widgets for a selected LAYER. */
static void sync_property_widgets_for_layer(flis_layer_t *lay) {
	gtk_widget_set_sensitive(g_panel->prop_frame, TRUE);
	gtk_widget_set_sensitive(g_panel->mask_frame, TRUE);

	gtk_editable_set_text(GTK_EDITABLE(g_panel->name_entry),
	                      lay->layer_name ? lay->layer_name : "");
	gtk_drop_down_set_selected(GTK_DROP_DOWN(g_panel->blend_dropdown),
	                            index_for_blend_mode(lay->blend_mode));
	gtk_adjustment_set_value(g_panel->opacity_adj, (double)(lay->opacity * 100.0f));

	/* Tint sub-frame sensitivity: only mono layers can be tinted. */
	const gboolean is_mono = (lay->fit && lay->fit->naxes[2] == 1);
	gtk_widget_set_sensitive(g_panel->tint_frame, is_mono);
	gtk_check_button_set_active(GTK_CHECK_BUTTON(g_panel->tint_check), lay->has_tint);
	if (lay->has_tint) {
		GdkRGBA rgba = {
			(float)lay->layer_tint.r,
			(float)lay->layer_tint.g,
			(float)lay->layer_tint.b,
			1.0f
		};
		gtk_color_dialog_button_set_rgba(GTK_COLOR_DIALOG_BUTTON(g_panel->tint_color_btn), &rgba);
	}

	/* Mask sub-frame */
	if (lay->lmask) {
		gchar *txt = g_strdup_printf(_("Layer mask (%zux%zu, %s)"),
			lay->lmask->w, lay->lmask->h,
			lay->lmask_active ? _("active") : _("inactive"));
		gtk_button_set_label(GTK_BUTTON(g_panel->mask_status_btn), txt);
		g_free(txt);
		gtk_button_set_label(GTK_BUTTON(g_panel->mask_toggle_btn), _("Remove…"));
		gtk_widget_set_sensitive(g_panel->mask_move_btn, TRUE);
	} else {
		gtk_button_set_label(GTK_BUTTON(g_panel->mask_status_btn), _("(no mask)"));
		gtk_button_set_label(GTK_BUTTON(g_panel->mask_toggle_btn), _("Add…"));
		gtk_widget_set_sensitive(g_panel->mask_move_btn, FALSE);
	}

	/* View-row visible only when both a processing mask and a layer
	 * mask exist — picks which one the global mask tab and the tint
	 * overlay on the image vports source from. */
	const gboolean has_proc_mask = (gfit && gfit->mask && gfit->mask->data);
	const gboolean has_lmask     = (lay->lmask && lay->lmask->data);
	const gboolean both_masks    = has_proc_mask && has_lmask;
	gtk_widget_set_visible(g_panel->mask_view_row, both_masks);

	/* Sync the radio state to com.uniq's view enum without firing
	 * the toggled handlers (which would loop back into the model). */
	if (both_masks && com.uniq) {
		GtkCheckButton *want = (com.uniq->flis_mask_view == 1)
		    ? GTK_CHECK_BUTTON(g_panel->mask_view_layer_radio)
		    : GTK_CHECK_BUTTON(g_panel->mask_view_proc_radio);
		if (!siril_toggle_get_active(GTK_WIDGET(want))) {
			g_signal_handlers_block_by_func(
			    g_panel->mask_view_proc_radio,
			    on_mask_view_radio_toggled, GINT_TO_POINTER(0));
			g_signal_handlers_block_by_func(
			    g_panel->mask_view_layer_radio,
			    on_mask_view_radio_toggled, GINT_TO_POINTER(1));
			gtk_check_button_set_active(want, TRUE);
			g_signal_handlers_unblock_by_func(
			    g_panel->mask_view_proc_radio,
			    on_mask_view_radio_toggled, GINT_TO_POINTER(0));
			g_signal_handlers_unblock_by_func(
			    g_panel->mask_view_layer_radio,
			    on_mask_view_radio_toggled, GINT_TO_POINTER(1));
		}
	}
}

/* Fill the property panel widgets for a selected GROUP.  Groups support
 * name / blend mode / opacity (the same three widgets as layers).
 * Tint and mask are layer-only concepts — disable / hide them. */
static void sync_property_widgets_for_group(flis_group_t *grp) {
	gtk_widget_set_sensitive(g_panel->prop_frame, TRUE);
	gtk_widget_set_sensitive(g_panel->mask_frame, FALSE);

	gtk_editable_set_text(GTK_EDITABLE(g_panel->name_entry),
	                      grp->name ? grp->name : "");
	gtk_drop_down_set_selected(GTK_DROP_DOWN(g_panel->blend_dropdown),
	                            index_for_blend_mode(grp->blend_mode));
	gtk_adjustment_set_value(g_panel->opacity_adj, (double)(grp->opacity * 100.0f));

	/* Tint sub-frame: not applicable to groups. */
	gtk_widget_set_sensitive(g_panel->tint_frame, FALSE);
	gtk_check_button_set_active(GTK_CHECK_BUTTON(g_panel->tint_check), FALSE);

	/* Mask sub-frame: groups don't carry their own lmask. */
	gtk_button_set_label(GTK_BUTTON(g_panel->mask_status_btn), _("(groups have no mask)"));
	gtk_widget_set_visible(g_panel->mask_view_row, FALSE);
}

/* Sync the property panel to current selection (layer, group, or none).
 * Guards against the cascade of widget "changed" signals firing
 * dispatch_op while we're updating display state. */
static void sync_property_widgets(flis_layer_t *lay) {
	const gboolean was_refreshing = g_panel->refreshing;
	g_panel->refreshing = TRUE;

	if (lay) {
		sync_property_widgets_for_layer(lay);
	} else {
		flis_group_t *grp = current_selected_group();
		if (grp) sync_property_widgets_for_group(grp);
		else     sync_property_widgets_empty();
	}

	g_panel->refreshing = was_refreshing;
}

/* =========================================================================
 * Operation dispatch helpers
 *
 * Every panel mutation routes through generic_layer_worker (§1.5).  The
 * hook is a tiny closure over (target_item_id, op_kind, payload) — when
 * com.script is set the worker short-circuits the undo gate, matching
 * the headless command-driven path.
 *
 * For property-only mutations (visibility, blend, opacity, tint), the
 * hook calls the matching flis_layer_set_* primitive.  The dispatch
 * helper sets invalidate_flags to the narrowest applicable
 * FLIS_INV_* combination (§3.4).
 * ========================================================================= */

enum op_kind {
	OP_SET_VISIBLE, OP_SET_LOCKED, OP_SET_NAME, OP_SET_BLEND,
	OP_SET_OPACITY, OP_SET_TINT, OP_CLEAR_TINT,
	OP_SET_LMASK_ACTIVE, OP_UNGROUP,
	OP_GROUP_SET_VISIBLE, OP_GROUP_SET_NAME, OP_GROUP_SET_BLEND, OP_GROUP_SET_OPACITY,
	OP_LAYER_REMOVE, OP_LAYER_DUPLICATE,
	OP_LAYER_MOVE_UP, OP_LAYER_MOVE_DOWN,
};

struct op_payload {
	destructor    destroy_fn;
	gint          target_id;
	enum op_kind  kind;
	gboolean      bool_v;
	gchar        *str_v;
	flis_blend_mode_t blend_v;
	gfloat        float_v;
	double        r, g, b;
	gint          dup_item_id;   /* out: new layer id set by OP_LAYER_DUPLICATE */
};

static void op_payload_free(gpointer p) {
	struct op_payload *op = p;
	if (!op) return;
	g_free(op->str_v);
	g_free(op);
}

static int op_hook(struct generic_layer_args *args) {
	struct op_payload *op = (struct op_payload *)args->user;
	const gboolean is_group_op =
		(op->kind == OP_GROUP_SET_VISIBLE ||
		 op->kind == OP_GROUP_SET_NAME    ||
		 op->kind == OP_GROUP_SET_BLEND   ||
		 op->kind == OP_GROUP_SET_OPACITY);
	flis_layer_t *lay = is_group_op ? NULL : flis_layer_get_by_id(op->target_id);
	flis_group_t *grp = is_group_op ? flis_group_get_by_id(op->target_id) : NULL;

	/* Remove frees the layer, so snapshot its name before running the op. */
	gchar *removed_name = (op->kind == OP_LAYER_REMOVE && lay && lay->layer_name)
	                      ? g_strdup(lay->layer_name) : NULL;

	int rc;
	switch (op->kind) {
		case OP_SET_VISIBLE:  rc = lay ? flis_layer_set_visible (lay, op->bool_v) : 1; break;
		case OP_SET_LOCKED:   rc = lay ? flis_layer_set_locked  (lay, op->bool_v) : 1; break;
		case OP_SET_NAME:     rc = lay ? flis_layer_set_name    (lay, op->str_v) : 1; break;
		case OP_SET_BLEND:    rc = lay ? flis_layer_set_blend_mode(lay, op->blend_v) : 1; break;
		case OP_SET_OPACITY:  rc = lay ? flis_layer_set_opacity (lay, op->float_v) : 1; break;
		case OP_SET_TINT:     rc = lay ? flis_layer_set_tint    (lay, op->r, op->g, op->b) : 1; break;
		case OP_CLEAR_TINT:   rc = lay ? flis_layer_clear_tint(lay) : 1; break;
		case OP_SET_LMASK_ACTIVE:
			if (!lay || !lay->lmask) { rc = 1; break; }
			lay->lmask_active = op->bool_v;
			flis_layer_touch_modified(lay);
			rc = 0; break;
		case OP_UNGROUP:      rc = lay ? flis_layer_set_group(lay, 0) : 1; break;
		case OP_GROUP_SET_VISIBLE: rc = grp ? flis_group_set_visible(grp, op->bool_v) : 1; break;
		case OP_GROUP_SET_NAME:    rc = grp ? flis_group_set_name(grp, op->str_v)    : 1; break;
		case OP_GROUP_SET_BLEND:   rc = grp ? flis_group_set_blend_mode(grp, op->blend_v) : 1; break;
		case OP_GROUP_SET_OPACITY: rc = grp ? flis_group_set_opacity(grp, op->float_v)   : 1; break;
		case OP_LAYER_REMOVE: rc = lay ? flis_layer_remove(lay) : 1; break;
		case OP_LAYER_DUPLICATE: {
			if (!lay) { rc = 1; break; }
			flis_layer_t *dup = flis_layer_duplicate(lay);
			op->dup_item_id = dup ? dup->item_id : 0;
			rc = dup ? 0 : 1; break;
		}
		case OP_LAYER_MOVE_UP:   rc = lay ? flis_layer_move_up(lay)   : 1; break;
		case OP_LAYER_MOVE_DOWN: rc = lay ? flis_layer_move_down(lay) : 1; break;
		default: rc = 1; break;
	}

	/* NDE provenance capture (sketch §13.2) — intent-site: this GUI funnel is
	 * reached only from user intent, whereas the low-level flis_layer_set_*
	 * setters are ALSO called by load and by test fixtures, so capture must
	 * NOT live in them.  Capture here, after the mutation succeeded.
	 * move_up/move_down/set_name/set_locked/set_tint/lmask-active/ungroup and
	 * the group ops are deliberately excluded here: move_up/down self-record
	 * in the backend; the rest are display attributes not recorded in phase 1
	 * (sketch §13.2 records visibility but not lock/tint/rename/position).
	 *
	 * Undo coupling for the three property commits: dispatch_op saved a
	 * props-only undo entry BEFORE the worker ran (gated on !com.script), so
	 * by now that entry is on top of the undo stack.  Couple it here with the
	 * same !com.script guard the undo save used (undo_save_flis_layer_props
	 * returns 0 WITHOUT pushing when com.script — sketch §13.3), mirroring
	 * generic_image_worker's post-hook tagging from the worker thread. */
	if (rc == 0) {
		switch (op->kind) {
			case OP_SET_OPACITY: {
				GString *kv = nde_kv_start();
				nde_kv_add_float(kv, "opacity", op->float_v);
				gint64 rid = nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER,
				                       op->target_id, nde_kv_end(kv), _("Set opacity"));
				if (!com.script) undo_tag_top_nde_record(rid);
				break;
			}
			case OP_SET_BLEND: {
				/* blend value is the frozen flis_blend_mode_t enum (int). */
				GString *kv = nde_kv_start();
				nde_kv_add_int(kv, "blend", op->blend_v);
				gint64 rid = nde_capture_structural("layer.set_blend", NDE_SCOPE_LAYER,
				                       op->target_id, nde_kv_end(kv), _("Set blend mode"));
				if (!com.script) undo_tag_top_nde_record(rid);
				break;
			}
			case OP_SET_VISIBLE: {
				GString *kv = nde_kv_start();
				nde_kv_add_bool(kv, "visible", op->bool_v);
				gint64 rid = nde_capture_structural("layer.set_visible", NDE_SCOPE_LAYER,
				                       op->target_id, nde_kv_end(kv), _("Set visibility"));
				if (!com.script) undo_tag_top_nde_record(rid);
				break;
			}
			case OP_LAYER_DUPLICATE: {
				GString *kv = nde_kv_start();
				nde_kv_add_str(kv, "name", "");
				nde_kv_add_int(kv, "src_item", op->target_id);
				nde_capture_structural("layer.duplicate", NDE_SCOPE_DOCUMENT,
				                       op->dup_item_id, nde_kv_end(kv), _("Duplicate layer"));
				break;
			}
			case OP_LAYER_REMOVE: {
				GString *kv = nde_kv_start();
				nde_kv_add_str(kv, "name", removed_name ? removed_name : "");
				nde_capture_structural("layer.remove", NDE_SCOPE_DOCUMENT,
				                       op->target_id, nde_kv_end(kv), _("Remove layer"));
				break;
			}
			default: break;
		}
	}
	g_free(removed_name);
	return rc;
}

static void dispatch_op(struct op_payload *op, const char *desc,
                         int invalidate_flags) {
	/* §C.1a: save undo BEFORE the mutation runs.  Pure-property layer
	 * ops use the props-only undo entry (no swap file, just a
	 * snapshot of the layer's pre-change props).  The undo function
	 * already short-circuits when com.script is set, so headless and
	 * Python script runs see no undo as designed.  Group property
	 * ops and structural ops (remove / duplicate / move) are not
	 * undoable here yet — TODO follow-up. */
	switch (op->kind) {
		case OP_SET_VISIBLE:
		case OP_SET_LOCKED:
		case OP_SET_NAME:
		case OP_SET_BLEND:
		case OP_SET_OPACITY:
		case OP_SET_TINT:
		case OP_CLEAR_TINT:
		case OP_SET_LMASK_ACTIVE: {
			flis_layer_t *lay = flis_layer_get_by_id(op->target_id);
			if (lay) undo_save_flis_layer_props(lay,
			                                    desc ? desc : "Layer op");
			break;
		}
		case OP_GROUP_SET_VISIBLE:
		case OP_GROUP_SET_NAME:
		case OP_GROUP_SET_BLEND:
		case OP_GROUP_SET_OPACITY:
			/* No group-props undo record type exists yet (TODO above);
			 * until then, at least tell the user instead of silently
			 * doing nothing on Ctrl-Z. */
			if (!com.script)
				siril_log_message(_("Note: group property changes cannot be undone yet\n"));
			break;
		default:
			break;
	}

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer        = NULL;     /* hook resolves layer from op->target_id */
	args->layer_hook   = op_hook;
	args->user         = op;
	args->description  = g_strdup(desc);
	args->verbose      = FALSE;
	args->invalidate_flags   = invalidate_flags;
	args->invalidate_item_id = op->target_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

/* Ungroup helper usable from code that sits above the op_payload
 * definition (the drag-and-drop edge handler). */
static void dispatch_ungroup(gint layer_id) {
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id  = layer_id;
	op->kind       = OP_UNGROUP;
	dispatch_op(op, _("Remove from group"), FLIS_INV_STACK);
}

/* =========================================================================
 * Signal handlers
 * ========================================================================= */

static void on_row_visible_toggled(GtkToggleButton *btn, gpointer u) {
	if (g_panel && g_panel->refreshing) return;
	gint encoded = GPOINTER_TO_INT(u);
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->bool_v = gtk_toggle_button_get_active(btn);
	if (encoded < 0) {
		op->target_id = -encoded;
		op->kind = OP_GROUP_SET_VISIBLE;
		dispatch_op(op, _("Group visibility"), FLIS_INV_STACK);
	} else {
		op->target_id = encoded;
		op->kind = OP_SET_VISIBLE;
		dispatch_op(op, _("Layer visibility"), FLIS_INV_STACK);
	}
}

static void on_row_lock_toggled(GtkToggleButton *btn, gpointer u) {
	if (g_panel && g_panel->refreshing) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = GPOINTER_TO_INT(u);
	op->kind = OP_SET_LOCKED;
	op->bool_v = gtk_toggle_button_get_active(btn);
	/* Locked is editorial state — no display invalidation needed but
	 * the panel needs a refresh; use LAYER_PROPS as the narrowest
	 * flag that still triggers the refresh cycle. */
	dispatch_op(op, _("Layer locked"), FLIS_INV_LAYER_PROPS);
}

/* Group header chevron toggle.  Pure UI state — flips
 * flis_group_t.collapsed and refreshes the panel.  No worker dispatch
 * (the collapse state is panel-only; the composite is unaffected). */
static void on_group_expander_clicked(GtkButton *btn, gpointer u) {
	(void)btn;
	if (g_panel && g_panel->refreshing) return;
	const gint gid = GPOINTER_TO_INT(u);
	flis_group_t *grp = flis_group_get_by_id(gid);
	if (!grp) return;
	grp->collapsed = !grp->collapsed;
	refresh_panel();
}

static void on_selection_changed(GtkSelectionModel *sel, guint pos, guint nitems, gpointer u) {
	(void)sel; (void)pos; (void)nitems; (void)u;
	if (g_panel && g_panel->refreshing) return;
	/* An un-applied opacity preview does not survive a selection
	 * change — put the outgoing item's opacity back. */
	opacity_edit_revert();
	flis_layer_t *lay = current_selected_layer();
	flis_group_t *grp = current_selected_group();
	if (lay)      { g_panel->selected_kind = FLIS_ROW_KIND_LAYER; g_panel->selected_item_id = lay->item_id; }
	else if (grp) { g_panel->selected_kind = FLIS_ROW_KIND_GROUP; g_panel->selected_item_id = grp->item_id; }
	else          { g_panel->selected_kind = -1;                  g_panel->selected_item_id = 0; }
	sync_property_widgets(lay);
	update_toolbar_sensitivity();   /* re-enable Remove/Duplicate/etc. */

	/* §4.5 checkoff: selecting a layer rebinds gfit (and uniq->fit,
	 * uniq->chans, uniq->active_layer) so tools that operate on gfit
	 * — asinh, GHS, MTF, histogram, etc. — act on the right layer.
	 * Composite cache is intentionally NOT rebuilt; the display still
	 * shows the composite, so the active-layer switch is invisible
	 * on screen until a tool touches gfit. */
	if (lay && com.uniq) {
		gint index = 0;
		gboolean found = FALSE;
		for (GSList *l = com.uniq->layers; l; l = l->next, index++) {
			if ((flis_layer_t *)l->data == lay) { found = TRUE; break; }
		}
		if (found && index != com.uniq->active_layer) {
			/* Preview/ROI-coherent switch: reverts a live preview on
			 * the outgoing layer and re-arms it on this one, and
			 * refreshes the ROI pixel cache (see
			 * flis_switch_active_layer_gui). */
			flis_switch_active_layer_gui(index);
			gui_iface.redraw_image(REDRAW_ALL);
			/* Refresh so the active-layer-row CSS class moves to the
			 * newly active row.  Selection is preserved across the
			 * rebuild by refresh_panel's snapshot/restore. */
			refresh_panel();
		}
	}
}

/* Property handlers */
static void on_name_focus_leave(GtkEventControllerFocus *c, gpointer u) {
	(void)c; (void)u;
	if (g_panel && g_panel->name_entry)
		on_name_activate(GTK_ENTRY(g_panel->name_entry), NULL);
}

static void on_name_activate(GtkEntry *e, gpointer u) {
	(void)u;
	if (g_panel && g_panel->refreshing) return;
	if (g_panel->selected_item_id == 0) return;
	/* Skip no-op renames — this also runs from the focus-leave
	 * controller, which fires on every focus change. */
	const char *txt = gtk_editable_get_text(GTK_EDITABLE(e));
	const char *cur = NULL;
	if (g_panel->selected_kind == FLIS_ROW_KIND_GROUP) {
		flis_group_t *grp = flis_group_get_by_id(g_panel->selected_item_id);
		cur = grp ? grp->name : NULL;
	} else {
		flis_layer_t *lay = flis_layer_get_by_id(g_panel->selected_item_id);
		cur = lay ? lay->layer_name : NULL;
	}
	if (!g_strcmp0(txt, cur ? cur : ""))
		return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = g_panel->selected_item_id;
	op->kind = (g_panel->selected_kind == FLIS_ROW_KIND_GROUP)
	             ? OP_GROUP_SET_NAME : OP_SET_NAME;
	op->str_v = g_strdup(gtk_editable_get_text(GTK_EDITABLE(e)));
	dispatch_op(op,
		(g_panel->selected_kind == FLIS_ROW_KIND_GROUP) ? _("Group name") : _("Layer name"),
		FLIS_INV_LAYER_PROPS);
}

static void on_blend_changed(GtkDropDown *dd, GParamSpec *p, gpointer u) {
	(void)p; (void)u;
	if (g_panel && g_panel->refreshing) return;
	if (g_panel->selected_item_id == 0) return;
	int idx = (int)gtk_drop_down_get_selected(dd);
	if (idx < 0 || idx >= N_BLENDS) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = g_panel->selected_item_id;
	op->kind = (g_panel->selected_kind == FLIS_ROW_KIND_GROUP)
	             ? OP_GROUP_SET_BLEND : OP_SET_BLEND;
	op->blend_v = blend_mode_for_index_table[idx];
	dispatch_op(op,
		(g_panel->selected_kind == FLIS_ROW_KIND_GROUP) ? _("Group blend mode") : _("Layer blend mode"),
		FLIS_INV_LAYER_PROPS);
}

/* Coalesced display refresh for per-motion drag ticks (opacity slider,
 * canvas layer drag).  Running the full notify_gfit_data_modified
 * pipeline (stats + histogram + CPU composite rebuild + remap) on every
 * motion event stalls the main thread on large canvases; queueing it at
 * idle priority collapses a burst of ticks into one rebuild once the
 * event queue drains, while an immediate paint keeps the GPU compose
 * path (which reads live property values per frame) fluid. */
static gint drag_refresh_pending = 0;
static gboolean drag_refresh_idle_cb(gpointer p) {
	(void)p;
	g_atomic_int_set(&drag_refresh_pending, 0);
	notify_gfit_data_modified();
	gui_iface.redraw_image(REDRAW_ALL);
	return G_SOURCE_REMOVE;
}
void flis_drag_tick_refresh(void) {
	gui_iface.redraw_image(REDRAW_ALL);
	if (!g_atomic_int_compare_and_exchange(&drag_refresh_pending, 0, 1))
		return;
	g_idle_add_full(G_PRIORITY_DEFAULT_IDLE, drag_refresh_idle_cb, NULL, NULL);
}

/* Opacity edits are "live preview" only: slider / spin-button ticks
 * write the property directly (plus the coalesced display refresh) and
 * never touch the worker or the undo stack.  The Apply button next to
 * the spin button commits the edit: it restores the pre-edit value and
 * dispatches a single worker op, so the undo snapshot taken by
 * dispatch_op records the true pre-edit state.  Selecting a different
 * layer or group without applying reverts the preview.
 *
 * (Historical note: dispatching a worker op per adjustment tick tripped
 * start_in_new_thread's busy guard — "The processing thread is busy" —
 * which spammed the log and dropped the refused ops; a per-tick commit
 * also produced one undo entry per tick.) */

/* Un-applied edit; id == 0 means none pending. */
static struct {
	gint     id;
	gboolean is_group;
	gfloat   start;      /* property value before the edit, for undo/revert */
} opacity_edit = { 0, FALSE, 0.0f };

static void opacity_edit_clear(void) {
	opacity_edit.id = 0;
	if (g_panel && g_panel->opacity_apply)
		gtk_widget_set_sensitive(g_panel->opacity_apply, FALSE);
}

/* Discard an un-applied edit, restoring the pre-edit value. */
static void opacity_edit_revert(void) {
	if (opacity_edit.id == 0) return;
	gboolean changed = FALSE;
	if (opacity_edit.is_group) {
		flis_group_t *grp = flis_group_get_by_id(opacity_edit.id);
		if (grp && grp->opacity != opacity_edit.start) {
			grp->opacity = opacity_edit.start;
			changed = TRUE;
		}
	} else {
		flis_layer_t *lay = flis_layer_get_by_id(opacity_edit.id);
		if (lay && lay->opacity != opacity_edit.start) {
			lay->opacity = opacity_edit.start;
			changed = TRUE;
		}
	}
	if (changed) {
		gui_iface.flis_display_invalidate(FLIS_INV_LAYER_PROPS, opacity_edit.id);
		flis_drag_tick_refresh();
	}
	opacity_edit_clear();
}

/* Commit payload captured at Apply-click time, so a subsequent
 * selection change (which reverts the *pending* edit) cannot disturb an
 * apply that is waiting out a busy worker. */
struct opacity_apply {
	gint     id;
	gboolean is_group;
	gfloat   start, final;
};

static gboolean opacity_apply_cb(gpointer p) {
	struct opacity_apply *a = (struct opacity_apply *)p;
	if (processing_is_job_active())
		return G_SOURCE_CONTINUE;   /* worker busy: retry, don't drop */
	/* Restore the pre-edit value so the worker-side undo snapshot
	 * records it; the worker re-applies the final value under the
	 * stack lock. */
	if (a->is_group) {
		flis_group_t *grp = flis_group_get_by_id(a->id);
		if (!grp) goto DONE;
		grp->opacity = a->start;
	} else {
		flis_layer_t *lay = flis_layer_get_by_id(a->id);
		if (!lay) goto DONE;
		lay->opacity = a->start;
	}
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = a->id;
	op->kind = a->is_group ? OP_GROUP_SET_OPACITY : OP_SET_OPACITY;
	op->float_v = a->final;
	dispatch_op(op,
		a->is_group ? _("Group opacity") : _("Layer opacity"),
		FLIS_INV_LAYER_PROPS);
DONE:
	g_free(a);
	return G_SOURCE_REMOVE;
}

static void on_opacity_apply_clicked(GtkButton *btn, gpointer u) {
	(void)btn; (void)u;
	if (opacity_edit.id == 0) return;
	gfloat final_v;
	if (opacity_edit.is_group) {
		flis_group_t *grp = flis_group_get_by_id(opacity_edit.id);
		if (!grp) { opacity_edit_clear(); return; }
		final_v = grp->opacity;
	} else {
		flis_layer_t *lay = flis_layer_get_by_id(opacity_edit.id);
		if (!lay) { opacity_edit_clear(); return; }
		final_v = lay->opacity;
	}
	if (final_v == opacity_edit.start) {
		opacity_edit_clear();
		return;
	}
	struct opacity_apply *a = g_new0(struct opacity_apply, 1);
	a->id       = opacity_edit.id;
	a->is_group = opacity_edit.is_group;
	a->start    = opacity_edit.start;
	a->final    = final_v;
	opacity_edit_clear();
	/* Immediate first attempt; falls back to polling only if the
	 * worker is busy right now. */
	if (opacity_apply_cb(a) == G_SOURCE_CONTINUE)
		g_timeout_add(100, opacity_apply_cb, a);
}

static void on_opacity_changed(GtkAdjustment *adj, gpointer u) {
	(void)u;
	if (!g_panel || g_panel->refreshing) return;
	if (g_panel->selected_item_id == 0) return;
	const gboolean is_group = (g_panel->selected_kind == FLIS_ROW_KIND_GROUP);
	const gint id = g_panel->selected_item_id;
	const gfloat v = (gfloat)(gtk_adjustment_get_value(adj) / 100.0);

	/* Shouldn't happen (selection change reverts first), but never let
	 * an edit of one item adopt another item's start value. */
	if (opacity_edit.id != 0 &&
	    (opacity_edit.id != id || opacity_edit.is_group != is_group))
		opacity_edit_revert();

	if (is_group) {
		flis_group_t *grp = flis_group_get_by_id(id);
		if (!grp) return;
		if (opacity_edit.id == 0)
			opacity_edit.start = grp->opacity;
		grp->opacity = v;
	} else {
		flis_layer_t *lay = flis_layer_get_by_id(id);
		if (!lay) return;
		if (opacity_edit.id == 0)
			opacity_edit.start = lay->opacity;
		lay->opacity = v;
	}
	opacity_edit.id = id;
	opacity_edit.is_group = is_group;
	gtk_widget_set_sensitive(g_panel->opacity_apply, v != opacity_edit.start);
	gui_iface.flis_display_invalidate(FLIS_INV_LAYER_PROPS, id);
	flis_drag_tick_refresh();
}

static void on_tint_check_toggled(GtkCheckButton *btn, gpointer u) {
	(void)u;
	if (g_panel && g_panel->refreshing) return;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	gboolean want = gtk_check_button_get_active(btn);
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = lay->item_id;
	if (want) {
		op->kind = OP_SET_TINT;
		/* Apply the colour the swatch is showing — enabling with the
		 * neutral broadcast while a colour is visibly selected was a
		 * "nothing happened" surprise. */
		const GdkRGBA *rgba = (g_panel && g_panel->tint_color_btn)
			? gtk_color_dialog_button_get_rgba(
				GTK_COLOR_DIALOG_BUTTON(g_panel->tint_color_btn))
			: NULL;
		if (rgba) {
			op->r = rgba->red;
			op->g = rgba->green;
			op->b = rgba->blue;
		} else {
			op->r = op->g = op->b = 1.0;
		}
	} else {
		op->kind = OP_CLEAR_TINT;
	}
	dispatch_op(op, _("Layer tint"), FLIS_INV_LAYER_PIXELS);
}

static void on_tint_color_chosen(GtkColorDialogButton *btn, GParamSpec *p, gpointer u) {
	(void)p; (void)u;
	if (g_panel && g_panel->refreshing) return;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	const GdkRGBA *rgba = gtk_color_dialog_button_get_rgba(btn);
	if (!rgba) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = lay->item_id;
	op->kind = OP_SET_TINT;
	op->r = rgba->red; op->g = rgba->green; op->b = rgba->blue;
	dispatch_op(op, _("Layer tint color"), FLIS_INV_LAYER_PIXELS);
}

/* Toolbar handlers — these are stubs that print a TODO until the
 * matching flis_* commands ship in §4.3.  The structural wiring is in
 * place; only the body of each handler awaits its command. */
/* GtkFileDialog completion for the Add Layer toolbar button.  Builds
 * the same flis_addlayer_args + generic_layer_args struct the
 * flis_addlayer command builds, then submits to generic_layer_worker
 * — guaranteed parity with the command path. */
static void on_add_file_chosen(GObject *src, GAsyncResult *res, gpointer ud) {
	(void)ud;
	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	GError *err = NULL;
	GFile *file = gtk_file_dialog_open_finish(fd, res, &err);
	if (!file) {
		if (err && !g_error_matches(err, GTK_DIALOG_ERROR, GTK_DIALOG_ERROR_DISMISSED))
			siril_log_error(_("FLIS: Add Layer: %s\n"), err->message);
		g_clear_error(&err);
		return;
	}
	gchar *path = g_file_get_path(file);
	g_object_unref(file);
	if (!path) return;

	struct flis_addlayer_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addlayer_args_free;
	payload->filename   = path;     /* takes ownership */
	payload->name       = NULL;     /* hook derives from basename */

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_addlayer_hook;
	args->user               = payload;
	args->description        = g_strdup(_("Add layer"));
	args->invalidate_flags   = FLIS_INV_ALL;
	args->invalidate_item_id = 0;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

static void on_add_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	if (!is_current_image_flis()) {
		siril_log_message(_("FLIS: Add Layer requires a FLIS image (use flis_promote on a plain FITS first)\n"));
		return;
	}
	GtkFileDialog *fd = gtk_file_dialog_new();
	gtk_file_dialog_set_title(fd, _("Add layer from FITS file"));
	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_set_name(filter, _("FITS files"));
	gtk_file_filter_add_pattern(filter, "*.fit");
	gtk_file_filter_add_pattern(filter, "*.fits");
	gtk_file_filter_add_pattern(filter, "*.fts");
	gtk_file_filter_add_pattern(filter, "*.FIT");
	gtk_file_filter_add_pattern(filter, "*.FITS");
	GListStore *filters = g_list_store_new(GTK_TYPE_FILE_FILTER);
	g_list_store_append(filters, filter);
	gtk_file_dialog_set_default_filter(fd, filter);
	gtk_file_dialog_set_filters(fd, G_LIST_MODEL(filters));
	g_object_unref(filters);
	g_object_unref(filter);
	gtk_file_dialog_open(fd, g_panel ? GTK_WINDOW(g_panel->window) : NULL,
	                      NULL, on_add_file_chosen, NULL);
	g_object_unref(fd);
}
/* Confirmation dialog callback for Remove Layer.  Remove is destructive
 * and not undoable on this branch (the §1.4 undo plumbing for layer
 * removal was deemed too costly given typical FLIS sizes), so a
 * GtkAlertDialog warns before dispatching. */
static void on_remove_confirm_resp(GObject *src, GAsyncResult *res, gpointer ud) {
	GtkAlertDialog *ad = GTK_ALERT_DIALOG(src);
	gint chosen = gtk_alert_dialog_choose_finish(ad, res, NULL);
	const gint id = GPOINTER_TO_INT(ud);
	/* Button 0 = Remove, 1 = Cancel.  Esc / close → -1. */
	if (chosen != 0) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = id;
	op->kind = OP_LAYER_REMOVE;
	dispatch_op(op, _("Remove layer"), FLIS_INV_ALL);
}

static void on_remove_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	GtkAlertDialog *ad = gtk_alert_dialog_new(
		_("Remove layer '%s'?"),
		lay->layer_name ? lay->layer_name : _("(unnamed)"));
	gtk_alert_dialog_set_detail(ad,
		_("This action cannot be undone.  The layer's pixel data and "
		  "mask will be discarded."));
	const char *buttons[] = { _("Remove"), _("Cancel"), NULL };
	gtk_alert_dialog_set_buttons(ad, buttons);
	gtk_alert_dialog_set_default_button(ad, 1);   /* Cancel is default */
	gtk_alert_dialog_set_cancel_button(ad, 1);
	gtk_alert_dialog_choose(ad,
		g_panel ? GTK_WINDOW(g_panel->window) : NULL,
		NULL, on_remove_confirm_resp, GINT_TO_POINTER(lay->item_id));
	g_object_unref(ad);
}
static void on_duplicate_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = lay->item_id;
	op->kind = OP_LAYER_DUPLICATE;
	dispatch_op(op, _("Duplicate layer"), FLIS_INV_ALL);
}
static void on_group_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	if (!is_current_image_flis()) return;
	/* No name prompt — the hook picks "Group N" with the smallest N
	 * not in use.  Users rename via the (future) group properties UI
	 * or by editing inline once group rows are made editable. */
	struct flis_addgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_addgroup_args_free;
	payload->name       = NULL;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_addgroup_hook;
	args->user               = payload;
	args->description        = g_strdup(_("Create group"));
	args->invalidate_flags   = FLIS_INV_STACK;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}
static void on_drag_toggled(GtkToggleButton *b, gpointer u) {
	(void)u;
	const gboolean want = gtk_toggle_button_get_active(b);
	if (want) {
		mouse_status = MOUSE_ACTION_FLIS_DRAG_LAYER;
		siril_log_message(_("FLIS: Canvas drag mode on — click and drag on the "
		                    "image to reposition the active layer\n"));
	} else {
		/* Defensive: if a drag was in flight when the user toggled off,
		 * end it cleanly without committing.  The current in-progress
		 * position stays — the user can toggle on again to continue or
		 * commit via flis_setposition. */
		gui.flis_layer_dragging = FALSE;
		/* Return to the default rectangular-select mouse mode, not to
		 * MOUSE_ACTION_NONE — leaving it at NONE leaves the user with
		 * no usable mouse action until they pick something else from
		 * the mouse-actions menu. */
		init_mouse();
		siril_log_message(_("FLIS: Canvas drag mode off\n"));
	}
}
/* Single-layer up/down — go through the slice-4 reorder hook with the
 * immediate visual neighbour as target.  Deliberately bypasses the
 * legacy flis_layer_move_up / _down primitives which jump past entire
 * groups (correct for "promote layer out of group" semantics but
 * wrong for "swap with immediate neighbour" which is what the panel
 * buttons should do — Photoshop / GIMP behaviour). */
static void move_relative(flis_layer_t *lay, gboolean direction_up) {
	if (!lay || !com.uniq) return;
	flis_layer_t *neighbour = NULL;
	for (GSList *l = com.uniq->layers; l; l = l->next) {
		flis_layer_t *cand = (flis_layer_t *)l->data;
		if (!cand || cand == lay) continue;
		if (direction_up) {
			if (cand->layer_order <= lay->layer_order) continue;
			/* lowest-order layer that's still above us */
			if (!neighbour || cand->layer_order < neighbour->layer_order)
				neighbour = cand;
		} else {
			if (cand->layer_order >= lay->layer_order) continue;
			/* highest-order layer that's still below us */
			if (!neighbour || cand->layer_order > neighbour->layer_order)
				neighbour = cand;
		}
	}
	if (!neighbour) return;  /* already at the top/bottom of the visual list */

	struct flis_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_reorder_args_free;
	payload->target_id   = neighbour->item_id;
	payload->place_above = direction_up;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_reorder_hook;
	args->user               = payload;
	args->description        = g_strdup(direction_up ? _("Move layer up")
	                                                  : _("Move layer down"));
	args->invalidate_flags   = FLIS_INV_STACK;
	args->invalidate_item_id = lay->item_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

/* Dispatch a "move group as a whole" op via flis_group_reorder_hook. */
static void move_group_relative(flis_group_t *grp, gboolean direction_up) {
	struct flis_group_reorder_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn   = flis_group_reorder_args_free;
	payload->direction_up = direction_up;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_group_reorder_hook;
	args->user               = payload;
	args->description        = g_strdup(direction_up ? _("Move group up")
	                                                  : _("Move group down"));
	args->invalidate_flags   = FLIS_INV_STACK;
	args->invalidate_item_id = grp->item_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

static void on_move_up_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (lay) { move_relative(lay, TRUE); return; }
	flis_group_t *grp = current_selected_group();
	if (grp) { move_group_relative(grp, TRUE); }
}
static void on_move_down_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (lay) { move_relative(lay, FALSE); return; }
	flis_group_t *grp = current_selected_group();
	if (grp) { move_group_relative(grp, FALSE); }
}
/* Mask sub-frame handlers (slice 2).
 *
 * status_btn click  → toggles lmask_active on the selected layer
 * toggle_btn click  → if layer has a mask: removes it via flis_clearmask_hook
 *                     if layer has no mask: opens FITS chooser → flis_setmask_hook
 * move_btn click    → stub (needs a layer chooser; deferred)
 *
 * Both worker dispatches mirror the flis_setmask / flis_clearmask
 * commands' args so panel and CLI paths produce identical state. */

static void dispatch_clearmask(flis_layer_t *lay) {
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_clearmask_hook;
	args->description        = g_strdup(_("Remove layer mask"));
	args->updates_lmask      = TRUE;
	args->invalidate_flags   = FLIS_INV_LAYER_PIXELS;
	args->invalidate_item_id = lay->item_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

static void on_mask_file_chosen(GObject *src, GAsyncResult *res, gpointer ud) {
	gint target_id = GPOINTER_TO_INT(ud);
	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	GError *err = NULL;
	GFile *file = gtk_file_dialog_open_finish(fd, res, &err);
	if (!file) {
		if (err && !g_error_matches(err, GTK_DIALOG_ERROR, GTK_DIALOG_ERROR_DISMISSED))
			siril_log_error(_("FLIS: Set Mask: %s\n"), err->message);
		g_clear_error(&err);
		return;
	}
	gchar *path = g_file_get_path(file);
	g_object_unref(file);
	if (!path) return;

	struct flis_setmask_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setmask_args_free;
	payload->filename   = path;     /* takes ownership */
	payload->bitpix     = com.pref.default_mask_bitpix > 0
	                       ? com.pref.default_mask_bitpix : 8;

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setmask_hook;
	args->user               = payload;
	args->description        = g_strdup(_("Set layer mask"));
	args->updates_lmask      = TRUE;
	args->invalidate_flags   = FLIS_INV_LAYER_PIXELS;
	args->invalidate_item_id = target_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

static void on_mask_status_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay || !lay->lmask) return;
	/* Dispatch through the worker like every other mutation: undo
	 * (lmask_active is part of the props snapshot), the stack writer
	 * lock, and invalidation all come from the shared path. */
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id  = lay->item_id;
	op->kind       = OP_SET_LMASK_ACTIVE;
	op->bool_v     = !lay->lmask_active;
	dispatch_op(op, _("Layer mask active"), FLIS_INV_LAYER_PIXELS);
}

static void on_mask_toggle_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	if (lay->lmask) {
		dispatch_clearmask(lay);
		return;
	}
	GtkFileDialog *fd = gtk_file_dialog_new();
	gtk_file_dialog_set_title(fd, _("Add layer mask from FITS file"));
	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_set_name(filter, _("FITS files"));
	gtk_file_filter_add_pattern(filter, "*.fit");
	gtk_file_filter_add_pattern(filter, "*.fits");
	gtk_file_filter_add_pattern(filter, "*.fts");
	gtk_file_filter_add_pattern(filter, "*.FIT");
	gtk_file_filter_add_pattern(filter, "*.FITS");
	GListStore *filters = g_list_store_new(GTK_TYPE_FILE_FILTER);
	g_list_store_append(filters, filter);
	gtk_file_dialog_set_default_filter(fd, filter);
	gtk_file_dialog_set_filters(fd, G_LIST_MODEL(filters));
	g_object_unref(filters);
	g_object_unref(filter);
	gtk_file_dialog_open(fd, g_panel ? GTK_WINDOW(g_panel->window) : NULL,
	                      NULL, on_mask_file_chosen,
	                      GINT_TO_POINTER(lay->item_id));
	g_object_unref(fd);
}

/* Mask Move dialog — same shape as Move-to-group: small modal with a
 * dropdown of layers OTHER than the source, plus OK/Cancel.  Dispatches
 * flis_movemask_hook on confirm.  flis_layer_move_lmask handles the
 * size-mismatch check itself. */
struct move_mask_ctx {
	GtkWidget *dialog;
	GtkWidget *dropdown;
	gint       from_id;
	GArray    *to_ids;       /* parallel to dropdown items: each is a layer item_id */
};

static void on_move_mask_ok(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct move_mask_ctx *ctx = ud;
	guint idx = gtk_drop_down_get_selected(GTK_DROP_DOWN(ctx->dropdown));
	if (idx >= ctx->to_ids->len) {
		gtk_window_destroy(GTK_WINDOW(ctx->dialog));
		g_array_unref(ctx->to_ids);
		g_free(ctx);
		return;
	}
	gint to_id = g_array_index(ctx->to_ids, gint, idx);

	struct flis_movemask_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn  = flis_movemask_args_free;
	payload->to_layer_id = to_id;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_movemask_hook;
	args->user               = payload;
	args->description        = g_strdup(_("Move layer mask"));
	args->updates_lmask      = TRUE;
	args->invalidate_flags   = FLIS_INV_LAYER_PIXELS;
	args->invalidate_item_id = ctx->from_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */

	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_array_unref(ctx->to_ids);
	g_free(ctx);
}

static void on_move_mask_cancel(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct move_mask_ctx *ctx = ud;
	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_array_unref(ctx->to_ids);
	g_free(ctx);
}

static void on_mask_move_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *src = current_selected_layer();
	if (!src) return;
	if (!src->lmask) {
		siril_log_message(_("FLIS: Move mask — selected layer has no mask\n"));
		return;
	}

	GtkWidget *dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Move layer mask"));
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(dialog),
	                              g_panel ? GTK_WINDOW(g_panel->window) : NULL);
	gtk_window_set_default_size(GTK_WINDOW(dialog), 320, 120);

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_widget_set_margin_start (box, 12);
	gtk_widget_set_margin_end   (box, 12);
	gtk_widget_set_margin_top   (box, 12);
	gtk_widget_set_margin_bottom(box, 12);
	gtk_window_set_child(GTK_WINDOW(dialog), box);

	gchar *prompt = g_strdup_printf(_("Move mask of '%s' to:"),
		src->layer_name ? src->layer_name : "");
	GtkWidget *label = gtk_label_new(prompt);
	g_free(prompt);
	gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
	gtk_box_append(GTK_BOX(box), label);

	GtkStringList *items = gtk_string_list_new(NULL);
	GArray *to_ids = g_array_new(FALSE, FALSE, sizeof(gint));
	for (GSList *l = com.uniq ? com.uniq->layers : NULL; l; l = l->next) {
		flis_layer_t *cand = (flis_layer_t *)l->data;
		if (!cand || cand == src) continue;     /* skip self */
		/* Size mismatch can be enforced at dispatch time but show all
		 * candidates so the user knows which to fix sizes for. */
		gtk_string_list_append(items, cand->layer_name ? cand->layer_name : _("(unnamed)"));
		g_array_append_val(to_ids, cand->item_id);
	}
	if (to_ids->len == 0) {
		GtkWidget *no_targets = gtk_label_new(_("(no other layers available)"));
		gtk_widget_add_css_class(no_targets, "dim-label");
		gtk_box_append(GTK_BOX(box), no_targets);
		g_array_unref(to_ids);
		g_object_unref(items);
		GtkWidget *close_btn = gtk_button_new_with_label(_("Close"));
		gtk_widget_set_halign(close_btn, GTK_ALIGN_END);
		g_signal_connect_swapped(close_btn, "clicked",
		                          G_CALLBACK(gtk_window_destroy), dialog);
		gtk_box_append(GTK_BOX(box), close_btn);
		gtk_window_present(GTK_WINDOW(dialog));
		return;
	}

	GtkWidget *dd = gtk_drop_down_new(G_LIST_MODEL(items), NULL);
	gtk_box_append(GTK_BOX(box), dd);

	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *cancel = gtk_button_new_with_label(_("Cancel"));
	GtkWidget *ok     = gtk_button_new_with_label(_("Move"));
	gtk_widget_add_css_class(ok, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), cancel);
	gtk_box_append(GTK_BOX(bbox), ok);
	gtk_box_append(GTK_BOX(box), bbox);

	struct move_mask_ctx *ctx = g_new0(struct move_mask_ctx, 1);
	ctx->dialog   = dialog;
	ctx->dropdown = dd;
	ctx->from_id  = src->item_id;
	ctx->to_ids   = to_ids;
	g_signal_connect(ok,     "clicked", G_CALLBACK(on_move_mask_ok),     ctx);
	g_signal_connect(cancel, "clicked", G_CALLBACK(on_move_mask_cancel), ctx);
	gtk_window_present(GTK_WINDOW(dialog));
}
static void on_mask_view_radio_toggled(GtkCheckButton *btn, gpointer u) {
	/* GINT_TO_POINTER(0) → Processing radio (view = 0)
	 * GINT_TO_POINTER(1) → Layer radio      (view = 1)
	 * Both radios in the same group emit "toggled" — handle the on
	 * transition only so we don't fire twice per click. */
	if (!siril_toggle_get_active(GTK_WIDGET(btn))) return;
	if (!is_current_image_flis() || !com.uniq) return;

	gint view = GPOINTER_TO_INT(u);
	if (com.uniq->flis_mask_view == view) return;
	com.uniq->flis_mask_view = view;

	/* Redraw both the mask tab and (if tint overlays are enabled) the
	 * image vports so the swap is immediately visible. */
	gui_iface.redraw_mask_idle(TRUE); // mask view swapped: tints are stale
	if (com.pref.gui.mask_tints_vports)
		gui_iface.redraw_image(REDRAW_ALL);
}

/* Context menu action handlers ---------------------------------- */

static int op_merge_down_hook(struct generic_layer_args *args) {
	flis_layer_t *lay = flis_layer_get_by_id(args->invalidate_item_id);
	if (!lay) return 1;
	return flis_merge_down_layer(lay);
}

/* Hook for "merge group": collapse every member into the bottom-most,
 * then remove the group so the surviving layer is ungrouped.  Members
 * are iterated top-down so each merge_down has a target underneath.
 * Group members must be at adjacent layer_order slots for the merge
 * sequence to be unambiguous; if a non-member layer is interleaved it
 * gets absorbed by mistake, so the user should reorder their stack
 * first.  Returns non-zero on any merge_down failure. */
static int op_merge_group_hook(struct generic_layer_args *args) {
	gint group_id = args->invalidate_item_id;
	flis_group_t *grp = flis_group_get_by_id(group_id);
	if (!grp) return 1;
	GSList *members = flis_group_get_layers(grp);  /* asc by layer_order */
	if (!members || !members->next) {
		/* 0 or 1 members — nothing to merge; just remove the group. */
		g_slist_free(members);
		flis_group_remove(grp);
		return 0;
	}
	/* Process from topmost member down to (and excluding) the bottom-
	 * most.  Reverse the list so we iterate top-first. */
	GSList *rev = g_slist_reverse(g_slist_copy(members));
	g_slist_free(members);
	int rv = 0;
	for (GSList *l = rev; l && l->next; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (flis_merge_down_layer(lay)) { rv = 1; break; }
	}
	g_slist_free(rev);
	if (!rv) flis_group_remove(grp);
	return rv;
}

static void on_ctx_merge_down(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	flis_group_t *grp = current_selected_group();
	if (grp) {
		GSList *members_tmp = flis_group_get_layers(grp);
		guint n_members = g_slist_length(members_tmp);
		g_slist_free(members_tmp);
		if (n_members < 2) {
			siril_log_message(
			    _("FLIS: Merge Down on group — group has fewer than 2 "
			      "members; nothing to merge\n"));
			return;
		}
		gchar *msg = g_strdup_printf(_(
		    "Merge all %u layers in group '%s' into the bottom-most "
		    "member and remove the group?\n\n"
		    "Undo for the merged layers is purged."),
		    n_members, grp->name ? grp->name : "?");
		gboolean ok = siril_confirm_dialog(_("Merge Group"), msg, _("Merge"));
		g_free(msg);
		if (!ok) return;
		struct generic_layer_args *args = calloc(1, sizeof(*args));
		args->layer = NULL;
		args->layer_hook = op_merge_group_hook;
		args->description = g_strdup(_("Merge group"));
		args->invalidate_flags   = FLIS_INV_ALL;
		args->invalidate_item_id = grp->item_id;
		if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
		return;
	}

	flis_layer_t *lay = current_selected_layer();
	if (!lay) {
		siril_log_message(_("FLIS: Merge Down — no layer or group selected\n"));
		return;
	}
	if (!siril_confirm_dialog(
	        _("Merge Down"),
	        _("Merge the active layer into the one beneath it?\n\n"
	          "Undo for both layers is purged after the merge."),
	        _("Merge"))) {
		return;
	}
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer = NULL;
	args->layer_hook = op_merge_down_hook;
	args->description = g_strdup(_("Merge Down"));
	args->invalidate_flags   = FLIS_INV_ALL;
	args->invalidate_item_id = lay->item_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

/* ---- Delete group dialog ---------------------------------------------
 *
 * Two-action choice: "Just remove group" (members become ungrouped,
 * layer pixel data preserved) vs "Delete group and members" (both
 * the group AND its member layers are destroyed).  Neither path is
 * undoable — group ops don't have an undo flavour yet. */

struct delete_group_ctx {
	GtkWidget *dialog;
	gint       group_id;
};

static int op_group_remove_hook(struct generic_layer_args *args) {
	flis_group_t *grp = flis_group_get_by_id(args->invalidate_item_id);
	return grp ? flis_group_remove(grp) : 1;
}

static int op_group_delete_with_layers_hook(struct generic_layer_args *args) {
	flis_group_t *grp = flis_group_get_by_id(args->invalidate_item_id);
	return grp ? flis_group_delete_with_layers(grp) : 1;
}

static void delete_group_dispatch(gint group_id,
                                  gboolean delete_members) {
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer              = NULL;
	args->layer_hook         = delete_members
	                           ? op_group_delete_with_layers_hook
	                           : op_group_remove_hook;
	args->user               = NULL;
	args->description        = g_strdup(delete_members
	                                    ? _("Delete group and members")
	                                    : _("Delete group (ungroup members)"));
	args->invalidate_flags   = FLIS_INV_ALL;
	args->invalidate_item_id = group_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

static void on_delete_group_just_remove(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct delete_group_ctx *ctx = ud;
	delete_group_dispatch(ctx->group_id, FALSE);
	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_free(ctx);
}

static void on_delete_group_with_members(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct delete_group_ctx *ctx = ud;
	flis_group_t *grp = flis_group_get_by_id(ctx->group_id);
	GSList *members = grp ? flis_group_get_layers(grp) : NULL;
	guint n_members = g_slist_length(members);
	g_slist_free(members);
	if (n_members > 0) {
		gchar *msg = g_strdup_printf(_(
		    "Permanently delete the group and all %u member layers?\n\n"
		    "Undo is not available for this operation."), n_members);
		gboolean ok = siril_confirm_dialog(_("Delete group and members"),
		                                    msg, _("Delete"));
		g_free(msg);
		if (!ok) return;  /* keep the chooser dialog open */
	}
	delete_group_dispatch(ctx->group_id, TRUE);
	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_free(ctx);
}

static void on_delete_group_cancel(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct delete_group_ctx *ctx = ud;
	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_free(ctx);
}

static void on_ctx_delete_group(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	flis_group_t *grp = current_selected_group();
	if (!grp) {
		siril_log_message(_("FLIS: Delete group — no group selected "
		                    "(select a group header row first)\n"));
		return;
	}

	struct delete_group_ctx *ctx = g_new0(struct delete_group_ctx, 1);
	ctx->group_id = grp->item_id;
	ctx->dialog   = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(ctx->dialog), _("Delete group"));
	gtk_window_set_modal(GTK_WINDOW(ctx->dialog), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(ctx->dialog),
	                             GTK_WINDOW(g_panel->window));
	gtk_window_set_default_size(GTK_WINDOW(ctx->dialog), 380, -1);

	GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_widget_set_margin_start (vbox, 12);
	gtk_widget_set_margin_end   (vbox, 12);
	gtk_widget_set_margin_top   (vbox, 12);
	gtk_widget_set_margin_bottom(vbox, 12);
	gtk_window_set_child(GTK_WINDOW(ctx->dialog), vbox);

	gchar *prompt = g_strdup_printf(
	    _("How should group '%s' be deleted?"),
	    grp->name ? grp->name : "?");
	GtkWidget *label = gtk_label_new(prompt);
	gtk_label_set_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
	gtk_box_append(GTK_BOX(vbox), label);
	g_free(prompt);

	GtkWidget *hint = gtk_label_new(_(
	    "• Just remove group — member layers are kept and become "
	    "ungrouped.\n"
	    "• Delete group and members — group AND every member layer "
	    "are permanently removed.\n\n"
	    "Neither action is undoable."));
	gtk_label_set_wrap(GTK_LABEL(hint), TRUE);
	gtk_label_set_xalign(GTK_LABEL(hint), 0.0f);
	gtk_widget_add_css_class(hint, "dim-label");
	gtk_box_append(GTK_BOX(vbox), hint);

	GtkWidget *btn_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(btn_row, GTK_ALIGN_END);
	GtkWidget *cancel_btn = gtk_button_new_with_label(_("Cancel"));
	GtkWidget *remove_btn = gtk_button_new_with_label(_("Just remove group"));
	GtkWidget *delete_btn = gtk_button_new_with_label(_("Delete group and members"));
	gtk_widget_add_css_class(delete_btn, "destructive-action");
	gtk_box_append(GTK_BOX(btn_row), cancel_btn);
	gtk_box_append(GTK_BOX(btn_row), remove_btn);
	gtk_box_append(GTK_BOX(btn_row), delete_btn);
	gtk_box_append(GTK_BOX(vbox), btn_row);

	g_signal_connect(cancel_btn, "clicked",
	                 G_CALLBACK(on_delete_group_cancel), ctx);
	g_signal_connect(remove_btn, "clicked",
	                 G_CALLBACK(on_delete_group_just_remove), ctx);
	g_signal_connect(delete_btn, "clicked",
	                 G_CALLBACK(on_delete_group_with_members), ctx);

	gtk_window_present(GTK_WINDOW(ctx->dialog));
}

static void on_ctx_flatten(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	if (!is_current_image_flis()) {
		siril_log_message(_("FLIS: Flatten — current image is not a FLIS\n"));
		return;
	}
	/* As destructive as Merge Down / Remove layer, which both confirm. */
	if (!siril_confirm_dialog(
	        _("Flatten Image"),
	        _("Flatten all visible layers into a single layer?\n\n"
	          "All other layers, masks and groups are removed, and the "
	          "undo history is purged."),
	        _("Flatten"))) {
		return;
	}
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer = NULL;
	args->layer_hook = flis_flatten_hook;
	args->description = g_strdup(_("Flatten Image"));
	args->invalidate_flags = FLIS_INV_ALL;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

/* Move-to-group dialog — a tiny modal window with a dropdown of
 * existing groups (plus "(none — remove from group)") and OK/Cancel.
 * On OK we dispatch flis_setgroup_hook with the chosen group_id. */
struct move_to_group_ctx {
	GtkWidget *dialog;
	GtkWidget *dropdown;
	gint       layer_id;
	GArray    *group_ids;  /* parallel to dropdown items: [0] = 0 (none), [n] = item_id */
};

static void on_move_to_group_ok(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct move_to_group_ctx *ctx = ud;
	guint idx = gtk_drop_down_get_selected(GTK_DROP_DOWN(ctx->dropdown));
	gint chosen_id = 0;
	if (idx < ctx->group_ids->len)
		chosen_id = g_array_index(ctx->group_ids, gint, idx);

	struct flis_setgroup_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_setgroup_args_free;
	payload->group_id   = chosen_id;
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_setgroup_hook;
	args->user               = payload;
	args->description        = g_strdup(_("Move layer to group"));
	args->invalidate_flags   = FLIS_INV_STACK;
	args->invalidate_item_id = ctx->layer_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */

	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_array_unref(ctx->group_ids);
	g_free(ctx);
}

static void on_move_to_group_cancel(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct move_to_group_ctx *ctx = ud;
	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_array_unref(ctx->group_ids);
	g_free(ctx);
}

static void on_ctx_move_to_group(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) {
		siril_log_message(_("FLIS: Move to group — no layer selected\n"));
		return;
	}

	GtkWidget *dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Move layer to group"));
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(dialog),
	                              g_panel ? GTK_WINDOW(g_panel->window) : NULL);
	gtk_window_set_default_size(GTK_WINDOW(dialog), 320, 120);

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_widget_set_margin_start (box, 12);
	gtk_widget_set_margin_end   (box, 12);
	gtk_widget_set_margin_top   (box, 12);
	gtk_widget_set_margin_bottom(box, 12);
	gtk_window_set_child(GTK_WINDOW(dialog), box);

	gchar *prompt = g_strdup_printf(_("Move layer '%s' to:"),
		lay->layer_name ? lay->layer_name : "");
	GtkWidget *label = gtk_label_new(prompt);
	g_free(prompt);
	gtk_label_set_xalign(GTK_LABEL(label), 0.0f);
	gtk_box_append(GTK_BOX(box), label);

	GtkStringList *items = gtk_string_list_new(NULL);
	gtk_string_list_append(items, _("(none — remove from group)"));
	GArray *group_ids = g_array_new(FALSE, FALSE, sizeof(gint));
	gint zero = 0;
	g_array_append_val(group_ids, zero);
	for (GSList *g = com.uniq ? com.uniq->groups : NULL; g; g = g->next) {
		flis_group_t *grp = (flis_group_t *)g->data;
		if (!grp) continue;
		gtk_string_list_append(items, grp->name ? grp->name : _("(unnamed)"));
		g_array_append_val(group_ids, grp->item_id);
	}
	GtkWidget *dd = gtk_drop_down_new(G_LIST_MODEL(items), NULL);
	/* Pre-select the layer's current group if any. */
	for (guint i = 0; i < group_ids->len; i++) {
		if (g_array_index(group_ids, gint, i) == lay->group_id) {
			gtk_drop_down_set_selected(GTK_DROP_DOWN(dd), i);
			break;
		}
	}
	gtk_box_append(GTK_BOX(box), dd);

	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *cancel = gtk_button_new_with_label(_("Cancel"));
	GtkWidget *ok     = gtk_button_new_with_label(_("OK"));
	gtk_widget_add_css_class(ok, "suggested-action");
	gtk_box_append(GTK_BOX(bbox), cancel);
	gtk_box_append(GTK_BOX(bbox), ok);
	gtk_box_append(GTK_BOX(box), bbox);

	struct move_to_group_ctx *ctx = g_new0(struct move_to_group_ctx, 1);
	ctx->dialog    = dialog;
	ctx->dropdown  = dd;
	ctx->layer_id  = lay->item_id;
	ctx->group_ids = group_ids;
	g_signal_connect(ok,     "clicked", G_CALLBACK(on_move_to_group_ok),     ctx);
	g_signal_connect(cancel, "clicked", G_CALLBACK(on_move_to_group_cancel), ctx);
	gtk_window_present(GTK_WINDOW(dialog));
}

/* Export-layer context menu action.  Opens a GtkFileDialog save-as
 * dialog; on confirm, dispatches flis_exportlayer_hook through the
 * worker so the panel and the flis_exportlayer command share the same
 * code path. */
static void on_export_file_chosen(GObject *src, GAsyncResult *res, gpointer ud) {
	gint target_id = GPOINTER_TO_INT(ud);
	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	GError *err = NULL;
	GFile *file = gtk_file_dialog_save_finish(fd, res, &err);
	if (!file) {
		if (err && !g_error_matches(err, GTK_DIALOG_ERROR, GTK_DIALOG_ERROR_DISMISSED))
			siril_log_error(_("FLIS: Export Layer: %s\n"), err->message);
		g_clear_error(&err);
		return;
	}
	gchar *path = g_file_get_path(file);
	g_object_unref(file);
	if (!path) return;

	struct flis_exportlayer_args *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = flis_exportlayer_args_free;
	payload->filename   = path;       /* takes ownership */
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer_hook         = flis_exportlayer_hook;
	args->user               = payload;
	args->description        = g_strdup(_("Export layer"));
	args->read_only          = TRUE;
	args->invalidate_item_id = target_id;
	if (!start_in_new_thread(generic_layer_worker, args))
		free_generic_layer_args(args);   /* busy refusal: nothing ran */
}

static void on_ctx_export_layer(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) {
		siril_log_message(_("FLIS: Export Layer — no layer selected\n"));
		return;
	}
	GtkFileDialog *fd = gtk_file_dialog_new();
	gtk_file_dialog_set_title(fd, _("Export layer as FITS file"));
	gchar *suggested = g_strdup_printf("%s.fit",
		lay->layer_name && *lay->layer_name ? lay->layer_name : "layer");
	gtk_file_dialog_set_initial_name(fd, suggested);
	g_free(suggested);
	gtk_file_dialog_save(fd, g_panel ? GTK_WINDOW(g_panel->window) : NULL,
	                      NULL, on_export_file_chosen,
	                      GINT_TO_POINTER(lay->item_id));
	g_object_unref(fd);
}

/* ---- Layers match (§5.7 — background neutralise) ---------------------
 *
 * Wraps flis_background_neutralise_layers in a generic_layer_worker
 * so progress / logging / completion idle all behave like other layer
 * ops.  When a group is currently selected in the panel we restrict
 * the operation to that group's layers (mirroring the GTK3 dialog's
 * implicit behaviour); otherwise every mono layer participates. */

struct ctx_layers_match_payload {
	destructor destroy_fn;
	GSList *subset;  /* NULL = all layers.  Owns the list shell only. */
};

static void ctx_layers_match_payload_free(void *p) {
	struct ctx_layers_match_payload *a = p;
	if (!a) return;
	if (a->subset) g_slist_free(a->subset);
	free(a);
}

static int op_layers_match_hook(struct generic_layer_args *args) {
	struct ctx_layers_match_payload *a =
	    (struct ctx_layers_match_payload *)args->user;
	return flis_background_neutralise_layers(a ? a->subset : NULL);
}

static void on_ctx_layers_match(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	if (!is_current_image_flis()) {
		siril_log_message(_("FLIS: Layers match — current image is not a FLIS\n"));
		return;
	}
	if (g_slist_length(com.uniq->layers) < 2) {
		siril_log_message(_("FLIS: Layers match — need at least two layers\n"));
		return;
	}
	if (!siril_confirm_dialog(
	        _("Layers match"),
	        _("Background-neutralise mono layers so their tinted "
	          "contributions sum to a neutral grey background. Layer "
	          "pixel data will be modified (undoable in one step)."),
	        _("Match"))) {
		return;
	}

	/* Build subset from the currently selected group, if any. */
	GSList *subset = NULL;
	flis_group_t *grp = current_selected_group();
	if (grp) subset = flis_group_get_layers(grp);

	GSList *snap_target = subset ? subset : com.uniq->layers;
	if (undo_save_flis_multi_layer(snap_target, _("Layers match"))) {
		siril_log_warning(_("Layers match: could not save undo state\n"));
	}

	struct ctx_layers_match_payload *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = ctx_layers_match_payload_free;
	payload->subset     = subset;  /* ownership transferred */

	struct generic_layer_args *gargs = calloc(1, sizeof(*gargs));
	gargs->layer             = NULL;       /* multi-layer op */
	gargs->layer_hook        = op_layers_match_hook;
	gargs->user              = payload;
	gargs->description       = g_strdup(_("Layers match"));
	gargs->invalidate_flags  = FLIS_INV_ALL;
	if (!start_in_new_thread(generic_layer_worker, gargs))
		free_generic_layer_args(gargs);  /* busy refusal: nothing ran */
}

/* ---- Register layers (§5.6) -----------------------------------------
 *
 * Tiny procedural dialog: reference layer dropdown + interpolation
 * dropdown + clamp toggle, then Cancel / Register.  On Register, the
 * primitive runs in a worker thread; a multi-layer undo snapshot is
 * taken first so the operation reverts in one step. */

struct ctx_register_args_dialog {
	GtkWidget *dialog;
	GtkWidget *method_dropdown;
	GtkWidget *ref_dropdown;
	GtkWidget *interp_dropdown;
	GtkWidget *clamp_toggle;
	GtkWidget *requirement_label;  /* dim label showing the current method's
	                                * selection requirement */
	GArray    *ref_ids;     /* parallel to ref_dropdown items: layer item_ids */
	gint       group_id;    /* 0 = register every layer; otherwise restrict
	                         * to this group's members.  Set from the panel
	                         * selection when the dialog opens; passed
	                         * through to op_register_layers_hook. */
};

struct ctx_register_payload {
	destructor   destroy_fn;
	gint         ref_item_id;       /* 0 = use current active */
	gint         group_id;          /* 0 = register all layers; otherwise
	                                 * restrict to this group's members */
	flis_reg_method_id   method_id;
	opencv_interpolation interp;
	gboolean     clamp;
};

static void ctx_register_payload_free(void *p) {
	free(p);
}

static int op_register_layers_hook(struct generic_layer_args *args) {
	struct ctx_register_payload *p = (struct ctx_register_payload *)args->user;
	flis_layer_t *ref = p->ref_item_id ? flis_layer_get_by_id(p->ref_item_id)
	                                   : flis_active_layer();
	selection_type sel;
	transformation_type tx;
	registration_function method = flis_register_resolve_method(
	    p->method_id, &sel, &tx);

	/* Group-scoped registration: pull just this group's members from
	 * com.uniq->layers so flis_register_layers builds its internal
	 * sequence over the subset only.  NULL → register all layers. */
	GSList *targets = NULL;
	gboolean owned = FALSE;
	if (p->group_id != 0) {
		flis_group_t *grp = flis_group_get_by_id(p->group_id);
		if (grp) {
			targets = flis_group_get_layers(grp);
			owned = TRUE;
		}
	}

	int rv = flis_register_layers(ref, targets, method, sel, tx,
	                              p->interp, p->clamp);
	if (owned) g_slist_free(targets);
	return rv;
}

static const opencv_interpolation REGISTER_INTERP_VALUES[] = {
	OPENCV_NEAREST, OPENCV_LINEAR, OPENCV_CUBIC, OPENCV_AREA, OPENCV_LANCZOS4
};
static const char *REGISTER_INTERP_NAMES[] = {
	N_("Nearest"), N_("Linear"), N_("Cubic"), N_("Area"), N_("Lanczos4")
};
#define REGISTER_INTERP_DEFAULT_IDX 4   /* Lanczos4 */

static const flis_reg_method_id REGISTER_METHOD_VALUES[] = {
	FLIS_REG_GLOBAL, FLIS_REG_2PASS, FLIS_REG_DFT, FLIS_REG_KOMBAT
};
static const char *REGISTER_METHOD_NAMES[] = {
	N_("Global star alignment (1-pass, recommended)"),
	N_("Global star alignment (2-pass)"),
	N_("Image pattern (DFT shift — needs a square selection)"),
	N_("KOMBAT pattern match (needs a selection)"),
};

/* Refresh the requirement-label below the method dropdown so the user
 * sees at a glance whether the chosen method needs a selection on the
 * image and whether the current com.selection satisfies it. */
static void update_register_requirement_hint(struct ctx_register_args_dialog *ctx) {
	guint idx = gtk_drop_down_get_selected(GTK_DROP_DOWN(ctx->method_dropdown));
	if (idx >= G_N_ELEMENTS(REGISTER_METHOD_VALUES)) idx = 0;
	selection_type sel;
	(void)flis_register_resolve_method(REGISTER_METHOD_VALUES[idx], &sel, NULL);
	const gboolean has_sel = (com.selection.w > 0 && com.selection.h > 0);
	const gboolean is_square = has_sel && (com.selection.w == com.selection.h);
	const char *txt = "";
	switch (sel) {
		case REQUIRES_NO_SELECTION:
			txt = _("No selection needed.");
			break;
		case REQUIRES_ANY_SELECTION:
			txt = has_sel
			    ? _("Selection ✓ — any size works.")
			    : _("Needs a selection — drag a rectangle on the image first.");
			break;
		case REQUIRES_SQUARED_SELECTION:
			txt = is_square
			    ? _("Selection ✓ — square, ready.")
			    : (has_sel
			        ? _("Selection must be SQUARE.  Adjust the rectangle to "
			            "equal width and height.")
			        : _("Needs a square selection — drag a square rectangle "
			            "on the image first."));
			break;
	}
	gtk_label_set_text(GTK_LABEL(ctx->requirement_label), txt);
}

static void on_register_method_changed(GObject *obj, GParamSpec *p,
                                       gpointer ud) {
	(void)obj; (void)p;
	update_register_requirement_hint((struct ctx_register_args_dialog *)ud);
}

static void on_register_dialog_ok(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct ctx_register_args_dialog *ctx = ud;

	guint ref_idx = gtk_drop_down_get_selected(GTK_DROP_DOWN(ctx->ref_dropdown));
	gint  ref_id = 0;
	if (ref_idx < ctx->ref_ids->len)
		ref_id = g_array_index(ctx->ref_ids, gint, ref_idx);

	guint method_idx = gtk_drop_down_get_selected(
	    GTK_DROP_DOWN(ctx->method_dropdown));
	if (method_idx >= G_N_ELEMENTS(REGISTER_METHOD_VALUES)) method_idx = 0;
	flis_reg_method_id method_id = REGISTER_METHOD_VALUES[method_idx];

	guint interp_idx = gtk_drop_down_get_selected(
	    GTK_DROP_DOWN(ctx->interp_dropdown));
	if (interp_idx >= G_N_ELEMENTS(REGISTER_INTERP_VALUES))
		interp_idx = REGISTER_INTERP_DEFAULT_IDX;

	gboolean clamp = siril_toggle_get_active(ctx->clamp_toggle);

	/* Up-front selection requirement check — duplicates the inside
	 * flis_register_layers check, but gives the user immediate feedback
	 * without closing the dialog or starting the worker. */
	selection_type sel_req;
	(void)flis_register_resolve_method(method_id, &sel_req, NULL);
	const gboolean has_sel    = (com.selection.w > 0 && com.selection.h > 0);
	const gboolean is_square  = has_sel && (com.selection.w == com.selection.h);
	if (sel_req == REQUIRES_ANY_SELECTION && !has_sel) {
		siril_log_error(_("Register layers: this method requires an "
		                  "image selection — drag a rectangle on the "
		                  "image first.\n"));
		return;
	}
	if (sel_req == REQUIRES_SQUARED_SELECTION && (!has_sel || !is_square)) {
		siril_log_error(_("Register layers: this method requires a SQUARE "
		                  "image selection — drag a square rectangle on "
		                  "the image first.\n"));
		return;
	}

	/* Snapshot every layer before the resampling pass. */
	if (undo_save_flis_multi_layer(com.uniq->layers, _("Register layers"))) {
		siril_log_warning(_("Register layers: could not save undo state\n"));
	}

	struct ctx_register_payload *payload = calloc(1, sizeof(*payload));
	payload->destroy_fn = ctx_register_payload_free;
	payload->ref_item_id = ref_id;
	payload->group_id    = ctx->group_id;
	payload->method_id   = method_id;
	payload->interp = REGISTER_INTERP_VALUES[interp_idx];
	payload->clamp = clamp;

	struct generic_layer_args *gargs = calloc(1, sizeof(*gargs));
	gargs->layer             = NULL;
	gargs->layer_hook        = op_register_layers_hook;
	gargs->user              = payload;
	gargs->description       = g_strdup(_("Register layers"));
	gargs->invalidate_flags  = FLIS_INV_ALL;
	if (!start_in_new_thread(generic_layer_worker, gargs))
		free_generic_layer_args(gargs);  /* busy refusal: nothing ran */

	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_array_unref(ctx->ref_ids);
	g_free(ctx);
}

static void on_register_dialog_cancel(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct ctx_register_args_dialog *ctx = ud;
	gtk_window_destroy(GTK_WINDOW(ctx->dialog));
	g_array_unref(ctx->ref_ids);
	g_free(ctx);
}

static void on_ctx_register_layers(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	if (!is_current_image_flis()) {
		siril_log_message(_("FLIS: Register layers — current image is not a FLIS\n"));
		return;
	}

	/* Scope: when a group is currently selected in the panel, register
	 * only the group's members.  Otherwise register every layer. */
	flis_group_t *sel_grp = current_selected_group();
	GSList *scope_layers = sel_grp ? flis_group_get_layers(sel_grp) : com.uniq->layers;
	const gboolean owned_scope = (sel_grp != NULL);
	const guint scope_n = g_slist_length(scope_layers);

	if (scope_n < 2) {
		if (owned_scope) g_slist_free(scope_layers);
		siril_log_message(_("FLIS: Register layers — need at least two "
		                    "layers in %s\n"),
		                  sel_grp ? _("the selected group") : _("the FLIS"));
		return;
	}

	struct ctx_register_args_dialog *ctx = g_new0(struct ctx_register_args_dialog, 1);
	ctx->ref_ids  = g_array_new(FALSE, FALSE, sizeof(gint));
	ctx->group_id = sel_grp ? sel_grp->item_id : 0;

	ctx->dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(ctx->dialog), _("Register layers"));
	gtk_window_set_modal(GTK_WINDOW(ctx->dialog), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(ctx->dialog),
	                             GTK_WINDOW(g_panel->window));
	gtk_window_set_default_size(GTK_WINDOW(ctx->dialog), 360, -1);

	GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_widget_set_margin_start (vbox, 12);
	gtk_widget_set_margin_end   (vbox, 12);
	gtk_widget_set_margin_top   (vbox, 12);
	gtk_widget_set_margin_bottom(vbox, 12);
	gtk_window_set_child(GTK_WINDOW(ctx->dialog), vbox);

	gtk_box_append(GTK_BOX(vbox),
	    gtk_label_new(_("Aligns the chosen layers to a reference, "
	                    "resampling pixel data and updating canvas offsets. "
	                    "Undoable in one step.")));

	/* Scope indicator — make it obvious which layers will be touched. */
	gchar *scope_txt = sel_grp
	    ? g_strdup_printf(_("Scope: %u layers in group '%s'"),
	                      scope_n, sel_grp->name ? sel_grp->name : "?")
	    : g_strdup_printf(_("Scope: all %u layers in the FLIS "
	                        "(select a group row in the panel before opening "
	                        "this dialog to register just that group)"),
	                      scope_n);
	GtkWidget *scope_label = gtk_label_new(scope_txt);
	gtk_label_set_xalign(GTK_LABEL(scope_label), 0.0f);
	gtk_label_set_wrap(GTK_LABEL(scope_label), TRUE);
	gtk_widget_add_css_class(scope_label, "dim-label");
	gtk_box_append(GTK_BOX(vbox), scope_label);
	g_free(scope_txt);

	/* Method dropdown — default to single-pass global star alignment so
	 * the dialog works on a freshly-loaded FLIS without any selection.
	 * DFT and KOMBAT need a selection; the requirement hint below the
	 * dropdown updates live as the user picks. */
	GtkWidget *method_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_box_append(GTK_BOX(method_row), gtk_label_new(_("Method:")));
	ctx->method_dropdown = gtk_drop_down_new(NULL, NULL);
	gtk_widget_set_hexpand(ctx->method_dropdown, TRUE);
	for (guint i = 0; i < G_N_ELEMENTS(REGISTER_METHOD_NAMES); i++)
		siril_drop_down_append_text(GTK_DROP_DOWN(ctx->method_dropdown),
		                            _(REGISTER_METHOD_NAMES[i]));
	gtk_drop_down_set_selected(GTK_DROP_DOWN(ctx->method_dropdown), 0);
	gtk_box_append(GTK_BOX(method_row), ctx->method_dropdown);
	gtk_box_append(GTK_BOX(vbox), method_row);

	ctx->requirement_label = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(ctx->requirement_label), 0.0f);
	gtk_label_set_wrap(GTK_LABEL(ctx->requirement_label), TRUE);
	gtk_widget_add_css_class(ctx->requirement_label, "dim-label");
	gtk_box_append(GTK_BOX(vbox), ctx->requirement_label);
	g_signal_connect(ctx->method_dropdown, "notify::selected",
	                 G_CALLBACK(on_register_method_changed), ctx);
	update_register_requirement_hint(ctx);

	/* Reference layer dropdown — only the layers in scope are listed
	 * (all layers, or the group's members when a group is selected). */
	GtkWidget *ref_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_box_append(GTK_BOX(ref_row), gtk_label_new(_("Reference:")));
	ctx->ref_dropdown = gtk_drop_down_new(NULL, NULL);
	gtk_widget_set_hexpand(ctx->ref_dropdown, TRUE);
	gint active_id = 0;
	{
		flis_layer_t *act = flis_active_layer();
		if (act) active_id = act->item_id;
	}
	guint active_idx = 0;
	guint idx = 0;
	for (GSList *l = scope_layers; l; l = l->next, idx++) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay) continue;
		const gchar *name = lay->layer_name ? lay->layer_name : "?";
		siril_drop_down_append_text(GTK_DROP_DOWN(ctx->ref_dropdown), name);
		g_array_append_val(ctx->ref_ids, lay->item_id);
		if (lay->item_id == active_id) active_idx = idx;
	}
	gtk_drop_down_set_selected(GTK_DROP_DOWN(ctx->ref_dropdown), active_idx);
	gtk_box_append(GTK_BOX(ref_row), ctx->ref_dropdown);
	gtk_box_append(GTK_BOX(vbox), ref_row);
	if (owned_scope) g_slist_free(scope_layers);
	scope_layers = NULL;

	/* Interpolation dropdown. */
	GtkWidget *interp_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_box_append(GTK_BOX(interp_row), gtk_label_new(_("Interpolation:")));
	ctx->interp_dropdown = gtk_drop_down_new(NULL, NULL);
	gtk_widget_set_hexpand(ctx->interp_dropdown, TRUE);
	for (guint i = 0; i < G_N_ELEMENTS(REGISTER_INTERP_NAMES); i++)
		siril_drop_down_append_text(GTK_DROP_DOWN(ctx->interp_dropdown),
		                            _(REGISTER_INTERP_NAMES[i]));
	gtk_drop_down_set_selected(GTK_DROP_DOWN(ctx->interp_dropdown),
	                           REGISTER_INTERP_DEFAULT_IDX);
	gtk_box_append(GTK_BOX(interp_row), ctx->interp_dropdown);
	gtk_box_append(GTK_BOX(vbox), interp_row);

	/* Clamp toggle (only meaningful for cubic / lanczos). */
	ctx->clamp_toggle = gtk_check_button_new_with_label(
	    _("Clamp cubic/Lanczos to suppress ringing"));
	siril_toggle_set_active(ctx->clamp_toggle, TRUE);
	gtk_box_append(GTK_BOX(vbox), ctx->clamp_toggle);

	/* Buttons row. */
	GtkWidget *btn_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(btn_row, GTK_ALIGN_END);
	GtkWidget *cancel_btn = gtk_button_new_with_label(_("Cancel"));
	GtkWidget *ok_btn     = gtk_button_new_with_label(_("Register"));
	gtk_widget_add_css_class(ok_btn, "suggested-action");
	gtk_box_append(GTK_BOX(btn_row), cancel_btn);
	gtk_box_append(GTK_BOX(btn_row), ok_btn);
	gtk_box_append(GTK_BOX(vbox), btn_row);

	g_signal_connect(cancel_btn, "clicked",
	                 G_CALLBACK(on_register_dialog_cancel), ctx);
	g_signal_connect(ok_btn, "clicked",
	                 G_CALLBACK(on_register_dialog_ok), ctx);

	gtk_window_present(GTK_WINDOW(ctx->dialog));
}

/* =========================================================================
 * Canvas properties dialog (§7 GUI)
 *
 * One modeless dialog holding every canvas-scoped operation: resize,
 * fit-to-layers, rotate, mirror.  The dialog stays open so the user can
 * chain several ops; the "Current" label refreshes after each one.
 * Each op saves a multi-layer-props undo entry before running so the
 * layer positions revert atomically (canvas dim change is not yet part
 * of the undo bundle — same caveat as the canvas commands per §7.10
 * deferred).
 * ========================================================================= */

/* Drag modes for the resize-preview drawingarea.  RESIZE_* corresponds
 * to a canvas-frame handle (compass directions).  MOVE drags the whole
 * canvas frame; ROTATE drags the dedicated rotation handle above the
 * (rotated) top edge.  NONE means no drag in progress / pointer landed
 * outside any actionable region. */
enum canvas_drag_mode {
	CDRAG_NONE = 0,
	CDRAG_MOVE,
	CDRAG_RESIZE_N, CDRAG_RESIZE_S, CDRAG_RESIZE_E, CDRAG_RESIZE_W,
	CDRAG_RESIZE_NE, CDRAG_RESIZE_NW, CDRAG_RESIZE_SE, CDRAG_RESIZE_SW,
	CDRAG_ROTATE,
};

struct canvas_dialog {
	GtkWidget *window;
	GtkWidget *current_label;
	GtkWidget *spin_w, *spin_h;
	GtkWidget *include_invisible;
	GtkWidget *spin_angle;

	/* Interactive resize preview ------------------------------------- */
	GtkWidget *preview_da;          /* GtkDrawingArea */

	/* Pending canvas geometry, in image coords.  The canvas top-left
	 * starts at (0,0) (same as the FLIS document); dragging the canvas
	 * frame shifts (cx,cy) so layers visually stay put.  Resize handles
	 * change (cw,ch) and possibly (cx,cy) when a top/left edge moves. */
	gint pending_cx, pending_cy;
	gint pending_cw, pending_ch;

	/* Pending rotation about to be applied, in degrees.  Visualised in
	 * the preview by rotating the canvas frame + layer outlines around
	 * the (pending) canvas centre; committed by the Rotate Apply
	 * button.  Sign matches flis_canvas_rotate (positive = CW in image
	 * y-down coords). */
	double pending_angle;

	/* Drag state ----------------------------------------------------- */
	enum canvas_drag_mode drag_mode;
	double drag_anchor_x, drag_anchor_y;     /* widget-space drag origin */
	gint   drag_start_cx, drag_start_cy;
	gint   drag_start_cw, drag_start_ch;
	/* Rotation drag: anchor angle (pointer→canvas-centre vector at the
	 * drag start, radians) and the pending_angle snapshot at that
	 * instant — we add (current_pointer_angle − anchor) to the start. */
	double drag_start_angle;          /* pending_angle at drag begin (deg) */
	double drag_anchor_rad;           /* pointer→centre at drag begin (rad) */

	gboolean syncing_spins;          /* re-entrancy guard for value-changed */
};

static void canvas_dialog_reset_pending(struct canvas_dialog *cd);

/* Singleton pointer to the currently-open canvas dialog (NULL when
 * the dialog isn't on screen).  Used by canvas_dialog_refresh_external
 * to repaint the preview when FLIS state changes underneath the dialog
 * — most importantly an undo / redo that reverts a canvas op. */
static struct canvas_dialog *g_canvas_dialog;

static void canvas_dialog_refresh_current(struct canvas_dialog *cd) {
	if (!cd || !cd->current_label) return;
	gchar *txt = g_strdup_printf(_("Current canvas: %u × %u"),
	    flis_canvas_rx(), flis_canvas_ry());
	gtk_label_set_text(GTK_LABEL(cd->current_label), txt);
	g_free(txt);
}

static void canvas_dialog_after_op(struct canvas_dialog *cd) {
	gui_iface.flis_invalidate_composite();
	gui_iface.flis_gui_update();
	notify_gfit_data_modified();
	/* notify_gfit_data_modified rebuilds the tile buffers but doesn't
	 * queue a paint of the image vports — the cvport widget keeps its
	 * last GtkSnapshot until something invalidates it.  Explicit
	 * redraw so the canvas op result is visible immediately rather
	 * than on the next incidental redraw (mouse move, focus change). */
	gui_iface.redraw_image(REDRAW_ALL);
	canvas_dialog_refresh_current(cd);
	/* Resync the spin defaults to the new canvas dims so a follow-up
	 * resize starts from the just-applied state.  Also reset the
	 * preview state (canvas back at (0,0), no pending drag). */
	canvas_dialog_reset_pending(cd);
	cd->syncing_spins = TRUE;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_w),
	                          (double)cd->pending_cw);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_h),
	                          (double)cd->pending_ch);
	if (cd->spin_angle)
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_angle), 0.0);
	cd->syncing_spins = FALSE;
	if (cd->preview_da)
		gtk_widget_queue_draw(cd->preview_da);
}

/* ---- Preview-area helpers --------------------------------------------
 *
 * Coordinate model: the document is laid out in image-space pixels.
 * Layers sit at (lay->position_x, lay->position_y) with size
 * (lay->fit->rx, lay->fit->ry).  The pending canvas sits at
 * (cd->pending_cx, cd->pending_cy) with size (cd->pending_cw,
 * cd->pending_ch).  The drawingarea fits the union of all those rects
 * (plus a margin) into its widget rect; widget→image and image→widget
 * conversions both use the same uniform scale + offset.
 * --------------------------------------------------------------------- */

#define CDRAG_HANDLE_PX  7.0   /* widget-space hit radius for edge / corner handles */
#define CROT_HANDLE_PX   8.0   /* widget-space hit radius for rotation handle */
#define CROT_OFFSET_PX  26.0   /* widget-space gap from top-centre of canvas to handle */

/* Image-space rotation of (*px,*py) around (cx,cy) by angle_rad.
 * Standard 2D rotation matrix — in image (y-down) coords, positive
 * angle = CW visually. */
static inline void rotate_pt(double cx, double cy, double angle_rad,
                              double *px, double *py) {
	double dx = *px - cx, dy = *py - cy;
	double c = cos(angle_rad), s = sin(angle_rad);
	*px = cx + dx * c - dy * s;
	*py = cy + dx * s + dy * c;
}

/* Centre of the pending canvas in image coords — also the pivot used
 * when applying the staged rotation to layer positions.  Conceptual
 * model: the canvas is the user's viewport, fixed and axis-aligned;
 * rotating the canvas by +θ is the same as rotating every layer's
 * position by -θ around this pivot.  Layer pixel content itself
 * never rotates (per the dialog's note). */
static inline void canvas_dialog_pivot(const struct canvas_dialog *cd,
                                        double *px, double *py) {
	*px = (double)cd->pending_cx + 0.5 * (double)cd->pending_cw;
	*py = (double)cd->pending_cy + 0.5 * (double)cd->pending_ch;
}

/* The 4 corners of a layer's outline in image coords, after the pending
 * canvas rotation has been applied.  Canvas-rotates-+θ is realised by
 * rotating every layer (position + pixel content) by -θ around the
 * canvas centre, so the corners of the original axis-aligned rect get
 * rotated together — producing a tilted quadrilateral the dialog can
 * draw as the layer's preview outline. */
static void canvas_dialog_rotated_layer_corners(const struct canvas_dialog *cd,
                                                 const flis_layer_t *lay,
                                                 double xs[4], double ys[4]) {
	double pcx, pcy;
	canvas_dialog_pivot(cd, &pcx, &pcy);
	double rad = -cd->pending_angle * M_PI / 180.0;
	double lx = lay->position_x, ly = lay->position_y;
	double lw = lay->fit->rx,    lh = lay->fit->ry;
	xs[0] = lx;        ys[0] = ly;
	xs[1] = lx + lw;   ys[1] = ly;
	xs[2] = lx + lw;   ys[2] = ly + lh;
	xs[3] = lx;        ys[3] = ly + lh;
	for (int i = 0; i < 4; i++)
		rotate_pt(pcx, pcy, rad, &xs[i], &ys[i]);
}

/* Union of the (axis-aligned) canvas frame and every layer's rotated
 * polygon corners.  Returns FALSE if the scene is empty. */
static gboolean canvas_dialog_scene_bbox(const struct canvas_dialog *cd,
                                          double *x0, double *y0,
                                          double *x1, double *y1) {
	if (!is_current_image_flis()) return FALSE;
	double bx0 = (double)cd->pending_cx;
	double by0 = (double)cd->pending_cy;
	double bx1 = bx0 + (double)cd->pending_cw;
	double by1 = by0 + (double)cd->pending_ch;
	for (GSList *l = com.uniq ? com.uniq->layers : NULL; l; l = l->next) {
		const flis_layer_t *lay = (const flis_layer_t *)l->data;
		if (!lay || !lay->fit) continue;
		double xs[4], ys[4];
		canvas_dialog_rotated_layer_corners(cd, lay, xs, ys);
		for (int i = 0; i < 4; i++) {
			if (xs[i] < bx0) bx0 = xs[i];
			if (ys[i] < by0) by0 = ys[i];
			if (xs[i] > bx1) bx1 = xs[i];
			if (ys[i] > by1) by1 = ys[i];
		}
	}
	*x0 = bx0; *y0 = by0; *x1 = bx1; *y1 = by1;
	return (bx1 > bx0) && (by1 > by0);
}

/* Uniform image→widget transform fitting the scene bbox into the widget
 * with a small margin.  Returned via out params; (sx,sy) is the offset
 * of image-origin in widget coords, scale is widget-pixels per image-px.
 *
 * Computes scale from the more constrained axis so aspect is preserved
 * and centres the scene in the unused axis. */
static void canvas_dialog_compute_transform(const struct canvas_dialog *cd,
                                             int da_w, int da_h,
                                             double *out_scale,
                                             double *out_ox, double *out_oy) {
	/* Margin sized to keep the rotation handle (tethered CROT_OFFSET_PX
	 * past the canvas top edge) inside the widget at angle=0; rotated
	 * angles can only be tighter because the rotated bbox grows. */
	const double margin = CROT_OFFSET_PX + 8.0;
	double sx0, sy0, sx1, sy1;
	if (!canvas_dialog_scene_bbox(cd, &sx0, &sy0, &sx1, &sy1) ||
	    da_w <= 2 * margin || da_h <= 2 * margin) {
		*out_scale = 1.0; *out_ox = 0.0; *out_oy = 0.0;
		return;
	}
	double scene_w = sx1 - sx0;
	double scene_h = sy1 - sy0;
	double avail_w = (double)da_w - 2.0 * margin;
	double avail_h = (double)da_h - 2.0 * margin;
	double s = MIN(avail_w / scene_w, avail_h / scene_h);
	if (s <= 0.0) s = 1.0;
	/* Centre. */
	double used_w = scene_w * s;
	double used_h = scene_h * s;
	*out_scale = s;
	*out_ox = margin + 0.5 * (avail_w - used_w) - sx0 * s;
	*out_oy = margin + 0.5 * (avail_h - used_h) - sy0 * s;
}

/* Convert a widget-space delta to an integer image-space delta given
 * the current transform.  Rounds to nearest. */
static gint widget_delta_to_image(double dpx, double scale) {
	if (scale <= 0.0) return 0;
	double v = dpx / scale;
	return (gint)(v >= 0.0 ? v + 0.5 : v - 0.5);
}

/* Position of the rotation handle in widget coords.  It hovers above the
 * (rotated) top-centre of the canvas frame, anchored to the canvas so
 * the user grabs the same point as the frame turns. */
static void canvas_dialog_rot_handle_pos(const struct canvas_dialog *cd,
                                          double scale, double ox, double oy,
                                          double *hx, double *hy) {
	double pcx_img, pcy_img;
	canvas_dialog_pivot(cd, &pcx_img, &pcy_img);
	double pcx_w = ox + pcx_img * scale;
	double pcy_w = oy + pcy_img * scale;
	double r_w   = 0.5 * (double)cd->pending_ch * scale + CROT_OFFSET_PX;
	double rad   = cd->pending_angle * M_PI / 180.0;
	/* "Up" in image-y-down coords is -y; rotating the unit -y vector
	 * by rad gives (sin rad, -cos rad) — that's the direction the
	 * handle sits from the canvas centre. */
	*hx = pcx_w + sin(rad) * r_w;
	*hy = pcy_w - cos(rad) * r_w;
}

/* Hit-test the pending canvas frame.  (wx,wy) are widget-space.  The
 * canvas frame is axis-aligned in the preview (it represents the
 * viewport, which never tilts), so the tests are straight axis-aligned
 * bounds.  The rotation handle takes priority because it can sit close
 * to the top edge at small canvases. */
static enum canvas_drag_mode
canvas_dialog_hit_test(const struct canvas_dialog *cd,
                       double wx, double wy,
                       double scale, double ox, double oy) {
	double hx, hy;
	canvas_dialog_rot_handle_pos(cd, scale, ox, oy, &hx, &hy);
	if (hypot(wx - hx, wy - hy) <= CROT_HANDLE_PX) return CDRAG_ROTATE;

	double cx0 = ox + (double)cd->pending_cx * scale;
	double cy0 = oy + (double)cd->pending_cy * scale;
	double cx1 = cx0 + (double)cd->pending_cw  * scale;
	double cy1 = cy0 + (double)cd->pending_ch  * scale;
	const double r = CDRAG_HANDLE_PX;

	gboolean near_l = fabs(wx - cx0) <= r && wy >= cy0 - r && wy <= cy1 + r;
	gboolean near_r = fabs(wx - cx1) <= r && wy >= cy0 - r && wy <= cy1 + r;
	gboolean near_t = fabs(wy - cy0) <= r && wx >= cx0 - r && wx <= cx1 + r;
	gboolean near_b = fabs(wy - cy1) <= r && wx >= cx0 - r && wx <= cx1 + r;
	if (near_t && near_l) return CDRAG_RESIZE_NW;
	if (near_t && near_r) return CDRAG_RESIZE_NE;
	if (near_b && near_l) return CDRAG_RESIZE_SW;
	if (near_b && near_r) return CDRAG_RESIZE_SE;
	if (near_l) return CDRAG_RESIZE_W;
	if (near_r) return CDRAG_RESIZE_E;
	if (near_t) return CDRAG_RESIZE_N;
	if (near_b) return CDRAG_RESIZE_S;
	if (wx >= cx0 && wx <= cx1 && wy >= cy0 && wy <= cy1) return CDRAG_MOVE;
	return CDRAG_NONE;
}

/* Draw callback.  The canvas is the user's viewport — always axis-
 * aligned — so the yellow frame, its resize handles, and the layer
 * outlines (which never tilt because layer pixels stay axis-aligned)
 * are all drawn as axis-aligned rectangles.  The rotation only moves
 * each layer's centre by −pending_angle around the canvas pivot,
 * showing the user where the layers will end up relative to the
 * fixed canvas frame.  The rotation handle is a blue dot tethered
 * from the canvas centre out into the rotation direction. */
static void on_canvas_preview_draw(GtkDrawingArea *area, cairo_t *cr,
                                    int width, int height, gpointer user_data) {
	(void)area;
	const struct canvas_dialog *cd = user_data;
	double scale, ox, oy;
	canvas_dialog_compute_transform(cd, width, height, &scale, &ox, &oy);
	double pcx_img, pcy_img;
	canvas_dialog_pivot(cd, &pcx_img, &pcy_img);

	/* Background. */
	cairo_set_source_rgb(cr, 0.12, 0.12, 0.13);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	/* Layers in light grey: each drawn as the 4-corner polygon you get
	 * from rotating the layer's axis-aligned rect by -pending_angle
	 * around the canvas centre.  This mirrors what Apply will produce
	 * (pixel content rotated, layer dims = rotated bbox, centre moved
	 * to the rotated position). */
	for (GSList *l = com.uniq ? com.uniq->layers : NULL; l; l = l->next) {
		const flis_layer_t *lay = (const flis_layer_t *)l->data;
		if (!lay || !lay->fit) continue;
		double xs[4], ys[4];
		canvas_dialog_rotated_layer_corners(cd, lay, xs, ys);
		cairo_move_to(cr, ox + xs[0] * scale, oy + ys[0] * scale);
		for (int i = 1; i < 4; i++)
			cairo_line_to(cr, ox + xs[i] * scale, oy + ys[i] * scale);
		cairo_close_path(cr);
		cairo_set_source_rgba(cr, 0.78, 0.78, 0.78, 0.18);
		cairo_fill_preserve(cr);
		cairo_set_source_rgba(cr, 0.85, 0.85, 0.85, 0.85);
		cairo_set_line_width(cr, 1.0);
		cairo_stroke(cr);
	}

	/* Canvas frame in yellow, axis-aligned. */
	double cx = ox + (double)cd->pending_cx * scale;
	double cy = oy + (double)cd->pending_cy * scale;
	double cw = (double)cd->pending_cw * scale;
	double ch = (double)cd->pending_ch * scale;
	cairo_set_source_rgba(cr, 1.0, 0.85, 0.10, 1.0);
	cairo_set_line_width(cr, 2.0);
	cairo_rectangle(cr, cx, cy, cw, ch);
	cairo_stroke(cr);

	/* Resize handles: 8 small filled squares on the axis-aligned frame
	 * (corners + edge midpoints). */
	const double hs = 4.0;
	const double mids_x[] = { cx, cx + cw * 0.5, cx + cw };
	const double mids_y[] = { cy, cy + ch * 0.5, cy + ch };
	cairo_set_source_rgb(cr, 1.0, 0.85, 0.10);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == 1 && j == 1) continue;  /* skip centre */
			cairo_rectangle(cr, mids_x[i] - hs, mids_y[j] - hs, 2 * hs, 2 * hs);
			cairo_fill(cr);
		}
	}

	/* Rotation handle: tether from canvas centre out to the angle
	 * indicator.  At pending_angle = 0 the handle sits straight above
	 * the top edge; positive angles swing it CW around the centre. */
	double pcx_w = ox + pcx_img * scale;
	double pcy_w = oy + pcy_img * scale;
	double rh_x, rh_y;
	canvas_dialog_rot_handle_pos(cd, scale, ox, oy, &rh_x, &rh_y);
	cairo_set_source_rgba(cr, 0.55, 0.85, 1.0, 0.85);
	cairo_set_line_width(cr, 1.5);
	cairo_move_to(cr, pcx_w, pcy_w);
	cairo_line_to(cr, rh_x, rh_y);
	cairo_stroke(cr);
	cairo_arc(cr, rh_x, rh_y, 5.0, 0.0, 2.0 * M_PI);
	cairo_fill(cr);
}

/* Push the pending (cw,ch) into the width/height spinbuttons without
 * triggering their value-changed handlers. */
static void canvas_dialog_sync_spins_from_pending(struct canvas_dialog *cd) {
	if (!cd->spin_w || !cd->spin_h) return;
	cd->syncing_spins = TRUE;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_w), (double)cd->pending_cw);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_h), (double)cd->pending_ch);
	cd->syncing_spins = FALSE;
}

/* Same, for the angle spinner. */
static void canvas_dialog_sync_angle_from_pending(struct canvas_dialog *cd) {
	if (!cd->spin_angle) return;
	cd->syncing_spins = TRUE;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_angle), cd->pending_angle);
	cd->syncing_spins = FALSE;
}

/* Spinbutton value-changed handler.  Pulls the new width/height into the
 * pending state and redraws the preview.  Width/height adjustments via
 * the spinners hold the canvas's top-left fixed (pending_cx/cy unchanged) —
 * this matches the "drag the SE handle" convention. */
static void on_canvas_size_spin_changed(GtkSpinButton *sb, gpointer ud) {
	struct canvas_dialog *cd = ud;
	if (cd->syncing_spins) return;
	cd->pending_cw = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(cd->spin_w));
	cd->pending_ch = (gint)gtk_spin_button_get_value(GTK_SPIN_BUTTON(cd->spin_h));
	if (cd->pending_cw < 1) cd->pending_cw = 1;
	if (cd->pending_ch < 1) cd->pending_ch = 1;
	(void)sb;
	if (cd->preview_da) gtk_widget_queue_draw(cd->preview_da);
}

/* Gesture drag callbacks ------------------------------------------------ */

static void on_canvas_preview_drag_begin(GtkGestureDrag *g,
                                          double start_x, double start_y,
                                          gpointer ud) {
	struct canvas_dialog *cd = ud;
	int da_w = gtk_widget_get_width (cd->preview_da);
	int da_h = gtk_widget_get_height(cd->preview_da);
	double scale, ox, oy;
	canvas_dialog_compute_transform(cd, da_w, da_h, &scale, &ox, &oy);
	cd->drag_mode = canvas_dialog_hit_test(cd, start_x, start_y, scale, ox, oy);
	cd->drag_anchor_x = start_x;
	cd->drag_anchor_y = start_y;
	cd->drag_start_cx = cd->pending_cx;
	cd->drag_start_cy = cd->pending_cy;
	cd->drag_start_cw = cd->pending_cw;
	cd->drag_start_ch = cd->pending_ch;
	if (cd->drag_mode == CDRAG_ROTATE) {
		/* Snapshot the pointer→pivot angle so the user's full drag arc
		 * maps 1:1 onto the rotation delta, regardless of where they
		 * grabbed the handle. */
		double pcx_img, pcy_img;
		canvas_dialog_pivot(cd, &pcx_img, &pcy_img);
		double pcx_w = ox + pcx_img * scale;
		double pcy_w = oy + pcy_img * scale;
		cd->drag_anchor_rad = atan2(start_y - pcy_w, start_x - pcx_w);
		cd->drag_start_angle = cd->pending_angle;
	}
	if (cd->drag_mode == CDRAG_NONE)
		gtk_gesture_set_state(GTK_GESTURE(g), GTK_EVENT_SEQUENCE_DENIED);
}

static void on_canvas_preview_drag_update(GtkGestureDrag *g,
                                           double offx, double offy,
                                           gpointer ud) {
	(void)g;
	struct canvas_dialog *cd = ud;
	if (cd->drag_mode == CDRAG_NONE) return;
	int da_w = gtk_widget_get_width (cd->preview_da);
	int da_h = gtk_widget_get_height(cd->preview_da);
	double scale, ox, oy;
	canvas_dialog_compute_transform(cd, da_w, da_h, &scale, &ox, &oy);

	if (cd->drag_mode == CDRAG_ROTATE) {
		double pcx_img, pcy_img;
		canvas_dialog_pivot(cd, &pcx_img, &pcy_img);
		double pcx_w = ox + pcx_img * scale;
		double pcy_w = oy + pcy_img * scale;
		double now_x = cd->drag_anchor_x + offx;
		double now_y = cd->drag_anchor_y + offy;
		double now_rad = atan2(now_y - pcy_w, now_x - pcx_w);
		double delta_deg = (now_rad - cd->drag_anchor_rad) * 180.0 / M_PI;
		cd->pending_angle = cd->drag_start_angle + delta_deg;
		/* Normalise to (-180, 180] so the spinner doesn't drift across
		 * many turns; the rotation op only cares about the mod-360
		 * residue anyway. */
		while (cd->pending_angle >  180.0) cd->pending_angle -= 360.0;
		while (cd->pending_angle <= -180.0) cd->pending_angle += 360.0;
		canvas_dialog_sync_angle_from_pending(cd);
		gtk_widget_queue_draw(cd->preview_da);
		return;
	}

	/* Move / resize: canvas frame is axis-aligned in the preview, so
	 * widget-space deltas map straight to image-space deltas. */
	gint dxi = widget_delta_to_image(offx, scale);
	gint dyi = widget_delta_to_image(offy, scale);

	gint nx = cd->drag_start_cx;
	gint ny = cd->drag_start_cy;
	gint nw = cd->drag_start_cw;
	gint nh = cd->drag_start_ch;

	switch (cd->drag_mode) {
	case CDRAG_MOVE:
		nx += dxi; ny += dyi;
		break;
	/* Edge / corner drags: a north or west edge shifts the origin AND
	 * adjusts the dimension by the opposite amount; south / east just
	 * change the dimension. */
	case CDRAG_RESIZE_N:  ny += dyi; nh -= dyi; break;
	case CDRAG_RESIZE_S:  nh += dyi;            break;
	case CDRAG_RESIZE_W:  nx += dxi; nw -= dxi; break;
	case CDRAG_RESIZE_E:  nw += dxi;            break;
	case CDRAG_RESIZE_NW: nx += dxi; nw -= dxi; ny += dyi; nh -= dyi; break;
	case CDRAG_RESIZE_NE: nw += dxi;            ny += dyi; nh -= dyi; break;
	case CDRAG_RESIZE_SW: nx += dxi; nw -= dxi; nh += dyi;            break;
	case CDRAG_RESIZE_SE: nw += dxi;            nh += dyi;            break;
	default: return;
	}
	/* Clamp dimensions; an edge drag that would invert the rect pins
	 * that edge at the opposite side instead of producing a 0/negative
	 * extent (which would also confuse subsequent hit-tests). */
	if (nw < 1) {
		if (cd->drag_mode == CDRAG_RESIZE_W || cd->drag_mode == CDRAG_RESIZE_NW ||
		    cd->drag_mode == CDRAG_RESIZE_SW)
			nx = cd->drag_start_cx + cd->drag_start_cw - 1;
		nw = 1;
	}
	if (nh < 1) {
		if (cd->drag_mode == CDRAG_RESIZE_N || cd->drag_mode == CDRAG_RESIZE_NW ||
		    cd->drag_mode == CDRAG_RESIZE_NE)
			ny = cd->drag_start_cy + cd->drag_start_ch - 1;
		nh = 1;
	}
	cd->pending_cx = nx;
	cd->pending_cy = ny;
	cd->pending_cw = nw;
	cd->pending_ch = nh;
	canvas_dialog_sync_spins_from_pending(cd);
	gtk_widget_queue_draw(cd->preview_da);
}

static void on_canvas_preview_drag_end(GtkGestureDrag *g,
                                        double offx, double offy,
                                        gpointer ud) {
	(void)g; (void)offx; (void)offy;
	struct canvas_dialog *cd = ud;
	cd->drag_mode = CDRAG_NONE;
}

/* Reset the pending resize state to match the current canvas.  Called
 * on dialog open and after every successful canvas op so the preview
 * starts from a clean slate (canvas at (0,0), layers unmoved). */
static void canvas_dialog_reset_pending(struct canvas_dialog *cd) {
	cd->pending_cx = 0;
	cd->pending_cy = 0;
	cd->pending_cw = (gint)flis_canvas_rx();
	cd->pending_ch = (gint)flis_canvas_ry();
	cd->pending_angle = 0.0;
	cd->drag_mode = CDRAG_NONE;
}

/* For every layer: resample its pixel data by -pending_angle around the
 * layer's own centre, then reposition the (now larger) layer so its
 * centre sits at the original-centre rotated around the canvas centre
 * by -pending_angle.  Net effect: each layer is rigidly rotated by
 * -pending_angle around the canvas centre — pixel content and position
 * stay in sync, so layers that previously tiled still tile after the
 * rotation.
 *
 * Sign convention: pending_angle (canvas rotate, +CW visually as the
 * handle traces) means layers must move and re-render at -pending_angle
 * visually (i.e. CCW for +pending_angle).  In both Siril's
 * verbose_rotate_image API (+ = CW visually, matching the menu's
 * "Rotate 90° clockwise" wired to angle = 90) AND our standard y-down
 * rotation matrix (+ = CW visually), the value we want to pass for the
 * direction "CCW visually" is therefore -pending_angle.  Both pixel
 * and position rotation use the same -pending_angle, keeping them in
 * lock-step.
 *
 * Returns 0 on success, non-zero (with an error logged) on any pixel-
 * rotation failure (e.g. OOM); the caller's undo step lets the user
 * revert. */
static int canvas_dialog_rotate_layers_pixels(double pivot_x, double pivot_y,
                                                double pending_angle_deg) {
	if (fabs(pending_angle_deg) < 1e-9 || !com.uniq || !com.uniq->layers)
		return 0;
	double rot_angle_deg = -pending_angle_deg;
	double rot_rad = rot_angle_deg * M_PI / 180.0;
	double c = cos(rot_rad), s = sin(rot_rad);
	for (GSList *l = com.uniq->layers; l; l = l->next) {
		flis_layer_t *lay = (flis_layer_t *)l->data;
		if (!lay || !lay->fit) continue;
		double old_lcx = (double)lay->position_x + 0.5 * (double)lay->fit->rx;
		double old_lcy = (double)lay->position_y + 0.5 * (double)lay->fit->ry;
		rectangle area = { 0, 0, (int)lay->fit->rx, (int)lay->fit->ry };
		/* uncropped (= 0) so the layer expands to the rotated bbox;
		 * clamp on so float layers don't wrap on overshoot. */
		if (verbose_rotate_image(lay->fit, area, rot_angle_deg,
		                         OPENCV_LANCZOS4, 0, TRUE) != 0) {
			siril_log_error(_("Canvas rotate: failed to rotate layer pixels.\n"));
			return 1;
		}
		double dx = old_lcx - pivot_x;
		double dy = old_lcy - pivot_y;
		double new_lcx = pivot_x + dx * c - dy * s;
		double new_lcy = pivot_y + dx * s + dy * c;
		lay->position_x = (gint)lround(new_lcx - 0.5 * (double)lay->fit->rx);
		lay->position_y = (gint)lround(new_lcy - 0.5 * (double)lay->fit->ry);
		flis_layer_touch_modified(lay);
	}
	gui_iface.flis_invalidate_composite();
	return 0;
}

/* Combined apply: resize (move + dim change) then resample each layer's
 * pixel data + reposition it around the NEW canvas centre, all under
 * one full-pixel undo entry so the user can undo back to the pre-Apply
 * state in a single step.
 *
 * The heavy work — full-pixel undo snapshots of every layer plus a
 * Lanczos4 resample per layer for rotation — runs on the processing
 * worker via generic_layer_worker: the GTK main thread stays live, the
 * worker's busy-guard refuses the op while another job runs, and the
 * hook executes under the FLIS stack writer lock like every other
 * layer mutation.  The default end idle plus the coalesced
 * flis_gui_update → canvas_dialog_refresh_external() reproduce what
 * canvas_dialog_after_op did for the synchronous path. */
struct canvas_apply_args {
	destructor destroy_fn;
	gint   cw, ch, cx, cy;
	double angle;
};

static void canvas_apply_args_free(gpointer p) {
	g_free(p);
}

static int canvas_apply_hook(struct generic_layer_args *args) {
	struct canvas_apply_args *a = (struct canvas_apply_args *)args->user;
	if (!a || !is_current_image_flis()) return 1;
	gboolean want_resize = (a->cx != 0 || a->cy != 0 ||
	                         a->cw != (gint)flis_canvas_rx() ||
	                         a->ch != (gint)flis_canvas_ry());
	gboolean want_rotate = (fabs(a->angle) >= 1e-9);
	if (!want_resize && !want_rotate) return 0;

	/* Full pixel undo when rotation is involved (each layer's pixel
	 * data changes); props-only is enough for resize alone. */
	if (want_rotate)
		undo_save_flis_multi_layer(com.uniq->layers, _("Adjust canvas"));
	else
		undo_save_flis_multi_layer_props(com.uniq->layers, _("Adjust canvas"));

	if (want_resize) {
		/* Canvas-moves-by-(cx,cy) ⇒ layer-shifts-by-the-opposite, so
		 * layers visually stay put.  flis_canvas_resize adds its dx/dy
		 * to each position. */
		if (flis_canvas_resize((guint)a->cw, (guint)a->ch,
		                       -a->cx, -a->cy) != 0)
			return 1;
	}
	if (want_rotate) {
		double pivot_x = 0.5 * (double)flis_canvas_rx();
		double pivot_y = 0.5 * (double)flis_canvas_ry();
		if (canvas_dialog_rotate_layers_pixels(pivot_x, pivot_y,
		                                        a->angle) != 0) {
			/* A mid-loop failure leaves some layers already rotated —
			 * refresh the display anyway so it doesn't show stale
			 * pixels (the end idle skips invalidation on failure). */
			gui_iface.flis_invalidate_composite();
			gui_iface.flis_gui_update();
			return 1;
		}
	}
	return 0;
}

static void on_canvas_apply_click(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct canvas_dialog *cd = ud;
	if (!is_current_image_flis()) return;
	if (cd->pending_cw <= 0 || cd->pending_ch <= 0) {
		siril_log_error(_("Canvas dims must be positive (got %dx%d)\n"),
		                cd->pending_cw, cd->pending_ch);
		return;
	}
	gboolean want_resize = (cd->pending_cx != 0 || cd->pending_cy != 0 ||
	                         cd->pending_cw != (gint)flis_canvas_rx() ||
	                         cd->pending_ch != (gint)flis_canvas_ry());
	gboolean want_rotate = (fabs(cd->pending_angle) >= 1e-9);
	if (!want_resize && !want_rotate) return;

	struct canvas_apply_args *payload = g_new0(struct canvas_apply_args, 1);
	payload->destroy_fn = canvas_apply_args_free;
	payload->cw    = cd->pending_cw;
	payload->ch    = cd->pending_ch;
	payload->cx    = cd->pending_cx;
	payload->cy    = cd->pending_cy;
	payload->angle = cd->pending_angle;

	struct generic_layer_args *args = calloc(1, sizeof(*args));
	if (!args) { canvas_apply_args_free(payload); return; }
	args->layer_hook       = canvas_apply_hook;
	args->user             = payload;
	args->description      = g_strdup(_("Adjust canvas"));
	args->invalidate_flags = FLIS_INV_ALL;
	if (!start_in_new_thread(generic_layer_worker, args)) {
		/* Busy-guard refusal: another job owns the worker.  Free the
		 * args; nothing was mutated and no undo entry was pushed. */
		free_generic_layer_args(args);
	}
}

static void on_canvas_fit_click(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct canvas_dialog *cd = ud;
	if (!is_current_image_flis()) return;
	gboolean include_invisible =
	    siril_toggle_get_active(cd->include_invisible);
	/* Instant op — runs on the main thread; writer lock against a
	 * concurrent worker walking the stack (M-F12). */
	flis_stack_writer_lock();
	undo_save_flis_multi_layer_props(com.uniq->layers, _("Fit canvas to layers"));
	int fit_ret = flis_canvas_fit_to_layers(include_invisible);
	flis_stack_writer_unlock();
	if (fit_ret == 0)
		canvas_dialog_after_op(cd);
}

/* Rotation is now staged into pending_angle and previewed live; the
 * quick-rotate buttons increment the pending angle, the spinner edits
 * it directly, and dragging the rotation handle drives it from the
 * drawing area.  Apply commits via flis_canvas_rotate. */
static void canvas_dialog_bump_angle(struct canvas_dialog *cd, double delta) {
	cd->pending_angle += delta;
	while (cd->pending_angle >  180.0) cd->pending_angle -= 360.0;
	while (cd->pending_angle <= -180.0) cd->pending_angle += 360.0;
	canvas_dialog_sync_angle_from_pending(cd);
	if (cd->preview_da) gtk_widget_queue_draw(cd->preview_da);
}
static void on_canvas_rotate_minus_90(GtkButton *b, gpointer ud) { (void)b; canvas_dialog_bump_angle(ud,  -90.0); }
static void on_canvas_rotate_180     (GtkButton *b, gpointer ud) { (void)b; canvas_dialog_bump_angle(ud, 180.0); }
static void on_canvas_rotate_plus_90 (GtkButton *b, gpointer ud) { (void)b; canvas_dialog_bump_angle(ud,  90.0); }

static void on_canvas_angle_spin_changed(GtkSpinButton *sb, gpointer ud) {
	struct canvas_dialog *cd = ud;
	(void)sb;
	if (cd->syncing_spins) return;
	cd->pending_angle = gtk_spin_button_get_value(GTK_SPIN_BUTTON(cd->spin_angle));
	if (cd->preview_da) gtk_widget_queue_draw(cd->preview_da);
}

static void on_canvas_mirrorx_click(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct canvas_dialog *cd = ud;
	if (!is_current_image_flis()) return;
	flis_stack_writer_lock();
	undo_save_flis_multi_layer_props(com.uniq->layers, _("Mirror canvas X"));
	int mx_ret = flis_canvas_mirrorx();
	flis_stack_writer_unlock();
	if (mx_ret == 0) canvas_dialog_after_op(cd);
}
static void on_canvas_mirrory_click(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct canvas_dialog *cd = ud;
	if (!is_current_image_flis()) return;
	flis_stack_writer_lock();
	undo_save_flis_multi_layer_props(com.uniq->layers, _("Mirror canvas Y"));
	int my_ret = flis_canvas_mirrory();
	flis_stack_writer_unlock();
	if (my_ret == 0) canvas_dialog_after_op(cd);
}

static void on_canvas_dialog_close(GtkButton *btn, gpointer ud) {
	(void)btn;
	struct canvas_dialog *cd = ud;
	/* cd ownership transfers to the window's "destroy" signal handler
	 * (canvas_dialog_window_destroyed), which clears the global pointer
	 * and frees the struct.  Just trigger the destroy here so X-button
	 * close and Close-button close take the same path. */
	gtk_window_destroy(GTK_WINDOW(cd->window));
}

/* Window-destroy hook.  Fires for any reason the dialog window goes
 * away (Close button, title-bar X, programmatic destroy).  Clears the
 * singleton pointer so future external refreshes know there's no
 * dialog to repaint, and frees the struct. */
static void canvas_dialog_window_destroyed(GtkWidget *win, gpointer ud) {
	(void)win;
	struct canvas_dialog *cd = ud;
	if (g_canvas_dialog == cd) g_canvas_dialog = NULL;
	g_free(cd);
}

/* External refresh entry point.  When FLIS state changes underneath
 * the dialog (e.g. undo / redo reverts a canvas op), reset the pending
 * preview to match the now-current canvas dims, resync the spinners,
 * and queue a redraw of the preview drawing area.  No-op when the
 * dialog isn't open. */
static void canvas_dialog_refresh_external(void) {
	struct canvas_dialog *cd = g_canvas_dialog;
	if (!cd || !is_current_image_flis()) return;
	canvas_dialog_reset_pending(cd);
	cd->syncing_spins = TRUE;
	if (cd->spin_w)
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_w),
		                          (double)cd->pending_cw);
	if (cd->spin_h)
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_h),
		                          (double)cd->pending_ch);
	if (cd->spin_angle)
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_angle), 0.0);
	cd->syncing_spins = FALSE;
	canvas_dialog_refresh_current(cd);
	if (cd->preview_da)
		gtk_widget_queue_draw(cd->preview_da);
}

/* Small helper to wrap a labelled section in a frame.  Keeps the
 * dialog layout compact while making each operation visually
 * distinct. */
static GtkWidget *canvas_section(const char *title, GtkWidget *body) {
	GtkWidget *frame = gtk_frame_new(title);
	gtk_widget_set_margin_start (body, 6);
	gtk_widget_set_margin_end   (body, 6);
	gtk_widget_set_margin_top   (body, 4);
	gtk_widget_set_margin_bottom(body, 4);
	gtk_frame_set_child(GTK_FRAME(frame), body);
	return frame;
}

static void on_canvas_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	if (!is_current_image_flis()) {
		siril_log_message(_("FLIS: Canvas properties — current image is not a FLIS\n"));
		return;
	}
	if (g_canvas_dialog) {
		/* Already open — bring it forward instead of opening a second
		 * one (which would race against the singleton refresh path). */
		gtk_window_present(GTK_WINDOW(g_canvas_dialog->window));
		return;
	}

	struct canvas_dialog *cd = g_new0(struct canvas_dialog, 1);
	g_canvas_dialog = cd;
	cd->window = gtk_window_new();
	g_signal_connect(cd->window, "destroy",
	                 G_CALLBACK(canvas_dialog_window_destroyed), cd);
	gtk_window_set_title(GTK_WINDOW(cd->window), _("Canvas properties"));
	gtk_window_set_transient_for(GTK_WINDOW(cd->window),
	                             GTK_WINDOW(g_panel->window));
	gtk_window_set_default_size(GTK_WINDOW(cd->window), 380, -1);

	GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start (vbox, 10);
	gtk_widget_set_margin_end   (vbox, 10);
	gtk_widget_set_margin_top   (vbox, 10);
	gtk_widget_set_margin_bottom(vbox, 10);
	gtk_window_set_child(GTK_WINDOW(cd->window), vbox);

	cd->current_label = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(cd->current_label), 0.0f);
	gtk_widget_add_css_class(cd->current_label, "dim-label");
	canvas_dialog_refresh_current(cd);
	gtk_box_append(GTK_BOX(vbox), cd->current_label);

	/* Resize section: interactive preview + numeric controls.  The
	 * drawingarea shows the canvas frame (yellow) over light-grey layer
	 * outlines; drag inside the frame to reposition the canvas relative
	 * to the layers, drag a handle to resize.  Numeric W/H spinners
	 * mirror the frame.  Layers visually stay put; the net shift is
	 * applied to every layer's position when Apply is pressed. */
	canvas_dialog_reset_pending(cd);
	{
		GtkWidget *vb = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);

		cd->preview_da = gtk_drawing_area_new();
		gtk_widget_set_size_request(cd->preview_da, 320, 200);
		gtk_widget_set_hexpand(cd->preview_da, TRUE);
		gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(cd->preview_da),
		    on_canvas_preview_draw, cd, NULL);
		GtkGesture *drag = gtk_gesture_drag_new();
		gtk_widget_add_controller(cd->preview_da, GTK_EVENT_CONTROLLER(drag));
		g_signal_connect(drag, "drag-begin",
		                 G_CALLBACK(on_canvas_preview_drag_begin),  cd);
		g_signal_connect(drag, "drag-update",
		                 G_CALLBACK(on_canvas_preview_drag_update), cd);
		g_signal_connect(drag, "drag-end",
		                 G_CALLBACK(on_canvas_preview_drag_end),    cd);
		gtk_box_append(GTK_BOX(vb), cd->preview_da);

		GtkWidget *row = gtk_grid_new();
		gtk_grid_set_row_spacing   (GTK_GRID(row), 4);
		gtk_grid_set_column_spacing(GTK_GRID(row), 6);
		gtk_grid_attach(GTK_GRID(row), gtk_label_new(_("Width")),  0, 0, 1, 1);
		cd->spin_w = gtk_spin_button_new_with_range(1, 100000, 1);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_w), (double)cd->pending_cw);
		gtk_grid_attach(GTK_GRID(row), cd->spin_w,                 1, 0, 1, 1);
		gtk_grid_attach(GTK_GRID(row), gtk_label_new(_("Height")), 2, 0, 1, 1);
		cd->spin_h = gtk_spin_button_new_with_range(1, 100000, 1);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_h), (double)cd->pending_ch);
		gtk_grid_attach(GTK_GRID(row), cd->spin_h,                 3, 0, 1, 1);
		g_signal_connect(cd->spin_w, "value-changed",
		                 G_CALLBACK(on_canvas_size_spin_changed), cd);
		g_signal_connect(cd->spin_h, "value-changed",
		                 G_CALLBACK(on_canvas_size_spin_changed), cd);
		gtk_box_append(GTK_BOX(vb), row);

		gtk_box_append(GTK_BOX(vbox),
		    canvas_section(_("Canvas size & position"), vb));
	}

	/* Fit-to-layers section */
	{
		GtkWidget *row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		cd->include_invisible = gtk_check_button_new_with_label(
		    _("Include hidden layers"));
		GtkWidget *btn_fit = gtk_button_new_with_label(_("Fit"));
		gtk_widget_set_hexpand(cd->include_invisible, TRUE);
		gtk_widget_set_halign(btn_fit, GTK_ALIGN_END);
		gtk_box_append(GTK_BOX(row), cd->include_invisible);
		gtk_box_append(GTK_BOX(row), btn_fit);
		g_signal_connect(btn_fit, "clicked",
		                 G_CALLBACK(on_canvas_fit_click), cd);
		gtk_box_append(GTK_BOX(vbox),
		    canvas_section(_("Fit canvas to layers"), row));
	}

	/* Rotate section */
	{
		GtkWidget *vb = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
		GtkWidget *quick = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		GtkWidget *b_m90 = gtk_button_new_with_label("−90°");
		GtkWidget *b_180 = gtk_button_new_with_label("180°");
		GtkWidget *b_p90 = gtk_button_new_with_label("+90°");
		gtk_box_append(GTK_BOX(quick), b_m90);
		gtk_box_append(GTK_BOX(quick), b_180);
		gtk_box_append(GTK_BOX(quick), b_p90);
		gtk_box_append(GTK_BOX(vb), quick);

		GtkWidget *free_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		gtk_box_append(GTK_BOX(free_row), gtk_label_new(_("Angle")));
		cd->spin_angle = gtk_spin_button_new_with_range(-180.0, 180.0, 1.0);
		gtk_spin_button_set_digits(GTK_SPIN_BUTTON(cd->spin_angle), 2);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cd->spin_angle), 0.0);
		gtk_widget_set_hexpand(cd->spin_angle, TRUE);
		gtk_box_append(GTK_BOX(free_row), cd->spin_angle);
		gtk_box_append(GTK_BOX(vb), free_row);

		g_signal_connect(b_m90, "clicked", G_CALLBACK(on_canvas_rotate_minus_90), cd);
		g_signal_connect(b_180, "clicked", G_CALLBACK(on_canvas_rotate_180),      cd);
		g_signal_connect(b_p90, "clicked", G_CALLBACK(on_canvas_rotate_plus_90),  cd);
		g_signal_connect(cd->spin_angle, "value-changed",
		                 G_CALLBACK(on_canvas_angle_spin_changed), cd);
		gtk_box_append(GTK_BOX(vbox),
		    canvas_section(_("Rotate canvas"), vb));
	}

	/* Mirror section */
	{
		GtkWidget *row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
		GtkWidget *btn_h = gtk_button_new_with_label(_("Flip horizontal"));
		GtkWidget *btn_v = gtk_button_new_with_label(_("Flip vertical"));
		gtk_box_append(GTK_BOX(row), btn_h);
		gtk_box_append(GTK_BOX(row), btn_v);
		g_signal_connect(btn_h, "clicked", G_CALLBACK(on_canvas_mirrory_click), cd);
		g_signal_connect(btn_v, "clicked", G_CALLBACK(on_canvas_mirrorx_click), cd);
		gtk_box_append(GTK_BOX(vbox),
		    canvas_section(_("Mirror canvas"), row));
	}

	/* Note + Apply + Close.  Apply commits the staged resize and the
	 * staged rotation together as one undo entry; Close is destructive
	 * for whatever's still in the preview. */
	GtkWidget *note = gtk_label_new(_(
	    "Resize / position changes update layer offsets only.  Rotation "
	    "is a full-pixel operation: every layer's pixel data is "
	    "resampled and its position updated together so the composition "
	    "stays coherent.  Repeated rotations are lossy due to "
	    "resampling — undo to revert."));
	gtk_label_set_wrap(GTK_LABEL(note), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(note), 50);
	gtk_widget_add_css_class(note, "dim-label");
	gtk_box_append(GTK_BOX(vbox), note);

	GtkWidget *action_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(action_row, GTK_ALIGN_END);
	GtkWidget *apply_btn = gtk_button_new_with_label(_("Apply"));
	GtkWidget *close_btn = gtk_button_new_with_label(_("Close"));
	gtk_box_append(GTK_BOX(action_row), apply_btn);
	gtk_box_append(GTK_BOX(action_row), close_btn);
	gtk_box_append(GTK_BOX(vbox), action_row);
	g_signal_connect(apply_btn, "clicked",
	                 G_CALLBACK(on_canvas_apply_click), cd);
	g_signal_connect(close_btn, "clicked",
	                 G_CALLBACK(on_canvas_dialog_close), cd);

	gtk_window_present(GTK_WINDOW(cd->window));
}
