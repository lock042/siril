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

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/gui_iface.h"
#include "io/image_format_flis.h"
#include "io/single_image.h"

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

/* =========================================================================
 * Panel singleton
 * ========================================================================= */

struct flis_panel {
	GtkWidget *window;
	GtkWidget *mode_label;          /* FITS / FLIS / FLIS·active group */

	/* Toolbar */
	GtkWidget *btn_add, *btn_remove, *btn_duplicate, *btn_group;
	GtkWidget *btn_drag, *btn_move_up, *btn_move_down;
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
	GtkAdjustment *opacity_adj;
	float      opacity_drag_start;  /* property snapshot for undo */
	gboolean   opacity_dragging;

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

	/* Refresh suppression: TRUE while we're programmatically syncing
	 * widget state to the model, so signal handlers don't loop back. */
	gboolean   refreshing;
	/* Pending idle refresh ID (0 = none).  Coalesces multiple
	 * flis_gui_update_from_idle calls during a worker burst. */
	guint      refresh_idle_id;
};

static struct flis_panel *g_panel;  /* NULL until first show */

/* Forward decls of the section builders / handlers. */
static void   build_panel(void);
static void   build_toolbar  (GtkWidget *box);
static void   build_list     (GtkWidget *box);
static void   build_property (GtkWidget *box);
static void   build_mask     (GtkWidget *box);
static GMenu *build_context_menu(void);
static void   register_panel_actions(void);
static void   refresh_panel(void);
static void   sync_property_widgets(flis_layer_t *lay);
static void   update_toolbar_sensitivity(void);
static void   ensure_flis_css(void);
static flis_layer_t *current_selected_layer(void);
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

static gboolean refresh_idle_cb(gpointer p) {
	(void)p;
	if (g_panel) {
		g_panel->refresh_idle_id = 0;
		if (gtk_widget_get_visible(g_panel->window))
			refresh_panel();
	}
	return G_SOURCE_REMOVE;
}

void flis_gui_update_from_idle(void) {
	if (!g_panel) return;                /* panel not yet built — nothing to refresh */
	if (g_panel->refresh_idle_id != 0) return;  /* already pending */
	g_panel->refresh_idle_id = g_idle_add(refresh_idle_cb, NULL);
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
	gtk_window_set_default_size(GTK_WINDOW(w), 360, 720);
	gtk_window_set_hide_on_close(GTK_WINDOW(w), TRUE);
	GtkWidget *main_w = lookup_widget("control_window");
	if (main_w)
		gtk_window_set_transient_for(GTK_WINDOW(w), GTK_WINDOW(main_w));

	GtkWidget *outer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start (outer, 6);
	gtk_widget_set_margin_end   (outer, 6);
	gtk_widget_set_margin_top   (outer, 6);
	gtk_widget_set_margin_bottom(outer, 6);
	gtk_window_set_child(GTK_WINDOW(w), outer);

	/* Mode indicator at the top. */
	g_panel->mode_label = gtk_label_new("FITS");
	gtk_label_set_xalign(GTK_LABEL(g_panel->mode_label), 0.0f);
	gtk_widget_add_css_class(g_panel->mode_label, "dim-label");
	gtk_box_append(GTK_BOX(outer), g_panel->mode_label);

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
static void on_drag_toggled      (GtkToggleButton *b, gpointer u);
static void on_move_up_clicked   (GtkButton *b, gpointer u);
static void on_move_down_clicked (GtkButton *b, gpointer u);

static void build_toolbar(GtkWidget *box) {
	GtkWidget *bar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
	gtk_widget_add_css_class(bar, "toolbar");

	g_panel->btn_add       = icon_button("list-add-symbolic",        _("Add layer"));
	g_panel->btn_remove    = icon_button("list-remove-symbolic",     _("Remove selected layer"));
	g_panel->btn_duplicate = icon_button("edit-copy-symbolic",       _("Duplicate selected layer"));
	g_panel->btn_group     = icon_button("folder-new-symbolic",      _("Create group"));
	g_panel->btn_drag      = icon_toggle("input-mouse-symbolic",     _("Drag layer in canvas"));
	g_panel->btn_move_up   = icon_button("go-up-symbolic",           _("Move layer up"));
	g_panel->btn_move_down = icon_button("go-down-symbolic",         _("Move layer down"));

	g_signal_connect(g_panel->btn_add,       "clicked", G_CALLBACK(on_add_clicked),       NULL);
	g_signal_connect(g_panel->btn_remove,    "clicked", G_CALLBACK(on_remove_clicked),    NULL);
	g_signal_connect(g_panel->btn_duplicate, "clicked", G_CALLBACK(on_duplicate_clicked), NULL);
	g_signal_connect(g_panel->btn_group,     "clicked", G_CALLBACK(on_group_clicked),     NULL);
	g_signal_connect(g_panel->btn_drag,      "toggled", G_CALLBACK(on_drag_toggled),      NULL);
	g_signal_connect(g_panel->btn_move_up,   "clicked", G_CALLBACK(on_move_up_clicked),   NULL);
	g_signal_connect(g_panel->btn_move_down, "clicked", G_CALLBACK(on_move_down_clicked), NULL);

	gtk_box_append(GTK_BOX(bar), g_panel->btn_add);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_remove);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_duplicate);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_group);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_drag);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_move_up);
	gtk_box_append(GTK_BOX(bar), g_panel->btn_move_down);

	/* Context menu button on the right. */
	GtkWidget *spacer = gtk_label_new(NULL);
	gtk_widget_set_hexpand(spacer, TRUE);
	gtk_box_append(GTK_BOX(bar), spacer);

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
static void on_group_expander_toggled(GtkToggleButton *btn, gpointer u);
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
	start_in_new_thread(generic_layer_worker, args);
	return TRUE;
}

static void on_row_setup(GtkListItemFactory *f, GtkListItem *item, gpointer u) {
	(void)f; (void)u;
	row_widgets_t *rw = g_new0(row_widgets_t, 1);

	rw->row_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_margin_start (rw->row_box, 4);
	gtk_widget_set_margin_end   (rw->row_box, 4);
	gtk_widget_set_margin_top   (rw->row_box, 2);
	gtk_widget_set_margin_bottom(rw->row_box, 2);

	rw->expander = gtk_toggle_button_new();
	gtk_button_set_icon_name(GTK_BUTTON(rw->expander), "pan-down-symbolic");
	gtk_widget_set_tooltip_text(rw->expander, _("Toggle group collapse"));
	gtk_widget_add_css_class(rw->expander, "flat");

	rw->visible_toggle = gtk_toggle_button_new();
	gtk_button_set_icon_name(GTK_BUTTON(rw->visible_toggle), "view-reveal-symbolic");
	gtk_widget_set_tooltip_text(rw->visible_toggle, _("Toggle visibility"));
	gtk_widget_add_css_class(rw->visible_toggle, "flat");

	rw->lock_toggle = gtk_toggle_button_new();
	gtk_button_set_icon_name(GTK_BUTTON(rw->lock_toggle), "changes-prevent-symbolic");
	gtk_widget_set_tooltip_text(rw->lock_toggle, _("Lock layer (prevent edits)"));
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
		gtk_label_set_text(GTK_LABEL(rw->name_label),
		                   lay->layer_name ? lay->layer_name : "(unnamed)");

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
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rw->expander), !grp->collapsed);
		gchar *markup = g_markup_printf_escaped(
			"\xF0\x9F\x97\x80 <b>%s</b>",   /* U+1F5C0 file folder */
			grp->name ? grp->name : _("(group)"));
		gtk_label_set_markup(GTK_LABEL(rw->name_label), markup);
		g_free(markup);
		gtk_label_set_text(GTK_LABEL(rw->kind_badge), "");

		rw->vis_handler_id = g_signal_connect(rw->visible_toggle, "toggled",
		                         G_CALLBACK(on_row_visible_toggled),
		                         GINT_TO_POINTER(-ri->item_id));  /* negative encodes group */
		rw->exp_handler_id = g_signal_connect(rw->expander, "toggled",
		                         G_CALLBACK(on_group_expander_toggled),
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
			src->group_id = 0;
			gui_iface.flis_display_invalidate(FLIS_INV_STACK, src_id);
			gui_iface.flis_gui_update();
			notify_gfit_data_modified();
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
	start_in_new_thread(generic_layer_worker, args);
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
static void on_blend_changed     (GtkDropDown *dd, GParamSpec *p, gpointer u);
static void on_opacity_changed   (GtkAdjustment *adj, gpointer u);
static void on_opacity_drag_begin(GtkGestureDrag *g, gdouble x, gdouble y, gpointer u);
static void on_opacity_drag_end  (GtkGestureDrag *g, gdouble dx, gdouble dy, gpointer u);
static void on_tint_check_toggled(GtkCheckButton *btn, gpointer u);
static void on_tint_color_chosen (GtkColorDialogButton *btn, GParamSpec *p, gpointer u);

/* 18 GSK-equivalent FLIS blend modes — the dropdown's index order maps
 * to the FLIS_BLEND_* values via blend_mode_for_index(). */
static const char *blend_names[] = {
	"Normal",       "Multiply",   "Screen",      "Overlay",
	"Soft Light",   "Hard Light", "Color Dodge", "Color Burn",
	"Darken",       "Lighten",    "Difference",  "Exclusion",
	"Hue",          "Saturation", "Color",       "Luminosity",
	"Chroma (LRGB)", "Pass-through",
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
	gtk_grid_set_row_spacing   (GTK_GRID(grid), 6);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
	gtk_widget_set_margin_start (grid, 8);
	gtk_widget_set_margin_end   (grid, 8);
	gtk_widget_set_margin_top   (grid, 8);
	gtk_widget_set_margin_bottom(grid, 8);

	int row = 0;

	/* Name */
	gtk_grid_attach(GTK_GRID(grid), gtk_label_new(_("Name")), 0, row, 1, 1);
	g_panel->name_entry = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(g_panel->name_entry), 32);
	gtk_widget_set_hexpand(g_panel->name_entry, TRUE);
	gtk_grid_attach(GTK_GRID(grid), g_panel->name_entry, 1, row, 2, 1);
	g_signal_connect(g_panel->name_entry, "activate",
	                 G_CALLBACK(on_name_activate), NULL);
	row++;

	/* Blend mode */
	gtk_grid_attach(GTK_GRID(grid), gtk_label_new(_("Blend")), 0, row, 1, 1);
	GtkStringList *bl = gtk_string_list_new(NULL);
	for (int i = 0; i < N_BLENDS; i++)
		gtk_string_list_append(bl, blend_names[i]);
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
	gtk_grid_attach(GTK_GRID(grid), g_panel->opacity_scale, 1, row, 1, 1);
	gtk_grid_attach(GTK_GRID(grid), g_panel->opacity_spin,  2, row, 1, 1);
	g_signal_connect(g_panel->opacity_adj, "value-changed",
	                 G_CALLBACK(on_opacity_changed), NULL);
	/* Drag-begin/end gestures on the scale for undo snapshot batching. */
	GtkGesture *drag = gtk_gesture_drag_new();
	gtk_widget_add_controller(g_panel->opacity_scale, GTK_EVENT_CONTROLLER(drag));
	g_signal_connect(drag, "drag-begin", G_CALLBACK(on_opacity_drag_begin), NULL);
	g_signal_connect(drag, "drag-end",   G_CALLBACK(on_opacity_drag_end),   NULL);
	row++;

	/* Tint frame */
	g_panel->tint_frame = gtk_frame_new(_("Tint"));
	GtkWidget *tbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_margin_start (tbox, 6);
	gtk_widget_set_margin_end   (tbox, 6);
	gtk_widget_set_margin_top   (tbox, 6);
	gtk_widget_set_margin_bottom(tbox, 6);
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
	GtkWidget *mbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_widget_set_margin_start (mbox, 6);
	gtk_widget_set_margin_end   (mbox, 6);
	gtk_widget_set_margin_top   (mbox, 6);
	gtk_widget_set_margin_bottom(mbox, 6);
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

/* ---- Context menu -------------------------------------------------- */

/* Stub action callback shared by all context menu items.  The body is
 * replaced piecewise as the wiring lands; for now each prints a
 * targeted "TODO: <verb>" log so users see why the menu item didn't
 * do what they expected.  Merge Down and Flatten Image have working
 * primitives (flis_merge_down_layer / flis_flatten_all) so they
 * dispatch through the worker directly. */
static void on_ctx_merge_down(GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_flatten   (GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_stub          (GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_move_to_group (GSimpleAction *a, GVariant *v, gpointer u);
static void on_ctx_export_layer  (GSimpleAction *a, GVariant *v, gpointer u);

static GMenu *build_context_menu(void) {
	GMenu *m = g_menu_new();
	g_menu_append(m, _("Export current layer as FITS…"), "win.flis-export-layer");
	g_menu_append(m, _("Register layers…"),              "win.flis-register-layers");
	g_menu_append(m, _("Layers match…"),                 "win.flis-layers-match");
	g_menu_append(m, _("Move layer to group…"),          "win.flis-move-to-group");
	g_menu_append(m, _("Merge Down"),                    "win.flis-merge-down");
	g_menu_append(m, _("Flatten Image"),                 "win.flis-flatten");
	return m;
}

static const GActionEntry flis_panel_actions[] = {
	{ "flis-export-layer",   on_ctx_export_layer, NULL, NULL, NULL },
	{ "flis-register-layers",on_ctx_stub,       NULL, NULL, NULL },
	{ "flis-layers-match",   on_ctx_stub,       NULL, NULL, NULL },
	{ "flis-move-to-group",  on_ctx_move_to_group, NULL, NULL, NULL },
	{ "flis-merge-down",     on_ctx_merge_down, NULL, NULL, NULL },
	{ "flis-flatten",        on_ctx_flatten,    NULL, NULL, NULL },
};

static void register_panel_actions(void) {
	GtkWidget *main_w = lookup_widget("control_window");
	if (!main_w || !G_IS_ACTION_MAP(main_w)) return;
	g_action_map_add_action_entries(G_ACTION_MAP(main_w),
	                                 flis_panel_actions,
	                                 G_N_ELEMENTS(flis_panel_actions),
	                                 NULL);
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

static void refresh_panel(void) {
	if (!g_panel) return;
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
}

static void sync_property_widgets(flis_layer_t *lay) {
	/* Programmatic widget state changes below — guard against the
	 * cascade of "changed" signals firing dispatch_op for what is
	 * really just a display refresh.  Caller may have already set
	 * the flag (refresh_panel does); set it locally just in case. */
	const gboolean was_refreshing = g_panel->refreshing;
	g_panel->refreshing = TRUE;

	gtk_widget_set_sensitive(g_panel->prop_frame, lay != NULL);
	gtk_widget_set_sensitive(g_panel->mask_frame, lay != NULL);
	if (!lay) {
		gtk_editable_set_text(GTK_EDITABLE(g_panel->name_entry), "");
		gtk_button_set_label(GTK_BUTTON(g_panel->mask_status_btn), _("(no mask)"));
		g_panel->refreshing = was_refreshing;
		return;
	}
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
	 * mask exist (chooses which one the global mask tab shows). */
	const gboolean both_masks = (lay->lmask != NULL && gfit && gfit->pdata[0] != NULL
	                              /* TODO: replace with real proc-mask test once mask plumbing arrives */
	                              && FALSE);
	gtk_widget_set_visible(g_panel->mask_view_row, both_masks);

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
	OP_GROUP_SET_VISIBLE,
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
};

static void op_payload_free(gpointer p) {
	struct op_payload *op = p;
	if (!op) return;
	g_free(op->str_v);
	g_free(op);
}

static int op_hook(struct generic_layer_args *args) {
	struct op_payload *op = (struct op_payload *)args->user;
	flis_layer_t *lay = (op->kind == OP_GROUP_SET_VISIBLE) ? NULL
	                    : flis_layer_get_by_id(op->target_id);
	flis_group_t *grp = (op->kind == OP_GROUP_SET_VISIBLE)
	                    ? flis_group_get_by_id(op->target_id) : NULL;
	switch (op->kind) {
		case OP_SET_VISIBLE:  return lay ? flis_layer_set_visible (lay, op->bool_v) : 1;
		case OP_SET_LOCKED:   return lay ? flis_layer_set_locked  (lay, op->bool_v) : 1;
		case OP_SET_NAME:     return lay ? flis_layer_set_name    (lay, op->str_v) : 1;
		case OP_SET_BLEND:    return lay ? flis_layer_set_blend_mode(lay, op->blend_v) : 1;
		case OP_SET_OPACITY:  return lay ? flis_layer_set_opacity (lay, op->float_v) : 1;
		case OP_SET_TINT:     return lay ? flis_layer_set_tint    (lay, op->r, op->g, op->b) : 1;
		case OP_CLEAR_TINT:   if (lay) { lay->has_tint = FALSE; return 0; } return 1;
		case OP_GROUP_SET_VISIBLE: return grp ? flis_group_set_visible(grp, op->bool_v) : 1;
		case OP_LAYER_REMOVE: return lay ? flis_layer_remove(lay) : 1;
		case OP_LAYER_DUPLICATE: {
			if (!lay) return 1;
			flis_layer_t *dup = flis_layer_duplicate(lay);
			return dup ? 0 : 1;
		}
		case OP_LAYER_MOVE_UP:   return lay ? flis_layer_move_up(lay)   : 1;
		case OP_LAYER_MOVE_DOWN: return lay ? flis_layer_move_down(lay) : 1;
	}
	return 1;
}

static void dispatch_op(struct op_payload *op, const char *desc,
                         int invalidate_flags) {
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer        = NULL;     /* hook resolves layer from op->target_id */
	args->layer_hook   = op_hook;
	args->user         = op;
	args->description  = g_strdup(desc);
	args->verbose      = FALSE;
	args->invalidate_flags   = invalidate_flags;
	args->invalidate_item_id = op->target_id;
	start_in_new_thread(generic_layer_worker, args);
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
static void on_group_expander_toggled(GtkToggleButton *btn, gpointer u) {
	if (g_panel && g_panel->refreshing) return;
	const gint gid = GPOINTER_TO_INT(u);
	flis_group_t *grp = flis_group_get_by_id(gid);
	if (!grp) return;
	grp->collapsed = !gtk_toggle_button_get_active(btn);  /* expanded = collapsed FALSE */
	refresh_panel();
}

static void on_selection_changed(GtkSelectionModel *sel, guint pos, guint nitems, gpointer u) {
	(void)sel; (void)pos; (void)nitems; (void)u;
	if (g_panel && g_panel->refreshing) return;
	flis_layer_t *lay = current_selected_layer();
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
			uniq_set_active_layer(com.uniq, index);
			gui_iface.redraw_image(REMAP_ALL);
			/* Refresh so the active-layer-row CSS class moves to the
			 * newly active row.  Selection is preserved across the
			 * rebuild by refresh_panel's snapshot/restore. */
			refresh_panel();
		}
	}
}

/* Property handlers */
static void on_name_activate(GtkEntry *e, gpointer u) {
	(void)u;
	if (g_panel && g_panel->refreshing) return;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = lay->item_id;
	op->kind = OP_SET_NAME;
	op->str_v = g_strdup(gtk_editable_get_text(GTK_EDITABLE(e)));
	dispatch_op(op, _("Layer name"), FLIS_INV_LAYER_PROPS);
}

static void on_blend_changed(GtkDropDown *dd, GParamSpec *p, gpointer u) {
	(void)p; (void)u;
	if (g_panel && g_panel->refreshing) return;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	int idx = (int)gtk_drop_down_get_selected(dd);
	if (idx < 0 || idx >= N_BLENDS) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = lay->item_id;
	op->kind = OP_SET_BLEND;
	op->blend_v = blend_mode_for_index_table[idx];
	dispatch_op(op, _("Layer blend mode"), FLIS_INV_LAYER_PROPS);
}

static void on_opacity_changed(GtkAdjustment *adj, gpointer u) {
	(void)u;
	if (g_panel && g_panel->refreshing) return;
	if (g_panel && g_panel->opacity_dragging) {
		/* While dragging, update layer state live but don't push undo
		 * per tick — the drag-end handler does that once.  We need a
		 * full notify_gfit_data_modified, not just a paint queue:
		 * the CPU composite cache + histogram + per-vport tile
		 * buffers all hold stale pixels until the composite is
		 * rebuilt.  This is per-tick during drag and could be heavy
		 * on huge FLIS images; if that becomes a problem we can
		 * throttle to every N ms — for now correctness wins. */
		flis_layer_t *lay = current_selected_layer();
		if (lay) lay->opacity = (gfloat)(gtk_adjustment_get_value(adj) / 100.0);
		gui_iface.flis_display_invalidate(FLIS_INV_LAYER_PROPS,
		                                   lay ? lay->item_id : 0);
		notify_gfit_data_modified();
		return;
	}
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = lay->item_id;
	op->kind = OP_SET_OPACITY;
	op->float_v = (gfloat)(gtk_adjustment_get_value(adj) / 100.0);
	dispatch_op(op, _("Layer opacity"), FLIS_INV_LAYER_PROPS);
}

static void on_opacity_drag_begin(GtkGestureDrag *g, gdouble x, gdouble y, gpointer u) {
	(void)g; (void)x; (void)y; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay || !g_panel) return;
	g_panel->opacity_drag_start = lay->opacity;
	g_panel->opacity_dragging   = TRUE;
}

static void on_opacity_drag_end(GtkGestureDrag *g, gdouble dx, gdouble dy, gpointer u) {
	(void)g; (void)dx; (void)dy; (void)u;
	if (!g_panel) return;
	g_panel->opacity_dragging = FALSE;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) return;
	const gfloat final_v = lay->opacity;
	if (final_v == g_panel->opacity_drag_start) return;  /* no net change */
	/* Stage the end-of-drag value through the worker so undo is recorded once. */
	lay->opacity = g_panel->opacity_drag_start;  /* reset so the worker sees a real change */
	struct op_payload *op = g_new0(struct op_payload, 1);
	op->destroy_fn = op_payload_free;
	op->target_id = lay->item_id;
	op->kind = OP_SET_OPACITY;
	op->float_v = final_v;
	dispatch_op(op, _("Layer opacity (drag)"), FLIS_INV_LAYER_PROPS);
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
		op->r = op->g = op->b = 1.0;   /* neutral broadcast until colour chosen */
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
	dispatch_op(op, _("Layer tint colour"), FLIS_INV_LAYER_PIXELS);
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
	start_in_new_thread(generic_layer_worker, args);
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
		lay->layer_name ? lay->layer_name : "(unnamed)");
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
	start_in_new_thread(generic_layer_worker, args);
}
static void on_drag_toggled(GtkToggleButton *b, gpointer u) {
	(void)u;
	siril_log_message(_("FLIS: Canvas drag mode %s — TODO\n"),
	                  gtk_toggle_button_get_active(b) ? _("on") : _("off"));
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
	start_in_new_thread(generic_layer_worker, args);
}

static void on_move_up_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (lay) { move_relative(lay, TRUE); return; }
	flis_group_t *grp = current_selected_group();
	if (grp) {
		/* Moving a group as a whole needs a batch reorder primitive
		 * — deferred to a follow-up.  For now, log and refuse so the
		 * user gets clear feedback instead of a partial move. */
		siril_log_message(_("FLIS: moving a group as a whole is not yet "
		                    "implemented; select an individual layer\n"));
	}
}
static void on_move_down_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (lay) { move_relative(lay, FALSE); return; }
	flis_group_t *grp = current_selected_group();
	if (grp) {
		siril_log_message(_("FLIS: moving a group as a whole is not yet "
		                    "implemented; select an individual layer\n"));
	}
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
	start_in_new_thread(generic_layer_worker, args);
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
	start_in_new_thread(generic_layer_worker, args);
}

static void on_mask_status_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay || !lay->lmask) return;
	/* Toggle active state directly on the layer struct — no primitive
	 * for this yet; flis_layer_set_lmask doesn't have an active-only
	 * mode.  The display invalidation does the right thing. */
	lay->lmask_active = !lay->lmask_active;
	gui_iface.flis_display_invalidate(FLIS_INV_LAYER_PIXELS, lay->item_id);
	gui_iface.redraw_image(REMAP_ALL);
	flis_gui_update_from_idle();
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

static void on_mask_move_clicked(GtkButton *b, gpointer u) {
	(void)b; (void)u;
	/* Move mask between layers — needs a layer-chooser dialog.
	 * Deferred until §4.3 slice 3 lands the Move-to-group infrastructure,
	 * which uses the same dialog pattern. */
	siril_log_message(_("FLIS: Move mask — not yet implemented\n"));
}
static void on_mask_view_radio_toggled(GtkCheckButton *btn, gpointer u) {
	(void)btn; (void)u;
	siril_log_message(_("FLIS: Mask view radio toggled — TODO\n"));
}

/* Context menu action handlers ---------------------------------- */

static int op_merge_down_hook(struct generic_layer_args *args) {
	flis_layer_t *lay = flis_layer_get_by_id(args->invalidate_item_id);
	if (!lay) return 1;
	return flis_merge_down_layer(lay);
}
static void on_ctx_merge_down(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	flis_layer_t *lay = current_selected_layer();
	if (!lay) {
		siril_log_message(_("FLIS: Merge Down — no layer selected\n"));
		return;
	}
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer = NULL;
	args->layer_hook = op_merge_down_hook;
	args->description = g_strdup(_("Merge Down"));
	args->invalidate_flags   = FLIS_INV_ALL;
	args->invalidate_item_id = lay->item_id;
	start_in_new_thread(generic_layer_worker, args);
}

static int op_flatten_hook(struct generic_layer_args *args) {
	(void)args;
	return flis_flatten_all();
}
static void on_ctx_flatten(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)a; (void)v; (void)u;
	if (!is_current_image_flis()) {
		siril_log_message(_("FLIS: Flatten — current image is not a FLIS\n"));
		return;
	}
	struct generic_layer_args *args = calloc(1, sizeof(*args));
	args->layer = NULL;
	args->layer_hook = op_flatten_hook;
	args->description = g_strdup(_("Flatten Image"));
	args->invalidate_flags = FLIS_INV_ALL;
	start_in_new_thread(generic_layer_worker, args);
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
	start_in_new_thread(generic_layer_worker, args);

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
		gtk_string_list_append(items, grp->name ? grp->name : "(unnamed)");
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
	start_in_new_thread(generic_layer_worker, args);
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

static void on_ctx_stub(GSimpleAction *a, GVariant *v, gpointer u) {
	(void)v; (void)u;
	const char *name = g_action_get_name(G_ACTION(a));
	siril_log_message(_("FLIS: %s — not yet implemented (needs §4.3 / §5)\n"), name);
}
