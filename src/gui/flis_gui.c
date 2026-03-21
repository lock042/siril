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
 * flis_gui.c — FLIS layers panel
 *
 * Provides a compact floating window for managing the FLIS layer stack.
 * Gracefully handles three states:
 *   1. No image loaded        → list shows placeholder, most controls disabled
 *   2. Plain FITS loaded      → list shows one "base layer" placeholder row;
 *                               Add button is enabled and promotes to FLIS
 *   3. FLIS image loaded      → full layer list and property controls
 *
 * Integration notes:
 *   • Call flis_gui_update() (or flis_gui_update_from_idle() from threads)
 *     after any operation that changes the layer stack.
 *   • flis_gui_open() is the public entry point; wire it to the menu item
 *     or toolbar button that shows the layers panel.
 *   • flis_gui_init() must be called once at startup (after the main UI has
 *     been loaded) to connect signals that cannot be expressed in Glade.
 *   • flis_layers.glade is processed alongside all other Siril .glade files
 *     at startup; widgets are accessed via lookup_widget() / lookup_gobject().
 */

#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/siril_app_dirs.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/single_image.h"
#include "core/undo.h"
#include "flis_gui.h"

/* =========================================================================
 * State
 * ========================================================================= */

/* TRUE while we are programmatically updating widget state, suppressing
 * change callbacks from firing back into our own logic. */
static gboolean      flis_updating   = FALSE;

/* The layer most recently selected in the list box. */
static flis_layer_t *flis_selected   = NULL;

/* =========================================================================
 * Blend mode index ↔ flis_blend_mode_t mapping
 * Must match the <items> order in flis_layers.glade exactly.
 * ========================================================================= */

static const flis_blend_mode_t blend_mode_map[] = {
    FLIS_BLEND_NORMAL,
    FLIS_BLEND_MULTIPLY,
    FLIS_BLEND_SCREEN,
    FLIS_BLEND_OVERLAY,
    FLIS_BLEND_SOFT_LIGHT,
    FLIS_BLEND_HARD_LIGHT,
    FLIS_BLEND_COLOR_DODGE,
    FLIS_BLEND_COLOR_BURN,
    FLIS_BLEND_DARKEN,
    FLIS_BLEND_LIGHTEN,
    FLIS_BLEND_DIFFERENCE,
    FLIS_BLEND_EXCLUSION,
    FLIS_BLEND_HUE,
    FLIS_BLEND_SATURATION,
    FLIS_BLEND_COLOR,
    FLIS_BLEND_LUMINOSITY,
};
static const int N_BLEND_MODES =
    (int)(sizeof(blend_mode_map) / sizeof(blend_mode_map[0]));

static int blend_mode_to_index(flis_blend_mode_t mode) {
    for (int i = 0; i < N_BLEND_MODES; i++)
        if (blend_mode_map[i] == mode) return i;
    return 0; /* default to Normal */
}

/* =========================================================================
 * Widget accessors
 *
 * All widgets were created by the startup UI loader; we reach them by name
 * via the same helpers used everywhere else in Siril.
 * ========================================================================= */

static inline GtkWidget *fw(const gchar *id) {
    return lookup_widget(id);
}

static inline GObject *fo(const gchar *id) {
    return lookup_gobject(id);
}

/* =========================================================================
 * Layer row construction
 *
 * Rows are built programmatically so we can attach the layer pointer and
 * connect per-row signals.  Each row is a GtkListBoxRow containing a
 * horizontal box with:
 *   [👁 visibility] [🔒 lock] [layer name — expands] [RGB/M badge] [⬛ mask]
 * ========================================================================= */

/* Internal callbacks for the per-row toggle buttons.
 * Declared here, defined below with other signal handlers. */
static void on_row_visibility_toggled(GtkToggleButton *btn, gpointer data);
static void on_row_lock_toggled(GtkToggleButton *btn, gpointer data);

/* Create an icon-only toggle button for use in a layer row. */
static GtkWidget *row_icon_toggle(const gchar *icon_on,
                                  const gchar *icon_off,
                                  gboolean     active,
                                  const gchar *tooltip,
                                  flis_layer_t *layer,
                                  GCallback    callback) {
    GtkWidget *btn  = gtk_toggle_button_new();
    GtkWidget *icon = gtk_image_new_from_icon_name(
        active ? icon_on : icon_off, GTK_ICON_SIZE_SMALL_TOOLBAR);

    gtk_button_set_image(GTK_BUTTON(btn), icon);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(btn), active);
    gtk_widget_set_tooltip_text(btn, tooltip);
    gtk_button_set_relief(GTK_BUTTON(btn), GTK_RELIEF_NONE);
    gtk_widget_set_focus_on_click(btn, FALSE);

    /* Attach both the layer pointer and the icon names so the callback
     * can flip the icon when the state changes. */
    g_object_set_data(G_OBJECT(btn), "flis-layer",    layer);
    g_object_set_data_full(G_OBJECT(btn), "icon-on",  g_strdup(icon_on),  g_free);
    g_object_set_data_full(G_OBJECT(btn), "icon-off", g_strdup(icon_off), g_free);
    g_signal_connect(btn, "toggled", callback, NULL);

    return btn;
}

static GtkWidget *flis_layer_row_new(flis_layer_t *layer) {
    GtkWidget *row  = gtk_list_box_row_new();
    GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
    gtk_container_set_border_width(GTK_CONTAINER(hbox), 2);

    /* Visibility toggle */
    GtkWidget *vis_btn = row_icon_toggle(
        "view-reveal-symbolic", "view-conceal-symbolic",
        layer->visible, _("Toggle layer visibility"), layer,
        G_CALLBACK(on_row_visibility_toggled));

    /* Lock toggle */
    GtkWidget *lock_btn = row_icon_toggle(
        "changes-prevent-symbolic", "changes-allow-symbolic",
        layer->locked, _("Toggle layer lock"), layer,
        G_CALLBACK(on_row_lock_toggled));

    /* Layer name label */
    GtkWidget *name_lbl = gtk_label_new(
        layer->layer_name ? layer->layer_name : _("Layer"));
    gtk_label_set_xalign(GTK_LABEL(name_lbl), 0.0f);
    gtk_widget_set_hexpand(name_lbl, TRUE);
    /* Store reference so flis_gui_update() can refresh it in-place if needed */
    g_object_set_data(G_OBJECT(row), "name-label", name_lbl);

    /* Colour model badge: "RGB" or "M" (mono) */
    gboolean mono = (layer->fit && layer->fit->naxes[2] == 1);
    GtkWidget *type_lbl = gtk_label_new(mono ? "M" : "RGB");
    gtk_style_context_add_class(
        gtk_widget_get_style_context(type_lbl), "dim-label");
    gtk_widget_set_margin_start(type_lbl, 4);
    gtk_widget_set_margin_end(type_lbl, 2);

    /* Layer mask indicator */
    GtkWidget *mask_lbl = gtk_label_new(layer->lmask ? "▪M" : "");
    gtk_widget_set_tooltip_text(mask_lbl,
        layer->lmask ? _("Layer mask present") : _("No layer mask"));
    gtk_widget_set_margin_end(mask_lbl, 4);
    g_object_set_data(G_OBJECT(row), "mask-label", mask_lbl);

    gtk_box_pack_start(GTK_BOX(hbox), vis_btn,  FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), lock_btn, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), name_lbl, TRUE,  TRUE,  0);
    gtk_box_pack_start(GTK_BOX(hbox), type_lbl, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), mask_lbl, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(row), hbox);
    gtk_widget_show_all(row);

    /* Attach the layer pointer for later retrieval in row-selected handler */
    g_object_set_data(G_OBJECT(row), "flis-layer", layer);

    return row;
}

/* =========================================================================
 * Properties panel update
 * ========================================================================= */

static void flis_properties_panel_update(flis_layer_t *layer) {
    flis_updating = TRUE;

    GtkWidget *props = fw("flis_props_box");

    if (!layer) {
        gtk_widget_set_sensitive(props, FALSE);
        flis_updating = FALSE;
        return;
    }

    gtk_widget_set_sensitive(props, TRUE);

    /* Name */
    gtk_entry_set_text(GTK_ENTRY(fw("flis_name_entry")),
                       layer->layer_name ? layer->layer_name : "");

    /* Blend mode */
    gtk_combo_box_set_active(GTK_COMBO_BOX(fw("flis_blend_combo")),
                             blend_mode_to_index(layer->blend_mode));

    /* Opacity (scale and spin share the same adjustment) */
    GtkAdjustment *adj = GTK_ADJUSTMENT(fo("flis_opacity_adj"));
    gtk_adjustment_set_value(adj, (gdouble)(layer->opacity * 100.f));

    /* Mask section */
    gboolean has_mask = (layer->lmask != NULL);
    gtk_button_set_label(GTK_BUTTON(fw("flis_mask_toggle_btn")),
                         has_mask ? _("Remove") : _("Add…"));
    gtk_label_set_text(GTK_LABEL(fw("flis_mask_status_label")),
                       has_mask ? _("Present") : _("None"));
    gtk_widget_set_sensitive(fw("flis_mask_move_btn"), has_mask);

    /* Tint section — always visible, but insensitive for RGB layers */
    gboolean mono = (layer->fit && layer->fit->naxes[2] == 1);
    gboolean has_tint = mono && layer->has_tint;
    gtk_widget_set_sensitive(fw("flis_tint_frame"), mono);
    gtk_toggle_button_set_active(
        GTK_TOGGLE_BUTTON(fw("flis_tint_check")), has_tint);
    gtk_widget_set_sensitive(fw("flis_tint_color_btn"), mono && has_tint);
    if (has_tint) {
        GdkRGBA colour = {
            layer->layer_tint.r,
            layer->layer_tint.g,
            layer->layer_tint.b,
            1.0
        };
        gtk_color_chooser_set_rgba(
            GTK_COLOR_CHOOSER(fw("flis_tint_color_btn")), &colour);
    }

    flis_updating = FALSE;
}

/* =========================================================================
 * Button sensitivity update
 *
 * Called after any stack change to reflect what is and isn't possible
 * for the currently selected layer.
 * ========================================================================= */

static void flis_toolbar_sensitivity_update(void) {
    gboolean have_image = (com.uniq != NULL && com.uniq->fit != NULL);
    gboolean is_flis    = is_current_image_flis();
    gboolean have_sel   = (flis_selected != NULL);
    gint     n          = flis_layer_count();

    /* Add is always available when an image is loaded */
    gtk_widget_set_sensitive(fw("flis_add_btn"), have_image);

    /* Operations on a selected layer */
    gtk_widget_set_sensitive(fw("flis_remove_btn"),
                             is_flis && have_sel && n > 1);
    /* Duplicate is available whenever an image is loaded: for a plain FITS
     * it promotes to FLIS first, then duplicates the base layer. */
    gtk_widget_set_sensitive(fw("flis_duplicate_btn"),
                             have_image && (have_sel || !is_flis));

    /* Move buttons: only active if there is a layer above/below */
    gboolean can_up = FALSE, can_down = FALSE;
    if (is_flis && have_sel) {
        gint idx = flis_layer_get_index(flis_selected);
        can_up   = (idx >= 0 && idx < n - 1);
        can_down = (idx > 0);
    }
    gtk_widget_set_sensitive(fw("flis_move_up_btn"),   can_up);
    gtk_widget_set_sensitive(fw("flis_move_down_btn"), can_down);

    /* Mode badge */
    gtk_label_set_text(GTK_LABEL(fw("flis_mode_label")),
                       is_flis ? "FLIS" : "FITS");
}

/* =========================================================================
 * Layer list rebuild
 *
 * Clears and repopulates the GtkListBox from com.uniq->layers (or shows
 * a placeholder when no layers are present).  After rebuilding, re-selects
 * the row corresponding to flis_selected if it still exists.
 * ========================================================================= */

static void flis_layers_list_rebuild(void) {
    GtkListBox *list = GTK_LIST_BOX(fw("flis_layer_list"));

    /* Remove all existing rows without triggering row-selected callbacks */
    flis_updating = TRUE;
    GList *children = gtk_container_get_children(GTK_CONTAINER(list));
    for (GList *c = children; c; c = c->next)
        gtk_container_remove(GTK_CONTAINER(list), GTK_WIDGET(c->data));
    g_list_free(children);
    flis_updating = FALSE;

    if (!is_current_image_flis()) {
        /* Plain FITS or no image: show a single informational placeholder row */
        GtkWidget *row  = gtk_list_box_row_new();
        GtkWidget *lbl  = gtk_label_new(
            com.uniq && com.uniq->fit
            ? _("Single layer — click + to convert to FLIS")
            : _("No image loaded"));
        gtk_widget_set_sensitive(row, FALSE);
        gtk_label_set_xalign(GTK_LABEL(lbl), 0.5f);
        gtk_widget_set_margin_top(lbl, 8);
        gtk_widget_set_margin_bottom(lbl, 8);
        gtk_style_context_add_class(
            gtk_widget_get_style_context(lbl), "dim-label");
        gtk_container_add(GTK_CONTAINER(row), lbl);
        gtk_widget_show_all(row);
        gtk_list_box_insert(list, row, -1);
        gtk_widget_show_all(GTK_WIDGET(list));

        flis_selected = NULL;
        flis_properties_panel_update(NULL);
        flis_toolbar_sensitivity_update();
        return;
    }

    /* Build rows in reverse layer_order (top of visual stack at top of list).
     * com.uniq->layers is sorted ascending by layer_order, so the last element
     * is the topmost layer. */
    gint n = flis_layer_count();
    GtkListBoxRow *row_to_select = NULL;

    for (gint i = n - 1; i >= 0; i--) {
        flis_layer_t *lay = (flis_layer_t *)g_slist_nth_data(
                                com.uniq->layers, (guint)i);
        if (!lay) continue;
        GtkWidget *row = flis_layer_row_new(lay);
        gtk_list_box_insert(list, row, -1);
        if (lay == flis_selected)
            row_to_select = GTK_LIST_BOX_ROW(row);
    }

    gtk_widget_show_all(GTK_WIDGET(list));

    /* Re-select the previously selected layer, or default to the topmost */
    flis_updating = TRUE;
    if (row_to_select) {
        gtk_list_box_select_row(list, row_to_select);
    } else {
        GtkListBoxRow *first = gtk_list_box_get_row_at_index(list, 0);
        if (first) {
            gtk_list_box_select_row(list, first);
            flis_selected = (flis_layer_t *)g_object_get_data(
                                G_OBJECT(first), "flis-layer");
        }
    }
    flis_updating = FALSE;

    flis_properties_panel_update(flis_selected);
    flis_toolbar_sensitivity_update();
}

/* =========================================================================
 * File helpers
 * ========================================================================= */

/* Open a GTK file chooser for FITS files. Returns a heap-allocated path
 * that the caller must g_free(), or NULL if the user cancelled. */
static gchar *flis_choose_fits_file(GtkWindow *parent,
                                    const gchar *title) {
    GtkWidget *dialog = gtk_file_chooser_dialog_new(
        title, parent, GTK_FILE_CHOOSER_ACTION_OPEN,
        _("Cancel"), GTK_RESPONSE_CANCEL,
        _("Open"),   GTK_RESPONSE_ACCEPT,
        NULL);

    GtkFileFilter *filter = gtk_file_filter_new();
    gtk_file_filter_set_name(filter, _("FITS images"));
    gtk_file_filter_add_pattern(filter, "*.fit");
    gtk_file_filter_add_pattern(filter, "*.fits");
    gtk_file_filter_add_pattern(filter, "*.fts");
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

    GtkFileFilter *all = gtk_file_filter_new();
    gtk_file_filter_set_name(all, _("All files"));
    gtk_file_filter_add_pattern(all, "*");
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), all);

    gchar *filename = NULL;
    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

    gtk_widget_destroy(dialog);
    return filename;
}

/* Return TRUE if the given file has FLIS=T in its primary HDU.
 * We need to catch this because readfits() would call load_flis() and
 * clobber the current layer stack. */
static gboolean file_is_flis(const gchar *filename) {
    fitsfile *fptr   = NULL;
    int       status = 0;
    int       flag   = 0;

    fits_open_diskfile(&fptr, filename, READONLY, &status);
    if (status) return FALSE;

    fits_movabs_hdu(fptr, 1, NULL, &status);
    status = 0;
    fits_read_key(fptr, TLOGICAL, "FLIS", &flag, NULL, &status);
    fits_close_file(fptr, &status);

    return flag != 0;
}

/* Load a FITS file into a heap-allocated fits* suitable for flis_layer_add().
 * Returns NULL on failure. */
static fits *flis_load_layer_file(const gchar *filename) {
    if (file_is_flis(filename)) {
        siril_log_color_message(
            _("Cannot use a FLIS file as a layer source.\n"), "salmon");
        return NULL;
    }

    fits *f = calloc(1, sizeof(fits));
    if (!f) { PRINT_ALLOC_ERR; return NULL; }

    if (readfits(filename, f, NULL, FALSE)) {
        siril_log_color_message(
            _("Failed to load layer file: %s\n"), "red", filename);
        clearfits(f);
        free(f);
        return NULL;
    }
    return f;
}

/* =========================================================================
 * Per-row signal handlers (visibility and lock toggle buttons)
 * ========================================================================= */

static void on_row_visibility_toggled(GtkToggleButton *btn,
                                      gpointer         data) {
    (void)data;
    if (flis_updating) return;

    flis_layer_t *layer = (flis_layer_t *)g_object_get_data(
                              G_OBJECT(btn), "flis-layer");
    if (!layer) return;

    gboolean active = gtk_toggle_button_get_active(btn);
    undo_save_flis_layer_props(layer, _("Layer visibility"));
    flis_layer_set_visible(layer, active);

    /* Flip the button icon */
    const gchar *icon_name = active
        ? (const gchar *)g_object_get_data(G_OBJECT(btn), "icon-on")
        : (const gchar *)g_object_get_data(G_OBJECT(btn), "icon-off");
    gtk_button_set_image(GTK_BUTTON(btn),
        gtk_image_new_from_icon_name(icon_name, GTK_ICON_SIZE_SMALL_TOOLBAR));

    flis_invalidate_composite();
    queue_redraw(REMAP_ALL);
}

static void on_row_lock_toggled(GtkToggleButton *btn,
                                gpointer         data) {
    (void)data;
    if (flis_updating) return;

    flis_layer_t *layer = (flis_layer_t *)g_object_get_data(
                              G_OBJECT(btn), "flis-layer");
    if (!layer) return;

    gboolean locked = gtk_toggle_button_get_active(btn);
    undo_save_flis_layer_props(layer, _("Layer lock"));
    flis_layer_set_locked(layer, locked);

    const gchar *icon_name = locked
        ? (const gchar *)g_object_get_data(G_OBJECT(btn), "icon-on")
        : (const gchar *)g_object_get_data(G_OBJECT(btn), "icon-off");
    gtk_button_set_image(GTK_BUTTON(btn),
        gtk_image_new_from_icon_name(icon_name, GTK_ICON_SIZE_SMALL_TOOLBAR));

    /* Lock state doesn't affect the composite */
    flis_toolbar_sensitivity_update();
}

/* =========================================================================
 * Layer list signal handlers
 * ========================================================================= */

G_MODULE_EXPORT void on_flis_layer_row_selected(GtkListBox    *box,
                                                GtkListBoxRow *row,
                                                gpointer       data) {
    (void)box; (void)data;
    if (flis_updating) return;

    if (!row) {
        flis_selected = NULL;
        flis_properties_panel_update(NULL);
        flis_toolbar_sensitivity_update();
        return;
    }

    flis_layer_t *layer = (flis_layer_t *)g_object_get_data(
                              G_OBJECT(row), "flis-layer");
    flis_selected = layer;

    if (layer && com.uniq) {
        gint idx = flis_layer_get_index(layer);
        if (idx >= 0)
            uniq_set_active_layer(com.uniq, idx);
    }

    flis_properties_panel_update(layer);
    flis_toolbar_sensitivity_update();

    /* Redraw so the display reflects the newly active layer.  The composite
     * itself does not change on a simple layer selection (all visible layers
     * are always composited), but this ensures the display is always current
     * and responsive when the user interacts with the panel. */
    if (is_current_image_flis() && !flis_updating)
        queue_redraw(REMAP_ALL);
}

/* =========================================================================
 * Toolbar signal handlers
 * ========================================================================= */

G_MODULE_EXPORT void on_flis_add_layer_clicked(GtkButton *btn,
                                               gpointer   data) {
    (void)btn; (void)data;

    gchar *filename = flis_choose_fits_file(
        GTK_WINDOW(lookup_widget("flis_layers_window")), _("Choose layer file"));
    if (!filename) return;

    fits *f = flis_load_layer_file(filename);
    if (!f) { g_free(filename); return; }

    /* Derive a sensible layer name from the filename */
    gchar *base = g_path_get_basename(filename);
    /* Strip extension: find last dot */
    gchar *dot = g_strrstr(base, ".");
    if (dot) *dot = '\0';

    /* If no FLIS is currently loaded, promote the existing FITS first */
    if (!is_current_image_flis()) {
        if (flis_promote_from_gfit(_("Background"))) {
            siril_log_color_message(
                _("FLIS: failed to convert image to FLIS\n"), "red");
            clearfits(f); free(f); g_free(base); g_free(filename);
            return;
        }
    }

    flis_layer_t *new_layer = flis_layer_add(f, base);
    g_free(base);
    g_free(filename);

    if (!new_layer) {
        siril_log_color_message(_("FLIS: failed to add layer\n"), "red");
        return;
    }

    flis_selected = new_layer;
//    flis_undo_notify_structural_change(); // Removed as part of layer structure undo implementation
    flis_invalidate_composite();
    flis_gui_update();
    queue_redraw(REMAP_ALL);
}

G_MODULE_EXPORT void on_flis_remove_layer_clicked(GtkButton *btn,
                                                  gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected) return;

    /* Confirmation dialog */
    GtkWidget *dlg = gtk_message_dialog_new(
        GTK_WINDOW(lookup_widget("flis_layers_window")),
        GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
        GTK_MESSAGE_QUESTION, GTK_BUTTONS_OK_CANCEL,
        _("Remove layer \"%s\"?"),
        flis_selected->layer_name ? flis_selected->layer_name : _("Layer"));
    gtk_message_dialog_format_secondary_text(
        GTK_MESSAGE_DIALOG(dlg),
        _("Pixel-level undo states for this layer will be discarded. "
          "States for other layers are preserved."));
    gint response = gtk_dialog_run(GTK_DIALOG(dlg));
    gtk_widget_destroy(dlg);
    if (response != GTK_RESPONSE_OK) return;

    flis_layer_t *to_remove = flis_selected;
    flis_selected = NULL;

    /* Purge undo states for this layer BEFORE destroying it so item_id is
     * still valid.  Other layers' states are preserved. */
    flis_undo_purge_layer(to_remove->item_id);

    if (flis_layer_remove(to_remove)) return; /* error already logged */

    flis_invalidate_composite();
    flis_gui_update();
    queue_redraw(REMAP_ALL);
}

G_MODULE_EXPORT void on_flis_duplicate_layer_clicked(GtkButton *btn,
                                                     gpointer   data) {
    (void)btn; (void)data;

    /* If a plain FITS is loaded, promote it to a single-layer FLIS first
     * so that flis_layer_duplicate() has a valid layer stack to work with.
     * This mirrors the behaviour of the Add button. */
    if (!is_current_image_flis()) {
        if (!com.uniq || !com.uniq->fit) return;
        if (flis_promote_from_gfit(_("Background"))) {
            siril_log_color_message(
                _("FLIS: failed to convert image to FLIS\n"), "red");
            return;
        }
        /* After promotion there is exactly one layer; select it so that
         * flis_layer_duplicate() below has a valid flis_selected. */
        flis_selected = flis_active_layer();
    }

    if (!flis_selected) return;

    flis_layer_t *dup = flis_layer_duplicate(flis_selected);
    if (!dup) return;

    flis_selected = dup;
//    flis_undo_notify_structural_change();
    flis_invalidate_composite();
    flis_gui_update();
    queue_redraw(REMAP_ALL);
}

G_MODULE_EXPORT void on_flis_move_up_clicked(GtkButton *btn,
                                             gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected) return;

    if (flis_layer_move_up(flis_selected) == 0) {
//        flis_undo_notify_structural_change();
        flis_invalidate_composite();
        flis_gui_update();
        queue_redraw(REMAP_ALL);
    }
}

G_MODULE_EXPORT void on_flis_move_down_clicked(GtkButton *btn,
                                               gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected) return;

    if (flis_layer_move_down(flis_selected) == 0) {
//        flis_undo_notify_structural_change();
        flis_invalidate_composite();
        flis_gui_update();
        queue_redraw(REMAP_ALL);
    }
}

/* =========================================================================
 * Properties panel signal handlers
 * ========================================================================= */

/* Apply the name entry to the selected layer (on Enter or focus-out) */
static void flis_apply_name_entry(void) {
    if (flis_updating || !flis_selected) return;

    const gchar *text = gtk_entry_get_text(
                            GTK_ENTRY(fw("flis_name_entry")));
    undo_save_flis_layer_props(flis_selected, _("Layer name"));
    if (flis_layer_set_name(flis_selected, text) == 0) {
        /* Update the name label in the list row without a full rebuild */
        GtkListBox    *list = GTK_LIST_BOX(fw("flis_layer_list"));
        GList *rows = gtk_container_get_children(GTK_CONTAINER(list));
        for (GList *r = rows; r; r = r->next) {
            GtkWidget *row = GTK_WIDGET(r->data);
            if (g_object_get_data(G_OBJECT(row), "flis-layer") == flis_selected) {
                GtkWidget *lbl = GTK_WIDGET(
                    g_object_get_data(G_OBJECT(row), "name-label"));
                if (lbl) gtk_label_set_text(GTK_LABEL(lbl), text);
                break;
            }
        }
        g_list_free(rows);
        /* Name change doesn't affect the composite pixels */
    }
}

G_MODULE_EXPORT void on_flis_name_activate(GtkEntry *entry, gpointer data) {
    (void)entry; (void)data;
    flis_apply_name_entry();
}

G_MODULE_EXPORT gboolean on_flis_name_focus_out(GtkWidget *widget,
                                                GdkEventFocus *event,
                                                gpointer data) {
    (void)widget; (void)event; (void)data;
    flis_apply_name_entry();
    return FALSE; /* allow default focus handling */
}

G_MODULE_EXPORT void on_flis_blend_mode_changed(GtkComboBox *combo,
                                                gpointer     data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    gint idx = gtk_combo_box_get_active(combo);
    if (idx < 0 || idx >= N_BLEND_MODES) return;

    undo_save_flis_layer_props(flis_selected, _("Blend mode"));
    flis_layer_set_blend_mode(flis_selected, blend_mode_map[idx]);
    flis_invalidate_composite();
    queue_redraw(REMAP_ALL);
}

/* Opacity: connect to the GtkAdjustment value-changed signal in code so
 * it fires exactly once regardless of which widget (scale or spin) was
 * used. Connected in flis_gui_open(). */
/* Opacity debouncing state. */
static gint               flis_opacity_drag_id       = FLIS_UNDO_LAYER_NONE;
static flis_layer_props_t flis_opacity_drag_snapshot;

static gboolean on_flis_opacity_scale_press(GtkWidget *widget,
                                            GdkEventButton *event,
                                            gpointer data) {
    (void)widget; (void)event; (void)data;
    if (!flis_selected) return FALSE;

    flis_opacity_drag_id                   = flis_selected->item_id;
    flis_opacity_drag_snapshot.blend_mode  = flis_selected->blend_mode;
    flis_opacity_drag_snapshot.opacity     = flis_selected->opacity;
    flis_opacity_drag_snapshot.visible     = flis_selected->visible;
    flis_opacity_drag_snapshot.locked      = flis_selected->locked;
    flis_opacity_drag_snapshot.has_tint    = flis_selected->has_tint;
    flis_opacity_drag_snapshot.tint        = flis_selected->layer_tint;
    g_strlcpy(flis_opacity_drag_snapshot.name,
              flis_selected->layer_name ? flis_selected->layer_name : "",
              sizeof(flis_opacity_drag_snapshot.name));
    return FALSE;
}

static gboolean on_flis_opacity_scale_release(GtkWidget *widget,
                                              GdkEventButton *event,
                                              gpointer data) {
    (void)widget; (void)event; (void)data;
    if (flis_opacity_drag_id == FLIS_UNDO_LAYER_NONE) return FALSE;

    if (flis_selected && flis_selected->item_id == flis_opacity_drag_id)
        undo_save_flis_layer_props_snapshot(flis_opacity_drag_id,
                                            &flis_opacity_drag_snapshot,
                                            _("Opacity"));
    flis_opacity_drag_id = FLIS_UNDO_LAYER_NONE;
    return FALSE;
}

static void on_flis_opacity_adj_changed(GtkAdjustment *adj, gpointer data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    gdouble pct = gtk_adjustment_get_value(adj);
    flis_layer_set_opacity(flis_selected, (gfloat)(pct / 100.0));
    flis_invalidate_composite();
    queue_redraw(REMAP_ALL);
}

static void on_flis_opacity_spin_commit(GtkSpinButton *spin, gpointer data) {
    (void)spin; (void)data;
    if (flis_updating || !flis_selected) return;
    if (flis_opacity_drag_id != FLIS_UNDO_LAYER_NONE) return;
    undo_save_flis_layer_props(flis_selected, _("Opacity"));
}

/* =========================================================================
 * Layer mask signal handlers
 * ========================================================================= */

G_MODULE_EXPORT void on_flis_mask_toggle_clicked(GtkButton *btn,
                                                 gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected) return;

    if (flis_selected->lmask) {
        /* Remove existing mask */
        GtkWidget *dlg = gtk_message_dialog_new(
            GTK_WINDOW(lookup_widget("flis_layers_window")),
            GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
            GTK_MESSAGE_QUESTION, GTK_BUTTONS_OK_CANCEL,
            _("Remove the layer mask from \"%s\"?"),
            flis_selected->layer_name ? flis_selected->layer_name : _("Layer"));
        gint r = gtk_dialog_run(GTK_DIALOG(dlg));
        gtk_widget_destroy(dlg);
        if (r != GTK_RESPONSE_OK) return;

        /* Save mask pixels before removal so undo can restore them */
        undo_save_flis_lmask(flis_selected, _("Remove layer mask"));
        flis_layer_remove_lmask(flis_selected);
    } else {
        /* Add mask from file */
        gchar *filename = flis_choose_fits_file(
            GTK_WINDOW(lookup_widget("flis_layers_window")),
            _("Choose greyscale mask file"));
        if (!filename) return;

        fits *mf = flis_load_layer_file(filename);
        g_free(filename);
        if (!mf) return;

        /* Convert to 8-bit greyscale layermask_t */
        guint w = mf->rx, h = mf->ry;
        size_t npix = (size_t)w * h;

        layermask_t *lm = calloc(1, sizeof(layermask_t));
        if (!lm) { PRINT_ALLOC_ERR; clearfits(mf); free(mf); return; }

        lm->w      = w;
        lm->h      = h;
        lm->bitpix = 8;
        lm->data   = malloc(npix);
        if (!lm->data) {
            PRINT_ALLOC_ERR;
            free(lm); clearfits(mf); free(mf); return;
        }

        /* Scale the first channel to 8-bit */
        uint8_t *dst = (uint8_t *)lm->data;
        if (mf->type == DATA_FLOAT && mf->fdata) {
            for (size_t i = 0; i < npix; i++)
                dst[i] = (uint8_t)(mf->fdata[i] * 255.f + 0.5f);
        } else if (mf->type == DATA_USHORT && mf->data) {
            for (size_t i = 0; i < npix; i++)
                dst[i] = (uint8_t)(mf->data[i] >> 8);
        } else {
            memset(dst, 255, npix); /* fallback: white (fully opaque) */
        }
        clearfits(mf); free(mf);

        /* Save "no mask" marker before adding, so undo removes the new mask */
        undo_save_flis_lmask(flis_selected, _("Add layer mask"));

        if (flis_layer_set_lmask(flis_selected, lm)) {
            /* Error (e.g. size mismatch) was already logged */
            layermask_free(lm);
            return;
        }
    }

    flis_invalidate_composite();
    flis_gui_update();
    queue_redraw(REMAP_ALL);
}

G_MODULE_EXPORT void on_flis_mask_move_clicked(GtkButton *btn,
                                               gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected || !flis_selected->lmask) return;

    gint n = flis_layer_count();
    if (n < 2) {
        siril_log_message(_("FLIS: no other layers to move the mask to.\n"));
        return;
    }

    /* Build a dialog with a combo listing compatible destination layers */
    GtkWidget *dlg = gtk_dialog_new_with_buttons(
        _("Move Layer Mask"),
        GTK_WINDOW(lookup_widget("flis_layers_window")),
        GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
        _("Cancel"), GTK_RESPONSE_CANCEL,
        _("Move"),   GTK_RESPONSE_OK,
        NULL);

    GtkWidget *content = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_container_set_border_width(GTK_CONTAINER(content), 12);

    GtkWidget *lbl = gtk_label_new(_("Destination layer:"));
    gtk_label_set_xalign(GTK_LABEL(lbl), 0.0f);
    gtk_box_pack_start(GTK_BOX(content), lbl, FALSE, FALSE, 4);

    GtkWidget *combo = gtk_combo_box_text_new();
    /* Layers are sorted ascending by layer_order; show highest first */
    for (gint i = n - 1; i >= 0; i--) {
        flis_layer_t *lay = (flis_layer_t *)g_slist_nth_data(
                                com.uniq->layers, (guint)i);
        if (!lay || lay == flis_selected) continue;

        /* Only show layers whose dimensions match the mask */
        if ((size_t)lay->fit->rx != flis_selected->lmask->w ||
            (size_t)lay->fit->ry != flis_selected->lmask->h) continue;

        gchar *label = g_strdup_printf("%s (%dx%d)",
            lay->layer_name ? lay->layer_name : _("Layer"),
            lay->fit->rx, lay->fit->ry);
        gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo),
                                  g_strdup_printf("%d", lay->item_id),
                                  label);
        g_free(label);
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    gtk_box_pack_start(GTK_BOX(content), combo, FALSE, FALSE, 4);
    gtk_widget_show_all(content);

    gint response = gtk_dialog_run(GTK_DIALOG(dlg));

    if (response == GTK_RESPONSE_OK) {
        const gchar *id_str =
            gtk_combo_box_get_active_id(GTK_COMBO_BOX(combo));
        if (id_str) {
            gint target_id = atoi(id_str);
            flis_layer_t *dest = flis_layer_get_by_id(target_id);
            if (dest) {
                /* Save both layers' mask states before the move.
                 * Order matters: undo is LIFO, so we save dest first and
                 * source second.  Undo then restores source first (mask
                 * goes back to source) then dest (mask removed from dest). */
                undo_save_flis_lmask(dest,          _("Move layer mask (dest)"));
                undo_save_flis_lmask(flis_selected, _("Move layer mask (source)"));
                flis_layer_move_lmask(flis_selected, dest);
                flis_invalidate_composite();
                flis_gui_update();
                queue_redraw(REMAP_ALL);
            }
        }
    }
    gtk_widget_destroy(dlg);
}

/* =========================================================================
 * Tint signal handlers
 * ========================================================================= */

/*
 * Returns TRUE if the FLIS composite will produce RGB output — i.e. if any
 * visible layer is either an RGB layer or a mono layer with a tint set.
 * A pure mono stack with no tints composites to mono; the moment any tint
 * is enabled the composite becomes RGB and the channel tabs must be updated.
 *
 * Called before and after a tint change to detect a colour-model transition
 * so that close_tab() / init_right_tab() are only triggered when necessary,
 * avoiding redundant GUI redraws.
 */
static gboolean flis_composite_will_be_rgb(void) {
    if (!com.uniq || !com.uniq->layers) return FALSE;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->fit || !lay->visible) continue;
        if (lay->fit->naxes[2] >= 3) return TRUE;   /* RGB layer */
        if (lay->has_tint)           return TRUE;   /* tinted mono */
    }
    return FALSE;
}

G_MODULE_EXPORT void on_flis_tint_toggled(GtkToggleButton *btn,
                                          gpointer         data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    gboolean was_rgb = flis_composite_will_be_rgb();

    undo_save_flis_layer_props(flis_selected, _("Tint"));
    gboolean enable = gtk_toggle_button_get_active(btn);
    gtk_widget_set_sensitive(fw("flis_tint_color_btn"), enable);

    if (!enable) {
        flis_layer_clear_tint(flis_selected);
    } else {
        GdkRGBA colour;
        gtk_color_chooser_get_rgba(
            GTK_COLOR_CHOOSER(fw("flis_tint_color_btn")), &colour);
        flis_layer_set_tint(flis_selected,
                            colour.red, colour.green, colour.blue);
    }

    gboolean now_rgb = flis_composite_will_be_rgb();
    if (now_rgb != was_rgb) {
        if (com.uniq)
            com.uniq->chans = now_rgb ? 3 : 1;
        close_tab(NULL);
        init_right_tab(NULL);
    }

    flis_invalidate_composite();
    queue_redraw(REMAP_ALL);
}

G_MODULE_EXPORT void on_flis_tint_color_set(GtkColorButton *btn,
                                            gpointer        data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    GdkRGBA colour;
    gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(btn), &colour);
    undo_save_flis_layer_props(flis_selected, _("Tint colour"));
    flis_layer_set_tint(flis_selected,
                        colour.red, colour.green, colour.blue);
    flis_invalidate_composite();
    queue_redraw(REMAP_ALL);
}

/* =========================================================================
 * Window signal handler
 * ========================================================================= */

G_MODULE_EXPORT gboolean on_flis_layers_window_delete_event(
    GtkWidget *widget, GdkEventAny *event, gpointer data) {
    (void)widget; (void)event; (void)data;
    /* Hide rather than destroy so state is preserved */
    gtk_widget_hide(lookup_widget("flis_layers_window"));
    return TRUE; /* prevent actual destruction */
}

/* =========================================================================
 * Public API
 * ========================================================================= */

static gboolean flis_gui_update_idle(gpointer data) {
    (void)data;
    flis_gui_update();
    return G_SOURCE_REMOVE;
}

void flis_gui_update_from_idle(void) {
    siril_add_idle(flis_gui_update_idle, NULL);
}

/**
 * flis_gui_update:
 *
 * Rebuilds the layer list and updates all property widgets to match the
 * current state of com.uniq.  Safe to call at any time — including when
 * no image is loaded and when the layers panel is not currently visible.
 * Must be called from the GTK main thread; use flis_gui_update_from_idle()
 * from worker threads.
 *
 * Call this from any code that opens or closes an image so the panel
 * stays in sync, e.g. from open_single_image_from_gfit() and
 * free_image_data_gui().
 */
void flis_gui_update(void) {
    if (!lookup_widget("flis_layers_window")) return;
    flis_layers_list_rebuild();
}

/**
 * flis_gui_init:
 *
 * One-time initialisation called once after the main UI has been loaded at
 * startup.  Connects the opacity GtkAdjustment signal in code so that it
 * fires exactly once per user gesture regardless of whether the GtkScale
 * or the GtkSpinButton triggered the change.  (Glade does not allow
 * connecting adjustment signals directly.)
 *
 * Call this from wherever Siril performs its post-load GUI setup.
 */
void flis_gui_init(void) {
    GObject *adj = lookup_gobject("flis_opacity_adj");
    if (adj)
        g_signal_connect(adj, "value-changed",
                         G_CALLBACK(on_flis_opacity_adj_changed), NULL);

    GtkWidget *scale = lookup_widget("flis_opacity_scale");
    if (scale) {
        g_signal_connect(scale, "button-press-event",
                         G_CALLBACK(on_flis_opacity_scale_press), NULL);
        g_signal_connect(scale, "button-release-event",
                         G_CALLBACK(on_flis_opacity_scale_release), NULL);
    }

    GtkWidget *spin = lookup_widget("flis_opacity_spin");
    if (spin)
        g_signal_connect(spin, "value-changed",
                         G_CALLBACK(on_flis_opacity_spin_commit), NULL);
}

void flis_gui_open(void) {
    GtkWidget *win = lookup_widget("flis_layers_window");
    if (!win) {
        siril_log_color_message(
            _("FLIS GUI: layers window not found — "
              "is flis_layers.glade included in the UI bundle?\n"), "red");
        return;
    }

    /* Set the transient parent every time we open the dialog.  This tells
     * the window manager which application window the dialog belongs to,
     * which (combined with keep-above=True in the UI file) ensures it
     * floats above the main window without covering other applications.
     * Doing this every time handles the case where the main window is
     * re-created (rare but possible). */
    GtkWindow *main_win = GTK_WINDOW(lookup_widget("control_window"));
    if (main_win)
        gtk_window_set_transient_for(GTK_WINDOW(win), main_win);

    /* keep-above must be set in code; GTK's UI loader rejects it as a
     * static property on GtkDialog. */
    gtk_window_set_keep_above(GTK_WINDOW(win), TRUE);

    gtk_widget_show_all(win);
    gtk_window_present(GTK_WINDOW(win));
    flis_layers_list_rebuild();
}
