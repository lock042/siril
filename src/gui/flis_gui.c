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
#include "io/sequence.h"
#include "algos/siril_wcs.h"
#include "registration/registration.h"
#include "gui/message_dialog.h"
#include "core/undo.h"
#include "flis_gui.h"
#include "gui/image_interactions.h"

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
    FLIS_BLEND_CHROMA,
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
void on_flis_mask_view_toggled(GtkToggleButton *btn, gpointer data);
void on_flis_lmask_active_toggle_clicked(GtkButton *btn, gpointer data);
void on_flis_move_layer_toggled(GtkToggleButton *btn, gpointer data);

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

    /* Layer mask indicator — shows "⊟" (mask present, bold if active) or
     * "○" (no mask). */
    GtkWidget *mask_lbl = gtk_label_new(NULL);
    if (layer->lmask && layer->lmask_active)
        gtk_label_set_markup(GTK_LABEL(mask_lbl), "<b>⊟</b>");
    else if (layer->lmask)
        gtk_label_set_text(GTK_LABEL(mask_lbl), "⊟");
    else
        gtk_label_set_text(GTK_LABEL(mask_lbl), "○");
    gtk_widget_set_tooltip_text(mask_lbl,
        layer->lmask
            ? (layer->lmask_active ? _("Layer mask present (active)")
                                   : _("Layer mask present (inactive)"))
            : _("No layer mask"));
    gtk_widget_set_margin_end(mask_lbl, 4);
    if (!layer->lmask) {
        gtk_style_context_add_class(
            gtk_widget_get_style_context(mask_lbl), "dim-label");
    }
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
    gboolean has_lmask    = (layer->lmask != NULL);
    gboolean has_procmask = (layer->fit && layer->fit->mask != NULL);
    gboolean has_both     = has_lmask && has_procmask;
    gtk_button_set_label(GTK_BUTTON(fw("flis_mask_toggle_btn")),
                         has_lmask ? _("Remove") : _("Add…"));
    /* Status button: "Active" (bold) / "Inactive" / "None"; sensitive when mask present */
    {
        GtkWidget *status_btn = fw("flis_mask_status_btn");
        if (has_lmask) {
            GtkLabel *lbl = GTK_LABEL(gtk_bin_get_child(GTK_BIN(status_btn)));
            if (layer->lmask_active)
                gtk_label_set_markup(lbl, _("<b>Active</b>"));
            else
                gtk_label_set_text(lbl, _("Inactive"));
        } else {
            gtk_button_set_label(GTK_BUTTON(status_btn), _("None"));
        }
        gtk_widget_set_sensitive(status_btn, has_lmask);
    }
    gtk_widget_set_sensitive(fw("flis_mask_move_btn"), has_lmask);

    /* Show the view-toggle row only when both mask types are present */
    gtk_widget_set_visible(fw("flis_mask_view_row"), has_both);
    if (has_both) {
        /* Sync radio buttons to current display state; block signals to avoid
         * triggering a redraw during the panel update. */
        GtkWidget *proc_radio  = fw("flis_mask_view_proc_radio");
        GtkWidget *layer_radio = fw("flis_mask_view_layer_radio");
        g_signal_handlers_block_by_func(proc_radio,  on_flis_mask_view_toggled, NULL);
        g_signal_handlers_block_by_func(layer_radio, on_flis_mask_view_toggled, NULL);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(proc_radio),
                                     !get_flis_show_layer_mask());
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(layer_radio),
                                     get_flis_show_layer_mask());
        g_signal_handlers_unblock_by_func(proc_radio,  on_flis_mask_view_toggled, NULL);
        g_signal_handlers_unblock_by_func(layer_radio, on_flis_mask_view_toggled, NULL);
    }

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

    /* Drag-to-move toggle — enabled whenever a FLIS layer is selected */
    gtk_widget_set_sensitive(fw("flis_drag_toggle_btn"), is_flis && have_sel);

    /* Layer actions menu — enabled whenever a FLIS layer is selected */
    gtk_widget_set_sensitive(fw("flis_layer_menu_btn"), is_flis && have_sel);

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

    /* Re-select the previously selected layer, or default to the topmost.
     * flis_updating suppresses on_flis_layer_row_selected during the
     * programmatic selection, so we must manually sync gfit afterwards. */
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

    /* Sync gfit with the newly selected layer.  on_flis_layer_row_selected
     * was suppressed by flis_updating above, so uniq_set_active_layer was
     * not called.  Without this, gfit keeps pointing at whichever layer was
     * active before the rebuild (typically the base layer from load_flis),
     * and every subsequent operation modifies the wrong layer regardless of
     * which one is visually selected in the panel. */
    if (flis_selected && com.uniq) {
        gint idx = flis_layer_get_index(flis_selected);
        if (idx >= 0)
            uniq_set_active_layer(com.uniq, idx);
    }

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

    /* Start in the current working directory */
    if (com.wd)
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), com.wd);

    gchar *filename = NULL;
    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

    gtk_widget_destroy(dialog);
    return filename;
}

/* Open a GTK file chooser for saving a FITS file.  suggested_name is the
 * pre-filled leaf filename (no directory).  Returns a heap-allocated path
 * the caller must g_free(), or NULL if the user cancelled. */
static gchar *flis_choose_fits_save(GtkWindow   *parent,
                                    const gchar *title,
                                    const gchar *suggested_name) {
    GtkWidget *dialog = gtk_file_chooser_dialog_new(
        title, parent, GTK_FILE_CHOOSER_ACTION_SAVE,
        _("Cancel"), GTK_RESPONSE_CANCEL,
        _("Save"),   GTK_RESPONSE_ACCEPT,
        NULL);
    gtk_file_chooser_set_do_overwrite_confirmation(
        GTK_FILE_CHOOSER(dialog), TRUE);

    GtkFileFilter *filter = gtk_file_filter_new();
    gtk_file_filter_set_name(filter, _("FITS images"));
    gtk_file_filter_add_pattern(filter, "*.fit");
    gtk_file_filter_add_pattern(filter, "*.fits");
    gtk_file_filter_add_pattern(filter, "*.fts");
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

    if (suggested_name)
        gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog),
                                          suggested_name);

    /* Start in the same directory as the current image if possible */
    if (com.uniq && com.uniq->filename) {
        gchar *dir = g_path_get_dirname(com.uniq->filename);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), dir);
        g_free(dir);
    }

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
 * Forward declarations
 * ========================================================================= */

/* Defined in Tint signal handlers section below */
static gboolean flis_composite_will_be_rgb(void);

/* =========================================================================
 * generic_layer_worker user-data structs and hook functions
 *
 * Each operation that modifies the layer stack or composite and needs a
 * queue_redraw() goes through generic_layer_worker so it does not race
 * with redraw().  The pattern mirrors generic_image_worker:
 *   1. GTK main thread: read widget state, save undo, allocate + fill args
 *   2. start_in_new_thread(generic_layer_worker, args)
 *   3. Worker thread: call layer_hook(args)
 *   4. GTK idle (end_generic_layer or custom): invalidate composite, redraw
 *
 * Operations that do NOT need the worker:
 *   - on_row_lock_toggled     — lock state doesn't affect composite/redraw
 *   - flis_apply_name_entry   — name doesn't affect composite pixels
 *   - on_flis_opacity_adj_changed — fires every drag tick; keep direct
 * ========================================================================= */

/* ---- Visibility ---- */
typedef struct {
    destructor  destroy_fn;
    gboolean    visible;
} flis_visibility_args_t;

static int flis_visibility_hook(struct generic_layer_args *args) {
    flis_visibility_args_t *a = args->user;
    flis_layer_set_visible(args->layer, a->visible);
    return 0;
}

/* ---- Blend mode ---- */
typedef struct {
    destructor         destroy_fn;
    flis_blend_mode_t  mode;
} flis_blend_args_t;

static int flis_blend_hook(struct generic_layer_args *args) {
    flis_blend_args_t *a = args->user;
    flis_layer_set_blend_mode(args->layer, a->mode);
    return 0;
}

/* ---- Move up / move down — no extra args needed ---- */
static int flis_move_up_hook(struct generic_layer_args *args) {
    return flis_layer_move_up(args->layer);
}

static int flis_move_down_hook(struct generic_layer_args *args) {
    return flis_layer_move_down(args->layer);
}

/* ---- Lmask toggle: add (lm != NULL) or remove (lm == NULL) ---- */
typedef struct {
    destructor    destroy_fn;
    layermask_t  *lm;   /* non-NULL = add, ownership taken by hook; NULL = remove */
} flis_lmask_toggle_args_t;

static void flis_lmask_toggle_destroy(void *p) {
    flis_lmask_toggle_args_t *a = p;
    /* lm is NULLed by the hook on success (ownership transferred).
     * Free it only if the hook never ran or failed. */
    if (a->lm) layermask_free(a->lm);
    free(a);
}

static int flis_lmask_toggle_hook(struct generic_layer_args *args) {
    flis_lmask_toggle_args_t *a = args->user;
    if (a->lm) {
        if (flis_layer_set_lmask(args->layer, a->lm)) return 1;
        a->lm = NULL;  /* ownership transferred to the layer */
    } else {
        flis_layer_remove_lmask(args->layer);
    }
    return 0;
}

/* ---- Lmask move ---- */
typedef struct {
    destructor    destroy_fn;
    flis_layer_t *dest;
} flis_lmask_move_args_t;

static int flis_lmask_move_hook(struct generic_layer_args *args) {
    flis_lmask_move_args_t *a = args->user;
    flis_layer_move_lmask(args->layer, a->dest);
    return 0;
}

/* ---- Tint toggle — needs custom end idle for tab management ---- */
typedef struct {
    destructor  destroy_fn;
    gboolean    enable;
    double      r, g, b;
    gboolean    was_rgb;
} flis_tint_toggle_args_t;

static int flis_tint_toggle_hook(struct generic_layer_args *args) {
    flis_tint_toggle_args_t *a = args->user;
    if (!a->enable)
        flis_layer_clear_tint(args->layer);
    else
        flis_layer_set_tint(args->layer, a->r, a->g, a->b);
    return 0;
}

/* Custom end idle for tint toggle: updates channel tabs when the composite
 * colour model flips between mono and RGB. */
static gboolean end_flis_tint_toggle(gpointer p) {
    struct generic_layer_args *args = p;
    flis_tint_toggle_args_t *a = args->user;
    stop_processing_thread();

    if (!args->retval && is_current_image_flis()) {
        flis_invalidate_composite();
        gboolean now_rgb = flis_composite_will_be_rgb();
        if (now_rgb != a->was_rgb) {
            if (com.uniq) com.uniq->chans = now_rgb ? 3 : 1;
            close_tab(NULL);
            init_right_tab(NULL);
        }
        flis_gui_update();
        queue_redraw(REMAP_ALL);
    }

    free_generic_layer_args(args);
    return FALSE;
}

/* ---- Tint colour ---- */
typedef struct {
    destructor  destroy_fn;
    double      r, g, b;
} flis_tint_color_args_t;

static int flis_tint_color_hook(struct generic_layer_args *args) {
    flis_tint_color_args_t *a = args->user;
    flis_layer_set_tint(args->layer, a->r, a->g, a->b);
    return 0;
}

/* Helper: allocate and launch a generic_layer_worker.
 * Takes ownership of user_data (freed via destroy_fn or plain free).
 * On thread-busy failure, frees args+user and returns FALSE. */
static gboolean flis_launch_layer_op(flis_layer_t *layer,
                                     int (*hook)(struct generic_layer_args *),
                                     void *user_data,
                                     const gchar *description,
                                     gboolean updates_lmask,
                                     gboolean (*idle_fn)(gpointer)) {
    struct generic_layer_args *args = calloc(1, sizeof(struct generic_layer_args));
    if (!args) { PRINT_ALLOC_ERR; return FALSE; }
    args->layer        = layer;
    args->layer_hook   = hook;
    args->user         = user_data;
    args->description  = (gchar *)description;
    args->verbose      = TRUE;
    args->updates_lmask= updates_lmask;
    args->idle_function= idle_fn;
    if (!start_in_new_thread(generic_layer_worker, args)) {
        free_generic_layer_args(args);
        return FALSE;
    }
    return TRUE;
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

    /* Flip the button icon immediately on the main thread (cosmetic only) */
    const gchar *icon_name = active
        ? (const gchar *)g_object_get_data(G_OBJECT(btn), "icon-on")
        : (const gchar *)g_object_get_data(G_OBJECT(btn), "icon-off");
    gtk_button_set_image(GTK_BUTTON(btn),
        gtk_image_new_from_icon_name(icon_name, GTK_ICON_SIZE_SMALL_TOOLBAR));

    undo_save_flis_layer_props(layer, _("Layer visibility"));

    flis_visibility_args_t *user = calloc(1, sizeof(flis_visibility_args_t));
    if (!user) { PRINT_ALLOC_ERR; return; }
    user->visible = active;

    flis_launch_layer_op(layer, flis_visibility_hook, user,
                         _("Layer visibility"), FALSE, NULL);
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


void on_flis_layer_row_selected(GtkListBox    *box,
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

    /* Auto-select which mask type to display for the newly active layer:
     * prefer processing mask if present; show layer mask only if it is the
     * sole mask available.  This resets any manual toggle the user made on
     * the previous layer, which is the least surprising behaviour. */
    if (layer) {
        gboolean has_proc  = (layer->fit && layer->fit->mask != NULL);
        gboolean has_lmask = (layer->lmask != NULL);
        set_flis_show_layer_mask(has_lmask && !has_proc);
    }

    /* Update mask tab visibility, title, and enable-check sensitivity. */
    show_or_hide_mask_tab();

    /* Populate the mask tab for the newly active layer. */
    queue_redraw_mask();

    /* Redraw so the display reflects the newly active layer. */
    if (is_current_image_flis() && !flis_updating)
        queue_redraw(REMAP_ALL);
}

/* =========================================================================
 * Toolbar signal handlers
 * ========================================================================= */

void on_flis_add_layer_clicked(GtkButton *btn,
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
        clearfits(f); free(f);
        return;
    }

    flis_selected = new_layer;
    flis_invalidate_composite();
    flis_gui_update();
    queue_redraw(REMAP_ALL);
}

void on_flis_remove_layer_clicked(GtkButton *btn,
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

void on_flis_duplicate_layer_clicked(GtkButton *btn,
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
    flis_invalidate_composite();
    flis_gui_update();
    queue_redraw(REMAP_ALL);
}

void on_flis_move_up_clicked(GtkButton *btn,
                                             gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected) return;

    gint idx = flis_layer_get_index(flis_selected);
    flis_layer_t *neighbor = (idx >= 0)
        ? (flis_layer_t *)g_slist_nth_data(com.uniq->layers, (guint)(idx + 1))
        : NULL;
    if (!neighbor) return;  /* already at top — nothing to do, no undo state */

    undo_save_flis_layer_reorder(flis_selected, neighbor, _("Move layer up"));
    flis_launch_layer_op(flis_selected, flis_move_up_hook, NULL,
                         _("Move layer up"), FALSE, NULL);
}

void on_flis_move_down_clicked(GtkButton *btn,
                                               gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected) return;

    gint idx = flis_layer_get_index(flis_selected);
    flis_layer_t *neighbor = (idx > 0)
        ? (flis_layer_t *)g_slist_nth_data(com.uniq->layers, (guint)(idx - 1))
        : NULL;
    if (!neighbor) return;  /* already at bottom */

    undo_save_flis_layer_reorder(flis_selected, neighbor, _("Move layer down"));
    flis_launch_layer_op(flis_selected, flis_move_down_hook, NULL,
                         _("Move layer down"), FALSE, NULL);
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

void on_flis_name_activate(GtkEntry *entry, gpointer data) {
    (void)entry; (void)data;
    flis_apply_name_entry();
}

gboolean on_flis_name_focus_out(GtkWidget *widget,
                                                GdkEventFocus *event,
                                                gpointer data) {
    (void)widget; (void)event; (void)data;
    flis_apply_name_entry();
    return FALSE; /* allow default focus handling */
}

void on_flis_blend_mode_changed(GtkComboBox *combo,
                                                gpointer     data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    gint idx = gtk_combo_box_get_active(combo);
    if (idx < 0 || idx >= N_BLEND_MODES) return;

    undo_save_flis_layer_props(flis_selected, _("Blend mode"));

    flis_blend_args_t *user = calloc(1, sizeof(flis_blend_args_t));
    if (!user) { PRINT_ALLOC_ERR; return; }
    user->mode = blend_mode_map[idx];

    flis_launch_layer_op(flis_selected, flis_blend_hook, user,
                         _("Blend mode"), FALSE, NULL);
}

/* Opacity debouncing state.
 * The GtkAdjustment fires value-changed on every drag tick; saving one undo
 * state per tick would flood the ring buffer.  Instead we snapshot the
 * pre-drag props on scale button-press and commit a single state on
 * button-release. */
static gint             flis_opacity_drag_id    = FLIS_UNDO_LAYER_NONE;
static flis_layer_props_t flis_opacity_drag_snapshot;

static gboolean on_flis_opacity_scale_press(GtkWidget *widget,
                                            GdkEventButton *event,
                                            gpointer data) {
    (void)widget; (void)event; (void)data;
    if (!flis_selected) return FALSE;

    /* Capture the full props at drag-start so we save one clean undo state */
    flis_opacity_drag_id              = flis_selected->item_id;
    flis_opacity_drag_snapshot.blend_mode = flis_selected->blend_mode;
    flis_opacity_drag_snapshot.opacity    = flis_selected->opacity;
    flis_opacity_drag_snapshot.visible    = flis_selected->visible;
    flis_opacity_drag_snapshot.locked     = flis_selected->locked;
    flis_opacity_drag_snapshot.has_tint   = flis_selected->has_tint;
    flis_opacity_drag_snapshot.tint       = flis_selected->layer_tint;
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

    /* Only commit if the layer is still the same one we started on */
    if (flis_selected && flis_selected->item_id == flis_opacity_drag_id)
        undo_save_flis_layer_props_snapshot(flis_opacity_drag_id,
                                            &flis_opacity_drag_snapshot,
                                            _("Opacity"));
    flis_opacity_drag_id = FLIS_UNDO_LAYER_NONE;
    return FALSE;
}

/* Opacity: connect to the GtkAdjustment value-changed signal in code so
 * it fires exactly once regardless of which widget (scale or spin) was
 * used.  Connected in flis_gui_init().
 * Note: undo state is NOT saved here — it is saved by the scale press/release
 * handlers above (one state per drag gesture).  Spinner changes save via
 * on_flis_opacity_spin_value_changed, connected in flis_gui_init(). */
static void on_flis_opacity_adj_changed(GtkAdjustment *adj, gpointer data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    gdouble pct = gtk_adjustment_get_value(adj);
    flis_layer_set_opacity(flis_selected, (gfloat)(pct / 100.0));
    flis_invalidate_composite();
    queue_redraw(REMAP_ALL);
}

/* Spinner-specific value-changed: fires on commit (Enter / focus-out), not
 * on every keypress, so saving a single undo state here is correct. */
static void on_flis_opacity_spin_commit(GtkSpinButton *spin, gpointer data) {
    (void)spin; (void)data;
    if (flis_updating || !flis_selected) return;
    /* The drag-snapshot path handles the scale; here we handle the spinner.
     * If a drag is in progress (drag_id is set), don't double-save. */
    if (flis_opacity_drag_id != FLIS_UNDO_LAYER_NONE) return;
    undo_save_flis_layer_props(flis_selected, _("Opacity"));
    /* Actual value change already applied by on_flis_opacity_adj_changed */
}

/* =========================================================================
 * Layer mask signal handlers
 * ========================================================================= */

void on_flis_mask_toggle_clicked(GtkButton *btn,
                                                 gpointer   data) {
    (void)btn; (void)data;
    if (!flis_selected) return;

    if (flis_selected->lmask) {
        /* Remove existing mask — confirm first */
        GtkWidget *dlg = gtk_message_dialog_new(
            GTK_WINDOW(lookup_widget("flis_layers_window")),
            GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
            GTK_MESSAGE_QUESTION, GTK_BUTTONS_OK_CANCEL,
            _("Remove the layer mask from \"%s\"?"),
            flis_selected->layer_name ? flis_selected->layer_name : _("Layer"));
        gint r = gtk_dialog_run(GTK_DIALOG(dlg));
        gtk_widget_destroy(dlg);
        if (r != GTK_RESPONSE_OK) return;

        undo_save_flis_lmask(flis_selected, _("Remove layer mask"));

        flis_lmask_toggle_args_t *user = calloc(1, sizeof(flis_lmask_toggle_args_t));
        if (!user) { PRINT_ALLOC_ERR; return; }
        user->destroy_fn = (destructor)flis_lmask_toggle_destroy;
        user->lm = NULL;  /* NULL = remove */

        flis_launch_layer_op(flis_selected, flis_lmask_toggle_hook, user,
                             _("Remove layer mask"), TRUE, NULL);
    } else {
        /* Add mask from file — file I/O and conversion happen on main thread
         * (fast, one-off), then the actual lmask attachment goes to the worker. */
        gchar *filename = flis_choose_fits_file(
            GTK_WINDOW(lookup_widget("flis_layers_window")),
            _("Choose greyscale mask file"));
        if (!filename) return;

        fits *mf = flis_load_layer_file(filename);
        g_free(filename);
        if (!mf) return;

        guint w = mf->rx, h = mf->ry;
        size_t npix = (size_t)w * h;

        layermask_t *lm = calloc(1, sizeof(layermask_t));
        if (!lm) { PRINT_ALLOC_ERR; clearfits(mf); free(mf); return; }

        lm->w = w; lm->h = h; lm->bitpix = 8;
        lm->data = malloc(npix);
        if (!lm->data) {
            PRINT_ALLOC_ERR; free(lm); clearfits(mf); free(mf); return;
        }

        uint8_t *dst = (uint8_t *)lm->data;
        if (mf->type == DATA_FLOAT && mf->fdata) {
            for (size_t i = 0; i < npix; i++) {
                float v = mf->fdata[i] * 255.f + 0.5f;
                dst[i] = v >= 255.f ? 255 : (v <= 0.f ? 0 : (uint8_t)v);
            }
        } else if (mf->type == DATA_USHORT && mf->data) {
            for (size_t i = 0; i < npix; i++)
                dst[i] = (uint8_t)(mf->data[i] >> 8);
        } else {
            memset(dst, 255, npix);
        }
        clearfits(mf); free(mf);

        undo_save_flis_lmask(flis_selected, _("Add layer mask"));

        flis_lmask_toggle_args_t *user = calloc(1, sizeof(flis_lmask_toggle_args_t));
        if (!user) { PRINT_ALLOC_ERR; layermask_free(lm); return; }
        user->destroy_fn = (destructor)flis_lmask_toggle_destroy;
        user->lm = lm;  /* non-NULL = add; ownership transfers to hook */

        flis_launch_layer_op(flis_selected, flis_lmask_toggle_hook, user,
                             _("Add layer mask"), TRUE, NULL);
    }
}

void on_flis_mask_move_clicked(GtkButton *btn,
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
                undo_save_flis_lmask_move(flis_selected, dest, _("Move layer mask"));

                flis_lmask_move_args_t *user = calloc(1, sizeof(flis_lmask_move_args_t));
                if (user) {
                    user->dest = dest;
                    flis_launch_layer_op(flis_selected, flis_lmask_move_hook, user,
                                         _("Move layer mask"), TRUE, NULL);
                }
            }
        }
    }
    gtk_widget_destroy(dlg);
}

/* =========================================================================
 * Mask view toggle (processing mask ↔ layer mask)
 * ========================================================================= */

void on_flis_mask_view_toggled(GtkToggleButton *btn,
                                                gpointer         data) {
    (void)data;
    if (flis_updating) return;

    GtkWidget *layer_radio = fw("flis_mask_view_layer_radio");
    gboolean show_layer = gtk_toggle_button_get_active(
                              GTK_TOGGLE_BUTTON(layer_radio));
    set_flis_show_layer_mask(show_layer);

    /* Update tab title and enable-check sensitivity */
    show_or_hide_mask_tab();
    /* Re-render the mask viewport and tint overlay */
    queue_redraw_mask();
}

/* Toggle layer mask active/inactive */
void on_flis_lmask_active_toggle_clicked(GtkButton *btn,
                                                          gpointer   data) {
    (void)btn;
    (void)data;
    if (flis_updating || !flis_selected || !flis_selected->lmask) return;

    undo_save_flis_layer_props(flis_selected,
                               flis_selected->lmask_active
                               ? _("Deactivate layer mask")
                               : _("Activate layer mask"));
    flis_selected->lmask_active = !flis_selected->lmask_active;
    flis_layer_touch_modified(flis_selected);

    /* Refresh UI and composite */
    flis_gui_update();
    redraw(REMAP_ALL);
}

/* flis_remove_selected_lmask — exported for use from callbacks.c */
void flis_remove_selected_lmask(void) {
    flis_layer_t *layer = flis_selected;
    if (!layer || !layer->lmask) return;

    GtkWidget *dlg = gtk_message_dialog_new(
        GTK_WINDOW(lookup_widget("flis_layers_window")),
        GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
        GTK_MESSAGE_QUESTION, GTK_BUTTONS_OK_CANCEL,
        _("Remove the layer mask from \"%s\"?"),
        layer->layer_name ? layer->layer_name : _("Layer"));
    gint r = gtk_dialog_run(GTK_DIALOG(dlg));
    gtk_widget_destroy(dlg);
    if (r != GTK_RESPONSE_OK) return;

    undo_save_flis_lmask(layer, _("Remove layer mask"));

    flis_lmask_toggle_args_t *user = calloc(1, sizeof(flis_lmask_toggle_args_t));
    if (!user) { PRINT_ALLOC_ERR; return; }
    user->destroy_fn = (destructor)flis_lmask_toggle_destroy;
    user->lm = NULL;  /* NULL = remove */

    flis_launch_layer_op(layer, flis_lmask_toggle_hook, user,
                         _("Remove layer mask"), TRUE, NULL);
}

/* flis_get_selected_layer_id — exported, used by masks_gui.c */
gint flis_get_selected_layer_id(void) {
    if (!flis_selected || !is_current_image_flis()) return 0;
    return flis_selected->item_id;
}

/* flis_populate_layer_combo — exported, used by masks_gui.c.
 * When filter_by_canvas_size is TRUE, only layers whose pixel dimensions
 * match gfit (the composited canvas) are included — masks are produced from
 * the canvas so target layers must be the same size.  Pass FALSE when the
 * mask source dimensions are not yet known (e.g. file-mode before a file
 * has been selected). */
void flis_populate_layer_combo(GtkComboBoxText *combo, gboolean filter_by_canvas_size) {
    if (!combo) return;
    gtk_combo_box_text_remove_all(combo);
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return;
    guint mask_rx = gfit->rx;
    guint mask_ry = gfit->ry;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->layer_name || !lay->fit) continue;
        if (filter_by_canvas_size &&
            ((guint)lay->fit->rx != mask_rx || (guint)lay->fit->ry != mask_ry))
            continue;
        gtk_combo_box_text_append_text(combo, lay->layer_name);
    }
}

/* flis_combo_select_active_layer — exported, used by masks_gui.c.
 * Selects the entry in @combo whose text matches the active layer's name.
 * No-op if the combo is empty or the active layer is not present in it
 * (e.g. filtered out by canvas-size). */
void flis_combo_select_active_layer(GtkComboBoxText *combo) {
    if (!combo) return;
    flis_layer_t *active = flis_active_layer();
    if (!active || !active->layer_name) {
        gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
        return;
    }
    GtkTreeModel *model = gtk_combo_box_get_model(GTK_COMBO_BOX(combo));
    if (!model) return;
    GtkTreeIter iter;
    gint idx = 0;
    if (gtk_tree_model_get_iter_first(model, &iter)) {
        do {
            gchar *text = NULL;
            gtk_tree_model_get(model, &iter, 0, &text, -1);
            gboolean match = text && g_strcmp0(text, active->layer_name) == 0;
            g_free(text);
            if (match) {
                gtk_combo_box_set_active(GTK_COMBO_BOX(combo), idx);
                return;
            }
            idx++;
        } while (gtk_tree_model_iter_next(model, &iter));
    }
    /* Active layer not in combo (filtered out): fall back to first entry */
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
}

/* =========================================================================
 * Layer context menu
 * ========================================================================= */

/* Right-click on the layer list: select the row under the cursor and pop
 * the context menu. */
static gboolean on_flis_layer_list_button_press(GtkWidget      *widget,
                                                 GdkEventButton *event,
                                                 gpointer        data) {
    (void)data;
    if (event->type != GDK_BUTTON_PRESS || event->button != 3)
        return GDK_EVENT_PROPAGATE;

    GtkListBox    *list = GTK_LIST_BOX(widget);
    GtkListBoxRow *row  = gtk_list_box_get_row_at_y(list, (gint)event->y);
    if (row)
        gtk_list_box_select_row(list, row);

    GtkMenu *menu = GTK_MENU(lookup_widget("flis_layer_context_menu"));
    if (menu && flis_selected) {
        /* Merge Down is only available when a layer exists below the selection */
        GtkWidget *merge_item = lookup_widget("flis_merge_down_item");
        if (merge_item) {
            gboolean can_merge = FALSE;
            if (flis_selected && com.uniq) {
                for (GSList *l = com.uniq->layers; l; l = l->next) {
                    flis_layer_t *lay = (flis_layer_t *)l->data;
                    if (lay && lay != flis_selected &&
                        lay->layer_order < flis_selected->layer_order) {
                        can_merge = TRUE;
                        break;
                    }
                }
            }
            gtk_widget_set_sensitive(merge_item, can_merge);
        }
        /* Sync the "Move Layer" toggle button to the current mouse mode */
        GtkToggleButton *drag_btn = GTK_TOGGLE_BUTTON(lookup_widget("flis_drag_toggle_btn"));
        if (drag_btn) {
            g_signal_handlers_block_by_func(drag_btn, on_flis_move_layer_toggled, NULL);
            gtk_toggle_button_set_active(drag_btn,
                mouse_status == MOUSE_ACTION_FLIS_DRAG_LAYER);
            g_signal_handlers_unblock_by_func(drag_btn, on_flis_move_layer_toggled, NULL);
        }
        gtk_menu_popup_at_pointer(menu, (GdkEvent *)event);
    }

    return GDK_EVENT_STOP;
}

void on_flis_merge_down_activate(GtkMenuItem *item, gpointer data) {
    (void)item; (void)data;
    if (!flis_selected || !is_current_image_flis()) return;

    if (!siril_confirm_dialog(_("Merge Down"),
            _("Merge the current layer onto the one below it?\n\n"
              "This operation cannot be undone."),
            _("Merge Down")))
        return;

    if (flis_merge_down_layer(flis_selected)) {
        siril_message_dialog(GTK_MESSAGE_ERROR, _("Merge Down"),
                             _("Merge Down failed."));
        return;
    }

    flis_gui_update();
    show_or_hide_mask_tab();
    redraw(REMAP_ALL);
}

void on_flis_background_neutralise_activate(GtkMenuItem *item, gpointer data) {
    (void)item; (void)data;
    if (!is_current_image_flis()) return;

    if (flis_background_neutralise()) {
        siril_message_dialog(GTK_MESSAGE_ERROR, _("Background Neutralise"),
                             _("Background neutralisation failed."));
        return;
    }

    flis_gui_update();
    redraw(REMAP_ALL);
}

void on_flis_flatten_activate(GtkMenuItem *item, gpointer data) {
    (void)item; (void)data;
    if (!is_current_image_flis()) return;

    if (flis_layer_count() <= 1) {
        siril_message_dialog(GTK_MESSAGE_INFO, _("Flatten Image"),
                             _("The image already has only one layer."));
        return;
    }

    if (!siril_confirm_dialog(_("Flatten Image"),
            _("Flatten all visible layers into a single layer?\n\n"
              "All masks will be removed and this operation cannot be undone."),
            _("Flatten")))
        return;

    if (flis_flatten_all()) {
        siril_message_dialog(GTK_MESSAGE_ERROR, _("Flatten Image"),
                             _("Flatten failed."));
        return;
    }

    flis_gui_update();
    show_or_hide_mask_tab();
    redraw(REMAP_ALL);
}

void on_flis_move_layer_toggled(GtkToggleButton *btn, gpointer data) {
    (void)data;
    gboolean active = gtk_toggle_button_get_active(btn);
    if (active) {
        mouse_status = MOUSE_ACTION_FLIS_DRAG_LAYER;
        set_cursor("grab");
    } else {
        mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
        gui.flis_layer_dragging = FALSE;
        set_cursor("crosshair");
    }
}

/*
 * Export the currently selected FLIS layer to a plain FITS file.
 *
 * The exported file gets:
 *   • pixel data from the layer itself (shared read-only via shallow copy)
 *   • FITS keywords (exposure, date, telescope, etc.) from the base layer,
 *     which is the authoritative metadata source for the whole FLIS stack
 *   • the base layer's ICC profile (borrowed; non-base layers carry no profile
 *     per the FLIS invariant)
 */
void on_flis_export_layer_activate(GtkMenuItem *item,
                                                    gpointer     data) {
    (void)item; (void)data;
    if (!flis_selected || !flis_selected->fit) return;

    flis_layer_t *layer  = flis_selected;
    GtkWindow    *parent = GTK_WINDOW(lookup_widget("flis_layers_window"));

    gchar *suggested = g_strdup_printf("%s.fits",
        (layer->layer_name && layer->layer_name[0])
            ? layer->layer_name : _("layer"));
    gchar *path = flis_choose_fits_save(parent,
                                        _("Export Layer as FITS"),
                                        suggested);
    g_free(suggested);
    if (!path) return;

    /* Shallow-copy the layer's fits struct so fdata/data/naxes/bitpix etc.
     * are correct.  We then overlay metadata from the base layer. */
    fits export_fit = *layer->fit;

    /* savefits creates a new fptr; null the borrowed one so there is no
     * risk of accidentally closing the layer's open file handle. */
    export_fit.fptr  = NULL;
    export_fit.stats = NULL;   /* belongs to layer->fit, not needed for save */

    /* Null out metadata fields that will be replaced by copy_fits_metadata()
     * so there is no confusion about ownership. */
    export_fit.keywords.date_obs = NULL;
    export_fit.keywords.wcslib   = NULL;
    export_fit.unknown_keys      = NULL;
    export_fit.icc_profile       = NULL;
    export_fit.color_managed     = FALSE;

    /* Copy FITS keywords from the base (profiled) layer — it carries the
     * authoritative header for the whole FLIS stack. */
    fits *base_fit = flis_get_profiled_fit();
    if (base_fit)
        copy_fits_metadata(base_fit, &export_fit);

    /* Borrow the base layer's ICC profile.  savefits only reads it (to
     * embed the serialised bytes), so this is safe without duplication. */
    if (base_fit && base_fit->icc_profile) {
        export_fit.icc_profile   = base_fit->icc_profile;
        export_fit.color_managed = base_fit->color_managed;
    }

    int ret = savefits(path, &export_fit);

    /* Free only what copy_fits_metadata() allocated.  Never touch fdata/data
     * — they belong to layer->fit.  icc_profile was borrowed, not duplicated. */
    if (export_fit.keywords.date_obs)
        g_date_time_unref(export_fit.keywords.date_obs);
    if (export_fit.keywords.wcslib)
        free_wcs(&export_fit);
    if (export_fit.unknown_keys)
        g_free(export_fit.unknown_keys);

    if (ret != 0)
        siril_log_color_message(
            _("FLIS: failed to export layer \"%s\"\n"), "red",
            layer->layer_name ? layer->layer_name : _("layer"));

    g_free(path);
}

/* =========================================================================
 * Layer registration
 * ========================================================================= */

/*
 * Build and run a small modal dialog asking for the registration method.
 * Framing is always FRAMING_MAX (maximum bounding box), so there is no
 * framing choice exposed to the user.
 *
 * Returns GTK_RESPONSE_ACCEPT if the user clicked Register, otherwise cancel.
 * Fills *method_out on accept.
 */
static gint flis_register_dialog(GtkWindow              *parent,
                                  struct registration_method **method_out) {
    /* Build the three methods we expose (same set as compositing.c). */
    struct registration_method *methods[3];
    methods[0] = new_reg_method(
        _("Deep Sky (global star alignment)"),
        &register_multi_step_global, REQUIRES_NO_SELECTION,
        REGTYPE_DEEPSKY);
    methods[1] = new_reg_method(
        _("Planetary (DFT pattern alignment)"),
        &register_shift_dft, REQUIRES_SQUARED_SELECTION,
        REGTYPE_PLANETARY);
    methods[2] = new_reg_method(
        _("Planetary (KOMBAT pattern alignment)"),
        &register_kombat, REQUIRES_ANY_SELECTION,
        REGTYPE_PLANETARY);

    GtkWidget *dlg = gtk_dialog_new_with_buttons(
        _("Register FLIS Layers"), parent,
        GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
        _("Cancel"),   GTK_RESPONSE_CANCEL,
        _("Register"), GTK_RESPONSE_ACCEPT,
        NULL);
    gtk_dialog_set_default_response(GTK_DIALOG(dlg), GTK_RESPONSE_ACCEPT);

    GtkWidget *content = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    GtkWidget *grid    = gtk_grid_new();
    gtk_grid_set_row_spacing(GTK_GRID(grid), 6);
    gtk_grid_set_column_spacing(GTK_GRID(grid), 12);
    gtk_widget_set_margin_start(grid, 12);
    gtk_widget_set_margin_end(grid, 12);
    gtk_widget_set_margin_top(grid, 12);
    gtk_widget_set_margin_bottom(grid, 12);
    gtk_container_add(GTK_CONTAINER(content), grid);

    /* Method row */
    GtkWidget *lbl_method = gtk_label_new(_("Method"));
    gtk_label_set_xalign(GTK_LABEL(lbl_method), 0.0);
    GtkWidget *method_combo = gtk_combo_box_text_new();
    for (int i = 0; i < 3; i++)
        gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(method_combo),
                                       methods[i]->name);
    gtk_combo_box_set_active(GTK_COMBO_BOX(method_combo), 0);
    gtk_widget_set_hexpand(method_combo, TRUE);
    gtk_grid_attach(GTK_GRID(grid), lbl_method,   0, 0, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), method_combo, 1, 0, 1, 1);

    gtk_widget_show_all(dlg);
    gint response = gtk_dialog_run(GTK_DIALOG(dlg));

    if (response == GTK_RESPONSE_ACCEPT) {
        int mi = gtk_combo_box_get_active(GTK_COMBO_BOX(method_combo));
        *method_out = methods[mi];
        /* free the two methods we won't use */
        for (int i = 0; i < 3; i++)
            if (i != mi) free(methods[i]);
    } else {
        for (int i = 0; i < 3; i++) free(methods[i]);
    }

    gtk_widget_destroy(dlg);
    return response;
}

/*
 * Register all FLIS layers against the currently selected layer (reference).
 *
 * Always uses FRAMING_MAX so each layer is transformed into its own minimum
 * bounding box.  After apply_reg the per-layer bounding-box offsets stored in
 * regparam[i].H are used to set position_x/y relative to the canvas layer.
 */
void on_flis_register_layers_activate(GtkMenuItem *item,
                                                       gpointer     data) {
    (void)item; (void)data;

    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return;

    gint n_layers = flis_layer_count();
    if (n_layers < 2) {
        siril_message_dialog(GTK_MESSAGE_WARNING, _("Register Layers"),
            _("At least two layers are required for registration."));
        return;
    }

    /* Show method dialog (framing is always FRAMING_MAX). */
    GtkWindow *parent = GTK_WINDOW(lookup_widget("flis_layers_window"));
    struct registration_method *method = NULL;
    if (flis_register_dialog(parent, &method) != GTK_RESPONSE_ACCEPT)
        return;

    /* The canvas (base) layer is the first layer in the stack.
     * The registration reference is the currently selected layer, or the
     * canvas layer if nothing is selected. */
    flis_layer_t *canvas_lay = (flis_layer_t *)com.uniq->layers->data;
    flis_layer_t *ref_lay    = flis_selected ? flis_selected : canvas_lay;

    /* Build an internal sequence.
     * Slot 0 = registration reference; remaining slots follow GSList order,
     * skipping the reference.  Track which slot the canvas layer lands in. */
    sequence *seq = create_internal_sequence(n_layers);
    seq->bitpix = ref_lay->fit->bitpix;
    seq->rx     = ref_lay->fit->rx;
    seq->ry     = ref_lay->fit->ry;

    const int ref_seq_idx = 0;
    int canvas_seq_idx = (ref_lay == canvas_lay) ? 0 : -1;

    /* Reference at slot 0 */
    internal_sequence_set(seq, 0, ref_lay->fit);
    seq->imgparam[0].rx = ref_lay->fit->rx;
    seq->imgparam[0].ry = ref_lay->fit->ry;

    int i = 1;
    gboolean is_variable = FALSE;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->fit) continue;
        if (lay == ref_lay) continue;
        if (lay == canvas_lay) canvas_seq_idx = i;
        if (lay->fit->rx != ref_lay->fit->rx || lay->fit->ry != ref_lay->fit->ry)
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
    regargs.interpolation        = OPENCV_LANCZOS4;
    regargs.output_scale         = 1.f;
    regargs.clamp                = TRUE;
    regargs.framing              = FRAMING_MAX;
    regargs.max_stars_candidates = MAX_STARS_FITTED;
    regargs.two_pass             = TRUE;
    regargs.percent_moved        = 0.50f;
    regargs.type = (method->method_ptr == register_shift_dft ||
                    method->method_ptr == register_kombat)
                   ? SHIFT_TRANSFORMATION : HOMOGRAPHY_TRANSFORMATION;
    get_the_registration_area(&regargs, method);

    set_cursor_waiting(TRUE);
    set_progress_bar_data(_("Registering FLIS layers…"), PROGRESS_RESET);
    com.run_thread = TRUE;

    /* Step 1: compute registration transforms. */
    int ret1 = method->method_ptr(&regargs);
    free(regargs.imgparam); regargs.imgparam = NULL;
    free(regargs.regparam); regargs.regparam = NULL;
    free(method);

    if (ret1) {
        free_sequence(seq, TRUE);
        com.run_thread = FALSE;
        set_cursor_waiting(FALSE);
        set_progress_bar_data(_("Registration failed."), PROGRESS_DONE);
        return;
    }

    /* Step 2: apply transforms with FRAMING_MAX.  Each layer is resampled into
     * its own minimum bounding box; regparam[i].H becomes the Hshift storing
     * the (xmin, ymin) offset of that bounding box in the global mosaic frame
     * (OpenCV convention: y increases downward). */
    int ret2 = register_apply_reg(&regargs);

    if (!ret2 && regargs.regparam) {
        /* Normalise offsets relative to the canvas layer so that the canvas
         * stays at position (0, 0) and all other layers are placed correctly
         * relative to it. */
        const double cx       = regargs.regparam[canvas_seq_idx].H.h02;
        const double cy       = regargs.regparam[canvas_seq_idx].H.h12;

        /* Canvas layer always sits at origin */
        canvas_lay->position_x = 0;
        canvas_lay->position_y = 0;

        /* All other layers: compute position from normalised Hshift offset.
         * position_y = canvas_H - dy_rel - layer_H  (FITS bottom-up convention)
         * where dy_rel = h12_i - cy  (y from top of canvas, OpenCV y-down) */
        int k = 1;  /* seq index counter for non-ref layers */
        for (GSList *l = com.uniq->layers; l; l = l->next) {
            flis_layer_t *lay = (flis_layer_t *)l->data;
            if (!lay || !lay->fit) continue;
            if (lay == canvas_lay) {
                if (lay != ref_lay) k++;
                continue;
            }
            int seq_idx = (lay == ref_lay) ? 0 : k++;
            double dx = regargs.regparam[seq_idx].H.h02 - cx;
            double dy = regargs.regparam[seq_idx].H.h12 - cy;
            lay->position_x = (gint) round(dx);
            lay->position_y = (gint) round(dy);
        }
    }

    free(regargs.imgparam); regargs.imgparam = NULL;
    free(regargs.regparam); regargs.regparam = NULL;

    free_sequence(seq, TRUE);
    com.run_thread = FALSE;
    set_cursor_waiting(FALSE);

    if (ret2) {
        set_progress_bar_data(_("Registration failed."), PROGRESS_DONE);
        return;
    }

    set_progress_bar_data(_("Layer registration complete."), PROGRESS_DONE);

    /* Layers have been transformed in-place; invalidate the composite. */
    flis_invalidate_composite();
    queue_redraw(REMAP_ALL);
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

void on_flis_tint_toggled(GtkToggleButton *btn,
                                          gpointer         data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    undo_save_flis_layer_props(flis_selected, _("Tint"));

    gboolean enable = gtk_toggle_button_get_active(btn);
    gtk_widget_set_sensitive(fw("flis_tint_color_btn"), enable);

    flis_tint_toggle_args_t *user = calloc(1, sizeof(flis_tint_toggle_args_t));
    if (!user) { PRINT_ALLOC_ERR; return; }
    user->enable  = enable;
    user->was_rgb = flis_composite_will_be_rgb();
    if (enable) {
        GdkRGBA colour;
        gtk_color_chooser_get_rgba(
            GTK_COLOR_CHOOSER(fw("flis_tint_color_btn")), &colour);
        user->r = colour.red;
        user->g = colour.green;
        user->b = colour.blue;
    }

    flis_launch_layer_op(flis_selected, flis_tint_toggle_hook, user,
                         _("Tint"), FALSE, end_flis_tint_toggle);
}

void on_flis_tint_color_set(GtkColorButton *btn,
                                            gpointer        data) {
    (void)data;
    if (flis_updating || !flis_selected) return;

    GdkRGBA colour;
    gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(btn), &colour);
    undo_save_flis_layer_props(flis_selected, _("Tint colour"));

    flis_tint_color_args_t *user = calloc(1, sizeof(flis_tint_color_args_t));
    if (!user) { PRINT_ALLOC_ERR; return; }
    user->r = colour.red;
    user->g = colour.green;
    user->b = colour.blue;

    flis_launch_layer_op(flis_selected, flis_tint_color_hook, user,
                         _("Tint colour"), FALSE, NULL);
}

/* =========================================================================
 * Window signal handler
 * ========================================================================= */

gboolean on_flis_layers_window_delete_event(
    GtkWidget *widget, GdkEventAny *event, gpointer data) {
    (void)widget; (void)event; (void)data;
    /* Hide rather than destroy so state is preserved */
    gtk_widget_hide(lookup_widget("flis_layers_window"));
    return TRUE; /* prevent actual destruction */
}

/* Fires every time the dialog is shown (gtk_widget_show / gtk_window_present).
 * This is the right place for setup that must happen each time the panel
 * becomes visible: setting the transient parent, keep-above, and rebuilding
 * the layer list so it is always current. */
void on_flis_layers_window_show(GtkWidget *widget,
                                                gpointer   data) {
    (void)data;
    GtkWindow *win     = GTK_WINDOW(widget);
    GtkWindow *main_win = GTK_WINDOW(lookup_widget("control_window"));

    /* Ensure the dialog is always parented to the main window.  This scopes
     * keep-above to the application rather than the whole desktop. */
    if (main_win)
        gtk_window_set_transient_for(win, main_win);

    /* keep-above must be set in code; GTK's UI loader rejects it as a
     * static property on GtkDialog. */
    gtk_window_set_keep_above(win, TRUE);

    /* Rebuild the layer list so it reflects the current image state. */
    flis_layers_list_rebuild();
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

    /* Opacity scale: capture pre-drag state on press, commit one undo state
     * on release — avoids flooding history with every drag tick. */
    GtkWidget *scale = lookup_widget("flis_opacity_scale");
    if (scale) {
        g_signal_connect(scale, "button-press-event",
                         G_CALLBACK(on_flis_opacity_scale_press), NULL);
        g_signal_connect(scale, "button-release-event",
                         G_CALLBACK(on_flis_opacity_scale_release), NULL);
    }

    /* Spinner commit: fires on Enter / focus-out, not on every keypress */
    GtkWidget *spin = lookup_widget("flis_opacity_spin");
    if (spin)
        g_signal_connect(spin, "value-changed",
                         G_CALLBACK(on_flis_opacity_spin_commit), NULL);

    /* Right-click on a layer row pops the context menu */
    GtkWidget *list = lookup_widget("flis_layer_list");
    if (list)
        g_signal_connect(list, "button-press-event",
                         G_CALLBACK(on_flis_layer_list_button_press), NULL);
}

void flis_gui_open(void) {
    GtkWidget *win = lookup_widget("flis_layers_window");
    if (!win) {
        siril_log_color_message(
            _("FLIS GUI: layers window not found — "
              "is flis_layers.ui included in the UI bundle?\n"), "red");
        return;
    }
    /* The show signal handler (on_flis_layers_window_show) takes care of
     * setting the transient parent, keep-above, and rebuilding the list. */
    gtk_widget_show(win);
    gtk_window_present(GTK_WINDOW(win));
}
