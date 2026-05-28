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

#include <math.h>

#include "core/initfile.h"
#include "core/siril_log.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/mouse_action_functions.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/utils.h"

static GtkScrolledWindow *mouse_scrolled_window = NULL;
static GtkScrolledWindow *scroll_scrolled_window = NULL;
static GtkSpinButton *mouse_speed_spin = NULL;
static GtkWidget *mouse_test_drawingarea = NULL;
static GtkWidget *mouse_test_frame_widget = NULL;
static GtkLabel *mouse_test_label = NULL;

/* Live-edited GTK4 model and view (per dialog instance). */
static GListStore       *mouse_action_store    = NULL;
static GtkColumnView    *mouse_action_columnview = NULL;
static GListStore       *scroll_action_store   = NULL;
static GtkColumnView    *scroll_action_columnview = NULL;

static void install_mouse_test_controllers(GtkWidget *area);

static void mouse_action_prefs_init_statics(void) {
	if (mouse_scrolled_window) return;
	mouse_scrolled_window = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "mouse_treeview_scrolled_window"));
	scroll_scrolled_window = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "scroll_treeview_scrolled_window"));
	mouse_speed_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mouse_speed_limit"));
	mouse_test_drawingarea = GTK_WIDGET(gtk_builder_get_object(gui.builder, "mouse_test_drawingarea"));
	mouse_test_frame_widget = GTK_WIDGET(gtk_builder_get_object(gui.builder, "mouse_test_frame"));
	mouse_test_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "mouse_test_area_label"));
}

/* Elements that require updating when adding new functions */

const mouse_function_metadata* mouse_metadata_array[] = {
	&main_action,
	&drag_action,
	&measure_action,
	&open_if_unloaded_action,
	&zoom_action,
	&photometry_box_action,
	&second_action,
	&save_on_click_action,
	&save_as_on_click_action,
	&snapshot_on_click_action,
	&undo_on_click_action,
	&redo_on_click_action,
	&findstar_on_click_action,
	&platesolve_on_click_action,
	&spcc_on_click_action,
	NULL
};

const scroll_function_metadata* scroll_metadata_array[] = {
	&scroll_zooms_action,
	&scroll_traverses_seq_action,
	&scroll_changes_tab_action,
	NULL
};

/* End of elements that require updating when adding new functions */

static const mouse_function_metadata *map_ref_to_metadata(const mouse_function_ref reference) {
	guint i = 0;
	while (mouse_metadata_array[i]) {
		if (reference == mouse_metadata_array[i]->reference)
			break;
		i++;
	}
	return (mouse_metadata_array[i]) ? mouse_metadata_array[i] : &null_action;
}

static const mouse_function_metadata *map_name_to_metadata(const gchar *name) {
	guint i = 0;
	while (mouse_metadata_array[i]) {
		if (!strcmp(name, mouse_metadata_array[i]->name))
			break;
		i++;
	}
	return (mouse_metadata_array[i]) ? mouse_metadata_array[i] : &null_action;
}

static const scroll_function_metadata *map_scroll_ref_to_metadata(const scroll_function_ref reference) {
	guint i = 0;
	while (scroll_metadata_array[i]) {
		if (reference == scroll_metadata_array[i]->reference)
			break;
		i++;
	}
	return (scroll_metadata_array[i]) ? scroll_metadata_array[i] : &scroll_null_action;
}

static const scroll_function_metadata *map_scroll_name_to_metadata(const gchar *name) {
	guint i = 0;
	while (scroll_metadata_array[i]) {
		if (!strcmp(name, scroll_metadata_array[i]->name))
			break;
		i++;
	}
	return (scroll_metadata_array[i]) ? scroll_metadata_array[i] : &scroll_null_action;
}

enum {
	NONE = 0,
	SHIFT = 1,
	CTRL = 2,
	CTRLSHIFT = 3
};

void on_config_mouse_buttons_clicked(GtkButton *dialog, gpointer user_data) {
	if (siril_confirm_dialog(_("Mouse Configuration"), _("Warning: this dialog allows you to reconfigure mouse behaviour so that it may no longer match the documentation. If you forget what settings you have changed, the dialog can also be used to view current mouse action assignments and also provides a button to reset them to default settings."), _("Proceed"))) {
		siril_open_dialog("mouse_actions_dialog");
	}
}

/* ====================================================================
 * GTK4 list models: per-row GObject types backing the GListStore
 * ==================================================================== */

#define SIRIL_TYPE_MOUSE_ACTION_ITEM   (siril_mouse_action_item_get_type())
G_DECLARE_FINAL_TYPE(SirilMouseActionItem, siril_mouse_action_item, SIRIL, MOUSE_ACTION_ITEM, GObject)

struct _SirilMouseActionItem {
	GObject parent_instance;
	gchar  *name;       /* mouse function name */
	guint   button;     /* 1..15 */
	gchar  *action;     /* CLICK_TEXT or DOUBLE_CLICK_TEXT */
	gchar  *modifier;   /* CTRL/SHIFT/CTRL_SHIFT/NO_MODIFIER */
	gchar  *tooltip;
	guint   reference;  /* MOUSE_REF_* */
};

G_DEFINE_TYPE(SirilMouseActionItem, siril_mouse_action_item, G_TYPE_OBJECT)

static void siril_mouse_action_item_finalize(GObject *obj) {
	SirilMouseActionItem *self = SIRIL_MOUSE_ACTION_ITEM(obj);
	g_clear_pointer(&self->name, g_free);
	g_clear_pointer(&self->action, g_free);
	g_clear_pointer(&self->modifier, g_free);
	g_clear_pointer(&self->tooltip, g_free);
	G_OBJECT_CLASS(siril_mouse_action_item_parent_class)->finalize(obj);
}

static void siril_mouse_action_item_class_init(SirilMouseActionItemClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_mouse_action_item_finalize;
}
static void siril_mouse_action_item_init(SirilMouseActionItem *self) { (void)self; }

static SirilMouseActionItem *siril_mouse_action_item_new(const gchar *name, guint button,
		const gchar *action, const gchar *modifier, const gchar *tooltip, guint reference) {
	SirilMouseActionItem *self = g_object_new(SIRIL_TYPE_MOUSE_ACTION_ITEM, NULL);
	self->name      = g_strdup(name ? name : NEW_ENTRY_TEXT);
	self->button    = button;
	self->action    = g_strdup(action ? action : CLICK_TEXT);
	self->modifier  = g_strdup(modifier ? modifier : NO_MODIFIER_TEXT);
	self->tooltip   = g_strdup(tooltip ? tooltip : "");
	self->reference = reference;
	return self;
}

#define SIRIL_TYPE_SCROLL_ACTION_ITEM   (siril_scroll_action_item_get_type())
G_DECLARE_FINAL_TYPE(SirilScrollActionItem, siril_scroll_action_item, SIRIL, SCROLL_ACTION_ITEM, GObject)

struct _SirilScrollActionItem {
	GObject parent_instance;
	gchar  *name;
	gchar  *direction;
	gchar  *modifier;
	gchar  *tooltip;
	guint   reference;
};

G_DEFINE_TYPE(SirilScrollActionItem, siril_scroll_action_item, G_TYPE_OBJECT)

static void siril_scroll_action_item_finalize(GObject *obj) {
	SirilScrollActionItem *self = SIRIL_SCROLL_ACTION_ITEM(obj);
	g_clear_pointer(&self->name, g_free);
	g_clear_pointer(&self->direction, g_free);
	g_clear_pointer(&self->modifier, g_free);
	g_clear_pointer(&self->tooltip, g_free);
	G_OBJECT_CLASS(siril_scroll_action_item_parent_class)->finalize(obj);
}

static void siril_scroll_action_item_class_init(SirilScrollActionItemClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_scroll_action_item_finalize;
}
static void siril_scroll_action_item_init(SirilScrollActionItem *self) { (void)self; }

static SirilScrollActionItem *siril_scroll_action_item_new(const gchar *name, const gchar *direction,
		const gchar *modifier, const gchar *tooltip, guint reference) {
	SirilScrollActionItem *self = g_object_new(SIRIL_TYPE_SCROLL_ACTION_ITEM, NULL);
	self->name      = g_strdup(name ? name : NEW_ENTRY_TEXT);
	self->direction = g_strdup(direction ? direction : SCROLL_VERTICAL_TEXT);
	self->modifier  = g_strdup(modifier ? modifier : NO_MODIFIER_TEXT);
	self->tooltip   = g_strdup(tooltip ? tooltip : "");
	self->reference = reference;
	return self;
}

/* ====================================================================
 * GtkSignalListItemFactory helpers: inline-editing cells.
 * Each cell hosts a permanent editor widget (GtkDropDown / GtkSpinButton)
 * bound to a property of the row item.  Edits propagate live to the
 * underlying GObject; "Apply" then reads the GListStore directly.
 * ==================================================================== */

static gint find_string_position(const char *const *table, gsize n, const gchar *target) {
	if (!target) return 0;
	for (gsize i = 0; i < n; i++) {
		if (table[i] && !g_strcmp0(table[i], target))
			return (gint)i;
	}
	return 0;
}

/* --- mouse: NAME (combo from mouse_metadata_array) ----------------- */

static GtkStringList *mouse_name_string_list(void) {
	GtkStringList *sl = gtk_string_list_new(NULL);
	for (guint i = 0; mouse_metadata_array[i]; i++)
		gtk_string_list_append(sl, mouse_metadata_array[i]->name);
	return sl;
}

static void mouse_name_dd_notify(GObject *dd, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkListItem *li = g_object_get_data(dd, "siril-listitem");
	if (!li) return;
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	if (!row) return;
	GObject *sel = gtk_drop_down_get_selected_item(GTK_DROP_DOWN(dd));
	if (!sel || !GTK_IS_STRING_OBJECT(sel)) return;
	const gchar *name = gtk_string_object_get_string(GTK_STRING_OBJECT(sel));
	const mouse_function_metadata *md = map_name_to_metadata(name);
	g_free(row->name);    row->name    = g_strdup(name);
	g_free(row->tooltip); row->tooltip = g_strdup(md->tooltip ? md->tooltip : "");
	row->reference = md->reference;
}

static void name_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_drop_down_new(G_LIST_MODEL(mouse_name_string_list()), NULL));
	gtk_list_item_set_child(li, GTK_WIDGET(dd));
	g_signal_connect(dd, "notify::selected", G_CALLBACK(mouse_name_dd_notify), NULL);
}

static void name_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", li);
	guint pos = 0;
	for (guint i = 0; mouse_metadata_array[i]; i++) {
		if (!g_strcmp0(row->name, mouse_metadata_array[i]->name)) { pos = i; break; }
	}
	g_signal_handlers_block_by_func(dd, mouse_name_dd_notify, NULL);
	gtk_drop_down_set_selected(dd, pos);
	g_signal_handlers_unblock_by_func(dd, mouse_name_dd_notify, NULL);
	gtk_widget_set_tooltip_text(GTK_WIDGET(dd), row->tooltip);
}

static void name_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", NULL);
}

/* --- mouse: BUTTON (spin 1..15) ------------------------------------ */

static void button_spin_changed(GtkSpinButton *btn, gpointer user_data) {
	(void)user_data;
	GtkListItem *li = g_object_get_data(G_OBJECT(btn), "siril-listitem");
	if (!li) return;
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	if (!row) return;
	row->button = (guint)gtk_spin_button_get_value_as_int(btn);
}

static void button_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkSpinButton *btn = GTK_SPIN_BUTTON(gtk_spin_button_new_with_range(1, 15, 1));
	gtk_spin_button_set_digits(btn, 0);
	gtk_list_item_set_child(li, GTK_WIDGET(btn));
	g_signal_connect(btn, "value-changed", G_CALLBACK(button_spin_changed), NULL);
}

static void button_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkSpinButton *btn = GTK_SPIN_BUTTON(gtk_list_item_get_child(li));
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", li);
	g_signal_handlers_block_by_func(btn, button_spin_changed, NULL);
	gtk_spin_button_set_value(btn, (gdouble)row->button);
	g_signal_handlers_unblock_by_func(btn, button_spin_changed, NULL);
}

static void button_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkSpinButton *btn = GTK_SPIN_BUTTON(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(btn), "siril-listitem", NULL);
}

/* --- mouse: ACTION (combo: Click / Double click) ------------------- */

static const char *const mouse_action_choices[] = { CLICK_TEXT, DOUBLE_CLICK_TEXT };
static const gsize mouse_action_choice_count = G_N_ELEMENTS(mouse_action_choices);

static GtkStringList *mouse_action_string_list(void) {
	GtkStringList *sl = gtk_string_list_new(NULL);
	for (gsize i = 0; i < mouse_action_choice_count; i++)
		gtk_string_list_append(sl, mouse_action_choices[i]);
	return sl;
}

static void action_dd_notify(GObject *dd, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkListItem *li = g_object_get_data(dd, "siril-listitem");
	if (!li) return;
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	if (!row) return;
	GObject *sel = gtk_drop_down_get_selected_item(GTK_DROP_DOWN(dd));
	if (!sel || !GTK_IS_STRING_OBJECT(sel)) return;
	g_free(row->action);
	row->action = g_strdup(gtk_string_object_get_string(GTK_STRING_OBJECT(sel)));
}

static void action_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_drop_down_new(G_LIST_MODEL(mouse_action_string_list()), NULL));
	gtk_list_item_set_child(li, GTK_WIDGET(dd));
	g_signal_connect(dd, "notify::selected", G_CALLBACK(action_dd_notify), NULL);
}

static void action_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", li);
	g_signal_handlers_block_by_func(dd, action_dd_notify, NULL);
	gtk_drop_down_set_selected(dd, find_string_position(mouse_action_choices, mouse_action_choice_count, row->action));
	g_signal_handlers_unblock_by_func(dd, action_dd_notify, NULL);
}

static void action_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", NULL);
}

/* --- modifier: combo (None / Shift / Ctrl / Ctrl+Shift) ------------ */

static const char *const modifier_choices[] = {
	NO_MODIFIER_TEXT, SHIFT_TEXT, CTRL_TEXT, CTRL_SHIFT_TEXT
};
static const gsize modifier_choice_count = G_N_ELEMENTS(modifier_choices);

static GtkStringList *modifier_string_list(void) {
	GtkStringList *sl = gtk_string_list_new(NULL);
	for (gsize i = 0; i < modifier_choice_count; i++)
		gtk_string_list_append(sl, modifier_choices[i]);
	return sl;
}

static void modifier_dd_notify_mouse(GObject *dd, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkListItem *li = g_object_get_data(dd, "siril-listitem");
	if (!li) return;
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	if (!row) return;
	GObject *sel = gtk_drop_down_get_selected_item(GTK_DROP_DOWN(dd));
	if (!sel || !GTK_IS_STRING_OBJECT(sel)) return;
	g_free(row->modifier);
	row->modifier = g_strdup(gtk_string_object_get_string(GTK_STRING_OBJECT(sel)));
}

static void modifier_mouse_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_drop_down_new(G_LIST_MODEL(modifier_string_list()), NULL));
	gtk_list_item_set_child(li, GTK_WIDGET(dd));
	g_signal_connect(dd, "notify::selected", G_CALLBACK(modifier_dd_notify_mouse), NULL);
}

static void modifier_mouse_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	SirilMouseActionItem *row = SIRIL_MOUSE_ACTION_ITEM(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", li);
	g_signal_handlers_block_by_func(dd, modifier_dd_notify_mouse, NULL);
	gtk_drop_down_set_selected(dd, find_string_position(modifier_choices, modifier_choice_count, row->modifier));
	g_signal_handlers_unblock_by_func(dd, modifier_dd_notify_mouse, NULL);
}

static void modifier_mouse_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", NULL);
}

/* --- scroll: NAME -------------------------------------------------- */

static GtkStringList *scroll_name_string_list(void) {
	GtkStringList *sl = gtk_string_list_new(NULL);
	for (guint i = 0; scroll_metadata_array[i]; i++)
		gtk_string_list_append(sl, scroll_metadata_array[i]->name);
	return sl;
}

static void scroll_name_dd_notify(GObject *dd, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkListItem *li = g_object_get_data(dd, "siril-listitem");
	if (!li) return;
	SirilScrollActionItem *row = SIRIL_SCROLL_ACTION_ITEM(gtk_list_item_get_item(li));
	if (!row) return;
	GObject *sel = gtk_drop_down_get_selected_item(GTK_DROP_DOWN(dd));
	if (!sel || !GTK_IS_STRING_OBJECT(sel)) return;
	const gchar *name = gtk_string_object_get_string(GTK_STRING_OBJECT(sel));
	const scroll_function_metadata *md = map_scroll_name_to_metadata(name);
	g_free(row->name);    row->name    = g_strdup(name);
	g_free(row->tooltip); row->tooltip = g_strdup(md->tooltip ? md->tooltip : "");
	row->reference = md->reference;
}

static void scroll_name_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_drop_down_new(G_LIST_MODEL(scroll_name_string_list()), NULL));
	gtk_list_item_set_child(li, GTK_WIDGET(dd));
	g_signal_connect(dd, "notify::selected", G_CALLBACK(scroll_name_dd_notify), NULL);
}

static void scroll_name_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	SirilScrollActionItem *row = SIRIL_SCROLL_ACTION_ITEM(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", li);
	guint pos = 0;
	for (guint i = 0; scroll_metadata_array[i]; i++) {
		if (!g_strcmp0(row->name, scroll_metadata_array[i]->name)) { pos = i; break; }
	}
	g_signal_handlers_block_by_func(dd, scroll_name_dd_notify, NULL);
	gtk_drop_down_set_selected(dd, pos);
	g_signal_handlers_unblock_by_func(dd, scroll_name_dd_notify, NULL);
	gtk_widget_set_tooltip_text(GTK_WIDGET(dd), row->tooltip);
}

static void scroll_name_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", NULL);
}

/* --- scroll: DIRECTION (Horizontal / Vertical) --------------------- */

static const char *const scroll_direction_choices[] = { SCROLL_HORIZ_TEXT, SCROLL_VERTICAL_TEXT };
static const gsize scroll_direction_count = G_N_ELEMENTS(scroll_direction_choices);

static GtkStringList *scroll_direction_string_list(void) {
	GtkStringList *sl = gtk_string_list_new(NULL);
	for (gsize i = 0; i < scroll_direction_count; i++)
		gtk_string_list_append(sl, scroll_direction_choices[i]);
	return sl;
}

static void scroll_direction_dd_notify(GObject *dd, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkListItem *li = g_object_get_data(dd, "siril-listitem");
	if (!li) return;
	SirilScrollActionItem *row = SIRIL_SCROLL_ACTION_ITEM(gtk_list_item_get_item(li));
	if (!row) return;
	GObject *sel = gtk_drop_down_get_selected_item(GTK_DROP_DOWN(dd));
	if (!sel || !GTK_IS_STRING_OBJECT(sel)) return;
	g_free(row->direction);
	row->direction = g_strdup(gtk_string_object_get_string(GTK_STRING_OBJECT(sel)));
}

static void scroll_direction_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_drop_down_new(G_LIST_MODEL(scroll_direction_string_list()), NULL));
	gtk_list_item_set_child(li, GTK_WIDGET(dd));
	g_signal_connect(dd, "notify::selected", G_CALLBACK(scroll_direction_dd_notify), NULL);
}

static void scroll_direction_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	SirilScrollActionItem *row = SIRIL_SCROLL_ACTION_ITEM(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", li);
	g_signal_handlers_block_by_func(dd, scroll_direction_dd_notify, NULL);
	gtk_drop_down_set_selected(dd, find_string_position(scroll_direction_choices, scroll_direction_count, row->direction));
	g_signal_handlers_unblock_by_func(dd, scroll_direction_dd_notify, NULL);
}

static void scroll_direction_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", NULL);
}

/* --- scroll: MODIFIER (shares choices with mouse) ------------------ */

static void modifier_dd_notify_scroll(GObject *dd, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkListItem *li = g_object_get_data(dd, "siril-listitem");
	if (!li) return;
	SirilScrollActionItem *row = SIRIL_SCROLL_ACTION_ITEM(gtk_list_item_get_item(li));
	if (!row) return;
	GObject *sel = gtk_drop_down_get_selected_item(GTK_DROP_DOWN(dd));
	if (!sel || !GTK_IS_STRING_OBJECT(sel)) return;
	g_free(row->modifier);
	row->modifier = g_strdup(gtk_string_object_get_string(GTK_STRING_OBJECT(sel)));
}

static void modifier_scroll_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_drop_down_new(G_LIST_MODEL(modifier_string_list()), NULL));
	gtk_list_item_set_child(li, GTK_WIDGET(dd));
	g_signal_connect(dd, "notify::selected", G_CALLBACK(modifier_dd_notify_scroll), NULL);
}

static void modifier_scroll_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	SirilScrollActionItem *row = SIRIL_SCROLL_ACTION_ITEM(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", li);
	g_signal_handlers_block_by_func(dd, modifier_dd_notify_scroll, NULL);
	gtk_drop_down_set_selected(dd, find_string_position(modifier_choices, modifier_choice_count, row->modifier));
	g_signal_handlers_unblock_by_func(dd, modifier_dd_notify_scroll, NULL);
}

static void modifier_scroll_unbind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkDropDown *dd = GTK_DROP_DOWN(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(dd), "siril-listitem", NULL);
}

/* ====================================================================
 * Build the GtkColumnView for mouse and scroll actions
 * ==================================================================== */

static GtkListItemFactory *make_factory(GCallback setup, GCallback bind, GCallback unbind) {
	GtkListItemFactory *f = gtk_signal_list_item_factory_new();
	g_signal_connect(f, "setup",   setup,   NULL);
	g_signal_connect(f, "bind",    bind,    NULL);
	if (unbind) g_signal_connect(f, "unbind", unbind, NULL);
	return f;
}

static gboolean fill_mouse_actions_list_idle(gpointer data) {
	(void) data;
	mouse_action_prefs_init_statics();
	GtkScrolledWindow *scrolled_window = mouse_scrolled_window;

	/* Drop existing column view if present so we can rebuild from scratch */
	gtk_scrolled_window_set_child(scrolled_window, NULL);
	if (mouse_action_store)    g_clear_object(&mouse_action_store);
	mouse_action_columnview = NULL;

	/* Build the GListStore from gui.mouse_actions */
	GListStore *store = g_list_store_new(SIRIL_TYPE_MOUSE_ACTION_ITEM);
	for (GSList *iterator = gui.mouse_actions; iterator; iterator = iterator->next) {
		mouse_action *action = (mouse_action *)iterator->data;
		const gchar *state_string = action->state & get_primary() && action->state & GDK_SHIFT_MASK ? CTRL_SHIFT_TEXT
				: action->state & get_primary() ? CTRL_TEXT
				: action->state & GDK_SHIFT_MASK ? SHIFT_TEXT : NO_MODIFIER_TEXT;
		const gchar *action_string = (action->type == GDK_BUTTON_PRESS) ? CLICK_TEXT : DOUBLE_CLICK_TEXT;
		SirilMouseActionItem *row = siril_mouse_action_item_new(
				action->data->name, action->button, action_string, state_string,
				action->data->tooltip, action->data->reference);
		g_list_store_append(store, row);
		g_object_unref(row);
	}
	mouse_action_store = store;

	/* Wrap in GtkSingleSelection and build the column view */
	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(store)));
	GtkColumnView *cv = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(sel)));
	gtk_widget_set_name(GTK_WIDGET(cv), "mouse_actions_columnview");

	GtkListItemFactory *fname = make_factory(G_CALLBACK(name_setup_cb),    G_CALLBACK(name_bind_cb),    G_CALLBACK(name_unbind_cb));
	GtkListItemFactory *fbtn  = make_factory(G_CALLBACK(button_setup_cb),  G_CALLBACK(button_bind_cb),  G_CALLBACK(button_unbind_cb));
	GtkListItemFactory *fact  = make_factory(G_CALLBACK(action_setup_cb),  G_CALLBACK(action_bind_cb),  G_CALLBACK(action_unbind_cb));
	GtkListItemFactory *fmod  = make_factory(G_CALLBACK(modifier_mouse_setup_cb), G_CALLBACK(modifier_mouse_bind_cb), G_CALLBACK(modifier_mouse_unbind_cb));

	GtkColumnViewColumn *c;
	c = gtk_column_view_column_new(N_("Name"),     fname); gtk_column_view_column_set_expand(c, TRUE); gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Button"),   fbtn);  gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Action"),   fact);  gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Modifier"), fmod);  gtk_column_view_append_column(cv, c); g_object_unref(c);

	gtk_scrolled_window_set_child(scrolled_window, GTK_WIDGET(cv));
	mouse_action_columnview = cv;
	/* No g_object_unref(sel): gtk_column_view_new() is transfer-full for
	 * the model, so cv consumed the ref we passed.  Unreffing here would
	 * destroy sel out from under cv and crash on the next rebuild. */

	return G_SOURCE_REMOVE;
}

void fill_mouse_actions_list(gboolean as_idle) {
	if (as_idle)
		g_idle_add(fill_mouse_actions_list_idle, NULL);
	else fill_mouse_actions_list_idle(NULL);
}

static gboolean validate_mouse_actions(GSList *list) {
	gboolean duplicate_found = FALSE, main_action_included = FALSE, doubleclick_conflict = FALSE;
	for (GSList *outer_iter = list; outer_iter != NULL; outer_iter = g_slist_next(outer_iter)) {
		mouse_action *outer_action = (mouse_action *)outer_iter->data;
		if (outer_action->data == &main_action) main_action_included = TRUE;
		for (GSList *inner_iter = g_slist_next(outer_iter); inner_iter != NULL; inner_iter = g_slist_next(inner_iter)) {
			mouse_action *inner_action = (mouse_action *)inner_iter->data;
			if (outer_action->type == inner_action->type &&
			    outer_action->state == inner_action->state &&
			    outer_action->button == inner_action->button) {
				duplicate_found = TRUE;
				siril_log_error(_("Duplicate mouse_actions found with type: %d, state: %d, button: %d\n"), outer_action->type, outer_action->state, outer_action->button);
			}
			if (outer_action->type == GDK_BUTTON_PRESS && !outer_action->data->doubleclick_compatible
			    && inner_action->type == GDK_DOUBLE_BUTTON_PRESS
			    && outer_action->button == inner_action->button
			    && outer_action->state == inner_action->state) {
				doubleclick_conflict = TRUE;
				siril_log_error(_("Double click action %s on button %u conflicts with single click action %s on the same button.\n"), inner_action->data->name, outer_action->button, outer_action->data->name);
			}
		}
	}
	if (!main_action_included) {
		siril_log_error(_("Main mouse action is not configured. This is essential for using Siril and must be assigned to a mouse action.\n"));
	}
	gboolean retval = main_action_included && !duplicate_found && !doubleclick_conflict;
	if (!retval) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid mouse action configuration"),
				_("There are conflicts between the actions you have configured. Please check the log for details and revise your action configuration."));
	}
	return retval;
}

/* Build a GSList<mouse_action*> from the live GListStore<SirilMouseActionItem*>. */
static void update_mouse_actions_from_store(GListStore *store) {
	GSList *new_mouse_actions = NULL;
	guint n = store ? g_list_model_get_n_items(G_LIST_MODEL(store)) : 0;
	for (guint i = 0; i < n; i++) {
		SirilMouseActionItem *row = g_list_model_get_item(G_LIST_MODEL(store), i);
		if (!g_strcmp0(row->name, NEW_ENTRY_TEXT)) {
			siril_log_warning(_("Warning: ignoring unconfigured action (\"" NEW_ENTRY_TEXT "\")\n"));
			g_object_unref(row);
			continue;
		}
		GdkEventType type = GDK_BUTTON_PRESS;
		if (!g_strcmp0(row->action, DOUBLE_CLICK_TEXT))
			type = GDK_DOUBLE_BUTTON_PRESS;
		GdkModifierType state = 0;
		if (g_strrstr(row->modifier, CTRL_TEXT))  state |= get_primary();
		if (g_strrstr(row->modifier, SHIFT_TEXT)) state |= GDK_SHIFT_MASK;
		const mouse_function_metadata *metadata = map_ref_to_metadata(row->reference);
		mouse_action *action = create_mouse_action(row->button, type, state, metadata);
		siril_log_debug("New action created: %s, ref %u, modifier %u, action %u\n",
				action->data->name, action->data->reference, action->state, action->type);
		new_mouse_actions = g_slist_append(new_mouse_actions, action);
		g_object_unref(row);
	}
	if (validate_mouse_actions(new_mouse_actions)) {
		if (gui.mouse_actions) g_slist_free_full(gui.mouse_actions, free);
		gui.mouse_actions = new_mouse_actions;
		fill_mouse_actions_list(FALSE);
		siril_log_info(_("Mouse configuration updated successfully.\n"));
	}
}

void on_mouse_actions_add_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	mouse_action_prefs_init_statics();
	if (!mouse_action_store) return;
	SirilMouseActionItem *row = siril_mouse_action_item_new(
			NEW_ENTRY_TEXT, 1, CLICK_TEXT, NO_MODIFIER_TEXT, "", MOUSE_REF_NULL);
	g_list_store_append(mouse_action_store, row);
	g_object_unref(row);
}

static void delete_selected_row_from_columnview(GtkColumnView *cv, GListStore *store) {
	if (!cv || !store) return;
	GtkSelectionModel *sm = gtk_column_view_get_model(cv);
	if (!GTK_IS_SINGLE_SELECTION(sm)) return;
	guint pos = gtk_single_selection_get_selected(GTK_SINGLE_SELECTION(sm));
	if (pos == GTK_INVALID_LIST_POSITION) return;
	g_list_store_remove(store, pos);
}

void on_mouse_actions_remove_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	mouse_action_prefs_init_statics();
	delete_selected_row_from_columnview(mouse_action_columnview, mouse_action_store);
}

/* ----- scroll list ------------------------------------------------- */

static gboolean fill_scroll_actions_list_idle(gpointer data) {
	(void) data;
	mouse_action_prefs_init_statics();
	GtkScrolledWindow *scrolled_window = scroll_scrolled_window;

	gtk_scrolled_window_set_child(scrolled_window, NULL);
	if (scroll_action_store) g_clear_object(&scroll_action_store);
	scroll_action_columnview = NULL;

	GListStore *store = g_list_store_new(SIRIL_TYPE_SCROLL_ACTION_ITEM);
	for (GSList *iterator = gui.scroll_actions; iterator; iterator = iterator->next) {
		scroll_action *action = (scroll_action *)iterator->data;
		const gchar *state_string = action->state & get_primary() && action->state & GDK_SHIFT_MASK ? CTRL_SHIFT_TEXT
				: action->state & get_primary() ? CTRL_TEXT
				: action->state & GDK_SHIFT_MASK ? SHIFT_TEXT : NO_MODIFIER_TEXT;
		const gchar *direction_string = action->direction == MOUSE_HORIZ_SCROLL ? SCROLL_HORIZ_TEXT : SCROLL_VERTICAL_TEXT;
		SirilScrollActionItem *row = siril_scroll_action_item_new(
				action->data->name, direction_string, state_string,
				action->data->tooltip, action->data->reference);
		g_list_store_append(store, row);
		g_object_unref(row);
	}
	scroll_action_store = store;

	GtkSingleSelection *sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(store)));
	GtkColumnView *cv = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(sel)));
	gtk_widget_set_name(GTK_WIDGET(cv), "scroll_actions_columnview");

	GtkListItemFactory *fname = make_factory(G_CALLBACK(scroll_name_setup_cb),      G_CALLBACK(scroll_name_bind_cb),      G_CALLBACK(scroll_name_unbind_cb));
	GtkListItemFactory *fdir  = make_factory(G_CALLBACK(scroll_direction_setup_cb), G_CALLBACK(scroll_direction_bind_cb), G_CALLBACK(scroll_direction_unbind_cb));
	GtkListItemFactory *fmod  = make_factory(G_CALLBACK(modifier_scroll_setup_cb),  G_CALLBACK(modifier_scroll_bind_cb),  G_CALLBACK(modifier_scroll_unbind_cb));

	GtkColumnViewColumn *c;
	c = gtk_column_view_column_new(N_("Name"),      fname); gtk_column_view_column_set_expand(c, TRUE); gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Direction"), fdir);  gtk_column_view_append_column(cv, c); g_object_unref(c);
	c = gtk_column_view_column_new(N_("Modifier"),  fmod);  gtk_column_view_append_column(cv, c); g_object_unref(c);

	gtk_scrolled_window_set_child(scrolled_window, GTK_WIDGET(cv));
	scroll_action_columnview = cv;
	/* No g_object_unref(sel): transferred to cv (see mouse list above). */

	return G_SOURCE_REMOVE;
}

void fill_scroll_actions_list(gboolean as_idle) {
	if (as_idle)
		g_idle_add(fill_scroll_actions_list_idle, NULL);
	else fill_scroll_actions_list_idle(NULL);
}

static gboolean validate_scroll_actions(GSList *list) {
	gboolean duplicate_found = FALSE;
	for (GSList *outer_iter = list; outer_iter != NULL; outer_iter = g_slist_next(outer_iter)) {
		scroll_action *outer_action = (scroll_action *)outer_iter->data;
		for (GSList *inner_iter = g_slist_next(outer_iter); inner_iter != NULL; inner_iter = g_slist_next(inner_iter)) {
			scroll_action *inner_action = (scroll_action *)inner_iter->data;
			if (outer_action->state == inner_action->state &&
			    outer_action->direction == inner_action->direction) {
				duplicate_found = TRUE;
				siril_log_error(_("Duplicate scroll_actions found with state: %d, direction: %d\n"), outer_action->state, outer_action->direction);
			}
		}
	}
	gboolean retval = !duplicate_found;
	if (!retval) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid mouse scroll action configuration"),
				_("There are conflicts between the actions you have configured. Please check the log for details and revise your action configuration."));
	}
	return retval;
}

static void update_scroll_actions_from_store(GListStore *store) {
	GSList *new_scroll_actions = NULL;
	guint n = store ? g_list_model_get_n_items(G_LIST_MODEL(store)) : 0;
	for (guint i = 0; i < n; i++) {
		SirilScrollActionItem *row = g_list_model_get_item(G_LIST_MODEL(store), i);
		if (!g_strcmp0(row->name, NEW_ENTRY_TEXT)) {
			siril_log_warning(_("Warning: ignoring unconfigured action (\"" NEW_ENTRY_TEXT "\")\n"));
			g_object_unref(row);
			continue;
		}
		SirilScrollDirection dir = MOUSE_HORIZ_SCROLL;
		if (!g_strcmp0(row->direction, SCROLL_VERTICAL_TEXT))
			dir = MOUSE_VERTICAL_SCROLL;
		GdkModifierType state = 0;
		if (g_strrstr(row->modifier, CTRL_TEXT))  state |= get_primary();
		if (g_strrstr(row->modifier, SHIFT_TEXT)) state |= GDK_SHIFT_MASK;
		const scroll_function_metadata *metadata = map_scroll_ref_to_metadata(row->reference);
		scroll_action *action = create_scroll_action(state, metadata, dir);
		siril_log_debug("New action created: %s, reference %u, modifier %u, action %u\n",
				action->data->name, action->data->reference, action->state, action->direction);
		new_scroll_actions = g_slist_append(new_scroll_actions, action);
		g_object_unref(row);
	}
	if (validate_scroll_actions(new_scroll_actions)) {
		if (gui.scroll_actions) g_slist_free_full(gui.scroll_actions, free);
		gui.scroll_actions = new_scroll_actions;
		fill_scroll_actions_list(FALSE);
		siril_log_info(_("Mouse scroll configuration updated successfully.\n"));
	}
}

void on_scroll_actions_add_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	mouse_action_prefs_init_statics();
	if (!scroll_action_store) return;
	SirilScrollActionItem *row = siril_scroll_action_item_new(
			NEW_ENTRY_TEXT, SCROLL_VERTICAL_TEXT, NO_MODIFIER_TEXT, "", SCROLL_REF_NULL);
	g_list_store_append(scroll_action_store, row);
	g_object_unref(row);
}

void on_scroll_actions_remove_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	mouse_action_prefs_init_statics();
	delete_selected_row_from_columnview(scroll_action_columnview, scroll_action_store);
}

void on_mouse_actions_reset_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	if (siril_confirm_dialog(_("Reset Mouse Actions"),
			_("This will reset all mouse button and scroll actions to the default settings. "
			  "It is not possible to undo this and you will need to reconfigure any changes. Are you sure?"),
			_("Reset"))) {
		if (gui.mouse_actions) {
			g_slist_free_full(gui.mouse_actions, free);
			gui.mouse_actions = NULL;
		}
		initialize_mouse_actions();
		fill_mouse_actions_list(FALSE);

		if (gui.scroll_actions) {
			g_slist_free_full(gui.scroll_actions, free);
			gui.scroll_actions = NULL;
		}
		initialize_scroll_actions();
		fill_scroll_actions_list(FALSE);
	}
}

void on_mouse_actions_apply_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	control_window_switch_to_tab(OUTPUT_LOGS);
	mouse_action_prefs_init_statics();
	com.pref.gui.mouse_speed_limit = gtk_spin_button_get_value(mouse_speed_spin);
	if (mouse_action_store)
		update_mouse_actions_from_store(mouse_action_store);
	if (scroll_action_store)
		update_scroll_actions_from_store(scroll_action_store);
	if (com.pref.gui.mouse_cfg.mouse_actions_array)
		g_slist_free_full(com.pref.gui.mouse_cfg.mouse_actions_array, g_free);
	com.pref.gui.mouse_cfg.mouse_actions_array = mouse_actions_list_to_config(gui.mouse_actions);
	if (com.pref.gui.mouse_cfg.scroll_actions_array)
		g_slist_free_full(com.pref.gui.mouse_cfg.scroll_actions_array, g_free);
	com.pref.gui.mouse_cfg.scroll_actions_array = scroll_actions_list_to_config(gui.scroll_actions);
	writeinitfile();
}

void on_mouse_actions_dialog_show(GtkWidget *dialog, gpointer user_data) {
	(void)dialog; (void)user_data;
	mouse_action_prefs_init_statics();

	siril_register_css_for_display("siril-mouse-test-frame",
		"frame.mouse-test-frame { border: 2px solid #aaa; border-radius: 4px; padding: 2px; }");
	gtk_widget_add_css_class(mouse_test_frame_widget, "mouse-test-frame");
	install_mouse_test_controllers(mouse_test_drawingarea);
	if (!gui.mouse_actions)  load_or_initialize_mouse_actions();
	if (!gui.scroll_actions) load_or_initialize_scroll_actions();

	fill_mouse_actions_list(FALSE);
	fill_scroll_actions_list(FALSE);
}

void on_mouse_actions_dialog_hide(GtkWidget *dialog, gpointer user_data) {
	(void)dialog; (void)user_data;
	siril_close_dialog("mouse_actions_dialog");
}

void on_mouse_actions_close_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	siril_close_dialog("mouse_actions_dialog");
}

gint timeout_ref = 0;

static gboolean reset_label_text(GtkLabel *label) {
    gtk_label_set_text(label, N_("Mouse button / scroll check"));
	timeout_ref = 0;
    return G_SOURCE_REMOVE;
}

static void update_test_label(const gchar *text) {
	if (timeout_ref) {
		g_source_remove(timeout_ref);
		timeout_ref = 0;
	}
	if (text)
		gtk_label_set_text(mouse_test_label, text);
	timeout_ref = g_timeout_add(1250, (GSourceFunc)reset_label_text, mouse_test_label);
}

static void on_mouse_test_pressed(GtkGestureClick *gesture, int n_press,
                                  double x, double y, gpointer user_data) {
	(void)x; (void)y; (void)user_data; (void)n_press;
	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	gchar *text = g_strdup_printf(_("Button %u"), button);
	update_test_label(text);
	g_free(text);
}

/* Report drag motion live so the user can verify that motion-while-
 * pressed is wired correctly on their platform — the macOS gesture-
 * routing bug fixed in the histo/plot/canvas areas went unnoticed
 * because there was no way to test drag delivery here.  drag-update
 * fires for every motion event after the drag threshold is crossed;
 * drag-end fires once on release with the final offset. */
static void on_mouse_test_drag_update(GtkGestureDrag *gesture,
                                      double offset_x, double offset_y, gpointer user_data) {
	(void)user_data;
	/* GtkGestureDrag fires drag-update even before the drag threshold
	 * is crossed (its threshold only gates drag-begin).  Ignore sub-
	 * pixel movement so a tap-and-release reads as "Button N" rather
	 * than "Drag (button N, +0, +0)". */
	if (fabs(offset_x) < 1.0 && fabs(offset_y) < 1.0) return;
	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	gchar *text = g_strdup_printf(_("Drag (button %u, %+.0f, %+.0f)"),
	                              button, offset_x, offset_y);
	update_test_label(text);
	g_free(text);
}

static void on_mouse_test_drag_end(GtkGestureDrag *gesture,
                                   double offset_x, double offset_y, gpointer user_data) {
	(void)user_data;
	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	gchar *text = g_strdup_printf(_("Drag end (button %u, %+.0f, %+.0f)"),
	                              button, offset_x, offset_y);
	update_test_label(text);
	g_free(text);
}

static gboolean on_mouse_test_scroll(GtkEventControllerScroll *ctrl,
                                     double dx, double dy, gpointer user_data) {
	(void)ctrl; (void)user_data;
	GdkEvent *event = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(ctrl));
	GdkScrollDirection dir = GDK_SCROLL_SMOOTH;
	if (event && gdk_event_get_event_type(event) == GDK_SCROLL)
		dir = gdk_scroll_event_get_direction(event);
	const gchar *text = NULL;
	switch (dir) {
		case GDK_SCROLL_UP:
		case GDK_SCROLL_DOWN:
			text = _("Scroll (vertical)");
			break;
		case GDK_SCROLL_LEFT:
		case GDK_SCROLL_RIGHT:
			text = _("Scroll (horizontal)");
			break;
		case GDK_SCROLL_SMOOTH:
		default:
			if (dy != 0.0)
				text = _("Smooth scroll (vertical)");
			else if (dx != 0.0)
				text = _("Smooth scroll (horizontal)");
			break;
	}
	if (text)
		update_test_label(text);
	return TRUE;
}

static void install_mouse_test_controllers(GtkWidget *area) {
	if (g_object_get_data(G_OBJECT(area), "siril-mouse-test-installed"))
		return;
	GtkGesture *click = gtk_gesture_click_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click), 0); /* any button */
	g_signal_connect(click, "pressed", G_CALLBACK(on_mouse_test_pressed), NULL);
	gtk_widget_add_controller(area, GTK_EVENT_CONTROLLER(click));

	/* Drag too, so the user can verify motion-while-pressed delivery —
	 * grouped with click so neither claim denies the other. */
	GtkGesture *drag = gtk_gesture_drag_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag), 0); /* any button */
	g_signal_connect(drag, "drag-update",
	                 G_CALLBACK(on_mouse_test_drag_update), NULL);
	g_signal_connect(drag, "drag-end",
	                 G_CALLBACK(on_mouse_test_drag_end),    NULL);
	gtk_widget_add_controller(area, GTK_EVENT_CONTROLLER(drag));
	gtk_gesture_group(click, drag);

	GtkEventController *scroll = gtk_event_controller_scroll_new(
		GTK_EVENT_CONTROLLER_SCROLL_BOTH_AXES);
	g_signal_connect(scroll, "scroll", G_CALLBACK(on_mouse_test_scroll), NULL);
	gtk_widget_add_controller(area, scroll);

	g_object_set_data(G_OBJECT(area), "siril-mouse-test-installed", GINT_TO_POINTER(1));
}

GSList *mouse_actions_list_to_config(GSList *actions) {
	GSList *descriptions = NULL;
	GSList *current = actions;
	while (current != NULL) {
		mouse_action *action = (mouse_action *)current->data;
		if (!action->data || action->data->reference == MOUSE_REF_NULL || action->data->reference >= MOUSE_REF_MAX)
			continue;
		gchar *description = g_strdup_printf(
			"%u,%u,%u,%u",
			action->button,
			(guint) action->type,
			(guint) action->state,
			(guint) action->data->reference
		);
		descriptions = g_slist_prepend(descriptions, description);
		current = g_slist_next(current);
	}
	descriptions = g_slist_reverse(descriptions);
	return descriptions;
}

GSList *mouse_actions_config_to_list(GSList *config) {
	GSList *actions = NULL;
	GSList *current = config;
	while (current != NULL) {
		gchar *input = (gchar*)current->data;
		gchar **tokens = g_strsplit(input, ",", -1);
		mouse_action *action = (mouse_action*) malloc(sizeof(mouse_action));
		if (tokens[0] && tokens[1] && tokens[2] && tokens[3]) {
			action->button = g_ascii_strtoull(tokens[0], NULL, 10);
			action->type   = g_ascii_strtoull(tokens[1], NULL, 10);
			action->state  = g_ascii_strtoull(tokens[2], NULL, 10);
			mouse_function_ref reference = (mouse_function_ref) g_ascii_strtoull(tokens[3], NULL, 10);
			if (reference == MOUSE_REF_NULL || reference >= MOUSE_REF_MAX || map_ref_to_metadata(reference) == &null_action) {
				siril_log_warning(_("Warning: when parsing mouse action config, config string parsed to an unknown function: \"%s\". Skipping..."), input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
			action->data = map_ref_to_metadata(reference);
			if (action->data == &null_action) {
				siril_log_warning(_("Warning: when parsing mouse action config, config string parsed to an unknown function: \"%s\". Skipping..."), input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
		} else {
			siril_log_warning(_("Warning: when parsing mouse action config, config string has incorrect format: \"%s\". Skipping..."), input);
			free(action);
			g_strfreev(tokens);
			current = g_slist_next(current);
			continue;
		}
		actions = g_slist_append(actions, action);
		g_strfreev(tokens);
		current = g_slist_next(current);
	}
	return actions;
}

GSList *scroll_actions_list_to_config(GSList *actions) {
	GSList *descriptions = NULL;
	GSList *current = actions;
	while (current != NULL) {
		scroll_action *action = (scroll_action *)current->data;
		if (!action->data || action->data->reference == SCROLL_REF_NULL || action->data->reference >= SCROLL_REF_MAX)
			continue;
		gchar *description = g_strdup_printf(
			"%u,%u,%u",
			(guint)action->direction,
			(guint)action->state,
			(guint) action->data->reference
		);
		descriptions = g_slist_prepend(descriptions, description);
		current = g_slist_next(current);
	}
	descriptions = g_slist_reverse(descriptions);
	return descriptions;
}

GSList *scroll_actions_config_to_list(GSList *config) {
	GSList *actions = NULL;
	GSList *current = config;
	while (current != NULL) {
		gchar *input = (gchar*)current->data;
		gchar **tokens = g_strsplit(input, ",", -1);
		scroll_action *action = (scroll_action*) malloc(sizeof(scroll_action));
		if (tokens[0] && tokens[1] && tokens[2]) {
			action->direction = g_ascii_strtoull(tokens[0], NULL, 10);
			action->state     = g_ascii_strtoull(tokens[1], NULL, 10);
			scroll_function_ref reference = (scroll_function_ref) g_ascii_strtoull(tokens[2], NULL, 10);
			if (reference == SCROLL_REF_NULL || reference >= SCROLL_REF_MAX || map_scroll_ref_to_metadata(reference) == &scroll_null_action) {
				siril_log_warning(_("Warning: when parsing mouse action config, config string parsed to an unknown function: \"%s\". Skipping..."), input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
			action->data = map_scroll_ref_to_metadata(reference);
			if (action->data == &scroll_null_action) {
				siril_log_warning(_("Warning: when parsing scroll action config, string parsed to an unknown function: \"%s\". Skipping...\n"), input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
		} else {
			siril_log_warning(_("Warning: when parsing mouse action config, string has incorrect format: \"%s\". Skipping...\n"), input);
			free(action);
			g_strfreev(tokens);
			current = g_slist_next(current);
			continue;
		}
		actions = g_slist_append(actions, action);
		g_strfreev(tokens);
		current = g_slist_next(current);
	}
	return actions;
}
