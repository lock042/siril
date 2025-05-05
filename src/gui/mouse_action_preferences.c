/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include "core/initfile.h"
#include "core/siril_log.h"
#include "gui/image_interactions.h"
#include "gui/mouse_action_functions.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"

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

// GtkTreeView column identifiers
enum {
	COLUMN_ACTION = 0,	// guint
	COLUMN_BUTTON,		// guint
	COLUMN_MODIFIER,	// guint
	COLUMN_FUNCTION,	// string
	COLUMN_TOOLTIP,		// string
	N_COLUMNS
};

enum {
	NONE = 0,
	SHIFT = 1,
	CTRL = 2,
	CTRLSHIFT = 3
};

void on_config_mouse_buttons_clicked(GtkDialog *dialog, gpointer user_data) {
	if (siril_confirm_dialog(_("Mouse Configuration"), _("Warning: this dialog allows you to reconfigure mouse behaviour so that it may no longer match the documentation. If you forget what settings you have changed, the dialog can also be used to view current mouse action assignments and also provides a button to reset them to default settings."), _("Proceed"))) {
		siril_open_dialog("mouse_actions_dialog");
	}
}

// Callback function for when the name cell is edited
static void name_cell_edited(GtkCellRendererText *renderer, const gchar *path, const gchar *new_text, GtkListStore *store) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter, path);
	if (valid) {
		// Update the value in the list store
		gtk_list_store_set(store, &iter, 0, new_text, -1);
		// Find the corresponding tooltip and unique reference and set those too
		const mouse_function_metadata *data = map_name_to_metadata(new_text);
		gtk_list_store_set(store, &iter, 4, data->tooltip, 5, data->reference, -1);
    }
}

// Callback function for when the button cell is edited
static void button_cell_edited(GtkCellRendererText *renderer, const gchar *path, const gchar *new_text, GtkListStore *store) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter, path);
	if (valid) {
		// Update the value in the list store
		gtk_list_store_set(store, &iter, 1, atoi(new_text), -1);
	}
}

// Callback function for when the state cell is edited
static void action_cell_edited(GtkCellRendererText *renderer, const gchar *path, const gchar *new_text, GtkListStore *store) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter, path);
	if (valid) {
		// Update the value in the list store
		gtk_list_store_set(store, &iter, 2, new_text, -1);
	}
}

// Callback function for when the state cell is edited
static void state_cell_edited(GtkCellRendererText *renderer, const gchar *path, const gchar *new_text, GtkListStore *store) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter, path);
	if (valid) {
		// Update the value in the list store
		gtk_list_store_set(store, &iter, 3, new_text, -1);
	}
}

static gboolean fill_mouse_actions_list_idle() {
	// Check if there is already a GtkTreeView child of the scrollable window
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("mouse_treeview_scrolled_window"));
	GtkWidget *existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "mouse_actions_treeview") == 0) {
		gtk_widget_destroy(existing_tree_view);
	}

	// Create a list store
	GtkListStore *store = gtk_list_store_new(6, G_TYPE_STRING, G_TYPE_INT, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_INT);

	// Populate the list store with mouse actions
	GSList *iterator;
	for (iterator = gui.mouse_actions; iterator; iterator = iterator->next) {
		mouse_action *action = (mouse_action *)iterator->data;
		GtkTreeIter iter;
		gtk_list_store_append(store, &iter);
		const gchar *state_string = action->state & get_primary() && action->state & GDK_SHIFT_MASK ? CTRL_SHIFT_TEXT : action->state & get_primary() ? CTRL_TEXT : action->state & GDK_SHIFT_MASK ? SHIFT_TEXT : NO_MODIFIER_TEXT;
		gtk_list_store_set(store, &iter, 0, action->data->name, 1, action->button, 2, (action->type == GDK_BUTTON_PRESS) ? CLICK_TEXT : DOUBLE_CLICK_TEXT, 3, state_string, 4, action->data->tooltip, 5, action->data->reference, -1);
	}

	// Create the tree view
	GtkWidget *tree_view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
	g_object_unref(store); // the tree view holds the reference now
    gtk_widget_set_name(tree_view, "mouse_actions_treeview"); // for the check next time this function is called

	// Create and add columns
	GtkCellRenderer *renderer;
	GtkTreeViewColumn *column;

	renderer = gtk_cell_renderer_combo_new();
	GtkListStore *state_store = gtk_list_store_new(1, G_TYPE_STRING);
	int i = 0;
	while (mouse_metadata_array[i]) {
		gtk_list_store_insert_with_values(state_store, NULL, -1, 0, mouse_metadata_array[i++]->name, -1);
	}
	g_object_set(renderer, "editable", TRUE, "model", state_store, "text-column", 0, NULL);
	g_signal_connect(renderer, "edited", G_CALLBACK(name_cell_edited), store); // Connect the edited signal
	column = gtk_tree_view_column_new_with_attributes(N_("Name"), renderer, "text", 0, NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	renderer = gtk_cell_renderer_spin_new();
	g_object_set(renderer, "adjustment", gtk_adjustment_new(1, 1, 15, 1, 0, 0), NULL);
	g_object_set(renderer, "editable", TRUE, "digits", 0, NULL); // Make the spin button editable and set digits to 0
	g_signal_connect(renderer, "edited", G_CALLBACK(button_cell_edited), store); // Connect the edited signal
	column = gtk_tree_view_column_new_with_attributes(N_("Button"), renderer, "text", 1, NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	renderer = gtk_cell_renderer_combo_new();
	state_store = gtk_list_store_new(1, G_TYPE_STRING);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, CLICK_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, DOUBLE_CLICK_TEXT, -1);
	g_object_set(renderer, "editable", TRUE, "model", state_store, "text-column", 0, NULL);
	g_signal_connect(renderer, "edited", G_CALLBACK(action_cell_edited), store); // Connect the edited signal
	column = gtk_tree_view_column_new_with_attributes(N_("Action"), renderer, "text", 2, NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	// Create a combo box renderer for the state column
	renderer = gtk_cell_renderer_combo_new();
	state_store = gtk_list_store_new(1, G_TYPE_STRING);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, CTRL_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, SHIFT_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, CTRL_SHIFT_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, NO_MODIFIER_TEXT, -1);
	g_object_set(renderer, "editable", TRUE, "model", state_store, "text-column", 0, NULL);
	g_signal_connect(renderer, "edited", G_CALLBACK(state_cell_edited), store); // Connect the edited signal
	column = gtk_tree_view_column_new_with_attributes(N_("Modifier"), renderer, "text", 3, NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	// Set tooltip column
	gtk_tree_view_set_tooltip_column(GTK_TREE_VIEW(tree_view), 4);

	// Add the tree view to the scrolled window
	gtk_container_add(GTK_CONTAINER(scrolled_window), tree_view);
	gtk_widget_show_all(tree_view);

	return FALSE;
}

/* called on mouse preference window loading.
 * It is executed safely in the GTK thread if as_idle is true. */
void fill_mouse_actions_list(gboolean as_idle) {
	if (as_idle)
		gdk_threads_add_idle(fill_mouse_actions_list_idle, NULL);
	else fill_mouse_actions_list_idle();
}

static gboolean validate_mouse_actions(GSList *list) {

	gboolean duplicate_found = FALSE, main_action_included = FALSE, doubleclick_conflict = FALSE;
	for (GSList *outer_iter = list; outer_iter != NULL; outer_iter = g_slist_next(outer_iter)) {
		mouse_action *outer_action = (mouse_action *)outer_iter->data;

		// Check the main mouse action is included
		if (outer_action->data == &main_action) {
			main_action_included = TRUE;
		}

		// Inner loop to compare the current mouse_action with subsequent mouse_actions
		for (GSList *inner_iter = g_slist_next(outer_iter); inner_iter != NULL; inner_iter = g_slist_next(inner_iter)) {
			mouse_action *inner_action = (mouse_action *)inner_iter->data;

			// Check there are no actions duplicating button, single / double click and modifiers
			if (outer_action->type == inner_action->type &&
						outer_action->state == inner_action->state &&
						outer_action->button == inner_action->button) {
				duplicate_found = TRUE;
				siril_log_color_message(_("Duplicate mouse_actions found with type: %d, state: %d, button: %d\n"), "red", outer_action->type, outer_action->state, outer_action->button);
			}
			// Check for actions assigned to double clicks where there is already an action
			// with lasting impact assigned to the same single click
			if (outer_action->type == GDK_BUTTON_PRESS && !outer_action->data->doubleclick_compatible && inner_action->type == GDK_DOUBLE_BUTTON_PRESS && outer_action->button == inner_action->button && outer_action->state == inner_action->state) {
				doubleclick_conflict = TRUE;
				siril_log_color_message(_("Double click action %s on button %u conflicts with single click action %s on the same button.\n"), "red", inner_action->data->name, outer_action->button, outer_action->data->name);
			}
		}
	}
	if (!main_action_included) {
		siril_log_color_message(_("Main mouse action is not configured. This is essential for using Siril and must be assigned to a mouse action.\n"), "red");
	}

	gboolean retval = main_action_included && !duplicate_found && !doubleclick_conflict;
	if (!retval) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid mouse action configuration"), _("There are conflicts between the actions you have configured. Please check the log for details and revise your action configuration."));
	}
	return retval;
}

void update_mouse_actions_from_treeview(GtkTreeView *tree_view) {
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_first(model, &iter);

	// Clear existing mouse actions list
	GSList *new_mouse_actions = NULL;

	// Iterate over rows in the tree view
	while (valid) {
		gchar *name;
		gint button;
		gchar *action_type;
		gchar *modifier;
		gint reference;

		// Get data from each column
		gtk_tree_model_get(model, &iter, 0, &name, 1, &button, 2, &action_type, 3, &modifier, 5, &reference, -1);

		// Convert action type from string to GdkEventType
		GdkEventType type = GDK_BUTTON_PRESS;
		if (strcmp(action_type, DOUBLE_CLICK_TEXT) == 0) {
			type = GDK_DOUBLE_BUTTON_PRESS;
		}

		// Convert modifier from string to GdkModifierType
		GdkModifierType state = 0;
		if (g_strrstr(modifier, CTRL_TEXT)) {
			state |= get_primary();
		}
		if (g_strrstr(modifier, SHIFT_TEXT)) {
			state |= GDK_SHIFT_MASK;
		}
		// Check for un-set new action
		if (!g_strcmp0(name, NEW_ENTRY_TEXT)) {
			g_free(name);
			g_free(action_type);
			g_free(modifier);
			siril_log_color_message(_("Warning: ignoring unconfigured action (\"" NEW_ENTRY_TEXT "\")\n"), "salmon");
			valid = gtk_tree_model_iter_next(model, &iter);
			continue;
		}

		// Convert name to metadata
		const mouse_function_metadata *metadata = map_ref_to_metadata(reference);

		// Create a new mouse action and add it to the list
		mouse_action *action = create_mouse_action(button, type, state, metadata);
		siril_debug_print("New action created: %s, ref %u, modifier %u, action %u\n", action->data->name, action->data->reference, action->state, action->type);

		new_mouse_actions = g_slist_append(new_mouse_actions, action);

		g_free(name);
		g_free(action_type);
		g_free(modifier);

		valid = gtk_tree_model_iter_next(model, &iter);
	}
	if (validate_mouse_actions(new_mouse_actions)) {
		if (gui.mouse_actions) {
			g_slist_free_full(gui.mouse_actions, free);
		}
		gui.mouse_actions = new_mouse_actions;
		fill_mouse_actions_list(FALSE);
		siril_log_color_message(_("Mouse configuration updated successfully.\n"), "green");
	}
}

void add_row_to_tree_view(GtkListStore *store) {
	GtkTreeIter iter;
	gtk_list_store_append(store, &iter);
	// Add default values for each column
	gtk_list_store_set(store, &iter,
					0, NEW_ENTRY_TEXT,
					1, 1,
					2, CLICK_TEXT,
					3, NO_MODIFIER_TEXT,
					4, "",
					5, MOUSE_REF_NULL,
					-1);
}

void on_mouse_actions_add_clicked(GtkButton *button, gpointer user_data) {
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("mouse_treeview_scrolled_window"));
	GtkWidget *existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "mouse_actions_treeview") == 0) {
		GtkTreeView *tree_view = GTK_TREE_VIEW(existing_tree_view);
		GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
		GtkListStore *store = GTK_LIST_STORE(model);
		add_row_to_tree_view(store);
	}
}

static void delete_selected_row(GtkTreeView *tree_view, GtkListStore *store) {
	GtkTreeSelection *selection = gtk_tree_view_get_selection(tree_view);
	GtkTreeModel *model;
	GtkTreeIter iter;
	if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
		gtk_list_store_remove(store, &iter);
	}
}

void on_mouse_actions_remove_clicked(GtkButton *button, gpointer user_data) {
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("mouse_treeview_scrolled_window"));
	GtkWidget *existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "mouse_actions_treeview") == 0) {
		GtkTreeView *tree_view = GTK_TREE_VIEW(existing_tree_view);
		GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
		GtkListStore *store = GTK_LIST_STORE(model);
		delete_selected_row(tree_view, store);
	}
}

// Callback function for when the name cell is edited
static void scroll_name_cell_edited(GtkCellRendererText *renderer, const gchar *path, const gchar *new_text, GtkListStore *store) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter, path);
	if (valid) {
		// Update the value in the list store
		gtk_list_store_set(store, &iter, 0, new_text, -1);
		// Find the corresponding tooltip and reference and set those too
		const scroll_function_metadata *data = map_scroll_name_to_metadata(new_text);
		gtk_list_store_set(store, &iter, 3, data->tooltip, 4, data->reference, -1);
    }
}

// Callback function for when the state cell is edited
static void scroll_direction_cell_edited(GtkCellRendererText *renderer, const gchar *path, const gchar *new_text, GtkListStore *store) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter, path);
	if (valid) {
		// Update the value in the list store
		gtk_list_store_set(store, &iter, 1, new_text, -1);
	}
}

// Callback function for when the state cell is edited
static void scroll_state_cell_edited(GtkCellRendererText *renderer, const gchar *path, const gchar *new_text, GtkListStore *store) {
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(store), &iter, path);
	if (valid) {
		// Update the value in the list store
		gtk_list_store_set(store, &iter, 2, new_text, -1);
	}
}

static gboolean fill_scroll_actions_list_idle() {
	// Check if there is already a GtkTreeView child of the scrollable window
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("scroll_treeview_scrolled_window"));
	GtkWidget *existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "scroll_actions_treeview") == 0) {
		gtk_widget_destroy(existing_tree_view);
	}

	// Create a list store
	GtkListStore *store = gtk_list_store_new(5, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_INT);

	// Populate the list store with mouse actions
	GSList *iterator;
	for (iterator = gui.scroll_actions; iterator; iterator = iterator->next) {
		scroll_action *action = (scroll_action *)iterator->data;
		GtkTreeIter iter;
		gtk_list_store_append(store, &iter);
		const gchar *state_string = action->state & get_primary() && action->state & GDK_SHIFT_MASK ? CTRL_SHIFT_TEXT : action->state & get_primary() ? CTRL_TEXT : action->state & GDK_SHIFT_MASK ? SHIFT_TEXT : NO_MODIFIER_TEXT;
		const gchar* direction_string = action->direction == MOUSE_HORIZ_SCROLL ? SCROLL_HORIZ_TEXT : SCROLL_VERTICAL_TEXT;
		gtk_list_store_set(store, &iter, 0, action->data->name, 1, direction_string, 2, state_string, 3, action->data->tooltip, 4, action->data->reference, -1);
	}

	// Create the tree view
	GtkWidget *tree_view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
	g_object_unref(store); // the tree view holds the reference now
    gtk_widget_set_name(tree_view, "scroll_actions_treeview"); // for the check next time this function is called

	// Create and add columns
	GtkCellRenderer *renderer;
	GtkTreeViewColumn *column;

	renderer = gtk_cell_renderer_combo_new();
    GtkListStore *state_store = gtk_list_store_new(1, G_TYPE_STRING);
	int i = 0;
	while (scroll_metadata_array[i]) {
		gtk_list_store_insert_with_values(state_store, NULL, -1, 0, scroll_metadata_array[i++]->name, -1);
	}
    g_object_set(renderer, "editable", TRUE, "model", state_store, "text-column", 0, NULL);
    g_signal_connect(renderer, "edited", G_CALLBACK(scroll_name_cell_edited), store); // Connect the edited signal
	column = gtk_tree_view_column_new_with_attributes(N_("Name"), renderer, "text", 0, NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	renderer = gtk_cell_renderer_combo_new();
	state_store = gtk_list_store_new(1, G_TYPE_STRING);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, SCROLL_HORIZ_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, SCROLL_VERTICAL_TEXT, -1);
	g_object_set(renderer, "editable", TRUE, "model", state_store, "text-column", 0, NULL);
	g_signal_connect(renderer, "edited", G_CALLBACK(scroll_direction_cell_edited), store); // Connect the edited signal
	column = gtk_tree_view_column_new_with_attributes(N_("Direction"), renderer, "text", 1, NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	// Create a combo box renderer for the state column
	renderer = gtk_cell_renderer_combo_new();
	state_store = gtk_list_store_new(1, G_TYPE_STRING);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, CTRL_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, SHIFT_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, CTRL_SHIFT_TEXT, -1);
	gtk_list_store_insert_with_values(state_store, NULL, -1, 0, NO_MODIFIER_TEXT, -1);
	g_object_set(renderer, "editable", TRUE, "model", state_store, "text-column", 0, NULL);
	g_signal_connect(renderer, "edited", G_CALLBACK(scroll_state_cell_edited), store); // Connect the edited signal
	column = gtk_tree_view_column_new_with_attributes(N_("Modifier"), renderer, "text", 2, NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(tree_view), column);

	// Set tooltip column
	gtk_tree_view_set_tooltip_column(GTK_TREE_VIEW(tree_view), 3);

	// Add the tree view to the scrolled window
	gtk_container_add(GTK_CONTAINER(scrolled_window), tree_view);
	gtk_widget_show_all(tree_view);

	return FALSE;
}

/* called on preference window loading.
 * It is executed safely in the GTK thread if as_idle is true. */
void fill_scroll_actions_list(gboolean as_idle) {

	if (as_idle)
		gdk_threads_add_idle(fill_scroll_actions_list_idle, NULL);
	else fill_scroll_actions_list_idle();
}

static gboolean validate_scroll_actions(GSList *list) {

	gboolean duplicate_found = FALSE;
	for (GSList *outer_iter = list; outer_iter != NULL; outer_iter = g_slist_next(outer_iter)) {
		scroll_action *outer_action = (scroll_action *)outer_iter->data;

		// Inner loop to compare the current mouse_action with subsequent mouse_actions
		for (GSList *inner_iter = g_slist_next(outer_iter); inner_iter != NULL; inner_iter = g_slist_next(inner_iter)) {
			scroll_action *inner_action = (scroll_action *)inner_iter->data;

			// Check there are no actions duplicating button, single / double click and modifiers
			if (outer_action->state == inner_action->state &&
						outer_action->direction == inner_action->direction) {
				duplicate_found = TRUE;
				siril_log_color_message(_("Duplicate scroll_actions found with state: %d, direction: %d\n"), "red", outer_action->state, outer_action->direction);
			}
		}
	}
	gboolean retval = !duplicate_found;
	if (!retval) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid mouse scroll action configuration"), _("There are conflicts between the actions you have configured. Please check the log for details and revise your action configuration."));
	}
	return retval;
}

void update_scroll_actions_from_treeview(GtkTreeView *tree_view) {
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GtkTreeIter iter;
	gboolean valid = gtk_tree_model_get_iter_first(model, &iter);

	// Clear existing mouse actions list
	GSList *new_scroll_actions = NULL;

	// Iterate over rows in the tree view
	while (valid) {
		gchar *name;
		gchar *direction;
		gchar *modifier;
		gint reference;

		// Get data from each column
		gtk_tree_model_get(model, &iter, 0, &name, 1, &direction, 2, &modifier, 4, &reference, -1);

		// Check for un-set new action
		if (!g_strcmp0(name, NEW_ENTRY_TEXT)) {
			g_free(name);
			g_free(direction);
			g_free(modifier);
			siril_log_color_message(_("Warning: ignoring unconfigured action (\"" NEW_ENTRY_TEXT "\")\n"), "salmon");
			valid = gtk_tree_model_iter_next(model, &iter);
			continue;
		}
		// Convert direction from string to SirilScrollDirection
		SirilScrollDirection dir = MOUSE_HORIZ_SCROLL;
		if (!strcmp(direction, SCROLL_VERTICAL_TEXT))
			dir = MOUSE_VERTICAL_SCROLL;

		// Convert modifier from string to GdkModifierType
		GdkModifierType state = 0;
		if (g_strrstr(modifier, CTRL_TEXT)) {
			state |= get_primary();
		}
		if (g_strrstr(modifier, SHIFT_TEXT)) {
			state |= GDK_SHIFT_MASK;
		}

		// Convert name to metadata
		const scroll_function_metadata *metadata = map_scroll_ref_to_metadata(reference);

		// Create a new mouse action and add it to the list
		scroll_action *action = create_scroll_action(state, metadata, dir);
		siril_debug_print("New action created: %s, reference %u, modifier %u, action %u\n", action->data->name, action->data->reference, action->state, action->direction);

		new_scroll_actions = g_slist_append(new_scroll_actions, action);

		g_free(name);
		g_free(direction);
		g_free(modifier);

		valid = gtk_tree_model_iter_next(model, &iter);
	}
	if (validate_scroll_actions(new_scroll_actions)) {
		if (gui.scroll_actions) {
			g_slist_free_full(gui.scroll_actions, free);
		}
		gui.scroll_actions = new_scroll_actions;
		fill_scroll_actions_list(FALSE);
		siril_log_color_message(_("Mouse scroll configuration updated successfully.\n"), "green");
	}
}

void add_row_to_scroll_tree_view(GtkListStore *store) {
	GtkTreeIter iter;
	gtk_list_store_append(store, &iter);
	// Add default values for each column
	gtk_list_store_set(store, &iter,
					0, NEW_ENTRY_TEXT,
					1, SCROLL_VERTICAL_TEXT,
					2, NO_MODIFIER_TEXT,
					3, "",
					-1);
}

void on_scroll_actions_add_clicked(GtkButton *button, gpointer user_data) {
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("scroll_treeview_scrolled_window"));
	GtkWidget *existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "scroll_actions_treeview") == 0) {
		GtkTreeView *tree_view = GTK_TREE_VIEW(existing_tree_view);
		GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
		GtkListStore *store = GTK_LIST_STORE(model);
		add_row_to_scroll_tree_view(store);
	}
}

void on_scroll_actions_remove_clicked(GtkButton *button, gpointer user_data) {
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("scroll_treeview_scrolled_window"));
	GtkWidget *existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "scroll_actions_treeview") == 0) {
		GtkTreeView *tree_view = GTK_TREE_VIEW(existing_tree_view);
		GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
		GtkListStore *store = GTK_LIST_STORE(model);
		delete_selected_row(tree_view, store);
	}
}

void on_mouse_actions_reset_clicked(GtkButton *button, gpointer user_data) {
	if (siril_confirm_dialog(_("Reset Mouse Actions"), _("This will reset all mouse button and scroll actions to the default settings. "
							   "It is not possible to undo this and you will need to reconfigure any changes. Are you sure?"), _("Reset"))) {
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
	control_window_switch_to_tab(OUTPUT_LOGS);
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("mouse_treeview_scrolled_window"));
	GtkWidget *existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	com.pref.gui.mouse_speed_limit = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_mouse_speed_limit")));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "mouse_actions_treeview") == 0) {
		update_mouse_actions_from_treeview(GTK_TREE_VIEW(existing_tree_view));
	}
	scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("scroll_treeview_scrolled_window"));
	existing_tree_view = gtk_bin_get_child(GTK_BIN(scrolled_window));
	if (existing_tree_view && strcmp(gtk_widget_get_name(existing_tree_view), "scroll_actions_treeview") == 0) {
		update_scroll_actions_from_treeview(GTK_TREE_VIEW(existing_tree_view));
	}
	if (com.pref.gui.mouse_cfg.mouse_actions_array)
		g_slist_free_full(com.pref.gui.mouse_cfg.mouse_actions_array, g_free);
	com.pref.gui.mouse_cfg.mouse_actions_array = mouse_actions_list_to_config(gui.mouse_actions);
	if (com.pref.gui.mouse_cfg.scroll_actions_array)
		g_slist_free_full(com.pref.gui.mouse_cfg.scroll_actions_array, g_free);
	com.pref.gui.mouse_cfg.scroll_actions_array = scroll_actions_list_to_config(gui.scroll_actions);
	writeinitfile();
}

void on_mouse_actions_dialog_show(GtkDialog *dialog, gpointer user_data) {
	GtkWidget *widget = lookup_widget("mouse_test_drawingarea");
	gtk_widget_add_events(widget, GDK_SCROLL_MASK | GDK_SMOOTH_SCROLL_MASK | GDK_BUTTON_PRESS_MASK);

	// Use CSS to draw a frame around the mouse test area
	widget = lookup_widget("mouse_test_frame");
	GtkCssProvider *provider = gtk_css_provider_new();
	gtk_css_provider_load_from_data(provider, "frame { border: 2px solid #aaa; border-radius: 4px; padding: 2px; }", -1, NULL);
	GtkStyleContext *context = gtk_widget_get_style_context(widget);
	gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(provider), GTK_STYLE_PROVIDER_PRIORITY_USER);

	if (!gui.mouse_actions) {
		load_or_initialize_mouse_actions();
	}
	if (!gui.scroll_actions) {
		load_or_initialize_scroll_actions();
	}

	fill_mouse_actions_list(FALSE);
	fill_scroll_actions_list(FALSE);
}

void on_mouse_actions_dialog_hide(GtkDialog *dialog, gpointer user_data) {
	siril_close_dialog("mouse_actions_dialog");
}

void on_mouse_actions_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("mouse_actions_dialog");
}

gint timeout_ref = 0;

// Function to reset label text after 5 seconds
static gboolean reset_label_text(GtkLabel *label) {
    gtk_label_set_text(label, N_("Mouse button / scroll check"));
	timeout_ref = 0;
    return G_SOURCE_REMOVE; // Remove the timeout source after executing once
}

gboolean on_mouse_test_drawingarea_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	GtkLabel *label = GTK_LABEL(lookup_widget("mouse_test_area_label"));
	gchar *text = NULL;
	gboolean refresh_callback = FALSE;
	if (event->type == GDK_SCROLL) {
		if (timeout_ref) {
			g_source_remove(timeout_ref);
			timeout_ref = 0;
			refresh_callback = TRUE;
		}
		GdkEventScroll *scrollevent = (GdkEventScroll*) event;
		switch (scrollevent->direction) {
			case GDK_SCROLL_DOWN:
			case GDK_SCROLL_UP:
				text = g_strdup_printf(_("Scroll (vertical)"));
				break;
			case GDK_SCROLL_LEFT:
			case GDK_SCROLL_RIGHT:
				text = g_strdup_printf(_("Scroll (horizontal)"));
				break;
			case GDK_SCROLL_SMOOTH:;
				point delta;
				gdk_event_get_scroll_deltas(event, &delta.x, &delta.y);
				if (delta.y != 0.0)
					text = g_strdup_printf(_("Smooth scroll (vertical)"));
				else if (delta.x != 0.0)
					text = g_strdup_printf(_("Smooth scroll (horizontal)"));
				break;
		}
	} else if (event->type == GDK_BUTTON_PRESS) {
		GdkEventButton *buttonevent = (GdkEventButton*) event;
		text = g_strdup_printf(_("Button %u"), buttonevent->button);
	}
	if (text != NULL) {
		gtk_label_set_text(label, text);
	}
	if (text != NULL || refresh_callback) {
		timeout_ref = g_timeout_add(1250, (GSourceFunc)reset_label_text, label);
	}
	g_free(text);

	return FALSE;
}

GSList *mouse_actions_list_to_config(GSList *actions) {
	GSList *descriptions = NULL;
	GSList *current = actions;
	while (current != NULL) {
		// Get the mouse_action from the current list node
		mouse_action *action = (mouse_action *)current->data;
		if (!action->data || action->data->reference == MOUSE_REF_NULL || action->data->reference >= MOUSE_REF_MAX)
			continue;
		// Format the string; ensure data and data->name are not NULL
		gchar *description = g_strdup_printf(
			"%u,%u,%u,%u",
			action->button,
			(guint) action->type,
			(guint) action->state,
			(guint) action->data->reference
		);
		// Store the string in the array
		descriptions = g_slist_prepend(descriptions, description);
		// Move to the next element in the list
		current = g_slist_next(current);
	}
	descriptions = g_slist_reverse(descriptions);
	return descriptions;
}

GSList *mouse_actions_config_to_list(GSList *config) {
	GSList *actions = NULL;
	GSList *current = config;
	while (current != NULL) {
		// Get the mouse_action from the current list node
		gchar *input = (gchar*)current->data;
		// Split the input string by commas
		gchar **tokens = g_strsplit(input, ",", -1);
		mouse_action *action = (mouse_action*) malloc(sizeof(mouse_action));
		if (tokens[0] && tokens[1] && tokens[2] && tokens[3]) {
			// Convert the first and second parts to guint
			action->button = g_ascii_strtoull(tokens[0], NULL, 10);
			action->type = g_ascii_strtoull(tokens[1], NULL, 10);
			action->state = g_ascii_strtoull(tokens[2], NULL, 10);
			// Assign the third part to the result string, trimming leading/trailing whitespace
			mouse_function_ref reference = (mouse_function_ref) g_ascii_strtoull(tokens[3], NULL, 10);
			if (reference == MOUSE_REF_NULL || reference >= MOUSE_REF_MAX || map_ref_to_metadata(reference) == &null_action) {
				siril_log_color_message(_("Warning: when parsing mouse action config, config string parsed to an unknown function: \"%s\". Skipping..."), "salmon", input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
			action->data = map_ref_to_metadata(reference);
			if (action->data == &null_action) {
				siril_log_color_message(_("Warning: when parsing mouse action config, config string parsed to an unknown function: \"%s\". Skipping..."), "salmon", input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
		} else {
			// Handle error: input format is incorrect
			siril_log_color_message(_("Warning: when parsing mouse action config, config string has incorrect format: \"%s\". Skipping..."), "salmon", input);
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
		// Get the scroll_action from the current list node
		scroll_action *action = (scroll_action *)current->data;
		if (!action->data || action->data->reference == SCROLL_REF_NULL || action->data->reference >= SCROLL_REF_MAX)
			continue;
		// Format the string; ensure data and data->name are not NULL
		gchar *description = g_strdup_printf(
			"%u,%u,%u",
			(guint)action->direction,
			(guint)action->state,
			(guint) action->data->reference
		);

		// Store the string in the array
		descriptions = g_slist_prepend(descriptions, description);

		// Move to the next element in the list
		current = g_slist_next(current);
	}
	descriptions = g_slist_reverse(descriptions);
	return descriptions;
}

GSList *scroll_actions_config_to_list(GSList *config) {
	GSList *actions = NULL;
	GSList *current = config;
	while (current != NULL) {
		// Get the mouse_action from the current list node
		gchar *input = (gchar*)current->data;
		// Split the input string by commas
		gchar **tokens = g_strsplit(input, ",", -1);
		scroll_action *action = (scroll_action*) malloc(sizeof(scroll_action));
		if (tokens[0] && tokens[1] && tokens[2]) {
			// Convert the first and second parts to guint
			action->direction = g_ascii_strtoull(tokens[0], NULL, 10);
			action->state = g_ascii_strtoull(tokens[1], NULL, 10);
			// Assign the third part to the result string, trimming leading/trailing whitespace
			scroll_function_ref reference = (scroll_function_ref) g_ascii_strtoull(tokens[2], NULL, 10);
			if (reference == SCROLL_REF_NULL || reference >= SCROLL_REF_MAX || map_scroll_ref_to_metadata(reference) == &scroll_null_action) {
				siril_log_color_message(_("Warning: when parsing mouse action config, config string parsed to an unknown function: \"%s\". Skipping..."), "salmon", input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
			action->data = map_scroll_ref_to_metadata(reference);
			if (action->data == &scroll_null_action) {
				siril_log_color_message(_("Warning: when parsing scroll action config, string parsed to an unknown function: \"%s\". Skipping...\n"), "salmon", input);
				free(action);
				g_strfreev(tokens);
				current = g_slist_next(current);
				continue;
			}
		} else {
			// Handle error: input format is incorrect
			siril_log_color_message(_("Warning: when parsing mouse action config, string has incorrect format: \"%s\". Skipping...\n"), "salmon", input);
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
