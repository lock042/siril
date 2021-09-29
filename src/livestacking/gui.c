#include "gui.h"
#include "livestacking.h"
#include <gtk/gtk.h>
#include "gui/utils.h"
#include "gui/callbacks.h"

/* for fullscreen and window management on activation,
 * see livestacking_action_activate() in core/siril_actions.c */

void livestacking_display(const char *str) {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("livest_label1"));
	gtk_label_set_text(label, str);
}

void on_livestacking_stop_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide((GtkWidget *)user_data);
}

void on_livestacking_player_hide(GtkWidget *widget, gpointer user_data) {
	GtkWindow *window;
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	gboolean is_fullscreen;

	window = GTK_WINDOW(lookup_widget("control_window"));

	GdkWindow *gdk_window = gtk_widget_get_window(GTK_WIDGET(window));
	is_fullscreen = gdk_window_get_state(gdk_window) & GDK_WINDOW_STATE_FULLSCREEN;

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);

	stop_live_stacking_engine();
}

void update_debayer_button_status(gboolean new_state) {
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("ls_debayer"));
	gtk_toggle_button_set_active(button, new_state);
}

gboolean livestacking_first_result_idle(gpointer p) {
	set_precision_switch(); // set precision on screen
	close_tab();	// b&w to rgb gui
	return FALSE;
}
