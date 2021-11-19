#include "gui.h"
#include "livestacking.h"
#include <gtk/gtk.h>
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "core/proto.h"

/* for fullscreen and window management on activation,
 * see livestacking_action_activate() in core/siril_actions.c */

// TODO: do in idle function
void livestacking_display(const char *str) {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("livest_label1"));
	gtk_label_set_text(label, str);
}

void livestacking_display_config(gboolean use_dark, transformation_type regtype) {
	static GtkLabel *conf_label = NULL;
	if (!conf_label)
		conf_label = GTK_LABEL(lookup_widget("ls_config_label"));
	gchar * txt = g_strdup_printf("%s dark and cosmetic correction, registering with %s transformation, stacking with sum", use_dark ? "Using" : "Not using", describe_transformation_type(regtype));
	gtk_label_set_text(conf_label, txt);
	g_free(txt);
}

// TODO: do in idle function
void livestacking_update_number_of_images(int nb, double total_exposure) {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("ls_cumul_label"));
	int secs = round_to_int(total_exposure);
	int time;
	const char *unit;
	if (secs >= 180) {
		time = secs / 60;
		unit = "minutes";
	} else {
		time = secs;
		unit = "seconds";
	}
	gchar *txt = g_strdup_printf("%d images stacked, %d %s of cumulated exposure", nb, time, unit);
	gtk_label_set_text(label, txt);
	g_free(txt);
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
