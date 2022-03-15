#include "gui.h"
#include "livestacking.h"
#include <gtk/gtk.h>
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "core/proto.h"

/* for fullscreen and window management on activation,
 * see livestacking_action_activate() in core/siril_actions.c */

struct _label_struct {
	GtkLabel *label;
	gchar *text;
	gboolean free_text;
};

static gboolean label_update_idle(gpointer ptr) {
	struct _label_struct *arg = (struct _label_struct *)ptr;
	gtk_label_set_text(arg->label, arg->text);
	if (arg->free_text)
		g_free(arg->text);
	free(arg);
	return FALSE;
}

void set_label(GtkLabel *label, gchar *text, gboolean free_after_display) {
	struct _label_struct *arg = malloc(sizeof(struct _label_struct));
	arg->label = label;
	arg->text = text;
	arg->free_text = free_after_display;
	gdk_threads_add_idle(label_update_idle, arg);
}

void livestacking_display(gchar *str, gboolean free_after_display) {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("livest_label1"));
	set_label(label, str, free_after_display);
}

void livestacking_display_config(gboolean use_dark, gboolean use_flat, transformation_type regtype) {
	static GtkLabel *conf_label = NULL;
	if (!conf_label)
		conf_label = GTK_LABEL(lookup_widget("ls_config_label"));
	gchar * txt = g_strdup_printf("%s dark and cosmetic correction, %s flat, registering with %s transformation, stacking with weighted mean",
			use_dark ? "Using" : "Not using",
			use_flat ? "Using" : "Not using",
			describe_transformation_type(regtype));
	gtk_label_set_text(conf_label, txt);
	g_free(txt);
}

void livestacking_update_number_of_images(int nb, double total_exposure, double noise) {
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
	gchar *txt;
	if (noise > 0.0)
		txt = g_strdup_printf("%d images stacked, %d %s of cumulated exposure, noise: %0.3f", nb, time, unit, noise);
	else txt = g_strdup_printf("%d images stacked, %d %s of cumulated exposure", nb, time, unit);
	set_label(label, txt, TRUE);
}

void on_livestacking_stop_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide((GtkWidget *)user_data);
}

void on_livestacking_pause_clicked(GtkButton *button, gpointer user_data) {
	pause_live_stacking_engine();
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

gboolean update_debayer_button_status_idle(gpointer new_state) {
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("ls_debayer"));
	gtk_toggle_button_set_active(button, GPOINTER_TO_INT(new_state));
	return FALSE;
}

void update_debayer_button_status(gboolean new_state) {
	gdk_threads_add_idle(update_debayer_button_status_idle, GINT_TO_POINTER(new_state));
}

gboolean livestacking_first_result_idle(gpointer p) {
	set_precision_switch(); // set precision on screen
	return FALSE;
}

static gboolean enable_debayer_idle(gpointer arg) {
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton"));
	gtk_toggle_button_set_active(button, GPOINTER_TO_INT(arg));
	return FALSE;
}

void enable_debayer(gboolean arg) {
	gdk_threads_add_idle(enable_debayer_idle, GINT_TO_POINTER(arg));
}

static gboolean end_image_loading(gpointer arg) {
	redraw(REMAP_ALL);
	return FALSE;
}

void complete_image_loading() {
	execute_idle_and_wait_for_it(end_image_loading, NULL);
}

