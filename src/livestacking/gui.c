#include "gui.h"
#include "livestacking.h"
#include <gtk/gtk.h>
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "core/siril_actions.h"
#include "core/proto.h"

static gchar *pause_play_button[] = {"media-playback-pause", "media-playback-start" };

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

void show_hide_toolbox() {
	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
	GAction *action_toolbar = g_action_map_lookup_action(G_ACTION_MAP(app_win), "hide-show-toolbar");
	g_action_activate(action_toolbar, NULL);
}

void force_unlinked_channels() {
	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
	GAction *action_chain = g_action_map_lookup_action(G_ACTION_MAP(app_win), "chain-chan");

	GVariant *state = g_action_get_state(G_ACTION(action_chain));
	gboolean need_to_enable = g_variant_get_boolean(state);

	if (need_to_enable)
		chain_channels_state_change(G_SIMPLE_ACTION(action_chain), g_variant_new_boolean(!need_to_enable), NULL);

	g_variant_unref(state);
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
	livestacking_update_number_of_images(0, 0.0, -1.0);
}

void livestacking_display_config(gboolean use_dark, gboolean use_flat, transformation_type regtype) {
	static GtkLabel *conf_label = NULL;
	if (!conf_label)
		conf_label = GTK_LABEL(lookup_widget("ls_config_label"));
	gchar * txt = g_strdup_printf(_("%s dark and cosmetic correction\n\n%s flat\n\n"
			"Registering with %s transformation\n\n"
			"Stacking with weighted mean"),
			use_dark ? _("Using") : _("Not using"),
			use_flat ? _("Using") : _("Not using"),
			describe_transformation_type(regtype));
	gtk_label_set_text(conf_label, txt);
	g_free(txt);
}

void livestacking_update_number_of_images(int nb, double total_exposure, double noise) {
	static GtkLabel *label_cumul = NULL, *label_stats = NULL;
	if (!label_cumul) {
		label_cumul = GTK_LABEL(lookup_widget("ls_cumul_label"));
		label_stats = GTK_LABEL(lookup_widget("ls_stats_label"));
	}
	int secs = round_to_int(total_exposure);
	int time;
	const char *unit;
	if (secs >= 180) {
		time = secs / 60;
		unit = _("minutes");
	} else {
		time = secs;
		unit = _("seconds");
	}
	gchar *txt;
	if (noise > 0.0) {
		txt = g_strdup_printf(_("noise: %0.3f"), noise);
		set_label(label_stats, txt, TRUE);
	}

	txt = g_strdup_printf(_("%d images stacked\n\n%d %s of cumulated exposure"), nb, time, unit);
	set_label(label_cumul, txt, TRUE);
}

void on_livestacking_player_hide(GtkWidget *widget, gpointer user_data) {
	stop_live_stacking_engine();
}

void on_livestacking_stop_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("livestacking_player"));
}

void on_livestacking_pause_clicked(GtkToolButton *button, gpointer user_data) {
	pause_live_stacking_engine();
	gtk_tool_button_set_icon_name(button, pause_play_button[get_paused_status()]);
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

