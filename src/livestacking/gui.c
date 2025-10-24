#include "gui.h"
#include "livestacking.h"
#include <gtk/gtk.h>
#include "io/image_format_fits.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "core/siril_actions.h"
#include "core/proto.h"
#include "core/preprocess.h"
#include "core/siril_log.h"

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
	siril_free(arg);
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

static void set_label(GtkLabel *label, gchar *text, gboolean free_after_display) {
	struct _label_struct *arg = siril_malloc(sizeof(struct _label_struct));
	arg->label = label;
	arg->text = text;
	arg->free_text = free_after_display;
	gdk_threads_add_idle(label_update_idle, arg);
}

void livestacking_display(gchar *str, gboolean free_after_display) {
	if (com.headless) {
		siril_log_message("%s\n", str);
		return;
	}
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("livest_label1"));
	set_label(label, str, free_after_display);
	livestacking_update_number_of_images(0, 0.0, -1.0, NULL);
}

static void update_icon(const gchar *name, gboolean is_loaded) {
	GtkImage *image = GTK_IMAGE(lookup_widget(name));
	if (is_loaded)
		gtk_image_set_from_icon_name(image, "gtk-yes", GTK_ICON_SIZE_LARGE_TOOLBAR);
	else
		gtk_image_set_from_icon_name(image, "gtk-no", GTK_ICON_SIZE_LARGE_TOOLBAR);
}

void livestacking_display_config(gboolean use_dark, gboolean use_flat, transformation_type regtype) {
	static GtkLabel *reg_conf_label = NULL;
	if (!reg_conf_label)
		reg_conf_label = GTK_LABEL(lookup_widget("ls_reg_config_label"));
	const char *desc = describe_transformation_type(regtype);
	update_icon("ls_reg_config", g_strcmp0(desc, "INVALID"));

	gchar *txt = g_strdup_printf(_("Registering with %s transformation"), desc);
	gtk_label_set_text(reg_conf_label, txt);
	g_free(txt);

	update_icon("ls_masterdark", use_dark);
	update_icon("ls_masterflat", use_flat);
}

void livestacking_update_number_of_images(int nb, double total_exposure, double noise, const char *process_time) {
	if (com.headless)
		return;
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
	if (noise > 0.0 && process_time) {
		txt = g_strdup_printf(_("Noise: %0.3f\n\nLast image processing time: %s"), noise, process_time);
		set_label(label_stats, txt, TRUE);
	} else if (noise > 0.0) {
		txt = g_strdup_printf(_("Noise: %0.3f"), noise);
		set_label(label_stats, txt, TRUE);
	} else {
		set_label(label_stats, _("Noise: N/A"), FALSE);
	}

	txt = g_strdup_printf(_("%d images stacked\n\n%d %s of cumulative exposure"), nb, time, unit);
	set_label(label_cumul, txt, TRUE);
}

void on_livestacking_playpause_clicked(GtkToolButton *button, gpointer user_data) {
	GtkWidget *label = lookup_widget("livest_label1");
	widget_set_class(label, "record", "");
	gtk_tool_button_set_icon_name(button, pause_play_button[get_paused_status()]);
	if (!livestacking_is_started()) {
		if (get_paused_status()) pause_live_stacking_engine();
		on_livestacking_start();
	} else {
		gtk_label_set_text(GTK_LABEL(label), _("Paused ..."));
		pause_live_stacking_engine();
	}
}

void on_livestacking_stop_clicked(GtkToolButton *button, gpointer user_data) {
	GtkWidget *label = lookup_widget("livest_label1");
	gtk_label_set_text(GTK_LABEL(label), _("Idle"));
	widget_set_class(label, "", "record");
	gtk_tool_button_set_icon_name(GTK_TOOL_BUTTON(lookup_widget("livestacking_playpause")),
			pause_play_button[1]);
	stop_live_stacking_engine();
}

void on_livestacking_player_hide(GtkWidget *widget, gpointer user_data) {
	on_livestacking_stop_clicked(NULL, NULL);
}

gboolean update_debayer_button_status_idle(gpointer new_state) {
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("ls_debayer"));
	gtk_toggle_button_set_active(button, GPOINTER_TO_INT(new_state));
	return FALSE;
}

void update_debayer_button_status(gboolean new_state) {
	siril_add_idle(update_debayer_button_status_idle, GINT_TO_POINTER(new_state));
}

gboolean livestacking_first_result_idle(gpointer p) {
	gui_function(set_precision_switch, NULL); // set precision on screen
	redraw(REMAP_ALL);
	return FALSE;
}

static gboolean enable_debayer_idle(gpointer arg) {
	GtkToggleButton *button = GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton"));
	gtk_toggle_button_set_active(button, GPOINTER_TO_INT(arg));
	return FALSE;
}

void enable_debayer(gboolean arg) {
	if (com.headless)
		com.pref.debayer.open_debayer = arg;
	else gdk_threads_add_idle(enable_debayer_idle, GINT_TO_POINTER(arg));
}

gboolean end_image_loading(gpointer arg) {
	redraw(REMAP_ALL);
	return FALSE;
}

void init_preprocessing_from_GUI() {
	/* copied from core/preprocess.c test_for_master_files() but modified to not check
	 * for some options and not check for image properties against gfit, which is not
	 * already loaded here. Error management is different too.
	 * Some options also come from on_prepro_button_clicked()
	 */
	struct preprocessing_data *prepro = siril_calloc(1, sizeof(struct preprocessing_data));

	/* checking for dark master and associated options */
	GtkToggleButton *tbutton = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		GtkEntry *entry = GTK_ENTRY(lookup_widget("darkname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			set_progress_bar_data(_("Opening dark image..."), PROGRESS_NONE);
			prepro->dark = siril_calloc(1, sizeof(fits));
			if (!readfits(filename, prepro->dark, NULL, FALSE)) {
				prepro->use_dark = TRUE;
			} else {
				livestacking_display(_("NOT USING DARK: cannot open the file"), FALSE);
				siril_free(prepro->dark);
				gtk_entry_set_text(entry, "");
				prepro->use_dark = FALSE;
			}
		}

		if (prepro->use_dark) {
			// cosmetic correction
			tbutton = GTK_TOGGLE_BUTTON(lookup_widget("cosmEnabledCheck"));
			prepro->use_cosmetic_correction = gtk_toggle_button_get_active(tbutton);

			if (prepro->use_cosmetic_correction) {
				GtkToggleButton *CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
				prepro->is_cfa = gtk_toggle_button_get_active(CFA);

				/* now we want to know which cosmetic correction was chosen */
				GtkStack *stack = GTK_STACK(lookup_widget("stack_cc"));
				GtkWidget *w = gtk_stack_get_visible_child(stack);
				if (w) {
					GValue value = G_VALUE_INIT;
					g_value_init(&value, G_TYPE_INT);
					gtk_container_child_get_property(GTK_CONTAINER(stack), w, "position", &value);
					gint position = g_value_get_int(&value);
					prepro->cc_from_dark = position == 0 ? TRUE : FALSE;
					g_value_unset(&value);
					if (prepro->cc_from_dark) { // cosmetic correction from masterdark
						tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkSigCold"));
						if (gtk_toggle_button_get_active(tbutton)) {
							GtkSpinButton *sigCold = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox"));
							prepro->sigma[0] = gtk_spin_button_get_value(sigCold);
						} else prepro->sigma[0] = -1.0;
						tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkSigHot"));
						if (gtk_toggle_button_get_active(tbutton)) {
							GtkSpinButton *sigHot = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox"));
							prepro->sigma[1] = gtk_spin_button_get_value(sigHot);
						} else prepro->sigma[1] = -1.0;
					} else {
						// TODO: cosmetic correction from bap pixel map is disabled for now
						// may need to take it out from the use_dark condition anyway
						prepro->use_cosmetic_correction = FALSE;
						livestacking_display(_("COSMETIC CORRECTION: can only use masterdark in livestacking mode - disabling"), FALSE);
						prepro->bad_pixel_map_file = NULL;
					}
				} else {  //fallback required?
					prepro->sigma[0] = -1.0;
					prepro->sigma[1] = -1.0;
					prepro->use_cosmetic_correction = FALSE;
					livestacking_display(_("COSMETIC CORRECTION: could not read the inputs - disabling"), FALSE);
				}
			}

			GtkToggleButton *fix_xtrans = GTK_TOGGLE_BUTTON(lookup_widget("fix_xtrans_af"));
			prepro->fix_xtrans = gtk_toggle_button_get_active(fix_xtrans);
		}
	}

	/* checking for flat master and associated options */
	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		GtkEntry *entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
			prepro->flat = siril_calloc(1, sizeof(fits));
			if (!readfits(filename, prepro->flat, NULL, !com.pref.force_16bit)) {
				prepro->use_flat = TRUE;
			} else {
				livestacking_display(_("NOT USING FLAT: cannot open the file"), FALSE);
				siril_free(prepro->flat);
				gtk_entry_set_text(entry, "");
				prepro->use_flat = FALSE;
			}

			if (prepro->use_flat) {
				GtkToggleButton *autobutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto_evaluate"));
				prepro->autolevel = gtk_toggle_button_get_active(autobutton);
				if (!prepro->autolevel) {
					GtkEntry *norm_entry = GTK_ENTRY(lookup_widget("entry_flat_norm"));
					prepro->normalisation = g_ascii_strtod(gtk_entry_get_text(norm_entry), NULL);
				}

				GtkToggleButton *CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
				prepro->is_cfa = gtk_toggle_button_get_active(CFA);

				GtkToggleButton *equalize_cfa = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
				prepro->equalize_cfa = gtk_toggle_button_get_active(equalize_cfa);
			}
		}
	}

	GtkToggleButton *use_32b_button = GTK_TOGGLE_BUTTON(lookup_widget("ls_32bits"));
	gboolean use_32_bits = gtk_toggle_button_get_active(use_32b_button);

	init_preprocessing_finalize(prepro, use_32_bits);

	GtkToggleButton *shift_reg_button = GTK_TOGGLE_BUTTON(lookup_widget("ls_shiftonly"));
	init_registration_finalize(gtk_toggle_button_get_active(shift_reg_button));
}

void on_livestacking_start() {
	init_preprocessing_from_GUI();
	if (start_livestacking(TRUE))
		stop_live_stacking_engine();
}

