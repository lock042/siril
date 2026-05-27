#include "livestacking/gui.h"
#include "livestacking/livestacking.h"
#include <gtk/gtk.h>
#include "io/image_format_fits.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_actions.h"
#include "core/proto.h"
#include "core/preprocess.h"
#include "core/siril_log.h"

static gchar *pause_play_button[] = {"media-playback-pause", "media-playback-start" };

/* for fullscreen and window management on activation,
 * see livestacking_action_activate() in gui/siril_actions.c */

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
	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window")));
	GAction *action_toolbar = g_action_map_lookup_action(G_ACTION_MAP(app_win), "hide-show-toolbar");
	g_action_activate(action_toolbar, NULL);
}

void force_unlinked_channels() {
	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window")));
	GAction *action_chain = g_action_map_lookup_action(G_ACTION_MAP(app_win), "chain-chan");

	GVariant *state = g_action_get_state(G_ACTION(action_chain));
	gboolean need_to_enable = g_variant_get_boolean(state);

	if (need_to_enable)
		chain_channels_state_change(G_SIMPLE_ACTION(action_chain), g_variant_new_boolean(!need_to_enable), NULL);

	g_variant_unref(state);
}

static void set_label(GtkLabel *label, gchar *text, gboolean free_after_display) {
	struct _label_struct *arg = malloc(sizeof(struct _label_struct));
	arg->label = label;
	arg->text = text;
	arg->free_text = free_after_display;
	g_idle_add(label_update_idle, arg);
}

void livestacking_display(gchar *str, gboolean free_after_display) {
	if (com.headless) {
		siril_log_message("%s\n", str);
		return;
	}
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "livest_label1")));
	set_label(label, str, free_after_display);
	livestacking_update_number_of_images(0, 0.0, -1.0, NULL);
}

static void update_icon(const gchar *name, gboolean is_loaded) {
	GtkImage *image = GTK_IMAGE(gtk_builder_get_object(gui.builder, name));
	if (is_loaded)
		gtk_image_set_from_icon_name(image, "emblem-ok-symbolic");
	else
		gtk_image_set_from_icon_name(image, "process-stop-symbolic");
}

void livestacking_display_config(gboolean use_dark, gboolean use_flat, transformation_type regtype) {
	static GtkLabel *reg_conf_label = NULL;
	if (!reg_conf_label)
		reg_conf_label = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "ls_reg_config_label")));
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
		label_cumul = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "ls_cumul_label")));
		label_stats = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "ls_stats_label")));
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

void on_livestacking_playpause_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "livest_label1"));
	widget_set_class(label, "record", "");
	gtk_button_set_icon_name(button, pause_play_button[get_paused_status()]);
	if (!livestacking_is_started()) {
		if (get_paused_status()) pause_live_stacking_engine();
		on_livestacking_start();
	} else {
		gtk_label_set_text(GTK_LABEL(label), _("Paused ..."));
		pause_live_stacking_engine();
	}
}

void on_livestacking_stop_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "livest_label1"));
	gtk_label_set_text(GTK_LABEL(label), _("Idle"));
	widget_set_class(label, "", "record");
	gtk_button_set_icon_name(GTK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "livestacking_playpause"))),
			pause_play_button[1]);
	stop_live_stacking_engine();
}

void on_livestacking_player_hide(GtkWidget *widget, gpointer user_data) {
	on_livestacking_stop_clicked(NULL, NULL);
}

gboolean update_debayer_button_status_idle(gpointer new_state) {
	GtkCheckButton *button = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "ls_debayer")));
	siril_toggle_set_active(GTK_WIDGET(button), GPOINTER_TO_INT(new_state));
	return FALSE;
}

void update_debayer_button_status(gboolean new_state) {
	siril_add_idle(update_debayer_button_status_idle, GINT_TO_POINTER(new_state));
}

gboolean livestacking_first_result_idle(gpointer p) {
	gui_function(set_precision_switch, NULL); // set precision on screen
	remap_all();
	gui_iface.redraw_image(REMAP_ALL);
	return FALSE;
}

static gboolean enable_debayer_idle(gpointer arg) {
	GtkCheckButton *button = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "demosaicingButton")));
	siril_toggle_set_active(GTK_WIDGET(button), GPOINTER_TO_INT(arg));
	return FALSE;
}

void enable_debayer(gboolean arg) {
	if (com.headless)
		com.pref.debayer.open_debayer = arg;
	else g_idle_add(enable_debayer_idle, GINT_TO_POINTER(arg));
}

gboolean end_image_loading(gpointer arg) {
	gui_iface.redraw_image(REMAP_ALL);
	return FALSE;
}

void init_preprocessing_from_GUI() {
	/* copied from core/preprocess.c test_for_master_files() but modified to not check
	 * for some options and not check for image properties against gfit, which is not
	 * already loaded here. Error management is different too.
	 * Some options also come from on_prepro_button_clicked()
	 */
	struct preprocessing_data *prepro = calloc(1, sizeof(struct preprocessing_data));

	/* checking for dark master and associated options */
	GtkCheckButton *tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "usedark_button")));
	if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
		const char *filename;
		GtkEntry *entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "darkname_entry")));
		filename = gtk_editable_get_text(GTK_EDITABLE(entry));
		if (filename[0] == '\0') {
			siril_toggle_set_active(GTK_WIDGET(tbutton), FALSE);
		} else {
			gui_iface.set_progress(PROGRESS_NONE, _("Opening dark image..."));
			prepro->dark = calloc(1, sizeof(fits));
			if (!readfits(filename, prepro->dark, NULL, FALSE)) {
				prepro->use_dark = TRUE;
			} else {
				livestacking_display(_("NOT USING DARK: cannot open the file"), FALSE);
				free(prepro->dark);
				gtk_editable_set_text(GTK_EDITABLE(entry), "");
				prepro->use_dark = FALSE;
			}
		}

		if (prepro->use_dark) {
			// cosmetic correction
			tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "cosmEnabledCheck")));
			prepro->use_cosmetic_correction = siril_toggle_get_active(GTK_WIDGET(tbutton));

			if (prepro->use_cosmetic_correction) {
				GtkCheckButton *CFA = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "cosmCFACheck")));
				prepro->is_cfa = siril_toggle_get_active(GTK_WIDGET(CFA));

				/* now we want to know which cosmetic correction was chosen */
				GtkStack *stack = GTK_STACK(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_cc")));
				GtkWidget *w = gtk_stack_get_visible_child(stack);
				if (w) {
					/* GTK4: gtk_stack_get_pages returns GtkSelectionModel*
					 * which implements GListModel — go through the iface
					 * cast for clean C typing. */
					GtkSelectionModel *pages_sel = gtk_stack_get_pages(stack);
					GListModel *pages = G_LIST_MODEL(pages_sel);
					guint n = g_list_model_get_n_items(pages);
					for (guint i = 0; i < n; i++) {
						GtkStackPage *page = g_list_model_get_item(pages, i);
						if (gtk_stack_page_get_child(page) == w) {
							prepro->cc_from_dark = (i == 0);
							g_object_unref(page);
							break;
						}
						g_object_unref(page);
					}
					g_object_unref(pages_sel);
					if (prepro->cc_from_dark) { // cosmetic correction from masterdark
						tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkSigCold")));
						if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
							GtkSpinButton *sigCold = GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spinSigCosmeColdBox")));
							prepro->sigma[0] = gtk_spin_button_get_value(sigCold);
						} else prepro->sigma[0] = -1.0;
						tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkSigHot")));
						if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
							GtkSpinButton *sigHot = GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spinSigCosmeHotBox")));
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

			GtkCheckButton *fix_xtrans = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "fix_xtrans_af")));
			prepro->fix_xtrans = siril_toggle_get_active(GTK_WIDGET(fix_xtrans));
		}
	}

	/* checking for flat master and associated options */
	tbutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "useflat_button")));
	if (siril_toggle_get_active(GTK_WIDGET(tbutton))) {
		const char *filename;
		GtkEntry *entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "flatname_entry")));
		filename = gtk_editable_get_text(GTK_EDITABLE(entry));
		if (filename[0] == '\0') {
			siril_toggle_set_active(GTK_WIDGET(tbutton), FALSE);
		} else {
			gui_iface.set_progress(PROGRESS_NONE, _("Opening flat image..."));
			prepro->flat = calloc(1, sizeof(fits));
			if (!readfits(filename, prepro->flat, NULL, !com.pref.force_16bit)) {
				prepro->use_flat = TRUE;
			} else {
				livestacking_display(_("NOT USING FLAT: cannot open the file"), FALSE);
				free(prepro->flat);
				gtk_editable_set_text(GTK_EDITABLE(entry), "");
				prepro->use_flat = FALSE;
			}

			if (prepro->use_flat) {
				GtkCheckButton *autobutton = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkbutton_auto_evaluate")));
				prepro->autolevel = siril_toggle_get_active(GTK_WIDGET(autobutton));
				if (!prepro->autolevel) {
					GtkEntry *norm_entry = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "entry_flat_norm")));
					prepro->normalisation = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(norm_entry)), NULL);
				}

				GtkCheckButton *CFA = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "cosmCFACheck")));
				prepro->is_cfa = siril_toggle_get_active(GTK_WIDGET(CFA));

				GtkCheckButton *equalize_cfa = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkbutton_equalize_cfa")));
				prepro->equalize_cfa = siril_toggle_get_active(GTK_WIDGET(equalize_cfa));
			}
		}
	}

	GtkCheckButton *use_32b_button = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "ls_32bits")));
	gboolean use_32_bits = siril_toggle_get_active(GTK_WIDGET(use_32b_button));

	init_preprocessing_finalize(prepro, use_32_bits);

	GtkCheckButton *shift_reg_button = GTK_CHECK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "ls_shiftonly")));
	init_registration_finalize(siril_toggle_get_active(GTK_WIDGET(shift_reg_button)));
}

void on_livestacking_start() {
	init_preprocessing_from_GUI();
	if (start_livestacking(TRUE))
		stop_live_stacking_engine();
}
