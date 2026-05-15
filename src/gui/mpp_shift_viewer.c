/*
 * Diagnostic shift viewer — Phase 9 slice 21.
 *
 * After a successful Register the cached mpp_run_t holds per-AP per-frame
 * shifts (run->shifts). This non-modal dialog lets the user pick a frame
 * and visualise those shifts as arrows on top of the AP overlay so
 * they can spot APs with anomalous or failed Stage B fits.
 *
 * Read-only: changes nothing about the run.
 */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/mpp_shift_viewer.h"
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_shift.h"

static GtkWidget       *dialog       = NULL;
static GtkSpinButton   *spin_frame   = NULL;
static GtkSpinButton   *spin_scale   = NULL;
static GtkAdjustment   *adj_frame    = NULL;
static GtkLabel        *label_stats  = NULL;

static void init_statics(void) {
	if (dialog) return;
	dialog       = GTK_WIDGET(gtk_builder_get_object(gui.builder, "mpp_shift_viewer_dialog"));
	spin_frame   = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_shift_viewer_frame"));
	spin_scale   = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_shift_viewer_scale"));
	adj_frame    = GTK_ADJUSTMENT (gtk_builder_get_object(gui.builder, "adj_mpp_shift_viewer_frame"));
	label_stats  = GTK_LABEL      (gtk_builder_get_object(gui.builder, "label_mpp_shift_viewer_stats"));
}

gboolean mpp_shift_viewer_is_open(void) {
	return dialog && gtk_widget_get_visible(dialog);
}

int mpp_shift_viewer_get_frame(void) {
	if (!spin_frame) return -1;
	return gtk_spin_button_get_value_as_int(spin_frame) - 1;   /* 1-based UI, 0-based internal */
}

double mpp_shift_viewer_get_scale(void) {
	if (!spin_scale) return 1.0;
	return gtk_spin_button_get_value(spin_scale);
}

void mpp_shift_viewer_update_button_sensitivity(void) {
	static GtkWidget *btn = NULL;
	if (!btn) {
		btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "button_mpp_view_shifts"));
		if (!btn) return;
	}
	const mpp_run_t *run = mpp_get_cached_run();
	gtk_widget_set_sensitive(btn,
		run != NULL && run->shifts != NULL && run->shifts->num_frames > 0);
}

/* Update the "frame X / Y: A APs ok, B failed" status line. */
static void refresh_stats_label(void) {
	if (!label_stats) return;
	const mpp_run_t *run = mpp_get_cached_run();
	if (!run || !run->shifts || !run->shifts->success) {
		gtk_label_set_text(label_stats, _("No shift data."));
		return;
	}
	const int f = mpp_shift_viewer_get_frame();
	const int M = run->shifts->num_aps;
	int ok = 0;
	if (f >= 0 && f < run->shifts->num_frames) {
		for (int a = 0; a < M; ++a)
			if (run->shifts->success[(size_t) f * M + a]) ++ok;
	}
	gchar *txt = g_strdup_printf(_("Frame %d / %d: %d / %d APs converged"),
	                             f + 1, run->shifts->num_frames, ok, M);
	gtk_label_set_text(label_stats, txt);
	g_free(txt);
}

void on_seqmpp_view_shifts_button_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	mpp_run_t *run = mpp_get_cached_run();
	if (!run || !run->shifts) {
		siril_log_color_message(_("Shift viewer: no Stage B shifts cached. "
		                          "Run Register first.\n"), "salmon");
		return;
	}
	init_statics();
	if (!dialog) return;
	if (adj_frame)
		gtk_adjustment_set_upper(adj_frame, (gdouble) run->shifts->num_frames);
	if (spin_frame) {
		gdouble v = gtk_spin_button_get_value(spin_frame);
		if (v < 1.0 || v > (gdouble) run->shifts->num_frames)
			gtk_spin_button_set_value(spin_frame, 1.0);
	}
	GtkWindow *parent = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	if (parent)
		gtk_window_set_transient_for(GTK_WINDOW(dialog), parent);
	gtk_widget_show(dialog);
}

gboolean on_mpp_shift_viewer_dialog_delete_event(GtkWidget *widget, GdkEvent *event,
                                                 gpointer user_data) {
	(void) event; (void) user_data;
	gtk_widget_hide(widget);
	return TRUE;
}

void on_mpp_shift_viewer_dialog_show(GtkWidget *widget, gpointer user_data) {
	(void) widget; (void) user_data;
	init_statics();
	refresh_stats_label();
	redraw(REDRAW_OVERLAY);
}

void on_mpp_shift_viewer_dialog_hide(GtkWidget *widget, gpointer user_data) {
	(void) widget; (void) user_data;
	redraw(REDRAW_OVERLAY);
}

void on_mpp_shift_viewer_frame_changed(GtkSpinButton *spin, gpointer user_data) {
	(void) spin; (void) user_data;
	refresh_stats_label();
	redraw(REDRAW_OVERLAY);
}

void on_mpp_shift_viewer_scale_changed(GtkSpinButton *spin, gpointer user_data) {
	(void) spin; (void) user_data;
	redraw(REDRAW_OVERLAY);
}

void on_mpp_shift_viewer_close_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (dialog) gtk_widget_hide(dialog);
}
