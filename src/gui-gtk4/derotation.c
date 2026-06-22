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

/*
 * Planetary derotation tool window. Reads the loaded sequence's per-frame
 * timestamps, proposes the capture-midpoint reference epoch, lets the user fit
 * the planet disk (auto-detect + numeric tweak, with a live outline drawn on
 * the image), and writes a <seqname>.derot sidecar via the shared builder. The
 * window does not stack: registration and stacking then run as usual and pick
 * up the sidecar automatically.
 */

#include <math.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/proto.h"
#include "io/sequence.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"   /* mouse_status, MOUSE_ACTION_FIT_DISK */
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/derotation.h"

#include "algos/planet_ephem.h"
#include "registration/mpp/mpp_derot_sidecar.h"
#include "registration/mpp/mpp_derot_build.h"

static GtkWidget *dialog = NULL;
static GtkDropDown *dd_body = NULL, *dd_system = NULL;
static GtkEntry *entry_epoch = NULL;
static GtkSpinButton *spin_cx = NULL, *spin_cy = NULL, *spin_radius = NULL,
                     *spin_rpol = NULL, *spin_pa = NULL, *spin_fps = NULL;
static GtkCheckButton *chk_parity = NULL;
static GtkLabel *status = NULL;
static gboolean window_open = FALSE;
static int g_drag_mode = 0;   /* 0 none, 1 move centre, 2 equatorial, 3 polar */

/* Known geometric flattening per body, used to default the polar radius. */
static double body_flattening(int body) {
	switch (body) {
		case 1:  return 0.09796;   /* Saturn */
		case 2:  return 0.00589;   /* Mars */
		default: return 0.06487;   /* Jupiter */
	}
}

/* Set the polar radius from the equatorial radius and the body's flattening. */
static void set_polar_default(void) {
	if (!spin_radius || !spin_rpol || !dd_body) return;
	const double req = gtk_spin_button_get_value(spin_radius);
	gtk_spin_button_set_value(spin_rpol,
	    req * (1.0 - body_flattening((int) gtk_drop_down_get_selected(dd_body))));
}

/* Default rotation system per body (matches the command): Jupiter -> II,
 * Saturn -> III, Mars -> I. Also refresh the polar-radius default. */
static void on_derot_body_changed(GObject *obj, GParamSpec *pspec, gpointer u) {
	(void) obj; (void) pspec; (void) u;
	if (!dd_system) return;
	const guint body = gtk_drop_down_get_selected(dd_body);
	const guint sys = (body == 0) ? 1u : (body == 1) ? 2u : 0u;   /* Jup II, Sat III, Mars I */
	gtk_drop_down_set_selected(dd_system, sys);
	set_polar_default();
	if (window_open) redraw(REDRAW_OVERLAY);
}

static void init_statics(void) {
	if (dialog) return;
	dialog      = GTK_WIDGET      (gtk_builder_get_object(gui.builder, "derotation_dialog"));
	dd_body     = GTK_DROP_DOWN   (gtk_builder_get_object(gui.builder, "derot_body"));
	dd_system   = GTK_DROP_DOWN   (gtk_builder_get_object(gui.builder, "derot_system"));
	entry_epoch = GTK_ENTRY       (gtk_builder_get_object(gui.builder, "derot_epoch"));
	spin_cx     = GTK_SPIN_BUTTON (gtk_builder_get_object(gui.builder, "derot_cx"));
	spin_cy     = GTK_SPIN_BUTTON (gtk_builder_get_object(gui.builder, "derot_cy"));
	spin_radius = GTK_SPIN_BUTTON (gtk_builder_get_object(gui.builder, "derot_radius"));
	spin_rpol   = GTK_SPIN_BUTTON (gtk_builder_get_object(gui.builder, "derot_rpol"));
	spin_pa     = GTK_SPIN_BUTTON (gtk_builder_get_object(gui.builder, "derot_pa"));
	spin_fps    = GTK_SPIN_BUTTON (gtk_builder_get_object(gui.builder, "derot_fps"));
	chk_parity  = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "derot_parity"));
	status      = GTK_LABEL       (gtk_builder_get_object(gui.builder, "derot_status"));
	if (dd_body)
		g_signal_connect(dd_body, "notify::selected",
		                 G_CALLBACK(on_derot_body_changed), NULL);
}

static void set_status(const char *msg) {
	if (status) gtk_label_set_text(status, msg);
}

/* Fill the epoch entry with the capture midpoint, and auto-detect the disk. */
static void populate_from_sequence(void) {
	if (!sequence_is_loaded() || com.seq.number <= 0) {
		set_status(_("Load a planetary sequence first."));
		return;
	}
	const int N = com.seq.number;
	double *jd = malloc((size_t) N * sizeof(double));
	if (jd && mpp_derot_frame_times(&com.seq, 0.0, NAN, jd)) {
		/* Build the midpoint datetime from the J2000 anchor (full JD 2451545.0):
		 * Julian_to_date_time() expects MJD, not the full JD that
		 * date_time_to_Julian()/the .derot use, so using it here would land the
		 * default epoch thousands of years off. */
		GDateTime *j2000 = g_date_time_new_utc(2000, 1, 1, 12, 0, 0);
		GDateTime *mid = j2000
		    ? g_date_time_add_seconds(j2000,
		          (mpp_derot_midpoint_epoch(jd, N) - 2451545.0) * 86400.0)
		    : NULL;
		if (j2000) g_date_time_unref(j2000);
		if (mid) {
			gchar *s = date_time_to_FITS_date(mid);
			if (s) { gtk_editable_set_text(GTK_EDITABLE(entry_epoch), s); g_free(s); }
			g_date_time_unref(mid);
		}
		set_status(_("Fit the disk to the planet, then compute."));
	} else {
		set_status(_("No per-frame timestamps — set the frame rate (fps)."));
	}
	free(jd);

	double cx, cy, r;
	if (gfit && mpp_derot_autodetect_disk(gfit, &cx, &cy, &r)) {
		gtk_spin_button_set_value(spin_cx, cx);
		gtk_spin_button_set_value(spin_cy, cy);
		gtk_spin_button_set_value(spin_radius, r);
		set_polar_default();
	}
}

void on_seqmpp_derotation_button_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (!sequence_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No sequence"),
		    _("Load the planetary sequence (video) before setting up derotation."));
		return;
	}
	init_statics();
	if (!dialog) return;
	GtkWindow *parent = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	if (parent)
		gtk_window_set_transient_for(GTK_WINDOW(dialog), parent);
	gtk_window_present(GTK_WINDOW(dialog));
}

void on_derotation_dialog_show(GtkWidget *widget, gpointer user_data) {
	(void) widget; (void) user_data;
	init_statics();
	window_open = TRUE;
	mouse_status = MOUSE_ACTION_FIT_DISK;   /* clicks on the image now fit the disc */
	on_derot_body_changed(NULL, NULL, NULL);   /* sync system to the body default */
	populate_from_sequence();
	redraw(REDRAW_OVERLAY);
}

static void leave_fit_mode(void) {
	window_open = FALSE;
	g_drag_mode = 0;
	if (mouse_status == MOUSE_ACTION_FIT_DISK)
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
}

gboolean on_derotation_dialog_close_request(GtkWindow *win, gpointer user_data) {
	(void) user_data;
	leave_fit_mode();
	gtk_widget_set_visible(GTK_WIDGET(win), FALSE);
	redraw(REDRAW_OVERLAY);
	return TRUE;   /* hide, don't destroy */
}

void on_derotation_close_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (dialog) {
		leave_fit_mode();
		gtk_widget_set_visible(dialog, FALSE);
		redraw(REDRAW_OVERLAY);
	}
}

void on_derotation_value_changed(GtkSpinButton *spin, gpointer user_data) {
	(void) spin; (void) user_data;
	if (window_open) redraw(REDRAW_OVERLAY);
}

void on_derotation_autofit_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	double cx, cy, r;
	if (gfit && mpp_derot_autodetect_disk(gfit, &cx, &cy, &r)) {
		gtk_spin_button_set_value(spin_cx, cx);
		gtk_spin_button_set_value(spin_cy, cy);
		gtk_spin_button_set_value(spin_radius, r);
		set_polar_default();
		set_status(_("Disk auto-detected — adjust if needed, then compute."));
	} else {
		set_status(_("Could not auto-detect a disk; set the centre and radius manually."));
	}
}

void on_derotation_compute_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	init_statics();
	if (!sequence_is_loaded() || com.seq.number <= 0 || !com.seq.seqname) {
		set_status(_("No sequence loaded."));
		return;
	}
	const int N = com.seq.number;
	const planet_body_t body = (planet_body_t) gtk_drop_down_get_selected(dd_body);
	const int system = (int) gtk_drop_down_get_selected(dd_system) + 1;
	const double cx = gtk_spin_button_get_value(spin_cx);
	const double cy = gtk_spin_button_get_value(spin_cy);
	const double radius = gtk_spin_button_get_value(spin_radius);
	const double rpol = gtk_spin_button_get_value(spin_rpol);
	const double pa = gtk_spin_button_get_value(spin_pa);
	const double parity = gtk_check_button_get_active(chk_parity) ? -1.0 : 1.0;
	const double fps = gtk_spin_button_get_value(spin_fps);

	if (radius < 1.0) { set_status(_("Set a valid disk radius.")); return; }

	double *jd = malloc((size_t) N * sizeof(double));
	if (!jd) { set_status(_("Out of memory.")); return; }
	if (!mpp_derot_frame_times(&com.seq, fps, NAN, jd)) {
		set_status(_("No per-frame timestamps — set the frame rate (fps)."));
		free(jd);
		return;
	}

	double epoch_jd;
	const char *eptxt = gtk_editable_get_text(GTK_EDITABLE(entry_epoch));
	GDateTime *ep = (eptxt && *eptxt) ? FITS_date_to_date_time((gchar *) eptxt) : NULL;
	if (ep) { epoch_jd = date_time_to_Julian(ep); g_date_time_unref(ep); }
	else      epoch_jd = mpp_derot_midpoint_epoch(jd, N);

	mpp_derot_t *d = mpp_derot_build(body, system, epoch_jd, jd, N,
	                                 (int) com.seq.ry, (int) com.seq.rx,
	                                 cx, cy, radius, pa, parity, NAN, NAN, NAN);
	free(jd);
	if (!d) { set_status(_("Ephemeris/plan build failed.")); return; }
	/* Let the user's polar radius set the disk oblateness (overriding the
	 * ephemeris default) so the fitted ellipse is used as drawn. */
	if (rpol > 0.0 && rpol < radius)
		d->flattening = 1.0 - rpol / radius;

	gchar *path = g_strdup_printf("%s.derot", com.seq.seqname);
	const mpp_status_t rc = mpp_derot_write(path, d);
	const double span_min = (d->jd[N - 1] - d->jd[0]) * 24.0 * 60.0;
	if (rc == MPP_OK) {
		gchar *msg = g_strdup_printf(_("Wrote %s — %d frames, %.2f min span. "
		                               "Now register and stack the sequence."),
		                             path, N, span_min);
		set_status(msg);
		siril_log_message(_("derotation: wrote %s (%d frames, %.2f min span, "
		                    "System %d)\n"), path, N, span_min, system);
		g_free(msg);
	} else {
		gchar *msg = g_strdup_printf(_("Failed to write %s."), path);
		set_status(msg);
		g_free(msg);
	}
	g_free(path);
	mpp_derot_free(d);
}

/* ---- overlay query ---- */

gboolean derotation_is_open(void) { return window_open; }

int derotation_get_body(void) {
	return (window_open && dd_body) ? (int) gtk_drop_down_get_selected(dd_body) : -1;
}

/* ---- interactive fitting (drag the disc to move, drag handles to resize) ---- */

/* Current fit in display coords (image pixels, y-down to match mouse events):
 * centre (cx, ycd) plus the equator (eu) and pole (pu) unit directions. */
static gboolean fit_display_geom(double *cx, double *ycd, double *req,
                                 double *rpol, double *eux, double *euy,
                                 double *pux, double *puy) {
	if (!window_open || !spin_cx || !gfit || gfit->ry <= 0) return FALSE;
	*cx  = gtk_spin_button_get_value(spin_cx);
	*ycd = (double) (gfit->ry - 1) - gtk_spin_button_get_value(spin_cy);
	*req  = gtk_spin_button_get_value(spin_radius);
	*rpol = gtk_spin_button_get_value(spin_rpol);
	const double pa = gtk_spin_button_get_value(spin_pa) * M_PI / 180.0;
	const double ex = gtk_check_button_get_active(chk_parity) ? -1.0 : 1.0;
	const double Px = ex * sin(pa), Py = -cos(pa);   /* pole unit (cairo/display) */
	*pux = Px; *puy = Py; *eux = -Py; *euy = Px;
	return TRUE;
}

/* Hit-test a display-space click: returns the drag mode, or 0 for a miss. */
int derotation_hit_test(double dx, double dy) {
	double cx, ycd, req, rpol, eux, euy, pux, puy;
	if (!fit_display_geom(&cx, &ycd, &req, &rpol, &eux, &euy, &pux, &puy))
		return 0;
	const double tol = fmax(8.0, req * 0.12);
	const double hx[4] = { cx + req * eux, cx - req * eux, cx + rpol * pux, cx - rpol * pux };
	const double hy[4] = { ycd + req * euy, ycd - req * euy, ycd + rpol * puy, ycd - rpol * puy };
	const int hmode[4] = { 2, 2, 3, 3 };
	for (int i = 0; i < 4; ++i)
		if (hypot(dx - hx[i], dy - hy[i]) <= tol) return hmode[i];
	if (hypot(dx - cx, dy - ycd) <= req) return 1;   /* inside the disc -> move */
	return 0;
}

void derotation_set_drag(int mode) { g_drag_mode = mode; }
int  derotation_get_drag(void)     { return g_drag_mode; }

/* Apply a drag to a display-space point, updating the spin controls (which
 * redraw the overlay). Resizes are centred (the centre stays put). */
void derotation_drag_to(double dx, double dy) {
	if (!window_open || !spin_cx || !gfit) return;
	if (g_drag_mode == 1) {
		gtk_spin_button_set_value(spin_cx, dx);
		gtk_spin_button_set_value(spin_cy, (double) (gfit->ry - 1) - dy);
		return;
	}
	const double cx = gtk_spin_button_get_value(spin_cx);
	const double ycd = (double) (gfit->ry - 1) - gtk_spin_button_get_value(spin_cy);
	const double r = hypot(dx - cx, dy - ycd);
	if (r < 1.0) return;
	if (g_drag_mode == 2) gtk_spin_button_set_value(spin_radius, r);
	else if (g_drag_mode == 3) gtk_spin_button_set_value(spin_rpol, r);
}

gboolean derotation_get_disk(double *cx, double *cy, double *radius, double *rpol,
                             double *pa_deg, gboolean *mirrored) {
	if (!window_open || !spin_cx) return FALSE;
	if (cx) *cx = gtk_spin_button_get_value(spin_cx);
	if (cy) *cy = gtk_spin_button_get_value(spin_cy);
	if (radius) *radius = gtk_spin_button_get_value(spin_radius);
	if (rpol) *rpol = gtk_spin_button_get_value(spin_rpol);
	if (pa_deg) *pa_deg = gtk_spin_button_get_value(spin_pa);
	if (mirrored) *mirrored = gtk_check_button_get_active(chk_parity);
	return TRUE;
}
