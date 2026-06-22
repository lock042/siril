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
#include "core/processing.h"
#include "core/processing_thread.h"
#include "io/sequence.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"   /* redraw, REDRAW_OVERLAY */
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
static int g_drag_mode = 0;   /* 0 none, 1 move centre, 2 equatorial, 3 polar, 4 rotate */

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

/* Auto-detect the disk into the spin controls, including the pole angle from
 * the fitted major axis. Returns TRUE on success. */
static gboolean apply_autodetect(void) {
	if (!gfit) return FALSE;
	double cx, cy, r, major_deg = 0.0;
	const int body = (int) gtk_drop_down_get_selected(dd_body);
	const double flat = body_flattening(body);
	const gboolean has_rings = (body == 1);   /* Saturn */
	if (!mpp_derot_autodetect_disk(gfit, flat, has_rings, &cx, &cy, &r, &major_deg))
		return FALSE;
	gtk_spin_button_set_value(spin_cx, cx);
	gtk_spin_button_set_value(spin_cy, cy);
	gtk_spin_button_set_value(spin_radius, r);
	set_polar_default();
	/* Pole = perpendicular to the major axis (the ring plane for Saturn, the
	 * equator for Jupiter). The dialog PA points the pole up at 0 and increases
	 * clockwise; from the y-up major-axis angle that is pa = -parity·major.
	 * North is placed up-ish — flip with the on-image handle if the visible
	 * pole is the southern one. */
	const double ex = gtk_check_button_get_active(chk_parity) ? -1.0 : 1.0;
	gtk_spin_button_set_value(spin_pa, -ex * major_deg);
	return TRUE;
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

	apply_autodetect();
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
	/* Disk fitting is layered on top of the normal mouse action (the press
	 * handler intercepts clicks that hit the disc/handles via derotation_is_open
	 * + hit-test); the main view stays interactive, so mouse_status is left
	 * untouched. */
	on_derot_body_changed(NULL, NULL, NULL);   /* sync system to the body default */
	populate_from_sequence();
	redraw(REDRAW_OVERLAY);
}

static void leave_fit_mode(void) {
	window_open = FALSE;
	g_drag_mode = 0;
}

/* Hide the tool window and clear the overlay. Crucially, hand active state back
 * to the main window — hiding a focused tool window without re-presenting the
 * parent leaves GTK's active-state accounting unbalanced (it warns "Broken
 * accounting of active state" at teardown). Mirrors reactivate_main_window(). */
static void hide_dialog(void) {
	leave_fit_mode();
	if (dialog) gtk_widget_set_visible(dialog, FALSE);
	redraw(REDRAW_OVERLAY);
	GtkWidget *main_win = GTK_WIDGET(gtk_builder_get_object(gui.builder, "control_window"));
	if (main_win && gtk_widget_get_visible(main_win))
		gtk_window_present(GTK_WINDOW(main_win));
}

gboolean on_derotation_dialog_close_request(GtkWindow *win, gpointer user_data) {
	(void) win; (void) user_data;
	hide_dialog();
	return TRUE;   /* hide, don't destroy */
}

void on_derotation_close_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	hide_dialog();
}

void on_derotation_value_changed(GtkSpinButton *spin, gpointer user_data) {
	(void) spin; (void) user_data;
	if (window_open) redraw(REDRAW_OVERLAY);
}

void on_derotation_autofit_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (apply_autodetect())
		set_status(_("Disk auto-detected — adjust if needed, then compute."));
	else
		set_status(_("Could not auto-detect a disk; set the centre and radius manually."));
}

/* Flip the pole 180°: the auto-fit aligns the rotation axis from the image but
 * can't tell which end is north, so this swaps it when the southern pole is the
 * visible one. */
void on_derotation_flipnorth_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (!spin_pa) return;
	double pa = gtk_spin_button_get_value(spin_pa) + 180.0;
	while (pa >  180.0) pa -= 360.0;
	while (pa < -180.0) pa += 360.0;
	gtk_spin_button_set_value(spin_pa, pa);
}

/* The ephemeris runs once per frame (thousands of evaluations) plus the SER
 * timestamp reads, so the compute is done off the GTK thread with the busy
 * cursor; the result goes to the log and, on success, the window closes. */
struct derot_compute_args {
	planet_body_t body;
	int system, N, ry, rx;
	double cx, cy, radius, rpol, pa, parity, fps;
	gchar *epoch_text;   /* user text, may be empty -> capture midpoint */
	gchar *seqname;
	/* results, filled by the worker */
	gboolean ok;
	gchar *path;
	double span_min;
	const char *err;     /* static string when !ok */
};

static void derot_compute_args_free(struct derot_compute_args *a) {
	if (!a) return;
	g_free(a->epoch_text);
	g_free(a->seqname);
	g_free(a->path);
	free(a);
}

/* Runs on the GTK thread once the worker has finished. */
static gboolean derot_compute_idle(gpointer p) {
	struct derot_compute_args *a = p;
	stop_processing_thread();
	if (a->ok) {
		siril_log_message(_("derotation: wrote %s — %d frames, %.2f min span, "
		                    "System %d. Now register and stack the sequence.\n"),
		                  a->path, a->N, a->span_min, a->system);
		control_window_switch_to_tab(OUTPUT_LOGS);
		hide_dialog();   /* success: close the tool window */
	} else {
		siril_log_error(_("derotation: %s\n"),
		                a->err ? a->err : _("compute failed"));
		control_window_switch_to_tab(OUTPUT_LOGS);
		set_status(a->err ? a->err : _("Derotation compute failed — see the log."));
	}
	set_cursor_waiting(FALSE);
	derot_compute_args_free(a);
	return FALSE;
}

/* Runs in the processing thread. */
static gpointer derot_compute_worker(gpointer p) {
	struct derot_compute_args *a = p;
	a->ok = FALSE;

	double *jd = malloc((size_t) a->N * sizeof(double));
	if (!jd) { a->err = _("Out of memory."); siril_add_idle(derot_compute_idle, a); return GINT_TO_POINTER(1); }
	if (!mpp_derot_frame_times(&com.seq, a->fps, NAN, jd)) {
		a->err = _("No per-frame timestamps — set the frame rate (fps).");
		free(jd);
		siril_add_idle(derot_compute_idle, a);
		return GINT_TO_POINTER(1);
	}

	double epoch_jd;
	GDateTime *ep = (a->epoch_text && *a->epoch_text)
	    ? FITS_date_to_date_time(a->epoch_text) : NULL;
	if (ep) { epoch_jd = date_time_to_Julian(ep); g_date_time_unref(ep); }
	else      epoch_jd = mpp_derot_midpoint_epoch(jd, a->N);

	mpp_derot_t *d = mpp_derot_build(a->body, a->system, epoch_jd, jd, a->N,
	                                 a->ry, a->rx, a->cx, a->cy, a->radius,
	                                 a->pa, a->parity, NAN, NAN, NAN);
	free(jd);
	if (!d) {
		a->err = _("Ephemeris/plan build failed.");
		siril_add_idle(derot_compute_idle, a);
		return GINT_TO_POINTER(1);
	}
	/* Let the user's polar radius set the disk oblateness (overriding the
	 * ephemeris default) so the fitted ellipse is used as drawn. */
	if (a->rpol > 0.0 && a->rpol < a->radius)
		d->flattening = 1.0 - a->rpol / a->radius;

	a->path = g_strdup_printf("%s.derot", a->seqname);
	a->span_min = (d->jd[a->N - 1] - d->jd[0]) * 24.0 * 60.0;
	if (mpp_derot_write(a->path, d) == MPP_OK) {
		a->ok = TRUE;
	} else {
		a->err = _("Failed to write the .derot sidecar.");
	}
	mpp_derot_free(d);
	siril_add_idle(derot_compute_idle, a);
	return GINT_TO_POINTER(a->ok ? 0 : 1);
}

void on_derotation_compute_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	init_statics();
	if (!sequence_is_loaded() || com.seq.number <= 0 || !com.seq.seqname) {
		set_status(_("No sequence loaded."));
		return;
	}
	const double radius = gtk_spin_button_get_value(spin_radius);
	if (radius < 1.0) { set_status(_("Set a valid disk radius.")); return; }

	struct derot_compute_args *a = calloc(1, sizeof(*a));
	if (!a) { set_status(_("Out of memory.")); return; }
	a->body       = (planet_body_t) gtk_drop_down_get_selected(dd_body);
	a->system     = (int) gtk_drop_down_get_selected(dd_system) + 1;
	a->cx         = gtk_spin_button_get_value(spin_cx);
	a->cy         = gtk_spin_button_get_value(spin_cy);
	a->radius     = radius;
	a->rpol       = gtk_spin_button_get_value(spin_rpol);
	a->pa         = gtk_spin_button_get_value(spin_pa);
	a->parity     = gtk_check_button_get_active(chk_parity) ? -1.0 : 1.0;
	a->fps        = gtk_spin_button_get_value(spin_fps);
	a->epoch_text = g_strdup(gtk_editable_get_text(GTK_EDITABLE(entry_epoch)));
	a->seqname    = g_strdup(com.seq.seqname);
	a->N          = com.seq.number;
	a->ry         = (int) com.seq.ry;
	a->rx         = (int) com.seq.rx;

	set_status(_("Computing derotation…"));
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(derot_compute_worker, a)) {
		set_cursor_waiting(FALSE);
		derot_compute_args_free(a);
		set_status(_("Could not start the derotation worker (another job is running?)."));
	}
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
gboolean derotation_rotate_handle(double *hx, double *hy, gboolean *inside) {
	double cx, ycd, req, rpol, eux, euy, pux, puy;
	if (!fit_display_geom(&cx, &ycd, &req, &rpol, &eux, &euy, &pux, &puy))
		return FALSE;
	double rh = rpol + fmax(16.0, rpol * 0.30);   /* normal: beyond the pole tip */
	double x = cx + rh * pux, y = ycd + rh * puy;
	gboolean in = FALSE;
	if (gfit && (x < 0.0 || x >= (double) gfit->rx
	          || y < 0.0 || y >= (double) gfit->ry)) {
		/* Off the image (tightly-cropped planet): pull the handle to the
		 * midpoint of the north axis, inside the disc, so it stays grabbable. */
		rh = 0.5 * rpol;
		x = cx + rh * pux; y = ycd + rh * puy;
		in = TRUE;
	}
	if (hx) *hx = x;
	if (hy) *hy = y;
	if (inside) *inside = in;
	return TRUE;
}

int derotation_hit_test(double dx, double dy) {
	double cx, ycd, req, rpol, eux, euy, pux, puy;
	if (!fit_display_geom(&cx, &ycd, &req, &rpol, &eux, &euy, &pux, &puy))
		return 0;
	const double tol = fmax(8.0, req * 0.12);
	/* Rotation handle (drag to spin the PA); position from the shared helper so
	 * it tracks the overlay, including the off-image fallback. */
	double rhx, rhy;
	if (derotation_rotate_handle(&rhx, &rhy, NULL)
	 && hypot(dx - rhx, dy - rhy) <= tol) return 4;
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
	if (g_drag_mode == 4) {
		const double vx = dx - cx, vy = dy - ycd;
		if (hypot(vx, vy) < 1.0) return;
		/* Pole unit (cairo) is (ex·sin pa, -cos pa); solve for the PA that
		 * points it at the cursor. ex is the image parity. */
		const double ex = gtk_check_button_get_active(chk_parity) ? -1.0 : 1.0;
		const double pa = atan2(vx * ex, -vy) * 180.0 / M_PI;
		gtk_spin_button_set_value(spin_pa, pa);
		return;
	}
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
