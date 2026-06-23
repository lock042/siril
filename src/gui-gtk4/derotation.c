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
#include "io/films.h"                /* film_struct->fps for AVI frame rate */
#include "io/single_image.h"        /* open_single_image */
#include "io/image_format_fits.h"   /* savefits, clearfits, set_right_extension */
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"   /* redraw, REDRAW_OVERLAY */
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/derotation.h"

#include "algos/planet_ephem.h"
#include "registration/mpp.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_derot_sidecar.h"
#include "registration/mpp/mpp_derot_build.h"
#include "registration/mpp/mpp_session.h"

static GtkWidget *dialog = NULL;
static GtkDropDown *dd_body = NULL, *dd_system = NULL;
static GtkEntry *entry_epoch = NULL;
static GtkSpinButton *spin_cx = NULL, *spin_cy = NULL, *spin_radius = NULL,
                     *spin_rpol = NULL, *spin_pa = NULL, *spin_fps = NULL;
static GtkCheckButton *chk_parity = NULL;
static GtkLabel *status = NULL;
static gboolean window_open = FALSE;
static int g_drag_mode = 0;   /* 0 none, 1 move centre, 2 equatorial, 3 polar, 4 rotate */

/* Multi-sequence session (Option B): the set of sequences to combine, each
 * tagged with a colour channel and its own disk fit, plus the designated
 * reference. Built up as the user fits one sequence at a time. */
static mpp_session_t *session = NULL;
static GtkDropDown *dd_channel = NULL;   /* channel tag for the next "Add" */
static GtkListBox *seqlist = NULL;       /* one row per session entry */
static GtkWidget *combine_frame = NULL;  /* multi-sequence panel (full mode only) */
/* The same dialog serves two entry points: the Registration-tab button opens it
 * in simple single-sequence mode (just fit + Compute & save, closes on save);
 * the Tools menu opens it in full multi-sequence mode (the combine panel and the
 * MPP parameter widgets, and Compute & save stays open). */
static gboolean multi_mode = FALSE;
/* Base name of the sequence whose fit currently sits in the spin controls, so
 * "Add" can refuse a sequence that hasn't been fitted (the spins would still
 * hold the previous sequence's disk). NULL until a fit is made. */
static gchar *fitted_for = NULL;

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

static gboolean apply_autodetect(void);
static void set_status(const char *msg);
static void on_derotation_seqlist_row_activated(GtkListBox *box, GtkListBoxRow *row, gpointer u);
static void mark_fit_valid(void);
static void update_session_epoch(void);
static void set_model_locked(gboolean locked);

/* Default rotation system per body (matches the command): Jupiter -> II,
 * Saturn -> III, Mars -> I. Also refresh the polar-radius default. */
static void on_derot_body_changed(GObject *obj, GParamSpec *pspec, gpointer u) {
	(void) pspec; (void) u;
	if (!dd_system) return;
	const guint body = gtk_drop_down_get_selected(dd_body);
	const guint sys = (body == 0) ? 1u : (body == 1) ? 2u : 0u;   /* Jup II, Sat III, Mars I */
	gtk_drop_down_set_selected(dd_system, sys);
	set_polar_default();
	/* A genuine change of planet (the notify signal, not the manual sync during
	 * window show) re-runs the auto-fit so the centre, size and orientation suit
	 * the new body — Saturn uses the ring sizing, Jupiter the minor axis, etc. */
	if (obj && window_open) apply_autodetect();
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
	dd_channel  = GTK_DROP_DOWN   (gtk_builder_get_object(gui.builder, "derot_channel"));
	seqlist     = GTK_LIST_BOX    (gtk_builder_get_object(gui.builder, "derot_seqlist"));
	combine_frame = GTK_WIDGET    (gtk_builder_get_object(gui.builder, "derot_combine_frame"));
	/* Default the MPP param widgets to mpp_config_defaults (the .ui sets the
	 * spin defaults; set the toggles/combo here). */
	{
		GtkBuilder *b = gui.builder;
		#define MDP_SET(id, on) do { GObject *o = gtk_builder_get_object(b, id); \
			if (o) gtk_check_button_set_active(GTK_CHECK_BUTTON(o), on); } while (0)
		MDP_SET("mdp_normalize", TRUE);
		MDP_SET("mdp_dewarp", TRUE);
		MDP_SET("mdp_seed", TRUE);
		MDP_SET("mdp_skip_failed", FALSE);
		#undef MDP_SET
	}
	/* One-shot colour is the common case, so make it the default tag. */
	if (dd_channel)
		gtk_drop_down_set_selected(dd_channel, MPP_CHAN_OSC);
	/* Double-click a row to switch to that sequence and restore its fit. */
	if (seqlist)
		g_signal_connect(seqlist, "row-activated",
		                 G_CALLBACK(on_derotation_seqlist_row_activated), NULL);
	/* Restore the remembered mirror state when the window is first created. */
	if (chk_parity)
		gtk_check_button_set_active(chk_parity, com.pref.gui.derot_mirror);
	if (dd_body)
		g_signal_connect(dd_body, "notify::selected",
		                 G_CALLBACK(on_derot_body_changed), NULL);
}

/* ---- multi-sequence session ---- */

static const char *channel_label(mpp_channel_t ch) {
	switch (ch) {
		case MPP_CHAN_R:   return "R";
		case MPP_CHAN_G:   return "G";
		case MPP_CHAN_B:   return "B";
		case MPP_CHAN_LUM: return "L";
		case MPP_CHAN_OSC: return "OSC";
		default:           return "mono";
	}
}

/* short filename suffix for a channel's stacked output */
static const char *channel_suffix(mpp_channel_t ch) {
	switch (ch) {
		case MPP_CHAN_R:   return "R";
		case MPP_CHAN_G:   return "G";
		case MPP_CHAN_B:   return "B";
		case MPP_CHAN_LUM: return "L";
		case MPP_CHAN_OSC: return "OSC";
		default:           return "mono";
	}
}

/* Rebuild the list box from the session. Each row shows the sequence name, its
 * channel tag, "(fit)" once a disk has been recorded, and a star on the
 * reference. The row index matches the session index. */
static void refresh_seqlist(void) {
	if (!seqlist) return;
	GtkWidget *child;
	while ((child = gtk_widget_get_first_child(GTK_WIDGET(seqlist))))
		gtk_list_box_remove(seqlist, child);
	if (!session) return;
	for (int i = 0; i < session->count; i++) {
		const mpp_session_seq_t *e = &session->seqs[i];
		gchar *base = g_path_get_basename(e->seqname);
		gchar *txt = g_strdup_printf("%s%s  [%s]%s",
		                             i == session->reference ? "★ " : "",
		                             base, channel_label(e->channel),
		                             e->has_fit ? "" : _("  — not fitted"));
		GtkWidget *row = gtk_label_new(txt);
		gtk_label_set_xalign(GTK_LABEL(row), 0.0);
		gtk_label_set_ellipsize(GTK_LABEL(row), PANGO_ELLIPSIZE_END);
		gtk_widget_set_tooltip_text(row, e->seqname);   /* full path on hover */
		gtk_widget_set_margin_start(row, 4);
		gtk_widget_set_margin_end(row, 4);
		gtk_list_box_append(seqlist, row);
		g_free(txt);
		g_free(base);
	}
}

/* Read the loaded sequence's capture span (UTC JD). Falls back to the dialog's
 * fps when there are no per-frame timestamps. Returns FALSE if neither works. */
static gboolean current_seq_span(double *first, double *last) {
	const double fps = spin_fps ? gtk_spin_button_get_value(spin_fps) : 0.0;
	return mpp_derot_sequence_span(&com.seq, fps, NAN, first, last);
}

/* A CFA-pattern signature for compatibility checks: the SER ColorID, or -1 for
 * mono / non-SER inputs. */
static int current_seq_bayer(void) {
	if (com.seq.type == SEQ_SER && com.seq.ser_file)
		return (int) com.seq.ser_file->color_id;
	return -1;
}

/* Record that the disk now in the spin controls was fitted for the loaded
 * sequence (gates "Add"; see fitted_for). */
static void mark_fit_valid(void) {
	g_free(fitted_for);
	fitted_for = (sequence_is_loaded() && com.seq.seqname)
	    ? g_strdup(com.seq.seqname) : NULL;
}

/* Lock the planet model (body + rotation system) once a session exists: all
 * sequences in one combine must be the same body, fixed by the reference. */
static void set_model_locked(gboolean locked) {
	if (dd_body)   gtk_widget_set_sensitive(GTK_WIDGET(dd_body), !locked);
	if (dd_system) gtk_widget_set_sensitive(GTK_WIDGET(dd_system), !locked);
}

/* Show the shared reference epoch (union-span midpoint over all session
 * sequences) in the epoch entry, so it is visible and editable. */
static void update_session_epoch(void) {
	if (!entry_epoch || !session || session->count <= 0) return;
	double *firsts = g_new(double, session->count);
	double *lasts  = g_new(double, session->count);
	for (int i = 0; i < session->count; i++) {
		firsts[i] = session->seqs[i].first_jd;
		lasts[i]  = session->seqs[i].last_jd;
	}
	const double epoch_jd = mpp_derot_union_epoch(firsts, lasts, session->count);
	g_free(firsts); g_free(lasts);

	GDateTime *j2000 = g_date_time_new_utc(2000, 1, 1, 12, 0, 0);
	GDateTime *mid = j2000
	    ? g_date_time_add_seconds(j2000, (epoch_jd - 2451545.0) * 86400.0) : NULL;
	if (j2000) g_date_time_unref(j2000);
	if (mid) {
		gchar *s = date_time_to_FITS_date(mid);
		if (s) { gtk_editable_set_text(GTK_EDITABLE(entry_epoch), s); g_free(s); }
		g_date_time_unref(mid);
	}
}

/* "Add current sequence": snapshot the live disk fit + channel into the session
 * (creating it on first use). A sequence already present is updated in place. */
void on_derotation_add_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	init_statics();
	if (!sequence_is_loaded() || com.seq.number <= 0 || !com.seq.seqname) {
		set_status(_("Load a sequence before adding it to the session."));
		return;
	}
	if (gtk_spin_button_get_value(spin_radius) < 1.0) {
		set_status(_("Fit the disk (radius) before adding the sequence."));
		return;
	}
	/* The spin controls must hold a fit made for THIS sequence — otherwise the
	 * previous sequence's disk would be recorded against it. */
	if (g_strcmp0(fitted_for, com.seq.seqname) != 0) {
		set_status(_("Fit this sequence's disk first (Auto-detect or adjust it on "
		             "the image)."));
		return;
	}

	const int bayer = current_seq_bayer();

	/* Basic compatibility: every sequence must share the reference's frame size
	 * and CFA layout, or they cannot be stacked into one canvas. */
	if (session && session->count > 0) {
		const mpp_session_seq_t *ref = &session->seqs[
		    (session->reference >= 0) ? session->reference : 0];
		if ((int) com.seq.rx != ref->frame_cols || (int) com.seq.ry != ref->frame_rows) {
			set_status(_("Incompatible dimensions — every sequence must match the "
			             "reference frame size."));
			return;
		}
		if (bayer != ref->bayer) {
			set_status(_("Incompatible Bayer pattern — it must match the reference "
			             "sequence."));
			return;
		}
	}

	if (!session) {
		const int body = (int) gtk_drop_down_get_selected(dd_body);
		const int sys  = (int) gtk_drop_down_get_selected(dd_system) + 1;
		session = mpp_session_new((planet_body_t) body, sys);
		if (!session) { set_status(_("Out of memory.")); return; }
	}

	double first = 0, last = 0;
	if (!current_seq_span(&first, &last)) {
		set_status(_("No per-frame timestamps — set the frame rate (fps) first."));
		return;
	}
	const mpp_channel_t ch = dd_channel
	    ? (mpp_channel_t) gtk_drop_down_get_selected(dd_channel) : MPP_CHAN_MONO;

	int idx = -1;
	for (int i = 0; i < session->count; i++)
		if (!g_strcmp0(session->seqs[i].seqname, com.seq.seqname)) { idx = i; break; }
	if (idx < 0) {
		idx = mpp_session_add(session, com.seq.seqname, ch, first, last, com.seq.number);
		if (idx < 0) { set_status(_("Could not add the sequence.")); return; }
	} else {
		mpp_session_set_channel(session, idx, ch);
	}
	mpp_session_set_fit(session, idx,
	                    gtk_spin_button_get_value(spin_cx),
	                    gtk_spin_button_get_value(spin_cy),
	                    gtk_spin_button_get_value(spin_radius),
	                    gtk_spin_button_get_value(spin_rpol),
	                    gtk_spin_button_get_value(spin_pa),
	                    gtk_check_button_get_active(chk_parity) ? -1.0 : 1.0,
	                    (int) com.seq.ry, (int) com.seq.rx, bayer,
	                    spin_fps ? gtk_spin_button_get_value(spin_fps) : 0.0);
	set_model_locked(TRUE);     /* the body is now fixed by the session */
	update_session_epoch();     /* show the shared epoch */
	refresh_seqlist();
	set_status(_("Sequence added. Load and fit the next, or combine."));
}

/* index of the currently selected session row, or -1 */
static int selected_session_index(void) {
	if (!seqlist) return -1;
	GtkListBoxRow *row = gtk_list_box_get_selected_row(seqlist);
	return row ? gtk_list_box_row_get_index(row) : -1;
}

void on_derotation_setref_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	const int idx = selected_session_index();
	if (!session || idx < 0) { set_status(_("Select a sequence in the list first.")); return; }
	mpp_session_set_reference(session, idx);
	refresh_seqlist();
	set_status(_("Reference sequence set — it defines the common output frame."));
}

void on_derotation_remove_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	const int idx = selected_session_index();
	if (!session || idx < 0) { set_status(_("Select a sequence in the list first.")); return; }
	mpp_session_remove(session, idx);
	if (session->count == 0)
		set_model_locked(FALSE);   /* empty session: planet selectable again */
	else
		update_session_epoch();
	refresh_seqlist();
}

/* Double-click a list row: load that sequence and restore its recorded fit into
 * the controls, so the user can review or re-fit it. */
static void on_derotation_seqlist_row_activated(GtkListBox *box, GtkListBoxRow *row,
                                                gpointer u) {
	(void) box; (void) u;
	if (!session || !row) return;
	const int idx = gtk_list_box_row_get_index(row);
	if (idx < 0 || idx >= session->count) return;
	const mpp_session_seq_t *e = &session->seqs[idx];

	/* switch the displayed sequence (no-op if it is already current) */
	if (g_strcmp0(com.seq.seqname, e->seqname) != 0)
		set_seq(e->seqname);

	/* restore the stored disk fit + channel into the controls */
	if (e->has_fit) {
		gtk_spin_button_set_value(spin_cx, e->cx);
		gtk_spin_button_set_value(spin_cy, e->cy);
		gtk_spin_button_set_value(spin_radius, e->r_eq);
		gtk_spin_button_set_value(spin_rpol, e->r_pol);
		gtk_spin_button_set_value(spin_pa, e->pa_deg);
		if (chk_parity) gtk_check_button_set_active(chk_parity, e->parity < 0.0);
		if (spin_fps && e->fps > 0.0) gtk_spin_button_set_value(spin_fps, e->fps);
	}
	if (dd_channel) gtk_drop_down_set_selected(dd_channel, e->channel);
	/* the fit now matches the loaded sequence */
	g_free(fitted_for);
	fitted_for = g_strdup(e->seqname);
	if (window_open) redraw(REDRAW_OVERLAY);
	set_status(_("Switched to the selected sequence; adjust and re-add if needed."));
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
	mark_fit_valid();   /* this fit belongs to the loaded sequence */
	return TRUE;
}

/* Fill the epoch entry with the capture midpoint, and auto-detect the disk.
 * Once a multi-sequence session exists the epoch is the shared union midpoint
 * (managed by update_session_epoch), so this single-sequence midpoint is only
 * applied when no session is in progress. */
static void populate_from_sequence(void) {
	if (!sequence_is_loaded() || com.seq.number <= 0) {
		set_status(_("Load a planetary sequence first."));
		return;
	}
	/* AVI has no per-frame timestamps; auto-fill the frame rate from the
	 * container so the epoch can be computed. SER carries timestamps, so its
	 * fps fallback is left at 0. */
#ifdef HAVE_FFMS2
	if (com.seq.type == SEQ_AVI && com.seq.film_file
	    && com.seq.film_file->fps > 0.0 && spin_fps)
		gtk_spin_button_set_value(spin_fps, com.seq.film_file->fps);
#endif
	const double fps = spin_fps ? gtk_spin_button_get_value(spin_fps) : 0.0;
	const gboolean has_session = session && session->count > 0;
	const int N = com.seq.number;
	double *jd = malloc((size_t) N * sizeof(double));
	if (jd && mpp_derot_frame_times(&com.seq, fps, NAN, jd)) {
		if (!has_session) {
			/* Build the midpoint datetime from the J2000 anchor (full JD
			 * 2451545.0): Julian_to_date_time() expects MJD, not the full JD
			 * that date_time_to_Julian()/the .derot use, so using it here would
			 * land the default epoch thousands of years off. */
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
		}
		set_status(has_session
		    ? _("Using the shared epoch — fit the disk, then add or combine.")
		    : _("Fit the disk to the planet, then compute."));
	} else {
		set_status(_("No per-frame timestamps — set the frame rate (fps)."));
	}
	free(jd);

	apply_autodetect();
}

/* A sequence was loaded while the tool is open: refresh the epoch + frame rate
 * for it and re-fit the disk. (Selecting a row during a future "pause" will be
 * handled separately; this is the ordinary load path.) */
void derotation_sequence_changed(void) {
	if (!window_open) return;
	init_statics();
	populate_from_sequence();
	redraw(REDRAW_OVERLAY);
}

/* Show/hide the MPP register/stack parameter column (full mode only).
 * Populated in Phase 2; a no-op until the param widgets exist. */
static void apply_params_visibility(gboolean show) {
	GtkWidget *col = GTK_WIDGET(gtk_builder_get_object(gui.builder, "derot_params_col"));
	if (col) gtk_widget_set_visible(col, show);
}

/* Show/hide the full-mode sections and set the title for the current mode. */
static void apply_dialog_mode(void) {
	if (combine_frame)
		gtk_widget_set_visible(combine_frame, multi_mode);
	apply_params_visibility(multi_mode);
	if (dialog)
		gtk_window_set_title(GTK_WINDOW(dialog),
		    multi_mode ? _("Multi-sequence planetary derotation")
		               : _("Planetary derotation"));
}

static void present_dialog(gboolean multi) {
	init_statics();
	if (!dialog) return;
	multi_mode = multi;
	apply_dialog_mode();
	GtkWindow *parent = GTK_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	if (parent)
		gtk_window_set_transient_for(GTK_WINDOW(dialog), parent);
	gtk_window_present(GTK_WINDOW(dialog));
}

/* Registration-tab "Derotation" button: simple single-sequence mode. */
void on_seqmpp_derotation_button_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	if (!sequence_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No sequence"),
		    _("Load the planetary sequence (video) before setting up derotation."));
		return;
	}
	present_dialog(FALSE);
}

/* Tools menu "Multi-sequence derotation...": full multi-sequence mode. Opens
 * with or without a sequence loaded — the user loads and adds the first (and
 * each subsequent) sequence from the main window while the tool stays open. */
void on_derotation_multi_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	present_dialog(TRUE);
}

void on_derotation_dialog_show(GtkWidget *widget, gpointer user_data) {
	(void) widget; (void) user_data;
	init_statics();
	window_open = TRUE;
	apply_dialog_mode();
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
	/* Any change to the disk controls counts as fitting the loaded sequence. */
	if (window_open) { mark_fit_valid(); redraw(REDRAW_OVERLAY); }
}

/* Mirror toggle: remember the state in preferences (persisted to the initfile)
 * and flip the overlay's east-west sense. */
void on_derotation_parity_toggled(GtkCheckButton *button, gpointer user_data) {
	(void) user_data;
	com.pref.gui.derot_mirror = gtk_check_button_get_active(button);
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
		/* Simple single-sequence mode closes on save (the old behaviour); the
		 * multi-sequence tool stays open so the user can keep adding sequences. */
		if (!multi_mode)
			hide_dialog();
		else
			set_status(_("Saved the .derot — add this sequence, or fit the next."));
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

/* ---- multi-sequence combine ---- */

struct derot_combine_args {
	int n;
	sequence **seqs;
	mpp_derot_t **derots;
	gboolean *is_com;
	mpp_channel_t *channels;
	int ref_index;
	gchar *base;
	mpp_config_t cfg;
	/* results */
	int saved, failed;
	gboolean composed;       /* an RGB image was assembled from R+G+B */
	gchar *rgb_name;
	gchar *load_name;        /* result to open on completion (base name, no ext) */
};

static void derot_combine_args_free(struct derot_combine_args *a) {
	if (!a) return;
	for (int i = 0; i < a->n; i++) {
		if (a->derots && a->derots[i]) mpp_derot_free(a->derots[i]);
		if (a->seqs && a->seqs[i] && !a->is_com[i]) free_sequence(a->seqs[i], TRUE);
	}
	g_free(a->seqs);
	g_free(a->derots);
	g_free(a->is_com);
	g_free(a->channels);
	g_free(a->base);
	g_free(a->rgb_name);
	g_free(a->load_name);
	free(a);
}

static gboolean derot_combine_idle(gpointer p) {
	struct derot_combine_args *a = p;
	stop_processing_thread();
	if (a->saved > 0 && a->composed)
		siril_log_message(_("Derotation combine: saved %d channel stack(s) and "
		                    "assembled %s%s.\n"),
		                  a->saved, a->rgb_name,
		                  a->failed ? _(" (some channels failed)") : "");
	else if (a->saved > 0)
		siril_log_message(_("Derotation combine: saved %d channel stack(s)%s. "
		                    "Assemble colour with the compositing tools.\n"),
		                  a->saved, a->failed ? _(" (some channels failed)") : "");
	else
		siril_log_error(_("Derotation combine: no channels stacked — see the log.\n"));
	control_window_switch_to_tab(OUTPUT_LOGS);
	set_status(a->saved > 0
	    ? (a->composed ? _("Combine done — colour image assembled (see the log).")
	                   : _("Combine done — see the log for the saved files."))
	    : _("Combine failed — see the log."));
	set_cursor_waiting(FALSE);
	/* Open the primary result so it is shown straight away. */
	if (a->saved > 0 && a->load_name) {
		gchar *path = set_right_extension(a->load_name);
		if (path) { open_single_image(path); g_free(path); }
	}
	derot_combine_args_free(a);
	return FALSE;
}

/* Assemble three co-registered mono channel stacks into one 3-layer image (Siril
 * planar layout: R plane, then G, then B). The channels already share the
 * reference canvas, so this is a straight per-plane copy — no alignment. Handles
 * both 16-bit and 32-bit float stacks. */
static gboolean compose_rgb(const fits *R, const fits *G, const fits *B, fits *out) {
	const int W = R->rx, H = R->ry;
	if (G->rx != W || G->ry != H || B->rx != W || B->ry != H) return FALSE;
	if (R->type != G->type || R->type != B->type) return FALSE;
	const size_t plane = (size_t) W * H;

	clearfits(out);
	out->rx = W; out->ry = H;
	out->naxes[0] = W; out->naxes[1] = H; out->naxes[2] = 3;
	out->naxis = 3;
	if (R->type == DATA_FLOAT) {
		out->fdata = malloc(plane * 3 * sizeof(float));
		if (!out->fdata) return FALSE;
		memcpy(out->fdata,             R->fdata, plane * sizeof(float));
		memcpy(out->fdata + plane,     G->fdata, plane * sizeof(float));
		memcpy(out->fdata + plane * 2, B->fdata, plane * sizeof(float));
		out->bitpix = FLOAT_IMG; out->orig_bitpix = FLOAT_IMG; out->type = DATA_FLOAT;
		out->fpdata[0] = out->fdata;
		out->fpdata[1] = out->fdata + plane;
		out->fpdata[2] = out->fdata + plane * 2;
	} else {
		out->data = malloc(plane * 3 * sizeof(WORD));
		if (!out->data) return FALSE;
		memcpy(out->data,             R->data, plane * sizeof(WORD));
		memcpy(out->data + plane,     G->data, plane * sizeof(WORD));
		memcpy(out->data + plane * 2, B->data, plane * sizeof(WORD));
		out->bitpix = USHORT_IMG; out->orig_bitpix = USHORT_IMG; out->type = DATA_USHORT;
		out->pdata[0] = out->data;
		out->pdata[1] = out->data + plane;
		out->pdata[2] = out->data + plane * 2;
	}
	return TRUE;
}

static gpointer derot_combine_worker(gpointer p) {
	struct derot_combine_args *a = p;
	a->saved = a->failed = 0;

	/* one stack per distinct channel, each landing in the reference canvas */
	static const mpp_channel_t order[] = {
		MPP_CHAN_R, MPP_CHAN_G, MPP_CHAN_B, MPP_CHAN_LUM,
		MPP_CHAN_OSC, MPP_CHAN_MONO
	};
	/* keep each channel's result so co-registered mono R/G/B can be composed */
	fits chan_fits[MPP_CHAN_COUNT];
	gboolean have[MPP_CHAN_COUNT];
	memset(chan_fits, 0, sizeof(chan_fits));
	memset(have, 0, sizeof(have));

	for (size_t k = 0; k < G_N_ELEMENTS(order); k++) {
		const mpp_channel_t ch = order[k];
		sequence **cs = g_new0(sequence *, a->n);
		mpp_derot_t **cd = g_new0(mpp_derot_t *, a->n);
		int m = 0;
		for (int i = 0; i < a->n; i++)
			if (a->channels[i] == ch) { cs[m] = a->seqs[i]; cd[m] = a->derots[i]; m++; }
		if (m == 0) { g_free(cs); g_free(cd); continue; }

		siril_log_message(_("Derotation combine: channel %s (%d sequence(s))…\n"),
		                  channel_label(ch), m);

		fits out = { 0 };
		const mpp_status_t st = mpp_multistack_to(
		    cs, cd, m, a->seqs[a->ref_index], a->derots[a->ref_index], &a->cfg, &out);
		if (st == MPP_OK) {
			gchar *name = (ch == MPP_CHAN_MONO)
			    ? g_strdup_printf("%s_derot_stack", a->base)
			    : g_strdup_printf("%s_%s", a->base, channel_suffix(ch));
			if (savefits(name, &out) == 0) {
				siril_log_message(_("Derotation combine: saved %s (%d sequence(s))\n"),
				                  name, m);
				a->saved++;
			} else {
				a->failed++;
			}
			g_free(name);
			chan_fits[ch] = out;   /* keep for compose (ownership transferred) */
			have[ch] = TRUE;
		} else {
			siril_log_error(_("Derotation combine: channel %s failed (code %d)\n"),
			                channel_label(ch), st);
			a->failed++;
			clearfits(&out);
		}
		g_free(cs);
		g_free(cd);
		if (!processing_should_continue()) break;   /* cancelled */
	}

	/* RGB compose: the three mono channels share the reference canvas, so a
	 * straight planar assembly yields a colour image (LRGB is left to the
	 * compositing tools when a separate Luminance was also stacked). */
	if (have[MPP_CHAN_R] && have[MPP_CHAN_G] && have[MPP_CHAN_B]) {
		fits rgb = { 0 };
		if (compose_rgb(&chan_fits[MPP_CHAN_R], &chan_fits[MPP_CHAN_G],
		                &chan_fits[MPP_CHAN_B], &rgb)) {
			a->rgb_name = g_strdup_printf("%s_rgb", a->base);
			if (savefits(a->rgb_name, &rgb) == 0) {
				siril_log_message(_("Derotation combine: assembled colour image %s\n"),
				                  a->rgb_name);
				a->composed = TRUE;
			} else {
				siril_log_error(_("Derotation combine: failed to save the colour image\n"));
			}
		}
		clearfits(&rgb);
	}

	/* Pick the result to open on completion: the assembled colour image if R+G+B
	 * were combined, otherwise the stack of the channel the reference sequence
	 * belongs to (which is the only channel in the single-channel case). */
	if (a->composed) {
		a->load_name = g_strdup(a->rgb_name);
	} else {
		const mpp_channel_t rch = a->channels[a->ref_index];
		if (have[rch])
			a->load_name = (rch == MPP_CHAN_MONO)
			    ? g_strdup_printf("%s_derot_stack", a->base)
			    : g_strdup_printf("%s_%s", a->base, channel_suffix(rch));
	}

	for (int ch = 0; ch < MPP_CHAN_COUNT; ch++)
		if (have[ch]) clearfits(&chan_fits[ch]);

	set_progress_bar_data(a->saved > 0 ? _("Derotation combine complete")
	                                   : _("Derotation combine failed"),
	                      a->saved > 0 ? PROGRESS_DONE : PROGRESS_RESET);
	siril_add_idle(derot_combine_idle, a);
	return GINT_TO_POINTER(a->saved > 0 ? 0 : 1);
}

/* Copy the MPP registration + stacking settings from their page widgets into
 * cfg, so the combine behaves exactly like the single-sequence register + stack
 * the user already tunes there — most importantly the per-AP frame selection
 * (stack/register percent), without which every frame is averaged and the
 * planet comes out soft. Reads gui.builder directly so it does not depend on
 * those pages having been opened. Fields whose widget is missing keep cfg's
 * existing (default) value. */
static void apply_gui_mpp_config(mpp_config_t *cfg) {
	GtkBuilder *b = gui.builder;
	if (!cfg || !b) return;
	#define MPP_SPIN_I(f, id) do { GObject *o = gtk_builder_get_object(b, id); \
		if (o) cfg->f = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(o)); } while (0)
	#define MPP_SPIN_D(f, id) do { GObject *o = gtk_builder_get_object(b, id); \
		if (o) cfg->f = gtk_spin_button_get_value(GTK_SPIN_BUTTON(o)); } while (0)
	#define MPP_TOGG(f, id) do { GObject *o = gtk_builder_get_object(b, id); \
		if (o) cfg->f = siril_toggle_get_active(GTK_WIDGET(o)); } while (0)
	/* registration parameters (the dialog's own widgets) */
	MPP_SPIN_I(alignment_points_half_box_width,       "mdp_half_box");
	MPP_SPIN_I(alignment_points_search_width,         "mdp_search_width");
	MPP_SPIN_I(align_frames_search_width,             "mdp_search_global");
	MPP_SPIN_I(align_frames_average_frame_percent,    "mdp_ref_percent");
	MPP_SPIN_I(alignment_points_brightness_threshold, "mdp_min_brightness");
	MPP_SPIN_I(alignment_points_contrast_threshold,   "mdp_min_contrast");
	MPP_SPIN_D(alignment_points_structure_threshold,  "mdp_min_structure");
	MPP_SPIN_I(alignment_points_frame_percent,        "mdp_reg_percent");
	MPP_TOGG(alignment_points_de_warp,                "mdp_dewarp");
	MPP_TOGG(frames_normalization,                    "mdp_normalize");
	MPP_TOGG(align_frames_seed_from_regdata,          "mdp_seed");
	/* Derotation only makes sense for a full-disc planet, so always use the
	 * Planet (centroid) global-alignment mode and treat it as a fast-changing
	 * object (per-AP frame selection from a sliding window). */
	cfg->align_frames_mode = MPP_ALIGN_PLANET;
	cfg->align_frames_fast_changing_object = TRUE;
	/* stacking parameters */
	MPP_SPIN_I(stack_frame_percent,                    "mdp_stack_percent");
	MPP_SPIN_I(stack_frame_number,                     "mdp_stack_frames");
	MPP_SPIN_D(stack_frames_background_fraction,        "mdp_bg_fraction");
	MPP_SPIN_D(stack_frames_background_blend_threshold, "mdp_bg_blend");
	MPP_TOGG(stack_skip_failed_aps,                     "mdp_skip_failed");
	{
		static const double scales[] = { 1.0, 1.5, 2.0, 3.0 };
		GObject *o = gtk_builder_get_object(b, "mdp_scale");
		if (o) {
			const guint i = gtk_drop_down_get_selected(GTK_DROP_DOWN(o));
			cfg->drizzle_scale = (i < G_N_ELEMENTS(scales)) ? scales[i] : 1.0;
		}
	}
	#undef MPP_SPIN_I
	#undef MPP_SPIN_D
	#undef MPP_TOGG
}

void on_derotation_combine_clicked(GtkButton *button, gpointer user_data) {
	(void) button; (void) user_data;
	init_statics();
	if (!session || session->count < 1) {
		set_status(_("Add at least one fitted sequence to the session first."));
		return;
	}
	for (int i = 0; i < session->count; i++)
		if (!session->seqs[i].has_fit) {
			set_status(_("Every sequence in the session must be fitted first."));
			return;
		}
	if (session->reference < 0 || session->reference >= session->count) {
		set_status(_("Select a reference sequence first."));
		return;
	}

	/* Shared reference epoch: the user's text if given, else the union midpoint
	 * of every sequence's span. */
	double epoch_jd;
	const char *etxt = entry_epoch ? gtk_editable_get_text(GTK_EDITABLE(entry_epoch)) : NULL;
	GDateTime *ep = (etxt && *etxt) ? FITS_date_to_date_time((gchar *) etxt) : NULL;
	if (ep) { epoch_jd = date_time_to_Julian(ep); g_date_time_unref(ep); }
	else {
		double *firsts = g_new(double, session->count);
		double *lasts  = g_new(double, session->count);
		for (int i = 0; i < session->count; i++) {
			firsts[i] = session->seqs[i].first_jd;
			lasts[i]  = session->seqs[i].last_jd;
		}
		epoch_jd = mpp_derot_union_epoch(firsts, lasts, session->count);
		g_free(firsts); g_free(lasts);
	}

	/* Load every sequence and build its derotation plan at the shared epoch
	 * (cheap; the heavy stacking runs off-thread below). */
	const int n = session->count;
	struct derot_combine_args *a = calloc(1, sizeof(*a));
	if (!a) { set_status(_("Out of memory.")); return; }
	a->n = n;
	a->seqs     = g_new0(sequence *, n);
	a->derots   = g_new0(mpp_derot_t *, n);
	a->is_com   = g_new0(gboolean, n);
	a->channels = g_new0(mpp_channel_t, n);
	a->ref_index = session->reference;
	a->base = g_strdup(session->seqs[session->reference].seqname);
	mpp_config_defaults(&a->cfg);
	apply_gui_mpp_config(&a->cfg);   /* use the user's register/stack settings */
	a->cfg.stack_method = MPP_STACK_WARP;
	/* The combine registers and stacks in one streamed pass, ranking every
	 * frame, so the per-AP selection IS the stack selection. Drive it from the
	 * Stacking-tab percent/number (the quality pass reads
	 * alignment_points_frame_*); otherwise the register percent (often 100)
	 * would stack every frame and average the detail away. */
	a->cfg.alignment_points_frame_number  = a->cfg.stack_frame_number;
	a->cfg.alignment_points_frame_percent = a->cfg.stack_frame_percent;

	gboolean ok = TRUE;
	for (int i = 0; i < n && ok; i++) {
		const mpp_session_seq_t *e = &session->seqs[i];
		a->channels[i] = e->channel;
		a->seqs[i] = load_sequence_force_debayer(e->seqname);
		if (!a->seqs[i]) {
			siril_log_error(_("Derotation combine: cannot load '%s'\n"), e->seqname);
			ok = FALSE; break;
		}
		a->is_com[i] = check_seq_is_comseq(a->seqs[i]);
		if (a->is_com[i]) { free_sequence(a->seqs[i], TRUE); a->seqs[i] = &com.seq; }

		double *jd = g_new(double, e->num_frames);
		if (!mpp_derot_frame_times(a->seqs[i], e->fps, NAN, jd)) {
			siril_log_error(_("Derotation combine: '%s' has no usable timestamps "
			                  "(set its fps)\n"), e->seqname);
			g_free(jd); ok = FALSE; break;
		}
		a->derots[i] = mpp_derot_build(session->body, session->system, epoch_jd,
		                               jd, e->num_frames, e->frame_rows, e->frame_cols,
		                               e->cx, e->cy, e->r_eq, e->pa_deg, e->parity,
		                               NAN, NAN, NAN);
		g_free(jd);
		if (!a->derots[i]) {
			siril_log_error(_("Derotation combine: plan build failed for '%s'\n"),
			                e->seqname);
			ok = FALSE; break;
		}
		/* honour the fitted polar radius (oblateness) as drawn */
		if (e->r_pol > 0.0 && e->r_pol < e->r_eq)
			a->derots[i]->flattening = 1.0 - e->r_pol / e->r_eq;
	}
	if (!ok) {
		derot_combine_args_free(a);
		set_status(_("Combine setup failed — see the log."));
		return;
	}

	set_status(_("Combining sequences…"));
	control_window_switch_to_tab(OUTPUT_LOGS);
	set_progress_bar_data(_("Combining sequences…"), PROGRESS_RESET);
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(derot_combine_worker, a)) {
		set_cursor_waiting(FALSE);
		set_progress_bar_data(_("Ready."), PROGRESS_RESET);
		derot_combine_args_free(a);
		set_status(_("Could not start the combine worker (another job running?)."));
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
