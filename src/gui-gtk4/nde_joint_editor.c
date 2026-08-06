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

/* See nde_joint_editor.h for the contract. */

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/nde_history.h"
#include "core/nde_joint.h"
#include "core/nde_replay.h"
#include "core/op_descriptor.h"
#include "gui-gtk4/nde_joint_editor.h"
#include "gui-gtk4/photometric_cc.h"
#include <gtk/gtk.h>

extern GtkWidget *lookup_widget(const gchar *widget_name);

typedef struct {
	gint64     record_id;
	gboolean   is_group;
	/* Owned typed params, one of the two set. */
	struct nde_joint_layers_match_data *match;
	struct nde_joint_group_calib_data  *calib;
	/* Editable widgets (group CC only; NULL otherwise). */
	GtkWidget *spin_white[4];   /* x y w h */
	GtkWidget *spin_black[4];
	GtkWidget *spin_low, *spin_high;
	GtkWidget *spin_mkw[3];
	GtkWidget *window;
} joint_editor;

static void joint_editor_free(gpointer p) {
	joint_editor *ed = p;
	if (!ed)
		return;
	if (ed->match)
		nde_joint_layers_match_data_free(ed->match);
	if (ed->calib)
		nde_joint_group_calib_data_free(ed->calib);
	g_free(ed);
}

static const nde_joint_participant *editor_parts(const joint_editor *ed, guint *n) {
	if (ed->match) {
		*n = ed->match->n;
		return ed->match->parts;
	}
	*n = ed->calib->n;
	return ed->calib->parts;
}

static void on_cancel(GtkButton *btn, gpointer user) {
	(void)btn;
	gtk_window_destroy(GTK_WINDOW(((joint_editor *)user)->window));
}

/* Apply = amend with the (possibly updated) operation parameters.  The core
 * recomputes the analysis and replays every participant; an unchanged blob
 * is exactly the "re-run" verb. */
static void on_apply(GtkButton *btn, gpointer user) {
	(void)btn;
	joint_editor *ed = user;
	const op_descriptor *op = op_descriptor_by_id(
			ed->is_group ? "flis.group_calibration" : "flis.layers_match");
	if (!op || !op->serialize)
		return;
	if (ed->calib && ed->spin_low) {
		ed->calib->white_sel.x = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_white[0]));
		ed->calib->white_sel.y = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_white[1]));
		ed->calib->white_sel.w = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_white[2]));
		ed->calib->white_sel.h = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_white[3]));
		ed->calib->black_sel.x = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_black[0]));
		ed->calib->black_sel.y = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_black[1]));
		ed->calib->black_sel.w = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_black[2]));
		ed->calib->black_sel.h = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ed->spin_black[3]));
		ed->calib->low  = gtk_spin_button_get_value(GTK_SPIN_BUTTON(ed->spin_low));
		ed->calib->high = gtk_spin_button_get_value(GTK_SPIN_BUTTON(ed->spin_high));
	} else if (ed->calib && ed->spin_mkw[0]) {
		for (int c = 0; c < 3; c++)
			ed->calib->manual_kw[c] =
					gtk_spin_button_get_value(GTK_SPIN_BUTTON(ed->spin_mkw[c]));
	}
	gchar *blob = op->serialize(ed->calib ? (gpointer)ed->calib
	                                      : (gpointer)ed->match);
	nde_amend_start(ed->record_id, blob);
	g_free(blob);
	gtk_window_destroy(GTK_WINDOW(ed->window));
}

static GtkWidget *rect_spin_row(GtkWidget *grid, int row, const char *label,
                                GtkWidget *spins[4], const rectangle *r) {
	GtkWidget *l = gtk_label_new(label);
	gtk_widget_set_halign(l, GTK_ALIGN_START);
	gtk_grid_attach(GTK_GRID(grid), l, 0, row, 1, 1);
	const int vals[4] = { r->x, r->y, r->w, r->h };
	for (int i = 0; i < 4; i++) {
		spins[i] = gtk_spin_button_new_with_range(0, 100000, 1);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spins[i]), vals[i]);
		gtk_grid_attach(GTK_GRID(grid), spins[i], 1 + i, row, 1, 1);
	}
	return l;
}

gboolean nde_joint_open_amend(gint64 record_id) {
	/* Typed params from the live record. */
	gchar *op_id = NULL, *params = NULL;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == record_id) {
			op_id = g_strdup(rec->op_id);
			params = g_strdup(rec->params);
			break;
		}
	}
	if (snap)
		g_ptr_array_unref(snap);
	if (!op_id || !params) {
		g_free(op_id);
		g_free(params);
		return FALSE;
	}
	const op_descriptor *op = op_descriptor_by_id(op_id);
	gpointer data = op && op->deserialize ? op->deserialize(params, op->version) : NULL;
	gboolean is_group = !g_strcmp0(op_id, "flis.group_calibration");
	g_free(params);
	if (!data) {
		g_free(op_id);
		return FALSE;   /* kv-grid fallback can still show the blob */
	}

	/* Group PCC/SPCC: the photometric parameters are the REAL dialog's to
	 * edit, pre-filled from the record's nested params — the same editing
	 * experience the single-image record gets.  Layers match (no parameters
	 * of its own) and group CC (a pair of selections) keep this window. */
	if (is_group) {
		struct nde_joint_group_calib_data *gd = data;
		if (gd->kind == 1 || gd->kind == 2) {
			nde_joint_group_calib_data_free(gd);
			g_free(op_id);
			return pcc_open_amend(record_id);
		}
	}

	joint_editor *ed = g_new0(joint_editor, 1);
	ed->record_id = record_id;
	ed->is_group = is_group;
	if (is_group)
		ed->calib = data;
	else
		ed->match = data;

	GtkWidget *win = gtk_window_new();
	ed->window = win;
	gchar *log_label = op->log_hook ? op->log_hook(data, 0) : NULL;
	gtk_window_set_title(GTK_WINDOW(win),
	                     log_label ? log_label : _("Joint calibration"));
	g_free(log_label);
	gtk_window_set_modal(GTK_WINDOW(win), TRUE);
	GtkWidget *main_w = lookup_widget("control_window");
	if (main_w)
		gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(main_w));

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_widget_set_margin_start(box, 12);
	gtk_widget_set_margin_end(box, 12);
	gtk_widget_set_margin_top(box, 12);
	gtk_widget_set_margin_bottom(box, 12);
	gtk_window_set_child(GTK_WINDOW(win), box);

	/* Participants + the factors the analysis last derived (read-only:
	 * they are RESULTS — re-running is how they change). */
	GtkWidget *pgrid = gtk_grid_new();
	gtk_grid_set_column_spacing(GTK_GRID(pgrid), 12);
	gtk_grid_set_row_spacing(GTK_GRID(pgrid), 2);
	guint n = 0;
	const nde_joint_participant *parts = editor_parts(ed, &n);
	for (guint i = 0; i < n; i++) {
		GtkWidget *name = gtk_label_new(parts[i].name && *parts[i].name ?
		                                parts[i].name : _("(unnamed layer)"));
		gtk_widget_set_halign(name, GTK_ALIGN_START);
		gchar *ftxt = parts[i].diag_offset != 0.0 ?
				g_strdup_printf("\303\227%.4f %+.5f", parts[i].diag_scale,
				                parts[i].diag_offset) :
				g_strdup_printf("\303\227%.4f", parts[i].diag_scale);
		GtkWidget *fac = gtk_label_new(ftxt);
		g_free(ftxt);
		gtk_widget_set_halign(fac, GTK_ALIGN_END);
		gtk_grid_attach(GTK_GRID(pgrid), name, 0, (int)i, 1, 1);
		gtk_grid_attach(GTK_GRID(pgrid), fac, 1, (int)i, 1, 1);
	}
	gtk_box_append(GTK_BOX(box), pgrid);

	/* The operation's own editable parameters, where the record has any. */
	if (ed->calib && ed->calib->kind == 0 && !ed->calib->manual) {
		GtkWidget *grid = gtk_grid_new();
		gtk_grid_set_column_spacing(GTK_GRID(grid), 6);
		gtk_grid_set_row_spacing(GTK_GRID(grid), 4);
		rect_spin_row(grid, 0, _("White reference area"), ed->spin_white,
		              &ed->calib->white_sel);
		rect_spin_row(grid, 1, _("Background area"), ed->spin_black,
		              &ed->calib->black_sel);
		GtkWidget *ll = gtk_label_new(_("Rejection limits"));
		gtk_widget_set_halign(ll, GTK_ALIGN_START);
		gtk_grid_attach(GTK_GRID(grid), ll, 0, 2, 1, 1);
		ed->spin_low = gtk_spin_button_new_with_range(0.0, 1.0, 0.001);
		gtk_spin_button_set_digits(GTK_SPIN_BUTTON(ed->spin_low), 3);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(ed->spin_low), ed->calib->low);
		gtk_grid_attach(GTK_GRID(grid), ed->spin_low, 1, 2, 1, 1);
		ed->spin_high = gtk_spin_button_new_with_range(0.0, 1.0, 0.001);
		gtk_spin_button_set_digits(GTK_SPIN_BUTTON(ed->spin_high), 3);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(ed->spin_high), ed->calib->high);
		gtk_grid_attach(GTK_GRID(grid), ed->spin_high, 2, 2, 1, 1);
		gtk_box_append(GTK_BOX(box), grid);
	} else if (ed->calib && ed->calib->kind == 0 && ed->calib->manual) {
		GtkWidget *grid = gtk_grid_new();
		gtk_grid_set_column_spacing(GTK_GRID(grid), 6);
		GtkWidget *ml = gtk_label_new(_("Manual white balance"));
		gtk_widget_set_halign(ml, GTK_ALIGN_START);
		gtk_grid_attach(GTK_GRID(grid), ml, 0, 0, 1, 1);
		for (int c = 0; c < 3; c++) {
			ed->spin_mkw[c] = gtk_spin_button_new_with_range(0.0, 10.0, 0.00001);
			gtk_spin_button_set_digits(GTK_SPIN_BUTTON(ed->spin_mkw[c]), 5);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(ed->spin_mkw[c]),
			                          ed->calib->manual_kw[c]);
			gtk_grid_attach(GTK_GRID(grid), ed->spin_mkw[c], 1 + c, 0, 1, 1);
		}
		gtk_box_append(GTK_BOX(box), grid);
	}

	GtkWidget *note = gtk_label_new(
			_("Re-running recomputes the analysis for every listed layer "
			  "against the current state of the composition.  Changes apply "
			  "on Re-run — there is no live preview.  The participating "
			  "layers are fixed; to change them, delete this step and run "
			  "the operation again."));
	gtk_label_set_wrap(GTK_LABEL(note), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(note), 52);
	gtk_widget_add_css_class(note, "dim-label");
	gtk_box_append(GTK_BOX(box), note);

	GtkWidget *btns = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_widget_set_halign(btns, GTK_ALIGN_END);
	GtkWidget *cancel = gtk_button_new_with_label(_("Cancel"));
	GtkWidget *apply = gtk_button_new_with_label(_("Re-run analysis"));
	gtk_widget_add_css_class(apply, "suggested-action");
	gtk_box_append(GTK_BOX(btns), cancel);
	gtk_box_append(GTK_BOX(btns), apply);
	gtk_box_append(GTK_BOX(box), btns);

	g_signal_connect(cancel, "clicked", G_CALLBACK(on_cancel), ed);
	g_signal_connect(apply, "clicked", G_CALLBACK(on_apply), ed);
	g_object_set_data_full(G_OBJECT(win), "joint-editor", ed, joint_editor_free);

	gtk_window_present(GTK_WINDOW(win));
	g_free(op_id);
	return TRUE;
}
