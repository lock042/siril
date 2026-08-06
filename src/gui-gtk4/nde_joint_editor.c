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
#include "gui-gtk4/utils.h"
#include "registration/flis_register.h"
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
	int op_version = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == record_id) {
			op_id = g_strdup(rec->op_id);
			params = g_strdup(rec->params);
			op_version = rec->op_version;
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
	/* The RECORD's version, not the descriptor's: a deserializer reads it as
	 * "the format this blob was written in" and branches on it (photometric_cc
	 * v1 and v2 have different contracts).  Handing it the current version
	 * makes every old record claim to be new. */
	gpointer data = op && op->deserialize ? op->deserialize(params, op_version) : NULL;
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

/* ======================================================================= */
/* flis.register                                                           */
/* ======================================================================= */

/* The live "Register layers" dialog (flis_gui.c) is not reusable here: it is
 * built around the CURRENT layer stack — it scopes itself to the selected
 * group, seeds its reference from flis_active_layer(), validates the method
 * against com.selection, saves a multi-layer undo state and dispatches a
 * worker.  An amend edits a RECORDED operation whose participants and
 * selection are fixed in the blob, and applies through nde_amend_start.  So
 * this is a dedicated window, in the shape nde_joint_open_amend uses, sharing
 * the value tables rather than the widgets. */

/* Mirrors the live dialog's tables so the two offer identical choices. */
static const flis_reg_method_id REG_METHOD_VALUES[] = {
	FLIS_REG_GLOBAL, FLIS_REG_2PASS, FLIS_REG_DFT, FLIS_REG_KOMBAT
};
static const char *REG_METHOD_NAMES[] = {
	N_("Global star alignment (1-pass, recommended)"),
	N_("Global star alignment (2-pass)"),
	N_("Image pattern (DFT shift — needs a square selection)"),
	N_("KOMBAT pattern match (needs a selection)"),
};
static const opencv_interpolation REG_INTERP_VALUES[] = {
	OPENCV_NEAREST, OPENCV_LINEAR, OPENCV_CUBIC, OPENCV_AREA, OPENCV_LANCZOS4
};
static const char *REG_INTERP_NAMES[] = {
	N_("Nearest"), N_("Linear"), N_("Cubic"), N_("Area"), N_("Lanczos4")
};

/* The transformation type is normally the method's own (DFT and KOMBAT solve
 * a shift and nothing else), but a record can carry any of these, so the
 * dropdown offers the range and the method change re-seeds it. */
static const transformation_type REG_TX_VALUES[] = {
	SHIFT_TRANSFORMATION, SIMILARITY_TRANSFORMATION, AFFINE_TRANSFORMATION,
	HOMOGRAPHY_TRANSFORMATION
};
static const char *REG_TX_NAMES[] = {
	N_("Shift"), N_("Similarity"), N_("Affine"), N_("Homography")
};

typedef struct {
	gint64 record_id;
	struct nde_joint_register_data *reg;   /* owned typed params */
	GtkWidget *method_dd, *tx_dd, *interp_dd, *ref_dd, *clamp_toggle;
	GArray    *ref_ids;                    /* parallel to ref_dd items */
	GtkWidget *window;
} register_editor;

static void register_editor_free(gpointer p) {
	register_editor *ed = p;
	if (!ed)
		return;
	if (ed->reg)
		nde_joint_register_data_free(ed->reg);
	if (ed->ref_ids)
		g_array_unref(ed->ref_ids);
	g_free(ed);
}

/* Index of @value in @values, or @fallback when the record carries something
 * the dropdown does not list (a blob from a newer build, or a hand-edited
 * one).  Never leaves the user staring at a control that misreports the
 * record. */
static guint index_of_int(const gint *values, guint n, gint value, guint fallback) {
	for (guint i = 0; i < n; i++)
		if (values[i] == value)
			return i;
	return fallback;
}

static void on_register_cancel(GtkButton *btn, gpointer user) {
	(void)btn;
	gtk_window_destroy(GTK_WINDOW(((register_editor *)user)->window));
}

/* Re-seed the transformation dropdown from the method's natural choice: DFT
 * and KOMBAT only ever solve a shift, so silently leaving HOMOGRAPHY selected
 * after switching to one of them would record a setting the solve ignores. */
static void on_register_method_changed(GObject *obj, GParamSpec *pspec,
                                       gpointer user) {
	(void)obj; (void)pspec;
	register_editor *ed = user;
	guint idx = gtk_drop_down_get_selected(GTK_DROP_DOWN(ed->method_dd));
	if (idx >= G_N_ELEMENTS(REG_METHOD_VALUES))
		return;
	transformation_type tx = HOMOGRAPHY_TRANSFORMATION;
	(void)flis_register_resolve_method(REG_METHOD_VALUES[idx], NULL, &tx);
	for (guint i = 0; i < G_N_ELEMENTS(REG_TX_VALUES); i++)
		if (REG_TX_VALUES[i] == tx) {
			gtk_drop_down_set_selected(GTK_DROP_DOWN(ed->tx_dd), i);
			break;
		}
}

/* Apply = amend with the record's own params, settings replaced and nothing
 * else touched.  Re-serializing the DESERIALIZED struct is what guarantees
 * that: the transforms, framing, positions, selection and canvas size go back
 * out exactly as they came in, and no kv key is ever hand-assembled. */
static void on_register_apply(GtkButton *btn, gpointer user) {
	(void)btn;
	register_editor *ed = user;
	const op_descriptor *op = op_descriptor_by_id("flis.register");
	if (!op || !op->serialize)
		return;

	guint i = gtk_drop_down_get_selected(GTK_DROP_DOWN(ed->method_dd));
	gint method = i < G_N_ELEMENTS(REG_METHOD_VALUES)
	            ? (gint)REG_METHOD_VALUES[i] : ed->reg->method;
	i = gtk_drop_down_get_selected(GTK_DROP_DOWN(ed->tx_dd));
	gint tx = i < G_N_ELEMENTS(REG_TX_VALUES)
	        ? (gint)REG_TX_VALUES[i] : ed->reg->tx_type;
	i = gtk_drop_down_get_selected(GTK_DROP_DOWN(ed->interp_dd));
	gint interp = i < G_N_ELEMENTS(REG_INTERP_VALUES)
	            ? (gint)REG_INTERP_VALUES[i] : ed->reg->interpolation;
	gboolean clamp = siril_toggle_get_active(ed->clamp_toggle);
	i = gtk_drop_down_get_selected(GTK_DROP_DOWN(ed->ref_dd));
	gint ref_item = i < ed->ref_ids->len ? g_array_index(ed->ref_ids, gint, i)
	                                     : ed->reg->ref_item;

	/* TRUE means a setting moved, so the stored transforms no longer answer
	 * the question and the participants' signatures have been poisoned —
	 * replay will re-solve.  FALSE means nothing moved: the blob is
	 * byte-identical and the amend is the plain "re-run" verb. */
	if (nde_joint_register_apply_settings(ed->reg, method, tx, interp, clamp,
	                                      ref_item))
		siril_log_message(_("Register layers: settings changed — the "
		                    "alignment will be solved again.\n"));

	gchar *blob = op->serialize(ed->reg);
	nde_amend_start(ed->record_id, blob);
	g_free(blob);
	gtk_window_destroy(GTK_WINDOW(ed->window));
}

/* A labelled dropdown row, appended to @grid at @row. */
static GtkWidget *register_dd_row(GtkWidget *grid, int row, const char *label,
                                  const char **names, guint n_names,
                                  guint selected) {
	GtkWidget *l = gtk_label_new(label);
	gtk_widget_set_halign(l, GTK_ALIGN_START);
	gtk_grid_attach(GTK_GRID(grid), l, 0, row, 1, 1);
	GtkWidget *dd = gtk_drop_down_new(NULL, NULL);
	gtk_widget_set_hexpand(dd, TRUE);
	for (guint i = 0; i < n_names; i++)
		siril_drop_down_append_text(GTK_DROP_DOWN(dd), _(names[i]));
	gtk_drop_down_set_selected(GTK_DROP_DOWN(dd), selected);
	gtk_grid_attach(GTK_GRID(grid), dd, 1, row, 1, 1);
	return dd;
}

gboolean nde_register_open_amend(gint64 record_id) {
	gchar *params = NULL;
	int op_version = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id == record_id) {
			params = g_strdup(rec->params);
			op_version = rec->op_version;
			break;
		}
	}
	if (snap)
		g_ptr_array_unref(snap);
	if (!params)
		return FALSE;

	const op_descriptor *op = op_descriptor_by_id("flis.register");
	/* The RECORD's version, not the descriptor's — see nde_joint_open_amend. */
	struct nde_joint_register_data *reg =
			op && op->deserialize ? op->deserialize(params, op_version) : NULL;
	g_free(params);
	if (!reg)
		return FALSE;   /* kv-grid fallback can still show the blob */

	register_editor *ed = g_new0(register_editor, 1);
	ed->record_id = record_id;
	ed->reg = reg;
	ed->ref_ids = g_array_new(FALSE, FALSE, sizeof(gint));

	GtkWidget *win = gtk_window_new();
	ed->window = win;
	gtk_window_set_title(GTK_WINDOW(win), _("Register layers"));
	gtk_window_set_modal(GTK_WINDOW(win), TRUE);
	GtkWidget *main_w = lookup_widget("control_window");
	if (main_w)
		gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(main_w));
	gtk_window_set_default_size(GTK_WINDOW(win), 420, -1);

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_widget_set_margin_start(box, 12);
	gtk_widget_set_margin_end(box, 12);
	gtk_widget_set_margin_top(box, 12);
	gtk_widget_set_margin_bottom(box, 12);
	gtk_window_set_child(GTK_WINDOW(win), box);

	GtkWidget *grid = gtk_grid_new();
	gtk_grid_set_column_spacing(GTK_GRID(grid), 12);
	gtk_grid_set_row_spacing(GTK_GRID(grid), 4);

	gint mvals[G_N_ELEMENTS(REG_METHOD_VALUES)];
	for (guint i = 0; i < G_N_ELEMENTS(REG_METHOD_VALUES); i++)
		mvals[i] = (gint)REG_METHOD_VALUES[i];
	ed->method_dd = register_dd_row(grid, 0, _("Method"), REG_METHOD_NAMES,
	                                G_N_ELEMENTS(REG_METHOD_NAMES),
	                                index_of_int(mvals,
	                                             G_N_ELEMENTS(REG_METHOD_VALUES),
	                                             reg->method, 0));

	gint tvals[G_N_ELEMENTS(REG_TX_VALUES)];
	for (guint i = 0; i < G_N_ELEMENTS(REG_TX_VALUES); i++)
		tvals[i] = (gint)REG_TX_VALUES[i];
	ed->tx_dd = register_dd_row(grid, 1, _("Transformation"), REG_TX_NAMES,
	                            G_N_ELEMENTS(REG_TX_NAMES),
	                            index_of_int(tvals, G_N_ELEMENTS(REG_TX_VALUES),
	                                         reg->tx_type,
	                                         G_N_ELEMENTS(REG_TX_VALUES) - 1));

	gint ivals[G_N_ELEMENTS(REG_INTERP_VALUES)];
	for (guint i = 0; i < G_N_ELEMENTS(REG_INTERP_VALUES); i++)
		ivals[i] = (gint)REG_INTERP_VALUES[i];
	ed->interp_dd = register_dd_row(grid, 2, _("Interpolation"),
	                                REG_INTERP_NAMES,
	                                G_N_ELEMENTS(REG_INTERP_NAMES),
	                                index_of_int(ivals,
	                                             G_N_ELEMENTS(REG_INTERP_VALUES),
	                                             reg->interpolation,
	                                             G_N_ELEMENTS(REG_INTERP_VALUES) - 1));

	/* Reference: the record's OWN participants, not the live layer list — the
	 * participating layers are fixed by the record and a reference outside
	 * them has no pixels to align to. */
	GtkWidget *rl = gtk_label_new(_("Reference layer"));
	gtk_widget_set_halign(rl, GTK_ALIGN_START);
	gtk_grid_attach(GTK_GRID(grid), rl, 0, 3, 1, 1);
	ed->ref_dd = gtk_drop_down_new(NULL, NULL);
	gtk_widget_set_hexpand(ed->ref_dd, TRUE);
	guint ref_idx = 0;
	for (guint k = 0; k < reg->n; k++) {
		const nde_joint_register_participant *part = &reg->parts[k];
		siril_drop_down_append_text(GTK_DROP_DOWN(ed->ref_dd),
		                            part->name && *part->name ? part->name
		                                                      : _("(unnamed layer)"));
		g_array_append_val(ed->ref_ids, reg->parts[k].item_id);
		if (part->item_id == reg->ref_item)
			ref_idx = k;
	}
	gtk_drop_down_set_selected(GTK_DROP_DOWN(ed->ref_dd), ref_idx);
	gtk_grid_attach(GTK_GRID(grid), ed->ref_dd, 1, 3, 1, 1);
	gtk_box_append(GTK_BOX(box), grid);

	ed->clamp_toggle = gtk_check_button_new_with_label(
			_("Clamp cubic/Lanczos to suppress ringing"));
	siril_toggle_set_active(ed->clamp_toggle, reg->clamp);
	gtk_box_append(GTK_BOX(box), ed->clamp_toggle);

	g_signal_connect(ed->method_dd, "notify::selected",
	                 G_CALLBACK(on_register_method_changed), ed);

	GtkWidget *note = gtk_label_new(
			_("Changing any setting here discards the stored alignment: the "
			  "layers are aligned again with the new settings when the step "
			  "replays.  Leaving them alone re-runs the step unchanged.  "
			  "Changes apply on OK — there is no live preview.  The "
			  "participating layers are fixed; to change them, delete this "
			  "step and register again."));
	gtk_label_set_wrap(GTK_LABEL(note), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(note), 52);
	gtk_widget_add_css_class(note, "dim-label");
	gtk_box_append(GTK_BOX(box), note);

	GtkWidget *btns = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_widget_set_halign(btns, GTK_ALIGN_END);
	GtkWidget *cancel = gtk_button_new_with_label(_("Cancel"));
	GtkWidget *apply = gtk_button_new_with_label(_("OK"));
	gtk_widget_add_css_class(apply, "suggested-action");
	gtk_box_append(GTK_BOX(btns), cancel);
	gtk_box_append(GTK_BOX(btns), apply);
	gtk_box_append(GTK_BOX(box), btns);

	g_signal_connect(cancel, "clicked", G_CALLBACK(on_register_cancel), ed);
	g_signal_connect(apply, "clicked", G_CALLBACK(on_register_apply), ed);
	g_object_set_data_full(G_OBJECT(win), "register-editor", ed,
	                       register_editor_free);

	gtk_window_present(GTK_WINDOW(win));
	return TRUE;
}
