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
 * PixelMath dialog — widget management, callbacks, expression builder, and
 * file-chooser logic.  The computation engine lives in
 * pixelMath/pixel_math_runner.c.
 */

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "core/initfile.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/dialog_preview.h"
#include "gui-gtk4/open_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/gui_state.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "pixelMath/pixel_math_runner.h"

/* ── Operator / function tables (GUI-side; used by tooltips and list views) */

typedef struct {
	const gchar *name;
	const gchar *prototype;
	const gchar *definition;
} _pm_op_func;

static const _pm_op_func functions[] = {
    { "abs",   "abs ( x )",                             N_("Absolute value of x.")                    },
    { "acos",  "acos ( x )",                            N_("Arc cosine of x.")                        },
    { "acosh", "acosh ( x )",                           N_("Hyperbolic arc cosine of x.")             },
    { "asin",  "asin ( x )",                            N_("Arc sine of x.")                          },
    { "asinh", "asinh ( x )",                           N_("Hyperbolic arc sine of x.")               },
    { "atan",  "atan ( x )",                            N_("Arc tangent of x.")                       },
    { "atan2", "atan2 ( y, x )",                        N_("Arc tangent of y/x.")                     },
    { "atanh", "atanh ( x )",                           N_("Hyperbolic arc tangent of x.")            },
    { "ceil",  "ceil ( x )",                            N_("Round x upwards to the nearest integer.") },
    { "cos",   "cos ( x )",                             N_("Cosine of x.")                            },
    { "cosh",  "cosh ( x )",                            N_("Hyperbolic cosine of x.")                 },
    { "e",     "e",                                     N_("The constant e=2.718282...")              },
    { "exp",   "exp ( x )",                             N_("Exponential function.")                   },
    { "fac",   "fac( x )",                              N_("Factorial function.")                     },
    { "iif",   "iif( cond, expr_true, expr_false )",    N_("Conditional function (or inline if function).\n"
                                        "\nReturns <i>expr_true</i> if <i>cond</i> evaluates to nonzero."
                                        "\nReturns <i>expr_false</i> if <i>cond</i> evaluates to zero.") },
    { "floor", "floor ( x )",                           N_("Highest integer less than or equal to x.")},
    { "ln",    "ln ( x )",                              N_("Natural logarithm of x.")                 },
    { "log",   "log ( x )",                             N_("Base-10 logarithm of x.")                 },
    { "log10", "log10 ( x )",                           N_("Base-10 logarithm of x.")                 },
    { "log2",  "log2 ( x )",                            N_("Base-2 logarithm of x.")                  },
    { "max",   "max ( x, y )",                          N_("Maximum function.")                       },
    { "min",   "min ( x, y )",                          N_("Minimum function.")                       },
    { "mtf",   "mtf ( m, x )",                          N_("Midtones Transfer Function (MTF) of x for a midtones balance parameter m in the [0, 1] range.") },
    { "ncr",   "ncr ( x, y )",                          N_("Combinations function.")                  },
    { "npr",   "npr ( x, y )",                          N_("Permutations function.")                  },
    { "pi",    "pi",                                    N_("The constant π=3.141592...")         },
    { "pow",   "pow ( x, y )",                          N_("Exponentiation function.")                },
    { "sign",  "sign ( x )",                            N_("Sign of x:\n"
                                                                        "\n\t +1 if x &gt; 0"
                                                                        "\n\t −1 if x &lt; 0"
                                                                        "\n\t  0 if x = 0.")  },
    { "sin",   "sin ( x )",                             N_("Sine of x.")                              },
    { "sinh",  "sinh ( x )",                            N_("Hyperbolic sine of x.")                   },
    { "sqrt",  "sqrt ( x )",                            N_("Square root of x.")                       },
    { "tan",   "tan ( x )",                             N_("Tangent of x.")                           },
    { "tanh",  "tanh ( x )",                            N_("Hyperbolic tangent of x.")                },
    { "trunc", "trunc ( x )",                           N_("Truncated integer part of x.")            }
};

static const _pm_op_func image_functions[] = {
	{ "adev",     "adev ( Image )",                    N_("Average absolute deviation of the image.")  },
	{ "bwmv",     "bwmv ( Image )",                    N_("Biweight midvariance of the image.")        },
	{ "height",   "height ( Image )",                  N_("Height in pixels of the specified image.")  },
	{ "mad",      "mad ( Image )",                     N_("Median absolute deviation of the image.")   },
	{ "max",      "max ( Image )",                     N_("Pixel maximum of the image.")               },
	{ "mdev",     "mdev ( Image )",                    N_("Median absolute deviation of the image.")   },
	{ "mean",     "mean ( Image )",                    N_("Mean of the image.")                        },
	{ "med",      "med ( Image )",                     N_("Median of the image.")                      },
	{ "median",   "median ( Image )",                  N_("Median of the image.")                      },
	{ "min",      "min ( Image )",                     N_("Pixel minimum of the image.")               },
	{ "noise",    "noise ( Image )",                   N_("Estimation of Gaussian noise in the image.")},
	{ "sdev",     "sdev ( Image )",                    N_("Standard deviation of the image.")          },
	{ "width",    "width ( Image )",                   N_("Width in pixels of the specified image.")   }
};

static const _pm_op_func operators[] = {
    { "~",   "~x",                             N_("Pixel Inversion operator.")                        },
    { "-",   "-x",                             N_("Unary Minus operator (sign change).")              },
    { "+",   "+x",                             N_("Unary Plus operator.")                             },
    { "!",   "!x",                             N_("Logical NOT operator.")                            },
    { "^",   "x ^ y",                          N_("Exponentiation operator.")                         },
    { "*",   "x * y",                          N_("Multiplication operator.")                         },
    { "/",   "x / y",                          N_("Division operator.")                               },
    { "%",   "x % y",                          N_("Modulus operator.")                                },
    { "+",   "x + y",                          N_("Addition operator.")                               },
    { "-",   "x - y",                          N_("Subtraction operator.")                            },
    { "<",   "x &lt; y",                       N_("Less Than relational operator.")                   },
    { "<=",  "x &lt;= y",                      N_("Less Than Or Equal relational operator.")          },
    { ">",   "x &gt; y",                       N_("Greater Than relational operator.")                },
    { ">=",  "x &gt;= y",                      N_("Greater Than Or Equal relational operator.")       },
    { "==",  "x == y",                         N_("Equal To relational operator.")                    },
    { "!=",  "x != y",                         N_("Not Equal To relational operator.")                },
    { "&&",  "x &amp;&amp; y",                 N_("Logical AND operator.")                            },
    { "||",  "x || y",                         N_("Logical OR operator.")                             }
};

#define MAX_FUNCTIONS      G_N_ELEMENTS(functions)
#define MAX_IMAGE_FUNCTIONS G_N_ELEMENTS(image_functions)
#define MAX_OPERATORS      G_N_ELEMENTS(operators)

/* ── Column indices for the image-variable tree view ─────────────────────── */

enum {
	COLUMN_IMAGE_NUM,
	COLUMN_IMAGE_PATH,
	COLUMN_IMAGE_CHAN,
	COLUMN_IMAGE_WIDTH,
	COLUMN_IMAGE_HEIGHT,
	N_COLUMNS
};

enum {
	COLUMN_NAME,
	COLUMN_INDEX
};

/* ── Static widget pointers ─────────────────────────────────────────────── */

/* Phase 11 list 4: presets — GtkStringList + GtkColumnView. */
static GtkStringList *pm_preset_store     = NULL;
static GtkColumnView *pm_preset_view      = NULL;
static GtkMultiSelection *pm_preset_sel   = NULL;
/* Phase 11 list 3: operators — GtkStringList + GtkColumnView. */
static GtkStringList *pm_op_store         = NULL;
static GtkColumnView *pm_op_view          = NULL;
/* Phase 11 list 2: functions — GtkStringList + GtkColumnView, with a
 * cached sorted (functions ++ image_functions) array used by the bind
 * callback for tooltips and by activate to read names. */
static GtkStringList *pm_fn_store         = NULL;
static GtkColumnView *pm_fn_view          = NULL;
static _pm_op_func   *pm_fn_sorted        = NULL;
static int            pm_fn_sorted_count  = 0;
/* Phase 11 list 1: variables — GListStore<SirilPmVarRow*> + GtkColumnView. */
static GListStore   *pm_var_store         = NULL;
static GtkColumnView *pm_var_view         = NULL;
static GtkMultiSelection *pm_var_sel      = NULL;
static GtkLabel     *pixel_math_status_bar = NULL;
static GtkLabel     *pixel_math_status_bar2 = NULL;
static GtkEntry     *pixel_math_entry_r = NULL;
static GtkEntry     *pixel_math_entry_g = NULL;
static GtkEntry     *pixel_math_entry_b = NULL;

static int entry_has_focus = 0;

/* ── Phase 11 list 1 (variables): per-row GObject ────────────────────────
 * 5 fields matching the original GtkListStore columns: var name, path,
 * channel count, width, height.  Variable name is editable; path is
 * read-only.
 */
#define SIRIL_TYPE_PM_VAR_ROW (siril_pm_var_row_get_type())
G_DECLARE_FINAL_TYPE(SirilPmVarRow, siril_pm_var_row, SIRIL, PM_VAR_ROW, GObject)
struct _SirilPmVarRow {
	GObject parent_instance;
	gchar *var;
	gchar *path;
	guint  chan;
	guint  width;
	guint  height;
};
G_DEFINE_TYPE(SirilPmVarRow, siril_pm_var_row, G_TYPE_OBJECT)
static void siril_pm_var_row_finalize(GObject *obj) {
	SirilPmVarRow *self = SIRIL_PM_VAR_ROW(obj);
	g_clear_pointer(&self->var,  g_free);
	g_clear_pointer(&self->path, g_free);
	G_OBJECT_CLASS(siril_pm_var_row_parent_class)->finalize(obj);
}
static void siril_pm_var_row_class_init(SirilPmVarRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_pm_var_row_finalize;
}
static void siril_pm_var_row_init(SirilPmVarRow *self) { (void)self; }

/* ── Phase 11 list 4 (presets) helpers ───────────────────────────────────
 * Defined ahead of init_widgets() so the column view can be wired up there.
 * The presets store is a GtkStringList; rows are GtkStringObject and bind
 * to a GtkLabel via a GtkSignalListItemFactory.
 */

static void pm_string_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}
static void pm_string_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	GObject *o = gtk_list_item_get_item(li);
	const gchar *s = (o && GTK_IS_STRING_OBJECT(o))
			? gtk_string_object_get_string(GTK_STRING_OBJECT(o)) : "";
	gtk_label_set_text(lbl, s);
}

/* Forward declarations: defined further down; needed by init_widgets. */
static GtkEntry *get_entry_with_focus(void);
static const gchar *get_preset_expr(int i);
static const gchar *get_operator_name(int i);
static const gchar *get_function_name(int i);
static _pm_op_func *concat_functions(void);
static const gchar *get_pixel_math_var_name(int i);
static gboolean check_for_variable_sanity(const gchar *new_text);
static void output_status_bar2(int width, int height, int channel);
static gboolean on_pm_preset_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data);
static void remove_presets_from_list(void);

static void on_pm_preset_activate(GtkColumnView *cv, guint position, gpointer user_data) {
	(void)cv; (void)user_data;
	const gchar *str = get_preset_expr((int)position);
	if (!str) return;
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	guint pos = gtk_editable_get_position(GTK_EDITABLE(entry));
	gtk_editable_delete_selection(GTK_EDITABLE(entry));
	gtk_entry_buffer_insert_text(buffer, pos, str, -1);
	gtk_widget_grab_focus(GTK_WIDGET(entry));
	gtk_editable_set_position(GTK_EDITABLE(entry), pos + strlen(str));
}

/* Operators list: text is the operator symbol; tooltip uses prototype +
 * definition from the static operators[] array (same order as the model). */
static void pm_op_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	GObject *o = gtk_list_item_get_item(li);
	const gchar *s = (o && GTK_IS_STRING_OBJECT(o))
			? gtk_string_object_get_string(GTK_STRING_OBJECT(o)) : "";
	gtk_label_set_text(lbl, s);
	guint pos = gtk_list_item_get_position(li);
	if (pos < MAX_OPERATORS) {
		gchar *tip = g_strdup_printf("<b>%s</b>\n\n%s",
				operators[pos].prototype, _(operators[pos].definition));
		gtk_widget_set_tooltip_markup(GTK_WIDGET(lbl), tip);
		g_free(tip);
	} else {
		gtk_widget_set_tooltip_text(GTK_WIDGET(lbl), NULL);
	}
}

static void on_pm_op_activate(GtkColumnView *cv, guint position, gpointer user_data) {
	(void)cv; (void)user_data;
	const gchar *str = get_operator_name((int)position);
	if (!str) return;
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	guint pos = gtk_editable_get_position(GTK_EDITABLE(entry));
	gtk_editable_delete_selection(GTK_EDITABLE(entry));
	gtk_entry_buffer_insert_text(buffer, pos, str, -1);
	gtk_widget_grab_focus(GTK_WIDGET(entry));
	gtk_editable_set_position(GTK_EDITABLE(entry), pos + strlen(str));
}

/* Functions list: bound to a GtkStringList of names from pm_fn_sorted[];
 * tooltip pulls prototype + definition from the same array.  The store
 * order matches the array order so position == array index. */
static void pm_fn_bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	GObject *o = gtk_list_item_get_item(li);
	const gchar *s = (o && GTK_IS_STRING_OBJECT(o))
			? gtk_string_object_get_string(GTK_STRING_OBJECT(o)) : "";
	gtk_label_set_text(lbl, s);
	guint pos = gtk_list_item_get_position(li);
	if (pm_fn_sorted && (int)pos < pm_fn_sorted_count) {
		gchar *tip = g_strdup_printf("<b>%s</b>\n\n%s",
				pm_fn_sorted[pos].prototype, _(pm_fn_sorted[pos].definition));
		gtk_widget_set_tooltip_markup(GTK_WIDGET(lbl), tip);
		g_free(tip);
	} else {
		gtk_widget_set_tooltip_text(GTK_WIDGET(lbl), NULL);
	}
}

static void on_pm_fn_activate(GtkColumnView *cv, guint position, gpointer user_data) {
	(void)cv; (void)user_data;
	const gchar *str = get_function_name((int)position);
	if (!str) return;
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	guint pos = gtk_editable_get_position(GTK_EDITABLE(entry));
	gtk_editable_delete_selection(GTK_EDITABLE(entry));
	gtk_entry_buffer_insert_text(buffer, pos, str, -1);
	gtk_widget_grab_focus(GTK_WIDGET(entry));
	gtk_editable_set_position(GTK_EDITABLE(entry), pos + strlen(str));
}

/* Variable list cell factories (Variable + Path columns). */
static void pm_var_setup_path_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}
static void pm_var_bind_path_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilPmVarRow *row = SIRIL_PM_VAR_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, row && row->path ? row->path : "");
}

/* Variable column: GtkEditableLabel; commits on focus-out / activate via the
 * "editing" property notify.  If the new text fails check_for_variable_sanity
 * the change is reverted. */
static void on_pm_var_editing_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkEditableLabel *lbl = GTK_EDITABLE_LABEL(obj);
	if (gtk_editable_label_get_editing(lbl))
		return;
	GtkListItem *li = g_object_get_data(G_OBJECT(lbl), "siril-listitem");
	if (!li) return;
	SirilPmVarRow *row = SIRIL_PM_VAR_ROW(gtk_list_item_get_item(li));
	if (!row) return;
	const gchar *new_text = gtk_editable_get_text(GTK_EDITABLE(lbl));
	if (g_strcmp0(new_text, row->var) == 0) return;
	if (check_for_variable_sanity(new_text)) {
		g_free(row->var);
		row->var = g_strdup(new_text);
	} else {
		/* revert */
		gtk_editable_set_text(GTK_EDITABLE(lbl), row->var ? row->var : "");
	}
}

static void pm_var_setup_var_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_editable_label_new("");
	gtk_widget_set_halign(lbl, GTK_ALIGN_FILL);
	gtk_list_item_set_child(li, lbl);
	g_signal_connect(lbl, "notify::editing",
			G_CALLBACK(on_pm_var_editing_changed), NULL);
}
static void pm_var_bind_var_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkEditableLabel *lbl = GTK_EDITABLE_LABEL(gtk_list_item_get_child(li));
	SirilPmVarRow *row = SIRIL_PM_VAR_ROW(gtk_list_item_get_item(li));
	g_object_set_data(G_OBJECT(lbl), "siril-listitem", li);
	gtk_editable_set_text(GTK_EDITABLE(lbl), row && row->var ? row->var : "");
}
static void pm_var_unbind_var_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkEditableLabel *lbl = GTK_EDITABLE_LABEL(gtk_list_item_get_child(li));
	g_object_set_data(G_OBJECT(lbl), "siril-listitem", NULL);
}

static void on_pm_var_activate(GtkColumnView *cv, guint position, gpointer user_data) {
	(void)cv; (void)user_data;
	const gchar *str = get_pixel_math_var_name((int)position);
	if (!str) return;
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	guint pos = gtk_editable_get_position(GTK_EDITABLE(entry));
	gtk_editable_delete_selection(GTK_EDITABLE(entry));
	gtk_entry_buffer_insert_text(buffer, pos, str, -1);
	gtk_widget_grab_focus(GTK_WIDGET(entry));
	gtk_editable_set_position(GTK_EDITABLE(entry), pos + strlen(str));
}

static void on_pm_var_selection_changed(GtkSelectionModel *sm, guint position, guint n_items, gpointer user_data) {
	(void)position; (void)n_items; (void)user_data;
	guint width = 0, height = 0, channel = 0;
	GtkBitset *bs = gtk_selection_model_get_selection(sm);
	if (gtk_bitset_get_size(bs) == 1) {
		guint p = gtk_bitset_get_nth(bs, 0);
		SirilPmVarRow *row = SIRIL_PM_VAR_ROW(g_list_model_get_item(G_LIST_MODEL(pm_var_store), p));
		if (row) {
			channel = row->chan;
			width   = row->width;
			height  = row->height;
			g_object_unref(row);
		}
	}
	gtk_bitset_unref(bs);
	output_status_bar2(width, height, channel);
}

/* ── Widget initialisation ──────────────────────────────────────────────── */

static void init_widgets() {
	if (!pm_var_store) {
		pixel_math_status_bar = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_status")));
		pixel_math_status_bar2 = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_status2")));
		pixel_math_entry_r = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_r")));
		pixel_math_entry_g = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_g")));
		pixel_math_entry_b = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_b")));
#if GTK_CHECK_VERSION(3, 22, 0)
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_functions"))), TRUE);
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_operators"))), TRUE);
#endif

		/* ---- Phase 11 list 1 (variables): build GtkColumnView in C ---- */
		GtkScrolledWindow *sw_var = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_variables"));
		pm_var_store = g_list_store_new(SIRIL_TYPE_PM_VAR_ROW);
		pm_var_sel = gtk_multi_selection_new(G_LIST_MODEL(g_object_ref(pm_var_store)));
		g_signal_connect(pm_var_sel, "selection-changed", G_CALLBACK(on_pm_var_selection_changed), NULL);
		pm_var_view = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(g_object_ref(pm_var_sel))));
		g_signal_connect(pm_var_view, "activate", G_CALLBACK(on_pm_var_activate), NULL);

		GtkSignalListItemFactory *facvar = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
		g_signal_connect(facvar, "setup",  G_CALLBACK(pm_var_setup_var_cb),  NULL);
		g_signal_connect(facvar, "bind",   G_CALLBACK(pm_var_bind_var_cb),   NULL);
		g_signal_connect(facvar, "unbind", G_CALLBACK(pm_var_unbind_var_cb), NULL);
		GtkColumnViewColumn *cvar = gtk_column_view_column_new(N_("Variable"), GTK_LIST_ITEM_FACTORY(facvar));
		gtk_column_view_column_set_resizable(cvar, TRUE);
		gtk_column_view_append_column(pm_var_view, cvar);
		g_object_unref(cvar);

		GtkSignalListItemFactory *facpath = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
		g_signal_connect(facpath, "setup", G_CALLBACK(pm_var_setup_path_cb), NULL);
		g_signal_connect(facpath, "bind",  G_CALLBACK(pm_var_bind_path_cb),  NULL);
		GtkColumnViewColumn *cpath = gtk_column_view_column_new(N_("Path"), GTK_LIST_ITEM_FACTORY(facpath));
		gtk_column_view_column_set_resizable(cpath, TRUE);
		gtk_column_view_column_set_expand(cpath, TRUE);
		gtk_column_view_append_column(pm_var_view, cpath);
		g_object_unref(cpath);

		gtk_scrolled_window_set_child(sw_var, GTK_WIDGET(pm_var_view));

		/* ---- Phase 11 list 4 (presets): build GtkColumnView in C ---- */
		GtkScrolledWindow *sw_pr = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_presets"));
		pm_preset_store = gtk_string_list_new(NULL);
		pm_preset_sel = gtk_multi_selection_new(G_LIST_MODEL(g_object_ref(pm_preset_store)));
		pm_preset_view = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(g_object_ref(pm_preset_sel))));
		gtk_column_view_set_show_column_separators(pm_preset_view, FALSE);
		g_signal_connect(pm_preset_view, "activate", G_CALLBACK(on_pm_preset_activate), NULL);

		GtkSignalListItemFactory *fac = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
		g_signal_connect(fac, "setup", G_CALLBACK(pm_string_setup_cb), NULL);
		g_signal_connect(fac, "bind",  G_CALLBACK(pm_string_bind_cb),  NULL);
		GtkColumnViewColumn *c = gtk_column_view_column_new(N_("expression"), GTK_LIST_ITEM_FACTORY(fac));
		gtk_column_view_column_set_expand(c, TRUE);
		gtk_column_view_append_column(pm_preset_view, c);
		g_object_unref(c);
		gtk_scrolled_window_set_child(sw_pr, GTK_WIDGET(pm_preset_view));

		/* GTK4: the .ui tooltip says "press the delete key to remove items".
		 * Wire that via GtkEventControllerKey on the preset column view. */
		GtkEventController *pres_key = gtk_event_controller_key_new();
		g_signal_connect(pres_key, "key-pressed",
		                 G_CALLBACK(on_pm_preset_key_pressed), NULL);
		gtk_widget_add_controller(GTK_WIDGET(pm_preset_view), pres_key);

		/* ---- Phase 11 list 3 (operators): build GtkColumnView in C ---- */
		GtkScrolledWindow *sw_op = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_operators"));
		pm_op_store = gtk_string_list_new(NULL);
		GtkSingleSelection *op_sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(pm_op_store)));
		gtk_single_selection_set_can_unselect(op_sel, TRUE);
		gtk_single_selection_set_autoselect(op_sel, FALSE);
		pm_op_view = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(op_sel)));
		g_signal_connect(pm_op_view, "activate", G_CALLBACK(on_pm_op_activate), NULL);

		GtkSignalListItemFactory *facop = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
		g_signal_connect(facop, "setup", G_CALLBACK(pm_string_setup_cb), NULL);
		g_signal_connect(facop, "bind",  G_CALLBACK(pm_op_bind_cb),      NULL);
		GtkColumnViewColumn *cop = gtk_column_view_column_new(N_("Operator"), GTK_LIST_ITEM_FACTORY(facop));
		gtk_column_view_column_set_expand(cop, TRUE);
		gtk_column_view_append_column(pm_op_view, cop);
		g_object_unref(cop);
		gtk_scrolled_window_set_child(sw_op, GTK_WIDGET(pm_op_view));

		/* ---- Phase 11 list 2 (functions): build GtkColumnView in C ---- */
		GtkScrolledWindow *sw_fn = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_functions"));
		pm_fn_store = gtk_string_list_new(NULL);
		GtkSingleSelection *fn_sel = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(pm_fn_store)));
		gtk_single_selection_set_can_unselect(fn_sel, TRUE);
		gtk_single_selection_set_autoselect(fn_sel, FALSE);
		pm_fn_view = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(fn_sel)));
		g_signal_connect(pm_fn_view, "activate", G_CALLBACK(on_pm_fn_activate), NULL);

		GtkSignalListItemFactory *facfn = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
		g_signal_connect(facfn, "setup", G_CALLBACK(pm_string_setup_cb), NULL);
		g_signal_connect(facfn, "bind",  G_CALLBACK(pm_fn_bind_cb),      NULL);
		GtkColumnViewColumn *cfn = gtk_column_view_column_new(N_("Function"), GTK_LIST_ITEM_FACTORY(facfn));
		gtk_column_view_column_set_expand(cfn, TRUE);
		gtk_column_view_append_column(pm_fn_view, cfn);
		g_object_unref(cfn);
		gtk_scrolled_window_set_child(sw_fn, GTK_WIDGET(pm_fn_view));
	}
	g_assert(pm_var_store);
	g_assert(pm_preset_store);
	g_assert(pm_op_store);
	g_assert(pm_fn_store);
	g_assert(pixel_math_status_bar);
	g_assert(pixel_math_status_bar2);
	g_assert(pixel_math_entry_r);
	g_assert(pixel_math_entry_g);
	g_assert(pixel_math_entry_b);
}

/* ── Function table helpers ─────────────────────────────────────────────── */

static int FnCompare_functions(const void *v1, const void *v2) {
	const _pm_op_func *i1 = v1;
	const _pm_op_func *i2 = v2;
	return (g_strcmp0(i1->name, i2->name));
}

static _pm_op_func *concat_functions() {
	_pm_op_func *all_functions = malloc((MAX_FUNCTIONS + MAX_IMAGE_FUNCTIONS) * sizeof(_pm_op_func));
	memcpy(all_functions, functions, MAX_FUNCTIONS * sizeof(_pm_op_func));
	memcpy(all_functions + MAX_FUNCTIONS, image_functions, MAX_IMAGE_FUNCTIONS * sizeof(_pm_op_func));
	qsort(all_functions, MAX_FUNCTIONS + MAX_IMAGE_FUNCTIONS, sizeof(_pm_op_func), FnCompare_functions);
	return all_functions;
}

/* ── Expression entry getters ───────────────────────────────────────────── */

static gchar* get_pixel_math_expression1() {
	init_widgets();
	return g_strdup(gtk_editable_get_text(GTK_EDITABLE(pixel_math_entry_r)));
}

static gchar* get_pixel_math_expression2() {
	init_widgets();
	return g_strdup(gtk_editable_get_text(GTK_EDITABLE(pixel_math_entry_g)));
}

static gchar* get_pixel_math_expression3() {
	init_widgets();
	return g_strdup(gtk_editable_get_text(GTK_EDITABLE(pixel_math_entry_b)));
}

/* ── Status bar helpers ─────────────────────────────────────────────────── */

static void output_status_bar2(int width, int height, int channel) {
	if (width != 0 && height != 0) {
		gchar *out = ngettext("%d x %d pixels (%d channel)", "%d x %d pixels (%d channels)", channel);
		out = g_strdup_printf(out, width, height, channel);
		gtk_label_set_text(pixel_math_status_bar2, out);
		g_free(out);
	} else {
		gtk_label_set_text(pixel_math_status_bar2, "");
	}
}

/* ── Variables-list dimension check (Phase 11) ─────────────────────────── */

static gboolean check_files_dimensions(guint *width, guint *height, guint *channel) {
	init_widgets();
	*channel = 0;
	*width = 0;
	*height = 0;
	guint n = pm_var_store ? g_list_model_get_n_items(G_LIST_MODEL(pm_var_store)) : 0;
	for (guint i = 0; i < n; i++) {
		SirilPmVarRow *row = SIRIL_PM_VAR_ROW(g_list_model_get_item(G_LIST_MODEL(pm_var_store), i));
		if (!row) continue;
		guint c = row->chan, w = row->width, h = row->height;
		g_object_unref(row);
		if (*width == 0) {
			*channel = c;
			*width = w;
			*height = h;
		} else if (w != *width || h != *height || c != *channel) {
			return FALSE;
		}
	}
	return TRUE;
}

/* ── Variables-list accessors (Phase 11) ───────────────────────────────── */

static const gchar *get_pixel_math_var_paths(int i) {
	init_widgets();
	if (!pm_var_store) return NULL;
	if (i < 0 || i >= (int)g_list_model_get_n_items(G_LIST_MODEL(pm_var_store))) return NULL;
	SirilPmVarRow *row = SIRIL_PM_VAR_ROW(g_list_model_get_item(G_LIST_MODEL(pm_var_store), i));
	const gchar *p = row ? row->path : NULL;
	if (row) g_object_unref(row);
	return p;
}

static int get_pixel_math_number_of_rows() {
	return pm_var_store ? (int)g_list_model_get_n_items(G_LIST_MODEL(pm_var_store)) : 0;
}

static int get_pixel_math_functions_number_of_rows() {
	if (pm_fn_store)
		return (int)g_list_model_get_n_items(G_LIST_MODEL(pm_fn_store));
	return 0;
}

static int get_pixel_math_operators_number_of_rows() {
	if (pm_op_store)
		return (int)g_list_model_get_n_items(G_LIST_MODEL(pm_op_store));
	return 0;
}

static int get_pixel_math_presets_number_of_rows() {
	if (pm_preset_store)
		return (int)g_list_model_get_n_items(G_LIST_MODEL(pm_preset_store));
	return 0;
}

static const gchar *get_pixel_math_var_name(int i) {
	init_widgets();
	if (!pm_var_store) return NULL;
	if (i < 0 || i >= (int)g_list_model_get_n_items(G_LIST_MODEL(pm_var_store))) return NULL;
	SirilPmVarRow *row = SIRIL_PM_VAR_ROW(g_list_model_get_item(G_LIST_MODEL(pm_var_store), i));
	const gchar *v = row ? row->var : NULL;
	if (row) g_object_unref(row);
	return v;
}

static const gchar *get_function_name(int i) {
	init_widgets();
	if (!pm_fn_store) return NULL;
	if (i < 0 || i >= (int)g_list_model_get_n_items(G_LIST_MODEL(pm_fn_store)))
		return NULL;
	return gtk_string_list_get_string(pm_fn_store, (guint)i);
}

static const gchar *get_operator_name(int i) {
	init_widgets();
	if (!pm_op_store) return NULL;
	if (i < 0 || i >= (int)g_list_model_get_n_items(G_LIST_MODEL(pm_op_store)))
		return NULL;
	return gtk_string_list_get_string(pm_op_store, (guint)i);
}

static const gchar *get_preset_expr(int i) {
	init_widgets();
	if (!pm_preset_store) return NULL;
	if (i < 0 || i >= (int)g_list_model_get_n_items(G_LIST_MODEL(pm_preset_store)))
		return NULL;
	return gtk_string_list_get_string(pm_preset_store, (guint)i);
}

/* ── Button/checkbox state readers ─────────────────────────────────────── */

static gboolean is_pm_use_rgb_button_checked() {
	return siril_toggle_get_active(GTK_WIDGET(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pm_use_rgb_button"))));
}

static gboolean is_pm_rescale_checked() {
	return siril_toggle_get_active(GTK_WIDGET(GTK_WIDGET(gtk_builder_get_object(gui.builder, "rescale_pm_button"))));
}

static gboolean is_cumulate_checked() {
	return siril_toggle_get_active(GTK_WIDGET(GTK_WIDGET(gtk_builder_get_object(gui.builder, "cumulate_pm_button"))));
}

static float get_min_rescale_value() {
	return (float)gtk_spin_button_get_value(GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_pm_low"))));
}

static float get_max_rescale_value() {
	return (float)gtk_spin_button_get_value(GTK_SPIN_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_pm_high"))));
}

/* ── Parameter substitution ─────────────────────────────────────────────── */

static gboolean is_op_or_null(const gchar c) {
	if (c == '\0') return TRUE;
	if (c == '(') return TRUE;
	if (c == ')') return TRUE;
	if (c == ',') return TRUE;
	for (int i = 0; i < (int)MAX_OPERATORS; i++) {
		const gchar op = operators[i].name[0];
		if (c == op) return TRUE;
	}
	return FALSE;
}

static guint siril_string_replace_parameter(GString *string, const gchar *find,
		const gchar *replace) {
	gsize f_len, r_len, pos;
	gchar *cur, *next;
	gint n = 0;

	g_return_val_if_fail(string != NULL, 0);
	g_return_val_if_fail(find != NULL, 0);
	g_return_val_if_fail(replace != NULL, 0);

	f_len = strlen(find);
	r_len = strlen(replace);
	cur = string->str;

	while ((next = strstr(cur, find)) != NULL) {
		pos = next - string->str;
		if ((pos == 0 && is_op_or_null(string->str[pos + f_len])) ||
				((pos > 0 && is_op_or_null(string->str[pos - 1])) &&
				(pos < strlen(string->str) && is_op_or_null(string->str[pos + f_len])))) {
			g_string_erase(string, pos, f_len);
			g_string_insert(string, pos, replace);
			cur = string->str + pos + r_len;
			n++;
			if (f_len == 0) {
				if (cur[0] == '\0')
					break;
				else
					cur++;
			}
		} else {
			cur = string->str + pos + f_len;
		}
	}

	return n;
}

static int parse_parameters(gchar **expression1, gchar **expression2, gchar **expression3) {
	GtkEntry *entry_param = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_param")));
	char *entry_text = g_strdup(gtk_editable_get_text(GTK_EDITABLE(entry_param)));

	if (entry_text == NULL) return 0;
	remove_spaces_from_str(entry_text);

	gchar **token = g_strsplit(entry_text, ",", -1);
	int nargs = g_strv_length(token);

	for (int i = 0; i < nargs; i++) {
		gchar **expr = g_strsplit(token[i], "=", -1);
		int n = g_strv_length(expr);
		if (n != 2) {
			g_strfreev(token);
			g_strfreev(expr);
			g_free(entry_text);
			return -1;
		}

		GString *string1 = g_string_new(*expression1);
		GString *string2 = g_string_new(*expression2);
		GString *string3 = g_string_new(*expression3);

		siril_string_replace_parameter(string1, expr[0], expr[1]);
		siril_string_replace_parameter(string2, expr[0], expr[1]);
		siril_string_replace_parameter(string3, expr[0], expr[1]);

		g_free(*expression1);
		g_free(*expression2);
		g_free(*expression3);

		*expression1 = g_string_free(string1, FALSE);
		*expression2 = g_string_free(string2, FALSE);
		*expression3 = g_string_free(string3, FALSE);

		g_strfreev(expr);
	}

	g_strfreev(token);
	g_free(entry_text);
	return 0;
}

/* ── Main evaluation entry point ────────────────────────────────────────── */

static int pixel_math_evaluate(gchar *expression1, gchar *expression2, gchar *expression3) {
	int nb_rows = 0;
	int width = -1, height = -1, channel = -1;
	int retval = 0;
	gboolean icc_warning_given = FALSE;
	gboolean single_rgb = is_pm_use_rgb_button_checked();
	gboolean rescale = is_pm_rescale_checked();
	gboolean do_sum = is_cumulate_checked();
	float min_val = get_min_rescale_value();
	float max_val = get_max_rescale_value();

	int max_images = pm_get_max_images();

	while (nb_rows < get_pixel_math_number_of_rows() && nb_rows < max_images) {
		const gchar *path = get_pixel_math_var_paths(nb_rows);
		if (readfits(path, &var_fit[nb_rows], NULL, TRUE)) {
			retval = 1;
			goto free_expressions;
		}

		if (nb_rows > 0) {
			if (!profiles_identical(var_fit[nb_rows].icc_profile, var_fit[0].icc_profile)) {
				if (!icc_warning_given) {
					siril_log_warning(_("ICC profiles are inconsistent. The output color profile will be based on the first layer to be loaded.\n"));
					icc_warning_given = TRUE;
				}
				if (var_fit[0].icc_profile)
					siril_log_warning(_("ICC profile of layer %d does not match the first image. Converting it to match.\n"), nb_rows + 1);
				else
					siril_log_warning(_("The first layer loaded had no color profile. All input layers will be treated as raw data.\n"));
				siril_colorspace_transform(&var_fit[nb_rows], var_fit[0].icc_profile);
			}
		}
		if (channel == -1) {
			width   = var_fit[nb_rows].rx;
			height  = var_fit[nb_rows].ry;
			channel = var_fit[nb_rows].naxes[2];
		}

		if (channel == 3 && !single_rgb) {
			gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Incompatible parameters"),
					_("3 channel images are incompatible with the \"Use single RGB/K expression\" unchecked."));
			retval = 1;
			goto free_expressions;
		}
		nb_rows++;
	}

	channel = single_rgb ? channel : 3;

	struct pixel_math_data *args = calloc(1, sizeof(struct pixel_math_data));

	if (parse_parameters(&expression1, &expression2, &expression3)) {
		gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Parameter error"), _("Parameter symbols could not be parsed."));
		retval = 1;
		goto free_expressions;
	}

	args->expression1 = g_strdup(expression1);
	if (!single_rgb) {
		args->expression2 = g_strdup(expression2);
		args->expression3 = g_strdup(expression3);
	}

	args->single_rgb = single_rgb;
	args->rescale    = rescale;
	args->do_sum     = do_sum;
	args->min        = min_val;
	args->max        = max_val;
	args->nb_rows    = nb_rows;
	args->ret        = 0;
	args->from_ui    = TRUE;
	args->has_gfit   = FALSE;

	args->varname = malloc(nb_rows * sizeof(gchar *));
	for (int i = 0; i < nb_rows; i++)
		args->varname[i] = g_strdup(get_pixel_math_var_name(i));

	if (replace_t_with_gfit(args) && gfit->naxes[2] > 1) {
		gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Incorrect formula"), _("RGB $T cannot be used in this context."));
		retval = 1;
		goto free_expressions;
	}

	fits *fit = NULL;
	if (args->has_gfit) {
		width   = gfit->rx;
		height  = gfit->ry;
		channel = args->single_rgb ? gfit->naxes[2] : 3;
		if (nb_rows > 0) {
			if (width != var_fit[0].naxes[0] || height != var_fit[0].naxes[1]) {
				gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Images have different size"),
					_("The image currently displayed must be the same size as the other images loaded into PixelMath."));
				retval = 1;
				goto free_expressions;
			}
		}
	}

	if (!nb_rows && !args->has_gfit) {
		gui_iface.message_dialog(SIRIL_MSG_ERROR, _("No images loaded"), _("You must load images first."));
		retval = 1;
		goto free_expressions;
	}

	if (new_fit_image(&fit, width, height, channel, com.pref.force_16bit ? DATA_USHORT : DATA_FLOAT)) {
		retval = 1;
		goto free_expressions;
	}

	args->fit = fit;

	if (!start_in_new_thread(apply_pixel_math_operation, args)) {
		g_free(args->expression1);
		g_free(args->expression2);
		g_free(args->expression3);
		free(args);
	}

free_expressions:
	g_free(expression1);
	g_free(expression2);
	g_free(expression3);
	return retval;
}

/* ── Variable validity / list management ────────────────────────────────── */

static gboolean allowed_char(char c) {
	return (g_ascii_isalnum(c) || c == '_');
}

static gboolean check_for_variable_sanity(const gchar *new_text) {
	gboolean abort = TRUE;
	const char *p = new_text;

	if (*p == '\0') return FALSE;
	while (*p) {
		if (!allowed_char(*p++)) return FALSE;
	}
	p = new_text;
	while (*p) {
		if (g_ascii_isalpha(*p++)) {
			abort = FALSE;
		}
	}
	if (abort) return FALSE;

	init_widgets();

	int nb_rows = 0;
	while (nb_rows < get_pixel_math_number_of_rows() && nb_rows < pm_get_max_images()) {
		const gchar *var = get_pixel_math_var_name(nb_rows);
		if (!g_strcmp0(var, new_text)) return FALSE;
		nb_rows++;
	}
	return TRUE;
}

static void add_image_to_variable_list(const gchar *path, const gchar *var, gchar *filter, guint chan, guint width, guint height) {
	const char *name = (filter[0] != '\0' && check_for_variable_sanity(filter)) ? filter : var;
	init_widgets();
	if (!pm_var_store) return;
	SirilPmVarRow *row = g_object_new(SIRIL_TYPE_PM_VAR_ROW, NULL);
	row->var    = g_strdup(name);
	row->path   = g_strdup(path);
	row->chan   = chan;
	row->width  = width;
	row->height = height;
	g_list_store_append(pm_var_store, row);
	g_object_unref(row);
}

static int search_for_free_index() {
	int max = pm_get_max_images();
	int i;
	for (i = 0; i < max; i++) {
		gboolean found = FALSE;
		for (int j = 0; j < max; j++) {
			const gchar *var = get_pixel_math_var_name(j);
			if (!var) break;
			if (!g_strcmp0(var, pm_get_variable_name(i))) {
				found = TRUE;
				break;
			}
		}
		if (!found)
			return i;
	}
	return i + 1;
}

/* Phase 14G.2: GtkFileChooserDialog → GtkFileDialog (multi-select). */
struct pm_select_image_ctx {
	int starting_nb;
};

static void on_pm_files_chosen(GObject *src, GAsyncResult *res, gpointer ud) {
	struct pm_select_image_ctx *ctx = (struct pm_select_image_ctx *)ud;
	int pos = ctx->starting_nb;
	g_free(ctx);

	GtkFileDialog *fd = GTK_FILE_DIALOG(src);
	GError *error = NULL;
	GListModel *files = gtk_file_dialog_open_multiple_finish(fd, res, &error);
	if (!files) {
		if (error) g_clear_error(&error);
		return;
	}
	guint n = g_list_model_get_n_items(files);
	guint width = 0, height = 0, channel = 0;
	for (guint i = 0; i < n; ++i) {
		GFile *gf = G_FILE(g_list_model_get_item(files, i));
		gchar *filename = gf ? g_file_get_path(gf) : NULL;
		if (gf) g_object_unref(gf);
		if (!filename) continue;

		gchar filter[FLEN_VALUE] = { 0 };
		fits f = { 0 };
		if (read_fits_metadata_from_path(filename, &f)) {
			siril_log_error(_("Could not open file: %s\n"), filename);
		} else if (check_files_dimensions(&width, &height, &channel)) {
			if (width != 0 && (channel != f.naxes[2] ||
					width != f.naxes[0] || height != f.naxes[1])) {
				gchar *name = g_path_get_basename(filename);
				gchar *str = g_strdup_printf("%s will not be added in the pixel math tool because its size is different from the other loaded images"
						" (width, height or number of channels).", name);
				gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Image must have same dimension"), str);
				g_free(name);
				g_free(str);
			} else {
				if (f.keywords.filter[0] != '\0')
					memcpy(filter, f.keywords.filter, FLEN_VALUE);
				int idx = search_for_free_index();
				if (idx >= pm_get_max_images()) {
					siril_log_error(_("Error: maximum variable index exceeded - too many variables!\n"));
					g_free(filename);
					clearfits(&f);
					break;
				}
				add_image_to_variable_list(filename, pm_get_variable_name(idx), filter, f.naxes[2], f.naxes[0], f.naxes[1]);
				pos++;
				if (pos == pm_get_max_images()) {
					g_free(filename);
					clearfits(&f);
					break;
				}
			}
		}
		g_free(filename);
		clearfits(&f);
	}
	g_object_unref(files);
}

static void select_image(int nb) {
	GtkFileDialog *fd = gtk_file_dialog_new();
	gtk_file_dialog_set_title(fd, "Open File");
	GFile *initial_folder = g_file_new_for_path(com.wd);
	gtk_file_dialog_set_initial_folder(fd, initial_folder);
	g_object_unref(initial_folder);

	GtkFileFilter *fitsf = gtk_file_filter_new();
	gtk_file_filter_set_name(fitsf, _("FITS Files (*.fit, *.fits, *.fts, *.fit.fz, *.fits.fz, *.fts.fz)"));
	{
		gchar **patterns = g_strsplit(FITS_EXTENSIONS, ";", -1);
		for (int i = 0; patterns[i]; ++i)
			gtk_file_filter_add_pattern(fitsf, patterns[i]);
		g_strfreev(patterns);
	}
	GListStore *fl = g_list_store_new(GTK_TYPE_FILE_FILTER);
	g_list_store_append(fl, fitsf);
	if (gui.file_ext_filter == TYPEFITS)
		gtk_file_dialog_set_default_filter(fd, fitsf);
	gtk_file_dialog_set_filters(fd, G_LIST_MODEL(fl));
	g_object_unref(fitsf);
	g_object_unref(fl);

	struct pm_select_image_ctx *ctx = g_new0(struct pm_select_image_ctx, 1);
	ctx->starting_nb = nb;

	gtk_file_dialog_open_multiple(fd,
		GTK_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_dialog"))),
		NULL, on_pm_files_chosen, ctx);
	g_object_unref(fd);
}

/* ── apply_pixel_math (collects widget state → calls evaluate) ──────────── */

static GtkEntry *get_entry_with_focus() {
	switch(entry_has_focus) {
	default:
	case 0: return pixel_math_entry_r;
	case 1: return pixel_math_entry_g;
	case 2: return pixel_math_entry_b;
	}
}

static void remove_selected_lines() {
	init_widgets();
	if (!pm_var_sel || !pm_var_store) return;
	GtkBitset *bs = gtk_selection_model_get_selection(GTK_SELECTION_MODEL(pm_var_sel));
	guint64 n = gtk_bitset_get_size(bs);
	for (guint64 i = n; i > 0; i--) {
		guint pos = gtk_bitset_get_nth(bs, (guint)(i - 1));
		g_list_store_remove(pm_var_store, pos);
	}
	gtk_bitset_unref(bs);
	gtk_selection_model_unselect_all(GTK_SELECTION_MODEL(pm_var_sel));
}

static void apply_pixel_math() {
	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	gchar *expression1 = get_pixel_math_expression1();
	remove_spaces_from_str(expression1);

	gchar *expression2 = get_pixel_math_expression2();
	remove_spaces_from_str(expression2);

	gchar *expression3 = get_pixel_math_expression3();
	remove_spaces_from_str(expression3);

	if (pixel_math_evaluate(expression1, expression2, expression3))
		siril_log_error(_("Error evaluating pixelmath expression.\n"));
}

/* ── Preset list management ─────────────────────────────────────────────── */

static void add_expr_to_tree(const gchar *expression) {
	init_widgets();
	if (pm_preset_store)
		gtk_string_list_append(pm_preset_store, expression);
}

static void save_presets_list() {
	g_slist_free_full(com.pref.gui.pm_presets, g_free);
	com.pref.gui.pm_presets = NULL;
	guint n = pm_preset_store ? g_list_model_get_n_items(G_LIST_MODEL(pm_preset_store)) : 0;
	for (guint i = 0; i < n; i++) {
		const gchar *expr = gtk_string_list_get_string(pm_preset_store, i);
		if (expr)
			com.pref.gui.pm_presets = g_slist_append(com.pref.gui.pm_presets, g_strdup(expr));
	}
	writeinitfile();
}

static void add_functions_to_list() {
	init_widgets();
	if (!pm_fn_store) return;
	if (!pm_fn_sorted) {
		pm_fn_sorted = concat_functions();
		pm_fn_sorted_count = (int)(MAX_FUNCTIONS + MAX_IMAGE_FUNCTIONS);
	}
	for (int i = 0; i < pm_fn_sorted_count; i++)
		gtk_string_list_append(pm_fn_store, pm_fn_sorted[i].name);
}

static void add_operators_to_list() {
	init_widgets();
	if (!pm_op_store) return;
	for (int i = 0; i < (int)MAX_OPERATORS; i++)
		gtk_string_list_append(pm_op_store, operators[i].name);
}

static void add_presets_to_list() {
	init_widgets();
	if (!pm_preset_store) return;
	for (GSList *l = com.pref.gui.pm_presets; l; l = l->next) {
		gtk_string_list_append(pm_preset_store, (const gchar *)l->data);
	}
}

/* Delete-key handler for the preset GtkColumnView. */
static gboolean on_pm_preset_key_pressed(GtkEventControllerKey *controller,
		guint keyval, guint keycode, GdkModifierType state, gpointer user_data) {
	(void)controller; (void)keycode; (void)state; (void)user_data;
	if (keyval == GDK_KEY_Delete || keyval == GDK_KEY_KP_Delete
			|| keyval == GDK_KEY_BackSpace) {
		remove_presets_from_list();
		return TRUE;
	}
	return FALSE;
}

static void remove_presets_from_list() {
	init_widgets();
	if (!pm_preset_sel || !pm_preset_store) return;
	GtkBitset *bs = gtk_selection_model_get_selection(GTK_SELECTION_MODEL(pm_preset_sel));
	guint64 n = gtk_bitset_get_size(bs);
	for (guint64 i = n; i > 0; i--) {
		guint pos = gtk_bitset_get_nth(bs, (guint)(i - 1));
		gtk_string_list_remove(pm_preset_store, pos);
	}
	gtk_bitset_unref(bs);
	gtk_selection_model_unselect_all(GTK_SELECTION_MODEL(pm_preset_sel));
	save_presets_list();
}

/* ── GTK signal callbacks ───────────────────────────────────────────────── */

void on_pm_use_rgb_button_toggled(GtkCheckButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_g")), !siril_toggle_get_active(GTK_WIDGET(button)));
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_b")), !siril_toggle_get_active(GTK_WIDGET(button)));
	gtk_label_set_text(GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_RGBK"))), siril_toggle_get_active(GTK_WIDGET(button)) ? _("RGB/K") : _("R"));
	if (siril_toggle_get_active(GTK_WIDGET(button))) entry_has_focus = 0;
}

gboolean on_pixel_math_entry_r_focus_in_event(GtkWidget *widget,
		GdkEvent *event, gpointer user_data) {
	entry_has_focus = 0;
	return FALSE;
}

gboolean on_pixel_math_entry_g_focus_in_event(GtkWidget *widget,
		GdkEvent *event, gpointer user_data) {
	entry_has_focus = 1;
	return FALSE;
}

gboolean on_pixel_math_entry_b_focus_in_event(GtkWidget *widget,
		GdkEvent *event, gpointer user_data) {
	entry_has_focus = 2;
	return FALSE;
}

void on_pixel_math_add_var_button_clicked(GtkButton *button, gpointer user_data) {
	init_widgets();
	int nb = get_pixel_math_number_of_rows();
	if (nb == pm_get_max_images()) {
		gui_iface.message_dialog(SIRIL_MSG_WARNING, _("Cannot load new image"),
				_("You've reached the maximum of loaded image file."));
	} else {
		select_image(nb);
	}
}

void on_pixel_math_remove_var_button_clicked(GtkButton *button, gpointer user_data) {
	remove_selected_lines();
}

void on_pixel_math_entry_activate(GtkEntry *entry, gpointer user_data) {
	apply_pixel_math();
}

void on_apply_pixel_math_clicked(GtkButton *button, gpointer user_data) {
	apply_pixel_math();
}

/* Phase 11 list 1: superseded by on_pm_var_activate (GtkColumnView::activate).
 * Stub kept for any stale .ui binding. */
void on_pixel_math_treeview_row_activated(GtkWidget *tree_view,
		gpointer path, gpointer column) {
	(void)tree_view; (void)path; (void)column;
}

/* Phase 11 list 2: tooltips are now set per-row inside pm_fn_bind_cb;
 * this old query-tooltip handler is no longer wired. */
gboolean query_tooltip_tree_view_cb(GtkWidget *widget, gint x, gint y,
		gboolean keyboard_tip, GtkTooltip *tooltip, gpointer data) {
	(void)widget; (void)x; (void)y; (void)keyboard_tip; (void)tooltip; (void)data;
	return FALSE;
}

/* Phase 11 list 3: tooltips are now set per-row inside pm_op_bind_cb;
 * this old query-tooltip handler is no longer wired. */
gboolean query_tooltip_op_tree_view_cb(GtkWidget *widget, gint x, gint y,
		gboolean keyboard_tip, GtkTooltip *tooltip, gpointer data) {
	(void)widget; (void)x; (void)y; (void)keyboard_tip; (void)tooltip; (void)data;
	return FALSE;
}

void on_pixel_math_dialog_show(GtkWidget *w, gpointer user_data) {
	if (!get_pixel_math_functions_number_of_rows())
		add_functions_to_list();
	if (!get_pixel_math_operators_number_of_rows())
		add_operators_to_list();
	if (!get_pixel_math_presets_number_of_rows())
		add_presets_to_list();
}

/* Phase 11 list 2: superseded by on_pm_fn_activate (GtkColumnView::activate).
 * Stub kept for any stale .ui binding. */
void on_pixel_math_treeview_functions_row_activated(GtkWidget *tree_view,
		gpointer path, gpointer column) {
	(void)tree_view; (void)path; (void)column;
}

/* Phase 11 list 3: superseded by on_pm_op_activate (GtkColumnView::activate).
 * Stub kept for any stale .ui binding. */
void on_pixel_math_treeview_operators_row_activated(GtkWidget *tree_view,
		gpointer path, gpointer column) {
	(void)tree_view; (void)path; (void)column;
}

void on_close_pixel_math_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("pixel_math_dialog");
}

/* TODO Phase8: convert to GtkEventController */
gboolean on_pixel_math_treeview_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	(void)widget; (void)event; (void)user_data;
	return FALSE;
}

/* Phase 11 list 1: editing now happens inline via GtkEditableLabel inside
 * pm_var_setup_var_cb / pm_var_bind_var_cb / on_pm_var_editing_changed.
 * Selection-change is wired to on_pm_var_selection_changed.  These three
 * old GtkCellRendererText handlers are no longer wired by .ui, but kept
 * as no-op stubs to keep the symbols available for stale builder caches. */
void on_cellrenderer_variables_edited(GtkWidget *renderer, char *path,
		char *new_text, gpointer user_data) {
	(void)renderer; (void)path; (void)new_text; (void)user_data;
}
void on_cellrenderer_variables_editing_started(GtkWidget *renderer,
		gpointer editable, char *path, gpointer user_data) {
	(void)renderer; (void)editable; (void)path; (void)user_data;
}
void on_cellrenderer_variables_editing_canceled(GtkWidget *renderer,
		gpointer user_data) {
	(void)renderer; (void)user_data;
}
void on_pixel_math_selection_changed(GObject *selection, gpointer user_data) {
	(void)selection; (void)user_data;
}

void on_pm_expr1_bt_clicked(GtkButton *button, gpointer user_data) {
	gchar *str = get_pixel_math_expression1();
	if (str) {
		g_strstrip(str);
		if (str[0] != '\0') {
			add_expr_to_tree(str);
			save_presets_list();
		}
		g_free(str);
	}
}

void on_pm_expr2_bt_clicked(GtkButton *button, gpointer user_data) {
	gchar *str = get_pixel_math_expression2();
	if (str) {
		g_strstrip(str);
		if (str[0] != '\0') {
			add_expr_to_tree(str);
			save_presets_list();
		}
		g_free(str);
	}
}

void on_pm_expr3_bt_clicked(GtkButton *button, gpointer user_data) {
	gchar *str = get_pixel_math_expression3();
	if (str) {
		g_strstrip(str);
		if (str[0] != '\0') {
			add_expr_to_tree(str);
			save_presets_list();
		}
		g_free(str);
	}
}

/* Phase 11 list 4: superseded by on_pm_preset_activate (GtkColumnView::activate).
 * Stub kept for any stale .ui binding. */
void on_pixel_math_treeview_presets_row_activated(GtkWidget *tree_view,
		gpointer path, gpointer column) {
	(void)tree_view; (void)path; (void)column;
}

/* TODO Phase8: convert to GtkEventController */
gboolean on_pixel_math_treeview_presets_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	(void)widget; (void)event; (void)user_data;
	return FALSE;
}

void on_rescale_pm_button_toggled(GtkCheckButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_pm_low")),  siril_toggle_get_active(GTK_WIDGET(button)));
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_pm_high")), siril_toggle_get_active(GTK_WIDGET(button)));
}
