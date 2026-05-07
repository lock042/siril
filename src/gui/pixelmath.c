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

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "core/initfile.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/dialog_preview.h"
#include "gui/open_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/gui_state.h"
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

static GtkListStore *pixel_math_list_store = NULL;
static GtkListStore *pixel_math_list_store_functions = NULL;
static GtkListStore *pixel_math_list_store_operators = NULL;
static GtkListStore *pixel_math_list_store_presets = NULL;
static GtkTreeView  *pixel_math_tree_view = NULL;
static GtkTreeView  *pixel_math_treeview_functions = NULL;
static GtkTreeView  *pixel_math_treeview_operators = NULL;
static GtkTreeView  *pixel_math_treeview_presets = NULL;
static GtkTreeModel *pixel_math_tree_model = NULL;
static GtkTreeModel *pixel_math_tree_model_functions = NULL;
static GtkTreeModel *pixel_math_tree_model_operators = NULL;
static GtkTreeModel *pixel_math_tree_model_presets = NULL;
static GtkLabel     *pixel_math_status_bar = NULL;
static GtkLabel     *pixel_math_status_bar2 = NULL;
static GtkEntry     *pixel_math_entry_r = NULL;
static GtkEntry     *pixel_math_entry_g = NULL;
static GtkEntry     *pixel_math_entry_b = NULL;

static int entry_has_focus = 0;

/* ── Widget initialisation ──────────────────────────────────────────────── */

static void init_widgets() {
	if (!pixel_math_tree_view) {
		pixel_math_tree_view = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview"));
		pixel_math_tree_model = gtk_tree_view_get_model(pixel_math_tree_view);
		pixel_math_list_store = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore"));
		pixel_math_status_bar = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_status")));
		pixel_math_status_bar2 = GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_status2")));
		pixel_math_entry_r = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_r")));
		pixel_math_entry_g = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_g")));
		pixel_math_entry_b = GTK_ENTRY(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_b")));
		pixel_math_treeview_functions = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview_functions"));
		pixel_math_treeview_operators = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview_operators"));
		pixel_math_treeview_presets = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview_presets"));
		pixel_math_tree_model_functions = gtk_tree_view_get_model(pixel_math_treeview_functions);
		pixel_math_tree_model_operators = gtk_tree_view_get_model(pixel_math_treeview_operators);
		pixel_math_tree_model_presets = gtk_tree_view_get_model(pixel_math_treeview_presets);
		pixel_math_list_store_functions = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_functions"));
		pixel_math_list_store_operators = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_operators"));
		pixel_math_list_store_presets = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_presets"));

#if GTK_CHECK_VERSION(3, 22, 0)
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_functions"))), TRUE);
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_scrolled_operators"))), TRUE);
#endif

		/* Due to a glade bug, this property is often removed, so we code it */
		GtkTreeSelection *selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "pixel_math_selection"));
		gtk_tree_selection_set_mode(selection, GTK_SELECTION_MULTIPLE);
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "pixel_math_presets_selection"));
		gtk_tree_selection_set_mode(selection, GTK_SELECTION_MULTIPLE);
	}
	g_assert(pixel_math_tree_view);
	g_assert(pixel_math_tree_model);
	g_assert(pixel_math_list_store);
	g_assert(pixel_math_list_store_functions);
	g_assert(pixel_math_list_store_operators);
	g_assert(pixel_math_list_store_presets);
	g_assert(pixel_math_status_bar);
	g_assert(pixel_math_status_bar2);
	g_assert(pixel_math_entry_r);
	g_assert(pixel_math_entry_g);
	g_assert(pixel_math_entry_b);
	g_assert(pixel_math_treeview_functions);
	g_assert(pixel_math_tree_model_functions);
	g_assert(pixel_math_treeview_operators);
	g_assert(pixel_math_tree_model_operators);
	g_assert(pixel_math_treeview_presets);
	g_assert(pixel_math_tree_model_presets);
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
	return g_strdup(gtk_entry_get_text(pixel_math_entry_r));
}

static gchar* get_pixel_math_expression2() {
	init_widgets();
	return g_strdup(gtk_entry_get_text(pixel_math_entry_g));
}

static gchar* get_pixel_math_expression3() {
	init_widgets();
	return g_strdup(gtk_entry_get_text(pixel_math_entry_b));
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

/* ── Tree-view dimension check ──────────────────────────────────────────── */

static gboolean check_files_dimensions(guint *width, guint *height, guint *channel) {
	init_widgets();
	GtkTreeIter iter;
	*channel = 0;
	*width = 0;
	*height = 0;
	gboolean valid = gtk_tree_model_get_iter_first(pixel_math_tree_model, &iter);

	while (valid) {
		guint c, w, h;
		gtk_tree_model_get(pixel_math_tree_model, &iter, COLUMN_IMAGE_CHAN, &c, COLUMN_IMAGE_WIDTH, &w, COLUMN_IMAGE_HEIGHT, &h, -1);
		if (*width == 0) {
			*channel = c;
			*width = w;
			*height = h;
		} else {
			if (w != *width || h != *height || c != *channel) {
				return FALSE;
			}
		}
		valid = gtk_tree_model_iter_next(pixel_math_tree_model, &iter);
	}
	return TRUE;
}

/* ── Tree-view accessors ────────────────────────────────────────────────── */

static const gchar *get_pixel_math_var_paths(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;
	init_widgets();
	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path);
	gtk_tree_model_get_value(pixel_math_tree_model, &iter, COLUMN_IMAGE_PATH, &value);
	return g_value_get_string(&value);
}

static int get_pixel_math_number_of_rows() {
	if (GTK_IS_TREE_MODEL(pixel_math_list_store))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store), NULL);
	else return 0;
}

static int get_pixel_math_functions_number_of_rows() {
	if (GTK_IS_TREE_MODEL(pixel_math_list_store_functions))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store_functions), NULL);
	else return 0;
}

static int get_pixel_math_operators_number_of_rows() {
	if (GTK_IS_TREE_MODEL(pixel_math_list_store_operators))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store_operators), NULL);
	else return 0;
}

static int get_pixel_math_presets_number_of_rows() {
	if (GTK_IS_TREE_MODEL(pixel_math_list_store_presets))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store_presets), NULL);
	else return 0;
}

static const gchar *get_pixel_math_var_name(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;
	init_widgets();
	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	if (gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path)) {
		gtk_tree_model_get_value(pixel_math_tree_model, &iter, COLUMN_IMAGE_NUM, &value);
		return g_value_get_string(&value);
	}
	return NULL;
}

static const gchar *get_function_name(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;
	init_widgets();
	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	if (gtk_tree_model_get_iter(pixel_math_tree_model_functions, &iter, path)) {
		gtk_tree_model_get_value(pixel_math_tree_model_functions, &iter, COLUMN_NAME, &value);
		return g_value_get_string(&value);
	}
	return NULL;
}

static const gchar *get_operator_name(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;
	init_widgets();
	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	if (gtk_tree_model_get_iter(pixel_math_tree_model_operators, &iter, path)) {
		gtk_tree_model_get_value(pixel_math_tree_model_operators, &iter, COLUMN_NAME, &value);
		return g_value_get_string(&value);
	}
	return NULL;
}

static const gchar *get_preset_expr(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;
	init_widgets();
	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	if (gtk_tree_model_get_iter(pixel_math_tree_model_presets, &iter, path)) {
		gtk_tree_model_get_value(pixel_math_tree_model_presets, &iter, 0, &value);
		return g_value_get_string(&value);
	}
	return NULL;
}

/* ── Button/checkbox state readers ─────────────────────────────────────── */

static gboolean is_pm_use_rgb_button_checked() {
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pm_use_rgb_button"))));
}

static gboolean is_pm_rescale_checked() {
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "rescale_pm_button"))));
}

static gboolean is_cumulate_checked() {
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "cumulate_pm_button"))));
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
	char *entry_text = g_strdup(gtk_entry_get_text(entry_param));

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
					siril_log_color_message(_("ICC profiles are inconsistent. The output color profile will be based on the first layer to be loaded.\n"), "salmon");
					icc_warning_given = TRUE;
				}
				if (var_fit[0].icc_profile)
					siril_log_color_message(_("ICC profile of layer %d does not match the first image. Converting it to match.\n"), "salmon", nb_rows + 1);
				else
					siril_log_color_message(_("The first layer loaded had no color profile. All input layers will be treated as raw data.\n"), "salmon");
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

static gboolean check_for_variable_sanity(char *new_text) {
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
	GtkTreeIter iter;
	const char *name;
	if (filter[0] != '\0' && check_for_variable_sanity(filter)) {
		name = filter;
	} else {
		name = var;
	}
	gtk_list_store_append(pixel_math_list_store, &iter);
	gtk_list_store_set(pixel_math_list_store, &iter,
			COLUMN_IMAGE_NUM, name,
			COLUMN_IMAGE_PATH, path,
			COLUMN_IMAGE_CHAN, chan,
			COLUMN_IMAGE_WIDTH, width,
			COLUMN_IMAGE_HEIGHT, height,
			-1);
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

static void select_image(int nb) {
	GtkWidget *dialog;
	fileChooserPreview *preview = NULL;
	GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
	gint res;

	dialog = gtk_file_chooser_dialog_new("Open File",
			GTK_WINDOW(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_dialog"))),
			action,
			_("_Cancel"), GTK_RESPONSE_CANCEL,
			_("_Open"),   GTK_RESPONSE_ACCEPT,
			NULL);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), com.wd);
	gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(dialog), TRUE);
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(dialog), FALSE);
	gtk_filter_add(GTK_FILE_CHOOSER(dialog), _("FITS Files (*.fit, *.fits, *.fts, *.fit.fz, *.fits.fz, *.fts.fz)"), FITS_EXTENSIONS, gui.file_ext_filter == TYPEFITS);
	siril_file_chooser_add_preview(GTK_FILE_CHOOSER(dialog), preview);

	res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		GSList *l;
		int pos = nb;
		guint width = 0, height = 0, channel = 0;
		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		GSList *filenames = siril_file_chooser_get_filenames(chooser);

		for (l = filenames; l; l = l->next) {
			char *filename = (char *)l->data;
			if (filename) {
				gchar filter[FLEN_VALUE] = { 0 };
				fits f = { 0 };
				if (read_fits_metadata_from_path(filename, &f)) {
					siril_log_color_message(_("Could not open file: %s\n"), "red", filename);
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
							siril_log_color_message(_("Error: maximum variable index exceeded - too many variables!\n"), "red");
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
		}
	}

	siril_preview_free(preview);
	gtk_widget_destroy(dialog);
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
	GtkTreeSelection *selection;
	GList *references, *list;

	init_widgets();
	selection = gtk_tree_view_get_selection(pixel_math_tree_view);
	references = get_row_references_of_selected_rows(selection, pixel_math_tree_model);
	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			if (gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path))
				gtk_list_store_remove(pixel_math_list_store, &iter);
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
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
		siril_log_color_message(_("Error evaluating pixelmath expression.\n"), "red");
}

/* ── Preset list management ─────────────────────────────────────────────── */

static void add_expr_to_tree(const gchar *expression) {
	GtkTreeIter iter;
	gtk_list_store_append(pixel_math_list_store_presets, &iter);
	gtk_list_store_set(pixel_math_list_store_presets, &iter, 0, expression, -1);
}

static gboolean foreach_func(GtkTreeModel *model, GtkTreePath *path,
		GtkTreeIter *iter, gpointer user_data) {
	gchar *expression;
	gtk_tree_model_get(model, iter, 0, &expression, -1);
	com.pref.gui.pm_presets = g_slist_prepend(com.pref.gui.pm_presets, expression);
	return FALSE;
}

static void save_presets_list() {
	g_slist_free_full(com.pref.gui.pm_presets, g_free);
	com.pref.gui.pm_presets = NULL;
	gtk_tree_model_foreach(GTK_TREE_MODEL(pixel_math_list_store_presets), foreach_func, com.pref.gui.pm_presets);
	com.pref.gui.pm_presets = g_slist_reverse(com.pref.gui.pm_presets);
	writeinitfile();
}

static void add_functions_to_list() {
	GtkTreeIter iter;
	_pm_op_func *all_functions = concat_functions();
	init_widgets();
	for (int i = 0; i < (int)(MAX_FUNCTIONS + MAX_IMAGE_FUNCTIONS); i++) {
		gtk_list_store_append(pixel_math_list_store_functions, &iter);
		gtk_list_store_set(pixel_math_list_store_functions, &iter, COLUMN_NAME, all_functions[i].name, COLUMN_INDEX, i, -1);
	}
	free(all_functions);
}

static void add_operators_to_list() {
	GtkTreeIter iter;
	init_widgets();
	for (int i = 0; i < (int)MAX_OPERATORS; i++) {
		gtk_list_store_append(pixel_math_list_store_operators, &iter);
		gtk_list_store_set(pixel_math_list_store_operators, &iter, COLUMN_NAME, operators[i].name, COLUMN_INDEX, i, -1);
	}
}

static void add_presets_to_list() {
	GtkTreeIter iter;
	init_widgets();
	com.pref.gui.pm_presets = g_slist_reverse(com.pref.gui.pm_presets);
	for (GSList *l = com.pref.gui.pm_presets; l; l = l->next) {
		gtk_list_store_append(pixel_math_list_store_presets, &iter);
		gtk_list_store_set(pixel_math_list_store_presets, &iter, 0, (gchar *)l->data, -1);
	}
}

static void remove_presets_from_list() {
	GtkTreeSelection *selection;
	GList *references, *list;
	init_widgets();
	selection = gtk_tree_view_get_selection(pixel_math_treeview_presets);
	references = get_row_references_of_selected_rows(selection, pixel_math_tree_model_presets);
	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			if (gtk_tree_model_get_iter(pixel_math_tree_model_presets, &iter, path))
				gtk_list_store_remove(pixel_math_list_store_presets, &iter);
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
	save_presets_list();
}

static gint get_real_index_from_index_in_list(GtkTreeModel *model, GtkTreeIter *iter) {
	gint real_index;
	gtk_tree_model_get(model, iter, COLUMN_INDEX, &real_index, -1);
	return real_index;
}

/* ── GTK signal callbacks ───────────────────────────────────────────────── */

void on_pm_use_rgb_button_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_g")), !gtk_toggle_button_get_active(button));
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_entry_b")), !gtk_toggle_button_get_active(button));
	gtk_label_set_text(GTK_LABEL(GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_RGBK"))), gtk_toggle_button_get_active(button) ? _("RGB/K") : _("R"));
	if (gtk_toggle_button_get_active(button)) entry_has_focus = 0;
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

void on_pixel_math_treeview_row_activated(GtkTreeView *tree_view,
		GtkTreePath *path, GtkTreeViewColumn *column) {
	init_widgets();
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	const gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_pixel_math_var_name(i[0]);

	if (str) {
		guint position = gtk_editable_get_position(GTK_EDITABLE(entry));
		gtk_editable_delete_selection(GTK_EDITABLE(entry));
		gtk_entry_buffer_insert_text(buffer, position, str, -1);
		gtk_widget_grab_focus(GTK_WIDGET(entry));
		gtk_editable_set_position(GTK_EDITABLE(entry), position + strlen(str));
	}
}

gboolean query_tooltip_tree_view_cb(GtkWidget *widget, gint x, gint y,
		gboolean keyboard_tip, GtkTooltip *tooltip, gpointer data) {
	GtkTreeIter iter;
	GtkTreeView *tree_view = GTK_TREE_VIEW(widget);
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GtkTreePath *path = NULL;
	_pm_op_func *all_functions = concat_functions();
	char buffer[512];

	if (!gtk_tree_view_get_tooltip_context(tree_view, &x, &y, keyboard_tip,
			&model, &path, &iter)) {
		free(all_functions);
		return FALSE;
	}

	gint real_index = get_real_index_from_index_in_list(model, &iter);
	g_snprintf(buffer, 511, "<b>%s</b>\n\n%s", all_functions[real_index].prototype, _(all_functions[real_index].definition));
	gtk_tooltip_set_markup(tooltip, buffer);
	gtk_tree_view_set_tooltip_row(tree_view, tooltip, path);
	gtk_tree_path_free(path);
	free(all_functions);
	return TRUE;
}

gboolean query_tooltip_op_tree_view_cb(GtkWidget *widget, gint x, gint y,
		gboolean keyboard_tip, GtkTooltip *tooltip, gpointer data) {
	GtkTreeIter iter;
	GtkTreeView *tree_view = GTK_TREE_VIEW(widget);
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	GtkTreePath *path = NULL;
	char buffer[512];

	if (!gtk_tree_view_get_tooltip_context(tree_view, &x, &y, keyboard_tip,
			&model, &path, &iter))
		return FALSE;

	gint real_index = get_real_index_from_index_in_list(model, &iter);
	g_snprintf(buffer, 511, "<b>%s</b>\n\n%s", operators[real_index].prototype, _(operators[real_index].definition));
	gtk_tooltip_set_markup(tooltip, buffer);
	gtk_tree_view_set_tooltip_row(tree_view, tooltip, path);
	gtk_tree_path_free(path);
	return TRUE;
}

void on_pixel_math_dialog_show(GtkWidget *w, gpointer user_data) {
	if (!get_pixel_math_functions_number_of_rows())
		add_functions_to_list();
	if (!get_pixel_math_operators_number_of_rows())
		add_operators_to_list();
	if (!get_pixel_math_presets_number_of_rows())
		add_presets_to_list();
}

void on_pixel_math_treeview_functions_row_activated(GtkTreeView *tree_view,
		GtkTreePath *path, GtkTreeViewColumn *column) {
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	const gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_function_name(i[0]);

	if (str) {
		guint position = gtk_editable_get_position(GTK_EDITABLE(entry));
		gtk_editable_delete_selection(GTK_EDITABLE(entry));
		gtk_entry_buffer_insert_text(buffer, position, str, -1);
		gtk_widget_grab_focus(GTK_WIDGET(entry));
		gtk_editable_set_position(GTK_EDITABLE(entry), position + strlen(str));
	}
}

void on_pixel_math_treeview_operators_row_activated(GtkTreeView *tree_view,
		GtkTreePath *path, GtkTreeViewColumn *column) {
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	const gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_operator_name(i[0]);

	if (str) {
		guint position = gtk_editable_get_position(GTK_EDITABLE(entry));
		gtk_editable_delete_selection(GTK_EDITABLE(entry));
		gtk_entry_buffer_insert_text(buffer, position, str, -1);
		gtk_widget_grab_focus(GTK_WIDGET(entry));
		gtk_editable_set_position(GTK_EDITABLE(entry), position + strlen(str));
	}
}

void on_close_pixel_math_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("pixel_math_dialog");
}

gboolean on_pixel_math_treeview_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {
		remove_selected_lines();
		return TRUE;
	}
	return FALSE;
}

void on_cellrenderer_variables_edited(GtkCellRendererText *renderer, char *path,
		char *new_text, gpointer user_data) {
	g_signal_handlers_unblock_by_func(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_treeview")), on_pixel_math_treeview_key_release_event, NULL);

	if (check_for_variable_sanity(new_text)) {
		GtkTreeIter iter;
		GValue value = G_VALUE_INIT;
		g_value_init(&value, G_TYPE_STRING);
		g_assert(G_VALUE_HOLDS_STRING(&value));
		gtk_tree_model_get_iter_from_string(pixel_math_tree_model, &iter, path);
		g_value_set_string(&value, new_text);
		gtk_list_store_set_value(pixel_math_list_store, &iter, COLUMN_IMAGE_NUM, &value);
		g_value_unset(&value);
	}
}

void on_cellrenderer_variables_editing_started(GtkCellRenderer *renderer,
		GtkCellEditable *editable, char *path, gpointer user_data) {
	g_signal_handlers_block_by_func(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_treeview")), on_pixel_math_treeview_key_release_event, NULL);
}

void on_cellrenderer_variables_editing_canceled(GtkCellRenderer *renderer,
		gpointer user_data) {
	g_signal_handlers_unblock_by_func(GTK_WIDGET(gtk_builder_get_object(gui.builder, "pixel_math_treeview")), on_pixel_math_treeview_key_release_event, NULL);
}

void on_pixel_math_selection_changed(GtkTreeSelection *selection, gpointer user_data) {
	GList *list;
	guint width = 0, height = 0, channel = 0;

	list = gtk_tree_selection_get_selected_rows(selection, &pixel_math_tree_model);
	if (g_list_length(list) == 1) {
		GtkTreeIter iter;
		gtk_tree_model_get_iter(pixel_math_tree_model, &iter, (GtkTreePath*)list->data);
		gtk_tree_model_get(pixel_math_tree_model, &iter, COLUMN_IMAGE_CHAN, &channel, COLUMN_IMAGE_WIDTH, &width, COLUMN_IMAGE_HEIGHT, &height, -1);
	}
	g_list_free_full(list, (GDestroyNotify)gtk_tree_path_free);
	output_status_bar2(width, height, channel);
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

void on_pixel_math_treeview_presets_row_activated(GtkTreeView *tree_view,
		GtkTreePath *path, GtkTreeViewColumn *column) {
	init_widgets();
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	const gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_preset_expr(i[0]);

	if (str) {
		guint position = gtk_editable_get_position(GTK_EDITABLE(entry));
		gtk_editable_delete_selection(GTK_EDITABLE(entry));
		gtk_entry_buffer_insert_text(buffer, position, str, -1);
		gtk_widget_grab_focus(GTK_WIDGET(entry));
		gtk_editable_set_position(GTK_EDITABLE(entry), position + strlen(str));
	}
}

gboolean on_pixel_math_treeview_presets_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {
		remove_presets_from_list();
		return TRUE;
	}
	return FALSE;
}

void on_rescale_pm_button_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_pm_low")),  gtk_toggle_button_get_active(button));
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "spin_pm_high")), gtk_toggle_button_get_active(button));
}
