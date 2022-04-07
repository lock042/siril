/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/dialog_preview.h"
#include "gui/image_display.h"
#include "gui/message_dialog.h"
#include "gui/open_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/histogram.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"

#include "tinyexpr.h"
#include "pixel_math_runner.h"

static const gchar *variables[] = {
		"I1",
		"I2",
		"I3",
		"I4",
		"I5",
		"I6",
		"I7",
		"I8",
		"I9",
		"I10"
};

static int entry_has_focus = 0;
#define MAX_IMAGES G_N_ELEMENTS(variables)

typedef struct {
	const gchar *name;
	const gchar *prototype;
	const gchar *definition;
} _pm_op_func;

static const _pm_op_func functions[] = {
    { "abs",   "abs ( x )",                             N_("Absolute value of x.")                    },
    { "acos",  "acos ( x )",                            N_("Arc cosine of x.")                        },
    { "asin",  "asin ( x )",                            N_("Arc sine of x.")                          },
    { "atan",  "atan ( x )",                            N_("Arc tangent of x.")                       },
    { "atan2", "atan2 ( y, x )",                        N_("Arc tangent of y/x.")                     },
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
    { "ncr",   "ncr ( x, y )",                          N_("Combinations function.")                  },
    { "npr",   "npr ( x, y )",                          N_("Permutations function.")                  },
    { "pi",    "pi",                                    N_("The constant \u03c0=3.141592...")         },
    { "pow",   "pow ( x, y )",                          N_("Exponentiation function.")                },
    { "sin",   "sin ( x )",                             N_("Sine of x.")                              },
    { "sinh",  "sinh ( x )",                            N_("Hyperbolic sine of x.")                   },
    { "sqrt",  "sqrt ( x )",                            N_("Square root of x.")                       },
    { "tan",   "tan ( x )",                             N_("Tangent of x.")                           },
    { "tanh",  "tanh ( x )",                            N_("Hyperbolic tangent of x.")                },
    { "trunc", "trunc ( x )",                           N_("Truncated integer part of x.")            }
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
    { ">",   "x > y",                          N_("Greater Than relational operator.")                },
    { ">=",  "x >= y",                         N_("Greater Than Or Equal relational operator.")       },
    { "==",  "x == y",                         N_("Equal To relational operator.")                    },
    { "!=",  "x != y",                         N_("Not Equal To relational operator.")                },
    { "&&",  "x &amp;&amp; y",                 N_("Logical AND operator.")                            },
    { "||",  "x || y",                         N_("Logical OR operator.")                             }
};

#define MAX_FUNCTIONS G_N_ELEMENTS(functions)
#define MAX_OPERATORS G_N_ELEMENTS(operators)


enum {
	COLUMN_IMAGE_NUM,		// string
	COLUMN_IMAGE_PATH,		// string
	COLUMN_IMAGE_CHAN,      // uint
	COLUMN_IMAGE_WIDTH,      // uint
	COLUMN_IMAGE_HEIGHT,      // uint
	N_COLUMNS
};

enum {
	COLUMN_NAME,		 // string
	COLUMN_INDEX
};

static fits var_fit[MAX_IMAGES] = { 0 };

static GtkListStore *pixel_math_list_store = NULL;
static GtkListStore *pixel_math_list_store_functions = NULL;
static GtkListStore *pixel_math_list_store_operators = NULL;
static GtkTreeView *pixel_math_tree_view = NULL;
static GtkTreeView *pixel_math_treeview_functions = NULL;
static GtkTreeView *pixel_math_treeview_operators = NULL;
static GtkTreeModel *pixel_math_tree_model = NULL;
static GtkTreeModel *pixel_math_tree_model_functions = NULL;
static GtkTreeModel *pixel_math_tree_model_operators = NULL;
static GtkLabel *pixel_math_status_bar = NULL;
static GtkLabel *pixel_math_status_bar2 = NULL;
static GtkEntry *pixel_math_entry_r = NULL;
static GtkEntry *pixel_math_entry_g = NULL;
static GtkEntry *pixel_math_entry_b = NULL;

static void init_widgets() {
	if (!pixel_math_tree_view) {
		pixel_math_tree_view = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview"));
		pixel_math_tree_model = gtk_tree_view_get_model(pixel_math_tree_view);
		pixel_math_list_store = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore"));
		pixel_math_status_bar = GTK_LABEL(lookup_widget("pixel_math_status"));
		pixel_math_status_bar2 = GTK_LABEL(lookup_widget("pixel_math_status2"));
		pixel_math_entry_r = GTK_ENTRY(lookup_widget("pixel_math_entry_r"));
		pixel_math_entry_g = GTK_ENTRY(lookup_widget("pixel_math_entry_g"));
		pixel_math_entry_b = GTK_ENTRY(lookup_widget("pixel_math_entry_b"));
		pixel_math_treeview_functions = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview_functions"));
		pixel_math_treeview_operators = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview_operators"));
		pixel_math_tree_model_functions = gtk_tree_view_get_model(pixel_math_treeview_functions);
		pixel_math_tree_model_operators = gtk_tree_view_get_model(pixel_math_treeview_operators);
		pixel_math_list_store_functions = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_functions"));
		pixel_math_list_store_operators = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_operators"));

#if GTK_CHECK_VERSION(3, 22, 0)
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(lookup_widget("pixel_math_scrolled_functions")), TRUE);
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(lookup_widget("pixel_math_scrolled_operators")), TRUE);
#endif

	}
	g_assert(pixel_math_tree_view);
	g_assert(pixel_math_tree_model);
	g_assert(pixel_math_list_store);
	g_assert(pixel_math_list_store_functions);
	g_assert(pixel_math_list_store_operators);
	g_assert(pixel_math_status_bar);
	g_assert(pixel_math_status_bar2);
	g_assert(pixel_math_entry_r);
	g_assert(pixel_math_entry_g);
	g_assert(pixel_math_entry_b);
	g_assert(pixel_math_treeview_functions);
	g_assert(pixel_math_tree_model_functions);
	g_assert(pixel_math_treeview_operators);
	g_assert(pixel_math_tree_model_operators);

}


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

static void output_status_bar(int status) {
	switch (status) {
	case 0:
		gtk_label_set_text(pixel_math_status_bar, "");
		break;
	default:
		gtk_label_set_text(pixel_math_status_bar, _("Syntax error"));
	}
}

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

static gboolean check_files_dimensions(guint *width, guint* height, guint *channel) {
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
		valid = gtk_tree_model_iter_next (pixel_math_tree_model, &iter);

	}
	return TRUE;
}

static gboolean end_pixel_math_operation(gpointer p) {
	struct pixel_math_data *args = (struct pixel_math_data *)p;
	stop_processing_thread();// can it be done here in case there is no thread?

	if (!args->ret) {
		if (sequence_is_loaded())
			close_sequence(FALSE);

		//	/* Create new image */
		com.seq.current = UNRELATED_IMAGE;
		com.uniq = calloc(1, sizeof(single));
		com.uniq->filename = strdup(_("new empty image"));
		com.uniq->fileexist = FALSE;
		com.uniq->nb_layers = args->fit->naxes[2];
		com.uniq->fit = args->fit;
		clearfits(&gfit);
		copyfits(args->fit, &gfit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		invalidate_gfit_histogram();
		open_single_image_from_gfit();
	}
	set_cursor_waiting(FALSE);
	g_free(args->expression1);
	if (args->single_rgb) {
		g_free(args->expression2);
		g_free(args->expression3);
	}
	output_status_bar(args->ret);

	free(args);
	return FALSE;
}

void on_pm_use_rgb_button_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("pixel_math_entry_g"), !gtk_toggle_button_get_active(button));
	gtk_widget_set_sensitive(lookup_widget("pixel_math_entry_b"), !gtk_toggle_button_get_active(button));
	gtk_label_set_text(GTK_LABEL(lookup_widget("label_RGBK")), gtk_toggle_button_get_active(button) ? _("RGB/K") : _("R"));
	if (gtk_toggle_button_get_active(button)) entry_has_focus = 0;
}

static const gchar *get_pixel_math_var_paths(int i) {
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;

	init_widgets();

	GtkTreePath *path = gtk_tree_path_new_from_indices(i, -1);
	gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path);
	gtk_tree_model_get_value(pixel_math_tree_model, &iter, COLUMN_IMAGE_PATH, &value);

	return g_value_get_string(&value);
}

static int get_pixel_math_number_of_rows(){
	if (GTK_IS_TREE_MODEL(pixel_math_list_store))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store), NULL);
	else return 0;
}

static int get_pixel_math_functions_number_of_rows(){
	if (GTK_IS_TREE_MODEL(pixel_math_list_store_functions))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store_functions), NULL);
	else return 0;
}

static int get_pixel_math_operators_number_of_rows(){
	if (GTK_IS_TREE_MODEL(pixel_math_list_store_operators))
		return gtk_tree_model_iter_n_children(GTK_TREE_MODEL(pixel_math_list_store_operators), NULL);
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

static gboolean is_pm_use_rgb_button_checked() {
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pm_use_rgb_button")));
}

static gpointer apply_pixel_math_operation(gpointer p) {
	struct pixel_math_data *args = (struct pixel_math_data *)p;

	fits *fit = args->fit;
	int nb_rows = args->nb_rows;
	const gchar *expression1 = args->expression1;
	const gchar *expression2 = args->expression2;
	const gchar *expression3 = args->expression3;

	int err;
	gboolean failed = FALSE;

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
	{
		te_variable *vars = malloc(nb_rows * sizeof(te_variable));
		double *x = malloc(nb_rows * sizeof(double));

		for (int i = 0; i < nb_rows; i++) {
			vars[i].name = get_pixel_math_var_name(i);
			vars[i].address = &x[i];
			vars[i].context = NULL;
			vars[i].type = 0;
		}
		te_expr *n1 = NULL, *n2 = NULL, *n3 = NULL;
		n1 = te_compile(expression1, vars, nb_rows, &err);
		if (!args->single_rgb) {
			n2 = te_compile(expression2, vars, nb_rows, &err);
			n3 = te_compile(expression3, vars, nb_rows, &err);
		}
		if (!n1 || (!args->single_rgb && (!n2 || !n3))) {
			failed = TRUE;
		} else {
			size_t px;
#ifdef _OPENMP
#pragma omp for private(px) schedule(static)
#endif
			for (px = 0; px < var_fit[0].naxes[0] * var_fit[0].naxes[1] * var_fit[0].naxes[2]; px++) {
				/* The variables can be changed here, and eval can be called as many
				 * times as you like. This is fairly efficient because the parsing has
				 * already been done. */
				for (int i = 0; i < nb_rows; i++) {
					x[i] = var_fit[i].fdata[px];
				}

				if (args->single_rgb) {
					fit->fdata[px] = (float) te_eval(n1);
				} else {
					fit->fpdata[RLAYER][px] = (float) te_eval(n1);
					fit->fpdata[GLAYER][px] = (float) te_eval(n2);
					fit->fpdata[BLAYER][px] = (float) te_eval(n3);
				}
			}

			te_free(n1);
			if (!args->single_rgb) {
				te_free(n2);
				te_free(n3);
			}
			args->ret = 0;
		}
		free(vars);
		free(x);
	}
	if (failed) {
		args->ret = err;
	}

	siril_add_idle(end_pixel_math_operation, args);
	return GINT_TO_POINTER(args->ret);
}

static gboolean is_op_or_null(const gchar c) {
	if (c == '\0') return TRUE;
	if (c == '(') return TRUE;
	if (c == ')') return TRUE;
	for (int i = 0; i < MAX_OPERATORS; i++) {
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
			/* Only match the empty string once at any given position, to
			 * avoid infinite loops */
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
	GtkEntry *entry_param = GTK_ENTRY(lookup_widget("pixel_math_entry_param"));

	char *entry_text = g_strdup(gtk_entry_get_text(entry_param));

	if (entry_text == NULL) return 0;
	/* first we remove spaces */
	remove_spaces_from_str(entry_text);

	/* all parameters are comma separated expressions */
	gchar **token = g_strsplit(entry_text, ",", -1);
	int nargs = g_strv_length(token);

	/* now we pare equality */
	for (int i = 0; i < nargs; i++) {
		gchar **expr = g_strsplit(token[i], "=", -1);
		int n = g_strv_length(expr);
		/* we want something like "expr[0]=expr[1]"
		 * so we need two tokens */
		if (n != 2) {
			g_strfreev(token);
			g_strfreev(expr);
			g_free(entry_text);
			return -1;
		}

		/* We copy original char* in a string structure */
		GString *string1 = g_string_new(*expression1);
		GString *string2 = g_string_new(*expression2);
		GString *string3 = g_string_new(*expression3);

		/* we replace old expression by new ones */
		siril_string_replace_parameter(string1, expr[0], expr[1]);
		siril_string_replace_parameter(string2, expr[0], expr[1]);
		siril_string_replace_parameter(string3, expr[0], expr[1]);

		/* copy string into char */
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

static int pixel_math_evaluate(gchar *expression1, gchar *expression2, gchar *expression3) {
	int nb_rows = 0;

	int width = -1;
	int height = -1;
	int channel = -1;

	gboolean single_rgb = is_pm_use_rgb_button_checked();

	while (nb_rows < get_pixel_math_number_of_rows() && nb_rows < MAX_IMAGES) {
		const gchar *path = get_pixel_math_var_paths(nb_rows);
		if (readfits(path, &var_fit[nb_rows], NULL, TRUE)) return -1;
		if (channel == - 1) {
			width = var_fit[nb_rows].rx;
			height = var_fit[nb_rows].ry;
			channel = var_fit[nb_rows].naxes[2];
		}

		if (channel == 3 && !single_rgb) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Incompatible parameters"),
					_("3 channel images are incompatible with the \"Use single RGB/K expression\" unchecked."));
			return 1;
		}
		nb_rows++;
	}

	if (!nb_rows) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No images loaded"), _("You must load images first."));
		return 1;
	}

	channel = single_rgb ? channel : 3;

	fits *fit = NULL;
	if (new_fit_image(&fit, width, height, channel, DATA_FLOAT))
		return 1;

	struct pixel_math_data *args = malloc(sizeof(struct pixel_math_data));

	if (parse_parameters(&expression1, &expression2, &expression3)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Parameter error"),
				_("Parameter symbols could not be parsed."));
		return 1;
	}

	args->expression1 = expression1;
	args->expression2 = single_rgb ? NULL : expression2;
	args->expression3 = single_rgb ? NULL : expression3;
	args->single_rgb = single_rgb;
	args->fit = fit;
	args->nb_rows = nb_rows;
	args->ret = 0;

	start_in_new_thread(apply_pixel_math_operation, args);

	return 0;
}

static gboolean allowed_char(char c) {
	return (g_ascii_isalnum(c) || c == '_');
}

static gboolean check_for_variable_sanity(char *new_text) {
	gboolean abort = TRUE;
	const char *p = new_text;

	if (*p == '\0') return FALSE;
	/* Exclude if non alphanum  (or _) variable */
	while (*p) {
		if (!allowed_char(*p++)) return FALSE;
	}

	/* no we exclude name that contain only digit */
	p = new_text;

	while (*p) {
		if (g_ascii_isalpha(*p++)) {
			/* at least one char is alphanum */
			abort = FALSE;
		}
	}

	if (abort) return FALSE;

	init_widgets();

	/* Exclude duplicate names */
	int nb_rows = 0;
	while (nb_rows < get_pixel_math_number_of_rows() && nb_rows < MAX_IMAGES) {
		const gchar *var = get_pixel_math_var_name(nb_rows);
		if (!g_strcmp0(var, new_text)) return FALSE;
		nb_rows++;
	}
	return TRUE;
}

/* Add an image to the list. */
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

static void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title,
		const gchar *pattern) {
	gchar **patterns;
	gint i;

	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, title);
	/* get the patterns */
	patterns = g_strsplit(pattern, ";", -1);
	for (i = 0; patterns[i] != NULL; i++)
		gtk_file_filter_add_pattern(f, patterns[i]);
	/* free the patterns */
	g_strfreev(patterns);
	gtk_file_chooser_add_filter(file_chooser, f);
	gtk_file_chooser_set_filter(file_chooser, f);
}

static int search_for_free_index() {
	int i;
	for (i = 0; i < MAX_IMAGES; i++) {
		gboolean found = FALSE;
		for (int j = 0; j < MAX_IMAGES; j++) {
			const gchar *var = get_pixel_math_var_name(j);
			if (!var) break;
			if (!g_strcmp0(var, variables[i])) {
				found = TRUE;
				break;
			}
		}
		if (!found) {
			return i;
		}
	}
	return i + 1;
}

static void select_image(int nb) {
	GtkWidget *dialog;
	fileChooserPreview *preview = NULL;
	GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
	gint res;

	dialog = gtk_file_chooser_dialog_new("Open File", GTK_WINDOW(lookup_widget("pixel_math_dialog")), action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Open"), GTK_RESPONSE_ACCEPT,	NULL);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), com.wd);
	gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(dialog), TRUE);
	gtk_filter_add(GTK_FILE_CHOOSER(dialog), _("FITS Files (*.fit, *.fits, *.fts, *.fits.fz)"), FITS_EXTENSIONS);
	siril_file_chooser_add_preview(GTK_FILE_CHOOSER(dialog), preview);

	res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		GSList *l;
		int pos = nb;
		guint width = 0;
		guint height = 0;
		guint channel = 0;

		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		GSList *filenames = gtk_file_chooser_get_filenames(chooser);

		for (l = filenames; l; l = l->next) {
			char *filename;

			filename = (char *) l->data;
			if (filename) {
				gchar filter[FLEN_VALUE] = { 0 };
				fits f = { 0 };
				read_fits_metadata_from_path(filename, &f);
				if (check_files_dimensions(&width, &height, &channel)) {
					if (width != 0 && (channel != f.naxes[2] ||
							width != f.naxes[0] ||
							height != f.naxes[1])) {
						gchar *f = g_path_get_basename(filename);
						gchar *str = g_strdup_printf("%s will not be added in the pixel math tool because its size is different from the other loaded images"
								" (width, height or number of channels).", f);
						siril_message_dialog(GTK_MESSAGE_ERROR, _("Image must have same dimension"), str);
						g_free(f);
						g_free(str);

					} else {
						if (f.filter[0] != '\0') {
							memcpy(filter, f.filter, FLEN_VALUE);
						}

						int idx = search_for_free_index();
						add_image_to_variable_list(filename, variables[idx], filter, f.naxes[2], f.naxes[0], f.naxes[1]);

						pos++;
						if (pos == MAX_IMAGES) {
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

static GtkEntry *get_entry_with_focus() {
	switch(entry_has_focus) {
	default:
	case 0:
		return pixel_math_entry_r;
	case 1:
		return pixel_math_entry_g;
	case 2:
		return pixel_math_entry_b;
	}
}

void on_pixel_math_add_var_button_clicked(GtkButton *button, gpointer user_data) {
	init_widgets();

	int nb = get_pixel_math_number_of_rows();
	if (nb == MAX_IMAGES) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Cannot load new image"),
				_("You've reached the maximum of loaded image file."));
	} else {
		select_image(nb);
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
			if (gtk_tree_model_get_iter(pixel_math_tree_model, &iter, path)) {
				gtk_list_store_remove(pixel_math_list_store, &iter);
			}
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
}

void on_pixel_math_remove_var_button_clicked(GtkButton *button, gpointer user_data) {
	remove_selected_lines();
}

static void apply_pixel_math() {
	gchar *expression1 = get_pixel_math_expression1();
	remove_spaces_from_str(expression1);

	gchar *expression2 = get_pixel_math_expression2();
	remove_spaces_from_str(expression2);

	gchar *expression3 = get_pixel_math_expression3();
	remove_spaces_from_str(expression3);

	pixel_math_evaluate(expression1, expression2, expression3);
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
	gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_pixel_math_var_name(i[0]);

	if (str) {
		guint position = gtk_editable_get_position(GTK_EDITABLE(entry));
		gtk_editable_delete_selection (GTK_EDITABLE(entry));

		gtk_entry_buffer_insert_text(buffer, position, str, -1);
		gtk_widget_grab_focus(GTK_WIDGET(entry));
		gtk_editable_set_position(GTK_EDITABLE(entry), position + strlen(str));
	}
}

static gint get_real_index_from_index_in_list(GtkTreeModel *model, GtkTreeIter *iter) {
	gint real_index;
	gtk_tree_model_get(model, iter, COLUMN_INDEX, &real_index, -1);
	return real_index;
}

gboolean query_tooltip_tree_view_cb(GtkWidget *widget, gint x, gint y,
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

	g_snprintf(buffer, 511, "<b>%s</b>\n\n%s", functions[real_index].prototype, _(functions[real_index].definition));
	gtk_tooltip_set_markup(tooltip, buffer);

	gtk_tree_view_set_tooltip_row(tree_view, tooltip, path);

	gtk_tree_path_free(path);

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

static void add_functions_to_list() {
	GtkTreeIter iter;

	init_widgets();
	for (int i = 0; i < MAX_FUNCTIONS; i++) {
		gtk_list_store_append(pixel_math_list_store_functions, &iter);
		gtk_list_store_set(pixel_math_list_store_functions, &iter,
				COLUMN_NAME, functions[i].name, COLUMN_INDEX, i,
				-1);
	}
}

static void add_operators_to_list() {
	GtkTreeIter iter;

	init_widgets();
	for (int i = 0; i < MAX_OPERATORS; i++) {
		gtk_list_store_append(pixel_math_list_store_operators, &iter);
		gtk_list_store_set(pixel_math_list_store_operators, &iter,
				COLUMN_NAME, operators[i].name, COLUMN_INDEX, i,
				-1);
	}
}

void on_pixel_math_dialog_show(GtkWidget *w, gpointer user_data) {
	if (!get_pixel_math_functions_number_of_rows())
		add_functions_to_list();
	if (!get_pixel_math_operators_number_of_rows())
		add_operators_to_list();
}

void on_pixel_math_treeview_functions_row_activated(GtkTreeView *tree_view,
		GtkTreePath *path, GtkTreeViewColumn *column) {
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_function_name(i[0]);

	if (str) {
		guint position = gtk_editable_get_position(GTK_EDITABLE(entry));
		gtk_editable_delete_selection (GTK_EDITABLE(entry));

		gtk_entry_buffer_insert_text(buffer, position, str, -1);
		gtk_widget_grab_focus(GTK_WIDGET(entry));
		gtk_editable_set_position(GTK_EDITABLE(entry), position + strlen(str));
	}
}

void on_pixel_math_treeview_operators_row_activated(GtkTreeView *tree_view,
		GtkTreePath *path, GtkTreeViewColumn *column) {
	GtkEntry *entry = get_entry_with_focus();
	GtkEntryBuffer *buffer = gtk_entry_get_buffer(entry);
	gint *i = gtk_tree_path_get_indices(path);
	const gchar *str = get_operator_name(i[0]);

	if (str) {
		guint position = gtk_editable_get_position(GTK_EDITABLE(entry));
		gtk_editable_delete_selection (GTK_EDITABLE(entry));

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

	g_signal_handlers_unblock_by_func(lookup_widget("pixel_math_treeview"), on_pixel_math_treeview_key_release_event, NULL);

	if (check_for_variable_sanity(new_text)) {
		/* Value looks fine, we copy it */
		GtkTreeIter iter;
		GValue value = G_VALUE_INIT;

		g_value_init(&value, G_TYPE_STRING);
		g_assert(G_VALUE_HOLDS_STRING(&value));

		gtk_tree_model_get_iter_from_string(pixel_math_tree_model, &iter, path);
		g_value_set_string(&value, new_text);

		gtk_list_store_set_value(pixel_math_list_store, &iter, COLUMN_IMAGE_NUM, &value);
		g_value_unset (&value);
	}
}

void on_cellrenderer_variables_editing_started(GtkCellRenderer *renderer,
		GtkCellEditable *editable, char *path, gpointer user_data) {
	g_signal_handlers_block_by_func(lookup_widget("pixel_math_treeview"), on_pixel_math_treeview_key_release_event, NULL);
}

void on_cellrenderer_variables_editing_canceled(GtkCellRenderer *renderer,
		gpointer user_data) {
	g_signal_handlers_unblock_by_func(lookup_widget("pixel_math_treeview"), on_pixel_math_treeview_key_release_event, NULL);
}

void on_pixel_math_selection_changed(GtkTreeSelection *selection, gpointer user_data) {
	GList *list;
	guint width = 0, height = 0, channel = 0;

	list = gtk_tree_selection_get_selected_rows(selection, &pixel_math_tree_model);
	if (g_list_length(list) == 1) {
		GtkTreeIter iter;
		gtk_tree_model_get_iter(pixel_math_tree_model, &iter, (GtkTreePath*) list->data);

		gtk_tree_model_get(pixel_math_tree_model, &iter, COLUMN_IMAGE_CHAN, &channel, COLUMN_IMAGE_WIDTH, &width, COLUMN_IMAGE_HEIGHT, &height, -1);
	}
	g_list_free_full(list, (GDestroyNotify) gtk_tree_path_free);
	output_status_bar2(width, height, channel);
}
