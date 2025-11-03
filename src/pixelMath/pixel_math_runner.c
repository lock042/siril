/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/initfile.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/dialog_preview.h"
#include "gui/message_dialog.h"
#include "gui/open_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/histogram.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"

#include "tinyexpr.h"
#include "pixel_math_runner.h"

#define T_CURRENT "gfit"

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
    { "pi",    "pi",                                    N_("The constant \u03c0=3.141592...")         },
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

#define MAX_FUNCTIONS G_N_ELEMENTS(functions)
#define MAX_IMAGE_FUNCTIONS G_N_ELEMENTS(image_functions)
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
static gboolean var_fit_mask[MAX_IMAGES] = { 0 };

static GtkListStore *pixel_math_list_store = NULL;
static GtkListStore *pixel_math_list_store_functions = NULL;
static GtkListStore *pixel_math_list_store_operators = NULL;
static GtkListStore *pixel_math_list_store_presets = NULL;
static GtkTreeView *pixel_math_tree_view = NULL;
static GtkTreeView *pixel_math_treeview_functions = NULL;
static GtkTreeView *pixel_math_treeview_operators = NULL;
static GtkTreeView *pixel_math_treeview_presets = NULL;
static GtkTreeModel *pixel_math_tree_model = NULL;
static GtkTreeModel *pixel_math_tree_model_functions = NULL;
static GtkTreeModel *pixel_math_tree_model_operators = NULL;
static GtkTreeModel *pixel_math_tree_model_presets = NULL;
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
		pixel_math_treeview_presets = GTK_TREE_VIEW(gtk_builder_get_object(gui.builder, "pixel_math_treeview_presets"));
		pixel_math_tree_model_functions = gtk_tree_view_get_model(pixel_math_treeview_functions);
		pixel_math_tree_model_operators = gtk_tree_view_get_model(pixel_math_treeview_operators);
		pixel_math_tree_model_presets = gtk_tree_view_get_model(pixel_math_treeview_presets);
		pixel_math_list_store_functions = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_functions"));
		pixel_math_list_store_operators = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_operators"));
		pixel_math_list_store_presets = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "pixel_math_liststore_presets"));

#if GTK_CHECK_VERSION(3, 22, 0)
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(lookup_widget("pixel_math_scrolled_functions")), TRUE);
		gtk_scrolled_window_set_propagate_natural_width(GTK_SCROLLED_WINDOW(lookup_widget("pixel_math_scrolled_operators")), TRUE);
#endif

		/* Due to a glade bug, this property is often removed, lets code it */
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
	stop_processing_thread();

	if (!args->ret) {
		/* write to gfit in the graphical thread */
		clearfits(&gfit);
		if (sequence_is_loaded())
			close_sequence(FALSE);
		invalidate_gfit_histogram();

		memcpy(&gfit, args->fit, sizeof(fits));
		icc_auto_assign(&gfit, ICC_ASSIGN_ON_COMPOSITION);
		com.seq.current = UNRELATED_IMAGE;
		create_uniq_from_gfit(strdup(_("Pixel Math result")), FALSE);
		gui_function(open_single_image_from_gfit, NULL);
	}
	else clearfits(args->fit);

	set_cursor_waiting(FALSE);
	if (args->from_ui)
		output_status_bar(args->ret);

	free(args->fit);
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

static int get_pixel_math_presets_number_of_rows(){
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

static gboolean is_pm_use_rgb_button_checked() {
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("pm_use_rgb_button")));
}

static gboolean is_pm_rescale_checked() {
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("rescale_pm_button")));
}

static gboolean is_cumulate_checked() {
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("cumulate_pm_button")));
}

static float get_min_rescale_value() {
	return (float) gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_pm_low")));
}

static float get_max_rescale_value() {
	return (float) gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_pm_high")));
}

static void update_metadata(fits *fit, gboolean do_sum) {
	fits **f = malloc((MAX_IMAGES + 1) * sizeof(fits *));
	int j = 0;
	for (int i = 0; i < MAX_IMAGES ; i++)
		if (var_fit[i].rx > 0 && var_fit_mask[i])
			f[j++] = &var_fit[i];
	f[j] = NULL;

	if (!f[0] && single_image_is_loaded() )
		// if no fit used (only constants),
		// we copy the metadata from gfit
		copy_fits_metadata(&gfit, fit);
	else
		merge_fits_headers_to_result2(fit, f, do_sum);
	update_fits_header(fit);
	free(f);
}

static gchar* parse_image_functions(gpointer p, int idx, int c) {
	struct pixel_math_data *args = (struct pixel_math_data*) p;
	gchar *expression;
	gchar **image = args->varname;
	int nb_images = args->nb_rows;
	switch (idx) {
	case 1:
		expression = args->expression1;
		break;
	case 2:
		expression = args->expression2;
		break;
	case 3:
		expression = args->expression3;
		break;
	default:
		return NULL;
	}
	if (!expression)
		return expression;

	gchar *result = g_strdup(expression);
	GRegex *regex = g_regex_new("(\\w+)\\((\\w+)\\)", 0, 0, NULL);

	gboolean replaced = TRUE;
	while (replaced) {
		replaced = FALSE;
		GMatchInfo *match_info;
		g_regex_match(regex, result, 0, &match_info);

		if (g_match_info_matches(match_info)) {
			gchar *function = g_match_info_fetch(match_info, 1);
			gchar *param = g_match_info_fetch(match_info, 2);
			gchar *full_match = g_match_info_fetch(match_info, 0);
			double median = 0.0, mean = 0.0, min = 0.0, max = 0.0, noise = 0.0, adev = 0.0, bwmv = 0.0, mad = 0.0, sdev = 0.0;
			double w = 0.0, h = 0.0;
			imstats *stats = NULL;

			if (g_strcmp0(param, T_CURRENT) == 0) {
				stats = statistics(NULL, -1, &gfit, c, NULL, STATS_MAIN, MULTI_THREADED);
				if (!stats) {
					g_free(full_match);
					g_free(function);
					g_free(param);
					g_match_info_free(match_info);
					g_regex_unref(regex);
					return result;
				}
				median = stats->median;
				mean = stats->mean;
				min = stats->min;
				max = stats->max;
				noise = stats->bgnoise;
				adev = stats->avgDev;
				bwmv = stats->sqrtbwmv * stats->sqrtbwmv;
				mad = stats->mad;
				sdev = stats->sigma;
				w = (double) gfit.rx;
				h = (double) gfit.ry;
				free_stats(stats);
			} else {
				for (int j = 0; j < nb_images; j++) {
					if (g_strcmp0(param, image[j]) == 0) {
						stats = statistics(NULL, -1, &var_fit[j], c, NULL, STATS_MAIN, MULTI_THREADED);
						if (!stats) {
							g_free(full_match);
							g_free(function);
							g_free(param);
							g_match_info_free(match_info);
							g_regex_unref(regex);
							return result;
						}
						median = stats->median;
						mean = stats->mean;
						min = stats->min;
						max = stats->max;
						noise = stats->bgnoise;
						adev = stats->avgDev;
						bwmv = stats->sqrtbwmv * stats->sqrtbwmv;
						mad = stats->mad;
						sdev = stats->sigma;
						w = (double) var_fit[j].rx;
						h = (double) var_fit[j].ry;
						free_stats(stats);
						break;
					}
				}
			}

			gchar *replace = NULL;
			if (!g_strcmp0(function, "mean")) {
				replace = g_strdup_printf("%g", mean);
			} else if (!g_strcmp0(function, "med") || !g_strcmp0(function, "median")) {
				replace = g_strdup_printf("%g", median);
			} else if (!g_strcmp0(function, "min")) {
				replace = g_strdup_printf("%g", min);
			} else if (!g_strcmp0(function, "max")) {
				replace = g_strdup_printf("%g", max);
			} else if (!g_strcmp0(function, "noise")) {
				replace = g_strdup_printf("%g", noise);
			} else if (!g_strcmp0(function, "adev")) {
				replace = g_strdup_printf("%g", adev);
			} else if (!g_strcmp0(function, "bwmv")) {
				replace = g_strdup_printf("%g", bwmv);
			} else if (!g_strcmp0(function, "mad") || !g_strcmp0(function, "mdev")) {
				replace = g_strdup_printf("%g", mad);
			} else if (!g_strcmp0(function, "sdev")) {
				replace = g_strdup_printf("%g", sdev);
			} else if (!g_strcmp0(function, "width") || !g_strcmp0(function, "w")) {
				replace = g_strdup_printf("%g", w);
			} else if (!g_strcmp0(function, "height") || !g_strcmp0(function, "h")) {
				replace = g_strdup_printf("%g", h);
			}

			if (replace) {
				gchar *temp = result;
				// Replace the specific matched string with the calculated value
				gchar **split = g_strsplit(result, full_match, 2);
				if (split[0] && split[1]) {
					result = g_strconcat(split[0], replace, split[1], NULL);
				} else {
					result = g_strdup(result);
				}
				g_strfreev(split);
				g_free(temp);
				g_free(replace);
				replaced = TRUE;
				siril_debug_print("Expression%d: %s\n", c, result);
			}

			g_free(full_match);
			g_free(function);
			g_free(param);
		}

		g_match_info_free(match_info);
	}

	g_regex_unref(regex);

	for (int j = 0; j < nb_images; j++) {
		const gchar *test = g_strrstr(result, image[j]);
		if (test) {
			var_fit_mask[j] = TRUE;
			siril_debug_print("found image name %s in the expression %s\n", image[j], result);
		}
	}

	return result;
}

gpointer apply_pixel_math_operation(gpointer p) {
	struct pixel_math_data *args = (struct pixel_math_data *)p;

	te_expr *n1 = NULL, *n2 = NULL, *n3 = NULL;
	fits *fit = args->fit;
	int nb_rows = args->nb_rows;
	gboolean failed = FALSE;
	args->ret = 0;
	float maximum = -FLT_MAX;
	float minimum = +FLT_MAX;
	long width, height, nchan;

	if (args->single_rgb && args->fit->naxes[2] > 1) {
		// No need to null check these two as they will be NULL if args->single_rgb is TRUE
		args->expression2 = g_strdup(args->expression1);
		args->expression3 = g_strdup(args->expression1);

		args->expression1 = parse_image_functions(args, 1, RLAYER);
		args->expression2 = parse_image_functions(args, 2, GLAYER);
		args->expression3 = parse_image_functions(args, 3, BLAYER);
	} else {
		args->expression1 = parse_image_functions(args, 1, RLAYER);
		args->expression2 = parse_image_functions(args, 2, RLAYER);
		args->expression3 = parse_image_functions(args, 3, RLAYER);
	}

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) firstprivate(n1,n2,n3)
#endif
	{
		// we build the expressions in parallel because tr_eval() is not thread-safe
		int k = 0;
		if (args->has_gfit) k = 1;
		te_variable *vars = malloc((nb_rows + k) * sizeof(te_variable));
		double *x = malloc((nb_rows + k) * sizeof(double));
		if (!vars || !x) {
			failed = TRUE;
			goto failure;
		}
		for (int i = 0; i < nb_rows; i++) {
			vars[i].name = args->varname[i];
			vars[i].address = &x[i];
			vars[i].context = NULL;
			vars[i].type = 0;
		}
		if (args->has_gfit) {
			vars[nb_rows].name = g_strdup(T_CURRENT);
			vars[nb_rows].address = &x[nb_rows];
			vars[nb_rows].context = NULL;
			vars[nb_rows].type = 0;
		}
		int err = 0;
		n1 = te_compile(args->expression1, vars, nb_rows + k, &err);
		if (!n1) {
#ifdef _OPENMP
			if (omp_get_thread_num() == 0)
#endif
				siril_log_color_message(_("Error in pixel math expression '%s' at character %d\n"), "red", args->expression1, err);
			failed = TRUE;
			goto failure;
		}

		if (args->expression2) {
			n2 = te_compile(args->expression2, vars, nb_rows + k, &err);
			if (!n2) {
#ifdef _OPENMP
				if (omp_get_thread_num() == 0)
#endif
					siril_log_color_message(_("Error in pixel math expression '%s' at character %d\n"), "red", args->expression2, err);
				failed = TRUE;
				goto failure;
			}

			n3 = te_compile(args->expression3, vars, nb_rows + k, &err);
			if (!n3) {
#ifdef _OPENMP
				if (omp_get_thread_num() == 0)
#endif
					siril_log_color_message(_("Error in pixel math expression '%s' at character %d\n"), "red", args->expression3, err);
				failed = TRUE;
				goto failure;
			}
		}
		if (args->has_gfit && nb_rows == 0) {
			width = gfit.naxes[0];
			height = gfit.naxes[1];
			nchan = gfit.naxes[2];
		} else {
			width = var_fit[0].naxes[0];
			height = var_fit[0].naxes[1];
			nchan = var_fit[0].naxes[2];
		}
		if (com.pref.force_16bit) {
#ifdef _OPENMP
#pragma omp for schedule(static) reduction(max:maximum) reduction(min:minimum)
#endif
			for (size_t px = 0; px < width * height * nchan; px++) {
				/* The variables can be changed here, and eval can be called as many
				 * times as you like. This is fairly efficient because the parsing has
				 * already been done. */
				for (int i = 0; i < nb_rows; i++) {
					x[i] = (double) var_fit[i].fdata[px];
				}
				if (args->has_gfit) {
					if (gfit.type == DATA_USHORT) {
						x[nb_rows] = (double) gfit.data[px] / USHRT_MAX_DOUBLE;
					} else {
						x[nb_rows] = (double) gfit.fdata[px];
					}
				}

				if (!args->single_rgb) { // in that case var_fit[0].naxes[2] == 1, but we built RGB
					fit->pdata[RLAYER][px] =  roundf_to_WORD((float) te_eval(n1) * USHRT_MAX_SINGLE);
					fit->pdata[GLAYER][px] =  roundf_to_WORD((float) te_eval(n2) * USHRT_MAX_SINGLE);
					fit->pdata[BLAYER][px] =  roundf_to_WORD((float) te_eval(n3) * USHRT_MAX_SINGLE);

					/* may not be used but (only if rescale) at least it is computed */
					maximum = max(maximum, max(fit->pdata[RLAYER][px], max(fit->pdata[GLAYER][px], fit->pdata[BLAYER][px])));
					minimum = min(minimum, min(fit->pdata[RLAYER][px], min(fit->pdata[GLAYER][px], fit->pdata[BLAYER][px])));
				} else {
					if (px < (var_fit[0].naxes[0] * var_fit[0].naxes[1])) {
						fit->data[px] =  roundf_to_WORD((float) te_eval(n1) * USHRT_MAX_SINGLE);
					} else if (px < 2 * (var_fit[0].naxes[0] * var_fit[0].naxes[1])) {
						fit->data[px] =  roundf_to_WORD((float) te_eval(n2) * USHRT_MAX_SINGLE);
					} else {
						fit->data[px] =  roundf_to_WORD((float) te_eval(n3) * USHRT_MAX_SINGLE);
					}

					/* may not be used but (only if rescale) at least it is computed */
					maximum = max(maximum, fit->data[px]);
					minimum = min(minimum, fit->data[px]);
				}
			}

		} else {
#ifdef _OPENMP
#pragma omp for schedule(static) reduction(max:maximum) reduction(min:minimum)
#endif
			for (size_t px = 0; px < width * height * nchan; px++) {
				/* The variables can be changed here, and eval can be called as many
				 * times as you like. This is fairly efficient because the parsing has
				 * already been done. */
				for (int i = 0; i < nb_rows; i++) {
					x[i] = var_fit[i].fdata[px];
				}
				if (args->has_gfit) {
					if (gfit.type == DATA_USHORT) {
						x[nb_rows] = gfit.data[px] / USHRT_MAX_DOUBLE;
					} else {
						x[nb_rows] = gfit.fdata[px];
					}
				}

				if (!args->single_rgb) { // in that case var_fit[0].naxes[2] == 1, but we built RGB
					fit->fpdata[RLAYER][px] = (float) te_eval(n1);
					fit->fpdata[GLAYER][px] = (float) te_eval(n2);
					fit->fpdata[BLAYER][px] = (float) te_eval(n3);

					/* may not be used but (only if rescale) at least it is computed */
					maximum = max(maximum, max(fit->fpdata[RLAYER][px], max(fit->fpdata[GLAYER][px], fit->fpdata[BLAYER][px])));
					minimum = min(minimum, min(fit->fpdata[RLAYER][px], min(fit->fpdata[GLAYER][px], fit->fpdata[BLAYER][px])));
				} else {
					if (px < (width * height)) {
						fit->fdata[px] = (float) te_eval(n1);
					} else if (px < 2 * (width * height)) {
						fit->fdata[px] = (float) te_eval(n2);
					} else {
						fit->fdata[px] = (float) te_eval(n3);
					}

					/* may not be used but (only if rescale) at least it is computed */
					maximum = max(maximum, fit->fdata[px]);
					minimum = min(minimum, fit->fdata[px]);
				}
			}
		}

failure: // failure before the eval loop
		te_free(n1);
		if (args->expression2) {
			te_free(n2);
			te_free(n3);
		}
		free(vars);
		free(x);
	} // end of parallel block

	if (args->rescale) {
		if (com.pref.force_16bit) {
			args->min *= USHRT_MAX_SINGLE;
			args->max *= USHRT_MAX_SINGLE;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (int i = 0; i < fit->rx * fit->ry * fit->naxes[2]; i++) {
				fit->data[i] = roundf_to_WORD((float)(args->max - args->min) * (float)(fit->data[i] - minimum) / (float)(maximum - minimum) + args->min);
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (int i = 0; i < fit->rx * fit->ry * fit->naxes[2]; i++) {
				fit->fdata[i] = (((args->max - args->min) * (fit->fdata[i] - minimum)) / (maximum - minimum)) + args->min;
			}
		}
	}

	if (failed)
		args->ret = 1;
	else update_metadata(args->fit, args->do_sum);

	/* free memory */
	g_free(args->expression1);
	if (args->expression2) {
		g_free(args->expression2);
		g_free(args->expression3);
	}
	for (int i = 0; i < args->nb_rows; i++)
		g_free(args->varname[i]);
	free(args->varname);
	free_pm_var(args->nb_rows);

	/* manage result and display */
	if (com.headless) {
		// no display or threading needed
		if (!failed) {
			clearfits(&gfit);
			memcpy(&gfit, args->fit, sizeof(fits));
			com.seq.current = UNRELATED_IMAGE;
			create_uniq_from_gfit(strdup(_("Pixel Math result")), FALSE);
		}
		else clearfits(args->fit);
		free(args->fit);
		free(args);
	}
	else {
		execute_idle_and_wait_for_it(end_pixel_math_operation, args);
	}
	return GINT_TO_POINTER((gint)failed);
}

static gboolean is_op_or_null(const gchar c) {
	if (c == '\0') return TRUE;
	if (c == '(') return TRUE;
	if (c == ')') return TRUE;
	if (c == ',') return TRUE;
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

	/* now we parse equality */
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


int load_pm_var(const gchar *var, int index, int *w, int *h, int *c) {
	if (index > MAX_IMAGES - 1) {
		siril_log_message(_("A maximum of %d images can be used in a single expression.\n"), MAX_IMAGES);
		return 1;
	}

	if (readfits(var, &var_fit[index], NULL, TRUE)) {
		*w = *h = *c = -1;
		return 1;
	}
	*w = var_fit[index].rx;
	*h = var_fit[index].ry;
	*c = var_fit[index].naxes[2];
	return 0;
}

void free_pm_var(int nb) {
	for (int i = 0; i < nb; i++) {
		clearfits(&var_fit[i]);
		var_fit_mask[i] = FALSE;
	}
}

static int replace_t_with_gfit(struct pixel_math_data *args) {
	int retval = 0;
	/* Check if single_rgb is FALSE and if $T is present in any expression */
	if (args->single_rgb == FALSE &&
		((args->expression1 && g_strstr_len(args->expression1, -1, "$T") != NULL) ||
		 (args->expression2 && g_strstr_len(args->expression2, -1, "$T") != NULL) ||
		 (args->expression3 && g_strstr_len(args->expression3, -1, "$T") != NULL))) {
		retval = 1;
	}

	/* If we reach here, either single_RGB is FALSE or no $T is present */
	if ((args->expression1 && g_strstr_len(args->expression1, -1, "$T") != NULL) ||
		(args->expression2 && g_strstr_len(args->expression2, -1, "$T") != NULL) ||
		(args->expression3 && g_strstr_len(args->expression3, -1, "$T") != NULL)) {

		const gchar *pattern = "\\$T";
		const gchar *replacement = T_CURRENT;

		args->has_gfit = TRUE;

		if (args->expression1 != NULL) {
			gchar *new_expression1 = g_regex_replace(g_regex_new(pattern, 0, 0, NULL), args->expression1, -1, 0, replacement, 0, NULL);
			g_free(args->expression1);
			args->expression1 = new_expression1;
		}

		if (args->expression2 != NULL) {
			gchar *new_expression2 = g_regex_replace(g_regex_new(pattern, 0, 0, NULL), args->expression2, -1, 0, replacement, 0, NULL);
			g_free(args->expression2);
			args->expression2 = new_expression2;
		}

		if (args->expression3 != NULL) {
			gchar *new_expression3 = g_regex_replace(g_regex_new(pattern, 0, 0, NULL), args->expression3, -1, 0, replacement, 0, NULL);
			g_free(args->expression3);
			args->expression3 = new_expression3;
		}
	}

	/* Success - no errors encountered */
	return retval;
}

static int pixel_math_evaluate(gchar *expression1, gchar *expression2, gchar *expression3) {
	int nb_rows = 0;

	int width = -1;
	int height = -1;
	int channel = -1;

	int retval = 0;

	gboolean icc_warning_given = FALSE;
	gboolean single_rgb = is_pm_use_rgb_button_checked();
	gboolean rescale = is_pm_rescale_checked();
	gboolean do_sum = is_cumulate_checked();
	float min = get_min_rescale_value();
	float max = get_max_rescale_value();

	while (nb_rows < get_pixel_math_number_of_rows() && nb_rows < MAX_IMAGES) {
		const gchar *path = get_pixel_math_var_paths(nb_rows);
		if (readfits(path, &var_fit[nb_rows], NULL, TRUE)) {
			retval = 1;
			goto free_expressions;
		}

		// Check ICC profiles are defined and the same
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
				// This also takes care of the stuation where a file doesn't have an
				// ICC profile - in that case it is assigned a profile to match.
				// Ultimately it is much better to use images that all have the same
				// color profile - when they don't we're reduced to guesswork about the
				// user's intention.
				siril_colorspace_transform(&var_fit[nb_rows], var_fit[0].icc_profile);
			}
		}
		// Check channels are compatible
		if (channel == -1) {
			width = var_fit[nb_rows].rx;
			height = var_fit[nb_rows].ry;
			channel = var_fit[nb_rows].naxes[2];
		}

		if (channel == 3 && !single_rgb) {
			queue_message_dialog(GTK_MESSAGE_ERROR, _("Incompatible parameters"),
					_("3 channel images are incompatible with the \"Use single RGB/K expression\" unchecked."));
			retval = 1;
			goto free_expressions;
		}
		nb_rows++;
	}

	channel = single_rgb ? channel : 3;

	struct pixel_math_data *args = calloc(1, sizeof(struct pixel_math_data));

	if (parse_parameters(&expression1, &expression2, &expression3)) {
		queue_message_dialog(GTK_MESSAGE_ERROR, _("Parameter error"), _("Parameter symbols could not be parsed."));
		retval = 1;
		goto free_expressions;
	}

	args->expression1 = g_strdup(expression1);
	if (!single_rgb) {
		args->expression2 = g_strdup(expression2);
		args->expression3 = g_strdup(expression3);
	}

	args->single_rgb = single_rgb;
	args->rescale = rescale;
	args->do_sum = do_sum;
	args->min = min;
	args->max = max;
	args->nb_rows = nb_rows;
	args->ret = 0;
	args->from_ui = TRUE;
	args->has_gfit = FALSE;

	args->varname = malloc(nb_rows * sizeof(gchar *));
	for (int i = 0; i < nb_rows; i++) {
		args->varname[i] = g_strdup(get_pixel_math_var_name(i));
	}

	if (replace_t_with_gfit(args) && gfit.naxes[2] > 1) {
		queue_message_dialog(GTK_MESSAGE_ERROR, _("Incorrect formula"), _("RGB $T cannot be used in this context."));
		retval = 1;
		goto free_expressions;
	}

	fits *fit = NULL;
	if (args->has_gfit) { // in the case where no images are loaded, at least gfit must be laded
		width = gfit.rx;
		height = gfit.ry;
		channel = args->single_rgb ? gfit.naxes[2] : 3;
		if (nb_rows > 0) {
			if (width != var_fit[0].naxes[0] ||
				height != var_fit[0].naxes[1]) {
					queue_message_dialog(GTK_MESSAGE_ERROR, _("Images have different size"),
						_("The image currently displayed must be the same size as the other images loaded into PixelMath."));
				retval = 1;
				goto free_expressions;
			}
		}
	}

	if (!nb_rows && !args->has_gfit) {
		queue_message_dialog(GTK_MESSAGE_ERROR, _("No images loaded"), _("You must load images first."));
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
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(dialog), FALSE);
	gtk_filter_add(GTK_FILE_CHOOSER(dialog), _("FITS Files (*.fit, *.fits, *.fts, *.fit.fz, *.fits.fz, *.fts.fz)"), FITS_EXTENSIONS, gui.file_ext_filter == TYPEFITS);
	siril_file_chooser_add_preview(GTK_FILE_CHOOSER(dialog), preview);

	res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		GSList *l;
		int pos = nb;
		guint width = 0;
		guint height = 0;
		guint channel = 0;

		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		GSList *filenames = siril_file_chooser_get_filenames(chooser);

		for (l = filenames; l; l = l->next) {
			char *filename;

			filename = (char *) l->data;
			if (filename) {
				gchar filter[FLEN_VALUE] = { 0 };
				fits f = { 0 };
				if (read_fits_metadata_from_path(filename, &f)) {
					siril_log_color_message(_("Could not open file: %s\n"), "red", filename);
				} else if (check_files_dimensions(&width, &height, &channel)) {
					if (width != 0 && (channel != f.naxes[2] ||
							width != f.naxes[0] ||
							height != f.naxes[1])) {
						gchar *name = g_path_get_basename(filename);
						gchar *str = g_strdup_printf("%s will not be added in the pixel math tool because its size is different from the other loaded images"
								" (width, height or number of channels).", name);
						siril_message_dialog(GTK_MESSAGE_ERROR, _("Image must have same dimension"), str);
						g_free(name);
						g_free(str);

					} else {
						if (f.keywords.filter[0] != '\0') {
							memcpy(filter, f.keywords.filter, FLEN_VALUE);
						}

						int idx = search_for_free_index();
						if (idx >= MAX_IMAGES) {
							siril_log_color_message(_("Error: maximum variable index exceeded - too many variables!\n"), "red");
							g_free(filename);
							clearfits(&f);
							break;
						}
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
	if (get_thread_run()) {
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

static void add_functions_to_list() {
	GtkTreeIter iter;

	_pm_op_func *all_functions = concat_functions();

	init_widgets();
	int i = 0;
	/* fill functions */
	for (i = 0; i < MAX_FUNCTIONS + MAX_IMAGE_FUNCTIONS; i++) {
		gtk_list_store_append(pixel_math_list_store_functions, &iter);
		gtk_list_store_set(pixel_math_list_store_functions, &iter, COLUMN_NAME, all_functions[i].name, COLUMN_INDEX, i, -1);
	}
	free(all_functions);
}

static void add_operators_to_list() {
	GtkTreeIter iter;

	init_widgets();
	for (int i = 0; i < MAX_OPERATORS; i++) {
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
	const gint *i = gtk_tree_path_get_indices(path);
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

	return FALSE; /* do not stop walking the store, call us with next row */
}

static void save_presets_list() {
	/* First we free the original list */
	g_slist_free_full(com.pref.gui.pm_presets, g_free);
	com.pref.gui.pm_presets = NULL;

	gtk_tree_model_foreach(GTK_TREE_MODEL(pixel_math_list_store_presets), foreach_func, com.pref.gui.pm_presets);
	com.pref.gui.pm_presets = g_slist_reverse(com.pref.gui.pm_presets);
	writeinitfile();
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
		gtk_editable_delete_selection (GTK_EDITABLE(entry));

		gtk_entry_buffer_insert_text(buffer, position, str, -1);
		gtk_widget_grab_focus(GTK_WIDGET(entry));
		gtk_editable_set_position(GTK_EDITABLE(entry), position + strlen(str));
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
			if (gtk_tree_model_get_iter(pixel_math_tree_model_presets, &iter, path)) {
				gtk_list_store_remove(pixel_math_list_store_presets, &iter);
			}
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
	save_presets_list();
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
	gtk_widget_set_sensitive(lookup_widget("spin_pm_low"), gtk_toggle_button_get_active(button));
	gtk_widget_set_sensitive(lookup_widget("spin_pm_high"), gtk_toggle_button_get_active(button));
}
