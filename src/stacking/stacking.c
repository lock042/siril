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

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <sys/stat.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/initfile.h"
#include "core/OS_utils.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "core/arithm.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "gui/sequence_list.h"
#include "gui/registration_preview.h"
#include "gui/stacking.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/ser.h"
#include "registration/registration.h"
#include "algos/noise.h"
#include "algos/sorting.h"
#include "algos/siril_wcs.h"
#include "stacking/sum.h"
#include "opencv/opencv.h"

#include "stacking.h"

static gboolean end_stacking(gpointer p);

gboolean evaluate_stacking_should_output_32bits(const stack_method method,
		sequence *seq, int nb_img_to_stack, gchar **err) {
	gchar *error = NULL;
	if (com.pref.force_16bit) {
		if (seq->bitpix == FLOAT_IMG) {
			error = _("Input sequence is in 32-bit format but preferences are set to 16-bit output format. "
					"Please, change your preference settings and retry.\n");
		}
		if (err) {
			*err = error;
		}
		return FALSE;
	}
	if (method == stack_summing_generic) {
		if (seq->bitpix == BYTE_IMG)
			return nb_img_to_stack > 256;
		return TRUE;
	}
	if (method == stack_mean_with_rejection) {
		return TRUE;
	}
	if (method == stack_median) {
		return TRUE;
	}
	return seq->bitpix == FLOAT_IMG; // for min or max, only use it if input is already float
}

/* the function that prepares the stacking and runs it */
void main_stack(struct stacking_args *args) {
	int nb_allowed_files;
	g_assert(args->ref_image >= 0 && args->ref_image < args->seq->number);

	/* first of all we need to check if we can process the files */
	if (args->seq->type == SEQ_REGULAR && args->method != stack_summing_generic) {
		if (!allow_to_open_files(args->nb_images_to_stack, &nb_allowed_files)) {
			siril_log_message(_("Your system does not allow one to open more than %d files at the same time. "
						"You may consider either to enhance this limit (the method depends of "
						"your Operating System) or to convert your FITS sequence into a SER "
						"sequence before stacking, or to stack with the \"sum\" method.\n"),
					nb_allowed_files);
			args->retval = -1;
			return;
		}
	}

	siril_log_message(args->description);
	if (args->use_32bit_output)
		siril_log_message(_("Stacking result will be stored as a 32-bit image\n"));

	// 1. normalization
	if (do_normalization(args)) // does nothing if NO_NORM
		return;
	// 2. up-scale
	if (upscale_sequence(args)) // does nothing if !args->upscale_at_stacking
		return;
	// 3. stack
	args->retval = args->method(args);

	// result is in args->result, not saved
	describe_stack_for_history(args, &args->result.history, FALSE, FALSE);
}

/* the function that runs the thread. */
gpointer stack_function_handler(gpointer p) {
	struct stacking_args *args = (struct stacking_args *)p;

	main_stack(args);

	// 4. save result and clean-up
	siril_add_idle(end_stacking, args);
	return GINT_TO_POINTER(args->retval);
}

void _show_summary(struct stacking_args *args) {
	siril_log_message(_("Integration of %d images on %d of the sequence:\n"), args->nb_images_to_stack, args->seq->number);

	/* Type of algorithm */
	if (args->method == &stack_mean_with_rejection) {
		siril_log_message(_("Pixel combination ......... average\n"));
	} else if (args->method == &stack_summing_generic) {
		siril_log_message(_("Pixel combination ......... normalized sum\n"));
	} else if (args->method == &stack_median) {
		siril_log_message(_("Pixel combination ......... median\n"));
	} else if (args->method == &stack_addmin) {
		siril_log_message(_("Pixel combination ......... minimum\n"));
	} else if (args->method == &stack_addmax) {
		siril_log_message(_("Pixel combination ......... maximum\n"));
	} else {
		siril_log_message(_("Pixel combination ......... none\n"));
	}

	/* Normalisation */
	const char *norm_str, *fast_norm = "", *overlap = "";
	if (args->method != &stack_mean_with_rejection &&
			args->method != &stack_median ) {
		norm_str = _("none");
	} else {
		switch (args->normalize) {
		default:
		case NO_NORM:
			norm_str = _("none");
			break;
		case ADDITIVE:
			norm_str = _("additive");
			break;
		case MULTIPLICATIVE:
			norm_str = _("multiplicative");
			break;
		case ADDITIVE_SCALING:
			norm_str = _("additive + scaling");
			break;
		case MULTIPLICATIVE_SCALING:
			norm_str = _("multiplicative + scaling");
			break;
		}
		if (args->normalize != NO_NORM && args->lite_norm)
			fast_norm = _(" (fast)");
		if (args->normalize != NO_NORM && args->overlap_norm)
			overlap = _(" (overlaps)");
	}
	siril_log_message(_("Input normalization ....... %s%s%s\n"), norm_str, fast_norm, overlap);

	if (args->output_norm)
		siril_log_message(_("Output normalization ...... enabled\n"));
	else siril_log_message(_("Output normalization ...... disabled\n"));

	/* Rejection */
	if (args->method != &stack_mean_with_rejection) {
		siril_log_message(_("Pixel rejection ........... none\n"));
	} else {
		const char *rej_str;
		switch (args->type_of_rejection) {
		default:
		case NO_REJEC:
			rej_str = _("none");
			break;
		case PERCENTILE:
			rej_str = _("percentile clipping");
			break;
		case SIGMA:
			rej_str = _("sigma clipping");
			break;
		case MAD:
			rej_str = _("MAD clipping");
			break;
		case SIGMEDIAN:
			rej_str = _("median sigma clipping");
			break;
		case WINSORIZED:
			rej_str = _("winsorized sigma clipping");
			break;
		case LINEARFIT:
			rej_str = _("linear fit clipping");
			break;
		case GESDT:
			rej_str = _("GESDT clipping");
			break;
		}
		siril_log_message(_("Pixel rejection ........... %s\n"), rej_str);
		if (args->type_of_rejection != NO_REJEC) {
			if (args->type_of_rejection == GESDT) {
				siril_log_message(_("Rejection parameters ...... outliers=%.3f significance=%.3f\n"),
						args->sig[0], args->sig[1]);
			} else {
				siril_log_message(_("Rejection parameters ...... low=%.3f high=%.3f\n"),
						args->sig[0], args->sig[1]);
			}
			if (!args->create_rejmaps)
				siril_log_message(_("Creating rejection maps ... no\n"));
			else {
				if (args->merge_lowhigh_rejmaps)
					siril_log_message(_("Creating rejection maps ... yes (merged)\n"));
				else	siril_log_message(_("Creating rejection maps ... yes\n"));
			}
		}
	}

	if (args->weighting_type == NOISE_WEIGHT)
		siril_log_message(_("Image weighting ........... from noise\n"));
	else if (args->weighting_type == NBSTACK_WEIGHT)
		siril_log_message(_("Image weighting ........... from image count\n"));
	else if (args->weighting_type == WFWHM_WEIGHT)
		siril_log_message(_("Image weighting ........... from weighted FWHM\n"));
	else if (args->weighting_type == NBSTARS_WEIGHT)
		siril_log_message(_("Image weighting ........... from star count\n"));
	else siril_log_message(_("Image weighting ........... disabled\n"));

	if (args->feather_dist > 0)
		siril_log_message(_("Feathering ................ over %4d pixels\n"), args->feather_dist);

	if (args->seq->nb_layers > 1) {
		if (args->equalizeRGB)
			siril_log_message(_("RGB equalization .......... enabled\n"));
		else siril_log_message(_("RGB equalization .......... disabled\n"));
	}
}

/* a short version of the above, for FITS header HISTORY */
void describe_stack_for_history(struct stacking_args *args, GSList **hist, gboolean for_rejmap, gboolean low_rejmap) {
	const char *stack_name;
	if (args->method == &stack_summing_generic)
		stack_name = "sum";
	else if (args->method == &stack_mean_with_rejection)
		stack_name = "mean";
	else if (args->method == &stack_median)
		stack_name = "median";
	else if (args->method == &stack_addmin)
		stack_name = "minimum";
	else if (args->method == &stack_addmax)
		stack_name = "maximum";
	else stack_name = "unknown";

	GString *str;
	if (for_rejmap) {
		g_assert(args->create_rejmaps);
		if (args->merge_lowhigh_rejmaps)
			str = g_string_new("merged low+high rejection map for a ");
		else {
			if (low_rejmap)
				str = g_string_new("low rejection map for a ");
			else str = g_string_new("high rejection map for a ");
		}
		g_string_append(str, stack_name);
		g_string_append(str, " stacking");
	} else {
		str = g_string_new(stack_name);
		g_string_append(str, " stacking");
	}

	/* Type of rejection */
	if (args->method != &stack_mean_with_rejection) {
		g_string_append(str, " without rejection");
	}
	else {
		switch (args->type_of_rejection) {
		default:
		case NO_REJEC:
			g_string_append(str, " without rejection");
			break;
		case PERCENTILE:
			g_string_append(str, " with percentile clipping rejection");
			break;
		case SIGMA:
			g_string_append(str, " with sigma clipping rejection");
			break;
		case MAD:
			g_string_append(str, " with MAD clipping rejection");
			break;
		case SIGMEDIAN:
			g_string_append(str, " with median sigma clipping rejection");
			break;
		case WINSORIZED:
			g_string_append(str, " with winsorized sigma clipping rejection");
			break;
		case LINEARFIT:
			g_string_append(str, " with linear fit clipping rejection");
			break;
		case GESDT:
			g_string_append(str, " with GESDT rejection");
			break;
		}
		if (args->type_of_rejection != NO_REJEC) {
			if (args->type_of_rejection == GESDT) {
				g_string_append_printf(str, " (outliers=%.3f significance=%.3f)",
						args->sig[0], args->sig[1]);
			} else {
				g_string_append_printf(str, " (low=%.3f high=%.3f)",
						args->sig[0], args->sig[1]);
			}
		}
	}

	/* Normalisation */
	if (args->method != &stack_mean_with_rejection &&
			args->method != &stack_median ) {
		g_string_append(str, ", unnormalized input");
	} else {
		switch (args->normalize) {
		default:
		case NO_NORM:
			g_string_append(str, ", unnormalized input");
			break;
		case ADDITIVE:
			g_string_append(str, ", additive normalized input");
			break;
		case MULTIPLICATIVE:
			g_string_append(str, ", multiplicative normalized input");
			break;
		case ADDITIVE_SCALING:
			g_string_append(str, ", additive+scaling normalized input");
			break;
		case MULTIPLICATIVE_SCALING:
			g_string_append(str, ", multiplicative+scaling normalized input");
			break;
		}
		if (args->normalize != NO_NORM && args->lite_norm)
			g_string_append(str, " (fast)");
		if (args->normalize != NO_NORM && args->overlap_norm)
			g_string_append(str, " (overlaps)");
	}

	if (args->output_norm)
		g_string_append(str, ", normalized output");
	else g_string_append(str, ", unnormalized output");

	if (args->weighting_type == NOISE_WEIGHT)
		g_string_append(str, ", image weighting from noise");
	else if (args->weighting_type == NBSTACK_WEIGHT)
		g_string_append(str, ", image weighting from image count");
	else if (args->weighting_type == WFWHM_WEIGHT)
		g_string_append(str, ", image weighting from weighted FWHM");
	else if (args->weighting_type == NBSTARS_WEIGHT)
		g_string_append(str, ", image weighting from star count");
	else g_string_append(str, ", no image weighting");

	if (args->feather_dist > 0)
		g_string_append_printf(str, ", feather distance %d pixels", args->feather_dist);

	if (args->seq->nb_layers > 1) {
		if (args->equalizeRGB)
			g_string_append(str, ", equalized RGB");
		else  g_string_append(str, ", unequalized RGB");
	}

	*hist = g_slist_append(*hist, g_string_free(str, FALSE));
}

void clean_end_stacking(struct stacking_args *args) {
	if (!args->retval)
		_show_summary(args);
	remove_tmp_upscaled_files(args);
}

/* because this idle function is called after one of many stacking method
 * functions, it contains all generic wrap-up stuff instead of only graphical
 * operations. */
static gboolean end_stacking(gpointer p) {
	struct timeval t_end;
	struct stacking_args *args = (struct stacking_args *)p;
	siril_debug_print("Ending stacking idle function, retval=%d\n", args->retval);
	stop_processing_thread();

	if (args->retval == ST_OK) {
		/* copy result to gfit if success */
		clearfits(&gfit);
		memcpy(&gfit, &args->result, sizeof(fits));
		if (!com.script)
			icc_auto_assign(&gfit, ICC_ASSIGN_ON_STACK);
		clear_stars_list(TRUE);
		/* check in com.seq, because args->seq may have been replaced */
		if (args->upscale_at_stacking)
			com.seq.current = SCALED_IMAGE;
		else com.seq.current = RESULT_IMAGE;
		/* Warning: the previous com.uniq is not freed, but calling
		 * close_single_image() will close everything before reopening it,
		 * which is quite slow */
		com.uniq = calloc(1, sizeof(single));
		com.uniq->comment = strdup(_("Stacking result image"));
		com.uniq->nb_layers = gfit.naxes[2];
		com.uniq->fit = &gfit;
		/* Giving summary if average rejection stacking */
		_show_summary(args);
		/* Giving noise estimation (new thread) */
		bgnoise_async(&gfit, TRUE);

		// updating the header string to parse the final name
		// and parse the name
		int status = PATHPARSE_ERR_OK;
		gchar *expression = g_strdup(args->output_filename);
		gchar *parsedname = update_header_and_parse(&gfit, expression, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);

		if (!parsedname || parsedname[0] == '\0') { // we cannot handout a NULL filename
			args->output_parsed_filename = g_strdup("unknown");
		} else {
			args->output_parsed_filename = g_strdup(parsedname);
		}
		g_free(parsedname);
		g_free(expression);

		/* save stacking result */
		if (args->output_parsed_filename != NULL && args->output_parsed_filename[0] != '\0') {
			int failed = 0;
			if (g_file_test(args->output_parsed_filename, G_FILE_TEST_EXISTS)) {
				failed = !args->output_overwrite;
				if (!failed) {
					if (g_unlink(args->output_parsed_filename) == -1)
						failed = 1;
					if (!failed && savefits(args->output_parsed_filename, &gfit))
						failed = 1;
					if (!failed) {
						com.uniq->filename = strdup(args->output_parsed_filename);
						com.uniq->fileexist = TRUE;
					}
				}
			}
			else {
				// output folder (if any) was already created by update_header_and_parse
				if (!savefits(args->output_parsed_filename, &gfit)) {
					com.uniq->filename = strdup(args->output_parsed_filename);
					com.uniq->fileexist = TRUE;
				} else {
					failed = 1;
				}
			}
			gfit.keywords.filename[0] = '\0'; // clear the reference to the original filename
			if (failed) {
				com.uniq->filename = strdup(_("Unsaved stacking result"));
				com.uniq->fileexist = FALSE;
			} else {
				if (args->create_rejmaps) {
					siril_log_message(_("Saving rejection maps\n"));
					if (args->merge_lowhigh_rejmaps) {
						char new_ext[30];
						sprintf(new_ext, "_low+high_rejmap%s", com.pref.ext);
						gchar *low_filename = replace_ext(args->output_parsed_filename, new_ext);
						soper_unscaled_div_ushort_to_float(args->rejmap_low, args->nb_images_to_stack);
						describe_stack_for_history(args, &args->rejmap_low->history, TRUE, FALSE);
						savefits(low_filename, args->rejmap_low);
						g_free(low_filename);
					} else {
						char new_ext[30];
						sprintf(new_ext, "_low_rejmap%s", com.pref.ext);
						gchar *low_filename = replace_ext(args->output_parsed_filename, new_ext);
						soper_unscaled_div_ushort_to_float(args->rejmap_low, args->nb_images_to_stack);
						describe_stack_for_history(args, &args->rejmap_low->history, TRUE, TRUE);
						savefits(low_filename, args->rejmap_low);
						g_free(low_filename);

						sprintf(new_ext, "_high_rejmap%s", com.pref.ext);
						gchar *high_filename = replace_ext(args->output_parsed_filename, new_ext);
						soper_unscaled_div_ushort_to_float(args->rejmap_high, args->nb_images_to_stack);
						describe_stack_for_history(args, &args->rejmap_high->history, TRUE, FALSE);
						savefits(high_filename, args->rejmap_high);
						g_free(high_filename);
					}
				}
			}

			display_filename();
			set_precision_switch(); // set precision on screen
		}

		initialize_display_mode();

		sliders_mode_set_state(gui.sliders);
		set_cutoff_sliders_max_values();
		set_sliders_value_to_gfit();
		adjust_cutoff_from_updated_gfit();	// computes min and max

		set_display_mode();

		/* update menus */
		update_MenuItem();

		redraw(REMAP_ALL);
		redraw_previews();
		sequence_list_change_current();
		update_stack_interface(TRUE);
		bgnoise_await();
	} else {
		clearfits(&args->result);
		siril_log_color_message(_("Stacking failed, please check the log to fix your issue.\n"), "red");
		if (args->retval == ST_ALLOC_ERROR) {
			siril_log_message(_("It looks like there is a memory allocation error, change memory settings and try to fix it.\n"));
		}
	}

	/* remove tmp files if exist (Upscale) */
	remove_tmp_upscaled_files(args);

	memset(&args->result, 0, sizeof(fits));
	if (args->create_rejmaps) {
		clearfits(args->rejmap_low);
		free(args->rejmap_low);
		if (!args->merge_lowhigh_rejmaps) {
			clearfits(args->rejmap_high);
			free(args->rejmap_high);
		}
	}
	set_cursor_waiting(FALSE);
	/* Do not display time for stack_summing_generic
	 * cause it uses the generic function that already
	 * displays the time
	 */
	if (args->method != &stack_summing_generic) {
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}
	stacking_args_deep_free(args);
	return FALSE;
}

int find_refimage_in_indices(const int *indices, int nb, int ref) {
	int i;
	for (i = 0; i < nb; i++) {
		if (indices[i] == ref)
			return i;
	}
	return -1;
}

void free_list_date(gpointer data) {
	DateEvent *item = (DateEvent *)data;
	g_date_time_unref(item->date_obs);
	g_slice_free(DateEvent, item);
}

DateEvent* new_date_item(GDateTime *dt, gdouble exposure) {
	DateEvent *item = g_slice_new(DateEvent);
	item->exposure = exposure;
	item->date_obs = dt;
	return item;
}

static gint list_date_compare(gconstpointer *a, gconstpointer *b) {
	const DateEvent *dt1 = (const DateEvent *) a;
	const DateEvent *dt2 = (const DateEvent *) b;

	return g_date_time_compare(dt1->date_obs, dt2->date_obs);
}

/* computes the observation date (beginning) and the start of exposure (same
 * but in julian date) and end of exposure (start of last shot + exposure time)
 */
void compute_date_time_keywords(GList *list_date, fits *fit) {
	if (!list_date)
		return;
	GDateTime *date_obs;
	gdouble start, end;
	/* First we want to sort the list */
	list_date = g_list_sort(list_date, (GCompareFunc) list_date_compare);

	/* go to the first stacked image and get needed values */
	list_date = g_list_first(list_date);
	date_obs = g_date_time_ref(((DateEvent *)list_date->data)->date_obs);
	start = date_time_to_Julian(((DateEvent *)list_date->data)->date_obs);

	/* go to the last stacked image and get needed values
	 * This time we need to add the exposure to the date_obs
	 * to exactly retrieve the end of the exposure
	 */
	list_date = g_list_last(list_date);
	gdouble last_exp = ((DateEvent *)list_date->data)->exposure;
	GDateTime *last_date = ((DateEvent *)list_date->data)->date_obs;
	GDateTime *corrected_last_date = g_date_time_add_seconds(last_date, (gdouble) last_exp);

	end = date_time_to_Julian(corrected_last_date);

	g_date_time_unref(corrected_last_date);

	/* we address the computed values to the keywords */
	fit->keywords.date_obs = date_obs;
	fit->keywords.expstart = start;
	fit->keywords.expend = end;
}

void init_stacking_args(struct stacking_args *args) {
	args->method = 0;
	args->seq = NULL;
	args->ref_image = -1;
	args->filtering_criterion = NULL;
	args->filtering_parameter = -1.;
	args->nb_images_to_stack = 0;
	args->image_indices = NULL;
	args->description = NULL;
	args->output_filename = NULL;
	args->output_parsed_filename = NULL;
	args->output_overwrite = FALSE;

	args->lite_norm = FALSE;
	args->normalize = NO_NORM;
	args->coeff = (norm_coeff){ 0 };
	args->force_norm = FALSE;
	args->output_norm = FALSE;
	args->use_32bit_output = FALSE;
	args->feather_dist = 0;
	args->overlap_norm = FALSE;
	args->reglayer = -1;
	args->equalizeRGB = FALSE;
	args->maximize_framing = FALSE;
	memset(args->offset, 0, 2 * sizeof(int));
	args->upscale_at_stacking = FALSE;

	args->type_of_rejection = NO_REJEC;
	memset(args->sig, 0, 2 * sizeof(float));
	args->critical_value = NULL;
	args->create_rejmaps = FALSE;
	args->merge_lowhigh_rejmaps = FALSE;
	args->rejmap_low = NULL;
	args->rejmap_high = NULL;

	args->weighting_type = NO_WEIGHT;
	args->weights = NULL;

	args->sd_calculator = NULL;
	args->mad_calculator = NULL;

	args->t_start = (struct timeval){ 0 };
	args->retval = 0;
	args->result = (fits){ 0 };
}

void stacking_args_deep_copy(struct stacking_args *from, struct stacking_args *to) {
	memcpy(to, from, sizeof(struct stacking_args));
	// sequence is not duplicated
	to->image_indices = malloc(from->nb_images_to_stack * sizeof(int));
	memcpy(to->image_indices, from->image_indices, from->nb_images_to_stack * sizeof(int));
	to->description = g_strdup(from->description);
	// output_filename is not duplicated, can be changed until the last minute
}

void stacking_args_deep_free(struct stacking_args *args) {
	free(args->image_indices);
	free(args->description);
	free(args->critical_value);
	free(args);
}

// used by sum, min and max stacking
// needs to take into account scale because those methods are called on the upscaled sequence
void compute_max_framing(struct stacking_args *args, int output_size[2], int offset[2]) {
	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;
	int nb_frames = args->nb_images_to_stack;
	double scale = (args->upscale_at_stacking) ? 2. : 1.;
	for (int i = 0; i < nb_frames; ++i) {
		int image_index = args->image_indices[i]; // image index in sequence
		regdata *regdat = args->seq->regparam[args->reglayer];
		int rx = (args->seq->is_variable) ? args->seq->imgparam[image_index].rx : args->seq->rx;
		int ry = (args->seq->is_variable) ? args->seq->imgparam[image_index].ry : args->seq->ry;
		xmin = (xmin > regdat[image_index].H.h02 * scale) ? regdat[image_index].H.h02 * scale : xmin;
		ymin = (ymin > regdat[image_index].H.h12 * scale) ? regdat[image_index].H.h12 * scale : ymin;
		xmax = (xmax < regdat[image_index].H.h02 * scale + rx) ? regdat[image_index].H.h02 * scale + rx : xmax;
		ymax = (ymax < regdat[image_index].H.h12 * scale + ry) ? regdat[image_index].H.h12 * scale + ry : ymax;
	}
	// using same formulas as in applyreg::compute_roi
	output_size[0] = (int)xmax - (int)xmin + 1;
	output_size[1] = (int)ymax - (int)ymin + 1;
	offset[0] = (int)xmin;
	offset[1] = -(int)ymax; // the stack is done with origin at bottom left but the shifts are computed from top right
	siril_debug_print("new size: %d %d\n", output_size[0], output_size[1]);
	siril_debug_print("new origin: %d %d\n", offset[0], offset[1]);
}
