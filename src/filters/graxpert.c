/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

// Comment out this #define before public release
//#define GRAXPERT_DEBUG

#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#include <gio/gwin32inputstream.h>
#else
#include <sys/types.h> // for waitpid(2)
#include <sys/wait.h> // for waitpid(2)
#include <gio/gunixinputstream.h>
#endif
#include <json-glib/json-glib.h>
#include "algos/background_extraction.h"
#include "core/siril.h"
#include "core/icc_profile.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_update.h"
#include "core/siril_log.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/siril_preview.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "filters/graxpert.h"

static gboolean verbose = TRUE;
static version_number graxpert_version = { 0 };

static void child_watch_cb(GPid pid, gint status, gpointer user_data) {
	siril_debug_print("GraXpert exited with status %d\n", status);
	g_spawn_close_pid(pid);
	// GraXpert has exited, reset the stored pid
	com.child_is_running = EXT_NONE;
#ifdef _WIN32
	com.childhandle = 0;		// For Windows, handle of a child process
#else
	com.childpid = 0;			// For other OSes, PID of a child process
#endif
}

static int exec_prog_graxpert(char **argv, gboolean graxpert_no_exit_report) {
	const gchar *progress_key = "Progress: ";
	gint child_stderr;
	GPid child_pid;
	g_autoptr(GError) error = NULL;
	int retval = -1;

	int index = 0;
	while (argv[index]) {
		fprintf(stdout, "%s ", argv[index++]);
	}
	fprintf(stdout, "\n");
	// g_spawn handles wchar so not need to convert
	set_progress_bar_data(_("Starting GraXpert..."), 0.0);
	g_spawn_async_with_pipes(NULL, argv, NULL,
			G_SPAWN_SEARCH_PATH |
			G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_DO_NOT_REAP_CHILD,
			NULL, NULL, &child_pid, NULL, NULL,
			&child_stderr, &error);

	int i = 0;
	while (argv[i])
		g_free(argv[i++]);

	if (error != NULL) {
		siril_log_color_message(_("Spawning GraXpert failed: %s\n"), "red", error->message);
		return retval;
	}
	g_child_watch_add(child_pid, child_watch_cb, NULL);
	com.child_is_running = EXT_GRAXPERT;
#ifdef _WIN32
	com.childhandle = child_pid;		// For Windows, handle of a child process
#else
	com.childpid = child_pid;			// For other OSes, PID of a child process
#endif

	GInputStream *stream = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(child_stderr), FALSE);
#else
	stream = g_unix_input_stream_new(child_stderr, FALSE);
#endif

	gboolean doprint = TRUE;
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
#ifdef GRAXPERT_DEBUG
	gchar *lastbuffer = NULL;
#endif
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL))) {
		gchar *arg = g_strstr_len(buffer, -1, progress_key);
		double value = -1.0;
		if (arg)
			value = g_ascii_strtod(arg + strlen(progress_key), NULL);
		if (value > 0.0 && value == value && verbose) {
			set_progress_bar_data(_("Running GraXpert"), value / 100.0);
		} else if (g_strrstr(buffer, "Finished")) {
			set_progress_bar_data(_("Done."), 1.0);
			retval = 0;
#ifdef GRAXPERT_DEBUG
			if (verbose)
				siril_log_message(_("GraXpert caught buffer: %s : exit successful\n"), buffer);
#endif
		} else if (doprint && verbose) {
			if (!g_strstr_len(buffer, -1, "ForkProcess") && !g_strstr_len(buffer, -1, "ai_version")) { // These are not useful messages
				gchar *print_from = g_strstr_len(buffer, -1, "INFO");
				if (print_from) {
					if (strlen(print_from) > 9) {
						print_from += 9;
						siril_log_message("GraXpert: %s\n", print_from);
					}
				}
			}
		}
#ifdef GRAXPERT_DEBUG
		lastbuffer = g_strdup(buffer);
#endif
		g_free(buffer);
	}
	// GraXpert has exited, reset the stored pid
#ifdef _WIN32
	com.childhandle = 0;		// For Windows, handle of a child process
#else
	com.childpid = 0;			// For other OSes, PID of a child process
#endif
	if (graxpert_no_exit_report) {
		siril_log_message(_("GraXpert GUI finished.\n"));
		set_progress_bar_data(_("Done."), 1.0);
		retval = 0; // No "Finished" log message is printed when closing the GUI, so we assume success
	}
#ifdef GRAXPERT_DEBUG
	if (retval)
		siril_log_message(_("GraXpert exit not caught, last buffer: %s\n"), lastbuffer);
	g_free(lastbuffer);
#endif
	g_object_unref(data_input);
	g_object_unref(stream);
	if (!g_close(child_stderr, &error))
		siril_debug_print("%s\n", error->message);
	return retval;
}

gboolean graxpert_executablecheck(gchar* executable, graxpert_operation operation) {
	const gchar *version_key = "version: ";
	char *test_argv[3] = { NULL };
	gint child_stderr;
	g_autoptr(GError) error = NULL;
	version_number null_version = { 0 };
	if (!memcmp(&graxpert_version, &null_version, sizeof(version_number))) {
		if (!executable || executable[0] == '\0') {
			return FALSE;
		}
		if (!g_file_test(executable, G_FILE_TEST_IS_EXECUTABLE)) {
			return FALSE; // It's not executable so return a zero version number
		}

		int nb = 0;
		test_argv[nb++] = executable;
		gchar *versionarg = g_strdup("-v");
		test_argv[nb++] = versionarg;
		// g_spawn handles wchar so not need to convert
		g_spawn_async_with_pipes(NULL, test_argv, NULL,
				G_SPAWN_SEARCH_PATH |
				G_SPAWN_LEAVE_DESCRIPTORS_OPEN,
				NULL, NULL, NULL, NULL, NULL,
				&child_stderr, &error);

		if (error != NULL) {
			siril_log_color_message(_("Spawning GraXpert failed during version check: %s\n"), "red", error->message);
			g_free(versionarg);
			return FALSE;
		}

		GInputStream *stream = NULL;
#ifdef _WIN32
		stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(child_stderr), FALSE);
#else
		stream = g_unix_input_stream_new(child_stderr, FALSE);
#endif
		gchar *buffer;
		gsize length = 0;
		GDataInputStream *data_input = g_data_input_stream_new(stream);
		gboolean done = FALSE;
		while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
						NULL, NULL)) && !done) {
			// Find the start of the version substring
			gchar *version_start = g_strstr_len(buffer, -1, version_key);
			if (version_start) {
				version_start += strlen(version_key); // Move past "version: "

				// Find the end of the version substring
				gchar *version_end = g_strstr_len(version_start, -1, " ");
				*version_end = '\0';

				// Extract the version substring
				graxpert_version = get_version_number_from_string(version_start);

				g_free(buffer);
				break;
			}
			g_free(buffer);
		}
		g_object_unref(data_input);
		g_object_unref(stream);
		g_free(versionarg);
		if (!g_close(child_stderr, &error)) {
			siril_debug_print("%s\n", error->message);
		}
	}

	if (compare_version(graxpert_version, (version_number) {.major_version = 3, .minor_version = 0, .micro_version = 0}) < 0) {
		siril_log_color_message(_("Error: GraXpert version is too old. You have version %d.%d.%d; at least version 3.0.0 is required.\n"), "red", graxpert_version.major_version, graxpert_version.minor_version, graxpert_version.micro_version);
		return FALSE;
	} else {
		return TRUE;
	}
}

gboolean save_graxpert_config(graxpert_data *args) {
	const gchar *filename = args->configfile;
	JsonBuilder *builder;
	JsonGenerator *generator;
	JsonNode *root;
	gchar *json_data;
	gsize json_length;
	gboolean success = FALSE;

	builder = json_builder_new();
	json_builder_begin_object(builder);

	// Add working directory
	json_builder_set_member_name(builder, "working_dir");
	json_builder_add_string_value(builder, com.wd);

	// Add width and height
	json_builder_set_member_name(builder, "width");
	json_builder_add_int_value(builder, args->fit->rx);
	json_builder_set_member_name(builder, "height");
	json_builder_add_int_value(builder, args->fit->ry);

	// Add background points
	if (args->bg_samples) {
		json_builder_set_member_name(builder, "background_points");
		json_builder_begin_array(builder);
		for (GSList *l = args->bg_samples; l != NULL; l = l->next) {
			background_sample *s = (background_sample *)l->data;
			json_builder_begin_array(builder);
			json_builder_add_int_value(builder, min(args->fit->rx - 1, round_to_int(s->position.x)));
			json_builder_add_int_value(builder, min(args->fit->ry - 1, round_to_int(s->position.y)));
			json_builder_add_int_value(builder, 1);
			json_builder_end_array(builder);
		}
		json_builder_end_array(builder);
	}

	// Add other options
	json_builder_set_member_name(builder, "bg_pts_option");
	json_builder_add_int_value(builder, args->bg_pts_option);

	json_builder_set_member_name(builder, "stretch_option");
	json_builder_add_string_value(builder,
		args->stretch_option == STRETCH_OPTION_10_BG_3_SIGMA ? "10% Bg, 3 sigma" :
		args->stretch_option == STRETCH_OPTION_15_BG_3_SIGMA ? "15% Bg, 3 sigma":
		args->stretch_option == STRETCH_OPTION_20_BG_3_SIGMA ? "20% Bg, 3 sigma":
		args->stretch_option == STRETCH_OPTION_30_BG_2_SIGMA ? "30% Bg, 2 sigma":
		"No Stretch");

	json_builder_set_member_name(builder, "bg_tol_option");
	json_builder_add_double_value(builder,
		args->bg_tol_option < -2.0 ? -2.0 : args->bg_tol_option > 10.0 ? 10.0 : args->bg_tol_option);

	json_builder_set_member_name(builder, "interpol_type_option");
	json_builder_add_string_value(builder, args->bg_algo == GRAXPERT_BG_KRIGING ? "Kriging" :
		args->bg_algo == GRAXPERT_BG_RBF ? "RBF" :
		args->bg_algo == GRAXPERT_BG_SPLINE ? "Splines" :
		"AI");

	json_builder_set_member_name(builder, "smoothing_option");
	json_builder_add_double_value(builder, args->bg_smoothing);

	if (args->ai_batch_size != -1) {
		json_builder_set_member_name(builder, "ai_batch_size");
		json_builder_add_int_value(builder, args->ai_batch_size);
	}

	if (args->sample_size != -1) {
		json_builder_set_member_name(builder, "sample_size");
		json_builder_add_int_value(builder, args->sample_size);
	}

	if (args->spline_order != -1) {
		json_builder_set_member_name(builder, "spline_order");
		json_builder_add_int_value(builder, args->spline_order);
	}

	json_builder_set_member_name(builder, "corr_type");
	json_builder_add_string_value(builder, args->bg_mode == GRAXPERT_SUBTRACTION ? "Subtraction" :
	"Division");

	json_builder_set_member_name(builder, "RBF_kernel");
	json_builder_add_string_value(builder, args->kernel == GRAXPERT_THIN_PLATE ? "thin_plate" :
	args->kernel == GRAXPERT_QUINTIC ? "quintic" :
	args->kernel == GRAXPERT_CUBIC ? "cubic" :
	"linear");

	json_builder_set_member_name(builder, "denoise_strength");
	json_builder_add_double_value(builder, args->denoise_strength);

	json_builder_set_member_name(builder, "saveas_option");
	json_builder_add_string_value(builder,
		args->fit->type == DATA_FLOAT ? "32 bit Fits" : "16 bit Fits");

	json_builder_end_object(builder);

	root = json_builder_get_root(builder);
	generator = json_generator_new();
	json_generator_set_root(generator, root);
	json_generator_set_pretty(generator, TRUE);

	json_data = json_generator_to_data(generator, &json_length);

	if (g_file_set_contents(filename, json_data, json_length, NULL)) {
		success = TRUE;
	}

	g_free(json_data);
	json_node_free(root);
	g_object_unref(generator);
	g_object_unref(builder);

	return success;
}

graxpert_data *new_graxpert_data() {
	graxpert_data *p = calloc(1, sizeof(graxpert_data));
	p->fit = &gfit;
	p->operation = GRAXPERT_BG;
	p->bg_smoothing = 0.5;
	p->bg_algo = GRAXPERT_BG_AI;
	p->bg_mode = GRAXPERT_SUBTRACTION;
	p->kernel = GRAXPERT_THIN_PLATE;
	p->stretch_option = STRETCH_OPTION_10_BG_3_SIGMA;
	p->sample_size = 25;
	p->spline_order = 3;
	p->bg_tol_option = 2;
	p->denoise_strength = 0.8;
	p->use_gpu = TRUE;
	p->ai_batch_size = 4;
	p->bg_pts_option = 15;
	return p;
}

void free_graxpert_data(graxpert_data *args) {
	g_free(args->path);
	g_free(args->configfile);
	free_background_sample_list(args->bg_samples);
	if (args->backup_icc)
		cmsCloseProfile(args->backup_icc);
	free(args);
}

static void open_graxpert_result(graxpert_data *args) {
	// Clean up config file if one was used
	if (args->configfile && g_unlink(args->configfile))
		siril_debug_print("Failed to remove GraXpert config file\n");

	// If successful, open the result image
	/* Note: we do this even when sequence working: it's a bit inefficient to read it in,
	 * delete it and save it again but it works for all sequence types (incl. ser and
	 * FITSEQ) whereas simply moving the file saved by GraXpert would only work for
	 * sequences of individual FITS files. */
	if (args->path) {
		disable_profile_check_verbose();
		if (args->previewing) {
			siril_log_message(_("Reading result from the temporary working file...\n"));
			if (readfits(args->path, &gui.roi.fit, NULL, (gui.roi.fit.type == DATA_FLOAT))) {
				siril_log_color_message(_("Error opening GraXpert result. Check the file %s\n"), "red", args->path);
				enable_profile_check_verbose();
				goto END_AND_RETURN;
			}
		} else {
			fits *result = calloc(1, sizeof(fits));
			siril_log_message(_("Reading result from the temporary working file...\n"));
			if (readfits(args->path, result, NULL, !com.pref.force_16bit)) {
				siril_log_color_message(_("Error opening GraXpert result. Check the file %s\n"), "red", args->path);
				enable_profile_check_verbose();
				goto END_AND_RETURN;
			}
			// Check the result dimensions and bit depth match what we expect
			// This should never fail, but we shouldn't trust external software too
			// much and a dimension mismatch would result in a crash
			if (args->fit->rx != result->rx || args->fit->ry != result->ry || args->fit->naxes[2] != result->naxes[2]
					|| args->fit->type != result->type) {
				siril_log_color_message(_("Error: the GraXpert image dimensions and bit depth "
						"do not match those of the original image. Please report this as a bug.\n"), "red");
				goto END_AND_RETURN;
			}
			// Swap the image data pointers. This way the data originally in args->fit
			// gets freed by the call to clearfits(result) and the result data becomes
			// owned by args->fit
			if (fits_swap_image_data(args->fit, result)) {
				siril_debug_print("Error, NULL pointer passed to fits_swap_image_data\n");
				goto END_AND_RETURN;
			}
			clearfits(result);
			free(result);
			if (args->fit == &gfit) {
				copy_gfit_to_backup();
				populate_roi();
			}
		}
		siril_log_message(_("Removing the temporary working file...\n"));
		if (g_unlink(args->path))
			siril_debug_print("Failed to unlink GraXpert working file\n");
		enable_profile_check_verbose();
	}
END_AND_RETURN:
	unlock_roi_mutex();
	return;
}

static gboolean end_graxpert(gpointer p) {
	stop_processing_thread();
	open_graxpert_result((graxpert_data *) p);
	free_graxpert_data((graxpert_data *) p);
	notify_gfit_modified();
	launch_clipboard_survey();
	return end_generic(NULL);
}

gpointer do_graxpert (gpointer p) {
	lock_roi_mutex();
	graxpert_data *args = (graxpert_data *) p;
	if (args->fit == &gfit)
		copy_backup_to_gfit();
	if (!args->previewing && !com.script) {
		gchar *text = NULL;
		switch (args->operation) {
			case GRAXPERT_BG:
				text = g_strdup_printf(_("GraXpert BG extraction, smoothness %.3f"), args->bg_smoothing);
				break;
			case GRAXPERT_DENOISE:
				text = g_strdup_printf(_("GraXpert denoising, strength %.3f"), args->denoise_strength);
				break;
			default:
				text = g_strdup(_("GraXpert operations using GUI"));
		}
		undo_save_state(&gfit, text);
		g_free(text);
	}
	char *my_argv[64] = { 0 };
	gchar *filename = NULL, *path = NULL, *outpath = NULL;
	int retval = 1;
	if (!graxpert_executablecheck(com.pref.graxpert_path, args->operation))
		goto ERROR_OR_FINISHED;

	// Configure input filename
	if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
		filename = g_path_get_basename(com.uniq->filename);
		gchar *temp = remove_ext_from_filename(filename);
		g_free(filename);
		filename = g_strdup_printf("graxpert_%s", temp);
		g_free(temp);
	} else {
		filename = g_strdup_printf("graxpert");
	}

	// We have to forcibly disable FITS compression as GraXpert cannot open compressed FITS (true as of GraXpert 3.0.2)
	// We also have to force the extension ".fits" to avoid GraXpert adding it onto other extensions
	gboolean pref_fitscomp = com.pref.comp.fits_enabled;
	gchar *backup_ext = g_strdup(com.pref.ext);
	g_free(com.pref.ext);
	com.pref.ext = g_strdup(".fits");
	com.pref.comp.fits_enabled = FALSE;
	gchar *temp = set_right_extension(filename);
	g_free(filename);
	filename = temp;
	g_free(com.pref.ext);
	com.pref.ext = backup_ext;

	path = g_build_filename(com.wd, filename, NULL);
	temp = g_strdup(path);
	outpath = remove_ext_from_filename(temp);
	g_free(temp);
	g_free(filename);
	// Save current image as input filename
	siril_log_message(_("Saving temporary working file...\n"));
	if (savefits(path, args->fit)) {
		siril_log_color_message(_("Error: failed to save temporary FITS\n"), "red");
		goto ERROR_OR_FINISHED;
	}

	com.pref.comp.fits_enabled = pref_fitscomp;

	// Configure GraXpert commandline
	int nb = 0;
	gboolean graxpert_no_exit_report = FALSE;
	gboolean is_gui = FALSE;
	my_argv[nb++] = g_strdup(com.pref.graxpert_path);
	if (args->operation == GRAXPERT_BG) {
		my_argv[nb++] = g_strdup("-cli");
		my_argv[nb++] = g_strdup("-cmd");
		my_argv[nb++] = g_strdup("background-extraction");
		if (args->bg_algo == GRAXPERT_BG_AI) {
			my_argv[nb++] = g_strdup("-correction");
			my_argv[nb++] = g_strdup_printf("%s", args->bg_mode == GRAXPERT_SUBTRACTION ? "Subtraction" : "Division");
			my_argv[nb++] = g_strdup("-smoothing");
			my_argv[nb++] = g_strdup_printf("%f", args->bg_smoothing);
			my_argv[nb++] = g_strdup("-gpu");
			my_argv[nb++] = g_strdup(args->use_gpu ? "true" : "false");
		} else {
			if (!args->bg_samples) {
				siril_log_color_message(_("Background samples must be computed for GraXpert RBF, Spline and Kriging methods\n"), "red");
				goto ERROR_OR_FINISHED;
			}
			my_argv[nb++] = g_strdup("-preferences_file");
			args->configfile = g_build_filename(com.wd, "siril-graxpert.pref", NULL);
			my_argv[nb++] = args->configfile;
			save_graxpert_config(args);
		}
		if (args->keep_bg)
			my_argv[nb++] = g_strdup("-bg");
		my_argv[nb++] = g_strdup("-output");
		my_argv[nb++] = g_strdup_printf("%s", outpath);
		graxpert_no_exit_report = TRUE;
	} else if (args->operation == GRAXPERT_DENOISE) {
		my_argv[nb++] = g_strdup("-cli");
		my_argv[nb++] = g_strdup("-cmd");
		my_argv[nb++] = g_strdup("denoising");
		my_argv[nb++] = g_strdup_printf("-output");
		my_argv[nb++] = g_strdup_printf("%s", outpath);
		my_argv[nb++] = g_strdup("-gpu");
		my_argv[nb++] = g_strdup(args->use_gpu ? "true" : "false");
		my_argv[nb++] = g_strdup("-strength");
		my_argv[nb++] = g_strdup_printf("%.2f", args->denoise_strength);
	} else if (args->operation == GRAXPERT_GUI) {
		siril_log_message(_("GraXpert GUI will open with the current image. When you have finished, save the "
							"image and return to Siril.\nYou will need to check and re-assign the ICC profile "
							"as GraXpert does not preserve ICC profiles embedded in FITS files.\n"));
		graxpert_no_exit_report = TRUE;
		is_gui = TRUE;
	} else {
		siril_log_message(_("Error: unknown GraXpert operation\n"));
		g_free(my_argv[0]);
		goto ERROR_OR_FINISHED;
	}
	my_argv[nb++] = g_strdup_printf("%s", path);

	// Save a copy of the current ICC profile, as GraXpert does not preserve these
	if (args->fit->icc_profile)
		args->backup_icc = copyICCProfile(args->fit->icc_profile);

	// Execute GraXpert
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	retval = exec_prog_graxpert(my_argv, graxpert_no_exit_report);
	gettimeofday(&t_end, NULL);
	if (verbose)
		show_time(t_start, t_end);
	if (retval) {
		if (g_unlink(path))
			siril_debug_print("Error removing temporary file %s\n", path);
		goto ERROR_OR_FINISHED;
	}

	if (!is_gui) {
		if (args->path)
			g_free(args->path);
		args->path = g_strdup(path);
	}
	g_free(path);

ERROR_OR_FINISHED:
	g_free(outpath);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	if (!args->seq && !com.script)
		siril_add_idle(end_graxpert, args); // this loads the result
	else
		open_graxpert_result(args);
	return GINT_TO_POINTER(retval);
}

/******** SEQUENCE FUNCTIONS *********/

static int graxpert_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
//	GraXpert cannot run in parallel as it fully utilizes the GPU / CPU. This function therefore
//	returns a maximum of 1 and all images will be processed in series.
	unsigned int MB_per_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	if (limit > 0) {
		required = MB_per_image;
		limit = MB_avail / required;
	}
	limit = (limit >= 1 ? 1 : 0);
	return limit;
}

static int graxpert_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 0;
	graxpert_data *data = (graxpert_data *) args->user;
	if (data->operation == GRAXPERT_BG && data->bg_algo != GRAXPERT_BG_AI) {
		const char *err;
		data->bg_samples = generate_samples(fit, data->bg_pts_option, data->bg_tol_option, data->sample_size, &err, MULTI_THREADED);
		if (!data->bg_samples) {
			siril_log_color_message(_("Failed to generate background samples for image: %s\n"), "red", _(err));
			return 1;
		}
	}
	data->fit = fit;
	siril_log_color_message(_("GraXpert: Processing image %d\n"), "green", o + 1);
	verbose = FALSE;
	do_graxpert(data);
	verbose = TRUE;
	free_background_sample_list(data->bg_samples);
	return ret;
}

void apply_graxpert_to_sequence(graxpert_data *args) {
	if (!graxpert_executablecheck(com.pref.graxpert_path, args->operation)) {
		return;
	}
	args->fit = NULL;
	struct generic_seq_args *seqargs = create_default_seqargs(args->seq);
	seqargs->seq = args->seq;
	seqargs->parallel = FALSE;
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = args->seq->selnum;
	seqargs->compute_mem_limits_hook = graxpert_compute_mem_limits;
	seqargs->prepare_hook = seq_prepare_hook;
	seqargs->image_hook = graxpert_image_hook;
	seqargs->description = _("GraXpert");
	seqargs->has_output = TRUE;
	seqargs->output_type = get_data_type(seqargs->seq->bitpix);
	if (args->operation == GRAXPERT_BG)
		seqargs->new_seq_prefix = strdup("gxbg_");
	else if (args->operation == GRAXPERT_DENOISE)
		seqargs->new_seq_prefix = strdup("gxnr_");
	else  {
		free_graxpert_data(args);
		free(seqargs);
		return;
	}
	seqargs->load_new_sequence = TRUE;
	seqargs->user = args;
	set_progress_bar_data(_("GraXpert: Processing..."), 0.);
	start_in_new_thread(generic_sequence_worker, seqargs);
}
