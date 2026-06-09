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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "core/siril.h"
#include "core/gui_iface.h"
#include "algos/statistics.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
/* gui_calls.h removed: no direct calls remain */
#include "core/processing.h"
#include "core/command_list.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/ser.h"
#include "livestacking/livestacking.h"

#include "command.h"
#include "command_line_processor.h"

const char *cmd_err_to_str(cmd_errors err) {
	switch (err) {
		case CMD_NOT_FOUND:
			return _("unknown command name");
		case CMD_NO_WAIT:
			return _("non-blocking command, no error");
		case CMD_NO_CWD:
			return _("current working directory is invalid");
		case CMD_NOT_SCRIPTABLE:
			return _("command cannot be used in a script");
		case CMD_WRONG_N_ARG:
			return _("invalid number of arguments");
		case CMD_ARG_ERROR:
			return _("invalid arguments");
		case CMD_SELECTION_ERROR:
			return _("invalid image area selection");
		case CMD_OK:
			return _("command succeeded");
		default:
		case CMD_GENERIC_ERROR:
			return _("generic error");
		case CMD_IMAGE_NOT_FOUND:
			return _("input image not found");
		case CMD_SEQUENCE_NOT_FOUND:
			return _("invalid input sequence");
		case CMD_INVALID_IMAGE:
			return _("invalid input image");
		case CMD_LOAD_IMAGE_FIRST:
			return _("load an image or sequence first");
		case CMD_ONLY_SINGLE_IMAGE:
			return _("command is only for a loaded single image");
		case CMD_NOT_FOR_SINGLE:
			return _("command is only for a sequence");
		case CMD_NOT_FOR_MONO:
			return _("command is not for monochrome images");
		case CMD_NOT_FOR_RGB:
			return _("command is not for color images");
		case CMD_FOR_CFA_IMAGE:
			return _("command is only for OSC CFA images");
		case CMD_FILE_NOT_FOUND:
			return _("file not found");
		case CMD_FOR_PLATE_SOLVED:
			return _("command is only for plate solved images");
		case CMD_NEED_INIT_FIRST:
			return _("command requires a preliminary step");
		case CMD_ALLOC_ERROR:
			return _("memory allocation error");
		case CMD_THREAD_RUNNING:
			return _("a processing is already running");
        case CMD_DIR_NOT_FOUND:
			return _("directory not found");
	}
}

void parse_line(char *myline, int len, int *nb) {
	int i = 0, wordnb = 0;
	char string_starter = '\0';	// quotes don't split words on spaces
	word[0] = NULL;

	do {
		while (i < len && isblank(myline[i]))
			i++;
		if (myline[i] == '"' || myline[i] == '\'')
			string_starter = myline[i++];
		if (myline[i] == '\0' || myline[i] == '\n' || myline[i] == '\r')
			break;
		word[wordnb++] = myline + i;	// the beginning of the word
		word[wordnb] = NULL;		// put next word to NULL
		do {
			i++;
			if (string_starter != '\0' && myline[i] == string_starter) {
				string_starter = '\0';
				break;
			}
		} while (i < len && (!isblank(myline[i]) || string_starter != '\0')
				&& myline[i] != '\r' && myline[i] != '\n');
		if (myline[i] == '\0')	// the end of the word and line (i == len)
			break;
		myline[i++] = '\0';		// the end of the word
	} while (wordnb < MAX_COMMAND_WORDS - 1);
	*nb = wordnb;
}

void remove_trailing_cr(char *str) {
	if (str == NULL || str[0] == '\0')
		return;
	int length = strlen(str);
	if (str[length - 1] == '\r')
		str[length - 1] = '\0';
}

int execute_command(int wordnb) {
	// search for the command in the list
	if (word[0] == NULL) return 1;
	int i = G_N_ELEMENTS(commands);
	while (g_ascii_strcasecmp(commands[--i].name, word[0])) {
		if (i == 0) {
			siril_log_message(_("Unknown command: '%s' or not implemented yet\n"), word[0]);
			return CMD_NOT_FOUND;
		}
	}

	// verify argument count
	if (wordnb - 1 < commands[i].nbarg) {
		siril_log_message(_("Usage: %s\n"), commands[i].usage);
		return CMD_WRONG_N_ARG;
	}

	// verify if command is scriptable
	if (com.script) {
		if (!commands[i].scriptable) {
			siril_log_warning(_("This command cannot be used in a script: %s\n"), commands[i].name);
			return CMD_NOT_SCRIPTABLE;
		}
	}

	/* verify prerequires */
	/* we check that (REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE) is inside the bitmask */
	if ((commands[i].prerequires & (REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE)) == (REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE)) {
		if (!single_image_is_loaded() && !sequence_is_loaded()) {
			return CMD_LOAD_IMAGE_FIRST;
		}
	} else if ((commands[i].prerequires & REQ_CMD_SINGLE_IMAGE) == REQ_CMD_SINGLE_IMAGE) {
		if (!single_image_is_loaded()) {
			return CMD_ONLY_SINGLE_IMAGE;
		}
	} else if ((commands[i].prerequires & REQ_CMD_SEQUENCE) == REQ_CMD_SEQUENCE) {
		if (!sequence_is_loaded()) {
			return CMD_NOT_FOR_SINGLE;
		}
	}

	if ((commands[i].prerequires & REQ_CMD_FOR_CFA) == REQ_CMD_FOR_CFA) {
		if (isrgb(gfit)) {
			return CMD_FOR_CFA_IMAGE;
		}
	} else if ((commands[i].prerequires & REQ_CMD_FOR_MONO) != 0) {
		if (isrgb(gfit)) {
			return CMD_NOT_FOR_RGB;
		}
	} else if ((commands[i].prerequires & REQ_CMD_FOR_RGB) != 0) {
		if (!isrgb(gfit)) {
			return CMD_NOT_FOR_MONO;
		}
	}

	if ((commands[i].prerequires & REQ_CMD_NO_THREAD) != 0) {
		if (processing_is_job_active()) {
			return CMD_THREAD_RUNNING;
		}
	}

	// process the command
	siril_log_warning(_("Running command: %s\n"), word[0]);
	fprintf(stdout, "%lu: running command %s\n", time(NULL), word[0]);
	int retval = commands[i].process(wordnb);
	if (gui_iface.roi_is_active())
		gui_iface.populate_roi();
	if (retval & CMD_NOTIFY_GFIT_MODIFIED) {
		waiting_for_thread(); // we can't proceed until the generic_image_Worker is done
		if (!com.python_script) {
			gfit_modified_update_gui();
		} else {
			invalidate_stats_from_fit(gfit);
			gui_iface.invalidate_histogram();
			gui_iface.execute_idle_sync(end_gfit_operation, NULL);
		}
		retval = retval & ~CMD_NOTIFY_GFIT_MODIFIED;
	}
	return (int) retval & ~CMD_NOTIFY_GFIT_MODIFIED;
}

/* update_log_icon, log_status_bar_idle_callback, display_command_on_status_bar,
 * clear_status_bar, end_script moved to gui/script_console.c */

int check_requires(gboolean *checked_requires, gboolean is_required) {
	int retval = CMD_OK;
	/* check for requires command */
	if (!g_ascii_strcasecmp(word[0], "requires")) {
		*checked_requires = TRUE;
	} else if (is_required && !*checked_requires) {
		siril_log_error(_("The \"requires\" command is missing at the top of the script file."
					" This command is needed to check script compatibility.\n"));
		retval = CMD_GENERIC_ERROR;
	}
	return retval;
}

int check_command_mode() {
	/* until we have a proper implementation of modes, we just forbid other
	 * commands to be run during live stacking */
	int retval = 0;
	if (livestacking_is_started()) {
		retval = g_ascii_strcasecmp(word[0], "livestack") &&
			g_ascii_strcasecmp(word[0], "stop_ls") &&
			g_ascii_strcasecmp(word[0], "exit");
		if (retval)
			siril_log_warning(_("This command cannot be run while live stacking is active, ignoring.\n"));

	}
	return retval;
}

gpointer execute_script(gpointer p) {
	GInputStream *input_stream = (GInputStream*) p;
	gboolean checked_requires = FALSE;
	gchar *buffer;
	int line = 0, retval = 0;
	int wordnb;
	int startmem, endmem;
	struct timeval t_start, t_end;

	com.script = TRUE;
	com.stop_script = FALSE;
	com.script_thread_exited = FALSE;

	gettimeofday(&t_start, NULL);

	/* Now we want to save the cwd in order to come back after
	 * script execution
	 */
	gchar *saved_cwd = g_strdup(com.wd);
	startmem = get_available_memory() / BYTES_IN_A_MB;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
			NULL, NULL))) {
		++line;
		if (com.stop_script) {
			retval = 1;
			g_free (buffer);
			break;
		}
		/* Displays comments */
		if (buffer[0] == '#') {
			siril_log_status("%s\n", buffer);
			g_free (buffer);
			continue;
		}

		/* in Windows case, remove trailing CR */
		remove_trailing_cr(buffer);
		g_strstrip(buffer);

		if (buffer[0] == '\0') {
			g_free (buffer);
			continue;
		}

		gui_iface.console_set_status(buffer, line);
		parse_line(buffer, length, &wordnb);
		if (check_requires(&checked_requires, com.pref.script_check_requires)) {
			g_free (buffer);
			break;
		}
		if (check_command_mode()) {
			g_free (buffer);
			continue;
		};

		retval = execute_command(wordnb);
		remove_child_from_children((GPid) -2); // remove the processing thread child
		// reference (speculative - not always necessary, but simplest to
		// call it every time just in case the command ran in the thread.)
		if (retval && retval != CMD_NO_WAIT) {
			siril_log_error(_("Error in line %d ('%s'): %s.\n"), line, buffer, cmd_err_to_str(retval));
			siril_log_message(_("Exiting batch processing.\n"));
			g_free (buffer);
			break;
		}
		if (retval != CMD_NO_WAIT && waiting_for_thread()) {
			retval = 1;
			g_free (buffer);
			break;	// abort script on command failure
		}
		endmem = get_available_memory() / BYTES_IN_A_MB;
		siril_log_debug("End of command %s, memory difference: %d MB\n", word[0], startmem - endmem);
		startmem = endmem;
		memset(word, 0, sizeof word);
		g_free (buffer);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);

	if (!com.headless) {
		com.script = FALSE;
		notify_gfit_data_modified();
		gui_iface.end_script_gui();
	}

	/* Now we want to restore the saved cwd */
	siril_change_dir(saved_cwd, NULL);
	writeinitfile();
	if (!retval) {
		siril_log_message(_("Script execution finished successfully.\n"));
		gettimeofday(&t_end, NULL);
		show_time_msg(t_start, t_end, _("Total execution time"));
	} else {
		char *msg = siril_log_error(_("Script execution failed.\n"));
		msg[strlen(msg) - 1] = '\0';
		gui_iface.set_progress(PROGRESS_DONE, msg);
	}
	g_free(saved_cwd);

	if (com.script_thread) {
		siril_log_debug("Script thread exiting\n");
		com.script_thread_exited = TRUE;
	}
	/* If called from the GUI, re-enable widgets blocked during the script */
	gui_iface.script_widgets_async(TRUE);
	return GINT_TO_POINTER(retval);
}

/* show_command_help_popup, on_command_key_press_event moved to gui/script_console.c */

int processcommand(const char *line, gboolean wait_for_completion) {
	int wordnb = 0;
	GError *error = NULL;
	int ret = 0;

	if (line[0] == '\0' || line[0] == '\n')
		return CMD_NOT_FOUND;
	if (line[0] == '@') { // case of files
		if (processing_is_job_active() || (get_script_thread_run() && !com.script_thread_exited)) {
			PRINT_ANOTHER_THREAD_RUNNING;
			return CMD_THREAD_RUNNING;
		}
		if (get_script_thread_run())
			wait_for_script_thread();

		/* Switch to console tab */
		gui_iface.show_panel("output_logs", TRUE);

		char filename[256];
		g_strlcpy(filename, line + 1, 250);
		expand_home_in_filename(filename, 256);

		GFile *file = g_file_new_for_path(filename);
		GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, &error);

		if (input_stream == NULL) {
			if (error != NULL) {
				g_clear_error(&error);
				siril_log_message(_("File [%s] does not exist\n"), filename);
			}

			g_object_unref(file);
			return CMD_FILE_NOT_FOUND;
		}
		/* Run the script */
		siril_log_message(_("Starting script %s\n"), filename);
		com.script_thread = g_thread_new("script", execute_script, input_stream);
		g_object_unref(file);
	} else {
		/* Switch to console tab */
		gui_iface.show_panel("output_logs", TRUE);

		gchar *myline = strdup(line);
		int len = strlen(line);
		if (len > 0)
			g_print("input command:%s\n", myline);

		parse_line(myline, len, &wordnb);
		ret = execute_command(wordnb);

		if (ret) {
			siril_log_error(_("Command execution failed: %s.\n"), cmd_err_to_str(ret));
			if (!(com.script || com.python_script) && !com.headless &&
				(ret == CMD_WRONG_N_ARG || ret == CMD_ARG_ERROR)) {
				gui_iface.show_command_help();
			}
			free(myline);
			return ret;
		}

		if (wait_for_completion && ret != CMD_NO_WAIT) {
			remove_child_from_children((GPid) -2);   /* tidy child list first    */
			gpointer job_retval = waiting_for_thread();
			if (GPOINTER_TO_INT(job_retval))
				ret = CMD_GENERIC_ERROR;
		}
		free(myline);
	}

	return ret;
}

// loads the sequence from com.wd
sequence *load_sequence(const char *name, char **get_filename) {
	gchar *file = NULL;
	gchar *altfile = NULL;
	if (name[0] == '.' && g_utf8_strlen(name, -1) == 1 && sequence_is_loaded())
		file = g_strdup(com.seq.seqname);
	else {
		file = g_strdup(name);
		if (!g_str_has_suffix(name, ".seq")) {
			str_append(&file, ".seq");
			if (!g_str_has_suffix(name, "_"))
				altfile = g_strdup_printf("%s_.seq", name);
		}
	}
	if (!is_readable_file(file) && (!altfile || !is_readable_file(altfile))) {
		if (check_seq()) {
			siril_log_error(_("No sequence `%s' found.\n"), name);
			g_free(file);
			g_free(altfile);
			return NULL;
		}
	}

	sequence *seq = NULL;
	if ((seq = readseqfile(file))) {
		if (get_filename) {
			*get_filename = file;
			file = NULL; // do not free
		}
	}
	else if (altfile && (seq = readseqfile(altfile))) {
		if (get_filename) {
			*get_filename = altfile;
			altfile = NULL; // do not free
		}
	}
	if (!seq)
		siril_log_error(_("Loading sequence `%s' failed.\n"), name);
	else {
		if (seq_check_basic_data(seq, FALSE) == -1) {
			free(seq);
			seq = NULL;
		}
		else if (seq->type == SEQ_SER)
			ser_display_info(seq->ser_file);
	}
	g_free(file);
	g_free(altfile);
	return seq;
}

/* callback functions */

/* COMPLETION_COLUMN, completion helpers, init_command, on_GtkCommandHelper moved to gui/script_console.c */

void log_several_lines(char *text) {
	char *line = text;
	do {
		char *eol = strchr(line, '\n');
		if (eol)
			*eol = '\0';
		siril_log_message("%s\n", line);
		if (!eol)
			break;
		line = eol+1;
	} while (line[0] != '\0');
}
