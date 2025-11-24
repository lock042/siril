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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "core/siril.h"
#include "algos/statistics.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/histogram.h"
#include "gui/script_menu.h"
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
			siril_log_message(_("This command cannot be used in a script: %s\n"), commands[i].name);
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
		if (get_thread_run()) {
			return CMD_THREAD_RUNNING;
		}
	}

	// process the command
	siril_log_color_message(_("Running command: %s\n"), "salmon", word[0]);
	fprintf(stdout, "%lu: running command %s\n", time(NULL), word[0]);
	int retval = commands[i].process(wordnb);
	if (gui.roi.active)
		populate_roi();
	if (retval & CMD_NOTIFY_GFIT_MODIFIED) {
		if (!com.python_script) {
			notify_gfit_modified();
		} else {
			waiting_for_thread();
			invalidate_stats_from_fit(gfit);
			invalidate_gfit_histogram();
			execute_idle_and_wait_for_it(end_gfit_operation, NULL);
		}
		retval = retval & ~CMD_NOTIFY_GFIT_MODIFIED;
	}
	return (int) retval;
}

static void update_log_icon(gboolean is_running) {
	GtkImage *image = GTK_IMAGE(lookup_widget("image_log"));
	if (is_running)
		gtk_image_set_from_icon_name(image, "gtk-yes", GTK_ICON_SIZE_LARGE_TOOLBAR);
	else
		gtk_image_set_from_icon_name(image, "gtk-no", GTK_ICON_SIZE_LARGE_TOOLBAR);
}

struct log_status_bar_idle_data {
	gchar *myline;
	int line;
};

static gboolean log_status_bar_idle_callback(gpointer p) {
	struct log_status_bar_idle_data *data = (struct log_status_bar_idle_data *) p;

	GtkStatusbar *statusbar_script = GTK_STATUSBAR(lookup_widget("statusbar_script"));
	gchar *status;
	gchar *newline;

	update_log_icon(TRUE);

	newline = g_strdup(data->myline);
	status = g_strdup_printf(_("Processing line %d: %s"), data->line, newline);

	gtk_statusbar_push(statusbar_script, 0, status);
	g_free(newline);
	g_free(status);
	g_free(data->myline);
	free(data);

	return FALSE;	// only run once
}

static void display_command_on_status_bar(int line, char *myline) {
	if (!com.headless) {
		struct log_status_bar_idle_data *data;

		data = malloc(sizeof(struct log_status_bar_idle_data));
		data->line = line;
		data->myline = myline ? g_strdup(myline) : NULL;
		gdk_threads_add_idle(log_status_bar_idle_callback, data);
	}
}

static void clear_status_bar() {
	GtkStatusbar *bar = GTK_STATUSBAR(lookup_widget("statusbar_script"));
	gtk_statusbar_remove_all(bar, 0);
	update_log_icon(FALSE);
}

static gboolean end_script(gpointer p) {
	/* GTK+ code is ignored during scripts, this is a good place to redraw everything */
	clear_status_bar();
	gui_function(set_GUI_CWD, NULL);
	gui_function(update_MenuItem, NULL);
	notify_gfit_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	update_zoom_label();
	update_display_fwhm();
	display_filename();
	gui_function(new_selection_zone, NULL);
	update_spinCPU(GINT_TO_POINTER(0));
	set_cursor_waiting(FALSE);
	return FALSE;
}

int check_requires(gboolean *checked_requires, gboolean is_required) {
	int retval = CMD_OK;
	/* check for requires command */
	if (!g_ascii_strcasecmp(word[0], "requires")) {
		*checked_requires = TRUE;
	} else if (is_required && !*checked_requires) {
		siril_log_color_message(_("The \"requires\" command is missing at the top of the script file."
					" This command is needed to check script compatibility.\n"), "red");
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
			siril_log_message(_("This command cannot be run while live stacking is active, ignoring.\n"));

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
			siril_log_color_message("%s\n", "blue", buffer);
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

		display_command_on_status_bar(line, buffer);
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
			siril_log_message(_("Error in line %d ('%s'): %s.\n"), line, buffer, cmd_err_to_str(retval));
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
		siril_debug_print("End of command %s, memory difference: %d MB\n", word[0], startmem - endmem);
		startmem = endmem;
		memset(word, 0, sizeof word);
		g_free (buffer);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);

	if (!com.headless) {
		com.script = FALSE;
		gui_function(end_script, NULL);
	}

	/* Now we want to restore the saved cwd */
	siril_change_dir(saved_cwd, NULL);
	writeinitfile();
	if (!retval) {
		siril_log_message(_("Script execution finished successfully.\n"));
		gettimeofday(&t_end, NULL);
		show_time_msg(t_start, t_end, _("Total execution time"));
	} else {
		char *msg = siril_log_message(_("Script execution failed.\n"));
		msg[strlen(msg) - 1] = '\0';
		set_progress_bar_data(msg, PROGRESS_DONE);
	}
	g_free(saved_cwd);

	if (com.script_thread) {
		siril_debug_print("Script thread exiting\n");
		com.script_thread_exited = TRUE;
	}
	/* If called from the GUI, re-enable widgets blocked during the script */
	siril_add_idle(script_widgets_idle, NULL);
	return GINT_TO_POINTER(retval);
}

static gboolean show_command_help_popup(gpointer user_data) {
	GtkEntry *entry = (GtkEntry*) user_data;
	gchar *helper = NULL;

	const gchar *text = gtk_entry_get_text(entry);
	if (*text == '\0') {
		helper = g_strdup(_("Please enter an existing command before hitting this button"));
	} else {
		command *current = commands;
		gchar **command_line = g_strsplit_set(text, " ", -1);
		while (current->process) {
			if (!g_ascii_strcasecmp(current->name, command_line[0])) {
				gchar **token;

				token = g_strsplit_set(current->usage, " \n", -1);
				GString *str = g_string_new(token[0]);
				str = g_string_prepend(str, "<span foreground=\"red\" size=\"larger\"><b>");
				str = g_string_append(str, "</b>");
				if (token[1] != NULL) {
					int i = 1;
					while (token[i]) {
						str = g_string_append(str, " ");
						if (!g_ascii_strcasecmp(current->name, token[i])) {
							str = g_string_append(str, "\n<b>");
						}
						str = g_string_append(str, token[i]);
						if (!g_ascii_strcasecmp(current->name, token[i])) {
							str = g_string_append(str, "</b>");
						}
						i++;
					}
				}
				str = g_string_append(str, "</span>\n\n\t");
				str = g_string_append(str, _(current->definition));
				str = g_string_append(str, "\n\n<b>");
				str = g_string_append(str, _("Can be used in a script: "));

				if (current->scriptable) {
					str = g_string_append(str, "<span foreground=\"green\">");
					str = g_string_append(str, _("YES"));
				} else {
					str = g_string_append(str, "<span foreground=\"red\">");
					str = g_string_append(str, _("NO"));
				}
				str = g_string_append(str, "</span></b>");
				helper = g_string_free(str, FALSE);
				g_strfreev(token);
				break;
			}
			current++;
		}
		g_strfreev(command_line);
	}
	if (!helper) {
		helper = g_strdup(_("No help for this command"));
	}

	GtkWidget *popover = popover_new(lookup_widget("command"), helper);
#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show(popover);
#endif
	g_free(helper);
	return FALSE;
}

/* handler for the single-line console */
#if GTK_CHECK_VERSION(3, 24, 24)
static gboolean on_command_key_press_event(GtkEventController *controller,
		guint keyval, guint keycode, GdkModifierType modifiers,
		GtkWidget *widget) {
#else
static gboolean on_command_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	guint keyval = event->keyval;
#endif
	int handled = 0;
	static GtkEntry *entry = NULL;
	if (!entry)
		entry = GTK_ENTRY(widget);
	GtkEditable *editable = GTK_EDITABLE(entry);
	int entrylength = 0;

	switch (keyval) {
	case GDK_KEY_Up:
		handled = 1;
		if (!gui.cmd_history)
			break;
		if (gui.cmd_hist_display > 0) {
			if (gui.cmd_history[gui.cmd_hist_display - 1])
				--gui.cmd_hist_display;
			// display previous entry
			gtk_entry_set_text(entry, gui.cmd_history[gui.cmd_hist_display]);
		} else if (gui.cmd_history[gui.cmd_hist_size - 1]) {
			// ring back, display previous
			gui.cmd_hist_display = gui.cmd_hist_size - 1;
			gtk_entry_set_text(entry, gui.cmd_history[gui.cmd_hist_display]);
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Down:
		handled = 1;
		if (!gui.cmd_history)
			break;
		if (gui.cmd_hist_display == gui.cmd_hist_current)
			break;
		if (gui.cmd_hist_display == gui.cmd_hist_size - 1) {
			if (gui.cmd_hist_current == 0) {
				// ring forward, end
				gtk_entry_set_text(entry, "");
				gui.cmd_hist_display++;
			} else if (gui.cmd_history[0]) {
				// ring forward, display next
				gui.cmd_hist_display = 0;
				gtk_entry_set_text(entry, gui.cmd_history[0]);
			}
		} else {
			if (gui.cmd_hist_display == gui.cmd_hist_current - 1) {
				// end
				gtk_entry_set_text(entry, "");
				gui.cmd_hist_display++;
			} else if (gui.cmd_history[gui.cmd_hist_display + 1]) {
				// display next
				gtk_entry_set_text(entry,
						gui.cmd_history[++gui.cmd_hist_display]);
			}
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Page_Up:
	case GDK_KEY_Page_Down:
		handled = 1;
		// go to first and last in history
		break;
	}
	return (handled == 1);
}

int processcommand(const char *line, gboolean wait_for_completion) {
	int wordnb = 0;
	GError *error = NULL;
	int ret = 0;

	if (line[0] == '\0' || line[0] == '\n')
		return CMD_NOT_FOUND;
	if (line[0] == '@') { // case of files
		if (get_thread_run() || (get_script_thread_run() && !com.script_thread_exited)) {
			PRINT_ANOTHER_THREAD_RUNNING;
			return CMD_THREAD_RUNNING;
		}
		if (get_script_thread_run())
			wait_for_script_thread();

		/* Switch to console tab */
		control_window_switch_to_tab(OUTPUT_LOGS);

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
		control_window_switch_to_tab(OUTPUT_LOGS);

		gchar *myline = strdup(line);
		int len = strlen(line);
		if (len > 0)
			g_print("input command:%s\n", myline);

		parse_line(myline, len, &wordnb);
		com.command = TRUE;
		ret = execute_command(wordnb);

		if (ret) {
			com.command = FALSE;
			siril_log_color_message(_("Command execution failed: %s.\n"), "red", cmd_err_to_str(ret));
			if (!(com.script || com.python_script) && !com.headless &&
				(ret == CMD_WRONG_N_ARG || ret == CMD_ARG_ERROR)) {
				gui_function(show_command_help_popup, GTK_ENTRY(lookup_widget("command")));
			}
			free(myline);
			return ret;
		}

		if (wait_for_completion && ret != CMD_NO_WAIT) {
			while (get_thread_run()) {
				if (waiting_for_thread()) {
					ret = CMD_GENERIC_ERROR;  // Command failed during execution
					break;
				}
				g_usleep(100000);  // Sleep for 100ms to avoid busy waiting
				remove_child_from_children((GPid) -2); // remove processing thread from list
			}
		}
		com.command = FALSE;
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
			siril_log_color_message(_("No sequence `%s' found.\n"), "red", name);
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
		siril_log_color_message(_("Loading sequence `%s' failed.\n"), "red", name);
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

#define COMPLETION_COLUMN 0

static gboolean on_match_selected(GtkEntryCompletion *widget, GtkTreeModel *model,
		GtkTreeIter *iter, gpointer user_data) {
	const gchar *cmd;
	GtkEditable *e = (GtkEditable *) gtk_entry_completion_get_entry(widget);
	gchar *s = gtk_editable_get_chars(e, 0, -1);
	gint cur_pos = gtk_editable_get_position(e);
	gint p = cur_pos;
	gchar *end;

	gtk_tree_model_get(model, iter, COMPLETION_COLUMN, &cmd, -1);

	end = s + cur_pos;
	gint del_end_pos = end - s + 1;

	gtk_editable_delete_text(e, 0, del_end_pos);
	gtk_editable_insert_text(e, cmd, -1, &p);
	gtk_editable_set_position(e, p);

	return TRUE;
}

static gboolean completion_match_func(GtkEntryCompletion *completion,
		const gchar *key, GtkTreeIter *iter, gpointer user_data) {
	gboolean res = FALSE;
	char *tag = NULL;

	if (*key == '\0') return FALSE;

	GtkTreeModel *model = gtk_entry_completion_get_model(completion);
	int column = gtk_entry_completion_get_text_column(completion);

	if (gtk_tree_model_get_column_type(model, column) != G_TYPE_STRING)
		return FALSE;

	gtk_tree_model_get(model, iter, column, &tag, -1);

	if (tag) {
		char *normalized = g_utf8_normalize(tag, -1, G_NORMALIZE_ALL);
		if (normalized) {
			char *casefold = g_utf8_casefold(normalized, -1);
			if (casefold) {
				res = g_strstr_len(casefold, -1, key) != NULL;
			}
			g_free(casefold);
		}
		g_free(normalized);
		g_free(tag);
	}

	return res;
}

static void init_completion_command() {
	GtkEntryCompletion *completion = gtk_entry_completion_new();
	GtkListStore *model = gtk_list_store_new(1, G_TYPE_STRING);
	GtkTreeIter iter;
	GtkEntry *entry = GTK_ENTRY(lookup_widget("command"));

	gtk_entry_completion_set_model(completion, GTK_TREE_MODEL(model));
	gtk_entry_completion_set_text_column(completion, COMPLETION_COLUMN);
	gtk_entry_completion_set_minimum_key_length(completion, 2);
	gtk_entry_completion_set_popup_completion(completion, TRUE);
	gtk_entry_completion_set_inline_completion(completion, TRUE);
	gtk_entry_completion_set_popup_single_match(completion, FALSE);
	gtk_entry_completion_set_match_func(completion, completion_match_func, NULL, NULL);
	gtk_entry_set_completion(entry, completion);
	g_signal_connect(G_OBJECT(completion), "match-selected", G_CALLBACK(on_match_selected), NULL);

	/* Populate the completion database. */
	command *current = commands;

	while (current->process) {
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter, COMPLETION_COLUMN, current->name, -1);
		current++;
	}
	g_object_unref(model);
}

static void init_controller_command() {
	GtkWidget *widget = lookup_widget("command");

#if GTK_CHECK_VERSION(3, 24, 24)
	GtkEventController *controller = gtk_event_controller_key_new(widget);
	g_signal_connect(controller, "key-pressed", G_CALLBACK(on_command_key_press_event), widget);
#else
	g_signal_connect(widget, "key-press-event", G_CALLBACK(on_command_key_press_event), NULL);
#endif
}

void init_command() {
	init_completion_command();
	init_controller_command();
}

void on_GtkCommandHelper_clicked(GtkButton *button, gpointer user_data) {
	show_command_help_popup((GtkEntry *)user_data);
}

/** Callbacks **/

/*
 * Command line history static function
 */

static void history_add_line(char *line) {
	if (!gui.cmd_history) {
		gui.cmd_hist_size = CMD_HISTORY_SIZE;
		gui.cmd_history = calloc(gui.cmd_hist_size, sizeof(const char*));
		gui.cmd_hist_current = 0;
		gui.cmd_hist_display = 0;
	}
	gui.cmd_history[gui.cmd_hist_current] = line;
	gui.cmd_hist_current++;
	// circle at the end
	if (gui.cmd_hist_current == gui.cmd_hist_size)
		gui.cmd_hist_current = 0;
	if (gui.cmd_history[gui.cmd_hist_current]) {
		free(gui.cmd_history[gui.cmd_hist_current]);
		gui.cmd_history[gui.cmd_hist_current] = NULL;
	}
	gui.cmd_hist_display = gui.cmd_hist_current;
}

void on_command_activate(GtkEntry *entry, gpointer user_data) {
	const gchar *text = gtk_entry_get_text(entry);
	history_add_line(strdup(text));
	if (!(processcommand(text, FALSE))) {
		gtk_entry_set_text(entry, "");
		gui_function(set_precision_switch, NULL);
	}
}

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
