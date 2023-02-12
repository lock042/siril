#ifndef SRC_CORE_COMMAND_LINE_PROCESSOR_H_
#define SRC_CORE_COMMAND_LINE_PROCESSOR_H_

#include <glib.h>

typedef enum {
	CMD_NOT_FOUND = -999,
	CMD_NO_WAIT = -667,
	CMD_NO_CWD = -500,
	CMD_NOT_SCRIPTABLE = -99,
	CMD_WRONG_N_ARG = -9,
	CMD_ARG_ERROR = -8,
	CMD_SELECTION_ERROR = -5,
	CMD_OK = 0,
	CMD_GENERIC_ERROR = 1,
	CMD_IMAGE_NOT_FOUND = 10,
	CMD_SEQUENCE_NOT_FOUND = 11,
	CMD_INVALID_IMAGE = 12,
	CMD_LOAD_IMAGE_FIRST = 20,
	CMD_ONLY_SINGLE_IMAGE = 30,
	CMD_NOT_FOR_SINGLE = 40,
	CMD_NOT_FOR_MONO = 50,
	CMD_NOT_FOR_RGB = 60,
	CMD_FOR_CFA_IMAGE = 70,
	CMD_FILE_NOT_FOUND = 80,
	CMD_FOR_PLATE_SOLVED = 90,
	CMD_NEED_INIT_FIRST = 100,
	CMD_ALLOC_ERROR = 666,
	CMD_THREAD_RUNNING = 999
} cmd_errors;

void parse_line(char *myline, int len, int *nb);
void remove_trailing_cr(char *str);
int execute_command(int wordnb);
int check_requires(gboolean *checked_requires, gboolean is_required);
void log_several_lines(char *text);

gpointer execute_script(gpointer p);
int processcommand(const char *line);
void init_command();
sequence *load_sequence(const char *name, char **get_filename);

#endif /* SRC_CORE_COMMAND_LINE_PROCESSOR_H_ */
