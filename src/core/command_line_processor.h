#ifndef SRC_CORE_COMMAND_LINE_PROCESSOR_H_
#define SRC_CORE_COMMAND_LINE_PROCESSOR_H_

#include <glib.h>

typedef enum {
	CMD_NOT_FOUND = 1,
	CMD_NO_WAIT = 1 << 1,
	CMD_NO_CWD = 1 << 2,
	CMD_NOT_SCRIPTABLE = 1 << 3,
	CMD_WRONG_N_ARG = 1 << 4,
	CMD_ARG_ERROR = 1 << 5,
	CMD_SELECTION_ERROR = 1 << 6,
	CMD_OK = 0,
	CMD_GENERIC_ERROR = 1 << 7,
	CMD_IMAGE_NOT_FOUND = 1 << 8,
	CMD_SEQUENCE_NOT_FOUND = 1 << 9,
	CMD_INVALID_IMAGE = 1 << 10,
	CMD_LOAD_IMAGE_FIRST = 1 << 11,
	CMD_ONLY_SINGLE_IMAGE = 1 << 12,
	CMD_NOT_FOR_SINGLE = 1 << 13,
	CMD_NOT_FOR_MONO = 1 << 14,
	CMD_NOT_FOR_RGB = 1 << 15,
	CMD_FOR_CFA_IMAGE = 1 << 16,
	CMD_FILE_NOT_FOUND = 1 << 17,
	CMD_FOR_PLATE_SOLVED = 1 << 18,
	CMD_NEED_INIT_FIRST = 1 << 19,
	CMD_ALLOC_ERROR = 1 << 20,
	CMD_THREAD_RUNNING = 1 << 21,
	CMD_DIR_NOT_FOUND = 1 << 22
} cmd_errors;

typedef enum {
	CMD_NOTIFY_GFIT_MODIFIED = 1 << 23
} cmd_status;

void parse_line(char *myline, int len, int *nb);
void remove_trailing_cr(char *str);
int execute_command(int wordnb);
int check_requires(gboolean *checked_requires, gboolean is_required);
void log_several_lines(char *text);

gpointer execute_script(gpointer p);
int processcommand(const char *line, gboolean wait_for_completion);
void init_command();
sequence *load_sequence(const char *name, char **get_filename);

#endif /* SRC_CORE_COMMAND_LINE_PROCESSOR_H_ */
