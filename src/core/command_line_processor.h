#ifndef SRC_CORE_COMMAND_LINE_PROCESSOR_H_
#define SRC_CORE_COMMAND_LINE_PROCESSOR_H_

#define PRINT_LOAD_IMAGE_FIRST siril_log_color_message(_("Load an image or a sequence first.\n"), "red")
#define PRINT_ONLY_SINGLE_IMAGE siril_log_color_message(_("Single image must be loaded, and this command cannot be applied on a sequence.\n"), "red")
#define PRINT_NOT_FOR_SINGLE siril_log_color_message(_("This command can only be used when a sequence is loaded.\n"), "red")
#define PRINT_NOT_FOR_MONO siril_log_color_message(_("This command cannot be applied on monochrome images.\n"), "red")
#define PRINT_NOT_FOR_RGB siril_log_color_message(_("This command cannot be applied on rgb images.\n"), "red")
#define PRINT_FOR_CFA_IMAGE PRINT_NOT_FOR_RGB; siril_log_color_message(_("Make sure your image is in CFA mode.\n"), "red")

enum {
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
	CMD_NEED_INIT_FIRST = 100,
	CMD_ALLOC_ERROR = 666,
	CMD_THREAD_RUNNING = 999
};

void parse_line(char *myline, int len, int *nb);
void remove_trailing_cr(char *str);
int execute_command(int wordnb);
int check_requires(gboolean *checked_requires);
void log_several_lines(char *text);

gpointer execute_script(gpointer p);
int	processcommand(const char *line);
void init_command();
sequence *load_sequence(const char *name, char **get_filename);

#endif /* SRC_CORE_COMMAND_LINE_PROCESSOR_H_ */
