#ifndef SRC_CORE_COMMAND_LINE_PROCESSOR_H_
#define SRC_CORE_COMMAND_LINE_PROCESSOR_H_

#define CMD_NO_WAIT -667

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
