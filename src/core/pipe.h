#ifndef _PIPE_H_
#define _PIPE_H_

/* These named pipe functions are not reentrant */

typedef enum {
	PIPE_LOG, PIPE_STATUS, PIPE_PROGRESS, PIPE_READY
} pipe_message;

typedef enum {
	PIPE_STARTING, PIPE_SUCCESS, PIPE_ERROR, PIPE_EXIT, PIPE_NA, PIPE_BUSY
} pipe_verb;

typedef enum {
	PIPE_READ_EMPTY = -2, PIPE_READ_ERROR, PIPE_READ_UNUSED, PIPE_READ_SUCCESS
} pipe_read_status;

/* not reentrant, arg needs the newline for logs but not for others */
int pipe_send_message(pipe_message msgtype, pipe_verb verb, const char *arg);

int pipe_start(char *r_path, char *w_path);
void pipe_stop();
void *read_pipe(void *p);

#endif
