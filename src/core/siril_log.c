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
#include "core/siril_date.h"
#include "core/command.h" // for process_clear()
#include "core/OS_utils.h"
#include "core/pipe.h"
#include "gui/progress_and_log.h"

/* This function writes a message on Siril's console/log. It is not thread safe.
 * There is a limit in number of characters that it is able to write in one call: 1023.
 * Return value is the string printed from arguments, or NULL if argument was empty or
 * only newline. It is an allocated string and must not be freed. It can be
 * reused until next call to this function.
 */
static char* siril_log_internal(const char* format, const char* color, va_list arglist) {
	static char *msg = NULL;

	if (msg == NULL) {
		msg = malloc(1024);
		msg[1023] = '\0';
	}

	vsnprintf(msg, 1023, format, arglist);

	if (msg == NULL || msg[0] == '\0')
		return NULL;

	if (msg[0] == '\n' && msg[1] == '\0') {
		fputc('\n', stdout);
		gui_log_message("\n", NULL);
		return NULL;
	}

	g_print("log: %s", msg);
	pipe_send_message(PIPE_LOG, PIPE_NA, msg);
	gui_log_message(msg, color);

	return msg;
}

char* siril_log_message(const char* format, ...) {
	va_list args;
	va_start(args, format);
	g_mutex_lock(&com.mutex);
	char *msg = siril_log_internal(format, NULL, args);
	g_mutex_unlock(&com.mutex);
	va_end(args);
	return msg;
}

char* siril_log_color_message(const char* format, const char* color, ...) {
	va_list args;
	va_start(args, color);
	g_mutex_lock(&com.mutex);
	char *msg = siril_log_internal(format, color, args);
	g_mutex_unlock(&com.mutex);
	va_end(args);
	return msg;
}

const char *format_time_diff(struct timeval t_start, struct timeval t_end) {
	static char str[32];
	double start = (double) (t_start.tv_sec + t_start.tv_usec / 1.0E6);
	double end = (double) (t_end.tv_sec + t_end.tv_usec / 1.0E6);
	double diff = end - start;
	if (diff < 0.0) {
		str[0] = '\0';
	} else {
		if (diff >= 3600.0) {
			int hour = (int) diff / 3600;
			int sec = (int) diff % 3600;
			int min = sec / 60;
			sec = sec % 60;
			sprintf(str, _("%d h %02d min %.2d s"), hour, min, sec);
		} else if (diff >= 60.0) {
			int min = (int) diff / 60;
			int sec = (int) diff % 60;
			sprintf(str, _("%d min %02d s"), min, sec);
		} else if (diff < 1.0) {
			double ms = diff * 1.0E3;
			char ms_str[32];
			g_snprintf(ms_str, sizeof(ms_str), "%.2lf", ms);
			sprintf(str, _("%s ms"), ms_str);
		} else {
			char diff_str[32];
			g_snprintf(diff_str, sizeof(diff_str), "%.2lf", diff);
			sprintf(str, _("%s s"), diff_str);
		}
	}
	return str;
}

void show_time_msg(struct timeval t_start, struct timeval t_end, const char *msg) {
	siril_log_color_message("%s: %s\n", "green", msg, format_time_diff(t_start, t_end));
}

void show_time(struct timeval t_start, struct timeval t_end) {
	show_time_msg(t_start, t_end, _("Execution time"));
}

void get_min_sec_from_timevals(struct timeval t_start, struct timeval t_end,
		int *min, int *sec) {
	double start, end, diff;
	start = (double)(t_start.tv_sec + t_start.tv_usec / 1.0E6);
	end = (double)(t_end.tv_sec + t_end.tv_usec / 1.0E6);
	diff = end - start;
	*min = (int)diff / 60;
	*sec = (int)diff % 60;
}

