/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include "core/siril_log.h"
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
static char* siril_log_internal(logpriority priority, const char* format, const char* color, va_list arglist) {
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
	if (priority >= com.log_threshold) {
		g_print(_("loglevel %d: %s"), (int) priority, msg);
		pipe_send_message(PIPE_LOG, PIPE_NA, msg);
	}
	if (priority >= gui.log_threshold)
		gui_log_message(msg, color);

	return msg;
}

char* siril_log_va(logpriority priority, const char* format, va_list args) {
	g_mutex_lock(&com.mutex);
	char *color;
	switch (priority) {
	// Case fallthrough is intentional
		case LOG_CRITICAL:
		case LOG_ERROR:
			color = "red";
			break;
		case LOG_WARNING:
			color = "salmon";
			break;
		case LOG_SUMMARY_INFO:
			color = "green";
			break;
		case LOG_DEBUG:
		case LOG_DEBUG_VERBOSE:
		case LOG_DEBUG_EXTRA_VERBOSE:
			color = "plum";
			break;
		case LOG_METATRON:
			color = "bold";
			break;
		default:
			color = NULL;
	}
	char *msg = siril_log_internal(priority, format, color, args);
//	free(color);
	g_mutex_unlock(&com.mutex);
	return msg;
}

char* siril_log(logpriority priority, const char* format, ...) {
	va_list args;
	va_start(args, format);
	char *msg = siril_log_va(priority, format, args);
	va_end(args);
	return msg;
}

/* Deprecated. Acts as siril_log(LOG_INFO...) and ideally all instances should
 * be migrated for consistency, but the function is kept available to make it
 * easier to merge MRs that are in-progress. In future new MRs should use
 * siril_log() with the correct log priority.
 */
char* siril_log_message(const char* format, ...) {
	va_list args;
	va_start(args, format);
	char *msg = siril_log_va(LOG_INFO, format, args);
	va_end(args);
	return msg;
}

/* Deprecated and ugly. Guesses the priority based on the color, but up to now
 * use of siril_log_color_message hasn't been consistent, so all instances ought
 * to be checked and converted to siril_log(). The function is kept available to
 * make it easier to merge MRs that are in-progress but all new MRs should use
 * siril_log()
 */
char* siril_log_color_message(const char* format, const char* color, ...) {
	va_list args;
	va_start(args, color);
//	g_mutex_lock(&com.mutex);
	logpriority priority;
	if (!strcmp(color, "red"))
		priority = LOG_ERROR;
	else if (!strcmp(color, "salmon"))
		priority = LOG_WARNING;
	else if (!strcmp(color, "green"))
		priority = LOG_SUMMARY_INFO;
	else
		priority = LOG_INFO;
//	char *msg = siril_log_internal(priority, format, color, args);
	char *msg = siril_log_va(priority, format, args);
//	g_mutex_unlock(&com.mutex);
	va_end(args);
	return msg;
}

void siril_debug_print(const char* format, ...) {
	// Don't waste time if we don't need to...
	if (com.headless && com.log_threshold > LOG_DEBUG)
		return;
	if (min(com.log_threshold, gui.log_threshold) > LOG_DEBUG)
		return;

	va_list args;
	va_start(args, format);
	(void) siril_log_va(LOG_DEBUG, format, args);
	va_end(args);
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
			sprintf(str, _("%.2lf ms"), ms);
		} else {
			sprintf(str, _("%.2lf s"), diff);
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

