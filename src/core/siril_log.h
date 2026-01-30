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

#ifndef _SIRIL_LOG_H
#define _SIRIL_LOG_H

#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif

char* siril_log_message(const char* format, ...);
char* siril_log_literal_message(const char* message);
char* siril_log_color_message(const char* format, const char* color, ...);
char* siril_log_literal_color_message(const char* message, const char* color);

void show_time(struct timeval, struct timeval);
void show_time_msg(struct timeval t_start, struct timeval t_end, const char *msg);
const char *format_time_diff(struct timeval t_start, struct timeval t_end);
void get_min_sec_from_timevals(struct timeval t_start, struct timeval t_end,
		int *min, int *sec);

#ifdef __cplusplus
}
#endif

#endif /* _SIRIL_LOG_H */
