/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_CORE_SIRILBACKTRACE_H_
#define SRC_CORE_SIRILBACKTRACE_H_

#include "siril.h"

typedef struct _SirilBacktraceAddressInfo SirilBacktraceAddressInfo;

struct _SirilBacktraceAddressInfo {
	gchar object_name[256];

	gchar symbol_name[256];
	guintptr symbol_address;

	gchar source_file[256];
	gint source_line;
};

void siril_backtrace_init(void);

gboolean siril_backtrace_start(void);
void siril_backtrace_stop(void);

SirilBacktrace* siril_backtrace_new(gboolean include_current_thread);
void siril_backtrace_free(SirilBacktrace *backtrace);

gint siril_backtrace_get_n_threads(SirilBacktrace *backtrace);
guintptr siril_backtrace_get_thread_id(SirilBacktrace *backtrace, gint thread);
const gchar* siril_backtrace_get_thread_name(SirilBacktrace *backtrace,
		gint thread);
gboolean siril_backtrace_is_thread_running(SirilBacktrace *backtrace,
		gint thread);

gint siril_backtrace_find_thread_by_id(SirilBacktrace *backtrace,
		guintptr thread_id, gint thread_hint);

gint siril_backtrace_get_n_frames(SirilBacktrace *backtrace, gint thread);
guintptr siril_backtrace_get_frame_address(SirilBacktrace *backtrace,
		gint thread, gint frame);

#endif /* SRC_CORE_SIRILBACKTRACE_H_ */
