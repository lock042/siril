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
#ifndef SRC_CORE_OS_UTILS_H_
#define SRC_CORE_OS_UTILS_H_

#include <glib.h>

#include "core/siril.h"

#ifdef __cplusplus
extern "C" {
#endif

gboolean is_space_disk_available(const gchar *disk);
gboolean update_displayed_memory(gpointer data);
int test_available_space(gint64 req_size);

guint64 get_available_memory();
int get_max_memory_in_MB();
void log_used_mem(gchar *when);

int get_available_cpu_cgroups();

void init_num_procs();

long get_pathmax(void);

GInputStream *siril_input_stream_from_stdin();

#ifdef _WIN32
gchar *get_special_folder(int csidl);
int ReconnectIO(int OpenNewConsole);
char* siril_real_path(const char *source);
gchar *get_siril_bundle_path();
gchar *find_executable_in_path(const char *exe_name, const char *path);
#endif

gboolean allow_to_open_files(int nb_frames, int *nb_allowed_file);

#ifdef __cplusplus
}
#endif

#endif /* SRC_CORE_OS_UTILS_H_ */
