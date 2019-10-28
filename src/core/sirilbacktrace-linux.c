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

#define _GNU_SOURCE

#include "config.h"

#include <gio/gio.h>

#include "sirilbacktrace-backend.h"

#ifdef SIRIL_BACKTRACE_BACKEND_LINUX

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <signal.h>
#include <execinfo.h>
#include <dlfcn.h>
#include <string.h>
#include <stdio.h>

#include "siril.h"

#ifdef HAVE_LIBBACKTRACE
#include <backtrace.h>
#endif

#ifdef HAVE_LIBUNWIND
#define UNW_LOCAL_ONLY
#include <libunwind.h>
#endif

#include "sirilbacktrace.h"

#define MAX_N_THREADS        256
#define MAX_N_FRAMES         256
#define MAX_THREAD_NAME_SIZE 32
#define N_SKIPPED_FRAMES     2
#define MAX_WAIT_TIME        (G_TIME_SPAN_SECOND / 20)
#define BACKTRACE_SIGNAL     SIGUSR1

typedef struct _SirilBacktraceThread SirilBacktraceThread;

struct _SirilBacktraceThread {
	pid_t tid;
	gchar name[MAX_THREAD_NAME_SIZE];
	gchar state;

	guintptr frames[MAX_N_FRAMES];
	gint n_frames;
};

struct _SirilBacktrace {
	SirilBacktraceThread *threads;
	gint n_threads;
};

/*  local function prototypes  */

static inline gint siril_backtrace_normalize_frame(SirilBacktrace *backtrace,
		gint thread, gint frame);

static gint siril_backtrace_enumerate_threads(gboolean include_current_thread,
		pid_t *threads, gint size);
static void siril_backtrace_read_thread_name(pid_t tid, gchar *name, gint size);
static gchar siril_backtrace_read_thread_state(pid_t tid);

static void siril_backtrace_signal_handler(gint signum);

/*  static variables  */

static GMutex mutex;
static gint n_initializations;
static gboolean initialized;
static struct sigaction orig_action;
static pid_t blacklisted_threads[MAX_N_THREADS];
static gint n_blacklisted_threads;
static SirilBacktrace *handler_backtrace;
static gint handler_n_remaining_threads;
static gint handler_lock;

#ifdef HAVE_LIBBACKTRACE
static struct backtrace_state *backtrace_state;
#endif

static const gchar *const blacklisted_thread_names[] =
		{ "gmain", "threaded-ml" };

/*  private functions  */

static inline gint siril_backtrace_normalize_frame(SirilBacktrace *backtrace,
		gint thread, gint frame) {
	if (frame >= 0)
		return frame + N_SKIPPED_FRAMES;
	else
		return backtrace->threads[thread].n_frames + frame;
}

static gint siril_backtrace_enumerate_threads(gboolean include_current_thread,
		pid_t *threads, gint size) {
	DIR *dir;
	struct dirent *dirent;
	pid_t tid;
	gint n_threads;

	dir = opendir("/proc/self/task");

	if (!dir)
		return 0;

	tid = syscall(SYS_gettid);

	n_threads = 0;

	while (n_threads < size && (dirent = readdir(dir))) {
		pid_t id = g_ascii_strtoull(dirent->d_name, NULL, 10);

		if (id) {
			if (!include_current_thread && id == tid)
				id = 0;
		}

		if (id) {
			gint i;

			for (i = 0; i < n_blacklisted_threads; i++) {
				if (id == blacklisted_threads[i]) {
					id = 0;

					break;
				}
			}
		}

		if (id)
			threads[n_threads++] = id;
	}

	closedir(dir);

	return n_threads;
}

static void siril_backtrace_read_thread_name(pid_t tid, gchar *name, gint size) {
	gchar filename[64];
	gint fd;

	if (size <= 0)
		return;

	name[0] = '\0';

	g_snprintf(filename, sizeof(filename), "/proc/self/task/%llu/comm",
			(unsigned long long) tid);

	fd = open(filename, O_RDONLY);

	if (fd >= 0) {
		gint n = read(fd, name, size);

		if (n > 0)
			name[n - 1] = '\0';

		close(fd);
	}
}

static gchar siril_backtrace_read_thread_state(pid_t tid) {
	gchar buffer[64];
	gint fd;
	gchar state = '\0';

	g_snprintf(buffer, sizeof(buffer), "/proc/self/task/%llu/stat",
			(unsigned long long) tid);

	fd = open(buffer, O_RDONLY);

	if (fd >= 0) {
		gint n = read(fd, buffer, sizeof(buffer));

		if (n > 0)
			buffer[n - 1] = '\0';

		sscanf(buffer, "%*d %*s %c", &state);

		close(fd);
	}

	return state;
}

static void siril_backtrace_signal_handler(gint signum) {
	SirilBacktrace *curr_backtrace;
	gint lock;

	do {
		lock = g_atomic_int_get(&handler_lock);

		if (lock < 0)
			continue;
	} while (!g_atomic_int_compare_and_exchange(&handler_lock, lock, lock + 1));

	curr_backtrace = g_atomic_pointer_get(&handler_backtrace);

	if (curr_backtrace) {
		pid_t tid = syscall(SYS_gettid);
		gint i;

		for (i = 0; i < curr_backtrace->n_threads; i++) {
			SirilBacktraceThread *thread = &curr_backtrace->threads[i];

			if (thread->tid == tid) {
				thread->n_frames = backtrace((gpointer*) thread->frames,
				MAX_N_FRAMES);

				g_atomic_int_dec_and_test(&handler_n_remaining_threads);

				break;
			}
		}
	}

	g_atomic_int_dec_and_test(&handler_lock);
}

/*  public functions  */

void siril_backtrace_init(void) {
#ifdef HAVE_LIBBACKTRACE
	backtrace_state = backtrace_create_state(NULL, 0, NULL, NULL);
#endif
}

gboolean siril_backtrace_start(void) {
	g_mutex_lock(&mutex);

	if (n_initializations == 0) {
		struct sigaction action = { };

		action.sa_handler = siril_backtrace_signal_handler;
		action.sa_flags = SA_RESTART;

		sigemptyset(&action.sa_mask);

		if (sigaction(BACKTRACE_SIGNAL, &action, &orig_action) == 0) {
			pid_t *threads;
			gint n_threads;
			gint i;

			n_blacklisted_threads = 0;

			threads = g_new(pid_t, MAX_N_THREADS);

			n_threads = siril_backtrace_enumerate_threads(TRUE, threads,
					MAX_N_THREADS);

			for (i = 0; i < n_threads; i++) {
				gchar name[MAX_THREAD_NAME_SIZE];
				gint j;

				siril_backtrace_read_thread_name(threads[i], name,
						MAX_THREAD_NAME_SIZE);

				for (j = 0; j < G_N_ELEMENTS(blacklisted_thread_names); j++) {
					if (!strcmp(name, blacklisted_thread_names[j])) {
						blacklisted_threads[n_blacklisted_threads++] =
								threads[i];
					}
				}
			}

			g_free(threads);

			initialized = TRUE;
		}
	}

	n_initializations++;

	g_mutex_unlock(&mutex);

	return initialized;
}

void siril_backtrace_stop(void) {
	g_return_if_fail(n_initializations > 0);

	g_mutex_lock(&mutex);

	n_initializations--;

	if (n_initializations == 0 && initialized) {
		if (sigaction(BACKTRACE_SIGNAL, &orig_action, NULL) < 0)
			g_warning("failed to restore origianl backtrace signal handler");

		initialized = FALSE;
	}

	g_mutex_unlock(&mutex);
}

SirilBacktrace *siril_backtrace_new(gboolean include_current_thread) {
	SirilBacktrace * backtrace;
	pid_t pid;
	pid_t *threads;
	gint n_threads;
	gint64 start_time;
	gint i;

	if (!initialized)
		return NULL;

	pid = getpid();

	threads = g_new(pid_t, MAX_N_THREADS);

	n_threads = siril_backtrace_enumerate_threads(include_current_thread,
			threads, MAX_N_THREADS);

	if (n_threads == 0) {
		g_free(threads);

		return NULL;
	}

	g_mutex_lock(&mutex);

	backtrace = g_slice_new(SirilBacktrace);

	backtrace->threads = g_new(SirilBacktraceThread, n_threads);
	backtrace->n_threads = n_threads;

	while (!g_atomic_int_compare_and_exchange(&handler_lock, 0, -1))
		;

	g_atomic_pointer_set(&handler_backtrace, backtrace);
	g_atomic_int_set(&handler_n_remaining_threads, n_threads);

	g_atomic_int_set(&handler_lock, 0);

	for (i = 0; i < n_threads; i++) {
		SirilBacktraceThread *thread = &backtrace->threads[i];

		thread->tid = threads[i];
		thread->n_frames = 0;

		siril_backtrace_read_thread_name(thread->tid, thread->name,
				MAX_THREAD_NAME_SIZE);

		thread->state = siril_backtrace_read_thread_state(thread->tid);

		syscall(SYS_tgkill, pid, threads[i], BACKTRACE_SIGNAL);
	}

	g_free(threads);

	start_time = g_get_monotonic_time();

	while (g_atomic_int_get (&handler_n_remaining_threads) > 0) {
		gint64 time = g_get_monotonic_time();

		if (time - start_time > MAX_WAIT_TIME)
			break;

		g_usleep(1000);
	}

	while (!g_atomic_int_compare_and_exchange(&handler_lock, 0, -1))
		;

	g_atomic_pointer_set(&handler_backtrace, NULL);

	g_atomic_int_set(&handler_lock, 0);

	g_mutex_unlock(&mutex);

	if (n_threads == 0) {
		siril_backtrace_free(backtrace);

		return NULL;
	}

	return backtrace;
}

void siril_backtrace_free(SirilBacktrace *backtrace) {
	if (!backtrace)
		return;

	g_free(backtrace->threads);

	g_slice_free(SirilBacktrace, backtrace);
}

gint siril_backtrace_get_n_threads(SirilBacktrace *backtrace) {
	g_return_val_if_fail(backtrace != NULL, 0);

	return backtrace->n_threads;
}

guintptr siril_backtrace_get_thread_id(SirilBacktrace *backtrace, gint thread) {
	g_return_val_if_fail(backtrace != NULL, 0);
	g_return_val_if_fail(thread >= 0 && thread < backtrace->n_threads, 0);

	return backtrace->threads[thread].tid;
}

const gchar*
siril_backtrace_get_thread_name(SirilBacktrace *backtrace, gint thread) {
	g_return_val_if_fail(backtrace != NULL, NULL);
	g_return_val_if_fail(thread >= 0 && thread < backtrace->n_threads, NULL);

	if (backtrace->threads[thread].name[0])
		return backtrace->threads[thread].name;
	else
		return NULL;
}

gboolean siril_backtrace_is_thread_running(SirilBacktrace *backtrace, gint thread) {
	g_return_val_if_fail(backtrace != NULL, FALSE);
	g_return_val_if_fail(thread >= 0 && thread < backtrace->n_threads, FALSE);

	return backtrace->threads[thread].state == 'R';
}

gint siril_backtrace_find_thread_by_id(SirilBacktrace *backtrace,
		guintptr thread_id, gint thread_hint) {
	pid_t tid = thread_id;
	gint i;

	g_return_val_if_fail(backtrace != NULL, -1);

	if (thread_hint >= 0 && thread_hint < backtrace->n_threads
			&& backtrace->threads[thread_hint].tid == tid) {
		return thread_hint;
	}

	for (i = 0; i < backtrace->n_threads; i++) {
		if (backtrace->threads[i].tid == tid)
			return i;
	}

	return -1;
}

gint siril_backtrace_get_n_frames(SirilBacktrace *backtrace, gint thread) {
	g_return_val_if_fail(backtrace != NULL, 0);
	g_return_val_if_fail(thread >= 0 && thread < backtrace->n_threads, 0);

	return MAX(backtrace->threads[thread].n_frames - N_SKIPPED_FRAMES, 0);
}

guintptr siril_backtrace_get_frame_address(SirilBacktrace *backtrace, gint thread,
		gint frame) {
	g_return_val_if_fail(backtrace != NULL, 0);
	g_return_val_if_fail(thread >= 0 && thread < backtrace->n_threads, 0);

	frame = siril_backtrace_normalize_frame(backtrace, thread, frame);

	g_return_val_if_fail(
			frame >= N_SKIPPED_FRAMES && frame < backtrace->threads[thread].n_frames,
			0);

	return backtrace->threads[thread].frames[frame];
}

#endif /* GIMP_BACKTRACE_BACKEND_LINUX */
