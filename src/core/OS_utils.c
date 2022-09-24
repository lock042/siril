/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#include <direct.h>
#include <shlobj.h>
#include <tchar.h>
#include <io.h>
#include <fcntl.h>
#include <gio/gwin32inputstream.h>
#else
#include <sys/resource.h>
#include <gio/gunixinputstream.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif
#if defined(__unix__) || defined(OS_OSX)
#include <sys/param.h>		// define or not BSD macro
#endif
#ifdef OS_OSX
#include <mach/task.h>
#include <mach/mach_init.h>
#include <mach/mach_types.h>
#include <mach/mach_host.h>
#include <sys/sysctl.h>
#include <mach/vm_statistics.h>
#endif
#ifdef HAVE_SYS_STATVFS_H
#include <sys/statvfs.h>
#endif
#if HAVE_SYS_VFS_H
#include <sys/vfs.h>
#elif HAVE_SYS_MOUNT_H
#if HAVE_SYS_PARAM_H
#include <sys/param.h>
#endif
#include <sys/mount.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"

#include "OS_utils.h"

/**
 * Find the space remaining in a directory, in bytes.
 * @param name the path of the directory to be tested
 * @return the disk space remaining in bytes, or a negative value if error
 */
#if HAVE_SYS_STATVFS_H
static gint64 find_space(const gchar *name) {
	struct statvfs st;
	gint64 available;
	if (statvfs (name, &st))
		return (gint64) -1;
	available = st.f_bavail;        // force 64 bits
	return available * st.f_frsize;
}
#elif (HAVE_SYS_VFS_H || HAVE_SYS_MOUNT_H)
static gint64 find_space(const gchar *name) {
	struct statfs st;
	gint64 available;
	if (statfs (name, &st))
		return (gint64) -1;
	available = st.f_bavail;        // force 64 bits
        return available * st.f_bsize;
}
#elif defined _WIN32
static gint64 find_space(const gchar *name) {
	ULARGE_INTEGER avail;
	gint64 sz;

	gchar *localdir = g_path_get_dirname(name);
	wchar_t *wdirname = g_utf8_to_utf16(localdir, -1, NULL, NULL, NULL);

	if (!GetDiskFreeSpaceExW(wdirname, &avail, NULL, NULL))
		sz = (gint64) -1;
	else
		sz = avail.QuadPart;

	g_free(localdir);
	g_free(wdirname);
	return sz;
}
#else
static gint64 find_space(const gchar *name) {
	return (gint64) -1;
}
#endif /*HAVE_SYS_STATVFS_H*/

/**
 * Compute the used RAM and returns the value in bytes.
 * @return
 */
static guint64 get_used_RAM_memory() {
#if defined(__linux__) || defined(__CYGWIN__)
	static gboolean initialized = FALSE;
	static long page_size;
	static gint fd = -1;
	gchar buffer[128];
	gint size;
	guint64 resident;
	guint64 shared;

	if (!initialized) {
		page_size = getpagesize();

		if (page_size > 0)
			fd = g_open("/proc/self/statm", O_RDONLY);

		initialized = TRUE;
	}

	if (fd < 0)
		return (guint64) 0;

	if (lseek(fd, 0, SEEK_SET))
		return (guint64) 0;

	size = read(fd, buffer, sizeof(buffer) - 1);

	if (size <= 0)
		return (guint64) 0;

	buffer[size] = '\0';

	if (sscanf(buffer, "%*u %" G_GUINT64_FORMAT " %" G_GUINT64_FORMAT, &resident, &shared) != 2)
		return (guint64) 0;

	return (guint64) (resident /*- shared*/) * page_size;
#elif defined(OS_OSX)
#ifndef TASK_VM_INFO_REV0_COUNT /* phys_footprint added in REV1 */
	struct mach_task_basic_info t_info;
	mach_msg_type_number_t t_info_count = MACH_TASK_BASIC_INFO_COUNT;

	if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count) != KERN_SUCCESS)
		return (guint64) 0;
	else
		return ((guint64) t_info.resident_size);
#else
	task_vm_info_data_t t_info;
	mach_msg_type_number_t t_info_count = TASK_VM_INFO_COUNT;

	if (task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&t_info, &t_info_count) != KERN_SUCCESS)
		return (guint64) 0;
	else
		return ((guint64) t_info.phys_footprint);
#endif
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). In fact, it could work with linux */
	struct rusage usage;

	getrusage(RUSAGE_SELF, &usage);
	return ((guint64) usage.ru_maxrss * 1024UL);
#elif defined(_WIN32) /* Windows */
	PROCESS_MEMORY_COUNTERS memCounter;

	if (GetProcessMemoryInfo(GetCurrentProcess(), &memCounter, sizeof(memCounter)))
		return (memCounter.WorkingSetSize);
	return (guint64) 0;
#else
	return (guint64) 0;
#endif
}

/**
 * Updates RAM memory used by siril, available free disk space
 * and displays information on the control window.
 * @return always return TRUE
 */
gboolean update_displayed_memory() {
	set_GUI_MEM(get_used_RAM_memory(), "labelmem");
	set_GUI_DiskSpace(find_space(com.wd), "labelFreeSpace");
	set_GUI_DiskSpace(find_space(com.pref.swap_dir), "free_mem_swap");
	return TRUE;
}

#define MAX_COMP_FREESPACE_RATIO 3

/**
 * Test if there is enough free disk space by returning the difference
 * in bytes between available free disk space and the size given as parameter
 * @param req_size available space to be tested
 * @return 0 if there is enough disk space, 1 otherwise, -1 on error.
 */
int test_available_space(gint64 req_size) {
	gint64 free_space = find_space(com.wd);
	int res = -1;
	if (free_space < 0) {
		siril_log_message(_("Error while computing available free disk space.\n"));
		return res;
	}
	if (req_size <= 0) {
		siril_log_message(_("Error in requested space disk.\n"));
		return res;
	}

	if (req_size > free_space) {
		char * msg;
		gchar *avail = g_format_size_full(free_space, G_FORMAT_SIZE_IEC_UNITS);
		gchar *required = g_format_size_full(req_size, G_FORMAT_SIZE_IEC_UNITS);
		gchar *missing = g_format_size_full(req_size - free_space, G_FORMAT_SIZE_IEC_UNITS);
		if (com.pref.comp.fits_enabled) {
			if (req_size / free_space < MAX_COMP_FREESPACE_RATIO) {
				msg = siril_log_message(_("Compression enabled: There may no be enough free disk space to perform this operation: "
						"%s available for %s needed (missing %s)\n"),
						avail, required, missing);
				queue_message_dialog(GTK_MESSAGE_WARNING, _("Compression enabled: There may not be enough free disk space to perform this operation"), msg);
			} else {
				msg = siril_log_message(_("Compression enabled: It is likely that there is not enough free disk space to perform this operation: "
						"%s available for %s needed (missing %s)\n"),
						avail, required, missing);
				queue_message_dialog(GTK_MESSAGE_WARNING, _("Compression enabled: It is likely that there is not enough free disk space to perform this operation"), msg);
			}
			res = 0;
		} else {
			msg = siril_log_message(_("Not enough free disk space to perform this operation: "
						"%s available for %s needed (missing %s)\n"),
						avail, required, missing);
			queue_message_dialog(GTK_MESSAGE_ERROR, _("Not enough disk space"), msg);
			res = 1;
		}
		g_free(avail);
		g_free(required);
		g_free(missing);
		return res;
	}
	siril_debug_print("Tested free space ok: %" G_GINT64_FORMAT " for %" G_GINT64_FORMAT " MB free\n",
			(gint64)(req_size / BYTES_IN_A_MB), (gint64)(free_space / BYTES_IN_A_MB));
	return 0;
}

#if defined(__linux__) || defined(__CYGWIN__)
/* read a value from a file that contains only an integer, like many procfs or sysfs files */
static int read_from_file(const char *filename, guint64 *value) {
	int fd = open(filename, O_RDONLY);
	if (fd == -1)
		return 1;
	char buf[30];
	int retval = read(fd, buf, 30);
	if (retval <= 0) {
		close(fd);
		return 1;
	}
	char *end;
	*value = strtoull(buf, &end, 10);
	if (end == buf)
		retval = 1;
	else retval = 0;
	close(fd);
	return retval;
}

static int read_2_from_file(const char *filename, guint64 *value1, guint64 *value2) {
	int fd = open(filename, O_RDONLY);
	if (fd == -1)
		return 1;
	char buf[50];
	int retval = read(fd, buf, 50);
	if (retval <= 0) {
		close(fd);
		return 1;
	}
	char *end;
	*value1 = strtoull(buf, &end, 10);
	if (end == buf)
		return 1;

	char *start = end;
	*value2 = strtoull(start, &end, 10);
	retval = (end == start);
	close(fd);
	return retval;
}

/* returns the directory in the cgroups filesystem that applies for the current
 * process for the requested module */
static gchar *find_cgroups_path(const char *module) {
	FILE *fd = g_fopen("/proc/self/cgroup", "r");
	if (!fd)
		return g_strdup("");
	char buf[2000];
	char *cgpath = NULL;
	while (!cgpath && fgets(buf, 2000, fd)) {
		remove_trailing_eol(buf);
		gchar **tokens = g_strsplit(buf, ":", -1);
		guint n = g_strv_length(tokens);
		if (n < 3) {
			siril_debug_print("malformed line in /proc/self/cgroup: %s\n", buf);
			continue;
		}
		if (atoi(tokens[0]) == 0 && tokens[1][0] == '\0') {
			// cgroups v2, only one entry
			siril_debug_print("cgroups v2 path: %s\n", tokens[2]);
			cgpath = g_strdup(tokens[2]);
		} else {
			gchar **controllers = g_strsplit(tokens[1], ",", -1);
			guint ncont = g_strv_length(controllers);
			for (int i = 0; i < ncont; i++) {
				if (!strcmp(controllers[i], module)) {
					siril_debug_print("cgroups v1 path: %s\n", tokens[2]);
					if (!strcmp("/", tokens[2]))
						cgpath = g_strdup("");
					else cgpath = g_strdup(tokens[2]);
					break;
				}
			}
			g_strfreev(controllers);
		}
		g_strfreev(tokens);
	}

	fclose(fd);
	if (!cgpath)
		cgpath = g_strdup("");
	return cgpath;
}

static int get_available_mem_cgroups(guint64 *amount) {
	/* useful files, cgroups v1:
	 * memory.usage_in_bytes	# show current usage for memory	(for the system?)
	 * memory.memsw.usage_in_bytes	# show current usage for memory+Swap (system?)
	 * memory.limit_in_bytes	# set/show limit of memory usage
	 * memory.soft_limit_in_bytes	# set/show soft limit of memory usage
	 *
	 * useful files, cgroups v2:
	 * memory.low			# soft limit, default is '0'
	 * memory.high			# hard limit, default is 'max'
	 * memory.max			# kill limit, default is 'max'
	 */
	static gchar *limits_filepath = NULL;
	static gboolean initialized = FALSE;
	// assuming cgroups fs is mounted on /sys/fs/cgroup, we could also
	// check the mount points but that's fairly common
	const char *limits_paths[] = {	// in order of files to try
		// v1
		"/sys/fs/cgroup/memory%s/memory.soft_limit_in_bytes",
		"/sys/fs/cgroup/memory%s/memory.limit_in_bytes",
		"/sys/fs/cgroup/memory/memory.soft_limit_in_bytes",
		"/sys/fs/cgroup/memory/memory.limit_in_bytes",
		// v2
		"/sys/fs/cgroup%s/memory.low",
		"/sys/fs/cgroup%s/memory.high",
		"/sys/fs/cgroup%s/memory.max"
	};

	if (initialized && !limits_filepath)
		return 1;
	// first, get the limit, currently it can be dynamic but it's slower to read the file every time
	guint64 limit = 0;
	if (!initialized) {
		gchar *cgroup_path = find_cgroups_path("memory");
		int nb_paths = sizeof limits_paths / sizeof(const char *);
		int source_file;
		for (source_file = 0; source_file < nb_paths; source_file++) {
			gchar *path = g_strdup_printf(limits_paths[source_file], cgroup_path);
			if (!read_from_file(path, &limit) && limit > (guint64)0 && limit < 0x7fffffffffff0000) {
				siril_debug_print("Found memory cgroups limit in %s\n", path);
				siril_log_message(_("Using cgroups limit on memory: %u MB\n"),
						(unsigned int)(limit / BYTES_IN_A_MB));
				limits_filepath = path;
				break;
			}
			g_free(path);
		}
		initialized = 1;
		g_free(cgroup_path);
		if (source_file == nb_paths) {
			siril_debug_print("no memory cgroup controller detected\n");
			// not using cgroups
			return 1;
		}
	}
	else {
		if (read_from_file(limits_filepath, &limit) || limit == (guint64)0) {
			siril_log_message(_("Error reading from %s, disabling cgroups memory limits\n"),
					limits_filepath);
			g_free(limits_filepath);
			limits_filepath = NULL;
		}
	}

	// then, get the current amount
	guint64 current = get_used_RAM_memory();

	siril_debug_print("current memory: %d, cgroup limit: %d MB\n",
			(int)(current / BYTES_IN_A_MB), (int)(limit / BYTES_IN_A_MB));
	if (limit < current)
		*amount = (guint64)1;	// 0 mean error in the caller
	else *amount = limit - current;
	return 0;
}

int get_available_cpu_cgroups() {
	/* computing a number of processors based on the cgroups throttles (cpu bandwidth)
	 * Using files from the cgroups filesystem mounted in /sys/fs/cgroup
	 * cgroups v1 (see https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/6/html/resource_management_guide/sec-cpu ):
	 * 	cpu/cpu.cfs_period_us	the period for which CPU resources should be reallocated
	 * 	cpu/cpu.cfs_quota_us	time for which all tasks in a cgroup can run during one period
	 *
	 * cgroups v2 (see https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/8/html/managing_monitoring_and_updating_the_kernel/using-cgroups-v2-to-control-distribution-of-cpu-time-for-applications_managing-monitoring-and-updating-the-kernel ):
	 *	cpu.max			containing 2 numbers, quota and period, space-separated
	 */
	guint64 period = 0, quota;
	gchar *cgroup_path = find_cgroups_path("cpu");
	gchar *v1periodpath = g_strdup_printf("/sys/fs/cgroup/cpu%s/cpu.cfs_period_us", cgroup_path);
	gchar *v1quotapath = g_strdup_printf("/sys/fs/cgroup/cpu%s/cpu.cfs_quota_us", cgroup_path);
	//siril_debug_print("trying cpu controller path: %s\n", v1periodpath);
	if (!read_from_file(v1periodpath, &period) &&
			!read_from_file(v1quotapath, &quota)) {
		siril_debug_print("found cgroups v1 cpu quota %"G_GUINT64_FORMAT" and period %"G_GUINT64_FORMAT" in cgroup %s\n", quota, period, cgroup_path);
	} else {
		/* retry with cgroups v1, but without the group name in the path */
		if (!read_from_file("/sys/fs/cgroup/cpu/cpu.cfs_period_us", &period) &&
				!read_from_file("/sys/fs/cgroup/cpu/cpu.cfs_quota_us", &quota)) {
			siril_debug_print("found cgroups v1 cpu quota %"G_GUINT64_FORMAT" and period %"G_GUINT64_FORMAT" in main controller\n", quota, period);
		} else {
			/* try with cgroups v2 */
			gchar *v2path = g_strdup_printf("/sys/fs/cgroup%s/cpu.max", cgroup_path);
			if (!read_2_from_file(v2path, &quota, &period)) {
				siril_debug_print("found cgroups v2 cpu quota %"G_GUINT64_FORMAT" and period %"G_GUINT64_FORMAT" in %s\n", quota, period, v2path);
			}
			else {
				siril_debug_print("no cgroups cpu bandwidth limitations found\n");
				period = 0;	// to be sure in case of partial read above
			}
			g_free(v2path);
		}
	}
	g_free(v1periodpath);
	g_free(v1quotapath);
	g_free(cgroup_path);
	if (period != 0)
		return round_to_int((double)quota / period);
	return 0;
}
#else
int get_available_cpu_cgroups() {
	return 0;
}
#endif

/**
 * Gets available memory in bytes.
 * @return available memory in Bytes, 0 if it fails.
 */
guint64 get_available_memory() {
#if defined(__linux__) || defined(__CYGWIN__)
	static super_bool initialized_cgroups = BOOL_NOT_SET;
	static gboolean initialized_meminfo = FALSE;

	/* first, we should check for cgroups limits */
	if (BOOL_FALSE != initialized_cgroups) {
		guint64 amount;
		if (!get_available_mem_cgroups(&amount)) {
			initialized_cgroups = BOOL_TRUE;
			return amount;
		}
		initialized_cgroups = BOOL_FALSE;
	}

	/* the code below gets the available memory on the OS in /proc/meminfo. */
	static gint64 last_check_time = 0;
	static gint fd;
	static guint64 available;
	static gboolean has_available = FALSE;
	gint64 time;

	if (!initialized_meminfo) {
		fd = g_open("/proc/meminfo", O_RDONLY);
		initialized_meminfo = TRUE;
	}
	if (fd < 0) {
		siril_debug_print("/proc/meminfo is unavailable\n");
		return (guint64) 0;
	}

	time = g_get_monotonic_time();
	if (time - last_check_time >= G_TIME_SPAN_SECOND) {
		last_check_time = time;
		has_available = FALSE;

		if (lseek(fd, 0, SEEK_SET))
			return (guint64) 0;

		gchar buffer[512];
		ssize_t size = read(fd, buffer, sizeof(buffer) - 1);
		if (size <= 0)
			return (guint64) 0;

		buffer[size] = '\0';

		char *str = strstr(buffer, "MemAvailable:");
		if (!str)
			return (guint64) 0;

		char *end;
		str += 13;
		available = strtoull(str, &end, 0);
		if (end == str)
			return (guint64) 0;
		str = end;

		for (; *str; str++) {
			if (*str == 'k') {
				available <<= 10;
				break;
			} else if (*str == 'M') {
				available <<= 20;
				break;
			}
		}
		if (*str == '\0')
			return (guint64) 0;

		has_available = TRUE;
	}

	if (!has_available)
		return (guint64) 0;

	return available;
#elif defined(OS_OSX)
	guint64 mem = (guint64) 0; /* this is the default value if we can't retrieve any values */
	vm_size_t page_size;
	mach_port_t mach_port;
	mach_msg_type_number_t count;
	vm_statistics64_data_t vm_stats;

	mach_port = mach_host_self();
	count = sizeof(vm_stats) / sizeof(natural_t);
	if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
			KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
					(host_info64_t)&vm_stats, &count))	{

		gint64 unused_memory = ((gint64)vm_stats.free_count +
				(gint64)vm_stats.inactive_count +
				(gint64)vm_stats.wire_count) * (gint64) page_size;

		mem = (guint64) (unused_memory);
	}
	return mem;
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). ----------- */
	/* FIXME: this is incorrect, this returns the available amount of memory at boot, not at runtime */
	static gboolean initialized = FALSE;
	static guint64 mem = (guint64) 0; /* this is the default value if we can't retrieve any values */
	if (!initialized) {
		FILE* fp = fopen("/var/run/dmesg.boot", "r");
		if (fp != NULL) {
			size_t bufsize = 1024 * sizeof(char);
			gchar *buf = g_new(gchar, bufsize);
			long value = -1L;
			while (getline(&buf, &bufsize, fp) >= 0) {
				if (strncmp(buf, "avail memory", 12) != 0)
					continue;
				sscanf(buf, "%*s%*s%*s%ld", &value);
				break;
			}
			fclose(fp);
			g_free(buf);
			if (value != -1L)
				mem = (guint64) (value * 1024UL);
		}
		initialized = TRUE;
	}
	return mem;
#elif defined(_WIN32) /* Windows */
	guint64 mem = (guint64) 0; /* this is the default value if we can't retrieve any values */
	MEMORYSTATUSEX memStatusEx = { 0 };
	memStatusEx.dwLength = sizeof(MEMORYSTATUSEX);
	if (GlobalMemoryStatusEx(&memStatusEx)) {
		mem = (guint64) (memStatusEx.ullAvailPhys);
	}
	return mem;
#else
	fprintf(stderr, "Siril failed to get available free RAM memory\n");
	return (guint64) 0;
#endif
}

/**
 * Get max memory depending on memory management mode
 * @return return the max memory
 */
int get_max_memory_in_MB() {
	int retval;
	switch (com.pref.mem_mode) {
		default:
		case RATIO:
			retval = round_to_int(com.pref.memory_ratio *
					(double)get_available_memory() / BYTES_IN_A_MB);
			break;
		case AMOUNT:
			retval = round_to_int(com.pref.memory_amount * 1024.0);
	}
	if (sizeof(void *) == 4 && retval > 1900) {
		siril_log_message(_("Limiting processing to 1900 MiB allocations (32-bit system)\n"));
		retval = 1900;
	}
	return retval;
}

/**
 *
 * @param filename
 * @param size
 */
#ifdef _WIN32
/* stolen from gimp which in turn stole it from glib 2.35 */
gchar* get_special_folder(int csidl) {
	HRESULT hr;
	LPITEMIDLIST pidl = NULL;
	BOOL b;
	gchar *retval = NULL;

	hr = SHGetSpecialFolderLocation(NULL, csidl, &pidl);
	if (hr == S_OK) {
		wchar_t path[MAX_PATH + 1];

		b = SHGetPathFromIDListW(pidl, path);
		if (b)
			retval = g_utf16_to_utf8(path, -1, NULL, NULL, NULL);
		CoTaskMemFree(pidl);
	}
	return retval;
}
#endif

/**
 * Check how many files a process can have open and try to extend the limit if possible.
 * The max files depends of the Operating System and of cfitsio (NMAXFILES)
 * @param nb_frames number of file processed
 * @param nb_allowed_file the maximum of file that can be opened
 * @return TRUE if the system can open all the files, FALSE otherwise
 */
gboolean allow_to_open_files(int nb_frames, int *nb_allowed_file) {
	int open_max, maxfile, MAX_NO_FILE_CFITSIO, MAX_NO_FILE;
	float version;

	/* get the limit of cfitsio */
	fits_get_version(&version);
	MAX_NO_FILE_CFITSIO = (version < 3.45f) ? 1000 : 10000;

	/* get the OS limit and extend it if possible */
#ifdef _WIN32
	MAX_NO_FILE = min(MAX_NO_FILE_CFITSIO, 2048);
	open_max = _getmaxstdio();
	if (open_max < MAX_NO_FILE) {
		/* extend the limit to 2048 if possible
		 * 2048 is the maximum on WINDOWS */
		_setmaxstdio(MAX_NO_FILE);
		open_max = _getmaxstdio();
	}
#else
	struct rlimit rlp;

/* we first set the limit to the CFITSIO limit */
	MAX_NO_FILE = MAX_NO_FILE_CFITSIO;
	if (getrlimit(RLIMIT_NOFILE, &rlp) == 0) {
		MAX_NO_FILE = (rlp.rlim_max == RLIM_INFINITY) ?
						MAX_NO_FILE_CFITSIO : rlp.rlim_max;

		if (rlp.rlim_cur != RLIM_INFINITY) {
			open_max = rlp.rlim_cur;
			MAX_NO_FILE = min(MAX_NO_FILE_CFITSIO, MAX_NO_FILE);
			if (open_max < MAX_NO_FILE) {
				rlp.rlim_cur = MAX_NO_FILE;
				/* extend the limit to NMAXFILES if possible */
				int retval = setrlimit(RLIMIT_NOFILE, &rlp);
				if (!retval) {
					getrlimit(RLIMIT_NOFILE, &rlp);
					open_max = rlp.rlim_cur;
				}
			}
		} else { // no soft limits
			open_max = MAX_NO_FILE;
		}
	} else {
		open_max = sysconf(_SC_OPEN_MAX); // if no success with getrlimit, try with sysconf
	}
#endif // _WIN32

	maxfile = min(open_max, MAX_NO_FILE);
	siril_debug_print("Maximum of files that will be opened=%d\n", maxfile);
	*nb_allowed_file = maxfile;

	return nb_frames < maxfile;
}

SirilWidget *siril_file_chooser_open(GtkWindow *parent, GtkFileChooserAction action) {
	gchar *title;
	SirilWidget *w;
	if (action == GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER) {
		title = g_strdup(_("Select Folder"));
	} else {
		title = g_strdup(_("Open File"));
	}
	w = gtk_file_chooser_dialog_new(title, parent, action, _("_Cancel"),
			GTK_RESPONSE_CANCEL, _("_Open"), GTK_RESPONSE_ACCEPT,
			NULL);
	g_free(title);
	return w;
}

SirilWidget *siril_file_chooser_add(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Add Files"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Add"), GTK_RESPONSE_ACCEPT,
			NULL);
}

SirilWidget *siril_file_chooser_save(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Save File"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);
}

gint siril_dialog_run(SirilWidget *widgetdialog) {
	return gtk_dialog_run(GTK_DIALOG(GTK_FILE_CHOOSER(widgetdialog)));
}

void siril_widget_destroy(SirilWidget *widgetdialog) {
	gtk_widget_destroy(widgetdialog);
}

GInputStream *siril_input_stream_from_stdin() {
	GInputStream *input_stream = NULL;
#ifdef _WIN32
	HANDLE handle = GetStdHandle(STD_INPUT_HANDLE);

	if (handle == INVALID_HANDLE_VALUE) {
		gchar *emsg = g_win32_error_message(GetLastError());
		g_printerr("Unable to acquire HANDLE for STDIN: %s\n", emsg);
		g_free(emsg);
		return NULL;
	}
	input_stream = g_win32_input_stream_new(handle, FALSE);
#else
	input_stream = g_unix_input_stream_new(fileno(stdin), FALSE);
#endif
	return input_stream;
}

#ifdef _WIN32
/* origin of sources: https://stackoverflow.com/questions/24171017/win32-console-application-that-can-open-windows */
int ReconnectIO(int OpenNewConsole) {
	int hConHandle;
	HANDLE lStdHandle;
	FILE *fp;
	int MadeConsole;

	MadeConsole = 0;
	if (!AttachConsole(ATTACH_PARENT_PROCESS)) {
		if (!OpenNewConsole)
			return 0;

		MadeConsole = 1;
		if (!AllocConsole())
			return 0;
	}

	// STDOUT to the console
	lStdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "w");
	*stdout = *fp;
	setvbuf( stdout, NULL, _IONBF, 0);

	// STDIN to the console
	lStdHandle = GetStdHandle(STD_INPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "r");
	*stdin = *fp;
	setvbuf( stdin, NULL, _IONBF, 0);

	// STDERR to the console
	lStdHandle = GetStdHandle(STD_ERROR_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "w");
	*stderr = *fp;
	setvbuf(stderr, NULL, _IONBF, 0);

	return MadeConsole;
}

char* siril_real_path(const char *source) {
	HANDLE hFile;
	DWORD maxchar = 2048;

	wchar_t *wsource = g_utf8_to_utf16(source, -1, NULL, NULL, NULL);
	if ( wsource == NULL ) {
		return NULL ;
	}

	if (!(GetFileAttributesW(wsource) & FILE_ATTRIBUTE_REPARSE_POINT)) { /* Ce n'est pas un lien symbolique , je sors */
		g_free(wsource);
		return NULL;
	}

	wchar_t *wFilePath = g_new(wchar_t, maxchar + 1);
	if (!wFilePath) {
		PRINT_ALLOC_ERR;
		g_free(wsource);
		return NULL;
	}
	wFilePath[0] = 0;

	hFile = CreateFileW(wsource, GENERIC_READ, FILE_SHARE_READ, NULL,
			OPEN_EXISTING, 0, NULL);
	if (hFile == INVALID_HANDLE_VALUE) {
		g_free(wFilePath);
		g_free(wsource);
		return NULL;
	}

	GetFinalPathNameByHandleW(hFile, wFilePath, maxchar, 0);

	gchar *gFilePath = g_utf16_to_utf8(wFilePath + 4, -1, NULL, NULL, NULL); // +4 = enleve les 4 caracteres du prefixe "//?/"
	g_free(wsource);
	g_free(wFilePath);
	CloseHandle(hFile);
	return gFilePath;
}
#endif

// for debug purposes
void log_used_mem(gchar *when) {
	guint64 used = get_used_RAM_memory();
	gchar *mem = g_format_size_full(used, G_FORMAT_SIZE_IEC_UNITS);
	siril_debug_print("Used memory %s: %s\n", when, mem);
	g_free(mem);
}
