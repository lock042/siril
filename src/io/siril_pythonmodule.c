// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include <glib.h>
#include <gio/gio.h>
#ifdef _WIN32
#include <windows.h>
#include <gio/gwin32inputstream.h>
#include "core/OS_utils.h"
#else
#include <gio/gunixinputstream.h>
#include <glib-unix.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#endif
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <limits.h>
#include <time.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/siril_update.h"
#include "core/siril_app_dirs.h"
#include "algos/siril_random.h"
#include "algos/background_extraction.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/siril_pythoncommands.h"
#include "io/siril_pythonmodule.h"
#include "io/siril_plot.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/siril_plot.h"
#include "gui/script_menu.h"
#include "gui/utils.h"

// 65k buffer is enough for any object except pixel data and things
// that could be an arbitrary length. For pixel data, FITS header,
// history, ICC profile and unknown keys we use a shared memory region
// (and really it's massively conservative to use shm for FITS
// header / keys and  history, and probably for 99% of ICC profiles
// anyone is likely to use with Siril, though there are some large
// ones such as scanner profiles)...
#define BUFFER_SIZE 65536
#define MAX_RETRIES 3
#define PIPE_NAME "\\\\.\\pipe\\mypipe"
#define SOCKET_PORT 12345

#ifdef _WIN32
#define PYTHON_EXE "python.exe"
#else
#define PYTHON_EXE "python3"
#endif

#define MODULE_DIR "python_module"
// Statics

// Forward declarations
static void cleanup_connection(Connection *conn);
static gboolean wait_for_client(Connection *conn);
static gboolean handle_client_communication(Connection *conn);

gboolean send_response(Connection* conn, uint8_t status, const void* data, uint32_t length) {
	ResponseHeader header = {
		.status = status,
		.length = GUINT32_TO_BE(length)  // Convert to network byte order
	};
#ifdef _WIN32
    DWORD bytes_written = 0;

    // Allocate a single buffer with header and data
    size_t total_size = sizeof(header) + (data && length > 0 ? length : 0);
    void* combined_buffer = malloc(total_size);
    if (!combined_buffer) {
        siril_log_message("Memory allocation failed for combined write\n");
        return FALSE;
    }

    // Copy header to the start of the buffer
    memcpy(combined_buffer, &header, sizeof(header));

    // Copy data after the header if present
    if (data && length > 0) {
        memcpy((char*)combined_buffer + sizeof(header), data, length);
    }

    // Single WriteFile call for atomic transfer
    if (!WriteFile(conn->pipe_handle, combined_buffer, total_size, &bytes_written, NULL) ||
        bytes_written != total_size) {
        siril_log_message("Failed to send response: %lu\n", GetLastError());
        free(combined_buffer);
        return FALSE;
    }

    free(combined_buffer);
#else
	ssize_t bytes_written;

	// Send header
	bytes_written = write(conn->client_fd, &header, sizeof(header));
	if (bytes_written != sizeof(header)) {
		siril_debug_print("Failed to send response header: %s\n", g_strerror(errno));
		return FALSE;
	}

	// Send data if present
	if (data && length > 0) {
		bytes_written = write(conn->client_fd, data, length);
		if (bytes_written != length) {
			siril_debug_print("Failed to send response data: %s\n", g_strerror(errno));
			return FALSE;
		}
	}
#endif
	return TRUE;
}

#ifdef _WIN32
static gboolean create_shared_memory_win32(const char* name, size_t size, win_shm_handle_t* handle) {
    handle->mapping = CreateFileMapping(
        INVALID_HANDLE_VALUE,    // Use paging file
        NULL,                    // Default security
        PAGE_READWRITE,          // Read/write access
        (DWORD)(size >> 32),     // High-order DWORD of size
        (DWORD)(size & 0xFFFFFFFF), // Low-order DWORD of size
        name);                   // Name of mapping object

    if (handle->mapping == NULL) {
        printf("Failed to create file mapping: %lu\n", GetLastError());
        return FALSE;
    }

    handle->ptr = MapViewOfFile(
        handle->mapping,         // Handle to mapping object
        FILE_MAP_ALL_ACCESS,     // Read/write permission
        0,                       // Offset high
        0,                       // Offset low
        size);           // Number of bytes to map

    if (handle->ptr == NULL) {
        CloseHandle(handle->mapping);
        printf("Failed to map view of file: %lu\n", GetLastError());
        return FALSE;
    }
    return TRUE;
}

gboolean siril_allocate_shm(void** shm_ptr_ptr,
							char* shm_name_ptr,
							size_t total_bytes,
							win_shm_handle_t *win_handle) {
	void *shm_ptr = NULL;
	snprintf(shm_name_ptr, 30, "%08x%08x%08x%04x", siril_random_int(), siril_random_int(), siril_random_int(), siril_random_int());
	*win_handle = (win_shm_handle_t){ NULL, NULL };
	if (!create_shared_memory_win32(shm_name_ptr, total_bytes, win_handle)) {
		return FALSE;
	}
	shm_ptr = win_handle->ptr;
	*shm_ptr_ptr = shm_ptr;
	return TRUE;
}

#else

gboolean siril_allocate_shm(void** shm_ptr_ptr,
							char* shm_name_ptr,
							size_t total_bytes,
							int *fd) {

	void *shm_ptr = NULL;
	snprintf(shm_name_ptr, 30, "/%08x%08x%08x%04x", siril_random_int(), siril_random_int(), siril_random_int(), siril_random_int());
	*fd = shm_open(shm_name_ptr, O_CREAT | O_RDWR | O_EXCL, S_IRUSR | S_IWUSR);
	if (*fd == -1) {
		siril_log_color_message(_("Invalid file descriptor after shm_open: %s\n"), "red", strerror(errno));
		return FALSE;
	}

	if (*fd < 0) {
		siril_log_color_message(_("Invalid file descriptor after shm_open\n"), "red");
		shm_unlink(shm_name_ptr);
		return FALSE;
	}

	// Round total_bytes up to page size
	long page_size = sysconf(_SC_PAGESIZE);
	if (page_size <= 0) {
		siril_log_color_message(_("Invalid page size reported\n"), "red");
		shm_unlink(shm_name_ptr);
		return FALSE;
	}
	off_t aligned_size = (total_bytes + page_size - 1) & ~(page_size - 1);
	printf("SHM allocation: Original size: %zu, Aligned size: %zu, Page size: %ld\n",
					  total_bytes, aligned_size, page_size);

	siril_debug_print("Truncating shm file to %lu bytes\n", total_bytes);

	// Truncate to ensure exact size
	if (ftruncate(*fd, aligned_size) == -1) {
		siril_log_color_message(_("Failed to set shared memory size (total_bytes: %ld): %s\n"), "red", aligned_size, strerror(errno));
		close(*fd);
		shm_unlink(shm_name_ptr);
		return FALSE;
	}

	// then mmap
	shm_ptr = mmap(NULL, (size_t) aligned_size, PROT_READ | PROT_WRITE,
				MAP_SHARED, *fd, 0);
	if (shm_ptr == MAP_FAILED) {
		siril_log_color_message(_("Failed to map shared memory: %s\n"), "red", strerror(errno));
		close(*fd);
		shm_unlink(shm_name_ptr);
		return FALSE;
	}

	*shm_ptr_ptr = shm_ptr;
	return TRUE;
}
#endif

// Initialize tracking system
static void init_shm_tracking(Connection *conn) {
	if (!conn->g_shm_initialized) {
		g_mutex_init(&conn->g_shm_mutex);
		conn->g_shm_initialized = TRUE;
	}
}

// Track new allocation
#ifdef _WIN32
static void track_shm_allocation(Connection *conn, const char* shm_name, void* shm_ptr, size_t size, win_shm_handle_t* handle) {
#else
static void track_shm_allocation(Connection *conn, const char* shm_name, void* shm_ptr, size_t size, int fd) {
#endif
	if (!conn->g_shm_initialized) init_shm_tracking(conn);

	// Using g_new0 and copying 255 chars ensures the last byte is 0
	shm_allocation_t* allocation = g_new0(shm_allocation_t, 1);
	memcpy(allocation->shm_name, shm_name, 255);

	allocation->shm_ptr = shm_ptr;
	allocation->size = size;
#ifdef _WIN32
	allocation->handle = *handle;
#else
	allocation->fd = fd;
#endif
	g_mutex_lock(&conn->g_shm_mutex);
	conn->g_shm_allocations = g_slist_append(conn->g_shm_allocations, allocation);
	g_mutex_unlock(&conn->g_shm_mutex);
}

// Helper function to match allocation by name
static gint find_allocation_by_name(gconstpointer a, gconstpointer b) {
	const shm_allocation_t* allocation = a;
	const char* name = b;
	return strcmp(allocation->shm_name, name);
}

// Cleanup allocation
void cleanup_shm_allocation(Connection *conn, const char* shm_name) {
	if (!conn->g_shm_initialized) return;

	g_mutex_lock(&conn->g_shm_mutex);

	GSList* link = g_slist_find_custom(conn->g_shm_allocations, shm_name, find_allocation_by_name);
	if (link) {
		shm_allocation_t* allocation = link->data;
#ifdef _WIN32
		UnmapViewOfFile(allocation->shm_ptr);
		CloseHandle(allocation->handle.mapping);
#else
		munmap(allocation->shm_ptr, allocation->size);
		close(allocation->fd);
		shm_unlink(allocation->shm_name);
#endif
		conn->g_shm_allocations = g_slist_remove(conn->g_shm_allocations, allocation);
		g_free(allocation);
	}

	g_mutex_unlock(&conn->g_shm_mutex);
}

static void cleanup_all_shm_allocations(Connection *conn) {
	if (!conn->g_shm_initialized) return;
	g_mutex_lock(&conn->g_shm_mutex);
	GSList* current = conn->g_shm_allocations;
	while (current) {
		shm_allocation_t* allocation = current->data;
#ifdef _WIN32
		UnmapViewOfFile(allocation->shm_ptr);
		CloseHandle(allocation->handle.mapping);
#else
		munmap(allocation->shm_ptr, allocation->size);
		close(allocation->fd);
		shm_unlink(allocation->shm_name);
#endif
		g_free(allocation);
		current = current->next;
	}
	g_slist_free(conn->g_shm_allocations);
	conn->g_shm_allocations = NULL;
	g_mutex_unlock(&conn->g_shm_mutex);
}

void cleanup_shm_resources(Connection *conn) {
	cleanup_all_shm_allocations(conn);
	if (conn->g_shm_initialized) {
		g_mutex_clear(&conn->g_shm_mutex);
		conn->g_shm_initialized = FALSE;
	}
}

// Handle a request for pixel data. We record the allocated SHM
// but leave clearup for another command
shared_memory_info_t* handle_pixeldata_request(Connection *conn, fits *fit, rectangle region) {
	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		const char* error_msg = _("Failed to retrieve pixel data - no image loaded");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}

	// Calculate total size of pixel data
	size_t total_bytes, row_bytes;
	if (fit->type == DATA_FLOAT) {
		row_bytes = region.w * sizeof(float);
		total_bytes = row_bytes * fit->naxes[2] * region.h;
	} else {
		row_bytes = region.w * sizeof(WORD);
		total_bytes = row_bytes * fit->naxes[2] * region.h;
	}

	// Generate unique name for shared memory and allocate it
	void* shm_ptr = NULL;
	char shm_name[256];
#ifdef _WIN32
	win_shm_handle_t win_handle;
	if (!siril_allocate_shm(&shm_ptr, shm_name, total_bytes, &win_handle)) {
		const char* error_msg = _("Failed to allocate shared memory");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}
#else
	char *shm_name_ptr = shm_name;
	int fd;
	if (!siril_allocate_shm(&shm_ptr, shm_name_ptr, total_bytes, &fd)) {
		const char* error_msg = _("Failed to allocate shared memory");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}
#endif

	// null check before memcpy
	if (shm_ptr == NULL) {
		const char* error_msg = _("Failed to allocate shared memory");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}

	// Copy data from gfit to shared memory
	size_t index = 0;
	int top = region.y + region.h;
	if (fit->type == DATA_FLOAT) {
		for (int chan = 0 ; chan < fit->naxes[2] ; chan++) {
			for (int i = region.y ; i < top ; i++) {
				memcpy((char*)shm_ptr + index, fit->fpdata[chan] + (i * fit->rx + region.x), region.w * sizeof(float));
				index += row_bytes;
			}
		}
	} else {
		for (int chan = 0 ; chan < fit->naxes[2] ; chan++) {
			for (int i = region.y ; i < top ; i++) {
				memcpy((char*)shm_ptr + index, fit->pdata[chan] + (i * fit->rx + region.x), region.w * sizeof(WORD));
				index += row_bytes;
			}
		}
	}

	// Track this allocation
	// Both sides have to munmap, close and unlink the shm for it to be deallocated entirely.
	// We can't do this immmediately as we don't know when the python side has opened it, so we
	// track it and wait for a CMD_RELEASE_SHM with the shm_name to release
#ifdef _WIN32
	track_shm_allocation(conn, shm_name, shm_ptr, total_bytes, &win_handle);
#else
	track_shm_allocation(conn, shm_name, shm_ptr, total_bytes, fd);
#endif

	// Prepare shared memory info structure
	shared_memory_info_t *info = calloc(1, sizeof(shared_memory_info_t));
	*info = (shared_memory_info_t) {
		.size = total_bytes,
		.data_type = (fit->type == DATA_FLOAT) ? 1 : 0,
		.width = region.w,
		.height = region.h,
		.channels = fit->naxes[2]
	};
	memcpy(info->shm_name, shm_name, strlen(shm_name)); // Safe copy with implicit null termination

	return info;
}

shared_memory_info_t* handle_rawdata_request(Connection *conn, void* data, size_t total_bytes) {
	if (total_bytes == 0) {
		const char* error_msg = _("Incorrect memory region specification");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}

	// Generate unique name for shared memory and allocate it
	void* shm_ptr = NULL;
	char shm_name[256];
#ifdef _WIN32
	win_shm_handle_t win_handle = { NULL, NULL };
	if (!siril_allocate_shm(&shm_ptr, shm_name, total_bytes, &win_handle)) {
		const char* error_msg = _("Failed to allocate shared memory");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}
#else
	int fd;
	if (!siril_allocate_shm(&shm_ptr, shm_name, total_bytes, &fd)) {
		const char* error_msg = _("Failed to allocate shared memory");
		if (send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}
#endif

	// null check before memcpy
	if (shm_ptr == NULL) {
		const char* error_msg = _("Failed to allocate shared memory");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return NULL;
	}

	// Copy data to shared memory
	if (data)
		memcpy(shm_ptr, data, total_bytes);
	else // for where Siril just requests the shm, provide it zeroed
		memset(shm_ptr, 0, total_bytes);

	// Track this allocation
#ifdef _WIN32
	track_shm_allocation(conn, shm_name, shm_ptr, total_bytes, &win_handle);
#else
	track_shm_allocation(conn, shm_name, shm_ptr, total_bytes, fd);
#endif

	// Prepare shared memory info structure
	shared_memory_info_t *info = calloc(1, sizeof(shared_memory_info_t));
	*info = (shared_memory_info_t) {
		.size = total_bytes,
		.data_type = 0,
		.width = 0,
		.height = 0,
		.channels = 0
	};
	memcpy(info->shm_name, shm_name, strlen(shm_name));

	// Send shared memory info to Python
	return info;
}

gboolean handle_set_pixeldata_request(Connection *conn, fits *fit, const char* payload, size_t payload_length) {
	if (!conn->thread_claimed) {
		const char* error_msg = _("Processing thread is not claimed: unable to update the current image. "
								"This is a script error: claim_thread() has either not been called or has failed, or "
								"the thread has been released too early.");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return FALSE;
	}

	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		const char* error_msg = _("No image or sequence loaded: set_pixel_data() can only be used to update a loaded image, not to create a new one");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return FALSE;
	}

	if (payload_length != sizeof(incoming_image_info_t)) {
		const char* error_msg = _("Invalid image info size");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return FALSE;
	}

	incoming_image_info_t* info = (incoming_image_info_t*)payload;
	info->width = GUINT32_FROM_BE(info->width);
	info->height = GUINT32_FROM_BE(info->height);
	info->channels = GUINT32_FROM_BE(info->channels);
	info->size = GUINT64_FROM_BE(info->size);
	if (info->size > get_available_memory() / 2) {
		const char* error_msg = _("Invalid image size: exceeds memory limit");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return FALSE;
	}
	info->data_type = GUINT32_FROM_BE(info->data_type);
	// Validate image dimensions and format
	if (info->width == 0 || info->height == 0 || info->channels == 0 ||
		info->channels > 3 || info->size == 0) {
		gchar size_str[32];
		g_snprintf(size_str, sizeof(size_str), "%" G_GUINT64_FORMAT, info->size);
		gchar *error_msg = g_strdup_printf(_("Invalid image dimensions or format: w = %u, h = %u, c = %u, size = %s"), info->width, info->height, info->channels, size_str);

		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		g_free(error_msg);
		return FALSE;
	}
	// Compute and sanitize ncpixels
	size_t ncpixels = info->width * info->height * info->channels;
	if (ncpixels * (info->data_type == 0 ? sizeof(WORD) : sizeof(float)) > get_available_memory() / 2) {
		const char* error_msg = _("Error: image dimensions exceed available memory");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return FALSE;
	}

	// Open shared memory
	void* shm_ptr = NULL;
	#ifdef _WIN32
		win_shm_handle_t win_handle = {NULL, NULL};
		HANDLE mapping = OpenFileMapping(FILE_MAP_READ, FALSE, info->shm_name);
		if (mapping == NULL) {
			const char* error_msg = "Failed to open shared memory mapping";
			if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
				siril_log_message("Error in send_response\n");
			return FALSE;
		}
		shm_ptr = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, info->size);
		if (shm_ptr == NULL) {
			CloseHandle(mapping);
			const char* error_msg = "Failed to map shared memory view";
			if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
				siril_log_message("Error in send_response\n");
			return FALSE;
		}
		win_handle.mapping = mapping;
		win_handle.ptr = shm_ptr;
	#else
		int fd = shm_open(info->shm_name, O_RDONLY, 0);
		if (fd == -1) {
			const char* error_msg = _("Failed to open shared memory");
			siril_debug_print("SHM ERROR: %s\n", error_msg);
			if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
				siril_log_message("Error in send_response\n");
			return FALSE;
		}
		shm_ptr = mmap(NULL, info->size, PROT_READ, MAP_SHARED, fd, 0);
		if (shm_ptr == MAP_FAILED) {
			close(fd);
			const char* error_msg = _("Failed to map shared memory");
			if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
				siril_log_message("Error in send_response\n");
			return FALSE;
		}
	#endif

	// Allocate new image buffer
	fit->type = DATA_UNSUPPORTED; // prevents a potential crash in on_drawingarea_motion_notify_event while we are messing with gfit
	free(fit->data);
	free(fit->fdata);
	fit->pdata[2] = fit->pdata[1] = fit->pdata[0] = fit->data = NULL;
	fit->fpdata[2] = fit->fpdata[1] = fit->fpdata[0] = fit->fdata = NULL;
	gboolean alloc_err = FALSE;
	if (info->data_type == 0) { // WORD data
		fit->pdata[2] = fit->pdata[1] = fit->pdata[0] = fit->data = calloc(ncpixels, sizeof(WORD));
		if (!fit->data) {
			alloc_err = TRUE;
		} else {
			for (int i = 0 ; i < info->channels ; i++) {
				fit->pdata[i] = fit->data + i * info->width * info->height;
			}
		}
	} else { // FLOAT data
		fit->fpdata[2] = fit->fpdata[1] = fit->fpdata[0] = fit->fdata = calloc(ncpixels, sizeof(float));
		if (!fit->fdata) {
			alloc_err = TRUE;
		} else {
			for (int i = 0 ; i < info->channels ; i++) {
				fit->fpdata[i] = fit->fdata + i * info->width * info->height;
			}
		}
	}
	if (alloc_err) {
		#ifdef _WIN32
			UnmapViewOfFile(shm_ptr);
			CloseHandle(win_handle.mapping);
		#else
			munmap(shm_ptr, info->size);
			close(fd);
		#endif
		const char* error_msg = _("Failed to allocate image buffer");
		if (!send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg)))
			siril_log_message("Error in send_response\n");
		return FALSE;
	}

	// Copy data from shared memory to gfit
	size_t total_bytes = ncpixels * (info->data_type == 1 ? sizeof(float) : sizeof(WORD));

	if (info->data_type == 0) {  // WORD data
		memcpy(fit->data, (char*) shm_ptr, total_bytes);
	} else {  // float data
		memcpy(fit->fdata, (char*) shm_ptr, total_bytes);
	}

	invalidate_stats_from_fit(fit);
	if (fit == &gfit) {
		// Update gfit metadata
		fit->type = info->data_type ? DATA_FLOAT : DATA_USHORT;
		fit->rx = fit->naxes[0] = info->width;
		fit->ry = fit->naxes[1] = info->height;
		fit->naxes[2] = info->channels;
		siril_debug_print("set_*_pixeldata: updating gfit\n");
		queue_redraw(REMAP_ALL);
	}
	// Cleanup shared memory
	#ifdef _WIN32
		UnmapViewOfFile(shm_ptr);
		CloseHandle(win_handle.mapping);
	#else
		munmap(shm_ptr, info->size);
		close(fd);
	#endif

	// In all cases we have now finished with the shm and closed and unlinked it.
	// On receipt of the response, python will also close and unlink the shm in its
	// finally: block.
	return send_response(conn, STATUS_OK, NULL, 0);
}

gboolean handle_plot_request(Connection* conn, const incoming_image_info_t* info) {
	// Open shared memory
	void* shm_ptr = NULL;
#ifdef _WIN32
	win_shm_handle_t win_handle = {NULL, NULL};
	HANDLE mapping = OpenFileMapping(FILE_MAP_READ, FALSE, info->shm_name);
	if (mapping == NULL) {
		const char* error_msg = "Failed to open shared memory mapping";
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
	shm_ptr = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, info->size);
	if (shm_ptr == NULL) {
		CloseHandle(mapping);
		const char* error_msg = "Failed to map shared memory view";
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
	win_handle.mapping = mapping;
	win_handle.ptr = shm_ptr;
#else
	int fd = shm_open(info->shm_name, O_RDONLY, 0);
	if (fd == -1) {
		const char* error_msg = _("Failed to open shared memory");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
	shm_ptr = mmap(NULL, info->size, PROT_READ, MAP_SHARED, fd, 0);
	if (shm_ptr == MAP_FAILED) {
		close(fd);
		const char* error_msg = _("Failed to map shared memory");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
#endif
	if (!shm_ptr) {
		const char* error_msg = _("Error: could not open shared memory");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	// Unpack the plot data
	siril_plot_data* plot_data = unpack_plot_data((uint8_t*) shm_ptr, info->size);

	// Cleanup shared memory
#ifdef _WIN32
	UnmapViewOfFile(shm_ptr);
	CloseHandle(win_handle.mapping);
#else
	munmap(shm_ptr, info->size);
	close(fd);
#endif

	// Plot the data in a siril_plot_window
	if (plot_data)
		siril_add_idle(create_new_siril_plot_window, plot_data);
	return send_response(conn, STATUS_OK, NULL, 0);
}

gboolean handle_set_bgsamples_request(Connection* conn, const incoming_image_info_t* info, gboolean show_samples, gboolean recalculate) {
	// Open shared memory
	void* shm_ptr = NULL;
#ifdef _WIN32
	win_shm_handle_t win_handle = {NULL, NULL};
	HANDLE mapping = OpenFileMapping(FILE_MAP_READ, FALSE, info->shm_name);
	if (mapping == NULL) {
		const char* error_msg = "Failed to open shared memory mapping";
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
	shm_ptr = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, info->size);
	if (shm_ptr == NULL) {
		CloseHandle(mapping);
		const char* error_msg = "Failed to map shared memory view";
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
	win_handle.mapping = mapping;
	win_handle.ptr = shm_ptr;
#else
	int fd = shm_open(info->shm_name, O_RDONLY, 0);
	if (fd == -1) {
		const char* error_msg = _("Failed to open shared memory");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
	shm_ptr = mmap(NULL, info->size, PROT_READ, MAP_SHARED, fd, 0);
	if (shm_ptr == MAP_FAILED) {
		close(fd);
		const char* error_msg = _("Failed to map shared memory");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}
#endif
	if (!shm_ptr) {
		const char* error_msg = _("Error: could not open shared memory");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	// Unpack the plot data
	size_t total_bytes = info->size;
	size_t nb_samples = total_bytes / sizeof(background_sample);
	background_sample* samples = (background_sample*) shm_ptr;

	// Build a list of the coords
	GSList *pts = NULL;

	// Reset the list in com.
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;

	// If we need to recalculate, build samples from the positions and call
	// add_background_samples
	if (recalculate) {
		for (int i = 0 ; i < nb_samples ; i++) {
			point* p = malloc(sizeof(point));
			memcpy(p, &samples[i].position, sizeof(point));
			pts = g_slist_append(pts, p);
		}
		com.grad_samples = add_background_samples(com.grad_samples, &gfit, pts);
	}
	// Otherwise, create copies of each individual sample and add them directly
	// to com.grad_samples
	else {
		for (int i = 0 ; i < nb_samples ; i++) {
			background_sample *s = malloc(sizeof(background_sample));
			memcpy(s, (background_sample*) (samples + i), sizeof(background_sample));
			// protect against zero sample size
			s->size = s->size ? s->size : SAMPLE_SIZE;
			com.grad_samples = g_slist_append(com.grad_samples, s);
		}
	}

	// Redraw if necessary
	if (show_samples && !com.headless) {
		queue_redraw(REDRAW_OVERLAY);
	}

	// Free the positions list
	g_slist_free_full(pts, free);

	// Cleanup shared memory
	#ifdef _WIN32
	UnmapViewOfFile(shm_ptr);
	CloseHandle(win_handle.mapping);
	#else
	munmap(shm_ptr, info->size);
	close(fd);
	#endif

	return send_response(conn, STATUS_OK, NULL, 0);
}

// Monitor stdout stream
static gpointer monitor_stream_stdout(GDataInputStream *data_input) {

	gsize length = 0;
	gchar *buffer;
	GError *error = NULL;

	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length, NULL, &error))) {
		siril_log_message("%s\n", buffer);
		g_free(buffer);
	}

	if (error) {
		siril_log_color_message(_("Error reading stdout: %s\n"), "red", error->message);
		g_error_free(error);
	}

	g_object_unref(data_input);
	return NULL;
}

// Monitor stderr stream
static gpointer monitor_stream_stderr(GDataInputStream *data_input) {

	gsize length = 0;
	gchar *buffer;
	GError *error = NULL;

	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length, NULL, &error))) {
//#ifdef __APPLE__
//		if (!g_strrstr(buffer, "resource_tracker"))
//#endif
			siril_log_color_message("%s\n", "red", buffer);
		g_free(buffer);
	}

	if (error) {
		siril_log_color_message(_("Error reading stderr: %s\n"), "red", error->message);
		g_error_free(error);
	}

	g_object_unref(data_input);
	return NULL;
}

static void cleanup_child_process(GPid pid, gint status, gpointer user_data) {
	Connection *conn = (Connection*) user_data;
	// Log the Python process exit status if needed
#ifdef G_OS_WIN32
	if (status == 0) {
		siril_debug_print("Python process (PID: %d) exited normally\n", pid);
	} else {
		siril_log_color_message(_("Python process (PID: %d) exited with status %d\n"), "salmon",
			pid, status);
	}
#else
	if (WIFEXITED(status)) {
		if (WEXITSTATUS(status) == 0)
			siril_debug_print("Python process (PID: %d) exited normally\n", pid);
		else
			siril_log_color_message(_("Python process (PID: %d) exited with status %d\n"), "salmon",
				pid, WEXITSTATUS(status));
	} else if (WIFSIGNALED(status)) {
		siril_log_color_message(_("Python process (PID: %d) terminated by signal %d\n"), "salmon",
				pid, WTERMSIG(status));
	}
#endif

	// If we had the python thread lock and failed to release it, release it now
	// and reset the busy cursor
	if (conn->thread_claimed) {
		com.python_claims_thread = FALSE;
		set_cursor_waiting(FALSE);
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	}

	// Clean up shared memory resources
	cleanup_shm_resources(conn);

	// Clean up the Connection
	free(conn);

	// Close the process handle
	g_spawn_close_pid(pid);
	remove_child_from_children(pid);

	// Check if it's OK to clear com.python_script
	check_python_flag();

	// Re-enable widgets
	gui_function(script_widgets_idle, NULL);
}

// Function to get Python version from venv (could be called during venv activation)
gchar* get_venv_python_version(const gchar* venv_path) {
	gchar* python_path;
#ifdef _WIN32
	python_path = g_build_filename(venv_path, "Scripts", PYTHON_EXE, NULL);
#else
	python_path = g_build_filename(venv_path, "bin", PYTHON_EXE, NULL);
#endif

	gchar* argv[] = { python_path, "-c", "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')", NULL };
	gchar* output = NULL;
	gint status;

	g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
				NULL, NULL, &output, NULL, &status, NULL);

	g_free(python_path);

	if (output) {
		g_strchomp(output);  // Remove trailing newline
		return output;
	}
	return NULL;
}

#ifdef _WIN32
static Connection* create_connection(const gchar *pipe_name) {
	Connection *conn = g_new0(Connection, 1);
	g_mutex_init(&conn->mutex);
	g_cond_init(&conn->condition);

	conn->pipe_handle = CreateNamedPipe(
		pipe_name,
		PIPE_ACCESS_DUPLEX,
		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
		1,
		BUFFER_SIZE,
		BUFFER_SIZE,
		0,
		NULL
	);

	if (conn->pipe_handle == INVALID_HANDLE_VALUE) {
		siril_debug_print("Failed to create pipe: %lu\n", GetLastError());
		g_free(conn);
		return NULL;
	}

	return conn;
}

static gboolean wait_for_client(Connection *conn) {
	g_mutex_lock(&conn->mutex);
	conn->is_connected = FALSE;
	g_mutex_unlock(&conn->mutex);

	BOOL result = ConnectNamedPipe(conn->pipe_handle, NULL);
	if (!result && GetLastError() != ERROR_PIPE_CONNECTED) {
		siril_debug_print("Failed to connect to client: %lu\n", GetLastError());
		return FALSE;
	}

	g_mutex_lock(&conn->mutex);
	conn->is_connected = TRUE;
	if (conn->client_connected_callback) {
		conn->client_connected_callback(conn->user_data);
	}
	g_mutex_unlock(&conn->mutex);

	return TRUE;
}

#else  // POSIX implementation

static Connection* create_connection(const gchar *socket_path) {
	Connection *conn = g_new0(Connection, 1);
	g_mutex_init(&conn->mutex);
	g_cond_init(&conn->condition);

	conn->server_fd = socket(AF_UNIX, SOCK_STREAM, 0);
	if (conn->server_fd == -1) {
		siril_debug_print("Failed to create socket: %s\n", g_strerror(errno));
		g_free(conn);
		return NULL;
	}

	struct sockaddr_un addr;
	memset(&addr, 0, sizeof(addr));
	addr.sun_family = AF_UNIX;
	strncpy(addr.sun_path, socket_path, sizeof(addr.sun_path) - 1);

	unlink(socket_path);  // Remove existing socket file if it exists

	if (bind(conn->server_fd, (struct sockaddr*)&addr, sizeof(addr)) == -1) {
		siril_debug_print("Failed to bind socket: %s\n", g_strerror(errno));
		close(conn->server_fd);
		g_free(conn);
		return NULL;
	}

	if (listen(conn->server_fd, 1) == -1) {
		siril_debug_print("Failed to listen on socket: %s\n", g_strerror(errno));
		close(conn->server_fd);
		unlink(socket_path);
		g_free(conn);
		return NULL;
	}

	conn->socket_path = g_strdup(socket_path);
	return conn;
}

static gboolean wait_for_client(Connection *conn) {
	g_mutex_lock(&conn->mutex);
	conn->is_connected = FALSE;
	g_mutex_unlock(&conn->mutex);

	conn->client_fd = accept(conn->server_fd, NULL, NULL);
	if (conn->client_fd == -1) {
		siril_debug_print("Failed to accept connection: %s\n", g_strerror(errno));
		return FALSE;
	}

	g_mutex_lock(&conn->mutex);
	conn->is_connected = TRUE;
	if (conn->client_connected_callback) {
		conn->client_connected_callback(conn->user_data);
	}
	g_mutex_unlock(&conn->mutex);

	return TRUE;
}

#endif  // _WIN32

/**
* handle_client_communication function with message processing
*/
static gboolean handle_client_communication(Connection *conn) {
	gchar buffer[BUFFER_SIZE];

#ifdef _WIN32
	DWORD bytes_read = 0;

	while (!conn->should_stop) {
		if (!ReadFile(conn->pipe_handle, buffer, BUFFER_SIZE, &bytes_read, NULL) ||
			bytes_read == 0) {
			g_mutex_lock(&conn->mutex);
			conn->is_connected = FALSE;
			if (conn->client_disconnected_callback) {
				conn->client_disconnected_callback(conn->user_data);
			}
			g_mutex_unlock(&conn->mutex);

			DisconnectNamedPipe(conn->pipe_handle);
			return TRUE;
		}

		process_connection(conn, buffer, bytes_read);
	}

#else
	gssize bytes_read;

	while (!conn->should_stop) {
		bytes_read = read(conn->client_fd, buffer, BUFFER_SIZE - 1);
		buffer[BUFFER_SIZE - 1] = '\0'; // Explicitly NULL terminate the buffer for safety

		if (bytes_read <= 0) {
			if (bytes_read < 0) {
				siril_debug_print("Error reading from socket: %s\n", g_strerror(errno));
			}

			g_mutex_lock(&conn->mutex);
			conn->is_connected = FALSE;
			if (conn->client_disconnected_callback) {
				conn->client_disconnected_callback(conn->user_data);
			}
			g_mutex_unlock(&conn->mutex);

			close(conn->client_fd);
			conn->client_fd = -1;
			return TRUE;
		}
		process_connection(conn, buffer, bytes_read);
	}
#endif

	return !conn->should_stop;
}

static void cleanup_connection(Connection *conn) {
	if (!conn) return;

	g_mutex_lock(&conn->mutex);
	conn->should_stop = TRUE;
	conn->is_connected = FALSE;
	g_mutex_unlock(&conn->mutex);

#ifdef _WIN32
	if (conn->pipe_handle != INVALID_HANDLE_VALUE) {
		CloseHandle(conn->pipe_handle);
	}
#else
	if (conn->client_fd > 0) {
		close(conn->client_fd);
	}
	if (conn->server_fd > 0) {
		close(conn->server_fd);
	}
	if (conn->socket_path) {
		unlink(conn->socket_path);
		g_free(conn->socket_path);
	}
#endif

	g_mutex_clear(&conn->mutex);
	g_cond_clear(&conn->condition);
	g_free(conn);
}

static gpointer connection_worker(gpointer data) {
	Connection *conn = (Connection*)data;
	siril_debug_print("Python communication initialized...\n");
	while (!conn->should_stop) {
		if (wait_for_client(conn)) {
			handle_client_communication(conn);
		} else {
			// Wait before retrying
			g_usleep(1000000);  // 1 second
		}
	}
	siril_log_message(_("Python communication worker finished...\n"));
	return NULL;
}

typedef struct {
	gchar *venv_path;
	gchar *python_version;
	GHashTable *env_vars;
} PythonVenvInfo;

static gchar* find_venv_bin_dir(const gchar *venv_path) {
	gchar *bin_dir = NULL;

#ifdef _WIN32
	// Try Scripts directory first
	bin_dir = g_build_filename(venv_path, "Scripts", NULL);
	if (!g_file_test(bin_dir, G_FILE_TEST_EXISTS)) {
		g_free(bin_dir);
		// Try bin directory as fallback
		bin_dir = g_build_filename(venv_path, "bin", NULL);
		if (!g_file_test(bin_dir, G_FILE_TEST_EXISTS)) {
			g_free(bin_dir);
			return NULL;
		}
	}
#else
	bin_dir = g_build_filename(venv_path, "bin", NULL);
	if (!g_file_test(bin_dir, G_FILE_TEST_EXISTS)) {
		g_free(bin_dir);
		return NULL;
	}
#endif

	return bin_dir;
}

static gchar* find_venv_python_exe(const gchar *venv_path, const gboolean verbose) {
	gchar *python_exe = NULL;

#ifdef _WIN32
	// Try Scripts directory first
	python_exe = g_build_filename(venv_path, "Scripts", PYTHON_EXE, NULL);
	if (!g_file_test(python_exe, G_FILE_TEST_EXISTS)) {
		g_free(python_exe);
		// Try bin directory as fallback
		python_exe = g_build_filename(venv_path, "bin", PYTHON_EXE, NULL);
		if (!g_file_test(python_exe, G_FILE_TEST_EXISTS)) {
			if (verbose) siril_debug_print("Error: python executable not found in the venv\n");
			g_free(python_exe);
			return NULL;
		}
	}
#else
	python_exe = g_build_filename(venv_path, "bin", PYTHON_EXE, NULL);
	if (!g_file_test(python_exe, G_FILE_TEST_EXISTS)) {
		if (verbose) siril_debug_print("Error: python executable not found in the venv\n");
		g_free(python_exe);
		return NULL;
	}
#endif

	return python_exe;
}

static gchar* get_setup_path(const gchar* module_dir) {
	return g_build_filename(module_dir, "pyproject.toml", NULL);
}

static version_number get_module_version(const gchar* filename, GError** error) {
	gchar* content = NULL;
	gsize length;
	gchar* version = NULL;
	version_number ver = { 0 };
	gboolean is_metadata;

	// Check if this is a METADATA file
	is_metadata = g_str_has_suffix(filename, "METADATA");

	// Read the entire file content
	if (!g_file_get_contents(filename, &content, &length, error)) {
		return ver;
	}

	// Create appropriate pattern based on file type
	const gchar* pattern = is_metadata ?
		"^Version:\\s*([^\\s]+)" :          // METADATA pattern
		"\\bversion\\s*=\\s*[\"']([^\"']+)[\"']"; // setup.py pattern

	GRegex* regex = g_regex_new(pattern,
							is_metadata ? 0 : G_REGEX_MULTILINE,
							0,
							error);
	if (regex == NULL) {
		g_free(content);
		return ver;
	}

	// Try to find the version match
	GMatchInfo* match_info;
	if (g_regex_match(regex, content, 0, &match_info)) {
		version = g_match_info_fetch(match_info, 1);
	} else {
		g_set_error(error,
				G_FILE_ERROR,
				G_FILE_ERROR_FAILED,
				"Version number not found in %s",
				is_metadata ? "METADATA" : "setup.py");
	}

	if (version) {
		ver = get_version_number_from_string(version);
	}

	// Clean up
	g_match_info_free(match_info);
	g_regex_unref(regex);
	g_free(content);
	g_free(version);

	return ver;
}

static version_number get_installed_module_version(const gchar* python_path, GError **error) {
	gchar *cmd = g_strdup(python_path);
	version_number ver = { 0 };
	gchar *stdout_data = NULL;
	gchar *stderr_data = NULL;
	gint exit_status;

	gchar *argv[] = { cmd, "-m", "pip", "show", "sirilpy", NULL };

	// Execute pip show command
	GError *spawn_error = NULL;
	if (!g_spawn_sync(NULL,  // Working directory (NULL = current)
					argv,
					NULL,  // Environment variables (NULL = inherit)
					G_SPAWN_SEARCH_PATH,
					NULL,  // Child setup function
					NULL,  // User data for child setup
					&stdout_data,
					&stderr_data,
					&exit_status,
					&spawn_error)) {
		g_set_error(error,
				G_FILE_ERROR,
				G_FILE_ERROR_FAILED,
				"Failed to execute pip: %s",
				spawn_error ? spawn_error->message : "Unknown error");
		g_clear_error(&spawn_error);
		g_free(cmd);
		return ver;
	}

	// Check if pip command was successful
	if (exit_status != 0) {
		g_set_error(error,
				G_FILE_ERROR,
				G_FILE_ERROR_FAILED,
				"pip command failed: %s",
				stderr_data ? stderr_data : "Unknown error");
		g_free(stdout_data);
		g_free(stderr_data);
		g_free(cmd);
		return ver;
	}

	// Create regex pattern for version
	GRegex *regex = g_regex_new("^Version:\\s*([^\\s]+)",
							G_REGEX_MULTILINE,
							0,
							error);
	if (regex == NULL) {
		g_free(stdout_data);
		g_free(stderr_data);
		g_free(cmd);
		return ver;
	}

	// Try to find the version match
	GMatchInfo *match_info;
	gchar *version = NULL;

	if (g_regex_match(regex, stdout_data, 0, &match_info)) {
		version = g_match_info_fetch(match_info, 1);
	} else {
		g_set_error(error,
				G_FILE_ERROR,
				G_FILE_ERROR_FAILED,
				"Version number not found in pip output");
	}

	if (version) {
		ver = get_version_number_from_string(version);
	}

	// Clean up
	g_match_info_free(match_info);
	g_regex_unref(regex);
	g_free(stdout_data);
	g_free(stderr_data);
	g_free(cmd);
	g_free(version);

	return ver;
}

gboolean copy_directory_recursive(const gchar *src_dir, const gchar *dest_dir, GError **error) {
	g_return_val_if_fail(src_dir != NULL && dest_dir != NULL, FALSE);

	GDir *dir = g_dir_open(src_dir, 0, error);
	if (!dir) {
		return FALSE;
	}

	GFile *src_file = g_file_new_for_path(src_dir);
	GFile *dest_file = g_file_new_for_path(dest_dir);
	gboolean success = TRUE;

	// Create destination directory if it doesn't exist
	if (!g_file_make_directory_with_parents(dest_file, NULL, error)) {
		if (!g_error_matches(*error, G_IO_ERROR, G_IO_ERROR_EXISTS)) {
			success = FALSE;
			goto cleanup_files;
		}
		g_clear_error(error);  // Clear the "already exists" error
	}

	const gchar *filename;
	while ((filename = g_dir_read_name(dir)) != NULL && success) {
		gchar *src_path = g_build_filename(src_dir, filename, NULL);
		gchar *dest_path = g_build_filename(dest_dir, filename, NULL);
		GFile *src_child = g_file_new_for_path(src_path);
		GFile *dest_child = g_file_new_for_path(dest_path);

		GFileType file_type = g_file_query_file_type(src_child,
													G_FILE_QUERY_INFO_NOFOLLOW_SYMLINKS,
													NULL);

		if (file_type == G_FILE_TYPE_DIRECTORY) {
			// Recursively copy subdirectory
			success = copy_directory_recursive(src_path, dest_path, error);
		} else if (file_type == G_FILE_TYPE_REGULAR) {
			// Copy file
			GFileCopyFlags flags = G_FILE_COPY_OVERWRITE |
								G_FILE_COPY_NOFOLLOW_SYMLINKS |
								G_FILE_COPY_ALL_METADATA;

			success = g_file_copy(src_child, dest_child, flags,
								NULL, NULL, NULL, error);
		}

		g_object_unref(src_child);
		g_object_unref(dest_child);
		g_free(src_path);
		g_free(dest_path);

		if (!success) {
			break;
		}
	}

cleanup_files:
	g_object_unref(src_file);
	g_object_unref(dest_file);
	g_dir_close(dir);

	return success;
}

gboolean delete_directory(const gchar *dir_path, GError **error) {
	GDir *dir;
	const gchar *name;
	gchar *full_path;
	GFile *file;
	gboolean success = TRUE;

	dir = g_dir_open(dir_path, 0, error);
	if (!dir) {
		return FALSE;
	}

	while ((name = g_dir_read_name(dir))) {
		full_path = g_build_filename(dir_path, name, NULL);

		if (g_file_test(full_path, G_FILE_TEST_IS_DIR)) {
			if (!delete_directory(full_path, error)) {
				success = FALSE;
				break;
			}
		} else {
			file = g_file_new_for_path(full_path);
			if (!g_file_delete(file, NULL, error)) {
				g_object_unref(file);
				success = FALSE;
				break;
			}
			g_object_unref(file);
		}

		g_free(full_path);
	}

	g_dir_close(dir);

	if (success) {
		file = g_file_new_for_path(dir_path);
		if (!g_file_delete(file, NULL, error)) {
			g_object_unref(file);
			return FALSE;
		}
		g_object_unref(file);
	}

	return success;
}

gboolean install_module_with_pip(const gchar* module_path, const gchar* user_module_path,
							const gchar* venv_path, GError** error) {
	g_return_val_if_fail(module_path != NULL, FALSE);
	g_return_val_if_fail(user_module_path != NULL, FALSE);
	g_return_val_if_fail(venv_path != NULL, FALSE);
	gchar *python_path = find_venv_python_exe(venv_path, TRUE);
	siril_debug_print("Python path: %s\n", python_path);
	g_return_val_if_fail(python_path != NULL, FALSE);
	gboolean needs_install = FALSE;
	gchar* module_setup_path = get_setup_path(module_path);
	GError* modver_error = NULL;
	version_number module_version = get_module_version(module_setup_path, &modver_error);
	if (modver_error) {
		g_propagate_error(error, modver_error);
		g_free(module_setup_path);
		g_free(python_path);
		return FALSE;
	}

	// Check if temp module directory exists
	if (!g_file_test(user_module_path, G_FILE_TEST_EXISTS)) {
		needs_install = TRUE;
	} else {
		GError* ver_error = NULL;
		version_number user_version = get_installed_module_version(python_path, &ver_error);
		if (ver_error) { // May just mean it's not installed: anyway we will try
			siril_debug_print("Module version check status (harmless): %s\n", ver_error->message);
			g_clear_error(&ver_error);
			needs_install = TRUE;
			// user_version is {0} so the version check will require us to install it
		}
		siril_debug_print("User version: %d.%d.%d\n", user_version.major_version, user_version.minor_version, user_version.micro_version);
		siril_debug_print("System version: %d.%d.%d\n", module_version.major_version, module_version.minor_version, module_version.micro_version);

		// Check if module version is higher than temp version
		if (compare_version(module_version, user_version) > 0) {
			needs_install = TRUE;
		} else if (compare_version(module_version, user_version) < 0) {
			// Downgrading: perhaps we have been trying a development branch of Siril
			// and are reverting to a stable branch or something. Attempt to uninstall
			// then reinstall: if not, we get drastic and yeet the entire venv then
			// rebuild it from scratch to force the downgrade
			gint uninstall_status;
			GError *uninstall_error = NULL;
			gchar *argv[] = {
				python_path,
				"-m",
				"pip",
				"uninstall",
				"-y",
				"sirilpy",
				NULL  // Array must be NULL-terminated
			};
			if (!g_spawn_sync(
					NULL,           // working_directory (NULL = inherit current)
					argv,           // argument vector
					NULL,           // inherit parent's environment
					G_SPAWN_DEFAULT, // flags
					NULL,           // child_setup function
					NULL,           // user_data for child_setup
					NULL,           // standard_output
					NULL,           // standard_error
					&uninstall_status,   // exit status
					&uninstall_error    // error
				)) {
				g_propagate_error(error, uninstall_error);
				GError* del_error = NULL;
				// Fallback, try to delete the entire venv and start again
				if (!delete_directory(user_module_path, &del_error)) {
					g_propagate_error(error, del_error);
					g_free(module_setup_path);
					g_free(python_path);
					return FALSE;
				}
			}
			needs_install = TRUE;
		}
	}

	g_free(module_setup_path);

	if (needs_install) {
		siril_log_message(_("Installing / updating python module in the background. This may take a few seconds...\n"));
		// Create user-owned directory and copy module
		if (!g_file_test(user_module_path, G_FILE_TEST_EXISTS)) {
			if (g_mkdir_with_parents(user_module_path, 0755) != 0) {
				g_set_error(error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
						"Failed to create directory: %s", user_module_path);
				g_free(python_path);
				return FALSE;
			}
		}

		GError* copy_error = NULL;
		if (!copy_directory_recursive(module_path, user_module_path, &copy_error)) {
			g_propagate_error(error, copy_error);
			g_free(python_path);
			return FALSE;
		}
		gchar *arg_module_path = g_strdup(user_module_path);
		// Install with pip
		gchar *argv[] = {
			python_path,
			"-m",
			"pip",
			"install",
			arg_module_path,
			NULL  // Array must be NULL-terminated
		};

		gint exit_status;
		GError *spawn_error = NULL;

		if (!g_spawn_sync(
				NULL,           // working_directory (NULL = inherit current)
				argv,           // argument vector
				NULL,           // inherit parent's environment
				G_SPAWN_DEFAULT, // flags
				NULL,           // child_setup function
				NULL,           // user_data for child_setup
				NULL,           // standard_output
				NULL,           // standard_error
				&exit_status,   // exit status
				&spawn_error    // error
			)) {
			g_propagate_error(error, spawn_error);
			g_free(arg_module_path);
			g_free(python_path);
			return FALSE;
		}
		g_free(arg_module_path);
		g_free(python_path);
		if (exit_status != 0) {
			g_set_error(error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
					"Pip installation failed with exit status: %d", exit_status);
			return FALSE;
		}
	}
	memcpy(&com.python_version, &module_version, sizeof(version_number));
	return TRUE;
}

static PythonVenvInfo* prepare_venv_environment(const gchar *venv_path) {
	PythonVenvInfo *info = g_new0(PythonVenvInfo, 1);
	info->venv_path = g_strdup(venv_path);
	info->env_vars = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);

	siril_log_message(_("Preparing python virtual environment: %s.\n"), venv_path);
	// Copy current environment
	gchar **current_env = g_get_environ();
	for (gchar **env = current_env; env && *env; env++) {
		gchar **parts = g_strsplit(*env, "=", 2);
		if (parts && parts[0] && parts[1]) {
			g_hash_table_insert(info->env_vars, g_strdup(parts[0]), g_strdup(parts[1]));
		}
		g_strfreev(parts);
	}
	g_strfreev(current_env);

	// Get Python version
	info->python_version = get_venv_python_version(venv_path);
	if (!info->python_version) {
		g_warning("Failed to determine Python version");
		goto cleanup;
	}

	// Set VIRTUAL_ENV
	g_hash_table_insert(info->env_vars, g_strdup("VIRTUAL_ENV"), g_strdup(venv_path));

	// Update PATH
	gchar *bin_dir = find_venv_bin_dir(venv_path);
	if (!bin_dir) {
		g_warning("Failed to locate virtual environment binary directory");
		goto cleanup;
	}
	gchar *old_path = g_hash_table_lookup(info->env_vars, "PATH");
	gchar *new_path = old_path ?
		g_strjoin(G_SEARCHPATH_SEPARATOR_S, bin_dir, old_path, NULL) :
		g_strdup(bin_dir);
	g_hash_table_insert(info->env_vars, g_strdup("PATH"), new_path);
	g_free(bin_dir);

	// Remove PYTHONPATH and PYTHONHOME to allow natural path discovery
	g_hash_table_remove(info->env_vars, "PYTHONPATH");
	g_hash_table_remove(info->env_vars, "PYTHONHOME");

	// Check the sirilpy python module is the latest version and install or
	// update it using the venv pip if not.
	siril_log_message(_("Checking the python module is up-to-date...\n"));
	gchar *module_path = g_build_filename(siril_get_system_data_dir(), MODULE_DIR, NULL);
	gchar *user_module_path = g_build_filename(g_get_user_data_dir(), "siril", ".python_module", NULL);
	GError *install_error = NULL;
	if (!install_module_with_pip(module_path, user_module_path, venv_path, &install_error)) {
		siril_log_color_message(_("Warning: unable to install or update the "
					"Siril python module.\n"), "salmon");
		g_warning("Failed to install Python module: %s",
				install_error ? install_error->message : "Unknown error");
		g_error_free(install_error);
	} else {
		siril_log_color_message(_("Python module is up-to-date\n"), "green");
	}
	g_free(module_path);

	return info;

cleanup:
	// No need to NULL check info, it can't be NULL here as it is never freed before
	// and if the g_new0 alloc fails the program will abort.
	g_free(info->venv_path);
	g_free(info->python_version);
	if (info->env_vars)
		g_hash_table_destroy(info->env_vars);
	g_free(info);
	return NULL;
}

//***********************************************************************************
// WARNING: the following function will IMMEDIATELY kill all running python scripts,
// delete the siril venv directory and rebuild it. Any call to this MUST be
// preceded by a siril_confirm_dialog(). On success a completely fresh venv will
// exist containing the current sirilpy python module, but any other modules that
// have been installed by scripts will require reinstallation.
//***********************************************************************************

void rebuild_venv() {
	gchar* venv_path = g_build_filename(g_get_user_data_dir(), "siril", "venv", NULL);
	GError *error = NULL;
	kill_all_python_scripts();
	delete_directory(venv_path, &error);
	initialize_python_venv_in_thread();
}

static gboolean check_or_create_venv(const gchar *project_path, GError **error) {

	gchar *venv_path = g_build_filename(project_path, "venv", NULL);
	siril_debug_print("venv path: %s\n", venv_path);
	gchar *python_exe = find_venv_python_exe(venv_path, FALSE);
	if (python_exe) {
		siril_debug_print("Found python executable in venv: %s\n", python_exe);
	} else {
		siril_debug_print("Did not find python executable in venv. Recreating the venv...\n");
#ifdef _WIN32
		// Check we aren't in a msys2 environment for the first init
		gchar **env = g_get_environ();
		const gchar *msys = g_environ_getenv(env, "MSYSTEM");
		g_strfreev(env);
		if (msys) {
			siril_log_color_message(_("Error: msys2 environment detected. Siril Python support cannot be correctly initialized.\n"), "red");
			siril_log_color_message(_("To complete the process, first make sure you have a Python installation (>=3.9) on your computer.\n"), "red");
			siril_log_color_message(_("Locate siril.exe (usually located in C:\\msys64\\mingw64\\bin) and start it from there.\n"), "red");
			siril_log_color_message(_("Next time you need to start siril, you can go back to starting it from msys2 terminal.\n"), "red");
			return FALSE;
		}
#endif
	}


	gboolean success = FALSE;
	GError *local_error = NULL;
	gchar *sys_python_exe = NULL;
	gchar **argv = NULL;

	// Check if venv exists
	if (!python_exe) {

#ifdef _WIN32
		gchar *bundle_python_exe = NULL;
		sys_python_exe = find_executable_in_path(PYTHON_EXE, NULL); // we want to find system python not mingw64 python
		if (sys_python_exe)
			printf("Python found in system: %s\n", sys_python_exe);
		if (!sys_python_exe) {
			const gchar *sirilrootpath = get_siril_bundle_path();
			printf("Siril bundle path: %s\n", sirilrootpath);
			bundle_python_exe = g_build_filename(sirilrootpath, "python", PYTHON_EXE, NULL);
			printf("Bundle python path: %s\n", bundle_python_exe);
			if (g_file_test(bundle_python_exe, G_FILE_TEST_IS_EXECUTABLE))
				printf("Python found in bundle: %s\n", bundle_python_exe);
			else {
				g_free(bundle_python_exe);
				bundle_python_exe = NULL;
			}
		}

		if (!sys_python_exe && !bundle_python_exe) {
			siril_log_color_message(_("No python installation found in the system or in the bundle, aborting\n"), "red");
			success = FALSE;
			goto cleanup;
		}
		if (!sys_python_exe) {
			sys_python_exe = g_strdup(bundle_python_exe);
			g_free(bundle_python_exe);
		}
		printf("Python executable: %s\n", sys_python_exe);
#else
		sys_python_exe = g_find_program_in_path(PYTHON_EXE);
#endif

		argv = g_new0(gchar*, 5);
		argv[0] = sys_python_exe;
		argv[1] = g_strdup("-m");
		argv[2] = g_strdup("venv");
		argv[3] = g_strdup(venv_path);
		argv[4] = NULL;
		siril_debug_print("Trying venv creation command: %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3]);
		gint exit_status;
		if (!g_spawn_sync(NULL, argv, NULL,
						G_SPAWN_SEARCH_PATH,
						NULL, NULL,
						NULL, NULL,
						&exit_status, &local_error)) {
			siril_log_color_message(_("Error in venv creation command: %s\n"), "red", local_error->message);
			g_propagate_error(error, local_error);
			success = FALSE;
			goto cleanup;
		}

		success = (exit_status == 0);
		if (!success) {
			g_set_error(error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
					"Failed to create virtual environment (exit status: %d)", exit_status);
			goto cleanup;
		}

cleanup:
		g_strfreev(argv);
	} else {
		success = TRUE;
		g_free(python_exe);
	}

	g_free(venv_path);
	return success;
}

gboolean python_venv_idle(gpointer user_data) {
//	g_thread_unref(com.python_init_thread);
	com.python_init_thread = NULL;
	return FALSE;
}

/*
 * Ensures that the venv is created and valid, and that the siril python module
 * is correctly installed in it along with its dependencies
 */
static gpointer initialize_python_venv(gpointer user_data) {
	GError *error = NULL;
	gchar* project_path = g_build_filename(g_get_user_data_dir(), "siril", NULL);

	// Check/create venv
	if (!check_or_create_venv(project_path, &error)) {
		siril_log_color_message(_("Failed to initialize Python virtual environment: %s\n"), "red",
				error ? error->message : "Unknown error");
		g_clear_error(&error);
		g_free(project_path);
		return GINT_TO_POINTER(1);
	}

	// Prepare venv environment
	gchar *venv_path = g_build_filename(project_path, "venv", NULL);
	PythonVenvInfo *venv_info = prepare_venv_environment(venv_path);
	if (!venv_info) {
		siril_log_color_message(_("Failed to prepare virtual environment\n"), "red");
		g_free(venv_path);
		g_free(project_path);
		return GINT_TO_POINTER(1);
	}

	// Set up environment for connection
	GHashTableIter iter;
	gpointer key, value;
	g_hash_table_iter_init(&iter, venv_info->env_vars);
	while (g_hash_table_iter_next(&iter, &key, &value)) {
		if (!g_setenv((const gchar*)key, (const gchar*)value, TRUE))
			siril_debug_print("Error in g_setenv: key = %s, value = %s\n", (const gchar*) key, (const gchar*) value);
	}

	// Clean up
	if (venv_info) {
		g_free(venv_info->venv_path);
		g_free(venv_info->python_version);
		g_hash_table_destroy(venv_info->env_vars);
		g_free(venv_info);
	}
	g_free(venv_path);
	g_free(project_path);
	if (!com.headless)
		gdk_threads_add_idle(python_venv_idle, NULL);
	else
		python_venv_idle(NULL);
	return GINT_TO_POINTER(0);
}

void initialize_python_venv_in_thread() {
	GError *error = NULL;
	com.python_init_thread = g_thread_try_new("initialize python venv", initialize_python_venv, NULL, &error);
	// We clean up the thread in python_venv_idle
}

void shutdown_python_communication(CommunicationState *commstate) {
	if (commstate->python_conn) {
		cleanup_connection(commstate->python_conn);
		commstate->python_conn = NULL;
	}

	if (commstate->worker_thread) {
		g_thread_join(commstate->worker_thread);
		commstate->worker_thread = NULL;
	}
}

void execute_python_script(gchar* script_name, gboolean from_file, gboolean sync, gchar** argv_script) {
	version_number none = { 0 };
	if (compare_version(none, com.python_version) >= 0) {
		if (com.python_init_thread) {
			g_thread_join(com.python_init_thread); // wait for python initialization to start
			com.python_init_thread = NULL;
		} else {
			siril_log_color_message(_("Error: python not ready yet. This may happen at first run "
					"if the python venv and module setup has not yet completed. Please wait a short "
					"time for a completion message in the log and try again.\n"), "red");
			return;
		}
	}

	// Generate a unique connection path for the pipe or socket for this script
	gchar *connection_path = NULL;
	gchar *uuid = g_uuid_string_random();
	#ifdef _WIN32
	connection_path = g_strdup_printf("\\\\.\\pipe\\%s", uuid);
	#else
	connection_path = g_strdup_printf("/tmp/%s.sock", uuid);
	#endif

	// Create connection
	CommunicationState commstate = {0};
	commstate.python_conn = create_connection(connection_path);
	if (!commstate.python_conn) {
		siril_log_color_message(_("Error: failed to create Python connection.\n"), "red");
		g_free(uuid);
		g_free(connection_path);
		return;
	}

	// Create worker thread
	commstate.worker_thread = g_thread_new("python-comm",
										connection_worker,
										commstate.python_conn);

	if (!commstate.python_conn) {
		siril_log_color_message(_("Error: Python connection not available.\n"), "red");
		return;
	}
	init_shm_tracking(commstate.python_conn);

	// Get base environment
	gchar** env = g_get_environ();

	// Handle virtual environment if active
	const gchar* venv_path = g_getenv("VIRTUAL_ENV");
	if (venv_path != NULL) {
		// Get the specific Python version for the site-packages path
		gchar* python_version = get_venv_python_version(venv_path);
		if (python_version) {
			gchar* site_packages = NULL;
#ifdef _WIN32
			site_packages = g_build_filename(venv_path, "Lib", "site-packages", NULL);
#else
			site_packages = g_build_filename(venv_path, "lib", g_strdup_printf("python%s", python_version), "site-packages", NULL);
#endif
			// Update PYTHONPATH to include site-packages
			const gchar* current_pythonpath = g_environ_getenv(env, "PYTHONPATH");
			if (current_pythonpath != NULL) {
				gchar* new_pythonpath = g_strjoin(G_SEARCHPATH_SEPARATOR_S, site_packages, current_pythonpath, NULL);
				env = g_environ_setenv(env, "PYTHONPATH", new_pythonpath, TRUE);
				g_free(new_pythonpath);
			} else {
				env = g_environ_setenv(env, "PYTHONPATH", site_packages, TRUE);
			}
			g_free(site_packages);
			g_free(python_version);
		}
	}

	// Set up connection information in environment
#ifdef _WIN32
	env = g_environ_setenv(env, "MY_PIPE", connection_path, TRUE);
#else
	env = g_environ_setenv(env, "MY_SOCKET", commstate.python_conn->socket_path, TRUE);
#endif
	// Set PYTHONUNBUFFERED in environment
	env = g_environ_setenv(env, "PYTHONUNBUFFERED", "1", TRUE);
	gchar *python_path = find_venv_python_exe(venv_path, TRUE);
	gboolean success = FALSE;
	gchar *working_dir = NULL;
	GError* error = NULL;
	GPid child_pid;
	gint stdout_fd, stderr_fd;
	if (!python_path) {
		siril_log_color_message(_("Error finding venv python path, unable to spawn python.\n"), "red");
	} else {
		// Clear any ROI that is set
		on_clear_roi();

		// Basic argv to spawn python to run the script
		GPtrArray* python_argv = g_ptr_array_new();
		g_ptr_array_add(python_argv, python_path);
		g_ptr_array_add(python_argv, "-u");  // Set unbuffered mode

		if (from_file) {
			g_ptr_array_add(python_argv, script_name);
		} else {
			g_ptr_array_add(python_argv, "-c");
			g_ptr_array_add(python_argv, script_name);
		}

		// Add delimiter and script arguments if script arguments are provided
		if (argv_script != NULL) {
			// Add all script arguments
			for (int i = 0; argv_script[i] != NULL; i++) {
				g_ptr_array_add(python_argv, argv_script[i]);
			}
		}

		// Null-terminate the array
		g_ptr_array_add(python_argv, NULL);

		// Use the GPtrArray for spawning
		GSpawnFlags spawn_flags = G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD;

		success = g_spawn_async_with_pipes(
			working_dir,
			(gchar**)python_argv->pdata,
			env,
			spawn_flags,
			NULL,
			NULL,
			&child_pid,
			NULL,
			&stdout_fd,
			&stderr_fd,
			&error
		);
		// Free the GPtrArray (but not its contents - the caller must free argv_script)
		g_ptr_array_free(python_argv, FALSE);
		g_free(python_path);

		if (success) {
			// Set the flag that a python script is running
			com.python_script = TRUE;
			siril_debug_print("***** com.python_script flag set\n");
			// Prepend this process to the list of child processes com.children
			child_info *child = g_malloc(sizeof(child_info));
			child->childpid = child_pid;
			child->program = EXT_PYTHON;
			gchar *script_basename = g_path_get_basename(script_name);
			child->name = g_strdup_printf("%s %s", PYTHON_EXE, from_file ? script_basename : "script");
			g_free(script_basename);
			child->datetime = g_date_time_new_now_local();
			com.children = g_slist_prepend(com.children, child);
		}
	}

	if (success) {
		if (sync) {
			// Cross-platform process waiting
#ifdef _WIN32
			// Use Windows-specific waiting
			HANDLE process_handle = (HANDLE) child_pid;
			if (process_handle != NULL) {
				WaitForSingleObject(process_handle, INFINITE);
				CloseHandle(process_handle);
			}
#else
			// Use POSIX waitpid
			gint status;
			waitpid(child_pid, &status, 0);
#endif
			// If we had the python thread lock and failed to release it, release it now
			if (commstate.python_conn->thread_claimed) {
				com.python_claims_thread = FALSE;
			}

			// Clean up shared memory resources
			cleanup_shm_resources(commstate.python_conn);

			// Clean up the Connection
			free(commstate.python_conn);

			// Close the process handle
			g_spawn_close_pid(child_pid);
			remove_child_from_children(child_pid);

			// Check if it's OK to clear com.python_script
			check_python_flag();

			// Re-enable widgets
			gui_function(script_widgets_idle, NULL);
		} else {
			// Set up child process monitoring: the callback will clean up any overlooked shm resources
			g_child_watch_add(child_pid, (GChildWatchFunc)cleanup_child_process, commstate.python_conn);
		}
	} else {
		// Clean up shared memory resources
		cleanup_shm_resources(commstate.python_conn);
		// Clean up the Connection
		free(commstate.python_conn);
		// Re-enable widgets
		gui_function(script_widgets_idle, NULL);
	}

	g_strfreev(env);
	g_free(script_name);

	if (!success && error) {
		siril_log_color_message(_("Failed to execute Python script: %s\n"), "red", error->message);
		g_error_free(error);
	}
	if (!success) {
		g_free(working_dir);
		return;
	}


	// Create input streams with appropriate flags
	GInputStream *stdout_stream = NULL;
	GInputStream *stderr_stream = NULL;

#ifdef _WIN32
	stdout_stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(stdout_fd), FALSE);
	stderr_stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(stderr_fd), FALSE);
#else
	stdout_stream = g_unix_input_stream_new(stdout_fd, TRUE);
	g_unix_set_fd_nonblocking(stdout_fd, TRUE, &error);
	stderr_stream = g_unix_input_stream_new(stderr_fd, TRUE);
	g_unix_set_fd_nonblocking(stderr_fd, TRUE, &error);
#endif

	// Create unbuffered data input streams
	GDataInputStream *stdout_data = g_data_input_stream_new(stdout_stream);
	GDataInputStream *stderr_data = g_data_input_stream_new(stderr_stream);

	// Set smaller buffer size for more responsive output
	g_object_set(stdout_data, "buffer-size", 1024, NULL);
	g_object_set(stderr_data, "buffer-size", 1024, NULL);

	// Start monitoring threads
	GThread *stdout_thread = g_thread_new("stdout-monitor",
		(GThreadFunc)monitor_stream_stdout,
		g_object_ref(stdout_data));
	GThread *stderr_thread = g_thread_new("stderr-monitor",
		(GThreadFunc)monitor_stream_stderr,
		g_object_ref(stderr_data));

	// Clean up references
	g_object_unref(stdout_stream);
	g_object_unref(stderr_stream);
	g_thread_unref(stdout_thread);
	g_thread_unref(stderr_thread);

	siril_debug_print("Python script launched asynchronously with PID %d\n", child_pid);
	g_free(working_dir);
}
