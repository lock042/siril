#include <glib.h>
#include <gio/gio.h>
#ifdef _WIN32
#include <windows.h>
#include <gio/gwin32inputstream.h>
#else
#include <gio/gunixinputstream.h>
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
#include <time.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "io/single_image.h"
#include "io/siril_pythoncommands.h"
#include "io/siril_pythonmodule.h"

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

// Statics
static CommunicationState commstate = {0};
static GSList* g_shm_allocations = NULL;
static GMutex g_shm_mutex;
static gboolean g_shm_initialized = FALSE;

// Forward declarations
static void cleanup_connection(Connection *conn);
static gboolean wait_for_client(Connection *conn);
static gboolean handle_client_communication(Connection *conn);

#ifdef _WIN32
static gboolean create_shared_memory_win32(const char* name, size_t size, win_shm_handle_t* handle) {
	handle->mapping = CreateFileMapping(
		INVALID_HANDLE_VALUE,    // Use paging file
		NULL,                    // Default security
		PAGE_READWRITE,          // Read/write access
		0,                       // Maximum object size (high-order DWORD)
		size,                    // Maximum object size (low-order DWORD)
		name);                   // Name of mapping object

	if (handle->mapping == NULL) {
		siril_debug_print("Failed to create file mapping: %lu\n", GetLastError());
		return FALSE;
	}

	handle->ptr = MapViewOfFile(
		handle->mapping,         // Handle to mapping object
		FILE_MAP_ALL_ACCESS,     // Read/write permission
		0,                       // Offset high
		0,                       // Offset low
		size);                   // Number of bytes to map

	if (handle->ptr == NULL) {
		CloseHandle(handle->mapping);
		siril_debug_print("Failed to map view of file: %lu\n", GetLastError());
		return FALSE;
	}

	return TRUE;
}
#endif

gboolean send_response(Connection* conn, uint8_t status, const void* data, uint32_t length) {
	ResponseHeader header = {
		.status = status,
		.length = GUINT32_TO_BE(length)  // Convert to network byte order
	};

#ifdef _WIN32
	DWORD bytes_written = 0;

	// Send header
	if (!WriteFile(conn->pipe_handle, &header, sizeof(header), &bytes_written, NULL) ||
		bytes_written != sizeof(header)) {
		siril_debug_print("Failed to send response header: %lu\n", GetLastError());
		return FALSE;
	}

	// Send data if present
	if (data && length > 0) {
		if (!WriteFile(conn->pipe_handle, data, length, &bytes_written, NULL) ||
			bytes_written != length) {
			siril_debug_print("Failed to send response data: %lu\n", GetLastError());
			return FALSE;
		}
	}

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
gboolean siril_allocate_shm(void** shm_ptr_ptr,
							char* shm_name_ptr,
							size_t total_bytes,
							win_shm_handle_t *win_handle) {
	void *shm_ptr = NULL;
	snprintf(shm_name_ptr, sizeof(shm_name_ptr), "siril_shm_%lu_%lu",
		(unsigned long)GetCurrentProcessId(),
		(unsigned long)time(NULL));
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
	snprintf(shm_name_ptr, 256, "/siril_shm_%d_%lu",
			getpid(), (unsigned long)time(NULL));

	*fd = shm_open(shm_name_ptr, O_CREAT | O_RDWR, 0600);
	if (*fd == -1) {
		siril_debug_print("Failed to create shared memory: %s\n", strerror(errno));
		return FALSE;
	}
	if (ftruncate(*fd, total_bytes) == -1) {
		siril_debug_print("Failed to set shared memory size: %s\n", strerror(errno));
		close(*fd);
		shm_unlink(shm_name_ptr);
		return FALSE;
	}
	shm_ptr = mmap(NULL, total_bytes, PROT_READ | PROT_WRITE,
				MAP_SHARED, *fd, 0);
	if (shm_ptr == MAP_FAILED) {
		siril_debug_print("Failed to map shared memory: %s\n", strerror(errno));
		close(*fd);
		shm_unlink(shm_name_ptr);
		return FALSE;
	}
	*shm_ptr_ptr = shm_ptr;
	return TRUE;
}
#endif

// Initialize tracking system
static void init_shm_tracking(void) {
    if (!g_shm_initialized) {
        g_mutex_init(&g_shm_mutex);
        g_shm_initialized = TRUE;
    }
}

// Track new allocation
#ifdef _WIN32
static void track_shm_allocation(const char* shm_name, void* shm_ptr, size_t size, win_shm_handle_t* handle) {
#else
static void track_shm_allocation(const char* shm_name, void* shm_ptr, size_t size, int fd) {
#endif
	if (!g_shm_initialized) init_shm_tracking();

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
    g_mutex_lock(&g_shm_mutex);
    g_shm_allocations = g_slist_append(g_shm_allocations, allocation);
    g_mutex_unlock(&g_shm_mutex);
}

// Helper function to match allocation by name
static gint find_allocation_by_name(gconstpointer a, gconstpointer b) {
    const shm_allocation_t* allocation = a;
    const char* name = b;
    return strcmp(allocation->shm_name, name);
}

// Cleanup allocation
void cleanup_shm_allocation(const char* shm_name) {
    if (!g_shm_initialized) return;

    g_mutex_lock(&g_shm_mutex);

    GSList* link = g_slist_find_custom(g_shm_allocations, shm_name, find_allocation_by_name);
    if (link) {
        shm_allocation_t* allocation = link->data;

        // Unmap memory
        munmap(allocation->shm_ptr, allocation->size);

#ifdef _WIN32
        CloseHandle(allocation->handle);
#else
        // Close file descriptor and then unlink shared memory
        close(allocation->fd);
        shm_unlink(allocation->shm_name);
#endif

        // Remove from list and free allocation structure
        g_shm_allocations = g_slist_remove(g_shm_allocations, allocation);
        g_free(allocation);
    }

    g_mutex_unlock(&g_shm_mutex);
}

// Cleanup all allocations (useful for shutdown)
static void cleanup_all_shm_allocations(void) {
    if (!g_shm_initialized) return;

    g_mutex_lock(&g_shm_mutex);

    GSList* current = g_shm_allocations;
    while (current) {
        shm_allocation_t* allocation = current->data;

        munmap(allocation->shm_ptr, allocation->size);

#ifdef _WIN32
        CloseHandle(allocation->handle);
#else
        close(allocation->fd);
        shm_unlink(allocation->shm_name);
#endif

        g_free(allocation);
        current = current->next;
    }

    g_slist_free(g_shm_allocations);
    g_shm_allocations = NULL;

    g_mutex_unlock(&g_shm_mutex);
}

void cleanup_shm_resources(void) {
    cleanup_all_shm_allocations();
    if (g_shm_initialized) {
        g_mutex_clear(&g_shm_mutex);
        g_shm_initialized = FALSE;
    }
}

// Handle a request for pixel data. We record the allocated SHM
// but leave clearup for another command
gboolean handle_pixeldata_request(Connection *conn, fits *fit, rectangle region) {
	if (!single_image_is_loaded()) {
		const char* error_msg = _("Failed to retrieve pixel data - no image loaded");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
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
	char *shm_name_ptr = shm_name;
#ifdef _WIN32
	win_shm_handle_t win_handle;
	if (!siril_allocate_shm(shm_ptr, shm_name, total_bytes, &win_handle))
		return FALSE;
#else
	int fd;
	if (!siril_allocate_shm(&shm_ptr, shm_name_ptr, total_bytes, &fd))
		return FALSE;
#endif

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
	track_shm_allocation(shm_name, shm_ptr, total_bytes, &win_handle);
#else
	track_shm_allocation(shm_name, shm_ptr, total_bytes, fd);
#endif

	// Prepare shared memory info structure
	shared_memory_info_t info = {
		.size = total_bytes,
		.data_type = (fit->type == DATA_FLOAT) ? 1 : 0,
		.width = region.w,
		.height = region.h,
		.channels = fit->naxes[2]
	};
	memset(info.shm_name, 0, sizeof(info.shm_name));  // Clear the buffer first
	memcpy(info.shm_name, shm_name, strlen(shm_name)); // Safe copy with implicit null termination

	// Send shared memory info to Python
	return send_response(conn, STATUS_OK, (const char*)&info, sizeof(info));
}

gboolean handle_rawdata_request(Connection *conn, void* data, size_t total_bytes) {
	if (data == NULL || total_bytes == 0) {
		const char* error_msg = _("Incorrect memory region specification");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	// Generate unique name for shared memory and allocate it
	void* shm_ptr = NULL;
	char shm_name[256];
	char *shm_name_ptr = shm_name;
#ifdef _WIN32
	win_shm_handle_t win_handle;
	if (!siril_allocate_shm(shm_ptr, shm_name, total_bytes, &win_handle))
		return FALSE;
#else
	int fd;
	if (!siril_allocate_shm(&shm_ptr, shm_name_ptr, total_bytes, &fd))
		return FALSE;
#endif

	// Copy data from gfit to shared memory
	memcpy((char*)shm_ptr, (char*) data, total_bytes);

	// Track this allocation
#ifdef _WIN32
	track_shm_allocation(shm_name, shm_ptr, total_bytes, &win_handle);
#else
	track_shm_allocation(shm_name, shm_ptr, total_bytes, fd);
#endif

	// Prepare shared memory info structure
	// We only care about the size and shm_name, but it's easier to use the same
	// struct as for handle_pixeldata_request
	shared_memory_info_t info = {
		.size = total_bytes,
		.data_type = 0,
		.width = 0,
		.height = 0,
		.channels = 0
	};
	memset(info.shm_name, 0, sizeof(info.shm_name));  // Clear the buffer first
	memcpy(info.shm_name, shm_name, strlen(shm_name)); // Safe copy with implicit null termination

	// Send shared memory info to Python
	return send_response(conn, STATUS_OK, (const char*)&info, sizeof(info));
}

gboolean handle_set_pixeldata_request(Connection *conn, fits *fit, const char* payload, size_t payload_length) {
	if (!single_image_is_loaded()) {
		const char* error_msg = _("No image loaded: set_pixel_data() can only be used to update a loaded image, not to create a new one");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	if (payload_length != sizeof(incoming_image_info_t)) {
		const char* error_msg = _("Invalid image info size");
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	incoming_image_info_t* info = (incoming_image_info_t*)payload;
	info->width = GUINT32_FROM_BE(info->width);
	info->height = GUINT32_FROM_BE(info->height);
	info->channels = GUINT32_FROM_BE(info->channels);
	info->size = GUINT64_FROM_BE(info->size);
	info->data_type = GUINT32_FROM_BE(info->data_type);
	// Validate image dimensions and format
	if (info->width == 0 || info->height == 0 || info->channels == 0 ||
		info->channels > 3 || info->size == 0) {
		gchar* error_msg = g_strdup_printf(_("Invalid image dimensions or format: w = %u, h = %u, c = %u, size = %lu"), info->width, info->height, info->channels, info->size);
		int retval = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	g_free(error_msg);
	return retval;
	}

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

	// Allocate new image buffer
	fit->type = DATA_UNSUPPORTED; // prevents a potential crash in on_drawingarea_motion_notify_event while we are messing with gfit
	free(fit->data);
	free(fit->fdata);
	fit->pdata[2] = fit->pdata[1] = fit->pdata[0] = fit->data = NULL;
	fit->fpdata[2] = fit->fpdata[1] = fit->fpdata[0] = fit->fdata = NULL;
	gboolean alloc_err = FALSE;
	if (info->data_type == 0) { // WORD data
		fit->pdata[2] = fit->pdata[1] = fit->pdata[0] = fit->data = calloc(info->width * info->height * info->channels, sizeof(WORD));
		if (!fit->data) {
			alloc_err = TRUE;
		} else {
			for (int i = 0 ; i < info->channels ; i++) {
				fit->pdata[i] = fit->data + i * info->width * info->height;
			}
		}
	} else { // FLOAT data
		fit->fpdata[2] = fit->fpdata[1] = fit->fpdata[0] = fit->fdata = calloc(info->width * info->height * info->channels, sizeof(float));
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
		return send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	// Copy data from shared memory to gfit
	size_t total_bytes = info->width * info->height * info->channels * (info->data_type == 1 ? sizeof(float) : sizeof(WORD));

	if (info->data_type == 0) {  // WORD data
		memcpy(fit->data, (char*) shm_ptr, total_bytes);
	} else {  // float data
		memcpy(fit->fdata, (char*) shm_ptr, total_bytes);
	}

	// Update gfit metadata
	fit->type = info->data_type ? DATA_FLOAT : DATA_USHORT;
	fit->rx = fit->naxes[0] = info->width;
	fit->ry = fit->naxes[1] = info->height;
	fit->naxes[2] = info->channels;

	notify_gfit_modified();

	// Cleanup shared memory
	#ifdef _WIN32
		UnmapViewOfFile(shm_ptr);
		CloseHandle(win_handle.mapping);
	#else
		munmap(shm_ptr, info->size);
		close(fd);
		shm_unlink(info->shm_name);  // Remove shared memory object
	#endif

	// In all cases we have now finished with the shm and closed and unlinked it.
	// On receipt of the response, python will also close and unlink the shm in its
	// finally: block.
	return send_response(conn, STATUS_OK, NULL, 0);
}

// Monitor stdout stream
static gpointer monitor_stream_stdout(GDataInputStream *data_input) {
	gsize length = 0;
	gchar *buffer;

	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length, NULL, NULL))) {
		siril_log_message("%s\n", buffer);
		g_free(buffer);
	}

	g_object_unref(data_input);
	return NULL;
}

// Monitor stderr stream
static gpointer monitor_stream_stderr(GDataInputStream *data_input) {
	gsize length = 0;
	gchar *buffer;

	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length, NULL, NULL))) {
		siril_log_color_message("%s\n", "red", buffer);
		g_free(buffer);
	}

	g_object_unref(data_input);
	return NULL;
}

static void cleanup_child_process(GPid pid, gint status, gpointer user_data) {
	// Log the Python process exit status if needed
	if (WIFEXITED(status)) {
		if (status == 0)
			siril_debug_print("Python process (PID: %d) exited normally\n", pid);
		else
			siril_log_color_message(_("Python process (PID: %d) exited with status %d\n"), "salmon",
				pid, WEXITSTATUS(status));
	} else if (WIFSIGNALED(status)) {
		siril_log_color_message(_("Python process (PID: %d) terminated by signal %d\n"), "salmon",
				pid, WTERMSIG(status));
	}

	// Clean up shared memory resources
	cleanup_shm_resources();

	// Close the process handle
	g_spawn_close_pid(pid);
}

void execute_python_script_async(gchar* script_name, gboolean from_file) {
	if (!commstate.python_conn) {
		siril_log_color_message(_("Error: Python connection not available.\n"), "red");
		return;
	}
	init_shm_tracking();
	// Prepare environment
	gchar** env = g_get_environ();
	const gchar* current_pythonpath = g_environ_getenv(env, "PYTHONPATH");

	// Set up module directory path
	gchar* module_dir = NULL; // replace with your module directory
	gchar* new_pythonpath = NULL;

	if (current_pythonpath != NULL) {
		new_pythonpath = g_strconcat(current_pythonpath,
								G_SEARCHPATH_SEPARATOR_S,
								module_dir,
								NULL);
	} else {
		new_pythonpath = g_strdup(module_dir);
	}

	env = g_environ_setenv(env, "PYTHONPATH", new_pythonpath, TRUE);
	g_free(new_pythonpath);

	// Set up connection information in environment
#ifdef _WIN32
	const gchar* pipe_name = "\\\\.\\pipe\\siril";
	env = g_environ_setenv(env, "MY_PIPE", pipe_name, TRUE);
#define PYTHON_EXE "python3.exe"
#else
	env = g_environ_setenv(env, "MY_SOCKET", commstate.python_conn->socket_path, TRUE);
#define PYTHON_EXE "python3"
#endif

	// Prepare command arguments
	gchar* python_argv[4];
	if (from_file) {
		python_argv[0] = PYTHON_EXE;
		python_argv[1] = script_name;
		python_argv[2] = NULL;
	} else {
		python_argv[0] = PYTHON_EXE;
		python_argv[1] = "-c";
		python_argv[2] = script_name;
		python_argv[3] = NULL;
	}

	// Set up process spawn
	GError* error = NULL;
	GPid child_pid;
	gint stdout_fd, stderr_fd;
	gchar* working_dir = g_strdup(com.wd);

	gboolean success = g_spawn_async_with_pipes(
		working_dir,
		python_argv,
		env,
		G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
		NULL,
		NULL,
		&child_pid,
		NULL,
		&stdout_fd,
		&stderr_fd,
		&error
	);

	g_strfreev(env);

	if (!success) {
		siril_log_color_message(_("Failed to execute Python script: %s\n"), "red", error->message);
		g_error_free(error);
		g_free(working_dir);
		return;
	}

	// Set up child process monitoring: the callback will clean up any overlooked shm resources
	g_child_watch_add(child_pid, (GChildWatchFunc)cleanup_child_process, NULL);

	// Create input streams for stdout and stderr
	GInputStream *stdout_stream = NULL;
	GInputStream *stderr_stream = NULL;

#ifdef _WIN32
	stdout_stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(stdout_fd), FALSE);
	stderr_stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(stderr_fd), FALSE);
#else
	stdout_stream = g_unix_input_stream_new(stdout_fd, FALSE);
	stderr_stream = g_unix_input_stream_new(stderr_fd, FALSE);
#endif

	// Create data input streams
	GDataInputStream *stdout_data = g_data_input_stream_new(stdout_stream);
	GDataInputStream *stderr_data = g_data_input_stream_new(stderr_stream);

	// Start monitoring threads
	GThread *stdout_thread = g_thread_new("stdout-monitor",
		(GThreadFunc)monitor_stream_stdout,
		g_object_ref(stdout_data));
	GThread *stderr_thread = g_thread_new("stderr-monitor",
		(GThreadFunc)monitor_stream_stderr,
		g_object_ref(stderr_data));

	// Clean up thread references
	g_thread_unref(stdout_thread);
	g_thread_unref(stderr_thread);

	siril_debug_print("Python script launched asynchronously with PID %d\n", child_pid);
	g_free(working_dir);
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
		bytes_read = read(conn->client_fd, buffer, BUFFER_SIZE);

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
	siril_log_message(_("Python communication initialized...\n"));
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

gboolean initialize_python_communication(const gchar *connection_path) {
	if (commstate.python_conn) {
		siril_debug_print("Python communication already initialized\n");
		return FALSE;
	}

#ifdef _WIN32
	commstate.python_conn = create_connection(connection_path);
#else
	commstate.python_conn = create_connection(connection_path);
#endif

	if (!commstate.python_conn) {
		return FALSE;
	}

	// Create worker thread for handling connections
	commstate.worker_thread = g_thread_new("python-comm", connection_worker, commstate.python_conn);

	return TRUE;
}

void shutdown_python_communication(void) {
	if (commstate.python_conn) {
		cleanup_connection(commstate.python_conn);
		commstate.python_conn = NULL;
	}

	if (commstate.worker_thread) {
		g_thread_join(commstate.worker_thread);
		commstate.worker_thread = NULL;
	}
}