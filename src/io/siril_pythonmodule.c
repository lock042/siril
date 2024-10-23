#include <glib.h>
#include <glib-unix.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
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

#define BUFFER_SIZE 65536
#define PIPE_NAME "\\\\.\\pipe\\mypipe"
#define SOCKET_PORT 12345

// Structure to hold process information
typedef struct {
	GPid pid;
	guint watch_source_id;
} ProcessInfo;

// Global list of active shared memory allocations
static GMutex shm_list_mutex;
static GList *active_shm_list = NULL;

// Global process info
static ProcessInfo *python_process = NULL;

// Function declarations
static void setup_io_channel(GIOChannel *channel, const char *name);
static gboolean io_watch_callback(GIOChannel *channel, GIOCondition condition, gpointer data);

// Add new shared memory allocation to tracking list
static void track_shm_allocation(const char *shm_name, void *shm_ptr, size_t size,
#ifdef _WIN32
							win_shm_handle_t *win_handle
#else
							int fd
#endif
) {
	shm_allocation_t *alloc = g_new(shm_allocation_t, 1);
	alloc->shm_name = g_strdup(shm_name);
	alloc->shm_ptr = shm_ptr;
	alloc->size = size;
#ifdef _WIN32
	alloc->win_handle = *win_handle;
#else
	alloc->fd = fd;
#endif

	g_mutex_lock(&shm_list_mutex);
	active_shm_list = g_list_append(active_shm_list, alloc);
	g_mutex_unlock(&shm_list_mutex);
}

// Free a specific shared memory allocation
static void free_shm_allocation(shm_allocation_t *alloc) {
	if (!alloc) return;

#ifdef _WIN32
	cleanup_shared_memory_win32(&alloc->win_handle);
#else
	munmap(alloc->shm_ptr, alloc->size);
	close(alloc->fd);
	shm_unlink(alloc->shm_name);
#endif

	g_free(alloc->shm_name);
	g_free(alloc);
}

// Find and remove a shared memory allocation by name
gboolean cleanup_shm_by_name(const char *shm_name) {
	gboolean found = FALSE;

	g_mutex_lock(&shm_list_mutex);
	GList *link = active_shm_list;
	while (link) {
		shm_allocation_t *alloc = link->data;
		if (strcmp(alloc->shm_name, shm_name) == 0) {
			free_shm_allocation(alloc);
			active_shm_list = g_list_delete_link(active_shm_list, link);
			found = TRUE;
			break;
		}
		link = link->next;
	}
	g_mutex_unlock(&shm_list_mutex);

	return found;
}

// Cleanup handler
static void cleanup_all_shm(GPid pid, gint status, gpointer user_data) {
	g_mutex_lock(&shm_list_mutex);

	if (active_shm_list) {
		g_warning("Cleaning up %d shared memory allocations at Python exit, these should have been cleaned before", g_list_length(active_shm_list));
		g_list_foreach(active_shm_list, (GFunc)free_shm_allocation, NULL);
		g_list_free(active_shm_list);
		active_shm_list = NULL;
	}

	g_mutex_unlock(&shm_list_mutex);

	// Remove the watch source
	if (python_process && python_process->watch_source_id > 0) {
		g_source_remove(python_process->watch_source_id);
		python_process->watch_source_id = 0;
	}
}

// Add cleanup function for program exit
static void cleanup_python_process(void) {
	if (python_process) {
		if (python_process->pid > 0) {
			// Try to terminate gracefully first
			kill(python_process->pid, SIGTERM);

			// Give it a short time to clean up
			struct timespec ts = { .tv_sec = 0, .tv_nsec = 500000000 }; // 500ms
			nanosleep(&ts, NULL);

			// Force kill if still running
			if (kill(python_process->pid, 0) == 0) {
				kill(python_process->pid, SIGKILL);
			}
		}

		if (python_process->watch_source_id > 0) {
			g_source_remove(python_process->watch_source_id);
		}

		g_free(python_process);
		python_process = NULL;
	}

	// Clean up any remaining shared memory
	cleanup_all_shm(0, 0, NULL);
}

// New function to initialize process monitoring
static void init_process_monitoring(GPid pid) {
	if (python_process) {
		// Clean up existing process info if any
		if (python_process->watch_source_id > 0) {
			g_source_remove(python_process->watch_source_id);
		}
		g_free(python_process);
	}

	python_process = g_new0(ProcessInfo, 1);
	python_process->pid = pid;

	// Set up process monitoring with full flags
	python_process->watch_source_id = g_child_watch_add_full(
		G_PRIORITY_DEFAULT,
		pid,
		cleanup_all_shm,
		NULL,
		NULL  // No need for destroy notify
	);
}

static gboolean handle_new_connection(GIOChannel* source, GIOCondition condition, gpointer data) {
	Connection* conn = (Connection*)data;

	if (condition & G_IO_IN) {
		int client_fd = accept(conn->server_fd, NULL, NULL);
		if (client_fd == -1) {
			g_warning("Failed to accept connection: %s", g_strerror(errno));
			return TRUE;  // Keep the source
		}

		// Close any existing client connection
		if (conn->channel) {
			g_io_channel_unref(conn->channel);
		}
		if (conn->socket_fd > 0) {
			close(conn->socket_fd);
		}

		// Set up the new client connection
		conn->socket_fd = client_fd;
		conn->channel = g_io_channel_unix_new(client_fd);

		// Set channel encoding to NULL for binary data
		GError *error = NULL;
		g_io_channel_set_encoding(conn->channel, NULL, &error);
		if (error) {
			g_warning("Failed to set channel encoding: %s", error->message);
			g_error_free(error);
		}

		// Set non-blocking mode
		g_io_channel_set_flags(conn->channel, G_IO_FLAG_NONBLOCK, &error);
		if (error) {
			g_warning("Failed to set non-blocking mode: %s", error->message);
			g_error_free(error);
		}

		// Notify that client has connected and channel is ready
		if (conn->client_connected_callback) {
			conn->client_connected_callback(conn);
		}
	}

	return TRUE;  // Return TRUE to keep the source
}

Connection* create_connection(void (*client_connected_callback)(Connection*)) {
	Connection* conn = g_new0(Connection, 1);
#ifdef _WIN32
	// Windows: Create a named pipe
	conn->pipe_handle = CreateNamedPipe(
		PIPE_NAME,
		PIPE_ACCESS_DUPLEX,
		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
		1,
		BUFFER_SIZE,
		BUFFER_SIZE,
		0,
		NULL
	);
	if (conn->pipe_handle == INVALID_HANDLE_VALUE) {
		g_error("Failed to create pipe");
		g_free(conn);
		return NULL;
	}
	conn->is_posix = FALSE;
	conn->channel = g_io_channel_win32_new_fd((int)conn->pipe_handle);
#else
	// POSIX: Create a socket
	int server_fd = socket(AF_UNIX, SOCK_STREAM, 0);
	if (server_fd == -1) {
		g_error("Failed to create socket");
		g_free(conn);
		return NULL;
	}

	struct sockaddr_un address;
	address.sun_family = AF_UNIX;
	strcpy(address.sun_path, "/tmp/siril.sock");  // Or whatever path you're using

	// Remove existing socket file if it exists
	unlink(address.sun_path);

	if (bind(server_fd, (struct sockaddr*)&address, sizeof(address)) < 0) {
		g_error("Failed to bind socket");
		close(server_fd);
		g_free(conn);
		return NULL;
	}
	if (listen(server_fd, 1) < 0) {
		g_error("Failed to listen on socket");
		close(server_fd);
		g_free(conn);
		return NULL;
	}

	conn->is_posix = TRUE;
	conn->server_path = g_strdup(address.sun_path);
	conn->server_fd = server_fd;
	conn->client_connected_callback = client_connected_callback;

	// Use the python MainContext
	GSource *source = g_io_create_watch(g_io_channel_unix_new(server_fd), G_IO_IN);
	g_source_set_callback(source, (GSourceFunc)handle_new_connection, conn, NULL);
	g_source_attach(source, com.python_context);
	g_source_unref(source);
#endif
	return conn;
}

/**
* Cleanup function
*/
void free_connection(Connection* conn) {
	if (!conn) return;

	if (conn->channel) {
		g_io_channel_unref(conn->channel);
	}

#ifdef _WIN32
	if (conn->pipe_handle != INVALID_HANDLE_VALUE) {
		CloseHandle(conn->pipe_handle);
	}
#else
	if (conn->socket_fd > 0) {
		close(conn->socket_fd);
	}
	if (conn->server_fd > 0) {
		close(conn->server_fd);
	}
	if (conn->server_channel) {
		g_io_channel_unref(conn->server_channel);
	}
#endif

	g_free(conn);
}

typedef struct {
	GIOChannel *channel;
	gboolean ready;
} ChannelData;

// Callback for when data is available on the channel
static gboolean channel_read_ready(GIOChannel *source, GIOCondition condition, gpointer data) {
	ChannelData *ch_data = (ChannelData *)data;

	if (condition & G_IO_IN) {
		// Data is ready to read
		ch_data->ready = TRUE;
		return FALSE; // Stop the watch
	} else if (condition & (G_IO_HUP | G_IO_ERR | G_IO_NVAL)) {
		// Handle errors or hangups
		ch_data->ready = FALSE;
		return FALSE; // Stop the watch
	}

	return TRUE; // Continue watching
}

gboolean receive_data_with_timeout(GIOChannel *channel, void *buffer, size_t *buffer_size, guint timeout_ms) {
	GError *error = NULL;
	gsize bytes_read = 0;
	gsize total_bytes_read = 0;
	GIOStatus status;

	// Set up the watch for the channel
	ChannelData ch_data = { .channel = channel, .ready = FALSE };
	GMainLoop *loop = g_main_loop_new(NULL, FALSE);

	while (total_bytes_read < *buffer_size) {
		ch_data.ready = FALSE;

		// Add watch for the channel
		guint watch_id = g_io_add_watch(channel, G_IO_IN | G_IO_ERR | G_IO_HUP,
									channel_read_ready, &ch_data);

		// Add timeout
		guint timeout_id = g_timeout_add(timeout_ms, (GSourceFunc)g_main_loop_quit, loop);

		// Run the loop
		g_main_loop_run(loop);

		// Clean up sources
		g_source_remove(watch_id);
		g_source_remove(timeout_id);

		if (!ch_data.ready) {
			g_main_loop_unref(loop);
			g_warning("Timeout or error occurred while waiting for data");
			return FALSE;
		}

		// Read available data
		status = g_io_channel_read_chars(channel,
									buffer + total_bytes_read,
									*buffer_size - total_bytes_read,
									&bytes_read,
									&error);

		if (status == G_IO_STATUS_ERROR) {
			g_warning("Error reading data from channel: %s", error->message);
			g_clear_error(&error);
			g_main_loop_unref(loop);
			return FALSE;
		} else if (status == G_IO_STATUS_EOF) {
			g_warning("End of file reached while reading from channel");
			g_main_loop_unref(loop);
			return FALSE;
		}

		total_bytes_read += bytes_read;

		// Break if we've read exactly what we needed
		if (total_bytes_read == *buffer_size) {
			break;
		}
	}

	g_main_loop_unref(loop);
	return TRUE;
}

/**
* Helper function to send a response back through the channel
*/
gboolean send_response(GIOChannel* channel, uint8_t status, const void* data, uint32_t length) {
	GError* error = NULL;
	gsize bytes_written = 0;
	ResponseHeader header = {
		.status = status,
		.length = GUINT32_TO_BE(length)  // Convert to network byte order
	};

	// Send header
	GIOStatus header_status = g_io_channel_write_chars(channel, (gchar*)&header,
													sizeof(header), &bytes_written, &error);
	if (header_status != G_IO_STATUS_NORMAL) {
		if (error) {
			g_warning("Error sending response header: %s", error->message);
			g_error_free(error);
		}
		return FALSE;
	}

	// Send data if present
	if (data && length > 0) {
		GIOStatus data_status = g_io_channel_write_chars(channel, data, length,
														&bytes_written, &error);
		if (data_status != G_IO_STATUS_NORMAL) {
			if (error) {
				g_warning("Error sending response data: %s", error->message);
				g_error_free(error);
			}
			return FALSE;
		}
	}

	// Flush the channel
	g_io_channel_flush(channel, &error);
	return TRUE;
}

static void cleanup_io_loop(void);
/**
* Function to poll the IOChannel and read data in chunks.
* Returns FALSE if the connection is closed and the watch should be removed.
*/
gboolean poll_io_channel(GIOChannel* channel, GIOCondition condition, gpointer user_data) {
	Connection* conn = (Connection*)user_data;
	GError* error = NULL;
	gsize bytes_read = 0;
	gchar buffer[BUFFER_SIZE];

	// Check for hangup or error conditions first
	if (condition & (G_IO_HUP | G_IO_ERR)) {
		siril_debug_print("Python client disconnected\n");
		goto cleanup_connection;
	}

	// Handle normal read condition
	if (condition & G_IO_IN) {
		GIOStatus status = g_io_channel_read_chars(channel, buffer, BUFFER_SIZE, &bytes_read, &error);

		switch (status) {
			case G_IO_STATUS_NORMAL:
				if (bytes_read > 0) {
					process_connection(conn, buffer, bytes_read);
				}
				break;

			case G_IO_STATUS_EOF:
				siril_debug_print("Python client closed connection\n");
				// Clean up and prepare for next connection
				goto cleanup_connection;

			case G_IO_STATUS_ERROR:
				if (error) {
					g_warning("Error reading from channel: %s", error->message);
					g_error_free(error);
				}
				return FALSE;  // Remove this watch

			case G_IO_STATUS_AGAIN:
				// Nothing to read right now, try again later
				break;
		}
	}

	return TRUE;  // Keep watching if we haven't returned FALSE above

cleanup_connection:
		// Clean up the current connection
		if (conn->channel) {
			g_io_channel_unref(conn->channel);
			conn->channel = NULL;
		}
		if (conn->socket_fd > 0) {
			close(conn->socket_fd);
			conn->socket_fd = 0;
		}
		cleanup_io_loop();
		// Return FALSE to remove this watch
		return FALSE;
}

#ifdef _WIN32
// Windows-specific shared memory handling
typedef struct {
	HANDLE mapping;
	void* ptr;
} win_shm_handle_t;

static gboolean create_shared_memory_win32(const char* name, size_t size, win_shm_handle_t* handle) {
	handle->mapping = CreateFileMapping(
		INVALID_HANDLE_VALUE,    // Use paging file
		NULL,                    // Default security
		PAGE_READWRITE,          // Read/write access
		0,                       // Maximum object size (high-order DWORD)
		size,                    // Maximum object size (low-order DWORD)
		name);                   // Name of mapping object

	if (handle->mapping == NULL) {
		g_warning("Failed to create file mapping: %lu", GetLastError());
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
		g_warning("Failed to map view of file: %lu", GetLastError());
		return FALSE;
	}

	return TRUE;
}

static void cleanup_shared_memory_win32(win_shm_handle_t* handle) {
	if (handle->ptr) {
		UnmapViewOfFile(handle->ptr);
	}
	if (handle->mapping) {
		CloseHandle(handle->mapping);
	}
}
#endif

// Handle a request for pixel data. We record the allocated SHM
// but leave clearup for another command
gboolean handle_pixeldata_request(Connection *conn, rectangle region) {
	if (!single_image_is_loaded()) {
		const char* error_msg = "Failed to retrieve pixel data - no image loaded";
		return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	// Calculate total size of pixel data
	size_t total_bytes, row_bytes;
	if (gfit.type == DATA_FLOAT) {
		row_bytes = region.w * sizeof(float);
		total_bytes = row_bytes * gfit.naxes[2] * region.h;
	} else {
		row_bytes = region.w * sizeof(WORD);
		total_bytes = row_bytes * gfit.naxes[2] * region.h;
	}

	// Generate unique name for shared memory
	char shm_name[256];
#ifdef _WIN32
	snprintf(shm_name, sizeof(shm_name), "siril_shm_%lu_%lu",
			(unsigned long)GetCurrentProcessId(),
			(unsigned long)time(NULL));
#else
	snprintf(shm_name, sizeof(shm_name), "/siril_shm_%d_%lu",
			getpid(), (unsigned long)time(NULL));
#endif

	void* shm_ptr = NULL;
#ifdef _WIN32
	win_shm_handle_t win_handle = {NULL, NULL};
	if (!create_shared_memory_win32(shm_name, total_bytes, &win_handle)) {
		return FALSE;
	}
	shm_ptr = win_handle.ptr;
#else
	int fd = shm_open(shm_name, O_CREAT | O_RDWR, 0600);
	if (fd == -1) {
		g_warning("Failed to create shared memory: %s", strerror(errno));
		return FALSE;
	}
	if (ftruncate(fd, total_bytes) == -1) {
		g_warning("Failed to set shared memory size: %s", strerror(errno));
		close(fd);
		shm_unlink(shm_name);
		return FALSE;
	}
	shm_ptr = mmap(NULL, total_bytes, PROT_READ | PROT_WRITE,
				MAP_SHARED, fd, 0);
	if (shm_ptr == MAP_FAILED) {
		g_warning("Failed to map shared memory: %s", strerror(errno));
		close(fd);
		shm_unlink(shm_name);
		return FALSE;
	}
#endif

	// Copy data from gfit to shared memory
	size_t index = 0;
	int top = region.y + region.h;
	if (gfit.type == DATA_FLOAT) {
		for (int chan = 0 ; chan < gfit.naxes[2] ; chan++) {
			for (int i = region.y ; i < top ; i++) {
				memcpy((char*)shm_ptr + index, gfit.fpdata[chan] + (i * gfit.rx + region.x), region.w * sizeof(float));
				index += row_bytes;
			}
		}
	} else {
		for (int chan = 0 ; chan < gfit.naxes[2] ; chan++) {
			for (int i = region.y ; i < top ; i++) {
				memcpy((char*)shm_ptr + index, gfit.fpdata[chan] + (i * gfit.rx + region.x), region.w * sizeof(WORD));
				index += row_bytes;
			}
		}
	}

	// Track this allocation
#ifdef _WIN32
	track_shm_allocation(shm_name, shm_ptr, total_bytes, &win_handle);
#else
	track_shm_allocation(shm_name, shm_ptr, total_bytes, fd);
#endif

	// Prepare shared memory info structure
	shared_memory_info_t info = {
		.size = total_bytes,
		.data_type = (gfit.type == DATA_FLOAT) ? 1 : 0,
		.width = region.w,
		.height = region.h,
		.channels = gfit.naxes[2]
	};
	strncpy(info.shm_name, shm_name, sizeof(info.shm_name) - 1);

	// Send shared memory info to Python
	return send_response(conn->channel, STATUS_OK, (const char*)&info, sizeof(info));
}

gboolean handle_set_pixeldata_request(Connection *conn, const char* payload, size_t payload_length) {
	if (!single_image_is_loaded()) {
		const char* error_msg = "No image loaded: set_pixel_data() can only be used to update a loaded image, not to create a new one";
		return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	if (payload_length != sizeof(incoming_image_info_t)) {
		const char* error_msg = "Invalid image info size";
		return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
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
		const char* error_msg = g_strdup_printf("Invalid image dimensions or format: w = %u, h = %u, c = %u, size = %lu", info->width, info->height, info->channels, info->size);
		return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	// Open shared memory
	void* shm_ptr = NULL;
	#ifdef _WIN32
		win_shm_handle_t win_handle = {NULL, NULL};
		HANDLE mapping = OpenFileMapping(FILE_MAP_READ, FALSE, info->shm_name);
		if (mapping == NULL) {
			const char* error_msg = "Failed to open shared memory mapping";
			return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
		}
		shm_ptr = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, info->size);
		if (shm_ptr == NULL) {
			CloseHandle(mapping);
			const char* error_msg = "Failed to map shared memory view";
			return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
		}
		win_handle.mapping = mapping;
		win_handle.ptr = shm_ptr;
	#else
		int fd = shm_open(info->shm_name, O_RDONLY, 0);
		if (fd == -1) {
			const char* error_msg = "Failed to open shared memory";
			return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
		}
		shm_ptr = mmap(NULL, info->size, PROT_READ, MAP_SHARED, fd, 0);
		if (shm_ptr == MAP_FAILED) {
			close(fd);
			const char* error_msg = "Failed to map shared memory";
			return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
		}
	#endif

	// Allocate new image buffer
	gfit.type = DATA_UNSUPPORTED; // prevents a potential crash in on_drawingarea_motion_notify_event while we are messing with gfit
	free(gfit.data);
	free(gfit.fdata);
	gfit.pdata[2] = gfit.pdata[1] = gfit.pdata[0] = gfit.data = NULL;
	gfit.fpdata[2] = gfit.fpdata[1] = gfit.fpdata[0] = gfit.fdata = NULL;
	gboolean alloc_err = FALSE;
	if (info->data_type == 0) {
		gfit.pdata[2] = gfit.pdata[1] = gfit.pdata[0] = gfit.data = calloc(info->width * info->height * info->channels, sizeof(WORD));
		if (!gfit.data) {
			alloc_err = TRUE;
		} else {
			for (int i = 0 ; i < info->channels ; i++) {
				gfit.pdata[i] = gfit.data + i * info->width * info->height;
			}
		}
	} else {
		gfit.fpdata[2] = gfit.fpdata[1] = gfit.fpdata[0] = gfit.fdata = calloc(info->width * info->height * info->channels, sizeof(float));
		if (!gfit.fdata) {
			alloc_err = TRUE;
		} else {
			for (int i = 0 ; i < info->channels ; i++) {
				gfit.fpdata[i] = gfit.fdata + i * info->width * info->height;
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
		const char* error_msg = "Failed to allocate image buffer";
		return send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
	}

	// Copy data from shared memory to gfit
	size_t total_bytes = info->width * info->height * info->channels * (info->data_type == 1 ? sizeof(float) : sizeof(WORD));

	if (info->data_type == 1) {  // float data
		memcpy(gfit.fdata, (char*) shm_ptr, total_bytes);
	} else {  // WORD data
		memcpy(gfit.data, (char*) shm_ptr, total_bytes);
	}

	// Update gfit metadata
	gfit.type = info->data_type ? DATA_FLOAT : DATA_USHORT;
	gfit.rx = gfit.naxes[0] = info->width;
	gfit.ry = gfit.naxes[1] = info->height;
	gfit.naxes[2] = info->channels;

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

	return send_response(conn->channel, STATUS_OK, NULL, 0);
}

void execute_python_script_async(gchar* script_name, gboolean from_file) {
	if (!com.python_conn) {
		siril_log_color_message(_("Error: python IO channel not available.\n"), "red");
		return;
	}

	// Initialize mutex if not already done
	static gsize initialization_done = 0;
	if (g_once_init_enter(&initialization_done)) {
		g_mutex_init(&shm_list_mutex);
		// Register cleanup handler for program exit
		atexit(cleanup_python_process);
		g_once_init_leave(&initialization_done, 1);
	}

	// Prepare environment
	gchar** env = g_get_environ();
	// Retrieve the current PYTHONPATH from the environment
	const gchar* current_pythonpath = g_environ_getenv(env, "PYTHONPATH");

	// If PYTHONPATH exists, append the new module directory; otherwise, just set it
	gchar* module_dir = NULL; // replace with somewhere we can put the module
	gchar* new_pythonpath = NULL;
	if (current_pythonpath != NULL) {
		new_pythonpath = g_strconcat(current_pythonpath, G_SEARCHPATH_SEPARATOR_S, module_dir, NULL);
	} else {
		new_pythonpath = g_strdup(module_dir);  // Just use the new module dir if PYTHONPATH is not set
	}

	// Set the new PYTHONPATH in the environment
	env = g_environ_setenv(env, "PYTHONPATH", new_pythonpath, TRUE);

	// Free the new PYTHONPATH string (env now has a copy of it)
	g_free(new_pythonpath);

	// Add or update MY_SOCKET environment variable
#ifdef _WIN32
	env = g_environ_setenv(env, "MY_PIPE", PIPE_NAME, TRUE);
#define PYTHON_EXE "python3.exe"
#else
	env = g_environ_setenv(env, "MY_SOCKET", com.python_conn->server_path, TRUE);
#define PYTHON_EXE "python3"
#endif

	// Prepare arguments
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

	GError* error = NULL;
	GPid pid;
	gint stdout_fd, stderr_fd;
	gchar* working_dir = g_strdup(com.wd);

	// Spawn process with DO_NOT_REAP_CHILD flag
	gboolean success = g_spawn_async_with_pipes(
		working_dir,
		python_argv,
		env,
		G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
		NULL,
		NULL,
		&pid,
		NULL,
		&stdout_fd,
		&stderr_fd,
		&error
	);

	if (!success) {
		g_error("Failed to execute Python script asynchronously: %s", error->message);
		g_error_free(error);
		g_free(working_dir);
		return;
	}

	// Set up process monitoring
	init_process_monitoring(pid);

	// Set up stdout monitoring
	GIOChannel *stdout_channel = g_io_channel_unix_new(stdout_fd);
	g_io_channel_set_encoding(stdout_channel, NULL, NULL);
	g_io_channel_set_flags(stdout_channel, G_IO_FLAG_NONBLOCK, NULL);
	setup_io_channel(stdout_channel, "stdout");

	// Set up stderr monitoring
	GIOChannel *stderr_channel = g_io_channel_unix_new(stderr_fd);
	g_io_channel_set_encoding(stderr_channel, NULL, NULL);
	g_io_channel_set_flags(stderr_channel, G_IO_FLAG_NONBLOCK, NULL);
	setup_io_channel(stderr_channel, "stderr");

	siril_log_message(_("Python script launched asynchronously with PID %d\n"), pid);

	g_free(working_dir);
}

// Helper function to set up I/O channels
static void setup_io_channel(GIOChannel *channel, const char *name) {
	g_io_add_watch(channel, G_IO_IN | G_IO_HUP, io_watch_callback, (gpointer)name);
	g_io_channel_unref(channel);  // Remove our reference, watch keeps its own
}

// Callback for I/O channel monitoring
static gboolean io_watch_callback(GIOChannel *channel, GIOCondition condition, gpointer data) {
	const char *name = (const char *)data;
	gchar *buffer = NULL;
	gsize bytes_read;
	GIOStatus status;

	if (condition & G_IO_HUP) {
		return FALSE;  // Remove source
	}

	if (condition & G_IO_IN) {
		status = g_io_channel_read_line(channel, &buffer, &bytes_read, NULL, NULL);
		if (status == G_IO_STATUS_NORMAL) {
			if (strcmp(name, "stderr") == 0) {
				siril_log_color_message(buffer, "red");
			} else {
				siril_log_message("%s", buffer);
			}
			g_free(buffer);
		}
	}

	return TRUE;  // Keep source
}

/**
* Main entry point.
*/
// Global variables for the IO context and loop
static GMainContext* io_context = NULL;
static GMainLoop* io_loop = NULL;

// Modified callback function that sets up the channel watch
static void on_python_client_connected(Connection* conn) {
	GSource* source = NULL;

	// Ensure IO loop is initialized
	if (!io_context) {
		io_context = g_main_context_new();
		io_loop = g_main_loop_new(io_context, FALSE);

		// Start the loop in a separate thread
		g_thread_new("io-loop", (GThreadFunc)g_main_loop_run, io_loop);
	}

	// Create an IO watch source
	source = g_io_create_watch(conn->channel, G_IO_IN | G_IO_ERR | G_IO_HUP);

	// Set the callback for the source
	g_source_set_callback(source, (GSourceFunc)poll_io_channel, conn, NULL);

	// Attach the source to our custom context
	g_source_attach(source, io_context);

	// We can unref the source now as the context holds a reference
	g_source_unref(source);

	siril_debug_print("Python IO channel opened\n");
}

// Cleanup function to be called when shutting down
static void cleanup_io_loop(void) {
	if (io_loop) {
		g_main_loop_quit(io_loop);
		g_main_loop_unref(io_loop);
		io_loop = NULL;
	}

	if (io_context) {
		g_main_context_unref(io_context);
		io_context = NULL;
	}
}

// Initialize the Python context and loop
gpointer open_python_channel(gpointer user_data) {
	// Create the python context
	com.python_context = g_main_context_new();
	com.python_loop = g_main_loop_new(com.python_context, FALSE);
	// Create the connection (pipe or socket)
	com.python_conn = create_connection(on_python_client_connected);
	if (!com.python_conn) {
		return GINT_TO_POINTER(1);
	}
	g_main_loop_run(com.python_loop);

	return GINT_TO_POINTER(0);
}

int release_python_channel() {
	if (com.python_conn) {
		if (com.python_conn->channel) {
			g_io_channel_shutdown(com.python_conn->channel, TRUE, NULL);
			g_io_channel_unref(com.python_conn->channel);
		}
		free_connection(com.python_conn);
		com.python_conn = NULL;
	}
	return 0;
}
