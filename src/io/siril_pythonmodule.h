#ifndef SRC_IO_SIRIL_PYTHONMODULE_H
#define SRC_IO_SIRIL_PYTHONMODULE_H

// Status codes for command responses
#include <sys/cdefs.h>
typedef enum {
	STATUS_OK = 0,
	STATUS_NONE = 1, // for "allowed to fail" commands e.g. those that may
			// legitimately return no data
	STATUS_ERROR = 0xFF
} StatusCode;

// Command definitions matching Python side
typedef enum {
	CMD_SEND_COMMAND = 1,
	CMD_LOG_MESSAGE = 2,
	CMD_UPDATE_PROGRESS = 3,
	CMD_GET_WORKING_DIRECTORY = 4,
	CMD_GET_FILENAME = 5,
	CMD_GET_DIMENSIONS = 6,
	CMD_GET_PIXELDATA = 7,
	CMD_GET_PIXELDATA_REGION = 8,
	CMD_RELEASE_SHM = 9,
	CMD_SET_PIXELDATA = 10,
	CMD_GET_IMAGE_STATS = 11,
	CMD_GET_KEYWORDS = 12,
	CMD_GET_ICC_PROFILE = 13,
	CMD_GET_FITS_HEADER = 14,
	CMD_GET_FITS_HISTORY = 15,
	CMD_GET_FITS_UNKNOWN_KEYS = 16,
	CMD_GET_IMAGE = 17,
	CMD_GET_PSFSTARS = 18,
	CMD_GET_SEQ_STATS = 19,
	CMD_GET_SEQ_REGDATA = 20,
	CMD_GET_SEQ_IMGDATA = 21,
	CMD_GET_SEQ_PIXELDATA = 22,
	CMD_GET_SEQ_IMAGE = 23,
	CMD_GET_SEQ = 24,
	CMD_GET_CONFIG = 25,
	CMD_GET_USERCONFIG_DIR = 26,
	CMD_GET_IS_IMAGE_LOADED = 27,
	CMD_GET_IS_SEQUENCE_LOADED = 28,
	CMD_GET_SELECTION = 29,
	CMD_SET_SELECTION = 30,
	CMD_ERROR = 0xFF
} CommandType;

// Config types matching python side
typedef enum {
	CONFIG_TYPE_BOOL = 0,
	CONFIG_TYPE_INT = 1,
	CONFIG_TYPE_DOUBLE = 2,
	CONFIG_TYPE_STR = 3,
	CONFIG_TYPE_STRDIR = 4,
	CONFIG_TYPE_STRLIST = 5
} config_type_t;

// Command header structure
typedef struct {
	uint8_t command;
	uint32_t length;
} __attribute__((packed)) CommandHeader;

// Response header structure
typedef struct {
	uint8_t status;
	uint32_t length;
} __attribute__((packed)) ResponseHeader;

// Structure to pass shared memory information to Python
typedef struct {
	size_t size;
	int data_type;  // 0 for WORD, 1 for float
	int width;
	int height;
	int channels;
	char shm_name[256];
} shared_memory_info_t;

#ifdef _WIN32
// Windows-specific shared memory handling
typedef struct win_shm_handle{
	void* mapping;
	void* ptr;
} win_shm_handle_t;
#endif

// Structure to track shared memory allocations
typedef struct {
	char shm_name[256];
	void* shm_ptr;
	size_t size;
#ifdef _WIN32
	win_shm_handle_t handle;
#else
	int fd;
#endif
} shm_allocation_t;

typedef struct {
	uint32_t width;
	uint32_t height;
	uint32_t channels;
	uint32_t data_type;  // 0 for WORD (uint16), 1 for float
	uint64_t size;       // Total size in bytes
	char shm_name[256];
} incoming_image_info_t;

typedef struct _Connection {
	gboolean is_connected;
	gboolean should_stop;
	GMutex mutex;
	GCond condition;
#ifdef _WIN32
	void *pipe_handle;
#else
	int server_fd;
	int client_fd;
	gchar *socket_path;
#endif
	void (*client_connected_callback)(gpointer);
	void (*client_disconnected_callback)(gpointer);
	gpointer user_data;
	GSList* g_shm_allocations;
	GMutex g_shm_mutex;
	gboolean g_shm_initialized;

} Connection;

typedef struct {
	Connection *python_conn;
	GThread *worker_thread;

} CommunicationState;

// Public functions
//gpointer open_python_channel(gpointer user_data);
//int release_python_channel();
void execute_python_script_async(gchar* script_name, gboolean from_file);
gboolean send_response(Connection *conn, uint8_t status, const void* data, uint32_t length);
gboolean handle_pixeldata_request(Connection *conn, fits *fit, rectangle region);
gboolean handle_set_pixeldata_request(Connection *conn, fits *fit, const char* payload, size_t payload_length);
void cleanup_shm_allocation(Connection *conn, const char* shm_name);
gboolean handle_rawdata_request(Connection *conn, void* data, size_t total_bytes);
void initialize_python_venv_in_thread();
void shutdown_python_communication(CommunicationState *commstate);

#endif
