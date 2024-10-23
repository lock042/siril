#ifndef SRC_IO_SIRIL_PYTHONMODULE_H
#define SRC_IO_SIRIL_PYTHONMODULE_H

// Status codes for command responses
typedef enum {
	STATUS_OK = 0,
	STATUS_ERROR = 0xFF
} StatusCode;

// Command definitions matching Python side
typedef enum {
	CMD_GET_DIMENSIONS = 1,
	CMD_GET_PIXELDATA = 2,
	CMD_GET_PIXELDATA_REGION = 3,
	CMD_NOTIFY_PROGRESS = 4,
	CMD_SEND_COMMAND = 5,
	CMD_LOG_MESSAGE = 6,
	CMD_GET_WORKING_DIRECTORY = 7,
	CMD_GET_FILENAME = 8,
	CMD_SET_PIXELDATA = 9,
	CMD_ERROR = 0xFF
} CommandType;

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

// Structure to track shared memory allocations
typedef struct {
    char *shm_name;
    void *shm_ptr;
    size_t size;
#ifdef _WIN32
    win_shm_handle_t win_handle;
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

// Public functions
gpointer open_python_channel(gpointer user_data);
int release_python_channel();
void execute_python_script_async(gchar* script_name, gboolean from_file);
gboolean send_response(GIOChannel* channel, uint8_t status, const void* data, uint32_t length);
gboolean handle_pixeldata_request(Connection *conn, rectangle region);
gboolean handle_set_pixeldata_request(Connection *conn, const char* payload, size_t payload_length);
gboolean cleanup_shm_by_name(const char *shm_name);
#endif
