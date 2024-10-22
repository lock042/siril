#include "core/siril.h"
#include "core/command_line_processor.h"
#include "core/siril_log.h"
#include "io/single_image.h"
#include "io/siril_pythonmodule.h"

/**
* Process the received connection data
*/
void process_connection(Connection* conn, const gchar* buffer, gsize length) {
	if (length < sizeof(CommandHeader)) {
		g_warning("Received incomplete command header");
		return;
	}

	// Get command header
	CommandHeader* header = (CommandHeader*)buffer;
	uint32_t payload_length = GUINT32_FROM_BE(header->length);  // Convert from network byte order

	// Verify we have complete message
	if (length < sizeof(CommandHeader) + payload_length) {
		g_warning("Received incomplete command payload, length = %lu; payload_length = %lu; CommandHeader = %lu, expected length = %lu", length, payload_length, sizeof(CommandHeader), sizeof(CommandHeader) + payload_length);
		return;
	}

	// Get payload
	const char* payload = buffer + sizeof(CommandHeader);
	gboolean success = FALSE;

	// Process commands
	switch (header->command) {
		case CMD_GET_DIMENSIONS: {
			// Assuming you have a function to get the current image dimensions
			gboolean result = single_image_is_loaded(); // Implement this function

			if (result) {
				// Prepare the response data: width (gfit.rx), height (gfit.ry), and channels (gfit.naxes[2])
				uint8_t response_data[12]; // 3 x 4 bytes for width, height, and channels

				// Convert the integers to BE format for consistency across the UNIX socket
				uint32_t width_BE = GUINT32_TO_BE(gfit.rx);
				uint32_t height_BE = GUINT32_TO_BE(gfit.ry);
				uint32_t channels_BE = GUINT32_TO_BE(gfit.naxes[2]);

				// Copy the packed data into the response buffer
				memcpy(response_data, &width_BE, sizeof(uint32_t));      // First 4 bytes: width
				memcpy(response_data + 4, &height_BE, sizeof(uint32_t)); // Next 4 bytes: height
				memcpy(response_data + 8, &channels_BE, sizeof(uint32_t)); // Final 4 bytes: channels

				// Send success response with dimensions
				success = send_response(conn->channel, STATUS_OK, response_data, sizeof(response_data));
			} else {
				// Handle error retrieving dimensions
				const char* error_msg = "Failed to retrieve image dimensions - no image loaded";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_PIXELDATA: {
				success = handle_pixeldata_request(conn);
				break;
		}

		case CMD_GET_PIXELDATA_REGION: {
			break;
		}

		case CMD_NOTIFY_PROGRESS: {
			// Clean up the specified shared memory
			if (payload_length == 256) {
				success = cleanup_shm_by_name(payload);
				if (!success) {
					g_warning("Failed to find or clean up shared memory: %s", payload);
				}
			} else {
				g_warning("Invalid shared memory name length: %u", payload_length);
			}
			break;
		}

		case CMD_SEND_COMMAND: {
			// Ensure null-terminated string for command
			char* cmd = g_strndup(payload, payload_length);
			int retval = processcommand(cmd, TRUE);
			g_free(cmd);

			// Send response based on command execution
			uint8_t status = (retval == CMD_OK) ? STATUS_OK : STATUS_ERROR;
			success = send_response(conn->channel, status, NULL, 0);
			break;
		}

		case CMD_LOG_MESSAGE: {
			// Ensure null-terminated string for log message
			char* log_msg = g_strndup(payload, payload_length);
			siril_log_message(log_msg);
			g_free(log_msg);

			// Send success response
			success = send_response(conn->channel, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_GET_WORKING_DIRECTORY: {
			// Ensure the working directory is available
			if (com.wd && strlen(com.wd) > 0) {
				// Send success response with the working directory string
				success = send_response(conn->channel, STATUS_OK, com.wd, strlen(com.wd));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = "Working directory not set";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_FILENAME: {
			// Ensure the image filename is available
			if (single_image_is_loaded() && com.uniq && strlen(com.uniq->filename) > 0) {
				// Send success response with the working directory string
				success = send_response(conn->channel, STATUS_OK, com.uniq->filename, strlen(com.uniq->filename));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = "Image not loaded";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		default:
			g_warning("Unknown command: %d", header->command);
			const char* error_msg = "Unknown command";
			send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
			break;
	}

	if (!success) {
		g_warning("Failed to send response");
	}
}
