#include "core/siril.h"
#include "algos/statistics.h"
#include "core/command_line_processor.h"
#include "core/siril_log.h"
#include "gui/progress_and_log.h"
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
		g_warning("Received incomplete command payload");
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
			rectangle region = {0, 0, gfit.rx, gfit.ry};
			success = handle_pixeldata_request(conn, region);
			break;
		}

		case CMD_GET_PIXELDATA_REGION: {
			if (payload_length == 16) {
				rectangle region_BE = *(rectangle*) payload;
				rectangle region = {GUINT32_FROM_BE(region_BE.x),
									GUINT32_FROM_BE(region_BE.y),
									GUINT32_FROM_BE(region_BE.w),
									GUINT32_FROM_BE(region_BE.h)};
				success = handle_pixeldata_request(conn, region);
			} else {
				g_warning("Unexpected payload length %u received for GET_PIXELDATA_REGION", payload_length);
			}
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

		case CMD_SET_PIXELDATA: {
			if (payload_length != sizeof(incoming_image_info_t)) {
				g_warning("Invalid payload length for SET_PIXELDATA: %u", payload_length);
				const char* error_msg = "Invalid payload length";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
			} else {
				success = handle_set_pixeldata_request(conn, payload, payload_length);
			}
			break;
		}

		case CMD_GET_IMAGE_STATS: {
			if (payload_length != sizeof(uint32_t)) {
				g_warning("Invalid payload length for GET_IMAGE_STATS: %u", payload_length);
				const char* error_msg = "Invalid payload length";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Get the channel number from payload (convert from network byte order)
			uint32_t channel_BE = *(uint32_t*)payload;
			uint32_t channel = GUINT32_FROM_BE(channel_BE);

			// Check if image is loaded and channel is valid
			if (!single_image_is_loaded() || channel >= gfit.naxes[2]) {
				const char* error_msg = "No image loaded or invalid channel";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			imstats **statsarray = NULL;
			// (Re)compute image stats if they are not all computed already
			if (!gfit.stats[channel] || gfit.stats[channel]->location == NULL_STATS) {
				statsarray = calloc(3, sizeof(imstats*));
				if (compute_all_channels_statistics_single_image(&gfit, STATS_EXTRA,
		MULTI_THREADED, statsarray)) {
					const char* error_msg = "Unable to compute image statistics";
					success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				} else {
					// Replace the stats in gfit with the new calculation
					full_stats_invalidation_from_fit(&gfit);
					gfit.stats = statsarray;
				}
			}

			imstats *stats = gfit.stats[channel];
			// Prepare response buffer with correct byte order
			// Size is 2 longs (16 bytes) + 12 doubles (96 bytes) = 112 bytes
			unsigned char response_buffer[112];
			unsigned char *ptr = response_buffer;

			// Convert and copy integer values (64-bit)
			int64_t total_BE = GINT64_TO_BE(stats->total);
			int64_t ngoodpix_BE = GINT64_TO_BE(stats->ngoodpix);
			memcpy(ptr, &total_BE, sizeof(int64_t));
			ptr += sizeof(int64_t);
			memcpy(ptr, &ngoodpix_BE, sizeof(int64_t));
			ptr += sizeof(int64_t);

			// Convert and copy double values ensuring big-endian byte order
			double values[] = {
				stats->mean, stats->median, stats->sigma, stats->avgDev,
				stats->mad, stats->sqrtbwmv, stats->location, stats->scale,
				stats->min, stats->max, stats->normValue, stats->bgnoise
			};

			for (int i = 0; i < 12; i++) {
				// Convert double to big-endian
				union {
					double d;
					unsigned char bytes[8];
				} convert;
				convert.d = values[i];

				// Swap bytes if on little-endian system
				#if G_BYTE_ORDER == G_LITTLE_ENDIAN
					for (int j = 0; j < 4; j++) {
						unsigned char temp = convert.bytes[j];
						convert.bytes[j] = convert.bytes[7-j];
						convert.bytes[7-j] = temp;
					}
				#endif

				memcpy(ptr, convert.bytes, sizeof(double));
				ptr += sizeof(double);
			}

			// Send the statistics data
			success = send_response(conn->channel, STATUS_OK, response_buffer, sizeof(response_buffer));
			break;
		}

		case CMD_UPDATE_PROGRESS: {
			if (payload_length < sizeof(float)) {
				g_warning("Invalid payload length for UPDATE_PROGRESS: %u", payload_length);
				const char* error_msg = "Invalid payload length";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Extract the progress value (first 4 bytes)
			float progress_BE;
			memcpy(&progress_BE, payload, sizeof(float));

			// Convert from network byte order
			union {
				float f;
				uint32_t i;
			} progress_convert;
			progress_convert.i = GUINT32_FROM_BE(*(uint32_t*)&progress_BE);
			float progress = progress_convert.f;

			// Get the message string (remaining bytes)
			const char* message = payload + sizeof(float);
			size_t message_length = payload_length - sizeof(float);

			// Create null-terminated copy of the message
			char* progress_msg = g_strndup(message, message_length);

			// Update the progress
			set_progress_bar_data(progress_msg, progress);

			// Clean up
			g_free(progress_msg);

			// Send success response
			success = send_response(conn->channel, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_GET_KEYWORDS: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t strings_size = FLEN_VALUE * 13;  // 13 string fields of FLEN_VALUE
			size_t numeric_size = sizeof(uint64_t) * 39; // 39 vars packed to 64-bit

			size_t total_size = strings_size + numeric_size;
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			// Copy string fields with fixed length
			#define COPY_STRING(str) \
				{ \
					memset(ptr, 0, FLEN_VALUE); \
					strncpy((char*)ptr, str, FLEN_VALUE - 1); \
					ptr += FLEN_VALUE; \
				}

			COPY_STRING(gfit.keywords.program);
			COPY_STRING(gfit.keywords.filename);
			COPY_STRING(gfit.keywords.row_order);
			COPY_STRING(gfit.keywords.filter);
			COPY_STRING(gfit.keywords.image_type);
			COPY_STRING(gfit.keywords.object);
			COPY_STRING(gfit.keywords.instrume);
			COPY_STRING(gfit.keywords.telescop);
			COPY_STRING(gfit.keywords.observer);
			COPY_STRING(gfit.keywords.sitelat_str);
			COPY_STRING(gfit.keywords.sitelong_str);
			COPY_STRING(gfit.keywords.bayer_pattern);
			COPY_STRING(gfit.keywords.focname);

			// Helper macro for network byte order conversion of different types
			#define COPY_BE(val, type) \
				{ \
					union { type v; uint64_t i; } conv; \
					conv.v = val; \
					uint64_t be = GUINT64_TO_BE(conv.i); \
					memcpy(ptr, &be, sizeof(type)); \
					ptr += sizeof(type); \
				}

			// Copy numeric values with proper byte order conversion. All
			// types shorter than 64bit are converted to 64bit types before
			// endianness conversion and transmission, to simplify the data
			COPY_BE(gfit.keywords.bscale, double);
			COPY_BE(gfit.keywords.bzero, double);
			COPY_BE((uint64_t) gfit.keywords.lo, uint64_t);
			COPY_BE((uint64_t) gfit.keywords.hi, uint64_t);
			COPY_BE((double) gfit.keywords.flo, double);
			COPY_BE((double) gfit.keywords.fhi, double);
			COPY_BE(gfit.keywords.data_max, double);
			COPY_BE(gfit.keywords.data_min, double);
			COPY_BE(gfit.keywords.pixel_size_x, double);
			COPY_BE(gfit.keywords.pixel_size_y, double);
			COPY_BE((uint64_t) gfit.keywords.binning_x, uint64_t);
			COPY_BE((uint64_t) gfit.keywords.binning_y, uint64_t);
			COPY_BE(gfit.keywords.expstart, double);
			COPY_BE(gfit.keywords.expend, double);
			COPY_BE(gfit.keywords.centalt, double);
			COPY_BE(gfit.keywords.centaz, double);
			COPY_BE(gfit.keywords.sitelat, double);
			COPY_BE(gfit.keywords.sitelong, double);
			COPY_BE(gfit.keywords.siteelev, double);
			COPY_BE((int64_t) gfit.keywords.bayer_xoffset, int64_t);
			COPY_BE((int64_t) gfit.keywords.bayer_yoffset, int64_t);
			COPY_BE(gfit.keywords.airmass, double);
			COPY_BE(gfit.keywords.focal_length, double);
			COPY_BE(gfit.keywords.flength, double);
			COPY_BE(gfit.keywords.iso_speed, double);
			COPY_BE(gfit.keywords.exposure, double);
			COPY_BE(gfit.keywords.aperture, double);
			COPY_BE(gfit.keywords.ccd_temp, double);
			COPY_BE(gfit.keywords.set_temp, double);
			COPY_BE(gfit.keywords.livetime, double);
			COPY_BE((uint64_t) gfit.keywords.stackcnt, uint64_t);
			COPY_BE(gfit.keywords.cvf, double);
			COPY_BE((int64_t) gfit.keywords.key_gain, int64_t);
			COPY_BE((int64_t) gfit.keywords.key_offset, int64_t);
			COPY_BE((int64_t) gfit.keywords.focuspos, int64_t);
			COPY_BE((int64_t) gfit.keywords.focussz, int64_t);
			COPY_BE(gfit.keywords.foctemp, double);

			// Convert GDateTime to Unix timestamp
			int64_t date_ts = gfit.keywords.date ? g_date_time_to_unix(gfit.keywords.date) : 0;
			int64_t date_obs_ts = gfit.keywords.date_obs ? g_date_time_to_unix(gfit.keywords.date_obs) : 0;
			COPY_BE(date_ts, int64_t);
			COPY_BE(date_obs_ts, int64_t);

			#undef COPY_STRING
			#undef COPY_BE

			success = send_response(conn->channel, STATUS_OK, response_buffer, total_size);
			g_free(response_buffer);
			break;
		}

		case CMD_GET_IMAGE: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn->channel, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t total_size = sizeof(uint64_t) * 11; // 11 vars packed to 64-bit

			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			// Helper macro for network byte order conversion of different types
			#define COPY_BE(val, type) \
				{ \
					union { type v; uint64_t i; } conv; \
					conv.v = val; \
					uint64_t be = GUINT64_TO_BE(conv.i); \
					memcpy(ptr, &be, sizeof(type)); \
					ptr += sizeof(type); \
				}

			// Copy numeric values with proper byte order conversion. All
			// types shorter than 64bit are converted to 64bit types before
			// endianness conversion and transmission, to simplify the data
			COPY_BE((int64_t) gfit.bitpix, int64_t);
			COPY_BE((int64_t) gfit.orig_bitpix, int64_t);
			COPY_BE((uint64_t) gfit.checksum, uint64_t);
			COPY_BE(gfit.mini, double);
			COPY_BE(gfit.maxi, double);
			COPY_BE((double) gfit.neg_ratio, double);
			COPY_BE((uint64_t) gfit.type, uint64_t);
			COPY_BE((uint64_t) gfit.top_down, uint64_t);
			COPY_BE((uint64_t) gfit.focalkey, uint64_t);
			COPY_BE((uint64_t) gfit.pixelkey, uint64_t);
			COPY_BE((uint64_t) gfit.color_managed, uint64_t);

			#undef COPY_BE

			success = send_response(conn->channel, STATUS_OK, response_buffer, total_size);
			g_free(response_buffer);
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
