#include "core/siril.h"
#include "algos/PSF.h"
#include "algos/statistics.h"
#include "core/command_line_processor.h"
#include "core/icc_profile.h"
#include "core/siril_log.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/siril_pythonmodule.h"
#include "siril_pythonmodule.h"

// Helper macros
#define COPY_FLEN_STRING(str) \
	{ \
		size_t len = FLEN_VALUE; \
		if ((ptr + len) - start_ptr > maxlen) { \
			fprintf(stderr, "Error: Exceeded max length for COPY_FLEN_STRING at %s\n", #str); \
			return 1; \
		} \
		memset(ptr, 0, len); \
		strncpy((char*)ptr, str, len - 1); \
		ptr += len; \
	}

#define COPY_STRING(str) \
	{ \
		size_t len = strlen(str) + 1; \
		if ((ptr + len) - start_ptr > maxlen) { \
			fprintf(stderr, "Error: Exceeded max length for COPY_STRING at %s\n", #str); \
			return 1; \
		} \
		memset(ptr, 0, len); \
		strncpy((char*)ptr, str, len - 1); \
		ptr += len; \
	}

#define COPY_BE(val, type) \
	{ \
		size_t len = sizeof(type); \
		if ((ptr + len) - start_ptr > maxlen) { \
			fprintf(stderr, "Error: Exceeded max length for COPY_BE at %s\n", #val); \
			return 1; \
		} \
		union { type v; uint64_t i; } conv; \
		conv.v = val; \
		uint64_t be = GUINT64_TO_BE(conv.i); \
		memcpy(ptr, &be, len); \
		ptr += len; \
	}


static int keywords_to_py(fits *fit, unsigned char *ptr, size_t maxlen) {
	if (!fit || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	// Convert GDateTime to Unix timestamp
	int64_t date_ts = fit->keywords.date ? g_date_time_to_unix(fit->keywords.date) : 0;
	int64_t date_obs_ts = fit->keywords.date_obs ? g_date_time_to_unix(fit->keywords.date_obs) : 0;
	// Copy string fields
	COPY_FLEN_STRING(fit->keywords.program);
	COPY_FLEN_STRING(fit->keywords.filename);
	COPY_FLEN_STRING(fit->keywords.row_order);
	COPY_FLEN_STRING(fit->keywords.filter);
	COPY_FLEN_STRING(fit->keywords.image_type);
	COPY_FLEN_STRING(fit->keywords.object);
	COPY_FLEN_STRING(fit->keywords.instrume);
	COPY_FLEN_STRING(fit->keywords.telescop);
	COPY_FLEN_STRING(fit->keywords.observer);
	COPY_FLEN_STRING(fit->keywords.sitelat_str);
	COPY_FLEN_STRING(fit->keywords.sitelong_str);
	COPY_FLEN_STRING(fit->keywords.bayer_pattern);
	COPY_FLEN_STRING(fit->keywords.focname);

	// Copy numeric values with proper byte order conversion. All
	// types shorter than 64bit are converted to 64bit types before
	// endianness conversion and transmission, to simplify the data
	COPY_BE(fit->keywords.bscale, double);
	COPY_BE(fit->keywords.bzero, double);
	COPY_BE((uint64_t) fit->keywords.lo, uint64_t);
	COPY_BE((uint64_t) fit->keywords.hi, uint64_t);
	COPY_BE((double) fit->keywords.flo, double);
	COPY_BE((double) fit->keywords.fhi, double);
	COPY_BE(fit->keywords.data_max, double);
	COPY_BE(fit->keywords.data_min, double);
	COPY_BE(fit->keywords.pixel_size_x, double);
	COPY_BE(fit->keywords.pixel_size_y, double);
	COPY_BE((uint64_t) fit->keywords.binning_x, uint64_t);
	COPY_BE((uint64_t) fit->keywords.binning_y, uint64_t);
	COPY_BE(fit->keywords.expstart, double);
	COPY_BE(fit->keywords.expend, double);
	COPY_BE(fit->keywords.centalt, double);
	COPY_BE(fit->keywords.centaz, double);
	COPY_BE(fit->keywords.sitelat, double);
	COPY_BE(fit->keywords.sitelong, double);
	COPY_BE(fit->keywords.siteelev, double);
	COPY_BE((int64_t) fit->keywords.bayer_xoffset, int64_t);
	COPY_BE((int64_t) fit->keywords.bayer_yoffset, int64_t);
	COPY_BE(fit->keywords.airmass, double);
	COPY_BE(fit->keywords.focal_length, double);
	COPY_BE(fit->keywords.flength, double);
	COPY_BE(fit->keywords.iso_speed, double);
	COPY_BE(fit->keywords.exposure, double);
	COPY_BE(fit->keywords.aperture, double);
	COPY_BE(fit->keywords.ccd_temp, double);
	COPY_BE(fit->keywords.set_temp, double);
	COPY_BE(fit->keywords.livetime, double);
	COPY_BE((uint64_t) fit->keywords.stackcnt, uint64_t);
	COPY_BE(fit->keywords.cvf, double);
	COPY_BE((int64_t) fit->keywords.key_gain, int64_t);
	COPY_BE((int64_t) fit->keywords.key_offset, int64_t);
	COPY_BE((int64_t) fit->keywords.focuspos, int64_t);
	COPY_BE((int64_t) fit->keywords.focussz, int64_t);
	COPY_BE(fit->keywords.foctemp, double);
	COPY_BE(date_ts, int64_t);
	COPY_BE(date_obs_ts, int64_t);
	return 0;
}

static int fits_to_py(fits *fit, unsigned char *ptr, size_t maxlen) {
	if (!fit || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	// Copy numeric values with proper byte order conversion. All
	// types shorter than 64bit are converted to 64bit types before
	// endianness conversion and transmission, to simplify the data
	COPY_BE((int64_t) fit->bitpix, int64_t);
	COPY_BE((int64_t) fit->orig_bitpix, int64_t);
	COPY_BE((uint64_t) fit->checksum, uint64_t);
	COPY_BE(fit->mini, double);
	COPY_BE(fit->maxi, double);
	COPY_BE((double) fit->neg_ratio, double);
	COPY_BE((uint64_t) fit->type, uint64_t);
	COPY_BE((uint64_t) fit->top_down, uint64_t);
	COPY_BE((uint64_t) fit->focalkey, uint64_t);
	COPY_BE((uint64_t) fit->pixelkey, uint64_t);
	COPY_BE((uint64_t) fit->color_managed, uint64_t);
	return 0;
}

static int homography_to_py(const Homography* H, unsigned char *ptr, size_t maxlen) {
	if (!H || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE((double) H->h00, double);
	COPY_BE((double) H->h01, double);
	COPY_BE((double) H->h02, double);
	COPY_BE((double) H->h10, double);
	COPY_BE((double) H->h11, double);
	COPY_BE((double) H->h12, double);
	COPY_BE((double) H->h20, double);
	COPY_BE((double) H->h21, double);
	COPY_BE((double) H->h22, double);
	COPY_BE((int64_t) H->pair_matched, int64_t);
	COPY_BE((int64_t) H->Inliers, int64_t);
	return 0;
}

static int regdata_to_py(const regdata *regparam, unsigned char *ptr, size_t maxlen) {
	if (!regparam || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE((double) regparam->fwhm, double);
	COPY_BE((double) regparam->weighted_fwhm, double);
	COPY_BE((double) regparam->roundness, double);
	COPY_BE(regparam->quality, double);
	COPY_BE((double) regparam->background_lvl, double);
	COPY_BE((int64_t) regparam->number_of_stars, int64_t);
	homography_to_py(&regparam->H, ptr, 11 * sizeof(double));
	return 0;
}

static int imgdata_to_py(const imgdata *imgparam, unsigned char* ptr, size_t maxlen) {
	if (!imgparam || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	int64_t date_obs_ts = imgparam->date_obs ? g_date_time_to_unix(imgparam->date_obs) : 0;
	COPY_BE((int64_t) imgparam->filenum, int64_t);
	COPY_BE((int64_t) imgparam->incl, int64_t);
	COPY_BE((int64_t) date_obs_ts, int64_t);
	COPY_BE((double) imgparam->airmass, double);
	COPY_BE((int64_t) imgparam->rx, int64_t);
	COPY_BE((int64_t) imgparam->ry, int64_t);
	return 0;
}

static int psfstar_to_py(const psf_star *data, unsigned char* ptr, const gchar *units, const gchar *starname, size_t maxlen) {
	if (!data || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE(data->B, double);
	COPY_BE(data->A, double);
	COPY_BE(data->x0, double);
	COPY_BE(data->y0, double);
	COPY_BE(data->sx, double);
	COPY_BE(data->sy, double);
	COPY_BE(data->fwhmx, double);
	COPY_BE(data->fwhmy, double);
	COPY_BE(data->fwhmx_arcsec, double);
	COPY_BE(data->fwhmy_arcsec, double);
	COPY_BE(data->angle, double);
	COPY_BE(data->rmse, double);
	COPY_BE(data->sat, double);
	COPY_BE((int64_t) data->R, int64_t);
	COPY_BE((int64_t) data->has_saturated, int64_t);
	COPY_BE(data->beta, double);
	COPY_BE((int64_t) data->profile, int64_t);
	COPY_BE(data->xpos, double);
	COPY_BE(data->ypos, double);
	COPY_BE(data->mag, double);
	COPY_BE(data->Bmag, double);
	COPY_BE(data->s_mag, double);
	COPY_BE(data->s_Bmag, double);
	COPY_BE(data->SNR, double);
	// photometry *phot not currently passed to python
	// gboolean phot_is_valid not currently passed to python
	COPY_BE(data->BV, double);
	COPY_BE(data->B_err, double);
	COPY_BE(data->A_err, double);
	COPY_BE(data->x_err, double);
	COPY_BE(data->y_err, double);
	COPY_BE(data->sx_err, double);
	COPY_BE(data->sy_err, double);
	COPY_BE(data->ang_err, double);
	COPY_BE(data->beta_err, double);
	COPY_BE((int64_t) data->layer, int64_t);
	COPY_BE(data->ra, double);
	COPY_BE(data->dec, double);
	COPY_STRING(units);
	COPY_STRING(starname);
	return 0;
}

static int seq_to_py(const sequence *seq, unsigned char* ptr, size_t maxlen) {
	if (!seq || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE((int64_t) seq->number, int64_t);
	COPY_BE((int64_t) seq->selnum, int64_t);
	COPY_BE((int64_t) seq->fixed, int64_t);
	COPY_BE((int64_t) seq->nb_layers, int64_t);
	COPY_BE((uint64_t) seq->rx, uint64_t);
	COPY_BE((uint64_t) seq->ry, uint64_t);
	COPY_BE((uint64_t) seq->is_variable, uint64_t);
	COPY_BE((int64_t) seq->bitpix, int64_t);
	COPY_BE((int64_t) seq->reference_image, int64_t);
	COPY_BE((int64_t) seq->beg, int64_t);
	COPY_BE((int64_t) seq->end, int64_t);
	COPY_BE(seq->exposure, double);
	COPY_BE((uint64_t) seq->fz, uint64_t);
	COPY_BE((int64_t) seq->type, int64_t);
	COPY_BE((uint64_t) seq->cfa_opened_monochrome, uint64_t);
	COPY_BE((int64_t) seq->current, int64_t);
	COPY_STRING(seq->seqname);
	// Registration preview coords are not passed to python
	// The dirty and invalid reg flags are not passed to python
	// The photometry data is not currently passed to python
	return 0;
}

static int imstats_to_py(const imstats *stats, unsigned char* ptr, size_t maxlen) {
	if (!stats || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE(stats->total, int64_t);
	COPY_BE(stats->ngoodpix, int64_t);
	COPY_BE(stats->mean, double);
	COPY_BE(stats->median, double);
	COPY_BE(stats->sigma, double);
	COPY_BE(stats->avgDev, double);
	COPY_BE(stats->mad, double);
	COPY_BE(stats->sqrtbwmv, double);
	COPY_BE(stats->location, double);
	COPY_BE(stats->scale, double);
	COPY_BE(stats->min, double);
	COPY_BE(stats->max, double);
	COPY_BE(stats->normValue, double);
	COPY_BE(stats->bgnoise, double);
	return 0;
}

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
	if (payload_length == -1) payload_length = 0;
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
				success = send_response(conn, STATUS_OK, response_data, sizeof(response_data));
			} else {
				// Handle error retrieving dimensions
				const char* error_msg = "Failed to retrieve image dimensions - no image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_PIXELDATA: {
			rectangle region = {0, 0, gfit.rx, gfit.ry};
			success = handle_pixeldata_request(conn, &gfit, region);
			break;
		}

		case CMD_GET_PIXELDATA_REGION: {
			if (payload_length == 16) {
				rectangle region_BE = *(rectangle*) payload;
				rectangle region = {GUINT32_FROM_BE(region_BE.x),
									GUINT32_FROM_BE(region_BE.y),
									GUINT32_FROM_BE(region_BE.w),
									GUINT32_FROM_BE(region_BE.h)};
				success = handle_pixeldata_request(conn, &gfit, region);
			} else {
				g_warning("Unexpected payload length %u received for GET_PIXELDATA_REGION", payload_length);
			}
			break;
		}

/*		case CMD_RELEASE_SHM: {
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
*/
		case CMD_SEND_COMMAND: {
			// Ensure null-terminated string for command
			char* cmd = g_strndup(payload, payload_length);
			int retval = processcommand(cmd, TRUE);
			g_free(cmd);

			// Send response based on command execution
			uint8_t status = (retval == CMD_OK) ? STATUS_OK : STATUS_ERROR;
			success = send_response(conn, status, NULL, 0);
			break;
		}

		case CMD_LOG_MESSAGE: {
			// Ensure null-terminated string for log message
			char* log_msg = g_strndup(payload, payload_length);
			siril_log_message(log_msg);
			g_free(log_msg);

			// Send success response
			success = send_response(conn, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_GET_WORKING_DIRECTORY: {
			// Ensure the working directory is available
			if (com.wd && strlen(com.wd) > 0) {
				// Send success response with the working directory string
				success = send_response(conn, STATUS_OK, com.wd, strlen(com.wd));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = "Working directory not set";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_FILENAME: {
			// Ensure the image filename is available
			if (single_image_is_loaded() && com.uniq && strlen(com.uniq->filename) > 0) {
				// Send success response with the working directory string
				success = send_response(conn, STATUS_OK, com.uniq->filename, strlen(com.uniq->filename));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = "Image not loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_SET_PIXELDATA: {
			if (payload_length != sizeof(incoming_image_info_t)) {
				g_warning("Invalid payload length for SET_PIXELDATA: %u", payload_length);
				const char* error_msg = "Invalid payload length";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			} else {
				success = handle_set_pixeldata_request(conn, &gfit, payload, payload_length);
			}
			break;
		}

		case CMD_GET_IMAGE_STATS: {
			if (payload_length != sizeof(uint32_t)) {
				g_warning("Invalid payload length for GET_IMAGE_STATS: %u", payload_length);
				const char* error_msg = "Invalid payload length";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Get the channel number from payload (convert from network byte order)
			uint32_t channel_BE = *(uint32_t*)payload;
			uint32_t channel = GUINT32_FROM_BE(channel_BE);

			// Check an image is loaded
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Check if channel is valid
			if (channel < 0 || channel >= gfit.naxes[2]) {
				const char* error_msg = "Invalid channel";
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}

			fits *fit = &gfit;

			// If there are no stats available we return NONE, the script can use
			// cmd("stat") to generate them if required
			if (!fit->stats || !fit->stats[channel]) {
				const char* error_message = "No stats";
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
				break;
			}

			// Prepare response buffer with correct byte order
			// Size is 2 longs (16 bytes) + 12 doubles (96 bytes) = 112 bytes
			size_t total_size = 14 * sizeof(double);
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			imstats *stats = fit->stats[channel];
			if (imstats_to_py(stats, ptr, total_size)) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_UPDATE_PROGRESS: {
			if (payload_length < sizeof(float)) {
				g_warning("Invalid payload length for UPDATE_PROGRESS: %u", payload_length);
				const char* error_msg = "Invalid payload length";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
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
			success = send_response(conn, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_GET_KEYWORDS: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t strings_size = FLEN_VALUE * 13;  // 13 string fields of FLEN_VALUE
			size_t numeric_size = sizeof(uint64_t) * 39; // 39 vars packed to 64-bit

			size_t total_size = strings_size + numeric_size;
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;
			if (keywords_to_py(&gfit, ptr, total_size)) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_REGDATA: {
			if (!sequence_is_loaded()) {
				const char* error_msg = "No sequence loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index, chan;
			if (payload_length == 8) {
				index = GUINT32_FROM_BE(*(int*) payload);
				chan = GUINT32_FROM_BE(*((int*) payload + 1));
			}
			if (payload_length != 8 || index < 0 || index >= com.seq.number || chan < 0 || chan > com.seq.nb_layers) {
				const char* error_msg = "Incorrect command arguments";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 6 * sizeof(double) + 11 * sizeof(double);
				// The representation of Homography is 11 x 64-bit vars
			unsigned char *response_buffer = g_try_malloc0(total_size);
			if (!response_buffer) {
				const char* error_msg = "Memory allocation failed";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			unsigned char *ptr = response_buffer;

			regdata *regparam = &com.seq.regparam[index][chan];
			if (regdata_to_py(regparam, ptr, total_size)) {
				const char* error_message = "No regdata available";
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_STATS: {
			if (!sequence_is_loaded()) {
				const char* error_msg = "No sequence loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index, chan;
			if (payload_length == 8) {
				index = GUINT32_FROM_BE(*(int*) payload);
				chan = GUINT32_FROM_BE(*((int*) payload + 1));
			}
			if (payload_length != 8 || index < 0 || index >= com.seq.number || chan < 0 || chan > com.seq.nb_layers) {
				const char* error_msg = "Incorrect command arguments";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 14 * sizeof(double);
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			imstats *stats = com.seq.stats[chan][index];

			if (!com.seq.stats[chan][index]) {
				const char* error_message = "No stats for this channel for this frame";
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
				break;
			}

			if (imstats_to_py(stats, ptr, total_size)) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_IMGDATA: {
			if (!sequence_is_loaded()) {
				const char* error_msg = "No sequence loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || index < 0 || index >= com.seq.number) {
				const char* error_msg = "Incorrect command argument";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 6 * sizeof(double);
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			imgdata *imgparam = &com.seq.imgparam[index];
			if (imgdata_to_py(imgparam, ptr, total_size)) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_PIXELDATA: {
			if (!sequence_is_loaded()) {
				const char* error_msg = "No sequence loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || index < 0 || index >= com.seq.number) {
				const char* error_msg = "Incorrect command argument";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame(&com.seq, index, fit, FALSE, MULTI_THREADED)) {
				free(fit);
				const char* error_msg = "Failed to read sequence frame";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			rectangle region = {0, 0, fit->rx, fit->ry};
			success = handle_pixeldata_request(conn, fit, region);
			clearfits(fit);
			free(fit);
			break;
		}

		case CMD_GET_SEQ_IMAGE: {
			if (!sequence_is_loaded()) {
				const char* error_msg = "No sequence loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || index < 0 || index >= com.seq.number) {
				const char* error_msg = "Incorrect command argument";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame_metadata(&com.seq, index, fit)) {
				free(fit);
				const char* error_msg = "Failed to read frame metadata";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t total_size = sizeof(uint64_t) * 11; // 11 vars packed to 64-bit

			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			int ret = fits_to_py(fit, ptr, total_size);
			clearfits(fit);
			free(fit);
			if (ret) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ: {
			if (!sequence_is_loaded()) {
				const char* error_msg = "No sequence loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			// Calculate size needed for the response
			size_t stringsize = strlen(com.seq.seqname) + 1;
			size_t varsize = sizeof(uint64_t) * 16;
			size_t total_size = varsize + stringsize;
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			if (seq_to_py(&com.seq, ptr, total_size)) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_PSFSTAR: {
			if (!com.stars) {
				const char* error_msg = "No stars list available";
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			// Count the number of stars in com.stars
			int nb_in_com_stars = 0;
			while (com.stars[nb_in_com_stars])
				nb_in_com_stars++;

			if (payload_length != 4 || index < 0) {
				const char* error_msg = "Incorrect command argument";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			if (index >= nb_in_com_stars) { // returning NONE isn't treated as an error, it means
				// we can iterate CMD_CET_PSFSTAR until NONE is returned
				const char* error_msg = "No star with this index";
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			psf_star *psf = com.stars[index];

			// Calculate size needed for the response
			size_t total_size = sizeof(uint64_t) * 36; // 36 variables packed to 64-bit
			gchar *units = NULL, *starname = NULL;
			if (psf->star_name) {
				total_size += strlen(psf->star_name) + 1;
				starname = g_strdup(psf->star_name);
			} else {
				starname = g_strdup("None");
				total_size += strlen(starname) + 1; // star_name may be NULL
			}
			if (psf->units) {
				total_size += strlen(psf->units) + 1;
				units = g_strdup(psf->units);
			} else {
				units = g_strdup("None");
				total_size += strlen(units) + 1;
			}

			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			int ret = psfstar_to_py(psf, ptr, units, starname, total_size);
			g_free(units);
			g_free(starname);
			if (ret) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_IMAGE: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t total_size = sizeof(uint64_t) * 11; // 11 vars packed to 64-bit

			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			if (fits_to_py(&gfit, ptr, total_size)) {
				const char* error_message = "Memory allocation error";
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				break;
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_ICC_PROFILE: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if (gfit.icc_profile == NULL) {
				const char* error_msg = "Image has no ICC profile";
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Prepare data
			guint32 profile_size;
			unsigned char* profile_data = get_icc_profile_data(gfit.icc_profile, &profile_size);

			success = handle_rawdata_request(conn, profile_data, profile_size);
			break;
		}

		case CMD_GET_FITS_HEADER: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = &gfit;
			if (fit->header == NULL) {
				const char* error_msg = "Image has no FITS header";
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Prepare data
			guint32 length = strlen(fit->header) + 1;
			printf("Header length: %u\n", length);
			success = handle_rawdata_request(conn, fit->header, length);
			break;
		}

		case CMD_GET_FITS_HISTORY: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = &gfit;
			if (fit->history == NULL) {
				const char* error_msg = "Image has no history entries";
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Prepare data
			size_t total_length = 0;
			for (GSList *item = fit->history; item != NULL; item = item->next) {
				gchar *str = (gchar *)item->data;
				total_length += strlen(str) + 1;  // +1 to account for the null terminator
			}
			gchar *buffer = malloc(total_length * sizeof(char));
			gchar *ptr = buffer;
			for (GSList *item = fit->history; item != NULL; item = item->next) {
				gchar *str = (gchar *)item->data;
				size_t len = strlen(str) + 1;
				memcpy(ptr, str, len * sizeof(char));
				ptr += len;
			}
			success = handle_rawdata_request(conn, buffer, total_length * sizeof(char));
			break;
		}

		case CMD_GET_FITS_UNKNOWN_KEYS: {
			if (!single_image_is_loaded()) {
				const char* error_msg = "No image loaded";
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = &gfit;
			if (fit->unknown_keys == NULL || fit->unknown_keys[0] == '\0') {
				const char* error_msg = "Image has no unknown keys";
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Prepare data
			guint32 length = strlen(fit->unknown_keys) + 1;

			success = handle_rawdata_request(conn, fit->unknown_keys, length * sizeof(char));
			break;
		}

		default:
			g_warning("Unknown command: %d", header->command);
			const char* error_msg = "Unknown command";
			send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			break;
	}

	if (!success) {
		g_warning("Failed to send response");
	}
}
