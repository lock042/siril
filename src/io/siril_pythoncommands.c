#include "core/siril.h"
#include "algos/PSF.h"
#include "algos/statistics.h"
#include "core/command_line_processor.h"
#include "core/siril_app_dirs.h"
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
			siril_debug_print("Error: Exceeded max length for COPY_FLEN_STRING at %s\n", #str); \
			return 1; \
		} \
		memset(ptr, 0, len); \
		memcpy((char*)ptr, str, len - 1); \
		((char*)ptr)[len - 1] = '\0';     /* Explicitly set null terminator */ \
		ptr += len; \
	}

#define COPY_STRING(str) \
	{ \
		size_t len = strlen(str) + 1; \
		if ((ptr + len) - start_ptr > maxlen) { \
			siril_debug_print("Error: Exceeded max length for COPY_STRING at %s\n", #str); \
			return 1; \
		} \
		memcpy((char*)ptr, str, len);     /* Copy including null terminator */ \
		ptr += len; \
	}

#define COPY_BE(val, type) \
	{ \
		size_t len = sizeof(type); \
		if ((ptr + len) - start_ptr > maxlen) { \
			siril_debug_print("Error: Exceeded max length for COPY_BE at %s\n", #val); \
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

static int psfstar_to_py(const psf_star *data, unsigned char* ptr, size_t maxlen) {
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

static gboolean get_config_value(const char* group, const char* key, config_type_t* type, void** value, size_t* value_size) {
	if (!group || !key || !type || !value || !value_size)
		return FALSE;

	// Get the value using the existing settings function
	const struct settings_access *desc = get_key_settings(group, key);
	if (!desc)
		return FALSE;

	// Example type determination based on key pattern or group
	// You'll need to adjust this logic based on your actual configuration system
	if (desc->type == STYPE_BOOL) {
		uint32_t* bool_val = g_malloc(sizeof(uint32_t));
		*bool_val = GUINT32_TO_BE(*(uint32_t*) desc->data);
		*type = CONFIG_TYPE_BOOL;
		*value = bool_val;
		*value_size = sizeof(uint32_t);
	}
	else if (desc->type == STYPE_INT) {
		gint32* int_val = g_malloc(sizeof(gint32));
		*int_val = GINT32_TO_BE(*(int*) desc->data);
		*type = CONFIG_TYPE_INT;
		*value = int_val;
		*value_size = sizeof(gint32);
	}
	else if (desc->type == STYPE_DOUBLE) {
		// Represented in network byte order
		size_t maxlen = sizeof(double);
		double* double_val = g_malloc(maxlen);
		unsigned char *ptr = (unsigned char*) double_val;
		unsigned char *start_ptr = ptr;
		COPY_BE(*(double*) desc->data, double);
		*type = CONFIG_TYPE_DOUBLE;
		*value = double_val;
		*value_size = sizeof(double);
	}
	else if (desc->type == STYPE_STR || desc->type == STYPE_STRDIR) {
		*type = CONFIG_TYPE_STR;
		*value = g_strdup(*(gchar**) desc->data);
		*value_size = strlen((gchar*) *value) + 1;
	}
	else if (desc->type == STYPE_STRLIST) {
		GSList *list = *((GSList**)desc->data);
		GSList *iter = list;
		size_t total_size = 0;
		while (iter) {
			g_strstrip((gchar*) iter->data); // Remove whitespace
			total_size += strlen((gchar*) iter->data) + 1;
			iter = iter->next;
		}
		total_size += 1; // Final null terminator

		char* list_val = g_malloc(total_size);
		char* ptr = list_val;
		iter = list;
		while (iter) {
			size_t len = strlen((gchar*) iter->data) + 1;
			memcpy(ptr, (gchar*) iter->data, len);
			ptr += len;
			iter = iter->next;
		}
		*ptr = '\0';

		*type = CONFIG_TYPE_STRLIST;
		*value = list_val;
		*value_size = total_size;
	}
	else {
		// Unknown type, report an error
		*type = CONFIG_TYPE_STR;
		*value = g_strdup(_("Unknown config type"));
		*value_size = strlen(*value) + 1;
	}

	return TRUE;
}

typedef struct {
    char shm_name[256];
} finished_shm_payload_t;

/**
* Process the received connection data
*/
void process_connection(Connection* conn, const gchar* buffer, gsize length) {
	if (length < sizeof(CommandHeader)) {
		siril_log_color_message(_("Received incomplete command header\n"), "red");
		return;
	}

	// Get command header
	CommandHeader* header = (CommandHeader*)buffer;
	uint32_t payload_length = GUINT32_FROM_BE(header->length);  // Convert from network byte order
	if (payload_length == -1) payload_length = 0;
	// Verify we have complete message
	if (length < sizeof(CommandHeader) + payload_length) {
		siril_log_color_message(_("Received incomplete command payload: length = %u, expected %u\n"), "red", length, payload_length);
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
				const char* error_msg = _("Failed to retrieve image dimensions - no image loaded");
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
				siril_debug_print(_("Unexpected payload length %u received for GET_PIXELDATA_REGION\n"), payload_length);
			}
			break;
		}

		case CMD_RELEASE_SHM: {
			if (payload_length >= sizeof(finished_shm_payload_t)) {
				finished_shm_payload_t* finished_payload = (finished_shm_payload_t*)payload;
				cleanup_shm_allocation(conn, finished_payload->shm_name);
				// Send acknowledgment
				success = send_response(conn, STATUS_OK, NULL, 0);
			} else {
				const char* error_msg = _("Invalid FINISHED_WITH_SHM payload");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_SEND_COMMAND: {
			// Ensure null-terminated string for command
			char* cmd = g_strndup(payload, payload_length);
			int retval = processcommand(cmd, TRUE);
			g_free(cmd);

			// Send response based on command execution
			if (retval == CMD_OK) {
				success = send_response(conn, STATUS_OK, NULL, 0);
			} else {
				const char* error_msg = _("Siril command error");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_IS_IMAGE_LOADED: {
			int32_t loaded = GINT32_TO_BE((int32_t) single_image_is_loaded());
			success = send_response(conn, STATUS_OK, &loaded, sizeof(int));
			break;
		}

		case CMD_GET_IS_SEQUENCE_LOADED: {
			int32_t loaded = GINT32_TO_BE((int32_t) sequence_is_loaded());
			success = send_response(conn, STATUS_OK, &loaded, sizeof(int));
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
				const char* error_msg = _("Working directory not set");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_USERCONFIG_DIR: {
				const char *configdir = siril_get_config_dir();
			// Ensure the config directory is available
			if (configdir && strlen(configdir) > 0) {
				// Send success response with the working directory string
				success = send_response(conn, STATUS_OK, configdir, strlen(configdir));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = _("Error: user config directory not set");
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
				const char* error_msg = _("Image not loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_SET_PIXELDATA: {
			if (payload_length != sizeof(incoming_image_info_t)) {
				siril_debug_print("Invalid payload length for SET_PIXELDATA: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			} else {
				success = handle_set_pixeldata_request(conn, &gfit, payload, payload_length);
			}
			break;
		}

		case CMD_GET_IMAGE_STATS: {
			if (payload_length != sizeof(uint32_t)) {
				siril_debug_print("Invalid payload length for GET_IMAGE_STATS: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Get the channel number from payload (convert from network byte order)
			uint32_t channel_BE = *(uint32_t*)payload;
			uint32_t channel = GUINT32_FROM_BE(channel_BE);

			// Check an image is loaded
			if (!single_image_is_loaded()) {
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Check if channel is valid
			if (channel >= gfit.naxes[2]) {
				const char* error_msg = _("Invalid channel");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}

			fits *fit = &gfit;

			// If there are no stats available we return NONE, the script can use
			// cmd("stat") to generate them if required
			if (!fit->stats || !fit->stats[channel]) {
				const char* error_message = _("No stats");
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
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_UPDATE_PROGRESS: {
			if (payload_length < sizeof(float)) {
				siril_debug_print("Invalid payload length for UPDATE_PROGRESS: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Use a union to safely handle the byte-order conversion
			union {
				float f;
				uint32_t i;
				unsigned char bytes[sizeof(float)];
			} progress_convert;

			// Copy the bytes from payload
			memcpy(progress_convert.bytes, payload, sizeof(float));

			// Convert from network byte order to host byte order
			progress_convert.i = GUINT32_FROM_BE(progress_convert.i);
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
				const char* error_msg = _("No image loaded");
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
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_REGDATA: {
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index, chan;
			if (payload_length == 8) {
				index = GUINT32_FROM_BE(*(int*) payload);
				chan = GUINT32_FROM_BE(*((int*) payload + 1));
			}
			if (payload_length != 8 || index < 0 || index >= com.seq.number || chan < 0 || chan > com.seq.nb_layers) {
				const char* error_msg = _("Incorrect command arguments");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 6 * sizeof(double) + 11 * sizeof(double);
				// The representation of Homography is 11 x 64-bit vars
			unsigned char *response_buffer = g_try_malloc0(total_size);
			if (!response_buffer) {
				const char* error_msg = _("Memory allocation failed");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			unsigned char *ptr = response_buffer;

			regdata *regparam = &com.seq.regparam[index][chan];
			if (regdata_to_py(regparam, ptr, total_size)) {
				const char* error_message = _("No regdata available");
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_STATS: {
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index, chan;
			if (payload_length == 8) {
				index = GUINT32_FROM_BE(*(int*) payload);
				chan = GUINT32_FROM_BE(*((int*) payload + 1));
			}
			if (payload_length != 8 || index < 0 || index >= com.seq.number || chan < 0 || chan > com.seq.nb_layers) {
				const char* error_msg = _("Incorrect command arguments");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 14 * sizeof(double);
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			imstats *stats = com.seq.stats[chan][index];

			if (!com.seq.stats[chan][index]) {
				g_free(response_buffer);
				const char* error_message = _("No stats for this channel for this frame");
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
				break;
			}

			if (imstats_to_py(stats, ptr, total_size)) {
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_IMGDATA: {
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || index < 0 || index >= com.seq.number) {
				const char* error_msg = _("Incorrect command argument");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 6 * sizeof(double);
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			imgdata *imgparam = &com.seq.imgparam[index];
			if (imgdata_to_py(imgparam, ptr, total_size)) {
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_PIXELDATA: {
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || index < 0 || index >= com.seq.number) {
				const char* error_msg = _("Incorrect command argument");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame(&com.seq, index, fit, FALSE, MULTI_THREADED)) {
				free(fit);
				const char* error_msg = _("Failed to read sequence frame");
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
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || index < 0 || index >= com.seq.number) {
				const char* error_msg = _("Incorrect command argument");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame_metadata(&com.seq, index, fit)) {
				free(fit);
				const char* error_msg = _("Failed to read frame metadata");
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
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ: {
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
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
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_PSFSTARS: {
			if (!com.stars || !com.stars[0]) {
				const char* error_msg = _("No stars list available");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}

			// Count the number of stars in com.stars
			int nb_in_com_stars = 0;
			while (com.stars[nb_in_com_stars])
				nb_in_com_stars++;

			const size_t psf_star_size = 36 * sizeof(double);
			const size_t total_size = nb_in_com_stars * psf_star_size;

			unsigned char* allstars = g_malloc0(total_size);
			if (!allstars) {
				const char* error_msg = _("Memory allocation failed");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			gboolean error_occurred = FALSE;
			for (int i = 0; i < nb_in_com_stars; i++) {
				if(!com.stars) {
					const char* error_msg = _("Stars array was cleared mid-process");
					error_occurred = TRUE;
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				psf_star *psf = com.stars[i];
				if (!psf) {
					error_occurred = TRUE;
					const char* error_msg = _("Unexpected null star entry");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}

				// Calculate the correct offset for this star's data
				unsigned char* ptr = allstars + (i * psf_star_size);

				if (psfstar_to_py(psf, ptr, psf_star_size)) {
					error_occurred = TRUE;
					const char* error_msg = _("Memory allocation error");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
			}

			if (!error_occurred) {
				success = handle_rawdata_request(conn, allstars, total_size);
			}

			g_free(allstars);
			break;
		}

		case CMD_GET_IMAGE: {
			if (!single_image_is_loaded()) {
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t total_size = sizeof(uint64_t) * 11; // 11 vars packed to 64-bit

			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			if (fits_to_py(&gfit, ptr, total_size)) {
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_ICC_PROFILE: {
			if (!single_image_is_loaded()) {
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if (gfit.icc_profile == NULL) {
				const char* error_msg = _("Image has no ICC profile");
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
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = &gfit;
			if (fit->header == NULL) {
				const char* error_msg = _("Image has no FITS header");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Prepare data
			guint32 length = strlen(fit->header) + 1;
			siril_debug_print("Header length: %u\n", length);
			success = handle_rawdata_request(conn, fit->header, length);
			break;
		}

		case CMD_GET_FITS_HISTORY: {
			if (!single_image_is_loaded()) {
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = &gfit;
			if (fit->history == NULL) {
				const char* error_msg = _("Image has no history entries");
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
			g_free(buffer);
			break;
		}

		case CMD_GET_FITS_UNKNOWN_KEYS: {
			if (!single_image_is_loaded()) {
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = &gfit;
			if (fit->unknown_keys == NULL || fit->unknown_keys[0] == '\0') {
				const char* error_msg = _("Image has no unknown keys");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Prepare data
			guint32 length = strlen(fit->unknown_keys) + 1;

			success = handle_rawdata_request(conn, fit->unknown_keys, length * sizeof(char));
			break;
		}

		case CMD_GET_CONFIG: {
			if (payload_length < 2) {  // Need at least two null-terminated strings
				const char* error_msg = _("Group and key must be provided");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Extract group and key from payload
			const char* group = (const char*)payload;
			size_t group_len = strlen(group);
			if (group_len >= payload_length - 1) {
				const char* error_msg = _("Malformed group/key data");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			const char* key = (const char*)payload + group_len + 1;

			config_type_t type;
			void* value;
			size_t value_size;

			if (!get_config_value(group, key, &type, &value, &value_size)) {
				const char* error_msg = _("Configuration key not found");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Allocate response buffer: 1 byte for type + value size
			size_t total_size = 1 + value_size;
			unsigned char* response_buffer = g_malloc(total_size);

			// Write type and value
			response_buffer[0] = (unsigned char)type;
			memcpy(response_buffer + 1, value, value_size);

			// Send response
			success = send_response(conn, STATUS_OK, response_buffer, total_size);

			// Clean up
			g_free(value);
			g_free(response_buffer);
			break;
		}

		default:
			siril_debug_print("Unknown command: %d\n", header->command);
			const char* error_msg = _("Unknown command");
			success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			break;
	}

	if (!success) {
		siril_debug_print("Failed to send response\n");
	}
}