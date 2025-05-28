#include <string.h>
#include <math.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "algos/geometry.h"
#include "gps_parser.h"


#if 0
struct _dvti_struct {
	BYTE preamble[4];
	BYTE overlay_version[2];
	BYTE firmware_version[2];
	BYTE start_timestamp[12];
	BYTE end_timestamp[12];
	BYTE reserved[4];
	uint32_t frame_number;
	uint16_t frame_width;
	uint16_t frame_height;
	BYTE xbinning;	// minus one
	BYTE ybinning;	// minus one
	BYTE reserved2[2];
	uint16_t signal_multiplier;
	uint16_t signal_divider;
	uint32_t exposure_seconds;
	uint32_t exposure_nanoseconds;
	uint16_t gain_mode;	// 0 for low gain
	uint16_t black_level;
	int32_t latitude;	// times 1e7, positive is N
	int32_t longitude;	// times 1e7, positive is E
	uint16_t position_accuracy;	// (meters)
	BYTE fix_state;
	BYTE number_of_stats;
	int32_t altitude;	// ASL (mm)
	BYTE camera_status_flags[4];
	BYTE camera_status_info[4];
	BYTE usb_phy_info[4];
	BYTE usb_link_info[4];
	BYTE sensor_temperature[2];
	BYTE gps_dop_data[10];
	BYTE reserved3[14];
	BYTE checksum[2];
};
#endif

struct __attribute__((packed)) _qhy_struct_raw {
	// big endian
	uint32_t sequence_number;
	BYTE reserved;
	uint16_t image_width;
	uint16_t image_height;
	uint32_t latitude;	// format is SDDMMMMMMM (sign, degrees, minutes)
	uint32_t longitude;
	BYTE start_flag;	// GPS status
	uint32_t start_seconds;
	BYTE start_microseconds[3];
	BYTE end_flag;	// GPS status
	uint32_t end_seconds;
	BYTE end_microseconds[3];
	BYTE now_flag;	// GPS status
	uint32_t now_seconds;
	BYTE now_microseconds[3];
	BYTE count_of_PPS[3];
};

static const char *receiver_status(BYTE flag, gboolean shifted) {
	BYTE receiver_flag = shifted ? flag : (flag & 48) >> 4;
	switch (receiver_flag) {
		case 0: return "just powered on";
		case 1: return "not locked";
		case 2: return "not locked but data valid";	// position and time
		case 3: return "locked and valid";
	}
	return "?";
}

static gboolean time_is_accurate(int shifted_flag) {
	return shifted_flag == 3 || shifted_flag == 2; // data->flag is already shifted
}

void print_qhy_data(struct _qhy_struct *qhy) {
	siril_log_message("sequence number: %u\n", qhy->sequence_number);
	siril_log_message("image dimensions: %hu x %hu\n", qhy->image_width, qhy->image_height);
	siril_log_message("lat/long: %f, %f\n", qhy->latitude, qhy->longitude);

	siril_log_message("start flag: receiver is %s\n",
			receiver_status(qhy->start_flag, qhy->flags_are_shifted));
	gchar *date_str = date_time_to_FITS_date(qhy->start);
	siril_log_message("start date: %s\n", date_str);
	g_free(date_str);

	siril_log_message("end flagd: receiver is %s\n",
			receiver_status(qhy->end_flag, qhy->flags_are_shifted));
	date_str = date_time_to_FITS_date(qhy->end);
	siril_log_message("end date: %s\n", date_str);
	g_free(date_str);

	GTimeSpan exposure = g_date_time_difference(qhy->end, qhy->start);
	double exp_sec = exposure / 1000000.0;
	siril_log_message("measured exposure (end - start): %f seconds\n", exp_sec);

	if (qhy->now) {
		siril_log_message("now flag: receiver is %s\n",
				receiver_status(qhy->now_flag, qhy->flags_are_shifted));
		date_str = date_time_to_FITS_date(qhy->now);
		siril_log_message("now date: %s\n", date_str);
		g_free(date_str);
	}

	if (qhy->count_of_PPS) {
		siril_log_message("PPS count: %u (%s)\n", qhy->count_of_PPS,
				qhy->count_of_PPS >= 10000500 ? "bad" : "ok");
	}
}

/* will read GPS metadata from QHY first pixel row, will not modify the image */
int parse_gps_image(fits *fit, struct _qhy_struct *qhy_header) {
	/* we read the metadata from the extra row in the image, not from text drawn in
	 * pixels, and we update the FITS header with it */
	if (fit->type != DATA_USHORT || !fit->data || fit->naxes[0] < 128) {
		siril_log_message("GPS metadata can only be extracted on RAW images\n");
		return -1;
	}
#if 0
	BYTE dvti_cookie8bit[] = { 0x00, 0x00, 0x54, 0x53 };
	BYTE dvti_cookie16bit[] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x54, 0x00, 0x53 };
	if ((fit->bitpix == BYTE_IMG && !memcmp(fit->data, dvti_cookie8bit, 4)) ||
			!memcmp(fit->data, dvti_cookie16bit, 8)) {
		// DVTI
		if (fit->bitpix == BYTE_IMG) {
			struct _dvti_struct header;
			memcpy(&header, fit->data, sizeof(struct _dvti_struct));
		} else {
			// relou
		}

		return 0;
	}
#endif

	struct _qhy_struct_raw raw_header = { 0 };
	memcpy(&raw_header, fit->data, sizeof raw_header);
	siril_debug_print("header size: %zd\n", sizeof raw_header);
	qhy_header->sequence_number = be32_to_cpu(raw_header.sequence_number);
	qhy_header->image_width = be16_to_cpu(raw_header.image_width);
	qhy_header->image_height = be16_to_cpu(raw_header.image_height);

	uint32_t lat = be32_to_cpu(raw_header.latitude);
	qhy_header->latitude = (lat % 1000000000) / 10000000;
	qhy_header->latitude += (lat % 10000000) / 6000000.0;
	if (lat > 1000000000) qhy_header->latitude = -qhy_header->latitude;

	uint32_t lon = be32_to_cpu(raw_header.longitude);
	qhy_header->longitude = (lon % 1000000000) / 1000000;
	qhy_header->longitude += (lon % 1000000) / 600000.0;
	if (lon > 1000000000) qhy_header->longitude = -qhy_header->longitude;
	siril_debug_print("lat: %u\tlon: %u\n", lat, lon);

	qhy_header->start_flag = raw_header.start_flag;
	uint32_t js_start = be32_to_cpu(raw_header.start_seconds);
	uint32_t us = be24_to_cpu(raw_header.start_microseconds) / 10;
	qhy_header->start = julian_sec_to_date(js_start, us);

	qhy_header->end_flag = raw_header.end_flag;
	uint32_t js_end = be32_to_cpu(raw_header.end_seconds);
	us = be24_to_cpu(raw_header.end_microseconds) / 10;
	qhy_header->end = julian_sec_to_date(js_end, us);

	qhy_header->now_flag = raw_header.now_flag;
	uint32_t js_now = be32_to_cpu(raw_header.now_seconds);
	us = be24_to_cpu(raw_header.now_microseconds) / 10;
	qhy_header->now = julian_sec_to_date(js_now, us);

	qhy_header->count_of_PPS = be24_to_cpu(raw_header.count_of_PPS);

	/* diagnose the validity of the metadata */
	GTimeSpan offset_s = 0;
	if (fit->keywords.date_obs) {
		GTimeSpan offset = g_date_time_difference(qhy_header->start, fit->keywords.date_obs);
		offset_s = offset / 1000000;	// rounded seconds
	}

	/* SharpCap uses the GPS time for DATE-OBS and the system clock as DATE-END, it handles
	 * loss of lock by setting: DATE-OBS= '0001-01-01T00:00:00.0000000'
	 */
	if (fit->keywords.date_obs && g_date_time_get_year(fit->keywords.date_obs) == 1 &&
			g_date_time_get_day_of_year(fit->keywords.date_obs) == 1 &&
			g_date_time_get_hour(fit->keywords.date_obs) == 0 &&
			g_date_time_get_minute(fit->keywords.date_obs) == 0 &&
			g_date_time_get_seconds(fit->keywords.date_obs) == 0.0)
		offset_s = 0;	// disable this check
	if (qhy_header->image_width > 12000 || qhy_header->image_height > 9000 ||
			qhy_header->longitude < -180.0 || qhy_header->longitude > 360.0 ||
			qhy_header->latitude < -90.0 || qhy_header->latitude > 90.0 ||
			js_end - js_start > round_to_int(fit->keywords.exposure) + 1 ||
			js_end - js_start < 0 ||
			offset_s < -900 || offset_s > 900) { // allow 15 minutes of clock offset
		siril_log_message("Extracted metadata seems incorrect, assuming this is not an image containing some\n");
		return 1;
	}
	return 0;
}

int siril_fits_open_diskfile_img(fitsfile **fptr, const char *filename, int iomode, int *status);

/* read GPS_* keywords from the FITS header, useful for global shutter GPS camera images
 * written by NINA (see the list below).
 * prerequisites: fit->fptr is still open, fit->naxes and fit->keywords.date_obs are initialized
 */
int parse_gps_from_header(fits *fit, const char *filename, struct _qhy_struct *qhy_header) {
	char date[FLEN_VALUE] = { 0 };
	int status = 0, s = 0, flag = 0;
	if (filename) {
		siril_debug_print("reading GPS data from header for %s\n", filename);
		siril_fits_open_diskfile_img(&(fit->fptr), filename, READONLY, &status);
	}
	fits_read_key(fit->fptr, TINT, "GPS_W", &qhy_header->image_width, NULL, &status);
	fits_read_key(fit->fptr, TINT, "GPS_H", &qhy_header->image_height, NULL, &status);
	fits_read_key(fit->fptr, TDOUBLE, "GPS_LAT", &qhy_header->latitude, NULL, &status);
	fits_read_key(fit->fptr, TDOUBLE, "GPS_LON", &qhy_header->longitude, NULL, &status);
	if (!fits_read_key(fit->fptr, TINT, "GPS_SFLG", &flag, NULL, &status))
		qhy_header->start_flag = flag;
	if (!fits_read_key(fit->fptr, TINT, "GPS_EFLG", &flag, NULL, &status))
		qhy_header->end_flag = flag;
	qhy_header->flags_are_shifted = TRUE;
	if (!fits_read_key(fit->fptr, TSTRING, "GPS_SUTC", &date, NULL, &status))
		qhy_header->start = FITS_date_to_date_time(date);
	if (!fits_read_key(fit->fptr, TSTRING, "GPS_EUTC", &date, NULL, &status))
		qhy_header->end = FITS_date_to_date_time(date);
	if (filename)
		fits_close_file(fit->fptr, &s);
	if (status || !qhy_header->start || !qhy_header->end) {
		siril_debug_print("Image does not have the expected GPS keywords\n");
		return -1;
	}

	/* sanity checks */
	GTimeSpan offset, offset_s = 0;
	if (fit->keywords.date_obs) {
		offset = g_date_time_difference(qhy_header->start, fit->keywords.date_obs);
		offset_s = offset / 1000000;	// rounded seconds
	}
	offset = g_date_time_difference(qhy_header->end, qhy_header->start);
	double calculated_exposure = offset / 1000000.0;

	/* SharpCap uses the GPS time for DATE-OBS and the system clock as DATE-END, it handles
	 * loss of lock by setting: DATE-OBS= '0001-01-01T00:00:00.0000000'
	 */
	if (fit->keywords.date_obs && g_date_time_get_year(fit->keywords.date_obs) == 1 &&
			g_date_time_get_day_of_year(fit->keywords.date_obs) == 1 &&
			g_date_time_get_hour(fit->keywords.date_obs) == 0 &&
			g_date_time_get_minute(fit->keywords.date_obs) == 0 &&
			g_date_time_get_seconds(fit->keywords.date_obs) == 0.0)
		offset_s = 0;	// disable this check
	if (qhy_header->image_width > 12000 || qhy_header->image_height > 9000 ||
			qhy_header->longitude < -180.0 || qhy_header->longitude > 360.0 ||
			qhy_header->latitude < -90.0 || qhy_header->latitude > 90.0 ||
			calculated_exposure > round_to_int(fit->keywords.exposure) + 1 ||
			calculated_exposure < 0 ||
			offset_s < -900 || offset_s > 900) { // allow 15 minutes of clock offset
		siril_log_message("Extracted metadata seems incorrect, ignoring\n");
		return 1;
	}
	if (calculated_exposure == 0 && fit->keywords.exposure > 0.001) {
		siril_log_message("Exposures don't match: %.3f from GPS, %f in header\n", calculated_exposure, fit->keywords.exposure);
		return 1;
	}
	if (!time_is_accurate(qhy_header->start_flag) || !time_is_accurate(qhy_header->end_flag)) {
		siril_log_message("Time is reported inaccurate, not using GPS data\n");
		return 1;
	}
	return 0;
}

/* global shutter metadata extracted by N.I.N.A.:
 * log: SWCREATE= 'N.I.N.A. 3.0.0.2007 (x64)' / Software that created this file
 * log: GPS_PPS =             10000500 / QHY pps
 * log: GPS_SEQ =              1890814 / QHY sequence nr
 * log: GPS_W   =                 1936 / QHY width
 * log: GPS_H   =                 1024 / QHY height
 * log: GPS_LAT =                   0. / latitude
 * log: GPS_LON =                   0. / longitude
 * log: GPS_SFLG=                    2 / QHY start_flag
 * log: GPS_SST = 'not locked but data valid' / QHY start_flag status
 * log: GPS_SSEC=            892554211 / [s] QHY start
 * log: GPS_SUS =               923882 / [us] QHY start
 * log: GPS_SUTC= '2024-01-21T11:43:31.923' / QHY start_time
 * log: GPS_EFLG=                    2 / QHY end_flag
 * log: GPS_EST = 'not locked but data valid' / QHY end_flag status
 * log: GPS_ESEC=            892554211 / [s] QHY end
 * log: GPS_EUS =               923903 / [us] QHY end
 * log: GPS_EUTC= '2024-01-21T11:43:31.923' / QHY end_time
 * log: GPS_NFLG=                    2 / QHY now_flag
 * log: GPS_NST = 'not locked but data valid' / QHY now_flag status
 * log: GPS_NSEC=            892554211 / [s] QHY now
 * log: GPS_NUS =               923903 / [us] QHY now
 * log: GPS_NUTC= '2024-01-21T11:43:31.923' / QHY now_time
 * log: GPS_EXP =                   21 / [us] QHY exposure
 */

void release_qhy_struct(struct _qhy_struct *data) {
	if (data->start)
		g_date_time_unref(data->start);
	if (data->end)
		g_date_time_unref(data->end);
	if (data->now)
		g_date_time_unref(data->now);
}

int update_fit_from_qhy_header(fits *fit, struct _qhy_struct *qhy_header) {
	gboolean rolling_shutter = (strstr(fit->keywords.instrume, "QHY268") ||
			strstr(fit->keywords.instrume, "QHY600")); // other exist but were not tested
	if (rolling_shutter)
		return 0;

	GTimeSpan exposure = g_date_time_difference(qhy_header->end, qhy_header->start);
	double real_exposure = exposure / 1000000.0;
	if (real_exposure > 900.0) {
		siril_log_message("Extracted metadata seems incorrect, assuming this is not an image containing some\n");
		return 1;
	}

	if (real_exposure != 0.0) {
		if (fabs((fit->keywords.exposure + real_exposure) / (2.0 * real_exposure) - 1.0) > 0.05)
			siril_log_message("Discrepancy in configured (%.3f) and actual exposure (%.4f)\n", fit->keywords.exposure, real_exposure);
	}
	fit->keywords.exposure = real_exposure;

	if (fit->keywords.date_obs && qhy_header->start) {
		gchar *date = date_time_to_date(fit->keywords.date_obs);
		double time1 = get_decimal_hours(fit->keywords.date_obs);
		double time2 = get_decimal_hours(qhy_header->start);
		siril_debug_print("D,%s,%f,%f,%f\n", date, time1, time2, (time1 - time2) * 3600.0);
		g_free(date);
		g_date_time_unref(fit->keywords.date_obs);
		fit->keywords.date_obs = g_date_time_ref(qhy_header->start);
	}
	else if (qhy_header->start) {
		fit->keywords.date_obs = g_date_time_ref(qhy_header->start);
	}

	siril_debug_print("DATE-OBS and EXPTIME overwritten by GPS metadata\n");
	fit->keywords.date_and_exp_from_gps = TRUE;
	fit->history = g_slist_append(fit->history, strdup("DATE-OBS and EXPTIME overwritten by GPS metadata"));
	return 0;
}

static gboolean readout_is_by_row_pairs(struct gps_rs_data *data) {	// = non-2CMS
	/* QHY cameras have several readout modes: Photographic DSO, High Gain, Extended Fullwell.
	 * As designed for color cameras, they read rows by pairs to keep colors aligned, but this
	 * is kept in monochrome camera. Some high-end cameras (268, 461, 600) support additional
	 * modes based on those that readout the sensor row by row, they are called 2CMS modes.
	 */
	return !strstr(data->readout_mode, "2CMS");
}

/* we currently only compute the timestamp for the row, the change in time for a pixel in column is much
 * less than 1 ms and the meaning of pixel_period is unclear
 * Y is in display coordinates
 */
GDateTime *get_timestamp_for_pixel(struct gps_rs_data *data, enum timestamp_type type, int x, int y) {
	if (!data) return NULL;
	if (!time_is_accurate(data->flag)) {
		siril_log_message("GPS timestamp will not be used because it's not locked (flag=%d)\n", data->flag);
		return NULL;
	}
	/* compute Y in bin1 FITS coordinates */
	if (data->top_down) // bottom-up images are displayed with the same coordinates as the sensor
		y = data->ry - y - 1;
	y += data->crop_offset_y;	// from siril crops, in binned coordinates
	y *= data->binning;
	siril_debug_print("FITS y: %d (ry: %d)\n", y, data->ry);

	/* compute offset */
	double y_offset;
	if (readout_is_by_row_pairs(data))
		y_offset = 2.0 * data->line_period * (y / 2) * 1e-3; // µs
	else y_offset = data->line_period * y * 1e-3; // µs

	double offset_seconds = (y_offset + data->end_offset0) * 1e-6;
	if (type == EXP_START)
		offset_seconds -= data->exposure;
	else if (type == EXP_MIDDLE)
		offset_seconds -= data->exposure * 0.5;
	/* result */
	return g_date_time_add_seconds(data->end_vsync_date, offset_seconds);
}

// to call before changing fit->ry
void apply_crop_to_gps_data(fits *fit, rectangle *bounds) {
	if (!fit->keywords.gps_data)
		return;
	fit->keywords.gps_data->crop_offset_x += bounds->x;
	// crop coordinates are in display coordinates, here we need FITS coordinates
	if (fit->keywords.gps_data->top_down)
		fit->keywords.gps_data->crop_offset_y += fit->ry - bounds->y - bounds->h;
	else fit->keywords.gps_data->crop_offset_y += bounds->y;
	fit->keywords.gps_data->ry = bounds->h;
}

// to call after changing fit->row_order
void apply_flip_to_gps_data(fits *fit) {
	if (!fit->keywords.gps_data)
		return;
	fit->keywords.gps_data->top_down = !g_strcmp0(fit->keywords.row_order, "TOP-DOWN");
	siril_debug_print("Image is now %s\n", fit->keywords.row_order);
}

// to call after changing fit->binning_x
void apply_binning_to_gps_data(fits *fit) {
	if (!fit->keywords.gps_data)
		return;
	if (fit->keywords.gps_data->crop_offset_x || fit->keywords.gps_data->crop_offset_y) {
		fit->keywords.gps_data->crop_offset_x = fit->keywords.gps_data->crop_offset_x * fit->keywords.gps_data->binning / fit->keywords.binning_x;
		fit->keywords.gps_data->crop_offset_y = fit->keywords.gps_data->crop_offset_y * fit->keywords.gps_data->binning / fit->keywords.binning_x;
	}
	siril_debug_print("updated binning value for GPS metadata from %d to %d\n",
			fit->keywords.gps_data->binning, fit->keywords.binning_x);
	fit->keywords.gps_data->binning = fit->keywords.binning_x;
	fit->keywords.gps_data->ry = fit->ry;
}

struct gps_rs_data *clone_gps_data(struct gps_rs_data *in) {
	if (!in) return NULL;
	struct gps_rs_data *out = malloc(sizeof(struct gps_rs_data));
	if (!out) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	memcpy(out, in, sizeof(struct gps_rs_data));
	if (in->end_vsync_date)
		out->end_vsync_date = g_date_time_ref(in->end_vsync_date);
	return out;
}



/* gets QHY GPS metadata from first image row and updates the header with extracted data */
int gps_extract_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
                rectangle *_, int threads) {
        int retval = 1;
        struct _qhy_struct qhy_header = { 0 };
        if (!parse_gps_image(fit, &qhy_header)) {
                int nb_rows_to_remove = GPOINTER_TO_INT(args->user);
                if (nb_rows_to_remove < 0)
                        nb_rows_to_remove = 6;
                retval = update_fit_from_qhy_header(fit, &qhy_header);
                release_qhy_struct(&qhy_header);
                // crop the metadata lines
                if (!retval) {
                        rectangle new_dimensions = { .x = 0, .y = 0,
                                .w = fit->rx, .h = fit->ry - nb_rows_to_remove };
                        if (g_strcmp0(fit->keywords.row_order, "TOP-DOWN"))
                                new_dimensions.y = nb_rows_to_remove;
                        retval = crop(fit, &new_dimensions);
                }
        }
        return retval;
}

