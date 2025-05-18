// Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
// Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
// Reference site is https://siril.org
// SPDX-License-Identifier: GPL-3.0-or-later

#include "core/siril.h"
#include "algos/background_extraction.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "algos/statistics.h"
#include "core/command_line_processor.h"
#include "core/siril_app_dirs.h"
#include "core/icc_profile.h"
#include "core/siril_log.h"
#include "core/proto.h"
#include "core/undo.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/user_polygons.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/siril_pythonmodule.h"
#include "filters/synthstar.h"

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

#define COPY_BE64(val, type) \
	{ \
		size_t len = sizeof(type); \
		if ((ptr + len) - start_ptr > maxlen) { \
			siril_debug_print("Error: Exceeded max length for COPY_BE64 at %s\n", #val); \
			return 1; \
		} \
		union { type v; uint64_t i; } conv; \
		conv.v = val; \
		uint64_t be = GUINT64_TO_BE(conv.i); \
		memcpy(ptr, &be, len); \
		ptr += len; \
	}

#define FROM_BE64_INTO(dest, val, type) \
do { \
	union { type v; uint64_t i; } conv; \
	memcpy(&conv.i, &val, sizeof(type)); \
	conv.i = GUINT64_FROM_BE(conv.i); \
	(dest) = conv.v; \
} while(0)

#define TO_BE64_INTO(dest, val, type) \
	do { \
		union { type v; uint64_t i; } conv; \
		conv.v = val; \
		conv.i = GUINT64_TO_BE(conv.i); \
		union { type v; uint64_t i; } result; \
		result.i = conv.i; \
		(dest) = result.v; \
	} while(0)

#define BOOL_FROM_BYTE(x) ((x) ? TRUE : FALSE)

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
	COPY_BE64(1.0, double); // aligned to the code in fits_keywords.c
	COPY_BE64(0.0, double); // aligned to the code in fits_keywords.c
	COPY_BE64((uint64_t) fit->keywords.lo, uint64_t);
	COPY_BE64((uint64_t) fit->keywords.hi, uint64_t);
	COPY_BE64((double) fit->keywords.flo, double);
	COPY_BE64((double) fit->keywords.fhi, double);
	COPY_BE64(fit->keywords.data_max, double);
	COPY_BE64(fit->keywords.data_min, double);
	COPY_BE64(fit->keywords.pixel_size_x, double);
	COPY_BE64(fit->keywords.pixel_size_y, double);
	COPY_BE64((uint64_t) fit->keywords.binning_x, uint64_t);
	COPY_BE64((uint64_t) fit->keywords.binning_y, uint64_t);
	COPY_BE64(fit->keywords.expstart, double);
	COPY_BE64(fit->keywords.expend, double);
	COPY_BE64(fit->keywords.centalt, double);
	COPY_BE64(fit->keywords.centaz, double);
	COPY_BE64(fit->keywords.sitelat, double);
	COPY_BE64(fit->keywords.sitelong, double);
	COPY_BE64(fit->keywords.siteelev, double);
	COPY_BE64((int64_t) fit->keywords.bayer_xoffset, int64_t);
	COPY_BE64((int64_t) fit->keywords.bayer_yoffset, int64_t);
	COPY_BE64(fit->keywords.airmass, double);
	COPY_BE64(fit->keywords.focal_length, double);
	COPY_BE64(fit->keywords.flength, double);
	COPY_BE64(fit->keywords.iso_speed, double);
	COPY_BE64(fit->keywords.exposure, double);
	COPY_BE64(fit->keywords.aperture, double);
	COPY_BE64(fit->keywords.ccd_temp, double);
	COPY_BE64(fit->keywords.set_temp, double);
	COPY_BE64(fit->keywords.livetime, double);
	COPY_BE64((uint64_t) fit->keywords.stackcnt, uint64_t);
	COPY_BE64(fit->keywords.cvf, double);
	COPY_BE64((int64_t) fit->keywords.key_gain, int64_t);
	COPY_BE64((int64_t) fit->keywords.key_offset, int64_t);
	COPY_BE64((int64_t) fit->keywords.focuspos, int64_t);
	COPY_BE64((int64_t) fit->keywords.focussz, int64_t);
	COPY_BE64(fit->keywords.foctemp, double);
	COPY_BE64(date_ts, int64_t);
	COPY_BE64(date_obs_ts, int64_t);
	return 0;
}

static int fits_to_py(fits *fit, unsigned char *ptr, size_t maxlen) {
	if (!fit || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	// Copy numeric values with proper byte order conversion. All
	// types shorter than 64bit are converted to 64bit types before
	// endianness conversion and transmission, to simplify the data
	COPY_BE64((int64_t) fit->rx, int64_t);
	COPY_BE64((int64_t) fit->ry, int64_t);
	COPY_BE64((int64_t) fit->naxes[2], int64_t);
	COPY_BE64((int64_t) fit->bitpix, int64_t);
	COPY_BE64((int64_t) fit->orig_bitpix, int64_t);
	COPY_BE64((uint64_t) fit->checksum, uint64_t);
	COPY_BE64(fit->mini, double);
	COPY_BE64(fit->maxi, double);
	COPY_BE64((double) fit->neg_ratio, double);
	COPY_BE64((uint64_t) fit->top_down, uint64_t);
	COPY_BE64((uint64_t) fit->focalkey, uint64_t);
	COPY_BE64((uint64_t) fit->pixelkey, uint64_t);
	COPY_BE64((uint64_t) fit->color_managed, uint64_t);
	return 0;
}

static int homography_to_py(const Homography* H, unsigned char *ptr, size_t maxlen) {
	if (!H || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE64((double) H->h00, double);
	COPY_BE64((double) H->h01, double);
	COPY_BE64((double) H->h02, double);
	COPY_BE64((double) H->h10, double);
	COPY_BE64((double) H->h11, double);
	COPY_BE64((double) H->h12, double);
	COPY_BE64((double) H->h20, double);
	COPY_BE64((double) H->h21, double);
	COPY_BE64((double) H->h22, double);
	COPY_BE64((int64_t) H->pair_matched, int64_t);
	COPY_BE64((int64_t) H->Inliers, int64_t);
	return 0;
}

static int sample_to_py(const background_sample *sample, unsigned char *ptr, size_t maxlen) {
	if (!sample || !ptr)
		return 1;
	unsigned char *start_ptr = ptr;

	COPY_BE64(sample->median[0], double);
	COPY_BE64(sample->median[1], double);
	COPY_BE64(sample->median[2], double);
	COPY_BE64(sample->mean, double);
	COPY_BE64(sample->min, double);
	COPY_BE64(sample->max, double);
	COPY_BE64((uint64_t) sample->size, uint64_t);
	COPY_BE64(sample->position.x, double);
	COPY_BE64(sample->position.y, double);
	COPY_BE64((uint64_t) sample->valid, uint64_t);
	return 0;
}

static int regdata_to_py(const regdata *regparam, unsigned char *ptr, size_t maxlen) {
	if (!regparam || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE64((double) regparam->fwhm, double);
	COPY_BE64((double) regparam->weighted_fwhm, double);
	COPY_BE64((double) regparam->roundness, double);
	COPY_BE64(regparam->quality, double);
	COPY_BE64((double) regparam->background_lvl, double);
	COPY_BE64((int64_t) regparam->number_of_stars, int64_t);
	homography_to_py(&regparam->H, ptr, 11 * sizeof(double));
	return 0;
}

static int imgdata_to_py(const imgdata *imgparam, unsigned char* ptr, size_t maxlen) {
	if (!imgparam || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	int64_t date_obs_ts = imgparam->date_obs ? g_date_time_to_unix(imgparam->date_obs) : 0;
	COPY_BE64((int64_t) imgparam->filenum, int64_t);
	COPY_BE64((int64_t) imgparam->incl, int64_t);
	COPY_BE64((int64_t) date_obs_ts, int64_t);
	COPY_BE64((double) imgparam->airmass, double);
	COPY_BE64((int64_t) imgparam->rx, int64_t);
	COPY_BE64((int64_t) imgparam->ry, int64_t);
	return 0;
}

static int psfstar_to_py(const psf_star *data, unsigned char* ptr, size_t maxlen) {
	if (!data || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE64(data->B, double);
	COPY_BE64(data->A, double);
	COPY_BE64(data->x0, double);
	COPY_BE64(data->y0, double);
	COPY_BE64(data->sx, double);
	COPY_BE64(data->sy, double);
	COPY_BE64(data->fwhmx, double);
	COPY_BE64(data->fwhmy, double);
	COPY_BE64(data->fwhmx_arcsec, double);
	COPY_BE64(data->fwhmy_arcsec, double);
	COPY_BE64(data->angle, double);
	COPY_BE64(data->rmse, double);
	COPY_BE64(data->sat, double);
	COPY_BE64((int64_t) data->R, int64_t);
	COPY_BE64((int64_t) data->has_saturated, int64_t);
	COPY_BE64(data->beta, double);
	COPY_BE64((int64_t) data->profile, int64_t);
	COPY_BE64(data->xpos, double);
	COPY_BE64(data->ypos, double);
	COPY_BE64(data->mag, double);
	COPY_BE64(data->Bmag, double);
	COPY_BE64(data->s_mag, double);
	COPY_BE64(data->s_Bmag, double);
	COPY_BE64(data->SNR, double);
	// photometry *phot not currently passed to python
	// gboolean phot_is_valid not currently passed to python
	COPY_BE64(data->BV, double);
	COPY_BE64(data->B_err, double);
	COPY_BE64(data->A_err, double);
	COPY_BE64(data->x_err, double);
	COPY_BE64(data->y_err, double);
	COPY_BE64(data->sx_err, double);
	COPY_BE64(data->sy_err, double);
	COPY_BE64(data->ang_err, double);
	COPY_BE64(data->beta_err, double);
	COPY_BE64((int64_t) data->layer, int64_t);
	COPY_BE64(data->ra, double);
	COPY_BE64(data->dec, double);
	return 0;
}

static int seq_to_py(const sequence *seq, unsigned char* ptr, size_t maxlen) {
	if (!seq || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE64((int64_t) seq->number, int64_t);
	COPY_BE64((int64_t) seq->selnum, int64_t);
	COPY_BE64((int64_t) seq->fixed, int64_t);
	COPY_BE64((int64_t) seq->nb_layers, int64_t);
	COPY_BE64((uint64_t) seq->rx, uint64_t);
	COPY_BE64((uint64_t) seq->ry, uint64_t);
	COPY_BE64((uint64_t) seq->is_variable, uint64_t);
	COPY_BE64((int64_t) seq->bitpix, int64_t);
	COPY_BE64((int64_t) seq->reference_image, int64_t);
	COPY_BE64((int64_t) seq->beg, int64_t);
	COPY_BE64((int64_t) seq->end, int64_t);
	COPY_BE64(seq->exposure, double);
	COPY_BE64((uint64_t) seq->fz, uint64_t);
	COPY_BE64((int64_t) seq->type, int64_t);
	COPY_BE64((uint64_t) seq->cfa_opened_monochrome, uint64_t);
	COPY_BE64((int64_t) seq->current, int64_t);
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

	COPY_BE64(stats->total, int64_t);
	COPY_BE64(stats->ngoodpix, int64_t);
	COPY_BE64(stats->mean, double);
	COPY_BE64(stats->median, double);
	COPY_BE64(stats->sigma, double);
	COPY_BE64(stats->avgDev, double);
	COPY_BE64(stats->mad, double);
	COPY_BE64(stats->sqrtbwmv, double);
	COPY_BE64(stats->location, double);
	COPY_BE64(stats->scale, double);
	COPY_BE64(stats->min, double);
	COPY_BE64(stats->max, double);
	COPY_BE64(stats->normValue, double);
	COPY_BE64(stats->bgnoise, double);
	return 0;
}

static const char* log_color_to_str(LogColor color) {
	switch (color) {
		case LOG_RED:
			return "red";
		case LOG_SALMON:
			return "salmon";
		case LOG_GREEN:
			return "green";
		case LOG_BLUE:
			return "blue";
		default:
			return NULL;
	}
}

static int distodata_to_py(const disto_params *disto, unsigned char* ptr, size_t maxlen) {
	if (!disto || !ptr)
		return 1;

	unsigned char *start_ptr = ptr;

	COPY_BE64((int64_t) disto->index, int64_t);
	COPY_BE64((double) disto->velocity.x, double);
	COPY_BE64((double) disto->velocity.y, double);
	if (disto->filename)
		COPY_STRING(disto->filename);
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
		COPY_BE64(*(double*) desc->data, double);
		*type = CONFIG_TYPE_DOUBLE;
		*value = double_val;
		*value_size = sizeof(double);
	}
	else if (desc->type == STYPE_STR || desc->type == STYPE_STRDIR) {
		*type = CONFIG_TYPE_STR;
		const gchar* const_val = *(gchar**) desc->data;
		if (const_val && const_val[0] != '\0') {
			*value = g_strdup(const_val);
			*value_size = strlen((gchar*) *value) + 1;
		} else {
			*value = g_strdup("(not set)");
			*value_size = strlen(*value);
		}
	}
	else if (desc->type == STYPE_STRLIST) {
		GSList *list = *((GSList**)desc->data);
		GSList *iter = list;
		size_t total_size = 0;
		while (iter && iter->data && ((gchar*)iter->data)[0] != '\0') {
			g_strstrip((gchar*) iter->data); // Remove whitespace
			total_size += strlen((gchar*) iter->data) + 1;
			iter = iter->next;
		}
		total_size += 1; // Final null terminator

		char* list_val = g_malloc(total_size);
		char* ptr = list_val;
		iter = list;
		while (iter && iter->data && ((gchar*)iter->data)[0] != '\0') {
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

siril_plot_data* unpack_plot_data(const uint8_t* buffer, size_t buffer_size) {
	size_t offset = 0;

	// Allocate the main plot data structure
	siril_plot_data* plot_data = init_siril_plot_data();
	if (!plot_data)
		return NULL;

	// We don't need to use the siril_plot_set_X functions here as we
	// know the plot_data is newly allocated and initialized

	// Unpack title (null-terminated string)
	plot_data->title = g_strdup((const char*)buffer + offset);
	offset += strlen(plot_data->title) + 1;

	// Unpack x and y axis labels
	plot_data->xlabel = g_strdup((const char*)buffer + offset);
	offset += strlen(plot_data->xlabel) + 1;

	plot_data->ylabel = g_strdup((const char*)buffer + offset);
	offset += strlen(plot_data->ylabel) + 1;

	// Unpack savename
	plot_data->savename = g_strdup((const char*)buffer + offset);
	offset += strlen(plot_data->savename) + 1;

	// Unpack show_legend (as a single byte)
	plot_data->show_legend = BOOL_FROM_BYTE(buffer[offset]);
	offset += sizeof(uint8_t);

	// Unpack number of series (network byte-order)
	uint32_t num_series;
	memcpy(&num_series, buffer + offset, sizeof(uint32_t));
	num_series = GUINT32_FROM_BE(num_series);
	offset += sizeof(uint32_t);

	gboolean datamin_set = BOOL_FROM_BYTE(buffer[offset]);
	offset += sizeof(uint8_t);
	if (datamin_set) {
		point datamin;
		double x_BE, y_BE;
		memcpy(&x_BE, buffer + offset, sizeof(double));
		offset += sizeof(double);
		FROM_BE64_INTO(datamin.x, x_BE, double);
		memcpy(&y_BE, buffer + offset, sizeof(double));
		offset += sizeof(double);
		FROM_BE64_INTO(datamin.y, y_BE, double);
		memcpy(&plot_data->datamin, &datamin, sizeof(point));
	}

	gboolean datamax_set = BOOL_FROM_BYTE(buffer[offset]);
	offset += sizeof(uint8_t);
	if (datamax_set) {
		point datamax;
		double x_BE, y_BE;
		memcpy(&x_BE, buffer + offset, sizeof(double));
		offset += sizeof(double);
		FROM_BE64_INTO(datamax.x, x_BE, double);
		memcpy(&y_BE, buffer + offset, sizeof(double));
		offset += sizeof(double);
		FROM_BE64_INTO(datamax.y, y_BE, double);
		memcpy(&plot_data->datamax, &datamax, sizeof(point));
	}

	// Unpack series data
	for (uint32_t series_idx = 0; series_idx < num_series; series_idx++) {
		// Read series label
		char* series_label = g_strdup((const char*)buffer + offset);
		offset += strlen(series_label) + 1;

		// Unpack with_errors (as a single byte)
		// This indicates if there are errorbar series or not
		gboolean with_errors = BOOL_FROM_BYTE(buffer[offset]);
		offset += sizeof(uint8_t);

		// Read number of points (network byte-order)
		uint32_t num_points;
		memcpy(&num_points, buffer + offset, sizeof(uint32_t));
		num_points = GUINT32_FROM_BE(num_points);
		offset += sizeof(uint32_t);

		// Read plot type (network byte-order)
		uint32_t plot_type;
		memcpy(&plot_type, buffer + offset, sizeof(uint32_t));
		plot_type = GUINT32_FROM_BE(plot_type);
		offset += sizeof(uint32_t);

		// Create a new dataseries and add it to plot_data
		double *xdata = malloc(num_points * sizeof(double));
		double *ydata = malloc(num_points * sizeof(double));
		double *nerror = with_errors ? malloc(num_points * sizeof(double)) : NULL;
		double *perror = with_errors ? malloc(num_points * sizeof(double)) : NULL;
		// Read coordinates (network byte-order)
		for (uint32_t point_idx = 0; point_idx < num_points; point_idx++) {
			double x, y, x_BE, y_BE, ne, pe, ne_BE, pe_BE;

			// Read raw bytes for x
			memcpy(&x_BE, buffer + offset, sizeof(double));
			offset += sizeof(double);
			FROM_BE64_INTO(x, x_BE, double);
			xdata[point_idx] = x;

			// Read raw bytes for y
			memcpy(&y_BE, buffer + offset, sizeof(double));
			offset += sizeof(double);
			FROM_BE64_INTO(y, y_BE, double);
			ydata[point_idx] = y;

			if (with_errors) {
				// Read raw bytes for negative error
				memcpy(&ne_BE, buffer + offset, sizeof(double));
				offset += sizeof(double);
				FROM_BE64_INTO(ne, ne_BE, double);
				nerror[point_idx] = ne;

				// Read raw bytes for positive error
				memcpy(&pe_BE, buffer + offset, sizeof(double));
				offset += sizeof(double);
				FROM_BE64_INTO(pe, pe_BE, double);
				perror[point_idx] = pe;
			}

		}

		// Add to plot list (assuming simple xy plot)
		siril_plot_add_xydata(plot_data, series_label, num_points, xdata, ydata, perror, nerror);
		siril_plot_set_nth_plot_type(plot_data, series_idx+1, (enum kplottype) plot_type);
		g_free(series_label);
		free(xdata);
		free(ydata);
		free(nerror);
		free(perror);
	}

	plot_data->plottype = KPLOT_LINES;  // Default plot type

	return plot_data;
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
	int32_t payload_length = GINT32_FROM_BE(header->length);  // Convert from network byte order
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
			gboolean result = single_image_is_loaded();

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

		case CMD_GET_SELECTION: {
			gboolean result = single_image_is_loaded();

			if (result) {
				// Prepare the response data: width (gfit.rx), height (gfit.ry), and channels (gfit.naxes[2])
				uint8_t response_data[16]; // 4 x 4 bytes for x,y, w, h

				// Convert the integers to BE format for consistency across the UNIX socket
				uint32_t x_BE = GUINT32_TO_BE(com.selection.x);
				uint32_t y_BE = GUINT32_TO_BE(com.selection.y);
				uint32_t w_BE = GUINT32_TO_BE(com.selection.w);
				uint32_t h_BE = GUINT32_TO_BE(com.selection.h);

				// Copy the packed data into the response buffer
				memcpy(response_data, &x_BE, sizeof(uint32_t));
				memcpy(response_data + 4, &y_BE, sizeof(uint32_t));
				memcpy(response_data + 8, &w_BE, sizeof(uint32_t));
				memcpy(response_data + 12, &h_BE, sizeof(uint32_t));

				// Send success response with dimensions
				success = send_response(conn, STATUS_OK, response_data, sizeof(response_data));
			} else {
				// Handle error retrieving dimensions
				const char* error_msg = _("Failed to retrieve image dimensions - no image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_ACTIVE_VPORT: {
			gboolean result = single_image_is_loaded();

			if (result) {
				// Prepare the response data: width (gfit.rx), height (gfit.ry), and channels (gfit.naxes[2])
				uint8_t response_data[4]; // 4 bytes for int

				// Convert the integers to BE format for consistency across the UNIX socket
				uint32_t vport_BE = GUINT32_TO_BE(gui.cvport);

				// Copy the packed data into the response buffer
				memcpy(response_data, &vport_BE, sizeof(uint32_t));

				// Send success response with dimensions
				success = send_response(conn, STATUS_OK, response_data, sizeof(response_data));
			} else {
				// Handle error retrieving dimensions
				const char* error_msg = _("Failed to retrieve current vport - no image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_SET_SELECTION: {
			gboolean result = single_image_is_loaded();
			if (result) {
				if (payload_length == 16) {
					rectangle region_BE = *(rectangle*) payload;
					rectangle selection = {GUINT32_FROM_BE(region_BE.x),
										GUINT32_FROM_BE(region_BE.y),
										GUINT32_FROM_BE(region_BE.w),
										GUINT32_FROM_BE(region_BE.h)};
					if (selection.x < 0 || selection.x + selection.w > gfit.rx - 1 ||
								selection.y < 0 || selection.y + selection.h > gfit.ry - 1) {
						const char* error_msg = _("Failed to set selection - selection exceeds image bounds");
						success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
						if (!success)
							siril_debug_print("Error in send_response\n");
					}
					memcpy(&com.selection, &selection, sizeof(rectangle));
					if (!com.headless)
						execute_idle_and_wait_for_it(new_selection_zone, NULL);
					success = send_response(conn, STATUS_OK, NULL, 0);
				}
			} else {
				// Handle error retrieving dimensions
				const char* error_msg = _("Failed to set selection - no image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_PIXELDATA: {
			if (payload_length == 1) {
				gboolean as_preview = (gboolean) (uint8_t) payload[0];
				rectangle region = {0, 0, gfit.rx, gfit.ry};
				shared_memory_info_t *info = handle_pixeldata_request(conn, &gfit, region, as_preview);
				// Send shared memory info to Python
				success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
				free(info);
			} else {
				siril_debug_print(_("Unexpected payload length %u received for GET_PIXELDATA\n"), payload_length);
			}
			break;
		}

		case CMD_GET_STAR_IN_SELECTION: {
			gboolean error_occurred = FALSE;
			int layer = 0;
			rectangle selection;
			memcpy(&selection, &com.selection, sizeof(rectangle));
			if (payload_length == 16 || payload_length == 20) {
				rectangle region_BE = *(rectangle*) payload;
				selection = (rectangle) {GUINT32_FROM_BE(region_BE.x),
									GUINT32_FROM_BE(region_BE.y),
									GUINT32_FROM_BE(region_BE.w),
									GUINT32_FROM_BE(region_BE.h)};
			} else {
				memcpy(&selection, &com.selection, sizeof(rectangle));
			}
			if (payload_length == 20) {
				layer = GUINT32_FROM_BE(*((int*) payload + 4));
			} else {
				layer = com.headless ? 0 : match_drawing_area_widget(gui.view[select_vport(gui.cvport)].drawarea, FALSE);
			}
			if (selection.x < 0 || selection.x + selection.w > gfit.rx - 1 ||
				selection.w < 5 || selection.w > 300 ||
				selection.h < 5 || selection.h > 300 ||
				selection.y < 0 || selection.y + selection.h > gfit.ry - 1) {
					const char* error_msg = _("Invalid selection");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
			if (layer < 0 || layer >= gfit.naxes[2]) {
					const char* error_msg = _("Invalid channel");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
			}
			starprofile profile = com.pref.starfinder_conf.profile;
			psf_error error = PSF_NO_ERR;
			psf_star *psf = psf_get_minimisation(&gfit, layer, &selection, FALSE, FALSE, NULL, TRUE, profile, &error);
			if (!psf || error) {
				free_psf(psf);
				error_occurred = TRUE;
				const char* error_msg = _("Failed to find a star");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Set xpos and ypos as these are not set by minimisation
			psf->xpos = selection.x + psf->x0;
			psf->ypos = selection.y + selection.h - psf->y0;
			// Set RA and dec if there is a plate solution
			if (gfit.keywords.wcslib) {
				double fx, fy;
				display_to_siril(psf->xpos, psf->ypos, &fx, &fy, gfit.ry);
				pix2wcs2(gfit.keywords.wcslib, fx, fy, &psf->ra, &psf->dec);
				// ra and dec = -1 is the error code
			}
			const size_t psf_star_size = 36 * sizeof(double);
			unsigned char* star = g_try_malloc0(psf_star_size);
			unsigned char* ptr = star;
			if (psfstar_to_py(psf, ptr, psf_star_size)) {
				error_occurred = TRUE;
				const char* error_msg = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if(!error_occurred) {
				success = send_response(conn, STATUS_OK, star, psf_star_size);
			}
			g_free(star);
			break;
		}

		case CMD_GET_STATS_FOR_SELECTION: {
			int layer = 0;
			rectangle selection = { 0 };
			if (payload_length == 16 || payload_length == 20) {
				rectangle region_BE = *(rectangle*) payload;
				selection = (rectangle) {GUINT32_FROM_BE(region_BE.x),
									GUINT32_FROM_BE(region_BE.y),
									GUINT32_FROM_BE(region_BE.w),
									GUINT32_FROM_BE(region_BE.h)};
				if (selection.x < 0 || selection.x + selection.w > gfit.rx - 1 ||
					selection.y < 0 || selection.y  + selection.h > gfit.ry - 1) {
					const char* error_msg = _("Invalid region: selection breaches image dimensions");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
			} else {
				memcpy(&selection, &com.selection, sizeof(rectangle));
			}
			if (payload_length == 20) {
				layer = GUINT32_FROM_BE(*((int*) payload + 4));
			} else {
				layer = com.headless ? 0 : match_drawing_area_widget(gui.view[select_vport(gui.cvport)].drawarea, FALSE);
			}
			// Check an image is loaded
			if (!single_image_is_loaded()) {
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Check if channel is valid
			if (layer >= gfit.naxes[2]) {
				const char* error_msg = _("Invalid channel");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			// Check selection is valid
			if (selection.x < 0 || selection.x + selection.w > gfit.rx - 1 ||
				selection.w * selection.h < 3 ||
				selection.y < 0 || selection.y + selection.h > gfit.ry - 1) {
					const char* error_msg = _("Invalid selection");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
			}
			// Compute stats
			imstats *stats = statistics(NULL, -1, &gfit, layer, &selection, STATS_MAIN, MULTI_THREADED);

			const size_t total_size = 14 * sizeof(double);
			unsigned char* response_buffer = g_try_malloc0(total_size);
			unsigned char* ptr = response_buffer;
			if (imstats_to_py(stats, ptr, total_size)) {
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			free_stats(stats);
			g_free(response_buffer);
			break;
		}

		case CMD_GET_PIXELDATA_REGION: {
			if (payload_length == 17) {
				gboolean as_preview = (gboolean) (uint8_t) payload[0];
				unsigned char* rectptr = (unsigned char*) payload + 1;
				rectangle region_BE = *(rectangle*) rectptr;
				rectangle region = {GUINT32_FROM_BE(region_BE.x),
									GUINT32_FROM_BE(region_BE.y),
									GUINT32_FROM_BE(region_BE.w),
									GUINT32_FROM_BE(region_BE.h)};
				shared_memory_info_t *info = handle_pixeldata_request(conn, &gfit, region, as_preview);
				success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
				free(info);
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
			com.python_command = TRUE;
			int retval = processcommand(cmd, TRUE);
			com.python_command = FALSE;
			g_free(cmd);

			// Send response based on command execution
			// We always return STATUS_OK (this is handled by _request_data)
			// but return the actual retval as the payload (this is handled
			// by cmd() and results in exception raising for anything except
			// CMD_OK or CMD_NO_WAIT)
			int32_t be_retval = GINT32_TO_BE(retval);
			success = send_response(conn, STATUS_OK, &be_retval, sizeof(int32_t));
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
			if (payload_length < 3) { // Color byte plus at least one char plus null termination
				success = send_response(conn, STATUS_ERROR, "Empty log message", 0);
				break;
			}

			// Extract the color from the first byte
			LogColor color = (LogColor) payload[0];

			// Create null-terminated string from remaining payload
			char* log_msg = g_strndup(payload + 1, payload_length - 1);
			siril_log_color_message(log_msg, log_color_to_str(color));  // Assuming function name updated to handle color
			g_free(log_msg);

			// Send success response
			success = send_response(conn, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_ERROR_MESSAGEBOX: // fallthrough intentional
		case CMD_ERROR_MESSAGEBOX_MODAL:
		case CMD_WARNING_MESSAGEBOX:
		case CMD_WARNING_MESSAGEBOX_MODAL:
		case CMD_INFO_MESSAGEBOX:
		case CMD_INFO_MESSAGEBOX_MODAL: {
			// Set the title and type
			GtkMessageType type;
			const gchar *title;
			gboolean modal = TRUE;
			switch (header->command) {
				case CMD_ERROR_MESSAGEBOX:
					modal = FALSE;
				case CMD_ERROR_MESSAGEBOX_MODAL:
					type = GTK_MESSAGE_ERROR;
					title = _("Error");
					break;
				case CMD_WARNING_MESSAGEBOX:
					modal = FALSE;
				case CMD_WARNING_MESSAGEBOX_MODAL:
					type = GTK_MESSAGE_WARNING;
					title = _("Warning");
					break;
				case CMD_INFO_MESSAGEBOX:
					modal = FALSE;
				case CMD_INFO_MESSAGEBOX_MODAL:
					type = GTK_MESSAGE_INFO;
					title = _("Information");
					break;
				default:
					type = GTK_MESSAGE_OTHER;
					title = _("Unknown dialog type");
			}
			// Ensure null-terminated string for log message
			char* log_msg = g_strndup(payload, payload_length);

			// If we are headess we can't use a siril_message_dialog, so we just print the message to the log
			if (com.headless) {
				if (type == GTK_MESSAGE_INFO)
					siril_log_message(log_msg);
				else if (type == GTK_MESSAGE_WARNING)
					siril_log_color_message(log_msg, "salmon");
				else if (type == GTK_MESSAGE_ERROR)
					siril_log_color_message(log_msg, "red");
				g_free(log_msg);
				success = send_response(conn, STATUS_OK, NULL, 0);
				break;
			}

			if (!modal) {
				queue_message_dialog(type, title, log_msg);
			} else {
				struct message_data *data = malloc(sizeof(struct message_data));
				data->type = type;
				data->title = strdup(title);
				data->text = strdup(log_msg);
				siril_debug_print("Executing modal dialog\n");
				if (!com.headless)
					execute_idle_and_wait_for_it(siril_message_dialog_idle, data);
			}
			g_free(log_msg);

			// Send success response
			success = send_response(conn, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_UNDO_SAVE_STATE: {
			if (single_image_is_loaded()) {
				if (com.headless) {
					const char* error_msg = _("Undo not supported in headless mode");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				// Ensure null-terminated string for undo message
				char* log_msg = g_strndup(payload, payload_length);
				undo_save_state(&gfit, log_msg);
				g_free(log_msg);

				// Send success response
				success = send_response(conn, STATUS_OK, NULL, 0);
			} else {
				// Handle error
				const char* error_msg = _("Image not loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
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
			gchar *configdir = g_build_path(G_DIR_SEPARATOR_S, siril_get_config_dir(), "siril", NULL);
			// Ensure the config directory is available
			if (configdir && strlen(configdir) > 0) {
				// Send success response with the working directory string
				success = send_response(conn, STATUS_OK, configdir, strlen(configdir));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = _("Error: user config directory not set");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			g_free(configdir);
			break;
		}

		case CMD_GET_USERDATA_DIR: {
			const gchar *uddir = siril_get_user_data_dir();
			// Ensure the config directory is available
			if (uddir && strlen(uddir) > 0) {
				// Send success response with the working directory string
				success = send_response(conn, STATUS_OK, uddir, strlen(uddir));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = _("Error: user data directory not set");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_SYSTEMDATA_DIR: {
			const gchar *sddir = siril_get_system_data_dir();
			// Ensure the config directory is available
			if (sddir && strlen(sddir) > 0) {
				// Send success response with the working directory string
				success = send_response(conn, STATUS_OK, sddir, strlen(sddir));
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = _("Error: system data directory not set");
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

		case CMD_SEQ_FRAME_SET_PIXELDATA: {
			if (!com.seq.seqname) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if (payload_length != 4 + sizeof(incoming_image_info_t)) {
				siril_debug_print("Invalid payload length for SET_PIXELDATA: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int32_t index = GUINT32_FROM_BE((int32_t) *payload);
			if (index >= com.seq.number) {
				const char* error_msg = _("Failed to load sequence frame");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			const char* info = payload + 4;

			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame(&com.seq, index, fit, FALSE, -1)) {
				const char* error_msg = _("Failed to load sequence frame");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				clearfits(fit);
				free(fit);
				break;
			}

			int out_index = -1;
			// Compute out_index for SER / FITSEQ
			if (com.seq.type == SEQ_SER || com.seq.type == SEQ_FITSEQ) {
				for (int temp_index = 0 ; temp_index <= index; temp_index++) {
					if (com.seq.imgparam[temp_index].incl) {
						out_index++;
					}
				}
				if (out_index == -1) {
					const char* error_msg = _("Failed to compute output index");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					clearfits(fit);
					free(fit);
					break;
				}
			}

			// Update the pixel data in the sequence frame fit
			success = handle_set_pixeldata_request(conn, fit, info, payload_length - 4);
			int writer_retval;
			// Write the sequence frame
			if (com.seq.type == SEQ_SER) {
				writer_retval = ser_write_frame_from_fit(com.seq.ser_file, fit, out_index);
			} else if (com.seq.type == SEQ_FITSEQ) {
				writer_retval = fitseq_write_image(com.seq.fitseq_file, fit, out_index);
			} else {
				char *dest = fit_sequence_get_image_filename_prefixed(&com.seq,
						"", index);
				fit->bitpix = fit->orig_bitpix;
				writer_retval = savefits(dest, fit);
				free(dest);
			}
			if (fit->rx != com.seq.rx || fit->ry != com.seq.ry) {
				// Mark the sequence as variable
				com.seq.is_variable = TRUE;
				// Update the imgparam rx and ry
				com.seq.imgparam[index].rx = fit->rx;
				com.seq.imgparam[index].ry = fit->ry;
				// Clean the sequence registration data, stats and selection as they will no longer be valid
				clean_sequence(&com.seq, TRUE, TRUE, TRUE);
			}
			clearfits(fit);
			free(fit);
			if (writer_retval) {
				siril_log_color_message(_("Error writing sequence frame %i from Python\n"), "red", index);
			}
			if (!com.headless && com.seq.current == index) {
				execute_idle_and_wait_for_it(seq_load_image_in_thread, &index);
			}
			break;
		}

		case CMD_PLOT: {
			if (payload_length != sizeof(incoming_image_info_t)) {
				siril_debug_print("Invalid payload length for PLOT: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			} else {
				incoming_image_info_t* info = (incoming_image_info_t*)payload;
				info->size = GUINT64_FROM_BE(info->size);
				success = handle_plot_request(conn, info);
			}
			break;
		}

		case CMD_SET_BGSAMPLES: {
			if (payload_length != sizeof(incoming_image_info_t)) {
				siril_debug_print("Invalid payload length for SET_BGSAMPLES: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			} else {
				incoming_image_info_t* info = (incoming_image_info_t*)payload;
				info->size = GUINT64_FROM_BE(info->size);
				info->data_type = GUINT32_FROM_BE(info->data_type);
				info->channels = GUINT32_FROM_BE(info->channels);
				gboolean show_samples = (gboolean) info->data_type;
				gboolean recalculate = (gboolean) info->channels;
				success = handle_set_bgsamples_request(conn, info, show_samples, recalculate);
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
			unsigned char *response_buffer = g_try_malloc0(total_size);
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
			unsigned char *response_buffer = g_try_malloc0(total_size);
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
			if (payload_length != 8 || index >= com.seq.number || chan < 0 || chan > com.seq.nb_layers) {
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
			if (!com.seq.regparam[chan]) {
				const char* error_message = _("No regdata available");
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
				g_free(response_buffer);
				break;
			}

			unsigned char *ptr = response_buffer;
			regdata *regparam = &com.seq.regparam[chan][index];
			if (regdata_to_py(regparam, ptr, total_size)) {
				const char* error_message = _("No regdata available");
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_GET_SEQ_FRAME_FILENAME: {
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
				const char* error_msg = _("Incorrect command arguments");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			char frame_filename[256];
			char *result = NULL;
			gchar *absolute_path = NULL;
			result = seq_get_image_filename(&com.seq, index, frame_filename);
			if (com.wd && result) {
				if (com.seq.type == SEQ_REGULAR) {
					absolute_path = g_build_filename(com.wd, frame_filename, NULL);
				} else {
					absolute_path = g_strdup(frame_filename);
				}
				success = send_response(conn, STATUS_OK, absolute_path, strlen(absolute_path));
			} else {
				const char* error_msg = _("Error building frame filename");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			g_free(absolute_path);
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
			if (payload_length != 8 || index >= com.seq.number || chan < 0 || chan > com.seq.nb_layers) {
				const char* error_msg = _("Incorrect command arguments");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 14 * sizeof(double);
			unsigned char *response_buffer = g_try_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			if (!com.seq.stats) {
				g_free(response_buffer);
				const char* error_message = _("No stats for this sequence");
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
				break;
			}
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
			if (payload_length != 4 || index >= com.seq.number) {
				const char* error_msg = _("Incorrect command argument");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for response
			size_t total_size = 6 * sizeof(double);
			unsigned char *response_buffer = g_try_malloc0(total_size);
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
			printf("payload length: %d\n", payload_length);
			if (payload_length != 21 && payload_length != 5) {
				const char* error_msg = _("Incorrect payload");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			gboolean as_preview = (uint8_t)payload[0] ? TRUE : FALSE;
			const char *indexptr = payload + 1;
			int index = GUINT32_FROM_BE(*(int*) indexptr);
			rectangle region = com.seq.is_variable ?
						(rectangle) {0, 0, com.seq.imgparam[index].rx, com.seq.imgparam[index].ry} :
						(rectangle) {0, 0, com.seq.rx, com.seq.ry};
			if (payload_length == 21) {
				// Use the provided rectangle instead of the full image
				const char *rectptr = payload + 4;
				region = *(rectangle*) rectptr;
				region.x = GUINT32_FROM_BE(region.x);
				region.y = GUINT32_FROM_BE(region.y);
				region.w = GUINT32_FROM_BE(region.w);
				region.h = GUINT32_FROM_BE(region.h);
			}
			if(enforce_area_in_image(&region, &com.seq, index)) {
				siril_log_message(_("Selection cropped to frame boundaries\n"));
			}
			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame(&com.seq, index, fit, FALSE, -1)) {
				clearfits(fit);
				free(fit);
				const char* error_msg = _("Failed to read sequence frame");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			rectangle nullrect = { 0 };
			if (!memcmp(&nullrect, &region, sizeof(rectangle))) {
				// If the sequence is variable, set the size to the size of this frame
				region = (rectangle) {0, 0, fit->rx, fit->ry};
			} else {
				if (region.x < 0 || region.y < 0 || region.x + region.w > fit->rx ||
									region.y + region.h > fit->ry) {
					clearfits(fit);
					free(fit);
					const char* error_msg = _("Invalid dimensions");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
			}
			shared_memory_info_t *info = handle_pixeldata_request(conn, fit, region, as_preview);
			success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
			free(info);
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
			gboolean with_pixels = TRUE;
			gboolean as_preview = FALSE;
			int index;
			if (payload_length == 6) {
				index = GUINT32_FROM_BE(*(int*) payload);
				const char* pixelbool = payload + 4;
				const char* previewbool = payload + 5;
				with_pixels = BOOL_FROM_BYTE(*pixelbool);
				as_preview = BOOL_FROM_BYTE(*previewbool);
			}
			if (payload_length != 6 || index >= com.seq.number) {
				const char* error_msg = _("Incorrect command argument");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame(&com.seq, index, fit, FALSE, -1)) {
				free(fit);
				const char* error_msg = _("Failed to read frame metadata");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t ffit_size = sizeof(uint64_t) * 13; // 14 vars packed to 64-bit
			size_t strings_size = FLEN_VALUE * 13;  // 13 string fields of FLEN_VALUE
			size_t numeric_size = sizeof(uint64_t) * 39; // 39 vars packed to 64-bit
			size_t total_size = ffit_size + strings_size + numeric_size;
			if (with_pixels) {
				size_t shminfo_size = sizeof(shared_memory_info_t);
				total_size += shminfo_size;
			}

			unsigned char *response_buffer = g_try_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			int ret = fits_to_py(fit, ptr, ffit_size);
			if (ret) {
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				goto CLEANUP;
			}

			ptr = response_buffer + ffit_size;
			ret = keywords_to_py(fit, ptr, (strings_size + numeric_size));
			if (ret) {
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
				goto CLEANUP;
			}

			if (with_pixels) {
				ptr += strings_size + numeric_size;
				rectangle region = (rectangle) {0, 0, fit->rx, fit->ry};
				shared_memory_info_t *info = handle_pixeldata_request(conn, fit, region, as_preview);
				if (!info) {
					const char* error_message = _("Memory allocation error");
					success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
					free(info);
					goto CLEANUP;
				}
				// Convert the values to BE
				TO_BE64_INTO(info->size, info->size, size_t);
				info->data_type = GUINT32_TO_BE(info->data_type);
				info->width = GUINT32_TO_BE(info->width);
				info->height = GUINT32_TO_BE(info->height);
				info->channels = GUINT32_TO_BE(info->channels);
				memcpy(ptr, info, sizeof(shared_memory_info_t));
				free(info);
			}
			success = send_response(conn, STATUS_OK, response_buffer, total_size);

CLEANUP:
			g_free(response_buffer);
			clearfits(fit);
			free(fit);
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
			unsigned char *response_buffer = g_try_malloc0(total_size);
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
			int nb_stars = starcount(com.stars);
			int channel = 1;
			psf_star **stars = NULL;
			gboolean stars_needs_freeing;
			if (nb_stars < 1) {
				image *input_image = NULL;
				input_image = calloc(1, sizeof(image));
				input_image->fit = &gfit;
				input_image->from_seq = NULL;
				input_image->index_in_seq = -1;
				if (gfit.naxes[2] == 1)
					channel = 0;
				stars = peaker(input_image, channel, &com.pref.starfinder_conf, &nb_stars,
						NULL, FALSE, FALSE, MAX_STARS, PSF_MOFFAT_BFREE, com.max_thread);
				free(input_image);
				stars_needs_freeing = TRUE;
			} else {
				stars = com.stars;
				stars_needs_freeing = FALSE;
			}

			// Count the number of stars in com.stars
			nb_stars = starcount(stars);
			if (nb_stars < 1) {
				const char* error_msg = _("No stars in image");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}

			const size_t psf_star_size = 36 * sizeof(double);
			const size_t total_size = nb_stars * psf_star_size;

			unsigned char* allstars = g_try_malloc0(total_size);
			if (!allstars) {
				const char* error_msg = _("Memory allocation failed");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			gboolean error_occurred = FALSE;
			for (int i = 0; i < nb_stars; i++) {
				if(!stars) {
					const char* error_msg = _("Stars array was cleared mid-process");
					error_occurred = TRUE;
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				psf_star *psf = stars[i];
				if (!psf) {
					error_occurred = TRUE;
					const char* error_msg = _("Unexpected null star entry");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					if (stars_needs_freeing) {
						free_psf_starstarstar(stars);
					}
					break;
				}

				// Calculate the correct offset for this star's data
				unsigned char* ptr = allstars + (i * psf_star_size);

				if (psfstar_to_py(psf, ptr, psf_star_size)) {
					error_occurred = TRUE;
					const char* error_msg = _("Memory allocation error");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					if (stars_needs_freeing) {
						free_psf_starstarstar(stars);
					}
					break;
				}
			}

			if (stars_needs_freeing) {
				free_psf_starstarstar(stars);
			}

			if (!error_occurred) {
				shared_memory_info_t *info = handle_rawdata_request(conn, allstars, total_size);
				success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
				free(info);
			}

			g_free(allstars);
			break;
		}

		case CMD_GET_BGSAMPLES: {
			int nb_samples = 0;
			if (com.grad_samples) {
				// Count the number of stars in com.grad_samples
				nb_samples = g_slist_length(com.grad_samples);
			}
			if (!nb_samples) {
				const char* error_msg = _("No background samples list available");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}

			const size_t sample_size = 8 * sizeof(double) + 2 * sizeof(uint64_t);
			const size_t total_size = nb_samples * sample_size;

			unsigned char* allsamples = g_try_malloc0(total_size);
			if (!allsamples) {
				const char* error_msg = _("Memory allocation failed");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			gboolean error_occurred = FALSE;
			int i = 0;
			for (GSList *iter = com.grad_samples; iter ; iter = iter->next) {
				if (i >= nb_samples) {
					const char* error_msg = _("Mismatch in samples count");
					error_occurred = TRUE;
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				background_sample *sample = (background_sample*) iter->data;
				if (!sample) {
					error_occurred = TRUE;
					const char* error_msg = _("Unexpected null background sample");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}

				// Calculate the correct offset for this star's data
				unsigned char* ptr = allsamples + (i * sample_size);

				if (sample_to_py(sample, ptr, sample_size)) {
					error_occurred = TRUE;
					const char* error_msg = _("Memory allocation error");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				i++;
			}

			if (!error_occurred) {
				shared_memory_info_t *info = handle_rawdata_request(conn, allsamples, total_size);
				success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
				free(info);
			}

			g_free(allsamples);
			break;
		}

		case CMD_GET_IMAGE: {
			if (!single_image_is_loaded()) {
				const char* error_msg = _("No image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}

			// Calculate size needed for the response
			size_t total_size = sizeof(uint64_t) * 13; // 14 vars packed to 64-bit

			unsigned char *response_buffer = g_try_malloc0(total_size);
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

			shared_memory_info_t *info = handle_rawdata_request(conn, profile_data, profile_size);
			success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
			free(info);
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
			shared_memory_info_t *info = handle_rawdata_request(conn, fit->header, length);
			success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
			free(info);
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
			shared_memory_info_t *info = handle_rawdata_request(conn, buffer, total_length * sizeof(char));
			success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
			free(info);
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

			shared_memory_info_t *info = handle_rawdata_request(conn, fit->unknown_keys, length * sizeof(char));
			success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
			free(info);
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
			unsigned char* response_buffer = g_try_malloc0(total_size);

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

		case CMD_PIX2WCS: {
			gboolean result = single_image_is_loaded() || sequence_is_loaded();
			if (result) {
				if (!has_wcs(&gfit)) {
					// Handle no WCS error
					const char* error_msg = _("Siril image is not plate solved");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				if (payload_length == 16) {
					double *DblPtrBE = (double*) payload;
					double x_BE = DblPtrBE[0];
					double y_BE = DblPtrBE[1];
					double x, y;
					FROM_BE64_INTO(x, x_BE, double);
					FROM_BE64_INTO(y, y_BE, double);
					double ra, dec, ra_BE, dec_BE;
					double fx, fy;
					display_to_siril(x, y, &fx, &fy, gfit.ry);
					pix2wcs2(gfit.keywords.wcslib, fx, fy, &ra, &dec);
					// ra and dec = -1 is the error code
					TO_BE64_INTO(ra_BE, ra, double);
					TO_BE64_INTO(dec_BE, dec, double);
					unsigned char* payload = g_try_malloc0(2 * sizeof(double));
					DblPtrBE = (double*) payload;
					DblPtrBE[0] = ra_BE;
					DblPtrBE[1] = dec_BE;
					success = send_response(conn, STATUS_OK, payload, 2 * sizeof(double));
					g_free(payload);
					break;
				}
			} else {
				// Handle error retrieving dimensions
				const char* error_msg = _("Failed to set selection - no image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_WCS2PIX: {
			gboolean result = single_image_is_loaded() || sequence_is_loaded();
			if (result) {
				if (!has_wcs(&gfit)) {
					// Handle no WCS error
					const char* error_msg = _("Siril image is not plate solved");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				if (payload_length == 16) {
					double *DblPtrBE = (double*) payload;
					double ra_BE = DblPtrBE[0];
					double dec_BE = DblPtrBE[1];
					double ra, dec;
					FROM_BE64_INTO(ra, ra_BE, double);
					FROM_BE64_INTO(dec, dec_BE, double);
					double x, y, fx, fy, x_BE, y_BE;
					wcs2pix(&gfit, ra, dec, &fx, &fy);
					siril_to_display(fx, fy, &x, &y, gfit.ry);
					TO_BE64_INTO(x_BE, x, double);
					TO_BE64_INTO(y_BE, y, double);
					unsigned char* payload = g_try_malloc0(2 * sizeof(double));
					DblPtrBE = (double*) payload;
					DblPtrBE[0] = x_BE;
					DblPtrBE[1] = y_BE;
					success = send_response(conn, STATUS_OK, payload, 2 * sizeof(double));
					g_free(payload);
					break;
				}
			} else {
				// Handle error retrieving dimensions
				const char* error_msg = _("Failed to set selection - no image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_BUNDLE_PATH: {
			// Return the value of get_siril_bundle_path()
#ifdef _WIN32
			gchar *path = get_siril_bundle_path();
			if (path) {
				// Send success response with the working directory string
				success = send_response(conn, STATUS_OK, path, strlen(path));
				g_free(path);
			} else {
				// Handle error retrieving the working directory
				const char* error_msg = _("Failed to retrieve bundle path");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
#else
			// Handle error: OS is not Windows
			const char* error_msg = _("_get_bundle_path() only applicable on Windows");
			success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
#endif
			break;
		}

		case CMD_CLAIM_THREAD: {
			int ret = claim_thread_for_python();
			if (ret == 1) {
				// Unable to claim the thread
				const char* error_msg = _("Thread is busy");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
			} else if (ret == 2) {
				// Unable to claim the thread
				const char* error_msg = _("Image processing dialog is open");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
			} else  if (ret == 0) {
				// Thread claimed, we can safely do gfit processing tasks
				conn->thread_claimed = TRUE;
				success = send_response(conn, STATUS_OK, NULL, 0);
			} else {
				// Unable to claim the thread
				const char* error_msg = _("Unknown error");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_RELEASE_THREAD: {
			if (conn->thread_claimed) {
				python_releases_thread();
				conn->thread_claimed = FALSE;
			}
			// Otherwise it's not ours to release. We fail silently here, because
			// release_thread() should be called in a finally: block
			success = send_response(conn, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_REQUEST_SHM: {
			uint64_t total_size;
			if (payload_length == 8) {
				total_size = GUINT64_FROM_BE(*(uint64_t*) payload);
				shared_memory_info_t *info = handle_rawdata_request(conn, NULL, total_size);
				success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
				free(info);

			} else {
				const char* error_msg = _("Incorrect payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_SET_SEQ_FRAME_INCL: {
			uint32_t index;
			gboolean incl;
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if (payload_length == 8) {
				// Extract the first 4 bytes as `index`
				memcpy(&index, payload, sizeof(uint32_t));
				index = GUINT32_FROM_BE(index);

				// Extract the second 4 bytes as `incl`
				uint32_t incl_encoded;
				memcpy(&incl_encoded, payload + sizeof(uint32_t), sizeof(uint32_t));
				incl_encoded = GUINT32_FROM_BE(incl_encoded);
				incl = (gboolean) incl_encoded;

				if (index >= com.seq.number) {
					const char* error_msg = _("Index is out of range");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
					break;
				}
				com.seq.imgparam[index].incl = incl;
				if (index == com.seq.current && !com.headless)
					redraw(REDRAW_OVERLAY);
				success = send_response(conn, STATUS_OK, NULL, 0);
			} else {
				const char* error_msg = _("Incorrect payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_SEQ_DISTODATA: {
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int chan;
			if (payload_length == 4) {
				chan = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || chan < 0 || chan > com.seq.nb_layers) {
				const char* error_msg = _("Incorrect command arguments");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if (com.seq.distoparam[chan].index == DISTO_UNDEF) {
				const char* error_message = _("No distodata available");
				success = send_response(conn, STATUS_NONE, error_message, strlen(error_message));
				break;
			}
			// Calculate size needed for the response
			size_t stringsize = 0;
			if (com.seq.distoparam[chan].filename)
				stringsize = strlen(com.seq.distoparam[chan].filename) + 1;
			size_t varsize = sizeof(int64_t) + 2 * sizeof(double);
			size_t total_size = varsize + stringsize;
			unsigned char *response_buffer = g_malloc0(total_size);
			unsigned char *ptr = response_buffer;

			if (distodata_to_py(&com.seq.distoparam[chan], ptr, total_size)) {
				const char* error_message = _("Memory allocation error");
				success = send_response(conn, STATUS_ERROR, error_message, strlen(error_message));
			} else {
				success = send_response(conn, STATUS_OK, response_buffer, total_size);
			}
			g_free(response_buffer);
			break;
		}

		case CMD_SET_IMAGE_HEADER: {
			if (payload_length != sizeof(incoming_image_info_t)) {
				siril_debug_print("Invalid payload length for SET_IMAGE_HEADER: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			} else {
				incoming_image_info_t* info = (incoming_image_info_t*)payload;
				info->size = GUINT64_FROM_BE(info->size);
				success = handle_set_image_header_request(conn, info);
			}
			break;
		}

		case CMD_ADD_USER_POLYGON: {
			UserPolygon *polygon = deserialize_polygon((const uint8_t*) payload, payload_length);
			if (!(single_image_is_loaded() || sequence_is_loaded())) {
				siril_debug_print("Failed to add user polygon\n");
				const char* error_msg = _("Failed to add user polygon: no image loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if (polygon) {
				int id = get_unused_polygon_id();
				polygon->id = id;
				gui.user_polygons = g_list_append(gui.user_polygons, polygon);
				redraw(REDRAW_OVERLAY);
				int id_be = GINT32_TO_BE(id);
				success = send_response(conn, STATUS_OK, &id_be, 4);
			} else {
				siril_debug_print("Failed to add user polygon\n");
				const char* error_msg = _("Failed to add user polygon");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_DELETE_USER_POLYGON: {
			if (payload_length == 4) {
				int32_t id = GINT32_FROM_BE(*(int*) payload);
				gboolean deleted = delete_user_polygon(id);
				if (!deleted) {
					siril_debug_print("Failed to delete user polygon with id %d\n", id);
					const char* error_msg = _("Invalid payload length");
					success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				}
				success = send_response(conn, STATUS_OK, NULL, 0);
			} else {
				siril_debug_print("Invalid payload length for DELETE_USER_POLYGON: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}
		case CMD_CLEAR_USER_POLYGONS: {
			clear_user_polygons();
			success = send_response(conn, STATUS_OK, NULL, 0);
			break;
		}

		case CMD_GET_USER_POLYGON: {
			if (payload_length == 4) {
				int32_t id = GINT32_FROM_BE(*(int*) payload);
				UserPolygon *polygon = find_polygon_by_id(id);
				if (!polygon) {
					siril_debug_print("Failed to find a user polygon with id %d\n", id);
					const char* error_msg = _("No polygon found matching id");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				}
				size_t polygon_size;
				uint8_t *serialized = serialize_polygon(polygon, &polygon_size);
				if (!serialized) {
					siril_debug_print("Failed to serialize the user polygon with id %d\n", id);
					const char* error_msg = _("Failed to serialize user polygon");
					success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				}
				shared_memory_info_t *info = handle_rawdata_request(conn, serialized, polygon_size);
				success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
				g_free(serialized);
				free(info);
			} else {
				siril_debug_print("Invalid payload length for GET_USER_POLYGON: %u\n", payload_length);
				const char* error_msg = _("Invalid payload length");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
			break;
		}

		case CMD_GET_USER_POLYGON_LIST: {
			size_t polygon_list_size;
			if (g_list_length(gui.user_polygons) == 0) {
				siril_debug_print("No user polygons defined\n");
				const char* error_msg = _("No user polygons to serialize");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
			}

			uint8_t *serialized = serialize_polygon_list(gui.user_polygons, &polygon_list_size);
			if (!serialized) {
				siril_debug_print("Failed to serialize the user polygon list\n");
				const char* error_msg = _("Failed to serialize user polygon list");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
			}

			shared_memory_info_t *info = handle_rawdata_request(conn, serialized, polygon_list_size);
			success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
			g_free(serialized);
			free(info);
			break;
		}

		case CMD_CONFIRM_MESSAGEBOX: {
			const char *title = (const char *)payload;
			if (!title || title[0] == '\0' || payload_length < strlen(title) + 2) {
				const char* error_msg = _("Argument error");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			const char *message = title + strlen(title) + 1;  // Move past the null terminator
			if (!message || message[0] == '\0' || payload_length < strlen(title) + strlen(message) + 3) {
				const char* error_msg = _("Argument error");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			const char *confirm_label = message + strlen(message) + 1; // Move past the next null terminator
			if (!confirm_label || confirm_label[0] == '\0' || payload_length < strlen(title) + strlen(message) + strlen(confirm_label) + 3) {
				const char* error_msg = _("Argument error");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				break;
			}
			uint8_t retval = siril_confirm_dialog((gchar*) title, (gchar*) message, (gchar*) confirm_label) ? 1 : 0;
			success = send_response(conn, STATUS_OK, (const char*)&retval, sizeof(uint8_t));
			break;
		}

		case CMD_GET_SEQ_FRAME_HEADER: {
			if (!sequence_is_loaded()) {
				const char* error_msg = _("No sequence loaded");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			int index;
			if (payload_length == 4) {
				index = GUINT32_FROM_BE(*(int*) payload);
			}
			if (payload_length != 4 || index >= com.seq.number) {
				const char* error_msg = _("Incorrect command argument");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			fits *fit = calloc(1, sizeof(fits));
			if (seq_read_frame_metadata(&com.seq, index, fit)) {
				clearfits(fit);
				free(fit);
				const char* error_msg = _("Failed to read frame metadata");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
				break;
			}
			if (fit->header == NULL) {
				const char* error_msg = _("Image has no FITS header");
				success = send_response(conn, STATUS_NONE, error_msg, strlen(error_msg));
				clearfits(fit);
				free(fit);
				break;
			}
			// Prepare data
			guint32 length = strlen(fit->header) + 1;
			shared_memory_info_t *info = handle_rawdata_request(conn, fit->header, length);
			success = send_response(conn, STATUS_OK, (const char*)info, sizeof(*info));
			free(info);
			clearfits(fit);
			free(fit);
			break;
		}

		case CMD_CREATE_NEW_SEQ: {
			const gchar* seqname = g_strndup(payload, payload_length);
			if (create_one_seq(seqname, SEQ_REGULAR)) {
				success = send_response(conn, STATUS_OK, NULL, 0);
			} else {
				const char* error_msg = _("Could not create the new sequence");
				success = send_response(conn, STATUS_ERROR, error_msg, strlen(error_msg));
			}
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
