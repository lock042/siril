/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/command.h"
#include "algos/demosaicing.h"
#include "algos/geometry.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "extraction.h"

/******************************************************************************
 * Note for maintainers: do not use the translation macro on the following    *
 * string. Color management relies on being able to detect "Extraction"       *
 * in FITS HISTORY header.                                                    *
 ******************************************************************************/
const gchar *extract_string = "Extraction";

static int cfa_extract_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer);

static void update_sampling_information(fits *fit, float factor) {
	clear_Bayer_information(fit);

	fit->keywords.pixel_size_x *= factor;
	fit->keywords.pixel_size_y *= factor;
}

void update_filter_information(fits *fit, char *filter, gboolean append) {
	gchar *filtername = NULL;
	if (append) {
		gchar *currfilter = g_strstrip(fit->keywords.filter);
		if (g_strcmp0(currfilter, filter) && !(currfilter[0] == '\0')) { // to deal with case of Ha filter (avoids Ha_Ha) or empty filter name
			filtername = g_strconcat(currfilter, "_", filter, NULL);
		} else filtername = g_strdup(filter);
	} else {
		filtername = g_strdup(filter);
	}
	strncpy(fit->keywords.filter, filtername, FLEN_VALUE - 1);
	g_free(filtername);
}

int extractHa_ushort(fits *in, fits *Ha, sensor_pattern pattern, extraction_scaling scaling) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Extract_Ha does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_USHORT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			WORD c0 = in->data[col + row * in->rx];
			WORD c1 = in->data[1 + col + row * in->rx];
			WORD c2 = in->data[col + (1 + row) * in->rx];
			WORD c3 = in->data[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c0) : c0;
				break;
			case BAYER_FILTER_BGGR:
				Ha->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c3) : c3;
				break;
			case BAYER_FILTER_GRBG:
				Ha->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c1) : c1;
				break;
			case BAYER_FILTER_GBRG:
				Ha->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c2) : c2;
				break;
			default:
				clearfits(Ha);
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}
	if (scaling == SCALING_HA_UP) {
		verbose_resize_gaussian(Ha, Ha->rx * 2, Ha->ry * 2, OPENCV_LANCZOS4, TRUE);
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha, scaling ? 1.f : 2.f);
	update_filter_information(Ha, "Ha", TRUE);
	Ha->history = g_slist_append(Ha->history, g_strdup(_("Ha channel")));

	return 0;
}

int extractHa_float(fits *in, fits *Ha, sensor_pattern pattern, extraction_scaling scaling) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Extract_Ha does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_FLOAT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			float c0 = in->fdata[col + row * in->rx];
			float c1 = in->fdata[1 + col + row * in->rx];
			float c2 = in->fdata[col + (1 + row) * in->rx];
			float c3 = in->fdata[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->fdata[j] = c0;
				break;
			case BAYER_FILTER_BGGR:
				Ha->fdata[j] = c3;
				break;
			case BAYER_FILTER_GRBG:
				Ha->fdata[j] = c1;
				break;
			case BAYER_FILTER_GBRG:
				Ha->fdata[j] = c2;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	if (scaling == SCALING_HA_UP) {
		verbose_resize_gaussian(Ha, Ha->rx * 2, Ha->ry * 2, OPENCV_LANCZOS4, TRUE);
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha, scaling ? 1.f : 2.f);
	update_filter_information(Ha, "Ha", TRUE);
	Ha->history = g_slist_append(Ha->history, g_strdup(_("Ha channel")));

	return 0;
}

sensor_pattern get_bayer_pattern(fits *fit) {
	/* Get Bayer informations from header if available */
	sensor_pattern tmp_pattern = com.pref.debayer.bayer_pattern;
	if (com.pref.debayer.use_bayer_header) {
		sensor_pattern bayer = get_cfa_pattern_index_from_string(fit->keywords.bayer_pattern);
		if (bayer <= BAYER_FILTER_MAX) {
			if (bayer != tmp_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "salmon");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"salmon", filter_pattern[bayer], filter_pattern[com.pref.debayer.bayer_pattern]);
					tmp_pattern = bayer;
				}
			}
		} else {
			siril_log_message(_("XTRANS pattern not supported for this feature.\n"));
			return 1;
		}
	}
	if (tmp_pattern >= BAYER_FILTER_MIN && tmp_pattern <= BAYER_FILTER_MAX) {
		siril_log_message(_("Filter Pattern: %s\n"),
				filter_pattern[tmp_pattern]);
	}

	adjust_Bayer_pattern(fit, &tmp_pattern);
	return tmp_pattern;
}

int extract_finalize_hook(struct generic_seq_args *args) {
	struct split_cfa_data *data = (struct split_cfa_data *) args->user;
	int retval = seq_finalize_hook(args);
	free(data);
	return retval;
}

int split_finalize_hook(struct generic_seq_args *args) {
	struct split_cfa_data *data = (struct split_cfa_data *) args->user;
	free(data);
	return 0;
}

int extractHa_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct split_cfa_data *data = (struct split_cfa_data *) args->user;
	int ret = 1;
	fits f_Ha = { 0 };
	sensor_pattern pattern = get_bayer_pattern(fit);

	if (fit->type == DATA_USHORT)
		ret = extractHa_ushort(fit, &f_Ha, pattern, data->scaling);
	else if (fit->type == DATA_FLOAT)
		ret = extractHa_float(fit, &f_Ha, pattern, data->scaling);
	else return 1;
	if (!ret) {
		clearfits(fit);
		memcpy(fit, &f_Ha, sizeof(fits));
	}
	return ret;
}

static int extract_prepare_hook(struct generic_seq_args *args) {
	int retval = seq_prepare_hook(args);
	if (!retval && args->new_ser) {
		retval = ser_reset_to_monochrome(args->new_ser);
	}
	return retval;
}

void apply_extractHa_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->seq = split_cfa_args->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->compute_mem_limits_hook = cfa_extract_compute_mem_limits;
	args->prepare_hook = extract_prepare_hook;
	args->finalize_hook = extract_finalize_hook;
	args->image_hook = extractHa_image_hook;
	args->description = _("Extract Ha");
	args->has_output = TRUE;
	args->new_seq_prefix = split_cfa_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = split_cfa_args;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

int extractGreen_ushort(fits *in, fits *green, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Extract_Green does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&green, width, height, 1, DATA_USHORT))
		return 1;

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			WORD c0, c1, c2, c3;
			switch(pattern) {
			case BAYER_FILTER_RGGB:
			case BAYER_FILTER_BGGR:
				c1 = in->data[1 + col + row * in->rx];
				c2 = in->data[col + (1 + row) * in->rx];
				green->data[j] = (c1 + c2) / 2;
				break;
			case BAYER_FILTER_GRBG:
			case BAYER_FILTER_GBRG:
				c0 = in->data[col + row * in->rx];
				c3 = in->data[1 + col + (1 + row) * in->rx];
				green->data[j] = (c0 + c3) / 2;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, green);
	update_sampling_information(green, 2.f);
	update_filter_information(green, "G", TRUE);
	green->history = g_slist_append(green->history, g_strdup_printf(_("%s: extract Green channel\n"), extract_string));

	return 0;
}

int extractGreen_float(fits *in, fits *green, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Extract_Green does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&green, width, height, 1, DATA_FLOAT))
		return 1;

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			float c0, c1, c2, c3;
			switch(pattern) {
			case BAYER_FILTER_RGGB:
			case BAYER_FILTER_BGGR:
				c1 = in->fdata[1 + col + row * in->rx];
				c2 = in->fdata[col + (1 + row) * in->rx];
				green->fdata[j] = (c1 + c2) * 0.5f;
				break;
			case BAYER_FILTER_GRBG:
			case BAYER_FILTER_GBRG:
				c0 = in->fdata[col + row * in->rx];
				c3 = in->fdata[1 + col + (1 + row) * in->rx];
				green->fdata[j] = (c0 + c3) * 0.5f;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, green);
	update_sampling_information(green, 2.f);
	update_filter_information(green, "G", TRUE);
	green->history = g_slist_append(green->history, g_strdup_printf(_("%s: extract Green channel\n"), extract_string));

	return 0;
}

int extractGreen_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 1;
	fits f_Ha = { 0 };
	sensor_pattern pattern = get_bayer_pattern(fit);

	if (fit->type == DATA_USHORT)
		ret = extractGreen_ushort(fit, &f_Ha, pattern);
	else if (fit->type == DATA_FLOAT)
		ret = extractGreen_float(fit, &f_Ha, pattern);
	else return 1;
	if (!ret) {
		clearfits(fit);
		memcpy(fit, &f_Ha, sizeof(fits));
	}
	return ret;
}

void apply_extractGreen_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->seq = split_cfa_args->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->compute_mem_limits_hook = cfa_extract_compute_mem_limits;
	args->prepare_hook = extract_prepare_hook;
	args->finalize_hook = extract_finalize_hook;
	args->image_hook = extractGreen_image_hook;
	args->description = _("Extract Green");
	args->has_output = TRUE;
	args->new_seq_prefix = split_cfa_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = split_cfa_args;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

struct _double_split {
	int index;
	fits *ha;
	fits *oiii;
};

int extractHaOIII_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 1;
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;
	sensor_pattern pattern = get_bayer_pattern(fit);

	/* Demosaic and store images for write */
	struct _double_split *double_data = malloc(sizeof(struct _double_split));
	double_data->ha = calloc(1, sizeof(fits));
	double_data->oiii = calloc(1, sizeof(fits));
	double_data->index = o;

	if (fit->type == DATA_USHORT) {
		ret = extractHaOIII_ushort(fit, double_data->ha, double_data->oiii, pattern, cfa_args->scaling, threads);
	}
	else if (fit->type == DATA_FLOAT) {
		ret = extractHaOIII_float(fit, double_data->ha, double_data->oiii, pattern, cfa_args->scaling, threads);
	}

	if (ret) {
		clearfits(double_data->ha);
		free(double_data->ha);
		clearfits(double_data->oiii);
		free(double_data->oiii);
		free(double_data);
	} else {
#ifdef _OPENMP
		omp_set_lock(&args->lock);
#endif
		cfa_args->processed_images = g_list_append(cfa_args->processed_images, double_data);
#ifdef _OPENMP
		omp_unset_lock(&args->lock);
#endif
		siril_debug_print("Ha-OIII: processed images added to the save list (%d)\n", o);
	}
	return ret;
}

static int dual_prepare(struct generic_seq_args *args) {
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;
	// we call the generic prepare twice with different prefixes
	args->new_seq_prefix = "Ha_"; // This is OK here, it does not get freed in the
	// end_generic_sequence later
	if (extract_prepare_hook(args))
		return 1;
	// but we copy the result between each call
	cfa_args->new_ser_ha = args->new_ser;
	cfa_args->new_fitseq_ha = args->new_fitseq;

	args->new_seq_prefix = "OIII_"; // This is OK here, it does not get freed in the
	// end_generic_sequence later
	if (extract_prepare_hook(args))
		return 1;
	cfa_args->new_ser_oiii = args->new_ser;
	cfa_args->new_fitseq_oiii = args->new_fitseq;

	args->new_seq_prefix = NULL;
	args->new_ser = NULL;
	args->new_fitseq = NULL;

	seqwriter_set_number_of_outputs(2);
	return 0;
}

static int dual_finalize(struct generic_seq_args *args) {
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;
	args->new_ser = cfa_args->new_ser_ha;
	args->new_fitseq = cfa_args->new_fitseq_ha;
	int retval = seq_finalize_hook(args);
	cfa_args->new_ser_ha = NULL;
	cfa_args->new_fitseq_ha = NULL;

	args->new_ser = cfa_args->new_ser_oiii;
	args->new_fitseq = cfa_args->new_fitseq_oiii;
	retval = seq_finalize_hook(args) || retval;
	cfa_args->new_ser_oiii = NULL;
	cfa_args->new_fitseq_oiii = NULL;
	seqwriter_set_number_of_outputs(1);
	free(cfa_args);
	return retval;
}

static int dual_save(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;
	struct _double_split *double_data = NULL;
	// images are passed from the image_hook to the save in a list, because
	// there are two, which is unsupported by the generic arguments
#ifdef _OPENMP
	omp_set_lock(&args->lock);
#endif
	GList *list = cfa_args->processed_images;
	while (list) {
		if (((struct _double_split *)list->data)->index == out_index) {
			double_data = list->data;
			break;
		}
		list = g_list_next(cfa_args->processed_images);
	}
	if (double_data)
		cfa_args->processed_images = g_list_remove(cfa_args->processed_images, double_data);
#ifdef _OPENMP
	omp_unset_lock(&args->lock);
#endif
	if (!double_data) {
		siril_log_color_message(_("Image %d not found for writing\n"), "red", in_index);
		return 1;
	}

	siril_debug_print("Ha-OIII: images to be saved (%d)\n", out_index);
	if (double_data->ha->naxes[0] == 0 || double_data->oiii->naxes[0] == 0) {
		siril_debug_print("empty data\n");
		return 1;
	}

	int retval1, retval2;
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		retval1 = ser_write_frame_from_fit(cfa_args->new_ser_ha, double_data->ha, out_index);
		retval2 = ser_write_frame_from_fit(cfa_args->new_ser_oiii, double_data->oiii, out_index);
		// the two fits are freed by the writing thread
		if (!retval1 && !retval2) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		retval1 = fitseq_write_image(cfa_args->new_fitseq_ha, double_data->ha, out_index);
		retval2 = fitseq_write_image(cfa_args->new_fitseq_oiii, double_data->oiii, out_index);
		// the two fits are freed by the writing thread
		if (!retval1 && !retval2) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else {
		char *dest = fit_sequence_get_image_filename_prefixed(args->seq, "Ha_", in_index);
		if (fit->type == DATA_USHORT) {
			retval1 = save1fits16(dest, double_data->ha, RLAYER);
		} else {
			retval1 = save1fits32(dest, double_data->ha, RLAYER);
		}
		free(dest);
		dest = fit_sequence_get_image_filename_prefixed(args->seq, "OIII_", in_index);
		if (fit->type == DATA_USHORT) {
			retval2 = save1fits16(dest, double_data->oiii, RLAYER);
		} else {
			retval2 = save1fits32(dest, double_data->oiii, RLAYER);
		}
		free(dest);
	}
	if (!cfa_args->new_fitseq_ha && !cfa_args->new_ser_ha) { // detect if there is a seqwriter
		clearfits(double_data->ha);
		free(double_data->ha);
		clearfits(double_data->oiii);
		free(double_data->oiii);
	}
	free(double_data);
	return retval1 || retval2;
}

void apply_extractHaOIII_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->seq = split_cfa_args->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->compute_mem_limits_hook = cfa_extract_compute_mem_limits;
	args->prepare_hook = dual_prepare;
	args->finalize_hook = dual_finalize;
	args->save_hook = dual_save;
	args->image_hook = extractHaOIII_image_hook;
	args->description = _("Extract Ha and OIII");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = 1.23;	// sqrt(1.5), for memory management
	args->new_seq_prefix = NULL;
	args->user = split_cfa_args;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

int split_cfa_ushort(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Split CFA does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&cfa0, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa1, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa2, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa3, width, height, 1, DATA_USHORT)) {
		return 1;
	}

	int j = 0;
	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			/* not c0, c1, c2 and c3 because of the read orientation */
			WORD c1 = in->data[col + row * in->rx];
			WORD c3 = in->data[1 + col + row * in->rx];
			WORD c0 = in->data[col + (1 + row) * in->rx];
			WORD c2 = in->data[1 + col + (1 + row) * in->rx];

			cfa0->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c0) : c0;
			cfa1->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c1) : c1;
			cfa2->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c2) : c2;
			cfa3->data[j] = (in->bitpix == 8) ? truncate_to_BYTE(c3) : c3;
			j++;
		}
	}

	copy_fits_metadata(in, cfa0);
	copy_fits_metadata(in, cfa1);
	copy_fits_metadata(in, cfa2);
	copy_fits_metadata(in, cfa3);

	/* we remove Bayer header because not needed now */
	clear_Bayer_information(cfa0);
	clear_Bayer_information(cfa1);
	clear_Bayer_information(cfa2);
	clear_Bayer_information(cfa3);

	/* Update history */
	cfa0->history = g_slist_append(cfa0->history, g_strdup(_("CFA0 Bayer channel\n")));
	cfa1->history = g_slist_append(cfa1->history, g_strdup(_("CFA1 Bayer channel\n")));
	cfa2->history = g_slist_append(cfa2->history, g_strdup(_("CFA2 Bayer channel\n")));
	cfa3->history = g_slist_append(cfa3->history, g_strdup(_("CFA3 Bayer channel\n")));

	return 0;
}

int split_cfa_float(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Split CFA does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&cfa0, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa1, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa2, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa3, width, height, 1, DATA_FLOAT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			/* not c0, c1, c2 and c3 because of the read orientation */
			float c1 = in->fdata[col + row * in->rx];
			float c3 = in->fdata[1 + col + row * in->rx];
			float c0 = in->fdata[col + (1 + row) * in->rx];
			float c2 = in->fdata[1 + col + (1 + row) * in->rx];

			cfa0->fdata[j] = c0;
			cfa1->fdata[j] = c1;
			cfa2->fdata[j] = c2;
			cfa3->fdata[j] = c3;
			j++;
		}
	}

	copy_fits_metadata(in, cfa0);
	copy_fits_metadata(in, cfa1);
	copy_fits_metadata(in, cfa2);
	copy_fits_metadata(in, cfa3);

	/* we remove Bayer header because not needed now */
	clear_Bayer_information(cfa0);
	clear_Bayer_information(cfa1);
	clear_Bayer_information(cfa2);
	clear_Bayer_information(cfa3);

	/* Update history */
	cfa0->history = g_slist_append(cfa0->history, g_strdup(_("CFA0 Bayer channel\n")));
	cfa1->history = g_slist_append(cfa1->history, g_strdup(_("CFA1 Bayer channel\n")));
	cfa2->history = g_slist_append(cfa2->history, g_strdup(_("CFA2 Bayer channel\n")));
	cfa3->history = g_slist_append(cfa3->history, g_strdup(_("CFA3 Bayer channel\n")));

	return 0;
}

int split_cfa_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 1;
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;

	fits f_cfa0 = { 0 }, f_cfa1 = { 0 }, f_cfa2 = { 0 }, f_cfa3 = { 0 };

	gchar *cfa0 = g_strdup_printf("%s0_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);
	gchar *cfa1 = g_strdup_printf("%s1_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);
	gchar *cfa2 = g_strdup_printf("%s2_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);
	gchar *cfa3 = g_strdup_printf("%s3_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);

	if (fit->type == DATA_USHORT) {
		if (!(ret = split_cfa_ushort(fit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits16(cfa0, &f_cfa0, RLAYER) ||
				save1fits16(cfa1, &f_cfa1, RLAYER) ||
				save1fits16(cfa2, &f_cfa2, RLAYER) ||
				save1fits16(cfa3, &f_cfa3, RLAYER);
		}
	}
	else if (fit->type == DATA_FLOAT) {
		if (!(ret = split_cfa_float(fit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits32(cfa0, &f_cfa0, RLAYER) ||
				save1fits32(cfa1, &f_cfa1, RLAYER) ||
				save1fits32(cfa2, &f_cfa2, RLAYER) ||
				save1fits32(cfa3, &f_cfa3, RLAYER);
		}
	}

	g_free(cfa0); g_free(cfa1);
	g_free(cfa2); g_free(cfa3);
	clearfits(&f_cfa0); clearfits(&f_cfa1);
	clearfits(&f_cfa2); clearfits(&f_cfa3);
	return ret;
}

static int cfa_extract_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_input_image, MB_per_output_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_input_image, NULL, &MB_avail);
	const struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;

	if (args->image_hook == extractHa_image_hook || args->image_hook == extractGreen_image_hook) {
		MB_per_output_image = min(1, MB_per_input_image / 4);
		required = 5 * MB_per_input_image / 4;
	}
	else if (args->image_hook == extractHaOIII_image_hook) {
		if (cfa_args->scaling == SCALING_NONE) {
			required = 3 * MB_per_input_image / 2;
			MB_per_output_image = 5 * MB_per_input_image / 4;
		} else {
			if (cfa_args->scaling == SCALING_HA_UP)
				MB_per_output_image = MB_per_input_image * 2;
			else MB_per_output_image = min(1, MB_per_input_image / 2);

			// Very slightly less for upscaling Ha but this is close enough
			required = 7 * MB_per_input_image / 2;
		}
	} else if (args->image_hook == split_cfa_image_hook) {
		required = 2 * MB_per_input_image;
		MB_per_output_image = MB_per_input_image;	// no writer for this split
	} else {
		required = MB_per_input_image;
		MB_per_output_image = MB_per_input_image;
		siril_log_color_message("unknown extraction type\n", "red");
	}

	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_output_image;
		}
		else limit = thread_limit;
	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
#ifdef _OPENMP
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_input_image, limit, for_writer ? "images" : "threads");
#else
		/* we still want the check of limit = 0 above */
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

void apply_split_cfa_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->compute_mem_limits_hook = cfa_extract_compute_mem_limits;
	args->image_hook = split_cfa_image_hook;
	args->description = _("Split CFA");
	args->new_seq_prefix = split_cfa_args->seqEntry;
	args->user = split_cfa_args;
	args->finalize_hook = split_finalize_hook;
	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

#define SQRTF_2 1.41421356f

int extractHaOIII_ushort(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern, extraction_scaling scaling, int threads) {
	int crop_x = in->rx, crop_y = in->ry;
	gboolean do_crop = FALSE;
	if (in->rx % 2) {
		crop_x--;
		do_crop = TRUE;
	}
	if (in->ry % 2) {
		crop_y--;
		do_crop = TRUE;
	}
	if (do_crop) {
		siril_log_message(_("Image dimensions must be even. Trimming by 1 pixel.\n"));
		rectangle bounds = {0, 0, crop_x, crop_y};
		crop(in, &bounds);
	}

	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Extract_HaOIII does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_USHORT) ||
			new_fit_image(&OIII, in->rx, in->ry, 1, DATA_USHORT)) {
		return 1;
	}
	// Loop through calculating the means of the 3 O-III photosite subchannels
	// Also populate the Ha fits data
	float g1 = 0.f, g2 = 0.f, b = 0.f;
	unsigned j = 0;
	unsigned error = 0;
	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			float c0 = (float) in->data[col + row * in->rx];
			float c1 = (float) in->data[1 + col + row * in->rx];
			float c2 = (float) in->data[col + (1 + row) * in->rx];
			float c3 = (float) in->data[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->data[j] = roundf_to_WORD(c0);
				g1 += c1;
				g2 += c2;
				b += c3;
				break;
			case BAYER_FILTER_BGGR:
				Ha->data[j] = roundf_to_WORD(c3);
				g1 += c1;
				g2 += c2;
				b += c0;
				break;
			case BAYER_FILTER_GRBG:
				Ha->data[j] = roundf_to_WORD(c1);
				g1 += c0;
				g2 += c3;
				b += c2;
				break;
			case BAYER_FILTER_GBRG:
				Ha->data[j] = roundf_to_WORD(c2);
				g1 += c0;
				g2 += c3;
				b += c1;
				break;
			default:
				printf("Should not happen.\n");
				error++;
			}
		j++;
		}
	}
	if (!error) {
		// g1, g2 and b are divided by j to give averages, work out the mean OIII level and calculate scaling ratios
		// for each of the 3 OIII subchannels to equalize them
		g1 /= j;
		g2 /= j;
		b /= j;
		float avgoiii = (g1 + g2 + b) / 3;
		float g1ratio = avgoiii / g1;
		float g2ratio = avgoiii / g2;
		float bratio = avgoiii / b;

		// Loop through to equalize the O-III photosite data and interpolate the O-III values at the Ha photosites
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
{
#pragma omp for simd schedule(static)
#endif
		for (int row = 0; row < in->ry - 1; row += 2) {
			for (int col = 0; col < in->rx - 1; col += 2) {
				switch(pattern) {
					case BAYER_FILTER_RGGB:
						OIII->data[1 + col + row * in->rx] = roundf_to_WORD(g1ratio * in->data[1 + col + row * in->rx]);
						OIII->data[col + (1 + row) * in->rx] = roundf_to_WORD(g2ratio * in->data[col + (1 + row) * in->rx]);
						OIII->data[1 + col + (1 + row) * in->rx] = roundf_to_WORD(bratio * in->data[1 + col + (1 + row) * in->rx]);
						break;
					case BAYER_FILTER_BGGR:
						OIII->data[1 + col + row * in->rx] = roundf_to_WORD(g1ratio * in->data[1 + col + row * in->rx]);
						OIII->data[col + (1 + row) * in->rx] = roundf_to_WORD(g2ratio * in->data[col + (1 + row) * in->rx]);
						OIII->data[col + row * in->rx] = roundf_to_WORD(bratio * in->data[col + row * in->rx]);
						break;
					case BAYER_FILTER_GRBG:
						OIII->data[col + row * in->rx] = roundf_to_WORD(g1ratio * in->data[col + row * in->rx]);
						OIII->data[1 + col + (1 + row) * in->rx] = roundf_to_WORD(g2ratio * in->data[1 + col + (1 + row) * in->rx]);
						OIII->data[col + (1 + row) * in->rx] = roundf_to_WORD(bratio * in->data[col + (1 + row) * in->rx]);
						break;
					case BAYER_FILTER_GBRG:
						OIII->data[col + row * in->rx] = roundf_to_WORD(g1ratio * in->data[col + row * in->rx]);
						OIII->data[1 + col + (1 + row) * in->rx] = roundf_to_WORD(g2ratio * in->data[1 + col + (1 + row) * in->rx]);
						OIII->data[1 + col + row * in->rx] = roundf_to_WORD(bratio * in->data[1 + col + row * in->rx]);
						break;
					default:
						printf("Should not happen.\n");
						error++;
				}
			}
		}
#ifdef _OPENMP
#pragma omp barrier
#endif
		// Separate loop for Ha site interpolation so this works for all Bayer patterns
		if (!error) {
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
			for (int row = 0; row < in->ry - 1; row += 2) {
				for (int col = 0; col < in->rx - 1; col += 2) {
					int HaIndex = 0;
					switch(pattern) {
						case BAYER_FILTER_RGGB:
							HaIndex = col + row * in->rx;
							break;
						case BAYER_FILTER_BGGR:
							HaIndex = 1 + col + (1 + row) * in->rx;
							break;
						case BAYER_FILTER_GRBG:
							HaIndex = 1 + col + row * in->rx;
							break;
						case BAYER_FILTER_GBRG:
							HaIndex = col + (1 + row) * in->rx;
							break;
						default:
							printf("Should not happen.\n");
							error++;
					}
					float interp = 0.f;
					float weight = 0.f;
					gboolean first_y = (HaIndex / in->rx == 0) ? TRUE : FALSE;
					gboolean last_y = (HaIndex / in->rx == in->ry - 1) ? TRUE : FALSE;
					gboolean first_x = (HaIndex % in->rx == 0) ? TRUE : FALSE;
					gboolean last_x = (HaIndex % in->rx == in->rx - 1) ? TRUE : FALSE;
					if (!first_y) {
						interp += OIII->data[HaIndex - in->rx] * SQRTF_2;
						weight += SQRTF_2;
						if (!first_x) {
							interp += (OIII->data[HaIndex - 1] * SQRTF_2);
							interp += OIII->data[(HaIndex - in->rx) - 1];
							weight += (1.f + SQRTF_2);
						}
						if (!last_x) {
							interp += OIII->data[(HaIndex - in->rx) + 1];
							interp += OIII->data[HaIndex + 1] * SQRTF_2;
							weight += (1.f + SQRTF_2);
						}
					} else { // first_y
						if (!first_x) {
							interp += OIII->data[HaIndex - 1] * SQRTF_2;
							weight += SQRTF_2;
						}
						if(!last_x) {
							interp += OIII->data[HaIndex + 1] * SQRTF_2;
							weight += SQRTF_2;
						}
					}
					if (!last_y) {
						interp += OIII->data[HaIndex + in->rx] * SQRTF_2;
						weight += SQRTF_2;
						if(!first_x) {
							interp += OIII->data[(HaIndex + in->rx) - 1];
							weight += 1.f;
						}
						if(!last_x) {
							interp += OIII->data[HaIndex + in->rx + 1];
							weight += 1.f;
						}
					}
					interp /= weight;
					OIII->data[HaIndex] = roundf_to_WORD(interp);
				}
			}
		}
#ifdef _OPENMP
}
#endif
	}
	if (error)
		return 1;
	// Scale images to match: either upsample Ha to match OIII, downsample OIII to match Ha
	// or do nothing. Hardcoded to upscale for now.

	float factorHa = 2.f, factorOIII = 1.f;
	switch (scaling) {
		case 1: // Upsample Ha to OIII size
			verbose_resize_gaussian(Ha, OIII->rx, OIII->ry, OPENCV_LANCZOS4, TRUE);
			factorHa = 1.f;
			break;
		case 2: // Downsample OIII to Ha size
			verbose_resize_gaussian(OIII, Ha->rx, Ha->ry, OPENCV_LANCZOS4, TRUE);
			factorOIII = 2.f;
			break;
		default:
			break;
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha, factorHa);
	update_filter_information(Ha, "Ha", TRUE);
	Ha->history = g_slist_append(Ha->history, g_strdup(_("Ha channel\n")));

	copy_fits_metadata(in, OIII);
	update_sampling_information(OIII, factorOIII);
	update_filter_information(OIII, "OIII", TRUE);
	OIII->history = g_slist_append(OIII->history, g_strdup(_("OIII channel\n")));

	return 0;
}

int extractHaOIII_float(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern, extraction_scaling scaling, int threads) {
	int crop_x = in->rx, crop_y = in->ry;
	gboolean do_crop = FALSE;
	if (in->rx % 2) {
		crop_x--;
		do_crop = TRUE;
	}
	if (in->ry % 2) {
		crop_y--;
		do_crop = TRUE;
	}
	if (do_crop) {
		siril_log_message(_("Image dimensions must be even. Trimming by 1 pixel.\n"));
		rectangle bounds = {0, 0, crop_x, crop_y};
		crop(in, &bounds);
	}

	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->keywords.bayer_pattern) > 4) {
		siril_log_message(_("Extract_HaOIII does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&OIII, in->rx, in->ry, 1, DATA_FLOAT)) {
		return 1;
	}
	// Loop through calculating the means of the 3 O-III photosite subchannels
	// Also populate the Ha fits data
	float g1 = 0.f, g2 = 0.f, b = 0.f;
	unsigned j = 0;
	unsigned error = 0;
	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			float c0 = in->fdata[col + row * in->rx];
			float c1 = in->fdata[1 + col + row * in->rx];
			float c2 = in->fdata[col + (1 + row) * in->rx];
			float c3 = in->fdata[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->fdata[j] = c0;
				g1 += c1;
				g2 += c2;
				b += c3;
				break;
			case BAYER_FILTER_BGGR:
				Ha->fdata[j] = c3;
				g1 += c1;
				g2 += c2;
				b += c0;
				break;
			case BAYER_FILTER_GRBG:
				Ha->fdata[j] = c1;
				g1 += c0;
				g2 += c3;
				b += c2;
				break;
			case BAYER_FILTER_GBRG:
				Ha->fdata[j] = c2;
				g1 += c0;
				g2 += c3;
				b += c1;
				break;
			default:
				printf("Should not happen.\n");
				error++;
			}
		j++;
		}
	}
	if (!error) {
		// g1, g2 and b are divided by j to give averages, work out the mean OIII level and calculate scaling ratios
		// for each of the 3 OIII subchannels to equalize them
		g1 /= j;
		g2 /= j;
		b /= j;
		float avgoiii = (g1 + g2 + b) / 3;
		float g1ratio = avgoiii / g1;
		float g2ratio = avgoiii / g2;
		float bratio = avgoiii / b;

		// Loop through to equalize the O-III photosite data and interpolate the O-III values at the Ha photosites
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
{
#pragma omp for simd schedule(static)
#endif
		for (int row = 0; row < in->ry - 1; row += 2) {
			for (int col = 0; col < in->rx - 1; col += 2) {
				switch(pattern) {
					case BAYER_FILTER_RGGB:
						OIII->fdata[1 + col + row * in->rx] = g1ratio * in->fdata[1 + col + row * in->rx];
						OIII->fdata[col + (1 + row) * in->rx] = g2ratio * in->fdata[col + (1 + row) * in->rx];
						OIII->fdata[1 + col + (1 + row) * in->rx] = bratio * in->fdata[1 + col + (1 + row) * in->rx];
						break;
					case BAYER_FILTER_BGGR:
						OIII->fdata[1 + col + row * in->rx] = g1ratio * in->fdata[1 + col + row * in->rx];
						OIII->fdata[col + (1 + row) * in->rx] = g2ratio * in->fdata[col + (1 + row) * in->rx];
						OIII->fdata[col + row * in->rx] = bratio * in->fdata[col + row * in->rx];
						break;
					case BAYER_FILTER_GRBG:
						OIII->fdata[col + row * in->rx] = g1ratio * in->fdata[col + row * in->rx];
						OIII->fdata[1 + col + (1 + row) * in->rx] = g2ratio * in->fdata[1 + col + (1 + row) * in->rx];
						OIII->fdata[col + (1 + row) * in->rx] = bratio * in->fdata[col + (1 + row) * in->rx];
						break;
					case BAYER_FILTER_GBRG:
						OIII->fdata[col + row * in->rx] = g1ratio * in->fdata[col + row * in->rx];
						OIII->fdata[1 + col + (1 + row) * in->rx] = g2ratio * in->fdata[1 + col + (1 + row) * in->rx];
						OIII->fdata[1 + col + row * in->rx] = bratio * in->fdata[1 + col + row * in->rx];
						break;
					default:
						printf("Should not happen.\n");
						error++;
				}
			}
		}
#ifdef _OPENMP
#pragma omp barrier
#endif
		// Separate loop for Ha site interpolation so this works for all Bayer patterns
		if (!error) {
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
			for (int row = 0; row < in->ry - 1; row += 2) {
				for (int col = 0; col < in->rx - 1; col += 2) {
					int HaIndex = 0;
					switch(pattern) {
						case BAYER_FILTER_RGGB:
							HaIndex = col + row * in->rx;
							break;
						case BAYER_FILTER_BGGR:
							HaIndex = 1 + col + (1 + row) * in->rx;
							break;
						case BAYER_FILTER_GRBG:
							HaIndex = 1 + col + row * in->rx;
							break;
						case BAYER_FILTER_GBRG:
							HaIndex = col + (1 + row) * in->rx;
							break;
						default:
							printf("Should not happen.\n");
							error++;
					}
					float interp = 0.f;
					float weight = 0.f;
					gboolean first_y = (HaIndex / in->rx == 0) ? TRUE : FALSE;
					gboolean last_y = (HaIndex / in->rx == in->ry - 1) ? TRUE : FALSE;
					gboolean first_x = (HaIndex % in->rx == 0) ? TRUE : FALSE;
					gboolean last_x = (HaIndex % in->rx == in->rx - 1) ? TRUE : FALSE;
					if (!first_y) {
						interp += OIII->fdata[HaIndex - in->rx] * SQRTF_2;
						weight += SQRTF_2;
						if (!first_x) {
							interp += (OIII->fdata[HaIndex - 1] * SQRTF_2);
							interp += OIII->fdata[HaIndex - in->rx - 1];
							weight += (1 + SQRTF_2);
						}
						if (!last_x) {
							interp += OIII->fdata[HaIndex - in->rx + 1];
							interp += OIII->fdata[HaIndex + 1] * SQRTF_2;
							weight += (1 + SQRTF_2);
						}
					} else { // first_y
						if (!first_x) {
							interp += OIII->fdata[HaIndex - 1] * SQRTF_2;
							weight += SQRTF_2;
						}
						if(!last_x) {
							interp += OIII->fdata[HaIndex+1] * SQRTF_2;
							weight += SQRTF_2;
						}
					}
					if (!last_y) {
						interp += OIII->fdata[HaIndex + in->rx] * SQRTF_2;
						weight += SQRTF_2;
						if(!first_x) {
							interp += OIII->fdata[HaIndex + in->rx - 1];
							weight += 1.f;
						}
						if(!last_x) {
							interp += OIII->fdata[HaIndex + in->rx + 1];
							weight += 1.f;
						}
					}
					interp /= weight;
					OIII->fdata[HaIndex] = interp;
				}
			}
		}
#ifdef _OPENMP
}
#endif
	}
	if (error)
		return 1;
	// Scale images to match: either upsample Ha to match OIII, downsample OIII to match Ha
	// or do nothing. Hardcoded to upscale for now.
	float factorHa = 2.f, factorOIII = 1.f;
	switch (scaling) {
		case SCALING_HA_UP: // Upsample Ha to OIII size
			verbose_resize_gaussian(Ha, OIII->rx, OIII->ry, OPENCV_LANCZOS4, TRUE);
			factorHa = 1.f;
			break;
		case SCALING_OIII_DOWN: // Downsample OIII to Ha size
			verbose_resize_gaussian(OIII, Ha->rx, Ha->ry, OPENCV_LANCZOS4, TRUE);
			factorOIII = 2.f;
			break;
		default:
			break;
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha, factorHa);
	update_filter_information(Ha, "Ha", TRUE);
	Ha->history = g_slist_append(Ha->history, g_strdup(_("Ha channel\n")));

	copy_fits_metadata(in, OIII);
	update_sampling_information(OIII, factorOIII);
	update_filter_information(OIII, "OIII", TRUE);
	OIII->history = g_slist_append(OIII->history, g_strdup(_("OIII channel\n")));

	return 0;
}
