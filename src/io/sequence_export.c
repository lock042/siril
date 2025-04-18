/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "sequence.h"
#include "ser.h"
#include "stacking/stacking.h"
#include "registration/registration.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/stacking.h"
#include "io/image_format_fits.h"
#include "io/Astro-TIFF.h"
#ifdef HAVE_FFMS2
#include "io/films.h"
#endif
#include "avi_pipp/avi_writer.h"
#ifdef HAVE_FFMPEG
#include "io/mp4_output.h"
#endif
#include "algos/geometry.h"

/* same order as in the combo box 'combo_export_preset' */
typedef enum {
	PRESET_RATIO,
	PRESET_FREE,
	PRESET_FULLHD,
	PRESET_4K,
	PRESET_8K
} export_presets;

struct exportseq_args {
	sequence *seq;
	seq_image_filter filtering_criterion;
	double filtering_parameter;

	char *basename;
	export_format output;
	gboolean normalize;

	gboolean tiff_compression; // compression of TIFF files

	int film_fps;		// has to be int for avi or ffmpeg
	int film_quality;	// [1, 5], for mp4 and webm

	gboolean resample;
	int32_t dest_width, dest_height;

	gboolean crop;
	rectangle crop_area;
};


/* Used for avi exporter, creates buffer as BGRBGR from ushort FITS */
static uint8_t *fits_to_uint8(fits *fit) {
	size_t i, j;
	size_t n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	int channel = fit->naxes[2];
	int step = (channel == 3 ? 2 : 0);

	uint8_t *data = malloc(n);
	if (fit->orig_bitpix == BYTE_IMG) {
		for (i = 0, j = 0; i < n; i += channel, j++) {
			data[i + step] = (BYTE)fit->pdata[RLAYER][j];
			if (channel > 1) {
				data[i + 1] = (BYTE)fit->pdata[GLAYER][j];
				data[i + 2 - step] = (BYTE)fit->pdata[BLAYER][j];
			}
		}
	} else {
		siril_debug_print("converting from ushort\n");
		for (i = 0, j = 0; i < n; i += channel, j++) {
			data[i + step] = (BYTE)(fit->pdata[RLAYER][j] >> 8);
			if (channel > 1) {
				data[i + 1] = (BYTE)(fit->pdata[GLAYER][j] >> 8);
				data[i + 2 - step] = (BYTE)(fit->pdata[BLAYER][j] >> 8);
			}
		}
	}
	return data;
}

static gpointer export_sequence(gpointer ptr) {
	int retval = 0, cur_nb = 0, j = 0;
	unsigned int out_width, out_height, in_width, in_height;
	uint8_t *data;
	fits *destfit = NULL;	// must be declared before any goto!
	char filename[256], dest[256];
	struct ser_struct *ser_file = NULL;
	fitseq *fitseq_file = NULL;
	GSList *timestamp = NULL;
	char *filter_descr;
	int nb_frames = 0;
	GDateTime *strTime;
	gboolean aborted = FALSE;
	cmsHPROFILE ref_icc = NULL;
#ifdef HAVE_FFMPEG
	struct mp4_struct *mp4_file = NULL;
#endif
	struct exportseq_args *args = (struct exportseq_args *)ptr;
	norm_coeff coeff = { 0 };

	int reglayer = get_registration_layer(args->seq);
	if (reglayer >= 0)
		siril_log_message(_("Using registration information from layer %d to export sequence\n"), reglayer);
	if (args->crop) {
		in_width  = args->crop_area.w;
		in_height = args->crop_area.h;
	} else {
		in_width  = args->seq->rx;
		in_height = args->seq->ry;
	}

	if (args->resample) {
		out_width = args->dest_width;
		out_height = args->dest_height;
		if (out_width == in_width && out_height == in_height)
			args->resample = FALSE;
	} else {
		out_width = in_width;
		out_height = in_height;
	}

	gboolean have_seqwriter = args->output == EXPORT_FITSEQ || args->output == EXPORT_SER;
	if (have_seqwriter)
		seqwriter_set_max_active_blocks(3);

	int output_bitpix = USHORT_IMG;

	/* possible output formats: FITS images, FITS cube, TIFF, SER, AVI, MP4, WEBM */
	// create the sequence file for single-file sequence formats
	fits ref = { 0 };
	int refindex = 0;
	switch (args->output) {
		case EXPORT_FITS:
			output_bitpix = args->seq->bitpix;
			refindex = sequence_find_refimage(args->seq);
			if (!seq_read_frame(args->seq, refindex, &ref, FALSE, -1) && ref.icc_profile) {
				ref_icc = copyICCProfile(ref.icc_profile);
				siril_log_message(_("Reference frame has an ICC profile. Will assign / convert other frames to to match.\n"));
			}
			clearfits(&ref);
			break;
		case EXPORT_FITSEQ:
			fitseq_file = calloc(1, sizeof(fitseq));
			snprintf(dest, 256, "%s%s", args->basename, com.pref.ext);
			if (fitseq_create_file(dest, fitseq_file, -1)) {
				free(fitseq_file);
				fitseq_file = NULL;
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			output_bitpix = args->seq->bitpix;
			break;
		case EXPORT_TIFF:
			output_bitpix = args->seq->bitpix;
			// limit to 16 bits for TIFF
			if (output_bitpix == FLOAT_IMG)
				output_bitpix = USHORT_IMG;
			break;

		case EXPORT_SER:
			ser_file = calloc(1, sizeof(struct ser_struct));
			snprintf(dest, 256, "%s.ser", args->basename);
			if (ser_create_file(dest, ser_file, TRUE, args->seq->ser_file)) {
				free(ser_file);
				ser_file = NULL;
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			output_bitpix = args->seq->bitpix;
			if (output_bitpix == FLOAT_IMG)
				output_bitpix = USHORT_IMG;
			break;

		case EXPORT_AVI:
			snprintf(dest, 256, "%s.avi", args->basename);
			int32_t avi_format;

			// Check if the sequence has an ICC profile. If so, we should convert to sRGB
			// as that's really the only suitable option here
			refindex = sequence_find_refimage(args->seq);
			if (!seq_read_frame(args->seq, refindex, &ref, FALSE, -1) && ref.icc_profile) {
				ref_icc = copyICCProfile(ref.icc_profile);
				siril_log_message(_("Reference frame has an ICC profile. Exporting as sRGB.\n"));
			}
			clearfits(&ref);

			if (args->seq->nb_layers == 1)
				avi_format = AVI_WRITER_INPUT_FORMAT_MONOCHROME;
			else avi_format = AVI_WRITER_INPUT_FORMAT_COLOUR;

			if (avi_file_create(dest, out_width, out_height, avi_format,
					AVI_WRITER_CODEC_DIB, args->film_fps)) {
				siril_log_color_message(_("AVI file `%s' could not be created\n"), "red", dest);
				retval = -1;
				goto free_and_reset_progress_bar;
			}

			output_bitpix = BYTE_IMG;
			break;

		case EXPORT_MP4:
		case EXPORT_MP4_H265:
		case EXPORT_WEBM_VP9:
#ifndef HAVE_FFMPEG
			siril_log_message(_("MP4 output is not supported because siril was not compiled with ffmpeg support.\n"));
			retval = -1;
			goto free_and_reset_progress_bar;
#else
			// Check if the sequence has an ICC profile. If so, we should convert to sRGB
			// as that's really the only suitable option here
			refindex = sequence_find_refimage(args->seq);
			if (!seq_read_frame(args->seq, refindex, &ref, FALSE, -1) && ref.icc_profile) {
				ref_icc = copyICCProfile(ref.icc_profile);
				siril_log_message(_("Reference frame has an ICC profile. Exporting as sRGB.\n"));
			}
			clearfits(&ref);

			/* resampling is managed by libswscale */
			snprintf(dest, 256, "%s.%s", args->basename,
					args->output == EXPORT_WEBM_VP9 ? "webm" : "mp4");

			if (in_width % 32 || out_height % 2 || out_width % 2) {
				siril_log_message(_("Film output needs to have a width that is a multiple of 32 and an even height, resizing selection.\n"));
				if (in_width % 32) in_width = (in_width / 32) * 32 + 32;
				if (in_height % 2) in_height++;
				if (args->crop) {
					args->crop_area.w = in_width;
					args->crop_area.h = in_height;
				} else {
					args->crop = TRUE;
					args->crop_area.x = 0;
					args->crop_area.y = 0;
					args->crop_area.w = in_width;
					args->crop_area.h = in_height;
				}
				compute_fitting_selection(&args->crop_area, 32, 2, 0);
				memcpy(&com.selection, &args->crop_area, sizeof(rectangle));
				fprintf(stdout, "final input area: %d,%d,\t%dx%d\n",
						args->crop_area.x, args->crop_area.y,
						args->crop_area.w, args->crop_area.h);
				in_width = args->crop_area.w;
				in_height = args->crop_area.h;
				if (!args->resample) {
					out_width = in_width;
					out_height = in_height;
				} else {
					if (out_width % 2) out_width++;
					if (out_height % 2) out_height++;
				}
			}

			mp4_file = mp4_create(dest, out_width, out_height, args->film_fps, args->seq->nb_layers, args->film_quality, in_width, in_height, args->output);
			if (!mp4_file) {
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			output_bitpix = BYTE_IMG;
			break;
#endif
	}

	if (output_bitpix == FLOAT_IMG)
		output_bitpix = com.pref.force_16bit ? USHORT_IMG : FLOAT_IMG;

	nb_frames = compute_nb_filtered_images(args->seq,
			args->filtering_criterion, args->filtering_parameter);
	filter_descr = describe_filter(args->seq, args->filtering_criterion,
			args->filtering_parameter);
	siril_log_message(filter_descr);
	g_free(filter_descr);

	if (args->normalize) {
		struct stacking_args stackargs = { 0 };
		stackargs.force_norm = FALSE;
		stackargs.seq = args->seq;
		stackargs.filtering_criterion = args->filtering_criterion;
		stackargs.filtering_parameter = args->filtering_parameter;
		stackargs.nb_images_to_stack = nb_frames;
		stackargs.normalize = ADDITIVE_SCALING;
		stackargs.reglayer = reglayer;
		stackargs.use_32bit_output = (output_bitpix == FLOAT_IMG);
		stackargs.ref_image = args->seq->reference_image;
		stackargs.equalizeRGB = FALSE;
		stackargs.lite_norm = FALSE;

		// build image indices used by normalization
		if (stack_fill_list_of_unfiltered_images(&stackargs))
			goto free_and_reset_progress_bar;

		do_normalization(&stackargs);
		coeff.offset = stackargs.coeff.offset;
		coeff.scale = stackargs.coeff.scale;
		for (int layer = 0; layer < args->seq->nb_layers; ++layer) {
			coeff.poffset[layer] = stackargs.coeff.poffset[layer];
			coeff.pscale[layer] = stackargs.coeff.pscale[layer];
		    for (int i = 0; i < nb_frames; ++i) {
				coeff.poffset[layer][i] = stackargs.coeff.poffset[layer][i];
				coeff.pscale[layer][i] = stackargs.coeff.pscale[layer][i];
			}
		}
		free(stackargs.coeff.mul);

		// and dispose them, because we don't need them anymore
		free(stackargs.image_indices);
	}

	long naxes[3];
	size_t nbpix = 0;
	set_progress_bar_data(NULL, PROGRESS_RESET);
	gboolean icc_msg_given = FALSE;
	for (int i = 0, skipped = 0; i < args->seq->number; ++i) {
		if (!get_thread_run()) {
			retval = -1;
			aborted = TRUE;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, i, args->filtering_parameter)) {
			siril_log_message(_("image %d is excluded from export\n"), i + 1);
			skipped++;
			continue;
		}

		if (have_seqwriter)
			seqwriter_wait_for_memory();

		if (!seq_get_image_filename(args->seq, i, filename)) {
			seqwriter_release_memory();
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		gchar *tmpmsg = g_strdup_printf(_("Processing image %s"), filename);
		set_progress_bar_data(tmpmsg, (double)cur_nb / (nb_frames == 0 ? 1 : (double)nb_frames));
		g_free(tmpmsg);

		/* we read the full frame */
		fits fit = { 0 };
		if (seq_read_frame(args->seq, i, &fit, FALSE, -1)) {
			seqwriter_release_memory();
			siril_log_message(_("Export: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		/* destfit is allocated to the full size. Data will be copied from fit,
		 * image buffers are duplicated. It will be cropped after the copy if
		 * needed */
		if (!nbpix || have_seqwriter) {
			if (!nbpix) {
				memcpy(naxes, fit.naxes, sizeof naxes);
				nbpix = fit.naxes[0] * fit.naxes[1];
			}
			else {
				if (memcmp(naxes, fit.naxes, sizeof naxes)) {
					siril_log_color_message(_("An image of the sequence doesn't have the same dimensions\n"), "red");
					retval = -3;
					clearfits(&fit);
					seqwriter_release_memory();
					goto free_and_reset_progress_bar;
				}
			}
			data_type dest_type = output_bitpix == FLOAT_IMG ? DATA_FLOAT : DATA_USHORT;
			destfit = NULL;	// otherwise it's overwritten by new_fit_image()
			if (new_fit_image(&destfit, fit.rx, fit.ry, fit.naxes[2], dest_type)) {
				retval = -1;
				clearfits(&fit);
				seqwriter_release_memory();
				goto free_and_reset_progress_bar;
			}
			destfit->bitpix = output_bitpix;
			destfit->orig_bitpix = fit.orig_bitpix;
		}
		else if (memcmp(naxes, fit.naxes, sizeof naxes)) {
			fprintf(stderr, "An image of the sequence doesn't have the same dimensions\n");
			retval = -3;
			clearfits(&fit);
			seqwriter_release_memory();
			goto free_and_reset_progress_bar;
		}
		else {
			/* we don't have a seq writer or it's not the first frame, we can reuse destfit */
			if (destfit->type == DATA_FLOAT) {
				memset(destfit->fdata, 0, nbpix * fit.naxes[2] * sizeof(float));
				if (args->crop) {
					/* reset destfit damaged by the crop function */
					if (fit.naxes[2] == 3) {
						destfit->fpdata[1] = destfit->fdata + nbpix;
						destfit->fpdata[2] = destfit->fdata + nbpix * 2;
					}
					destfit->rx = destfit->naxes[0] = fit.rx;
					destfit->ry = destfit->naxes[1] = fit.ry;
				}
			} else {
				memset(destfit->data, 0, nbpix * fit.naxes[2] * sizeof(WORD));
				if (args->crop) {
					/* reset destfit damaged by the crop function */
					if (fit.naxes[2] == 3) {
						destfit->pdata[1] = destfit->data + nbpix;
						destfit->pdata[2] = destfit->data + nbpix * 2;
					}
					destfit->rx = destfit->naxes[0] = fit.rx;
					destfit->ry = destfit->naxes[1] = fit.ry;
				}
			}
		}

		// Copy the ICC profile from fit if available
		if (destfit->icc_profile) {
			cmsCloseProfile(destfit->icc_profile);
			destfit->icc_profile = copyICCProfile(fit.icc_profile);
		} else {
		// Otherwise see if the reference frame has an ICC profile, and if so copy that
			if (ref_icc) {
				destfit->icc_profile = copyICCProfile(ref_icc);
			}
		}
		color_manage(destfit, destfit->icc_profile != NULL);


		// we copy the header
		copy_fits_metadata(&fit, destfit);

		int shiftx, shifty;
		/* load registration data for current image */
		if (reglayer != -1 && args->seq->regparam[reglayer]) {
			double dx, dy;
			translation_from_H(args->seq->regparam[reglayer][i].H, &dx, &dy);
			shiftx = round_to_int(dx);
			shifty = round_to_int(dy);
		} else {
			shiftx = 0;
			shifty = 0;
		}

		/* fill the image with shifted data and normalization */
		for (int layer = 0; layer < fit.naxes[2]; ++layer) {
			for (int y = 0; y < fit.ry; ++y) {
				for (int x = 0; x < fit.rx; ++x) {
					int nx = x + shiftx;
					int ny = y + shifty;
					if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
						if (fit.type == DATA_USHORT) {
							WORD pixel = fit.pdata[layer][x + y * fit.rx];
							if (args->normalize) {
								double tmp = (double) pixel;
								if (pixel > 0) { // do not offset null pixels
									tmp *= coeff.pscale[layer][j];
									tmp -= (coeff.poffset[layer][j]);
									pixel = round_to_WORD(tmp);
								}
							}
							destfit->pdata[layer][nx + ny * fit.rx] = pixel;
						}
						else if (fit.type == DATA_FLOAT) {
							float pixel = fit.fpdata[layer][x + y * fit.rx];
							if (args->normalize) {
								if (pixel != 0.f) { // do not offset null pixels
									pixel *= (float) coeff.pscale[layer][j];
									pixel -= (float) coeff.poffset[layer][j];
								}
							}
							if (destfit->type == DATA_FLOAT) {
								destfit->fpdata[layer][nx + ny * fit.rx] = pixel;
							} else {
								destfit->pdata[layer][nx + ny * fit.rx] = roundf_to_WORD(pixel * USHRT_MAX_SINGLE);
							}
						} else {
							retval = -1;
							clearfits(&fit);
							seqwriter_release_memory();
							goto free_and_reset_progress_bar;
						}
						// for 8 bit output, destfit will be transformed later with a linear scale
					}
				}
			}
		}
		j++;
		clearfits(&fit);

		if (args->crop) {
			crop(destfit, &args->crop_area);
		}

		switch (args->output) {
			case EXPORT_FITS:
				if (ref_icc) {
					siril_colorspace_transform(destfit, ref_icc);
				} else {
					if (destfit->icc_profile && !icc_msg_given) {
						siril_log_message(_("Info: this frame has an ICC profile but the reference frame does not. Profile will be preserved...\n"));
						icc_msg_given = TRUE;
					}
				}
				snprintf(dest, 255, "%s%05d%s", args->basename, i + 1, com.pref.ext);
				retval = savefits(dest, destfit);
				break;
			case EXPORT_FITSEQ:
				retval = fitseq_write_image(fitseq_file, destfit, i - skipped);
				break;
#ifdef HAVE_LIBTIFF
			case EXPORT_TIFF:
				snprintf(dest, 255, "%s%05d", args->basename, i + 1);
				gchar *astro_tiff = AstroTiff_build_header(destfit);
				retval = savetif(dest, destfit, 16, astro_tiff, com.pref.copyright, args->tiff_compression, TRUE, TRUE);
				g_free(astro_tiff);
				break;
#endif
			case EXPORT_SER:
				if (destfit->keywords.date_obs) {
					strTime = g_date_time_ref(destfit->keywords.date_obs);
					timestamp = g_slist_append(timestamp, strTime);
				}
				retval = ser_write_frame_from_fit(ser_file, destfit, i - skipped);
				break;
			case EXPORT_AVI:
				if (ref_icc) {
					siril_colorspace_transform(destfit, srgb_trc());
				}
				data = fits_to_uint8(destfit);
				retval = avi_file_write_frame(0, data);
				break;
#ifdef HAVE_FFMPEG
			case EXPORT_MP4:
			case EXPORT_MP4_H265:
			case EXPORT_WEBM_VP9:
				// an equivalent to fits_to_uint8 is called in there (fill_rgb_image)...
				if (ref_icc) {
					siril_colorspace_transform(destfit, srgb_trc());
				}
				retval = mp4_add_frame(mp4_file, destfit);
				break;
#endif
			default:
				fprintf(stderr, "Case not handled, should not happen.\n");
				seqwriter_release_memory();
				goto free_and_reset_progress_bar;
		}
		if (retval) {
			seqwriter_release_memory();
			goto free_and_reset_progress_bar;
		}
		cur_nb++;
	}

free_and_reset_progress_bar:
	cmsCloseProfile(ref_icc);
	ref_icc = NULL;
	if (destfit && !have_seqwriter)
		clearfits(destfit);
	if (args->normalize) {
		free(coeff.offset);
		free(coeff.scale);
	}
	// close the sequence file for single-file sequence formats
	switch (args->output) {
		case EXPORT_FITSEQ:
			if (fitseq_file) {
				fitseq_close_file(fitseq_file);
				free(fitseq_file);
			}
			break;
		case EXPORT_SER:
			if (ser_file) {
				if (timestamp)
					ser_convertTimeStamp(ser_file, timestamp);
				ser_write_and_close(ser_file);
				free(ser_file);
			}
			g_slist_free_full(timestamp, (GDestroyNotify) g_date_time_unref);
			break;
		case EXPORT_AVI:
			avi_file_close(0);
			break;
		case EXPORT_MP4:
		case EXPORT_MP4_H265:
		case EXPORT_WEBM_VP9:
#ifdef HAVE_FFMPEG
			if (mp4_file) {
				mp4_close(mp4_file, aborted);
				free(mp4_file);
			}
#endif
			break;
		default:
			break;
	}

	if (retval) {
		set_progress_bar_data(_("Sequence export failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Sequence export failed\n"));
	}
	else {
		set_progress_bar_data(_("Sequence export succeeded."), PROGRESS_RESET);
		siril_log_message(_("Sequence export succeeded.\n"));
	}
	if (aborted) { //disposing of the file if Stop button was hit
		switch (args->output){
			case EXPORT_FITSEQ:
			case EXPORT_SER:
			case EXPORT_AVI:
			case EXPORT_MP4:
			case EXPORT_MP4_H265:
			case EXPORT_WEBM_VP9:
				if (!(dest[0] == '\0')) {
					GFile *file = g_file_new_build_filename(dest, NULL);
					g_autoptr(GError) local_error = NULL;
					if (!g_file_delete(file, NULL, &local_error)) {
						g_warning("Failed to delete %s: %s", g_file_peek_path(file), local_error->message);
					}
					g_object_unref(file);
				}
				break;
			default:
				break;
		}
	}

	free(args->basename);
	free(args);
	args = NULL;
	siril_add_idle(end_generic, args);
	return NULL;
}

void on_buttonExportSeq_clicked(GtkButton *button, gpointer user_data) {
	int selected = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboExport")));
	GtkEntry *entry = GTK_ENTRY(lookup_widget("entryExportSeq"));
	const char *bname = gtk_entry_get_text(entry);
	gboolean normalize = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("exportNormalize")));

	if (bname[0] == '\0') {
		widget_set_class(GTK_WIDGET(entry), "warning", "");
		return;
	}
	if (selected == -1) return;

	// checking regdata is absent, or if present, is only shift
	if (!test_regdata_is_valid_and_shift(&com.seq, get_registration_layer(&com.seq))) {
		int confirm = siril_confirm_dialog(_("Registration data found"),
			_("Export has detected registration data with more than simple shifts.\n"
			"Normally, you should apply existing registration before exporting."),
			_("Export anyway"));
		if (!confirm)
			return;
	}

	struct exportseq_args *args = calloc(1, sizeof(struct exportseq_args));
	args->seq = &com.seq;
	get_sequence_filtering_from_gui(&args->filtering_criterion, &args->filtering_parameter);
	args->basename = g_str_to_ascii(bname, NULL);
	args->output = (export_format)selected;
	args->normalize = normalize;
	args->resample = FALSE;
	args->crop = com.selection.w && com.selection.h;
	if (args->crop)
		memcpy(&args->crop_area, &com.selection, sizeof(rectangle));

	if (args->output == EXPORT_AVI || args->output == EXPORT_MP4 || args->output == EXPORT_MP4_H265 || args->output == EXPORT_WEBM_VP9) {
		GtkEntry *fpsEntry = GTK_ENTRY(lookup_widget("entryAviFps"));
		args->film_fps = round_to_int(g_ascii_strtod(gtk_entry_get_text(fpsEntry), NULL));
		if (args->film_fps <= 0) args->film_fps = 1;
	}
	if (args->output == EXPORT_MP4 || args->output == EXPORT_MP4_H265 || args->output == EXPORT_WEBM_VP9) {
		GtkAdjustment *adjQual = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder,"adjustmentqualscale"));
		GtkToggleButton *checkResize = GTK_TOGGLE_BUTTON(lookup_widget("checkAviResize"));
		args->film_quality = (int)gtk_adjustment_get_value(adjQual);
		args->resample = gtk_toggle_button_get_active(checkResize);
		if (args->resample) {
			GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));
			GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));
			args->dest_width = g_ascii_strtoll(gtk_entry_get_text(widthEntry), NULL, 10);
			args->dest_height = g_ascii_strtoll(gtk_entry_get_text(heightEntry), NULL, 10);
			if (args->dest_height == 0 || args->dest_width == 0) {
				siril_log_message(_("Width or height cannot be null. Not resizing.\n"));
				gtk_toggle_button_set_active(checkResize, FALSE);
				args->resample = FALSE;
			} else if (args->dest_height == args->seq->ry && args->dest_width == args->seq->rx) {
				gtk_toggle_button_set_active(checkResize, FALSE);
				args->resample = FALSE;
			}
		}
	}
	else if (args->output == EXPORT_FITS || args->output == EXPORT_TIFF) {
		// add a trailing '_' for multiple-files sequences
		args->basename = format_basename(args->basename, TRUE);
		if (args->output == EXPORT_TIFF) {
#ifdef HAVE_LIBTIFF
			args->tiff_compression = get_tiff_compression();
#endif
		}
	}
	// Display a useful warning because I always forget to remove selection
	if (args->crop) {
		gboolean confirm = siril_confirm_dialog(_("Export cropped sequence?"),
				_("An active selection was detected. The exported sequence will only contain data within the drawn selection. You can confirm the crop or cancel it. "
						"If you choose to click on cancel, the exported sequence will contain all data."), _("Confirm Crop"));
		args->crop = confirm;
	}

	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(export_sequence, args)) {
		g_free(args->basename);
		free(args);
	}
}

void on_comboExport_changed(GtkComboBox *box, gpointer user_data) {
	GtkWidget *avi_options = lookup_widget("boxAviOptions");
	GtkWidget *checkAviResize = lookup_widget("checkAviResize");
	GtkWidget *quality = lookup_widget("exportQualScale");
	int output_type = gtk_combo_box_get_active(box);
	gtk_widget_set_visible(avi_options, output_type >= EXPORT_AVI);
	gtk_widget_set_visible(quality, output_type >= EXPORT_MP4);
	gtk_widget_set_sensitive(checkAviResize, output_type >= EXPORT_MP4);
}

void on_checkAviResize_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *heightEntry = lookup_widget("entryAviHeight");
	GtkWidget *widthEntry = lookup_widget("entryAviWidth");
	GtkWidget *combo_export_preset = lookup_widget("combo_export_preset");
	gtk_widget_set_sensitive(heightEntry, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(widthEntry, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(combo_export_preset, gtk_toggle_button_get_active(togglebutton));
}

void update_export_crop_label() {
	static GtkLabel *label = NULL;
	if (!label)
		label = GTK_LABEL(lookup_widget("exportLabel"));
	if (com.selection.w && com.selection.h)
		gtk_label_set_text(label, _("Cropping to selection"));
	else gtk_label_set_text(label, _("Select area to crop"));
}

void on_entryExportSeq_changed(GtkEditable *editable, gpointer user_data){
	gchar *name = (gchar *)gtk_entry_get_text(GTK_ENTRY(editable));
	if (*name != 0) {
		if (check_if_seq_exist(name, !g_str_has_suffix(name, ".ser"))) {
			set_icon_entry(GTK_ENTRY(editable), "dialog-warning");
		} else {
			set_icon_entry(GTK_ENTRY(editable), NULL);
		}
	} else {
		set_icon_entry(GTK_ENTRY(editable), NULL);
	}
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data);

void on_entryAviWidth_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_height;
	GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.ry / (double) com.seq.rx;
	width = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(editable)),NULL);
	height = ratio * width;
	c_height = g_strdup_printf("%d", (int)(height));

	g_signal_handlers_block_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	gtk_entry_set_text(heightEntry, c_height);
	g_signal_handlers_unblock_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	g_free(c_height);
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_width;
	GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.rx / (double) com.seq.ry;
	height = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(editable)), NULL);
	width = ratio * height;
	c_width = g_strdup_printf("%d", (int)(width));

	g_signal_handlers_block_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	gtk_entry_set_text(widthEntry, c_width);
	g_signal_handlers_unblock_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	g_free(c_width);
}

static gboolean g_signal_handlers_is_blocked_by_func(gpointer instance, GFunc func, gpointer data) {
	return g_signal_handler_find(instance,
			G_SIGNAL_MATCH_FUNC | G_SIGNAL_MATCH_DATA
					| G_SIGNAL_MATCH_UNBLOCKED, 0, 0, NULL, func, data) == 0;
}

void on_combo_export_preset_changed(GtkComboBox *box, gpointer user_data) {
	int preset = gtk_combo_box_get_active(box);
	GtkWidget *widthEntry = lookup_widget("entryAviWidth");
	GtkWidget *heightEntry = lookup_widget("entryAviHeight");

	switch(preset) {
	default:
	case PRESET_RATIO:
		if (g_signal_handlers_is_blocked_by_func(GTK_ENTRY(widthEntry), (GFunc) on_entryAviWidth_changed, NULL)) {
			g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
			g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		}
		gtk_widget_set_sensitive(widthEntry, TRUE);
		gtk_widget_set_sensitive(heightEntry, TRUE);
		break;
	case PRESET_FREE:
		if (!g_signal_handlers_is_blocked_by_func(GTK_ENTRY(widthEntry), (GFunc) on_entryAviWidth_changed, NULL)) {
			g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
			g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		}
		gtk_widget_set_sensitive(widthEntry, TRUE);
		gtk_widget_set_sensitive(heightEntry, TRUE);
		break;
	case PRESET_FULLHD:
		g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		gtk_entry_set_text(GTK_ENTRY(widthEntry), "1920");
		gtk_entry_set_text(GTK_ENTRY(heightEntry), "1080");
		gtk_widget_set_sensitive(widthEntry, FALSE);
		gtk_widget_set_sensitive(heightEntry, FALSE);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		break;
	case PRESET_4K:
		g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		gtk_entry_set_text(GTK_ENTRY(widthEntry), "3840");
		gtk_entry_set_text(GTK_ENTRY(heightEntry), "2160");
		gtk_widget_set_sensitive(widthEntry, FALSE);
		gtk_widget_set_sensitive(heightEntry, FALSE);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		break;
	case PRESET_8K:
		g_signal_handlers_block_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_block_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
		gtk_entry_set_text(GTK_ENTRY(widthEntry), "7680");
		gtk_entry_set_text(GTK_ENTRY(heightEntry), "4320");
		gtk_widget_set_sensitive(widthEntry, FALSE);
		gtk_widget_set_sensitive(heightEntry, FALSE);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(widthEntry), on_entryAviWidth_changed, NULL);
		g_signal_handlers_unblock_by_func(GTK_ENTRY(heightEntry), on_entryAviHeight_changed, NULL);
	}
}
