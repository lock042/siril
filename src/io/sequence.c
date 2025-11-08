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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <libgen.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "core/siril_log.h"
#include "io/conversion.h"
#include "gui/utils.h"
#include "gui/cut.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/plot.h"
#include "gui/registration.h"
#include "ser.h"
#include "fits_sequence.h"
#ifdef HAVE_FFMS2
#include "films.h"
#endif
#include "single_image.h"
#include "image_format_fits.h"
#include "gui/histogram.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"	// clear_stars_list
#include "gui/sequence_list.h"
#include "gui/registration_preview.h"
#include "gui/stacking.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"
#include "algos/demosaicing.h"
#include "registration/registration.h"
#include "stacking/stacking.h"	// for update_stack_interface
#include "opencv/opencv.h"

#include "sequence.h"


/* com.seq is a static struct containing the sequence currently selected by the
 * user from the interface. It may change to be a pointer to any sequence
 * someday, until then, the seqname is NULL when no sequence is loaded and the
 * number of images in the sequence is also negative.
 * com.uniq represents information about an image opened and displayed outside
 * a sequence, for example from the load command, the open menu, or the result
 * of a stacking operation.
 * com.seq.number is used to provide a relationship between a possibly loaded
 * sequence and the single image. A single image can be loaded without
 * unloading the sequence. This information could as well be moved to
 * com.status if com.seq becomes a pointer. Three constants have been declared
 * in siril.h to explicit this relationship: RESULT_IMAGE, UNRELATED_IMAGE and
 * SCALED_IMAGE. They are mostly used to understand what to do to display
 * single images when a sequence is loaded or not.
 */

static void fillSeqAviExport() {
	char width[6], height[6];
	GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));
	GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));

	g_snprintf(width, sizeof(width), "%d", com.seq.rx);
	g_snprintf(height, sizeof(width), "%d", com.seq.ry);
	gtk_entry_set_text(widthEntry, width);
	gtk_entry_set_text(heightEntry, height);
	if (com.seq.type == SEQ_SER) {
		GtkEntry *entryAviFps = GTK_ENTRY(lookup_widget("entryAviFps"));

		if (com.seq.ser_file != NULL) {
			char fps[7];

			if (com.seq.ser_file->fps <= 0.0) {
				g_snprintf(fps, sizeof(fps), "25.000");
			} else {
				g_snprintf(fps, sizeof(fps), "%2.3lf", com.seq.ser_file->fps);
			}
			gtk_entry_set_text(entryAviFps, fps);
		}
	}
}

static sequence *check_seq_one_file(const char* name, gboolean check_for_fitseq);

gboolean populate_seqcombo(gpointer user_data) {
	const gchar *realname = (const gchar*) user_data;
	control_window_switch_to_tab(IMAGE_SEQ);
	GtkComboBoxText *combo_box_text = GTK_COMBO_BOX_TEXT(lookup_widget("sequence_list_combobox"));
	gtk_combo_box_text_remove_all(combo_box_text);
	gchar *rname = g_path_get_basename(realname);
	gtk_combo_box_text_append(combo_box_text, 0, rname);
	g_signal_handlers_block_by_func(GTK_COMBO_BOX(combo_box_text), on_seqproc_entry_changed, NULL);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo_box_text), 0);
	g_signal_handlers_unblock_by_func(GTK_COMBO_BOX(combo_box_text), on_seqproc_entry_changed, NULL);
	g_free(rname);
	return FALSE;
}

/* normalizes sequence name
 * takes a string and
 * - removes the extension if known
 * - appends _ at the end if required and add_underscore is TRUE
 * returns a newly allocated string to be freed with free
 */
char *normalize_seqname(char *name, gboolean add_underscore) {
	char *file_no_ext;
	if (g_str_has_suffix(name, ".seq") || g_str_has_suffix(name, ".fit") || g_str_has_suffix(name, ".fits") ||
	g_str_has_suffix(name, ".fts") || g_str_has_suffix(name, ".ser")) {
		file_no_ext = remove_ext_from_filename(name);
	} else {
		file_no_ext = strdup(name);
	}
	gboolean needs_underscore = add_underscore && !g_str_has_suffix(name, "_");
	gchar *outname = g_strdup_printf("%s%s", file_no_ext, needs_underscore ? "_" : "");
	free(file_no_ext);
	return outname;
}

/* when opening a file outside the main sequence loading system and that file
 * is a sequence (SER/AVI), this function is called to load this sequence. */
int read_single_sequence(char *realname, image_type imagetype) {
	int retval = 0, len;
	gchar *dirname = g_path_get_dirname(realname);
	if (!siril_change_dir(dirname, NULL)) {
		writeinitfile();
		gui_function(set_GUI_CWD, NULL);
	}
	g_free(dirname);

	sequence *new_seq = check_seq_one_file(realname, TRUE); // it's not the real .seq read
	if (!new_seq)
		return 1;
	free_sequence(new_seq, TRUE);

	char *name = strdup(realname);
	const char *ext;
#ifdef HAVE_FFMS2
	const char *film_ext;
#endif
	switch (imagetype) {
	case TYPESER:
		name[strlen(name) - 1] = 'q';
		break;
	case TYPEFITS:
		ext = get_filename_ext(realname);
		assert(ext);
		len = strlen(ext);
		strncpy(name + strlen(name) - len, "seq", len);
		break;
#ifdef HAVE_FFMS2
	case TYPEAVI:
		film_ext = get_filename_ext(realname);
		assert(film_ext);
		len = strlen(film_ext);
		strncpy(name + strlen(name) - len, "seq", len);
		break;
#endif
		default:
			retval = 1;
	}
	gchar *fname = g_path_get_basename(name);
	if (!set_seq(fname) && !com.script) {
		/* if it loads, make it selected and only element in the list of sequences */
		gui_function(populate_seqcombo, realname);
	}
	else retval = 1;
	g_free(fname);
	free(name);
	return retval;
}

/* Find sequences in CWD and create .seq files.
 * In the current working directory, looks for sequences of fits files or files
 * already representing sequences like SER and AVI formats and builds the
 * corresponding sequence files.
 * Called when changing wd with name == NULL or when an explicit root name is
 * given in the GUI or when searching for sequences.
 * force clears the stats in the seqfile.
 */
int check_seq() {
	int curidx, fixed;
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	sequence **sequences;
	int i, nb_seq = 0, max_seq = 10;

	siril_log_color_message(_("Checking sequences in the directory: %s.\n"), "blue", com.wd);

	if (!com.wd) {
		siril_log_message(_("Current working directory is not set, aborting.\n"));
		return 1;
	}
	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL) {
		fprintf (stderr, "check_seq: %s\n", error->message);
		g_clear_error(&error);
		g_free(com.wd);
		com.wd = NULL;
		return 1;
	}

	sequences = malloc(sizeof(sequence *) * max_seq);
	if (!sequences) {
		PRINT_ALLOC_ERR;
		g_dir_close(dir);
		return 1;
	}
	set_progress_bar_data(NULL, PROGRESS_PULSATE);

	while ((file = g_dir_read_name(dir)) != NULL) {
		sequence *new_seq;
		int fnlen = strlen(file);
		if (fnlen < 4) continue;
		const char *ext = get_filename_ext(file);
		if (!ext) continue;

		gboolean is_fz = g_str_has_suffix(ext, ".fz");
		const gchar *com_ext = get_com_ext(is_fz);

		if ((new_seq = check_seq_one_file(file, FALSE))) {
			sequences[nb_seq] = new_seq;
			nb_seq++;
		} else if (!strcasecmp(ext, com_ext + 1)) {
			char *basename = NULL;
			if (!get_index_and_basename(file, &basename, &curidx, &fixed, com_ext)) {
				int current_seq = -1;
				/* search in known sequences if we already have it */
				for (i = 0; i < nb_seq; i++) {
					if (!strcmp(sequences[i]->seqname, basename)) {
						current_seq = i;
						free(basename);
						break;
					}
				}
				/* not found */
				if (current_seq < 0) {
					new_seq = calloc(1, sizeof(sequence));
					initialize_sequence(new_seq, TRUE);
					new_seq->seqname = basename;
					new_seq->beg = INT_MAX;
					new_seq->end = 0;
					new_seq->fixed = fixed;
					new_seq->fz = is_fz;
					sequences[nb_seq] = new_seq;
					current_seq = nb_seq;
					nb_seq++;
					siril_debug_print("Found a sequence (number %d) with base name"
							" \"%s\", looking for first and last indexes.\n",
							nb_seq, basename);
				}
				if (curidx < sequences[current_seq]->beg)
					sequences[current_seq]->beg = curidx;
				if (curidx > sequences[current_seq]->end)
					sequences[current_seq]->end = curidx;
				if (fixed > sequences[current_seq]->fixed)
					sequences[current_seq]->fixed = fixed;
			}
			else if ((new_seq = check_seq_one_file(file, TRUE))) {
				sequences[nb_seq] = new_seq;
				nb_seq++;
			}
		}
		if (nb_seq == max_seq) {
			max_seq *= 2;
			sequence **tmp = realloc(sequences, sizeof(sequence *) * max_seq);
			if (tmp)
				sequences = tmp;
			else {
				PRINT_ALLOC_ERR;
				break;
			}
		}
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);
	g_dir_close(dir);

	if (nb_seq > 0) {
		int retval = 1;
		for (i = 0; i < nb_seq; i++) {
			if (sequences[i]->beg == sequences[i]->end) {
				/* maybe it's a FITSEQ? */
				char *name = malloc(strlen(sequences[i]->seqname) + 20);
				if (!name) {
					PRINT_ALLOC_ERR;
					free_sequence(sequences[i], TRUE);
					continue;
				}
				sequence *new_seq;
				if (get_possible_image_filename(sequences[i],
							sequences[i]->beg, name) &&
						(new_seq = check_seq_one_file(name, TRUE))) {
					free_sequence(sequences[i], TRUE);
					sequences[i] = new_seq;
					free(name);
				}
				else {
					free_sequence(sequences[i], TRUE);
					free(name);
					continue;
				}
			}
			siril_debug_print(_("sequence %d, found: %d to %d\n"),
					i + 1, sequences[i]->beg, sequences[i]->end);
			if (!buildseqfile(sequences[i], 0) && retval) {
				retval = 0;	// at least one succeeded to be created
			}
			free_sequence(sequences[i], TRUE);
		}
		free(sequences);
		return retval;
	}
	free(sequences);
	siril_log_message(_("No sequence found, verify working directory or "
				"change FITS extension in settings (current is %s)\n"), com.pref.ext);
	return 1;	// no sequence found
}

/* Creates a .seq file for regular FITS sequence passed in argument */
static sequence *create_one_regular_seq(const char *seqname) {

	gchar *abs_path = g_canonicalize_filename(seqname, com.wd);
	gchar *search_folder = g_path_get_dirname(abs_path);
	gchar *filename = g_path_get_basename(abs_path);
	g_free(abs_path);

	char *root = remove_ext_from_filename(filename);
	g_free(filename);
	int fixed = 5; // TODO: isn't it defined somewhere else?
	const gchar *ext = get_com_ext(com.pref.comp.fits_enabled);

	GError* error = NULL;
	GDir* dir = g_dir_open(search_folder, 0, &error);
	g_free(search_folder);
	if (error) {
		siril_log_color_message(_("Error opening directory: %s\n"), "red", error->message);
		g_error_free(error);
		free(root);
		return NULL;
    }
	const gchar* pattern = g_strdup_printf("^%s(\\d{%d})\\%s$", 
										root,
										fixed,
										ext);
	GRegex* regex = g_regex_new(pattern, 0, 0, &error);
 
	sequence *new_seq = calloc(1, sizeof(sequence));
	initialize_sequence(new_seq, TRUE);
	new_seq->seqname = root; // move it as we don't need it any more in this fn
	new_seq->beg = INT_MAX;
	new_seq->end = 0;
	new_seq->number = 0;
	new_seq->type = SEQ_REGULAR;
	new_seq->fz = com.pref.comp.fits_enabled;
	new_seq->fixed = fixed;
	int n = 0;

	const gchar *newfits;
	while ((newfits = g_dir_read_name(dir))) {
		GMatchInfo* match_info;
		if (g_regex_match(regex, newfits, 0, &match_info)) {
			const gchar* number_str = g_match_info_fetch(match_info, 1);
			gint number = atoi(number_str);
			if (number < new_seq->beg)
				new_seq->beg = number;
			if (number > new_seq->end)
				new_seq->end = number;
			n++;
		}
		g_match_info_free(match_info);
	}
	if (n < 2) {
		siril_log_color_message(_("Cannot create sequence %s. Need at least 2 frames to be usable in Siril.\n"), "salmon", seqname);
		free_sequence(new_seq, TRUE);
		return NULL;
	}
	if (buildseqfile(new_seq, 1)) {
		free_sequence(new_seq, TRUE);
		new_seq = NULL;
	}
	return new_seq;
}

/* Creates a .seq file for one-file sequence passed in argument */
static sequence *check_seq_one_file(const char* name, gboolean check_for_fitseq) {
	if (!com.wd) {
		siril_log_message(_("Current working directory is not set, aborting.\n"));
		return NULL;
	}
	int fnlen = strlen(name);
	const char *ext = get_filename_ext(name);
	sequence *new_seq = NULL;

	if (!strcasecmp(ext, "ser")) {
		struct ser_struct *ser_file = malloc(sizeof(struct ser_struct));
		ser_init_struct(ser_file);
		int ret = ser_open_file(name, ser_file);
		if (ret) {
			return NULL;
		}

		new_seq = calloc(1, sizeof(sequence));
		initialize_sequence(new_seq, TRUE);
		new_seq->seqname = g_strndup(name, fnlen - 4);
		new_seq->beg = 0;
		new_seq->end = ser_file->frame_count - 1;
		new_seq->number = ser_file->frame_count;
		new_seq->type = SEQ_SER;
		new_seq->ser_file = ser_file;
		siril_debug_print("Found a SER sequence\n");
	}
#ifdef HAVE_FFMS2
	else if (!check_for_film_extensions(ext)) {
		struct film_struct *film_file = malloc(sizeof(struct film_struct));
		if (film_open_file(name, film_file)) {
			free(film_file);
			return NULL;
		}
		new_seq = calloc(1, sizeof(sequence));
		initialize_sequence(new_seq, TRUE);
		int len = strlen(ext);
		new_seq->seqname = g_strndup(name, fnlen-len-1);
		new_seq->beg = 0;
		new_seq->end = film_file->frame_count-1;
		new_seq->number = film_file->frame_count;
		new_seq->type = SEQ_AVI;
		new_seq->film_file = film_file;
		siril_debug_print("Found a AVI sequence\n");
	}
#endif
	else if (check_for_fitseq && TYPEFITS == get_type_for_extension(ext) && fitseq_is_fitseq(name, NULL)) {
		gboolean is_fz = g_str_has_suffix(ext, ".fz");
		const gchar *com_ext = get_com_ext(is_fz);

		/* set the configured extention to the extension of the file, otherwise reading will fail */
		if (strcasecmp(ext, com_ext + 1)) {
			g_free(com.pref.ext);
			com.pref.ext = g_strdup_printf(".%s", ext);
			if (is_fz) com.pref.ext[strlen(com.pref.ext) - 2] = '\0';
		}

		fitseq *fitseq_file = malloc(sizeof(fitseq));
		fitseq_init_struct(fitseq_file);
		if (fitseq_open(name, fitseq_file, READONLY)) {
			free(fitseq_file);
			return NULL;
		}
		new_seq = calloc(1, sizeof(sequence));
		initialize_sequence(new_seq, TRUE);
		new_seq->seqname = g_strndup(name, fnlen - strlen(com_ext));
		new_seq->beg = 0;
		new_seq->end = fitseq_file->frame_count - 1;
		new_seq->number = fitseq_file->frame_count;
		new_seq->type = SEQ_FITSEQ;
		new_seq->fitseq_file = fitseq_file;
		new_seq->fz = is_fz;
		siril_debug_print("Found a FITS sequence\n");
	}

	if (new_seq) {
		if (new_seq->beg != new_seq->end) {
			if (buildseqfile(new_seq, 0)) {
				free_sequence(new_seq, TRUE);
				new_seq = NULL;
			}
		} else if (new_seq->type == SEQ_SER) {
			siril_log_color_message(_("Cannot load SER sequence. Need at least 2 frames to be usable in Siril. "
					"Please convert the SER file into FITS file.\n"), "salmon");
		}
	}
	return new_seq;
}

// get the number of layers and image size for a new sequence
// if load_ref_into_gfit is true, the image is kept into gfit if data loading was
// required, and 1 is returned when it required loading
int seq_check_basic_data(sequence *seq, gboolean load_ref_into_gfit) {
	if (seq->nb_layers == -1 || seq->rx == 0) {	// not init yet, first loading of the sequence
		int image_to_load = sequence_find_refimage(seq);
		fits tmpfit = { 0 }, *fit;

		if (load_ref_into_gfit) {
			clearfits(gfit);
			fit = gfit;
		} else {
			fit = &tmpfit;
			memset(fit, 0, sizeof(fits));
		}

		if (load_ref_into_gfit) {
			if (seq_read_frame(seq, image_to_load, fit, FALSE, -1)) {
				fprintf(stderr, "could not load first image from sequence\n");
				return -1;
			}
		} else {
			if (seq_read_frame_metadata(seq, image_to_load, fit)) {
				fprintf(stderr, "could not load first image from sequence\n");
				return -1;
			}
		}

		/* initialize sequence-related runtime data */
		seq->rx = fit->rx;
		seq->ry = fit->ry;
		seq->bitpix = fit->orig_bitpix;	// for partial read
		fprintf(stdout, "bitpix for the sequence is set as %d\n", seq->bitpix);
		if (seq->nb_layers == -1 || seq->nb_layers != fit->naxes[2]) {
			seq->nb_layers = fit->naxes[2];
			seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
			if (!seq->regparam) {
				PRINT_ALLOC_ERR;
				clearfits(fit);
				return 1;
			}
			seq->distoparam = calloc(seq->nb_layers, sizeof(disto_params));
		}
		seq->needs_saving = TRUE;

		if (load_ref_into_gfit) {
			seq->current = image_to_load;
		} else {
			clearfits(fit);
		}
		return 1;
	}
	return 0;
}

// generates the .seq file 
gboolean create_one_seq(const char *seqname, sequence_type seqtype) {
	sequence *seq = NULL;
	switch (seqtype) {
		case SEQ_REGULAR:
			seq = create_one_regular_seq(seqname);
			break;
		case SEQ_FITSEQ:
		case SEQ_SER:;
			char *root = remove_ext_from_filename(seqname);
			const gchar *ext;
			if (seqtype == SEQ_FITSEQ)
				ext = get_com_ext(com.pref.comp.fits_enabled);
			else
				ext = ".ser";
			const gchar *filename = g_strdup_printf("%s%s", root, ext);
			g_free(root);
			seq = check_seq_one_file(filename, seqtype == SEQ_FITSEQ);
			break;
		default:
			return FALSE;
	}
	if (!seq)
		return FALSE;
	free_sequence(seq, TRUE);
	return TRUE;
}

static gboolean set_seq_gui(gpointer user_data) {
	sequence *seq = (sequence *) user_data;
	init_layers_hi_and_lo_values(MIPSLOHI); // set some hi and lo values in seq->layers,
	set_cutoff_sliders_max_values();// update min and max values for contrast sliders
	set_cutoff_sliders_values();	// update values for contrast sliders for this image
	int layer = set_layers_for_registration();	// set layers in the combo box for registration
	update_seqlist(layer);
	fill_sequence_list(seq, max(layer, 0), FALSE);// display list of files in the sequence on active layer if regdata exists
	set_output_filename_to_sequence_name();
	sliders_mode_set_state(gui.sliders);
	initialize_display_mode();
	update_zoom_label();
	reset_plot(); // reset all plots
	reset_3stars();

	/* initialize image-related runtime data */
	set_display_mode();		// display the display mode in the combo box
	display_filename();		// display filename in gray window
	gui_function(set_precision_switch, NULL); // set precision on screen
	adjust_refimage(seq->current);	// check or uncheck reference image checkbox
	update_prepro_interface(seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ || seq->type == SEQ_SER); // enable or not the preprobutton
	update_reg_interface(FALSE);	// change the registration prereq message
	update_stack_interface(FALSE);	// get stacking info and enable the Go button, already done in set_layers_for_registration
	adjust_reginfo();		// change registration displayed/editable values
	update_gfit_histogram_if_needed();
	adjust_sellabel();
	fillSeqAviExport();	// fill GtkEntry of export box

	/* update menus */
	update_MenuItem(NULL);
	/* update parameters in GUI */
	set_GUI_CAMERA();

	/* redraw and display image */
	gui_function(close_tab, NULL);	//close Green and Blue Tab if a 1-layer sequence is loaded
	gui_function(init_right_tab, NULL);

	redraw(REMAP_ALL);
	drawPlot();
	return FALSE;
}

static void free_cbbt_layers() {
	GtkComboBoxText *cbbt_layers = GTK_COMBO_BOX_TEXT(lookup_widget("comboboxreglayer"));
	gtk_combo_box_text_remove_all(cbbt_layers);
}

/* load a sequence and initializes everything that relates */
gboolean set_seq(gpointer user_data){
	const char *name = (const char*) user_data;
	sequence *seq = NULL;
	char *basename;

	if ((seq = readseqfile(name)) == NULL) {
		fprintf(stderr, "could not load sequence %s\n", name);
		return TRUE;
	}
	free_image_data();
	close_sequence(TRUE);

#ifdef HAVE_FFMS2
	int convert = (int)((com.headless));
	if (!com.headless) {
		if (seq->type == SEQ_AVI) {
			convert = siril_confirm_dialog(_("Deprecated sequence"),
					_("Film sequences are now deprecated in Siril: some features are disabled and others may crash."
							" We strongly encourage you to convert this sequence into a SER file."
							" SER file format is a simple image sequence format, similar to uncompressed films."), _("Convert to SER"));
		}
	} else {
		siril_log_color_message(_("Warning: deprecated sequence. Film sequences are now deprecated "
			"in Siril: some features are disabled and others may crash. Continuing, but "
							"we strongly encourage you to convert this sequence into a SER file."
							"SER file format is a simple image sequence format, similar to uncompressed films.\n"), "salmon");
	}
	if (convert) {
		close_sequence(FALSE);
		convert_single_film_to_ser(seq);
		return FALSE;
	}
#endif
	int retval = seq_check_basic_data(seq, TRUE);
	if (retval == -1) {
		free_sequence(seq, TRUE);
		return TRUE;
	}
	if (retval == 0) {
		int image_to_load = sequence_find_refimage(seq);
		if (seq_read_frame(seq, image_to_load, gfit, FALSE, -1)) {
			siril_log_color_message(_("could not load reference image from sequence\n"), "red");
			free_sequence(seq, TRUE);
			return TRUE;
		}
		seq->current = image_to_load;
	}
	if (seq->type == SEQ_SER)
		ser_display_info(seq->ser_file);

	basename = g_path_get_basename(seq->seqname);
	siril_log_message(_("Sequence loaded: %s (%d->%d)\n"),
			basename, seq->beg, seq->end);
	g_free(basename);

	/* Sequence is stored in com.seq for now */
	memcpy(&com.seq, seq, sizeof(sequence));
	update_gain_from_gfit();

	if (!com.script && !com.headless) {
		execute_idle_and_wait_for_it(set_seq_gui, seq);
	}

	free(seq);
	return FALSE;
}

/* Load image number index from the sequence and display it.
 * if load_it is true, dest is assumed to be gfit
 * TODO: cut that method in two, with an internal func taking a filename and a fits
 */
// This function is OK to have GUI calls in it as it is only ever called from GUI functions
int seq_load_image(sequence *seq, int index, gboolean load_it) {
	gboolean do_refresh_annotations = com.found_object != NULL;
	if (!single_image_is_loaded())
		save_stats_from_fit(gfit, seq, seq->current);
	on_clear_roi(); // Always clear a ROI when changing images
	cleanup_annotation_catalogues(FALSE);
	clear_stars_list(TRUE);
	invalidate_gfit_histogram();
	undo_flush();
	close_single_image();
	clearfits(gfit);
	if (seq->current == SCALED_IMAGE) {
		gfit->rx = seq->rx;
		gfit->ry = seq->ry;
	}
	seq->current = index;

	if (load_it && !com.script) {
		set_cursor_waiting(TRUE);
		if (seq_read_frame(seq, index, gfit, FALSE, -1)) {
			set_cursor_waiting(FALSE);
			return 1;
		}
		set_fwhm_star_as_star_list(seq);// display the fwhm star if possible
		if (gui.sliders != USER) {
			init_layers_hi_and_lo_values(gui.sliders);
			sliders_mode_set_state(gui.sliders);
			set_cutoff_sliders_max_values();// update min and max values for contrast sliders
			set_cutoff_sliders_values();	// update values for contrast sliders for this image
			set_display_mode();		// display the display mode in the combo box
		}
		if (do_refresh_annotations)
			refresh_found_objects();
		redraw(REMAP_ALL);
		if (seq->is_variable)
			clear_previews();
		else
			gui_function(redraw_previews, NULL);		// redraw registration preview areas
		display_filename();		// display filename in gray window
		gui_function(set_precision_switch, NULL); // set precision on screen
		adjust_reginfo();		// change registration displayed/editable values
		update_display_fwhm();
		update_gfit_histogram_if_needed();
		set_cursor_waiting(FALSE);
		reset_3stars();
	}

	gui_function(update_MenuItem, NULL);		// initialize menu gui
	set_GUI_CAMERA();		// update image information
	sequence_list_change_current();
	adjust_refimage(index);	// check or uncheck reference image checkbox

	return 0;
}

// Used by the python interface to ensure seq_load_image is run in the GUI thread
gboolean seq_load_image_in_thread(gpointer user_data) {
	int index = *(int*) user_data;
	seq_load_image(&com.seq, index, TRUE);
	return FALSE;
}

/**
 * Computes size of an opened sequence in bytes for a passed number of frames.
 * For SER or films, it returns the size of the file.
 * For FITS sequences, the reference image's size is used as baseline.
 * Unsupported for internal sequences.
 * @param seq input sequence
 * @param nb_frames number of frames to compute the size of the sequence of
 * @return the size of the sequence in bytes, or -1 if an error happened.
 */
int64_t seq_compute_size(sequence *seq, int nb_frames, data_type depth) {
	int64_t frame_size, size = -1LL;
#ifdef HAVE_FFMS2
	GStatBuf sts;
#endif

	switch(seq->type) {
	case SEQ_SER:
		size = ser_compute_file_size(seq->ser_file, nb_frames);
		break;
	case SEQ_REGULAR:
	case SEQ_FITSEQ:
		frame_size = (int64_t) seq->rx * seq->ry * seq->nb_layers;
		if (depth == DATA_USHORT)
			frame_size *= sizeof(WORD);
		else if (depth == DATA_FLOAT)
			frame_size *= sizeof(float);
		frame_size += FITS_DOUBLE_BLOC_SIZE; // FITS double HDU size
		size = frame_size * nb_frames;
		break;
#ifdef HAVE_FFMS2
	case SEQ_AVI:
		if (g_stat(seq->film_file->filename, &sts) == 0) {
			// this is a close approximation
			frame_size = sts.st_size / seq->film_file->frame_count;
			size = nb_frames * frame_size;
		}
		break;
#endif
	default:
		fprintf(stderr, "Failure: computing sequence size on internal sequence is unsupported\n");
	}
	return size;
}

/**
 * Check if a sequence with a basename 'name' or a full name 'name' already exists
 * @param name either the base name of the sequence or its full name in case of
 * single file sequence
 * @param name_is_base TRUE is the name is a base name
 * @return TRUE if the name already exists, FALSE otherwise
 */
gboolean check_if_seq_exist(gchar *name, gboolean name_is_base) {
	gchar *path;
	if (name_is_base) {
		gchar *path_;

		gchar *seq = g_strdup_printf("%s.seq", name);
		gchar *seq_ = g_strdup_printf("%s_.seq", name);
		path = g_build_filename(com.wd, seq, NULL);
		path_ = g_build_filename(com.wd, seq_, NULL);
		g_free(seq);
		g_free(seq_);
		gboolean retval = is_readable_file(path);
		if (!retval) {
			retval = is_readable_file(path_);
		}
		g_free(path);
		g_free(path_);
		return retval;
	} else {
		path = g_build_filename(com.wd, name, NULL);
		gboolean retval = is_readable_file(path);
		g_free(path);
		return retval;
	}
}

/*****************************************************************************
 *              SEQUENCE FUNCTIONS FOR NON-OPENED SEQUENCES                  *
 * **************************************************************************/

/* Get the filename of an image in a sequence.
 * Return value is the same as the name_buf argument, which must be
 * pre-allocated to at least 256 characters. If sequence has no file names, a
 * description like image "42 from awesome_mars.ser" is made. */
char *seq_get_image_filename(sequence *seq, int index, char *name_buf) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return fit_sequence_get_image_filename(seq, index, name_buf, TRUE);
		case SEQ_SER:
			if (!name_buf || index < 0 || index > seq->end) {
				return NULL;
			}
			snprintf(name_buf, 255, "%s_%d.ser", seq->seqname,  index);
			return name_buf;
		case SEQ_FITSEQ:
			if (!name_buf || index < 0 || index > seq->end) {
				return NULL;
			}
			snprintf(name_buf, 255, "%s_%d%s", seq->seqname,  index, get_com_ext(seq->fz));
			return name_buf;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			if (!name_buf || index < 0 || index > seq->end) {
				return NULL;
			}
			snprintf(name_buf, 255, "%s_%d", seq->seqname, index);
			return name_buf;
#endif
		case SEQ_INTERNAL:
			snprintf(name_buf, 255, "%s_%d", seq->seqname, index);
			return name_buf;
	}
	return NULL;
}

gchar *get_image_filename_no_ext(sequence *seq, int idx) {
        if (!seq) {
                // check loaded image
                if (com.uniq)
                        return g_strdup(com.uniq->filename);
                seq = &com.seq;
                idx = com.seq.current;
        }
        char root[256];
        if (!fit_sequence_get_image_filename(seq, idx, root, FALSE)) {
                return NULL;
        }
	// TO COMPLETE FOR OTHER TYPES OF SEQUENCE
        return g_strdup(root);
}

/* Read an entire image from a sequence, inside a pre-allocated fits.
 * Opens the file, reads data, closes the file.
 */
int seq_read_frame(sequence *seq, int index, fits *dest, gboolean force_float, int thread_id) {
	char filename[256];
	assert(index < seq->number);
	switch (seq->type) {
		case SEQ_REGULAR:
			fit_sequence_get_image_filename_checkext(seq, index, filename);
			int ret = readfits(filename, dest, NULL, force_float);
			if (ret) {
				char *base = remove_all_ext_from_filename(filename);
				ret = readfits(base, dest, NULL, force_float);
				free(base);
				if (ret) {
					siril_log_message(_("Could not load image %d from sequence %s\n"),
							index, seq->seqname);
					return 1;
				}
			}
			break;
		case SEQ_SER:
			assert(seq->ser_file);
			if (ser_read_frame(seq->ser_file, index, dest, force_float, com.pref.debayer.open_debayer)) {
				siril_log_message(_("Could not load frame %d from SER sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			break;
		case SEQ_FITSEQ:
			assert(seq->fitseq_file);
			if (fitseq_read_frame(seq->fitseq_file, index, dest, force_float, thread_id)) {
				siril_log_message(_("Could not load frame %d from FITS sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			break;

#ifdef HAVE_FFMS2
		case SEQ_AVI:
			assert(seq->film_file);
			if (film_read_frame(seq->film_file, index, dest)) {
				siril_log_message(_("Could not load frame %d from AVI sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			// should dest->maxi be set to 255 here?
			break;
#endif
		case SEQ_INTERNAL:
			assert(seq->internal_fits);
			// copyfits copies ICC profile so internal sequences do retain ICC profiles
			copyfits(seq->internal_fits[index], dest, CP_FORMAT, -1);
			copy_fits_metadata(seq->internal_fits[index], dest);
			if (seq->internal_fits[index]->type == DATA_FLOAT) {
				dest->fdata = seq->internal_fits[index]->fdata;
				dest->fpdata[0] = seq->internal_fits[index]->fpdata[0];
				dest->fpdata[1] = seq->internal_fits[index]->fpdata[1];
				dest->fpdata[2] = seq->internal_fits[index]->fpdata[2];
			}
			else if (seq->internal_fits[index]->type == DATA_USHORT) {
				dest->data = seq->internal_fits[index]->data;
				dest->pdata[0] = seq->internal_fits[index]->pdata[0];
				dest->pdata[1] = seq->internal_fits[index]->pdata[1];
				dest->pdata[2] = seq->internal_fits[index]->pdata[2];
			}
			else return 1;
			break;
	}
	if (seq->nb_layers > 0 &&  seq->nb_layers != dest->naxes[2]) {
		siril_log_color_message(_("Image #%d: number of layers (%d) is not consistent with sequence (%d), aborting\n"), "red",
			index, dest->naxes[2], seq->nb_layers);
		return 1;
	}

	full_stats_invalidation_from_fit(dest);
	copy_seq_stats_to_fit(seq, index, dest);
	seq->imgparam[index].rx = dest->rx;
	seq->imgparam[index].ry = dest->ry;
	if (seq->rx != 0 && seq->ry != 0 && (dest->rx != seq->rx || dest->ry != seq->ry)) {
		siril_debug_print("sequence detected as containing images of different sizes\n");
		seq->is_variable = TRUE;
	}
	return 0;
}

/* same as seq_read_frame above, but creates an image the size of the selection
 * rectangle only. layer is set to the layer number in the read partial frame.
 * The partial image result is only one-channel deep, so it cannot be used to
 * have a partial RGB image. */
int seq_read_frame_part(sequence *seq, int layer, int index, fits *dest, const rectangle *area, gboolean do_photometry, int thread_id) {
	char filename[256];
#ifdef HAVE_FFMS2
	fits tmp_fit;
#endif
	switch (seq->type) {
		case SEQ_REGULAR:
			fit_sequence_get_image_filename_checkext(seq, index, filename);
			if (readfits_partial(filename, layer, dest, area, do_photometry)) {
				siril_log_message(_("Could not load partial image %d from sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			break;
		case SEQ_SER:
			assert(seq->ser_file);
			if (ser_read_opened_partial_fits(seq->ser_file, layer, index, dest, area)) {
				siril_log_message(_("Could not load frame %d from SER sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			break;
		case SEQ_FITSEQ:
			assert(seq->fitseq_file);
			if (fitseq_read_partial_fits(seq->fitseq_file, layer, index, dest, area, do_photometry, thread_id)) {
				siril_log_message(_("Could not load partial image %d from sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			break;

#ifdef HAVE_FFMS2
		case SEQ_AVI:
			assert(seq->film_file);
			memset(&tmp_fit, 0, sizeof(fits));
			if (film_read_frame(seq->film_file, index, &tmp_fit)) {
				siril_log_message(_("Could not load frame %d from AVI sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			extract_region_from_fits(&tmp_fit, layer, dest, area);
			clearfits(&tmp_fit);
			break;
#endif
		case SEQ_INTERNAL:
			assert(seq->internal_fits);
			extract_region_from_fits(seq->internal_fits[index], 0, dest, area);
			break;
	}
	return 0;
}

// not thread-safe
// gets image naxes and bitpix
int seq_read_frame_metadata(sequence *seq, int index, fits *dest) {
	assert(index < seq->number);
	char filename[256];
	switch (seq->type) {
		case SEQ_REGULAR:
			fit_sequence_get_image_filename_checkext(seq, index, filename);
			if (read_fits_metadata_from_path(filename, dest)) {
				siril_log_message(_("Could not load image %d from sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			break;
		case SEQ_SER:
			assert(seq->ser_file);
			if (ser_metadata_as_fits(seq->ser_file, dest)) {
				siril_log_message(_("Could not load frame %d from SER sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			break;
		case SEQ_FITSEQ:
			assert(seq->fitseq_file);
			if (seq->fitseq_file->thread_fptr) {
#ifdef _OPENMP
				int thread_id = omp_get_thread_num();
				dest->fptr = seq->fitseq_file->thread_fptr[thread_id];
				int status = 0;
				fits_movabs_hdu(dest->fptr, seq->fitseq_file->hdu_index[index], NULL, &status);
				if (status) {
					siril_log_message(_("Could not seek frame %d from FITS sequence %s. Error status: %d\n"),
							index, seq->seqname, status);
					return 1;
				}
				if (read_fits_metadata(dest)) {
					siril_log_message(_("Could not load frame %d from FITS sequence %s\n"),
							index, seq->seqname);
					return 1;
				}
#else
				return 1;
#endif
			} else {
				dest->fptr = seq->fitseq_file->fptr;
				if (fitseq_set_current_frame(seq->fitseq_file, index) ||
						read_fits_metadata(dest)) {
					siril_log_message(_("Could not load frame %d from FITS sequence %s\n"),
							index, seq->seqname);
					return 1;
				}
			}
			break;

#ifdef HAVE_FFMS2
		case SEQ_AVI:
			assert(seq->film_file);
			// TODO: do a metadata-only read in films
			if (film_read_frame(seq->film_file, index, dest)) {
				siril_log_message(_("Could not load frame %d from AVI sequence %s\n"),
						index, seq->seqname);
				return 1;
			}
			// should dest->maxi be set to 255 here?
			break;
#endif
		case SEQ_INTERNAL:
			assert(seq->internal_fits);
			copyfits(seq->internal_fits[index], dest, CP_FORMAT, -1);
			break;
	}
	seq->imgparam[index].rx = dest->rx;
	seq->imgparam[index].ry = dest->ry;
	return 0;
}

/*****************************************************************************
 *                 SEQUENCE FUNCTIONS FOR OPENED SEQUENCES                   *
 * **************************************************************************/

/* locks cannot be probed to see if they are init or not, so we have to keep
 * all of them in the same state, which is initialized if the array is non-nul. */
static int _allocate_sequence_locks(sequence *seq) {
#ifdef _OPENMP
	if (!seq->fd_lock) {
		int i;
		seq->fd_lock = malloc(seq->number * sizeof(omp_lock_t));
		if (!seq->fd_lock) {
			PRINT_ALLOC_ERR;
			return 1;
		}

		for (i=0; i<seq->number; i++)
			omp_init_lock(&seq->fd_lock[i]);
	}
#endif
	return 0;
}

/* open image for future intensive operations (read only) */
int seq_open_image(sequence *seq, int index) {
	int status = 0;
	char filename[256];
	switch (seq->type) {
		case SEQ_REGULAR:
			if (!seq->fptr) {
				seq->fptr = calloc(seq->number, sizeof(fitsfile *));
				if (!seq->fptr) {
					PRINT_ALLOC_ERR;
					return 1;
				}
			}
			if (_allocate_sequence_locks(seq))
				return 1;

			fit_sequence_get_image_filename_checkext(seq, index, filename);
			siril_fits_open_diskfile_img(&seq->fptr[index], filename, READONLY, &status);
			if (status) {
				fits_report_error(stderr, status);
				return status;
			}
			/* should we check image parameters here? such as bitpix or naxis */
			break;
		case SEQ_SER:
			assert(seq->ser_file->file == NULL);
			break;
		case SEQ_FITSEQ:
			assert(seq->fitseq_file->fptr == NULL);
			break;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			siril_log_message(_("This operation is not supported on AVI sequences (seq_open_image)\n"));
			return 1;
#endif
		case SEQ_INTERNAL:
			siril_log_message(_("This operation is not supported on internal sequences (seq_open_image)\n"));
			return 1;
	}
	return 0;
}

/* close opened images, only useful for regular FITS sequences */
void seq_close_image(sequence *seq, int index) {
	int status = 0;
	switch (seq->type) {
		case SEQ_REGULAR:
			if (seq->fptr && seq->fptr[index]) {
				fits_close_file(seq->fptr[index], &status);
				seq->fptr[index] = NULL;
			}
			break;
		default:
			break;
	}
}

/* read a region in a layer of an opened file from a sequence.
 * The buffer must have been allocated to the size of the area, with type of
 * float if seq->bitpix is FLOAT_IMG, with 16-bit type otherwise
 * Used only by median and rejection stacking.
 */
int seq_opened_read_region(sequence *seq, int layer, int index, void *buffer, const rectangle *area, int thread_id) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return read_opened_fits_partial(seq, layer, index, buffer, area);
		case SEQ_SER:
			return ser_read_opened_partial(seq->ser_file, layer, index, buffer, area);
		case SEQ_FITSEQ:
			return fitseq_read_partial(seq->fitseq_file, layer, index, buffer, area, thread_id);
		default:
			break;
	}
	return 0;
}


/*****************************************************************************
 *                         SEQUENCE DATA MANAGEMENT                          *
 * **************************************************************************/

/* if FWHM was calculated on the sequence, a minimisation exists for all
 * images, and when switching to a new image, it should be set as the only item
 * in the star list, in order to be displayed.
 * A special care is required in PSF_list.c:clear_stars_list(), to not free this data. */
static void set_fwhm_star_as_star_list_with_layer(sequence *seq, int layer) {
	assert(seq->regparam);
	/* we chose here the first layer that has been allocated, which doesn't
	 * mean it contains data for all images. Handle with care. */
	if (seq->regparam && layer >= 0 && layer < seq->nb_layers
			&& seq->regparam[layer] && seq->current >= 0
			&& seq->regparam[layer][seq->current].fwhm_data && !com.stars) {
		com.stars = new_fitted_stars(1);
		com.stars[0] = seq->regparam[layer][seq->current].fwhm_data;
		com.stars[1] = NULL;
		// this is freed in PSF_list.c:clear_stars_list()
		com.star_is_seqdata = TRUE;
	}
}

// cannot be called in the worker thread
void set_fwhm_star_as_star_list(sequence *seq) {
	int layer = get_registration_layer(seq);
	set_fwhm_star_as_star_list_with_layer(seq, layer);
}

/* Rebuilds the file name of an image in a sequence.
 * The file name is stored in name_buffer, which must be allocated 256 bytes
 * The index is the index in the sequence, not the number appearing in the file name
 * Return value: NULL on error, name_buffer on success.
 */
char *fit_sequence_get_image_filename(sequence *seq, int index, char *name_buffer, gboolean add_fits_ext) {
	char format[20];
	const gchar *com_ext = get_com_ext(seq->fz);

	if (index < 0 || index > seq->number || name_buffer == NULL)
		return NULL;
	if (seq->fixed <= 1) {
		sprintf(format, "%%s%%d");
	} else {
		sprintf(format, "%%s%%.%dd", seq->fixed);
	}
	if (add_fits_ext)
		strncat(format, com_ext, 19); // static char* length 20, leave 1 char for the NULL termination
	snprintf(name_buffer, 255, format, seq->seqname, seq->imgparam[index].filenum);
	name_buffer[255] = '\0';

	return name_buffer;
}

char *fit_sequence_get_image_filename_checkext(sequence *seq, int index, char *name_buffer) {
	char format[20];
	const gchar *default_ext;

	if (index < 0 || index > seq->number || name_buffer == NULL)
		return NULL;

	// Use cached extension if available, otherwise get default
	if (seq->cached_ext != NULL) {
		default_ext = seq->cached_ext;
	} else {
		default_ext = get_com_ext(seq->fz);
	}

	if (seq->fixed <= 1) {
		sprintf(format, "%%s%%d");
	} else {
		sprintf(format, "%%s%%.%dd", seq->fixed);
	}

	// List of extensions to check
	const char *extensions[] = {".fit", ".fits", ".fts", ".FIT", ".FITS", ".FTS"};
	int num_extensions = 6;

	// Build base filename without extension
	snprintf(name_buffer, 255, format, seq->seqname, seq->imgparam[index].filenum);

	// First, try default_ext (com_ext or cached_ext)
	char test_path[256];
	snprintf(test_path, 255, "%s%s", name_buffer, default_ext);

	if (g_file_test(test_path, G_FILE_TEST_EXISTS)) {
		strncpy(name_buffer, test_path, 255);
		name_buffer[255] = '\0';

		// Cache the extension if not already cached
		if (seq->cached_ext == NULL) {
			seq->cached_ext = strdup(default_ext);
		}

		return name_buffer;
	}

	// If default_ext didn't match, try other extensions
	for (int i = 0; i < num_extensions; i++) {
		// Skip if this extension matches default_ext (already checked)
		if (seq->fz) {
			// For compressed files, check if extension + .fz matches default_ext
			char temp_ext[20];
			snprintf(temp_ext, 19, "%s%s", extensions[i], ".fz");
			if (strcmp(temp_ext, default_ext) == 0) continue;
		} else {
			// For uncompressed files, check if extension matches default_ext
			if (strcmp(extensions[i], default_ext) == 0) continue;
		}

		snprintf(test_path, 255, "%s%s", name_buffer, extensions[i]);

		// If seq->fz is TRUE, add .fz suffix
		if (seq->fz) {
			strncat(test_path, ".fz", 255 - strlen(test_path) - 1);
		}

		// Check if file exists
		if (g_file_test(test_path, G_FILE_TEST_EXISTS)) {
			strncpy(name_buffer, test_path, 255);
			name_buffer[255] = '\0';

			// Cache the found extension if not already cached
			if (seq->cached_ext == NULL) {
				if (seq->fz) {
					// Cache extension with .fz suffix
					char full_ext[20];
					snprintf(full_ext, 19, "%s.fz", extensions[i]);
					seq->cached_ext = strdup(full_ext);
				} else {
					seq->cached_ext = strdup(extensions[i]);
				}
			}

			return name_buffer;
		}
	}

	// If no file found, fall back to default behavior (use default_ext)
	snprintf(name_buffer, 255, format, seq->seqname, seq->imgparam[index].filenum);
	strncat(name_buffer, default_ext, 255 - strlen(name_buffer) - 1);

	name_buffer[255] = '\0';
	return name_buffer;
}

char *fit_sequence_get_image_filename_prefixed(sequence *seq, const char *prefix, int index) {
	char format[16];
	const gchar *com_ext = get_com_ext(seq->fz);
	gchar *basename = g_path_get_basename(seq->seqname);
	GString *str = g_string_sized_new(70);

	sprintf(format, "%%s%%s%%0%dd%%s", seq->fixed);
	g_string_printf(str, format, prefix, basename, seq->imgparam[index].filenum, com_ext);
	g_free(basename);

	return g_string_free(str, FALSE);
}

/* Returns a filename for an image that could be in a sequence, but the sequence structure
 * has not been fully initialized yet. Only beg, end, fixed and seqname are used.
 */
char *get_possible_image_filename(sequence *seq, int image_number, char *name_buffer) {
	char format[20];
	const gchar *com_ext = get_com_ext(seq->fz);

	if (image_number < seq->beg || image_number > seq->end || name_buffer == NULL)
		return NULL;
	if (seq->fixed <= 1){
		sprintf(format, "%%s%%d%s", com_ext);
	} else {
		sprintf(format, "%%s%%.%dd%s", seq->fixed, com_ext);
	}

	sprintf(name_buffer, format, seq->seqname, image_number);
	return name_buffer;
}

/* splits a filename in a base name and an index number, if the file name ends with .fit
 * it also computes the fixed length if there are zeros in the index */
int	get_index_and_basename(const char *filename, char **basename, int *index, int *fixed, const gchar *com_ext){
	char *buffer;
	int i, fnlen, first_zero, digit_idx;

	*index = -1;		// error values
	*fixed = 0;
	first_zero = -1;
	*basename = NULL;
	fnlen = strlen(filename);
	if (fnlen < strlen(com_ext) + 2) return -1;
	if (!g_str_has_suffix(filename, com_ext)) return -1;
	i = fnlen - strlen(com_ext) - 1;
	if (!isdigit(filename[i])) return -1;
	digit_idx = i;

	buffer = strdup(filename);
	buffer[fnlen - strlen(com_ext)] = '\0';		// for g_ascii_strtoll()
	do {
		if (buffer[i] == '0' && first_zero < 0)
			first_zero = i;
		if (buffer[i] != '0' && first_zero > 0)
			first_zero = -1;
		i--;
	} while (i >= 0 && isdigit(buffer[i]));
	i++;
	if (i == 0) {
		free(buffer);
		return -1;	// no base name, only number
	}
	if (first_zero >= 0)
		*fixed = digit_idx - i + 1;
	//else *fixed = 0;
	*index = g_ascii_strtoll(buffer+i, NULL, 10);
	if (*basename == NULL) {	// don't copy it if we already have it
		*basename = malloc(i * sizeof(char) + 1);
		strncpy(*basename, buffer, i);
		(*basename)[i] = '\0';
	}
	//fprintf(stdout, "from filename %s, base name is %s, index is %d\n", filename, *basename, *index);
	free(buffer);
	return 0;
}

void remove_prefixed_sequence_files(sequence *seq, const char *prefix) {
	int i, len;
	gchar *basename, *seqname;
	if (!prefix || prefix[0] == '\0')
		return;
	basename = g_path_get_basename(seq->seqname);
	len = strlen(basename) + 5 + strlen(prefix);
	seqname = malloc(len);
	g_snprintf(seqname, len, "%s%s.seq", prefix, basename);
	siril_debug_print("Removing %s\n", seqname);
	if (g_unlink(seqname))
		siril_debug_print("g_unlink() failed\n"); // removing the seqfile
	free(seqname);
	g_free(basename);

	switch (seq->type) {
	default:
	case SEQ_REGULAR:
		for (i = 0; i < seq->number; i++) {
			// TODO: use com.cache_upscaled and the current sequence
			// filter to leave the images to be up-scaled.
			char *filename = fit_sequence_get_image_filename_prefixed(
					seq, prefix, i);
			siril_debug_print("Removing %s\n", filename);
			if (g_unlink(filename))
				siril_debug_print("g_unlink() failed\n");
			free(filename);
		}
		break;
	case SEQ_SER:
	case SEQ_FITSEQ:
		if (seq->type == SEQ_SER)
			basename = seq->ser_file->filename;
		else basename = seq->fitseq_file->filename;
		len = strlen(basename) + strlen(prefix) + 1;
		seqname = malloc(len);
		g_snprintf(seqname, len, "%s%s", prefix, basename);
		siril_debug_print("Removing %s\n", seqname);
		if (g_unlink(seqname))
			siril_debug_print("g_unlink() failed\n");
		free(seqname);
		break;
	}
}

void remove_prefixed_star_files(sequence *seq, const char *prefix) {
	for (int i = 0; i < seq->number; i++) {
		const gchar *star_filename = get_sequence_cache_filename(seq, i, "cache", "lst", prefix);
		siril_debug_print("Removing %s\n", star_filename);
		if (g_unlink(star_filename))
			siril_debug_print("g_unlink() failed\n");
	}
}

void remove_prefixed_drizzle_files(sequence *seq, const char *prefix) {
	for (int i = 0; i < seq->number; i++) {
		const gchar *drizzle_filename = get_sequence_cache_filename(seq, i, "drizztmp", "fit", prefix);
		siril_debug_print("Removing %s\n", drizzle_filename);
		if (g_unlink(drizzle_filename))
			siril_debug_print("g_unlink() failed\n");
	}
}

/* sets default values for the sequence */
void initialize_sequence(sequence *seq, gboolean is_zeroed) {
	int i;
	if (!is_zeroed) {
		memset(seq, 0, sizeof(sequence));
	}
	seq->nb_layers = -1;		// uninit value
	seq->reference_image = -1;	// uninit value
	seq->reference_star = -1;	// uninit value
	seq->type = SEQ_REGULAR;
	for (i=0; i<PREVIEW_NB; i++) {
		seq->previewX[i] = -1;
		seq->previewY[i] = -1;
	}
}

/* call this to close a sequence. Second arg must be FALSE for com.seq
 * WARNING: the data is not reset to NULL, if seq is to be reused,
 * initialize_sequence() must be called on it right after free_sequence()
 * (= do it for com.seq) */
void free_sequence(sequence *seq, gboolean free_seq_too) {
	if (seq == NULL) return;
	siril_debug_print("free_sequence(%s)\n", seq->seqname ? seq->seqname : "null name");
	int layer, j;

	// free regparam
	if (seq->nb_layers > 0 && seq->regparam) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->regparam[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->regparam[layer][j].fwhm_data &&
							(!seq->photometry[0] ||
							 seq->regparam[layer][j].fwhm_data != seq->photometry[0][j])) // avoid double free
						free_psf(seq->regparam[layer][j].fwhm_data);
				}
				free(seq->regparam[layer]);
			}
		}
		free(seq->regparam);
	}
	if (seq->distoparam) {
		g_free(seq->distoparam->filename);
	}
	// free stats
	if (seq->nb_layers > 0 && seq->stats) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->stats[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->stats[layer][j])
						free_stats(seq->stats[layer][j]);
				}
				free(seq->stats[layer]);
			}
		}
		free(seq->stats);
	}
	// free backup regparam
	if (seq->nb_layers > 0 && seq->regparam_bkp) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->regparam_bkp[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->regparam_bkp[layer][j].fwhm_data &&
							(!seq->photometry[0] ||
							 seq->regparam_bkp[layer][j].fwhm_data != seq->photometry[0][j])) // avoid double free
						free(seq->regparam_bkp[layer][j].fwhm_data);
				}
				free(seq->regparam_bkp[layer]);
			}
		}
		free(seq->regparam_bkp);
	}
	// free backup stats
	if (seq->nb_layers > 0 && seq->stats_bkp) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			if (seq->stats_bkp[layer]) {
				for (j = 0; j < seq->number; j++) {
					if (seq->stats_bkp[layer][j])
						free_stats(seq->stats_bkp[layer][j]);
				}
				free(seq->stats_bkp[layer]);
			}
		}
		free(seq->stats_bkp);
	}

	for (j = 0; j < seq->number; j++) {
		if (seq->fptr && seq->fptr[j]) {
			int status = 0;
			fits_close_file(seq->fptr[j], &status);
		}
		if (seq->imgparam) {
			if (seq->imgparam[j].date_obs) {
				g_date_time_unref(seq->imgparam[j].date_obs);
			}
		}
	}
	if (seq->seqname)	free(seq->seqname);
	if (seq->imgparam)	free(seq->imgparam);
	if (seq->fptr)		free(seq->fptr);

#ifdef _OPENMP
	if (seq->fd_lock) {
		for (j=0; j<seq->number; j++) {
			omp_destroy_lock(&seq->fd_lock[j]);
		}
		free(seq->fd_lock);
	}
#endif

	if (seq->ser_file) {
		ser_close_file(seq->ser_file);	// frees the data too
		free(seq->ser_file);
	}
#ifdef HAVE_FFMS2
	if (seq->film_file) {
		film_close_file(seq->film_file);	// frees the data too
		free(seq->film_file);
	}
#endif
	if (seq->fitseq_file) {
		fitseq_close_file(seq->fitseq_file);
		free(seq->fitseq_file);
	}

	if (seq->internal_fits) {
		// Compositing still uses references to the images in the sequence
		//for (j=0; j<seq->number; j++)
		//	clearfits(seq->internal_fits[j]);
		free(seq->internal_fits);
	}
	/* Here this is a bit tricky. An internal sequence is a single image. So some
	 * processes like RGB alignment could free sequences and load it again: we need
	 * to keep undo history.
	 * In the case of a standard sequence, loading a new sequence MUST remove all
	 * undo history.
	 */
	if (seq->type != SEQ_INTERNAL)
		undo_flush();

	for (j = 0; j < MAX_SEQPSF && seq->photometry[j]; j++) {
		free_photometry_set(seq, j);
	}
	free(seq->cached_ext);
	if (free_seq_too)	free(seq);
}

void free_photometry_set(sequence *seq, int set) {
	for (int j = 0; j < seq->number; j++) {
		if (seq->photometry[set][j])
			free_psf(seq->photometry[set][j]);
	}
	free(seq->photometry[set]);
	seq->photometry[set] = NULL;
}

gboolean sequence_is_loaded() {
	return (com.seq.seqname != NULL && com.seq.imgparam != NULL);
}

gboolean check_seq_is_comseq(const sequence *seq) {
	if (sequence_is_loaded() && !g_strcmp0(com.seq.seqname, seq->seqname))
		return TRUE;
	return FALSE;
}

gboolean check_seq_is_variable(const sequence *seq) {
	if(!seq || !seq->imgparam)
		return FALSE;
	int rx = seq->imgparam[0].rx;
	int ry = seq->imgparam[0].ry;
	for (int i = 1; i < seq->number; i++) {
		if (seq->imgparam[i].rx != rx || seq->imgparam[i].ry != ry)
			return TRUE;
	}
	return FALSE;
}


gboolean close_sequence_idle(gpointer data) {
	fprintf(stdout, "closing sequence idle\n");
	free_cbbt_layers();
	clear_sequence_list();
	clear_stars_list(TRUE);
	reset_3stars();
	clear_previews();
	free_reference_image();
	update_stack_interface(TRUE);
	adjust_sellabel();
	update_seqlist(-1);
	free_cut_args(&gui.cut);
	initialize_cut_struct(&gui.cut);
	/* unselect the sequence in the sequence list if it's not the one
	 * being loaded */
	if (!data) { // loading_sequence_from_combo
		GtkComboBox *seqcombo = GTK_COMBO_BOX(lookup_widget("sequence_list_combobox"));
		gtk_combo_box_set_active(seqcombo, -1);
	}
	return FALSE;
}

static void close_sequence_gui(gboolean loading_sequence_from_combo) {
	if (com.script || com.python_command)
		execute_idle_and_wait_for_it(close_sequence_idle,
				GINT_TO_POINTER(loading_sequence_from_combo));
	else close_sequence_idle(GINT_TO_POINTER(loading_sequence_from_combo));
}

/* Close the com.seq sequence */
void close_sequence(int loading_sequence_from_combo) {
	if (sequence_is_loaded()) {
		fprintf(stdout, "MODE: closing sequence\n");
		siril_log_message(_("Closing sequence %s\n"), com.seq.seqname);
		if (!single_image_is_loaded())
			save_stats_from_fit(gfit, &com.seq, com.seq.current);
		if (com.seq.needs_saving)
			writeseqfile(&com.seq);

		free_sequence(&com.seq, FALSE);
		initialize_sequence(&com.seq, FALSE);

		if (!com.headless)
			close_sequence_gui(loading_sequence_from_combo);
	}
}

/* if no reference image has been set, return the index of an image that is
 * selected in the sequence, the best of the first registration data found if
 * any, the first otherwise */
int sequence_find_refimage(sequence *seq) {
	if (seq->reference_image != -1)
		return seq->reference_image;
	if (seq->type == SEQ_INTERNAL)
		return 1; // green channel
	int layer, image, best = -1;
	for (layer = 0; layer < seq->nb_layers; layer++) {
		if (seq->regparam && seq->regparam[layer]) {
			gboolean use_fwhm;
			double best_val;
			if (seq->regparam[layer][0].fwhm > 0.0) {
				use_fwhm = TRUE;
				best_val = 1000000.0;
			} else if (seq->regparam[layer][0].quality > 0.0) {
				use_fwhm = FALSE;
				best_val = 0.0;
			}
			else continue;

			for (image = 0; image < seq->number; image++) {
				if (!seq->imgparam[image].incl)
					continue;
				if (use_fwhm && seq->regparam[layer][image].fwhm > 0 &&
						seq->regparam[layer][image].fwhm < best_val) {
					best_val = seq->regparam[layer][image].fwhm;
					best = image;
				} else if (seq->regparam[layer][image].quality > 0 &&
						seq->regparam[layer][image].quality > best_val) {
					best_val = seq->regparam[layer][image].quality;
					best = image;
				}
			}
		}
	}

	if (best == -1 && seq->selnum > 0) {
		for (image = 0; image < seq->number; image++) {
			if (seq->imgparam[image].incl) {
				best = image;
				break;
			}
		}
	}

	if (best == -1) best = 0;	// the first anyway if no regdata and no selected

	return best;
}

/* requires seq->nb_layers and seq->number to be already set */
void check_or_allocate_regparam(sequence *seq, int layer) {
	assert(layer < seq->nb_layers);
	if (!seq->regparam && seq->nb_layers > 0) {
		seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
		seq->distoparam = calloc(seq->nb_layers, sizeof(disto_params));
	}
	if (seq->regparam && !seq->regparam[layer] && seq->number > 0) {
		seq->regparam[layer] = calloc(seq->number, sizeof(regdata));
		for (int i = 0; i < seq->number; i++) {
			cvGetEye(&seq->regparam[layer][i].H);
		}
	}
}

/* assign shift values for registration data of a sequence, depending on its type and sign */
void set_shifts(sequence *seq, int frame, int layer, double shiftx, double shifty, gboolean data_is_top_down) {
	if (seq->regparam[layer]) {
		seq->regparam[layer][frame].H = H_from_translation(shiftx,
				data_is_top_down ? -shifty : shifty);
	}
}

void cum_shifts(regdata *regparam, int frame, double shiftx, double shifty) {
	if (regparam) {
		Homography Hshift = { 0 }, res = { 0 };
		cvGetEye(&Hshift);
		Hshift.h02 = shiftx;
		Hshift.h12 = shifty;
		cvMultH(Hshift, regparam[frame].H, &res);
		regparam[frame].H = res;
	}
}

// Checks that the number of degrees of freedoms is not more than shift
// returns FALSE if not
gboolean test_regdata_is_valid_and_shift(sequence *seq, int reglayer) {
	if (!layer_has_registration(seq, reglayer)) return TRUE;
	transformation_type regmin, regmax;
	guess_transform_from_seq(seq, reglayer, &regmin, &regmax, FALSE);
	if (regmax > SHIFT_TRANSFORMATION)
		return FALSE;
	else if (regmax == SHIFT_TRANSFORMATION) {
		siril_log_color_message(_("This operation will use registration data of layer %d\n"),
				"salmon", reglayer);
	}
	return TRUE;
}

/* internal sequence are a set of 1-layer images already loaded elsewhere, and
 * directly referenced as fits *.
 * This is used in LRGV composition.
 * The returned sequence does not contain any reference to files, and thus has
 * to be populated with internal_sequence_set() */
sequence *create_internal_sequence(int size) {
	int i;
	sequence *seq = calloc(1, sizeof(sequence));
	initialize_sequence(seq, TRUE);
	seq->type = SEQ_INTERNAL;
	seq->number = size;
	seq->selnum = size;
	seq->nb_layers = 1;
	seq->internal_fits = calloc(size, sizeof(fits *));
	seq->seqname = strdup(_("internal sequence"));
	seq->imgparam = calloc(size, sizeof(imgdata));
	for (i = 0; i < size; i++) {
		seq->imgparam[i].filenum = i;
		seq->imgparam[i].incl = 1;
		seq->imgparam[i].date_obs = NULL;
	}
	check_or_allocate_regparam(seq, 0);
	return seq;
}

void internal_sequence_set(sequence *seq, int index, fits *fit) {
	assert(seq);
	assert(seq->internal_fits);
	assert(index < seq->number);
	seq->internal_fits[index] = fit;
}

fits *internal_sequence_get(sequence *seq, int index) {
	if (index > seq->number)
		return NULL;
	return seq->internal_fits[index];
}

// find index of the fit argument in the sequence
int internal_sequence_find_index(sequence *seq, const fits *fit) {
	int i;
	assert(seq);
	assert(seq->internal_fits);
	for (i = 0; i < seq->number; i++) {
		if (fit == seq->internal_fits[i])
			return i;
	}
	return -1;
}

// check if the passed sequence is used as a color sequence. It can be a CFA
// sequence explicitly demoisaiced too, which returns true.
gboolean sequence_is_rgb(sequence *seq) {
	switch (seq->type) {
		case SEQ_REGULAR:
			return seq->nb_layers == 3;
		case SEQ_SER:
			return (seq->ser_file->color_id != SER_MONO && com.pref.debayer.open_debayer) ||
				seq->ser_file->color_id == SER_RGB ||
				seq->ser_file->color_id == SER_BGR;
		case SEQ_FITSEQ:
			return seq->fitseq_file->naxes[2] == 3;
		default:
			return TRUE;
	}
}

/* Ensures that an area does not derive off-image.
 * Verifies coordinates of the center and moves it inside the image if the area crosses the bounds.
 */
gboolean enforce_area_in_image(rectangle *area, sequence *seq, int index) {
	gboolean has_crossed = FALSE;
	// need to check against current image size in case not the same as seq->rx/ry
	int rx = (seq->is_variable) ? seq->imgparam[index].rx : seq->rx;
	int ry = (seq->is_variable) ? seq->imgparam[index].ry : seq->ry;
	if (area->x < 0) {
		area->x = 0;
		has_crossed = TRUE;
	}
	if (area->y < 0) {
		area->y = 0;
		has_crossed = TRUE;
	}
	if (area->x + area->w > rx) {
		area->x = rx - area->w;
		has_crossed = TRUE;
	}
	if (area->y + area->h > ry) {
		area->y = ry - area->h;
		has_crossed = TRUE;
	}
	return has_crossed;
}

/********************************************************************
 *                                             __                   *
 *                   ___  ___  __ _ _ __  ___ / _|                  *
 *                  / __|/ _ \/ _` | '_ \/ __| |_                   *
 *                  \__ \  __/ (_| | |_) \__ \  _|                  *
 *                  |___/\___|\__, | .__/|___/_|                    *
 *                               |_|_|                              *
 ********************************************************************/

/* Computes FWHM for a sequence image.
 * area is the area from which fit was extracted from the full frame.
 * when the framing is set to follow star, args->area is centered on the found star
 */
int seqpsf_image_hook(struct generic_seq_args *args, int out_index, int index, fits *fit, rectangle *area, int threads) {
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;
	struct seqpsf_data *data = calloc(1, sizeof(struct seqpsf_data));
	if (!data) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	data->image_index = index;

	/* Backup the original pointer to fit. If there is a Bayer pattern we need
	 * to interpolate non-green pixels, so make a copy we can work on. */
	fits *orig_fit = fit;
	gboolean handle_cfa = fit_is_cfa(fit);
	if (handle_cfa) {
		fit = calloc(1, sizeof(fits));
		copyfits(orig_fit, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		interpolate_nongreen(fit);
	}
#if BAYER_DEBUG
	const gchar *green_filename = get_sequence_cache_filename(args->seq, index, "fit", "green_");
	savefits(green_filename, fit);
#endif

	rectangle psfarea = { .x = 0, .y = 0, .w = fit->rx, .h = fit->ry };
	psf_error error;
	struct phot_config *ps = NULL;
	if (spsfargs->for_photometry)
		ps = phot_set_adjusted_for_image(fit);
	data->psf = psf_get_minimisation(fit, 0, &psfarea, spsfargs->for_photometry, spsfargs->init_from_center, ps, TRUE, com.pref.starfinder_conf.profile, &error);
	free(ps);
	if (data->psf) {
		/* for photometry ? */
		if (spsfargs->for_photometry) {
			if (data->psf->s_mag > 9.0 || !data->psf->phot_is_valid) {
				siril_log_color_message(_("Photometry analysis failed for image %d (%s)\n"),
						"salmon", index, psf_error_to_string(error));
			}
		}
		// TODO: should we check for error or not?
		// for 3 stars reg, even though it did not converge, we still want to 
		// have the data to make the registration even with a poor fit (or not?)
		data->psf->xpos = data->psf->x0 + area->x;
		if (fit->top_down)
			data->psf->ypos = data->psf->y0 + area->y;
		else data->psf->ypos = area->y + area->h - data->psf->y0;

		/* let's move args->area to center it on the star */
		if (spsfargs->framing == FOLLOW_STAR_FRAME) {
			args->area.x = round_to_int(data->psf->xpos - args->area.w*0.5);
			args->area.y = round_to_int(data->psf->ypos - args->area.h*0.5);
			siril_debug_print("moving area to %d, %d\n", args->area.x, args->area.y);
		}

		if (!args->seq->imgparam[index].date_obs && fit->keywords.date_obs)
			args->seq->imgparam[index].date_obs = g_date_time_ref(fit->keywords.date_obs);
		data->exposure = fit->keywords.exposure;

		args->seq->imgparam[index].airmass = fit->keywords.airmass;
	}
	else {
		if (spsfargs->framing == FOLLOW_STAR_FRAME) {
			siril_log_color_message(_("No star found in the area image %d around %d,%d:"
						" error %s (use a larger area?)\n"),
					"red", index, area->x, area->y, psf_error_to_string(error));
		} else {
			siril_log_color_message(_("No star found in the area image %d around %d,%d:"
					" error %s (use 'follow star' option?)\n"),
				"red", index, area->x, area->y, psf_error_to_string(error));
		}
	}
	if (handle_cfa) {
		// Get rid of the temporary copy and restore the original frame fits
		// now that we have computed the actual registration data
		clearfits(fit);
		free(fit);
		fit = orig_fit;
	}
#ifdef _OPENMP
	omp_set_lock(&args->lock);
#endif
	spsfargs->list = g_slist_prepend(spsfargs->list, data);
#ifdef _OPENMP
	omp_unset_lock(&args->lock);
#endif
	return !data->psf;
}

static void write_regdata(sequence *seq, int layer, GSList *list, gboolean duplicate_for_regdata) {
	check_or_allocate_regparam(seq, layer);
	double xref = 0., yref = 0.;
	gboolean ref_set = FALSE;
	GSList *iterator;
	for (iterator = list; iterator; iterator = iterator->next) {
		struct seqpsf_data *data = iterator->data;
		seq->regparam[layer][data->image_index].fwhm_data =
			duplicate_for_regdata ? duplicate_psf(data->psf) : data->psf;
		if (data->psf) {
			seq->regparam[layer][data->image_index].fwhm = data->psf->fwhmx;
			seq->regparam[layer][data->image_index].roundness =
				data->psf->fwhmy / data->psf->fwhmx;
			seq->regparam[layer][data->image_index].weighted_fwhm = data->psf->fwhmx;
			seq->regparam[layer][data->image_index].background_lvl = data->psf->B;
			seq->regparam[layer][data->image_index].number_of_stars = 1;
			// we compute shift wrt to the first image which has data (even if not the ref)
			// the shifts will then be correctly computed as the Homography matrices are recomposed
			if (!ref_set) {
				xref = data->psf->xpos;
				yref = data->psf->ypos;
				ref_set = TRUE;
			}
			double shiftx = data->psf->xpos - xref;
			double shifty = data->psf->ypos - yref;
			seq->regparam[layer][data->image_index].H = H_from_translation(-shiftx, shifty);
		}
	}
	seq->needs_saving = TRUE;
}

int seqpsf_finalize_hook(struct generic_seq_args *args) {
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;
	sequence *seq = args->seq;
	int photometry_index = 0;
	gboolean displayed_warning = FALSE;

	if (args->retval)
		return 0;

	if (!spsfargs->for_photometry) {
		if (spsfargs->allow_use_as_regdata == BOOL_TRUE) {
			write_regdata(seq, args->layer_for_partial, spsfargs->list, FALSE);
		}
		return 0;
	}

	int i;
	for (i = 0; i < MAX_SEQPSF && seq->photometry[i]; i++);
	if (i == MAX_SEQPSF) {
		free_photometry_set(seq, 1);
		i = 1;
	}
	else seq->photometry[i+1] = NULL;
	seq->photometry[i] = calloc(seq->number, sizeof(psf_star *));
	photometry_index = i;

	for (GSList *iterator = spsfargs->list; iterator; iterator = iterator->next) {
		struct seqpsf_data *data = iterator->data;

		/* check exposure consistency */
		if (seq->exposure > 0.0 && data->psf && seq->exposure != data->exposure &&
				!displayed_warning) {
			siril_log_color_message(_("Star analysis does not give consistent results when exposure changes across the sequence.\n"), "red");
			displayed_warning = TRUE;
		}
		seq->exposure = data->exposure;

		// for photometry use: store data in seq->photometry
		seq->photometry[photometry_index][data->image_index] = data->psf;
	}

	if (args->already_in_a_thread || com.script || com.python_script) { // the idle won't be called
		// printing results ordered, the list isn't
		gboolean first = TRUE;
		for (int j = 0; j < seq->number; j++) {
			if (seq->photometry[photometry_index][j]) {
				psf_star *psf = seq->photometry[photometry_index][j];
				if (first) {
					siril_log_message(_("Photometry for star at %.1f, %.1f in image %d\n"), psf->xpos, psf->ypos, j);
					siril_debug_print("image_index magnitude error fwhm amplitude background\n");
					first = FALSE;
				}
				siril_debug_print("%d %f %f %f %f %f\n", j, psf->mag, psf->s_mag, psf->fwhmx, psf->A, psf->B);
			}
		}

		/* the idle below won't be called, we free data here */
		if (spsfargs->list)
			g_slist_free_full(spsfargs->list, free);
		memset(&com.selection, 0, sizeof(rectangle)); // we don't call delete_selected_area to avoid its idle when running python scripts
		free(spsfargs);
		args->user = NULL;
	}

	return 0;
}

// only does something if allow_use_as_regdata != false or GUI can be used
gboolean end_seqpsf(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *)p;
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;
	sequence *seq = args->seq;
	int layer = args->layer_for_partial;
	gboolean write_to_regdata, duplicate_for_regdata;

	if (args->retval)
		goto proper_ending;

	if (!spsfargs->for_photometry) {
		write_to_regdata = TRUE;
		duplicate_for_regdata = FALSE;
	} else {
		// for photometry data was saved in seq.photometry before, so we duplicate
		duplicate_for_regdata = TRUE;
		if (spsfargs->allow_use_as_regdata == BOOL_FALSE)
			write_to_regdata = FALSE;
		else if (spsfargs->allow_use_as_regdata == BOOL_TRUE)
			write_to_regdata = TRUE;
		else {
			gboolean has_any_regdata = FALSE;
			for (int i = 0; i < seq->nb_layers; i++)
				has_any_regdata = has_any_regdata || seq->regparam[i];
			if (!seq->regparam[layer] && has_any_regdata) {
				write_to_regdata = siril_confirm_dialog(_("No registration data stored for this layer"),
						_("Some registration data was found for another layer.\n"
							"Do you want to save the psf data as registration data for this layer?"), _("Save"));
			}
			else write_to_regdata = !seq->regparam[layer];
		}
	}

	if (write_to_regdata) {
		write_regdata(seq, layer, spsfargs->list, duplicate_for_regdata);
	}

	if (seq->needs_saving)
		writeseqfile(seq);

	// GUI things
	if (seq == &com.seq) {
		set_fwhm_star_as_star_list_with_layer(seq, layer);

		/* do here all GUI-related items, because it runs in the main thread.
		 * Most of these things are already done in end_register_idle
		 * in case seqpsf is called for registration. */
		if (seq->type != SEQ_INTERNAL) {
			// update the list in the GUI
			update_seqlist(layer);
			fill_sequence_list(seq, layer, FALSE);
		}
		set_layers_for_registration();	// update display of available reg data
		drawPlot();
		notify_new_photometry();	// switch to and update plot tab
		redraw(REDRAW_OVERLAY);
	}

proper_ending:
	if (spsfargs->list)
		g_slist_free_full(spsfargs->list, free);

	if (seq == &com.seq)
		adjust_sellabel();

	if (!check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);

	free(spsfargs);
	free(args);
	return end_generic(NULL);
}

/* process PSF for the given sequence, on the given layer, the area of the
 * image selection (com.selection), as a threaded operation or not.
 * expects seq->current to be valid
 */
int seqpsf(sequence *seq, int layer, gboolean for_registration,
		gboolean init_from_center, gboolean regall, framing_mode framing,
		gboolean run_in_thread, gboolean no_GUI) {

	if (framing == REGISTERED_FRAME && !layer_has_usable_registration(seq, layer))
		framing = ORIGINAL_FRAME;

	if (com.selection.w <= 0 || com.selection.h <= 0){
		siril_log_message(_("Select an area first\n"));
		return 1;
	}

	struct generic_seq_args *args = create_default_seqargs(seq);
	struct seqpsf_args *spsfargs = calloc(1, sizeof(struct seqpsf_args));

	spsfargs->for_photometry = !for_registration;
	if (!no_GUI)
		spsfargs->allow_use_as_regdata = BOOL_NOT_SET;
	else spsfargs->allow_use_as_regdata = for_registration ? BOOL_TRUE : BOOL_FALSE;
	spsfargs->framing = framing;
	spsfargs->list = NULL;	// GSList init is NULL
	spsfargs->init_from_center = init_from_center;

	args->partial_image = TRUE;
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	if (framing == REGISTERED_FRAME) {
		if (seq->reference_image < 0) seq->reference_image = sequence_find_refimage(seq);
		if (guess_transform_from_H(seq->regparam[layer][seq->reference_image].H) == NULL_TRANSFORMATION) {
			siril_log_color_message(_("The reference image has a null matrix and was not previously registered. Please select another one.\n"), "red");
			free(args);
			free(spsfargs);
			return 1;
		}
		// transform selection back from current to ref frame coordinates
		if (seq->current != seq->reference_image) {
			if (guess_transform_from_H(seq->regparam[layer][seq->current].H) == NULL_TRANSFORMATION) {
				siril_log_color_message(_("The current image has a null matrix and was not previously registered. Please load another one to select the star.\n"), "red");
				free(args);
				free(spsfargs);
				return 1;
			}
			selection_H_transform(&args->area, seq->regparam[layer][seq->current].H, seq->regparam[layer][seq->reference_image].H);
			if (args->area.x < 0 || args->area.x > seq->rx - args->area.w ||
					args->area.y < 0 || args->area.y > seq->ry - args->area.h) {
				siril_log_color_message(_("This area is outside of the reference image. Please select the reference image to select another star.\n"), "red");
				free(args);
				free(spsfargs);
				return 1;
			}
		}
	}

	if (framing == FOLLOW_STAR_FRAME)
		siril_log_color_message(_("The sequence analysis of the PSF will use a sliding selection area centred on the previous found star; this disables parallel processing.\n"), "salmon");
	else if (framing == REGISTERED_FRAME)
		siril_log_color_message(_("The sequence analysis of the PSF will use registration data to move the selection area for each image; this is compatible with parallel processing.\n"), "salmon");
	args->layer_for_partial = layer;
	args->regdata_for_partial = framing == REGISTERED_FRAME;
	args->get_photometry_data_for_partial = !for_registration;
	if (!regall) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = seq->selnum;
	}
	args->image_hook = seqpsf_image_hook;
	args->finalize_hook = seqpsf_finalize_hook;
	args->idle_function = end_seqpsf;
	args->stop_on_error = FALSE;
	args->description = _("PSF on area");
	args->user = spsfargs;
	args->already_in_a_thread = !run_in_thread;
	args->parallel = framing != FOLLOW_STAR_FRAME;

	if (run_in_thread) {
		if (!start_in_new_thread(generic_sequence_worker, args)) {
			free(spsfargs);
			free_generic_seq_args(args, TRUE);
			return 1;
		}
		return 0;
	} else {
		int retval = GPOINTER_TO_INT(generic_sequence_worker(args));
		free(args);
		return retval;
	}
}

void free_reference_image() {
	fprintf(stdout, "Purging previously saved reference frame data.\n");
	if (gui.refimage_regbuffer) {
		free(gui.refimage_regbuffer);
		gui.refimage_regbuffer = NULL;
	}
	if (gui.refimage_surface) {
		cairo_surface_destroy(gui.refimage_surface);
		gui.refimage_surface = NULL;
	}
	if (com.seq.reference_image == -1)
		enable_view_reference_checkbox(FALSE);
}

/* returns the number of images of the sequence that can fit into memory based
 * on the configured memory ratio */
static int compute_nb_images_fit_memory_from_dimensions(int rx, int ry, int nb_layers, data_type type, double factor, gboolean force_float, unsigned int *MB_per_orig_image, unsigned int *MB_per_scaled_image, unsigned int *max_mem_MB) {
	int max_memory_MB = get_max_memory_in_MB();

	if (factor > 3.0) {
		siril_debug_print("Info: image scaling is very large! (> 3.0)\n");
	}
	uint64_t newx = round_to_int((double) rx * factor);
	uint64_t newy = round_to_int((double) ry * factor);
	uint64_t memory_per_orig_image = (uint64_t) rx * ry * nb_layers;
	uint64_t memory_per_scaled_image = newx * newy * nb_layers;
	if (force_float || type == DATA_FLOAT) {
		memory_per_orig_image *= sizeof(float);
		memory_per_scaled_image *= sizeof(float);
	} else {
		memory_per_orig_image *= sizeof(WORD);
		memory_per_scaled_image *= sizeof(WORD);
	}
	unsigned int memory_per_orig_image_MB = memory_per_orig_image / BYTES_IN_A_MB;
	unsigned int memory_per_scaled_image_MB = memory_per_scaled_image / BYTES_IN_A_MB;
	if (memory_per_scaled_image_MB == 0) // in theory we should only do this if factor > 0. But just case we make a division by memory_per_scaled_image_MB... we'll keep this here for now
		memory_per_scaled_image_MB = 1;
	if (memory_per_orig_image_MB == 0)
		memory_per_orig_image_MB = 1;
	fprintf(stdout, "Memory per image: %u MB. Max memory: %d MB\n", memory_per_scaled_image_MB, max_memory_MB);
	if (MB_per_orig_image)
		*MB_per_orig_image = memory_per_orig_image_MB;
	if (MB_per_scaled_image)
		*MB_per_scaled_image = memory_per_scaled_image_MB;
	if (max_mem_MB)
		*max_mem_MB = max_memory_MB;
	return max_memory_MB / memory_per_scaled_image_MB;
}

size_t get_max_seq_dimension(sequence *seq, int *rx, int *ry) {
	*rx = seq->rx;
	*ry = seq->ry;
	size_t maxdim = 0;
	if (!seq->is_variable || !seq->imgparam)
		return seq->rx * seq->ry;
	for (int i = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl)
			continue;
		size_t dim = seq->imgparam[i].rx * seq->imgparam[i].ry;
		if (dim > maxdim) {
			maxdim = dim;
			*rx = seq->imgparam[i].rx;
			*ry = seq->imgparam[i].ry;
		}
	}
	return maxdim;
}

/* returns the number of images of the sequence that can fit into memory based
 * on the configured memory ratio
 * If the sequence is of variable size, we use the size of the largest image
*/
int compute_nb_images_fit_memory(sequence *seq, double factor, gboolean force_float, unsigned int *MB_per_orig_image, unsigned int *MB_per_scaled_image, unsigned int *max_mem_MB) {
	int rx = 0, ry = 0;
	get_max_seq_dimension(seq, &rx, &ry);
	return compute_nb_images_fit_memory_from_dimensions(rx, ry, seq->nb_layers, get_data_type(seq->bitpix), factor, force_float, MB_per_orig_image, MB_per_scaled_image, max_mem_MB);
}

/* returns the number of images of the sequence that can fit into memory based
 * on the configured memory ratio */
int compute_nb_images_fit_memory_from_fit(fits *fit, double factor, gboolean force_float, unsigned int *MB_per_orig_image, unsigned int *MB_per_scaled_image, unsigned int *max_mem_MB) {
	return compute_nb_images_fit_memory_from_dimensions(fit->rx, fit->ry, fit->naxes[2], fit->type, factor, force_float, MB_per_orig_image, MB_per_scaled_image, max_mem_MB);
}

void fix_selnum(sequence *seq, gboolean warn) {
	int nbsel = 0;
	for (int i = 0; i < seq->number; i++)
		if (seq->imgparam[i].incl)
			nbsel++;

	if (nbsel != seq->selnum) {
		if (warn)
			siril_log_message(_("Fixing the selection number in the .seq file (%d) to the actual value (%d) (not saved)\n"), seq->selnum, nbsel);
		seq->selnum = nbsel;
	}
}

gboolean sequence_ref_has_wcs(sequence *seq) {
	int refidx = sequence_find_refimage(seq);
	fits ref = { 0 };
	if (seq_read_frame_metadata(seq, refidx, &ref)) {
		siril_log_message(_("Could not load reference image\n"));
		return FALSE;
	}
	gboolean ret = has_wcs(&ref);
	clearfits(&ref);
	return ret;
}

struct wcsprm *get_wcs_ref(sequence *seq) {
	if (!seq)
		return NULL;
	struct wcsprm *wcsref = NULL;
	int refimage = sequence_find_refimage(seq);
	if (check_seq_is_comseq(seq) && seq->current == refimage && has_wcs(gfit)) { // we are in GUI
		wcsref = wcs_deepcopy(gfit->keywords.wcslib, NULL);
	} else { // we are in GUI with another image loaded or we are in script or headless, loading the seq has loaded the ref image, we check if it has wcs info
		fits ref = { 0 };
		if (seq_read_frame_metadata(seq, refimage, &ref)) {
			siril_log_message(_("Could not load reference image\n"));
			return FALSE;
		}
		if (has_wcs(&ref))
			wcsref = wcs_deepcopy(ref.keywords.wcslib, NULL);
		clearfits(&ref);
	}
	return wcsref;
}

gboolean sequence_drifts(sequence *seq, int reglayer, int threshold) {
	if (!seq->regparam || !seq->regparam[reglayer]) {
		siril_log_message(_("Sequence drift could not be checked as sequence has no regdata on layer %d\n"), reglayer);
		return FALSE;
	}
	double orig_x = (double)(seq->rx / 2.);
	double orig_y = (double)(seq->ry / 2.);
	for (int i = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl)
			continue;
		double x = orig_x, y = orig_y;
		cvTransfPoint(&x, &y, seq->regparam[reglayer][i].H, seq->regparam[reglayer][seq->reference_image].H, 1.);
		double dist = sqrt((x - orig_x) * (x - orig_x) + (y - orig_y) * (y - orig_y));
		if (dist > threshold) {
			siril_log_color_message(_("Warning: the sequence appears to have heavy drifted images (%d pixels for image %d), photometry will probably not be reliable. Check the sequence and exclude some images\n"), "salmon", (int)dist, i);
			return TRUE;
		}
	}
	siril_debug_print("no heavy drift detected\n");
	return FALSE;
}

void clean_sequence(sequence *seq, gboolean cleanreg, gboolean cleanstat, gboolean cleansel) {
	if (cleanreg && seq->regparam) {
		for (int i = 0; i < seq->nb_layers; i++) {
			if (seq->regparam[i]) {
				free(seq->regparam[i]);
				seq->regparam[i] = NULL;
				siril_log_message(_("Registration data cleared for layer %d\n"), i);
			}
		}
		// remove_prefixed_star_files(seq, "");
	}
	if (cleanreg && seq->regparam_bkp) {
		for (int i = 0; i < seq->nb_layers; i++) {
			if (seq->regparam_bkp[i]) {
				free(seq->regparam_bkp[i]);
				seq->regparam_bkp[i] = NULL;
				siril_log_message(_("Registration (back-up) data cleared for layer %d\n"), i);
			}
		}
	}
	if (cleanstat && seq->stats) {
		for (int i = 0; i < seq->nb_layers; i++) {
			clear_stats(seq, i);
			siril_log_message(_("Statistics data cleared for layer %d\n"), i);
		}
	}
	if (cleanstat && seq->stats_bkp) {
		for (int i = 0; i < seq->nb_layers; i++) {
			clear_stats_bkp(seq, i);
			siril_log_message(_("Statistics data (back-up) cleared for layer %d\n"), i);
		}
	}
	if (cleansel && seq->imgparam) {
		for (int i = 0; i < seq->number; i++) {
			seq->imgparam[i].incl = SEQUENCE_DEFAULT_INCLUDE;
		}
		// unsetting ref image
		seq->reference_image = -1;
		fix_selnum(seq, TRUE);
	}
	writeseqfile(seq);
}

// returns CACHE_NEWER if cache_filename last modification is more recent than image (for FITS) or
// sequence (for FITSEQ or SER). This is used for mask files (*.msk) and star files (*.lst)
// returns CACHE_NOT_FOUND if cache_filename does not exist
// returns CACHE_OLDER if cache_filename is older than image or sequence
// for FITS and lst files, we accept 30s difference because FITS may be updated by platesolving which also caches stars
// for FITSEQ and SER, we don't accept this as original sequence is not altered by platesolving
// we cannot use st_mtime because it is not reliable on all systems
cache_status check_cachefile_date(sequence *seq, int index, const gchar *cache_filename) {
	if (seq->type == SEQ_INTERNAL)
		return CACHE_NOT_FOUND; // internal sequences don't have filenames or cache files

	if (!g_file_test(cache_filename, G_FILE_TEST_EXISTS))
		return CACHE_NOT_FOUND;

	struct stat imgfileInfo, cachefileInfo;
	// if sequence is FITS, we check individual img file date vs cachefile date
	if (seq->type == SEQ_REGULAR) {
		char last_char = cache_filename[strlen(cache_filename) - 1];
		int margin = (last_char == 't') ? 30 : 0; // we use a margin for lst files but not for msk files
		char img_filename[256];
		if (!fit_sequence_get_image_filename_checkext(seq, index, img_filename) ||
				!g_file_test(img_filename, G_FILE_TEST_EXISTS) ||
				stat(img_filename, &imgfileInfo) ||
				stat(cache_filename, &cachefileInfo))
			return CACHE_NOT_FOUND;
		if (cachefileInfo.st_mtime < imgfileInfo.st_mtime - margin) {
			siril_debug_print("%s is older than %s, removing\n", cache_filename, img_filename);
			if (g_unlink(cache_filename))
				siril_debug_print(_("Removed outdated cache file %s failed\n"), cache_filename);
			return CACHE_OLDER;
		}
		return CACHE_NEWER;
	}
	// else, we check the sequence date vs cachefile date
	gchar *seqname;
	if (seq->type == SEQ_SER)
		seqname = seq->ser_file->filename;
	else seqname = seq->fitseq_file->filename;
	if (stat(seqname, &imgfileInfo) || stat(cache_filename, &cachefileInfo))
		return CACHE_NOT_FOUND;
	if (cachefileInfo.st_mtime < imgfileInfo.st_mtime) {
		siril_debug_print("%s is older than %s, removing\n", cache_filename, seqname);
		if (g_unlink(cache_filename)) {
			siril_debug_print(_("Removed outdated cache file %s failed\n"), cache_filename);
		}
		return CACHE_OLDER;
	}
	return CACHE_NEWER;
}

gchar *get_sequence_cache_filename(sequence *seq, int index, const gchar *cachefolder, const gchar *ext, const gchar *prefix) {
	char root[256];
	if (!fit_sequence_get_image_filename(seq, index, root, FALSE)) {
		return NULL;
	}
	gchar *base_root = g_path_get_basename(root);
	gchar *cache_filename = NULL;
	if (prefix)
		cache_filename = g_strdup_printf("%s%s.%s", prefix, base_root, ext);
	else
		cache_filename = g_strdup_printf("%s.%s", base_root, ext);
	if (g_strcmp0(ext, "fit") == 0) {
		gchar *tmp_filename = set_right_extension(cache_filename);
		g_free(cache_filename);
		cache_filename = tmp_filename;
	}
	gchar *cache_path = g_build_path(G_DIR_SEPARATOR_S, com.wd, cachefolder, cache_filename, NULL);
	g_free(cache_filename);
	g_free(base_root);
	return cache_path;
}
