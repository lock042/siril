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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

#include "core/siril.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "io/fits_sequence.h"
#include "io/ser.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "gui/progress_and_log.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#ifdef HAVE_FFMS2
#include "io/films.h"
#endif

/* seqfile version history *
 * no version up to 0.9.9
 * version 1 introduced roundness in regdata, 0.9.9
 * version 2 allowed regdata to be stored for CFA SER sequences, 0.9.11
 * version 3 introduced new weighted fwhm criteria, 0.99.0
 * version 4 introduced variable size sequences, extended registration data (incl. H), 1.1.0
 * version 5:
 * 	- removed upscale at stacking (U card) 1.3.4
 *  - added D* cards containing distortion and astrometry information 1.3.4  - see enum disto_source
 *  - added overlap statistics in the O* cards:
 *  	=> ON i j areai.x areai.y areaj.x areaj.y areai.w areai.h Nij medij medji madij madji locij locji scaij scji
 * 		=> with N the layer number and i,j the ith and jth images of the sequence
 */
#define CURRENT_SEQFILE_VERSION 5	// to increment on format change

/* File format (lines starting with # are comments, lines that are (for all
 * something) need to be in all in sequence of this only type of line):
 *
 * S sequence_name beg number selnum fixed reference_image [version]
 * L nb_layers
 * (for all images) I filenum incl [width,height] [stats+] <- stats added at some point, removed in 0.9.9
 * (for all layers (x)) Rx regparam+
 * TS | TA | TF (type for ser or film (avi) or fits)
 * U up-scale_ratio -> discarded in v5
 * (for all images (y) and layers (x)) Mx-y stats+
 */

/* name is sequence filename, with or without .seq extension
 * It should always be used with seq_check_basic_data() because on first loading
 * of a .seq that was created from scan of the filesystem, number of layers and
 * image size are unknown and some properties of the sequence are null or unset.
 * Returns NULL if the sequence could not be loaded.
 */
sequence * readseqfile(const char *name){
	char line[512], *scanformat;
	char filename[512], *seqfilename;
	int i, nb_tokens, allocated = 0, current_layer = -1, image;
	int to_backup = 0, version = -1;
	FILE *seqfile;
	sequence *seq;
	imstats *stats;
	regdata *regparam;

	if (!name) return NULL;
	fprintf(stdout, "Reading sequence file `%s'.\n", name);

	if(!g_str_has_suffix(name, ".seq")){
		seqfilename = malloc(strlen(name) + 6);	/* 6 stands for a max length of 4 + '.' + '\0' */
		sprintf(seqfilename, "%s.seq", name);
	} else {
		seqfilename = strdup(name);
	}

	if ((seqfile = g_fopen(seqfilename, "r")) == NULL) {
		fprintf(stderr, "Reading sequence failed, file cannot be opened: %s.\n", seqfilename);
		free(seqfilename);
		return NULL;
	}

	seq = calloc(1, sizeof(sequence));
	initialize_sequence(seq, TRUE);
	i = 0;
	while (fgets(line, 511, seqfile)) {
		switch (line[0]) {
			case '#':
				continue;
			case 'S':
				/* The double quote as sequence name is a sequence with no name.
				 * Such sequences don't exist anymore. */
				assert(line[2] != '"');
				if (line[2] == '\'')	/* new format, quoted string */
					scanformat = "'%511[^']' %d %d %d %d %d %d %d %d";
				else scanformat = "%511s %d %d %d %d %d %d %d %d";

				if(sscanf(line+2, scanformat,
							filename, &seq->beg, &seq->number,
							&seq->selnum, &seq->fixed,
							&seq->reference_image, &version, &seq->is_variable, &seq->fz) < 6 ||
						allocated != 0){
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				if (seq->number == 0) {
					fprintf(stderr, "readseqfile: sequence is empty?\n");
					goto error;
				}
				if (version > CURRENT_SEQFILE_VERSION)
					siril_log_message(_("This sequence file was created by a version of "
								"siril that is newer than this one, it may not "
								"be loaded as expected\n"),
							"salmon");
				/* for now, only the R* line is not supported in the previous version */
				seq->seqname = strdup(filename);
				seq->imgparam = calloc(seq->number, sizeof(imgdata));
				allocated = 1;
				break;

			case 'L':
				/* for now, the L line stores the number of layers for each image. */
				if (line[1] == ' ') {
					int nbl_backup = seq->nb_layers;
					if (sscanf(line+2, "%d", &seq->nb_layers) != 1) {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
						goto error;
					}
					/* seq->nb_layers can be -1 when the sequence has not been
					 * opened for the first time, but it may already have been
					 * set in SER opening below, so we keep the backup in this
					 * case */
					if (nbl_backup > 0 && ser_is_cfa(seq->ser_file)) {
						if (com.pref.debayer.open_debayer)
							seq->nb_layers = nbl_backup;
						else {
							seq->nb_layers = 1;
						}
					}
					// else if nbl_backup is 3 but opening debayer is not
					// enabled, we keep 1 in the nb_layers, which will be set in
					// the seq_check_basic_data() call later
					if (seq->nb_layers >= 1) {
						seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
						if (ser_is_cfa(seq->ser_file))
							seq->regparam_bkp = calloc(3, sizeof(regdata*));
						seq->distoparam = calloc(seq->nb_layers, sizeof(disto_params));
					}
				} else if (line[1] >= '0' && line[1] <= '9') {
					/* in the future, wavelength and name of each layer will be added here */
				}
				break;

			case 'I':
				/* First sequence file format was I filenum and incl.
				 * A later sequence file format added the stats to this.
				 * The current file format comes back to the first,
				 * moving stats to the M-line. */
				stats = NULL;
				if (!seq->imgparam) {
					fprintf(stderr, "readseqfile: sequence file format error, missing S line\n");
					goto error;
				}

				if (version <= 3) {
					allocate_stats(&stats);
					nb_tokens = sscanf(line + 2,
							"%d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
							&(seq->imgparam[i].filenum),
							&(seq->imgparam[i].incl),
							&(stats->mean),
							&(stats->median),
							&(stats->sigma),
							&(stats->avgDev),
							&(stats->mad),
							&(stats->sqrtbwmv),
							&(stats->location),
							&(stats->scale),
							&(stats->min),
							&(stats->max));
					if (nb_tokens == 12) {
						add_stats_to_seq(seq, i, 0, stats);
						free_stats(stats);	// we unreference it here
					} else {
						free_stats(stats);
						if (nb_tokens != 2) {
							fprintf(stderr, "readseqfile: sequence file format error: %s\n", line);
							goto error;
						}
					}
				} else {
					// v4: width and height for variable sequences
					nb_tokens = sscanf(line + 2, "%d %d %d,%d",
							&(seq->imgparam[i].filenum),
							&(seq->imgparam[i].incl),
							&(seq->imgparam[i].rx),
							&(seq->imgparam[i].ry));
					if ((nb_tokens != 4 && seq->is_variable) ||
							(nb_tokens != 2 && !seq->is_variable)) {
						fprintf(stderr, "readseqfile: sequence file format error: %s\n", line);
						goto error;
					}
				}
				++i;
				break;
			case 'D': // Distortion data - from version 5 onwards
				current_layer = line[1] - '0';
				if (current_layer < 0 || current_layer > seq->nb_layers) {
					fprintf(stderr, "readseqfile: sequence file bad distortion layer: %s\n", line);
					goto error;
				}
				int index;
				char buf0[256], buf1[256], buf2[256];
				nb_tokens = sscanf(line + 3, "%d %s %s %s\n",
							&index,
							buf0, buf1, buf2);
				if (nb_tokens < 1 || nb_tokens > 4) {
					fprintf(stderr, "readseqfile: sequence file bad distortion param: %s\n", line);
					goto error;
				}
				if (!seq->distoparam)
					seq->distoparam = calloc(seq->nb_layers, sizeof(disto_params));
				seq->distoparam[current_layer].index = index;
				if (index == DISTO_FILE || index == DISTO_MASTER) {
					if (nb_tokens == 1) {
						fprintf(stderr, "readseqfile: sequence file bad distortion param: %s\n", line);
						goto error;
					}
					seq->distoparam[current_layer].filename = g_strdup(buf0);
				}
				if (index == DISTO_FILE_COMET) {
					if (nb_tokens < 3) {
						fprintf(stderr, "readseqfile: sequence file bad distortion param: %s\n", line);
						goto error;
					}
					seq->distoparam[current_layer].velocity.x = (float)g_strtod(buf0, NULL);
					seq->distoparam[current_layer].velocity.y = (float)g_strtod(buf1, NULL);
					if (nb_tokens > 3) {
						seq->distoparam[current_layer].filename = g_strdup(buf2);
					}
				}
				++i;
				break;
			case 'R':
				/* registration info */
				if (line[1] == '*') {
					/* these are registration data for the CFA channel, the
					 * star is a way to differentiate stats belonging to
					 * CFA and those belonging to the demosaiced red
					 * channel, both would have layer number 0 otherwise */
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.pref.debayer.open_debayer) {
						siril_debug_print("- using CFA registration info\n");
						to_backup = 0;
					} else {
						siril_debug_print("- backing up CFA registration info\n");
						to_backup = 1;
					}
					current_layer = 0;
				}
				else {
					to_backup = 0;
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.pref.debayer.open_debayer) {
						to_backup = 1;
						siril_debug_print("- stats: backing up demosaiced registration info\n");
					}
					current_layer = line[1] - '0';
				}

				if (current_layer < 0 || current_layer > 9 ||
						(!seq->cfa_opened_monochrome && current_layer >= seq->nb_layers) ||
						(seq->cfa_opened_monochrome && current_layer >= 3)) {
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}

				if (to_backup) {
					if (!seq->regparam_bkp) {
						fprintf(stderr, "readseqfile: sequence type probably changed from CFA to MONO, invalid file\n");
						goto error;
					}
					regparam = seq->regparam_bkp[current_layer];
				} else {
					if (!seq->regparam) {
						fprintf(stderr, "readseqfile: file contains registration data but not the basic information\n");
						goto error;
					}
					regparam = seq->regparam[current_layer];
				}

				if (!regparam) {
					regparam = calloc(seq->number, sizeof(regdata));
					if (!regparam) {
						PRINT_ALLOC_ERR;
						goto error;
					}
					i = 0;	// one line per image, starting with 0
					// reassign, because we didn't use a pointer
					if (to_backup)
						seq->regparam_bkp[current_layer] = regparam;
					else seq->regparam[current_layer] = regparam;
				}
				if (!seq->distoparam)
					seq->distoparam = calloc(seq->nb_layers, sizeof(disto_params));
				if (i >= seq->number) {
					fprintf(stderr, "\nreadseqfile: out of array bounds in reg info!\n\n");
					goto error;
				}
				if (version < 1) {
					float rot_centre_x, rot_centre_y, angle, shiftx, shifty;
					nb_tokens = sscanf(line+3, "%f %f %g %g %g %g %lg",
							&shiftx,
							&shifty,
							&rot_centre_x, &rot_centre_y,
							&angle,
							&(regparam[i].fwhm),
							&(regparam[i].quality));
					if (nb_tokens != 7) {
						if (nb_tokens == 3) {
							// old format, with quality as third token
							regparam[i].quality = rot_centre_x;
							// the rest is already zero due to the calloc
						} else {
							fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
							goto error;
						}
					}
					regparam[i].H = H_from_translation(shiftx, shifty);
				} else if (version <= 2) { // include version 1
					// version 2 with roundness instead of weird things
					float shiftx, shifty;
					if (sscanf(line+3, "%f %f %g %g %lg",
								&shiftx,
								&shifty,
								&(regparam[i].fwhm),
								&(regparam[i].roundness),
								&(regparam[i].quality)) != 5) {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
						goto error;
					}
					regparam[i].H = H_from_translation(shiftx, shifty);
				} else if (version == 3) {
					// version 3 with weighted_fwhm
					float shiftx, shifty;
					if (sscanf(line+3, "%f %f %g %g %g %lg",
								&shiftx,
								&shifty,
								&(regparam[i].fwhm),
								&(regparam[i].weighted_fwhm),
								&(regparam[i].roundness),
								&(regparam[i].quality)) != 6) {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
						goto error;
					}
					regparam[i].H = H_from_translation(shiftx, shifty);
				}
				else {
					// version 4 without shifts and with homography matrix
					// version 5
					if (sscanf(line+3, "%g %g %g %lg %g %d H %lg %lg %lg %lg %lg %lg %lg %lg %lg",
								&(regparam[i].fwhm),
								&(regparam[i].weighted_fwhm),
								&(regparam[i].roundness),
								&(regparam[i].quality),
								&(regparam[i].background_lvl),
								&(regparam[i].number_of_stars),
								&(regparam[i].H.h00),
								&(regparam[i].H.h01),
								&(regparam[i].H.h02),
								&(regparam[i].H.h10),
								&(regparam[i].H.h11),
								&(regparam[i].H.h12),
								&(regparam[i].H.h20),
								&(regparam[i].H.h21),
								&(regparam[i].H.h22)) != 15) {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
						goto error;
					}
				}
				++i;
				break;

			case 'T':
				/* sequence type (several files or a single file) */
				if (line[1] == 'S') {
					seq->type = SEQ_SER;
#ifdef HAVE_FFMS2
					seq->ext = "ser";
#endif
					if (seq->ser_file) break;
					seq->ser_file = malloc(sizeof(struct ser_struct));
					ser_init_struct(seq->ser_file);
					seqfilename[strlen(seqfilename)-1] = 'r';
					if (ser_open_file(seqfilename, seq->ser_file)) {
						free(seq->ser_file);
						seq->ser_file = NULL;
						goto error;
					}

					if (ser_is_cfa(seq->ser_file)) {
						if (!com.pref.debayer.open_debayer) {
							// we set this flag instead of relying on the
							// com.debayer.open_debayer flag which varies
							// as the user changes the GUI
							seq->cfa_opened_monochrome = TRUE;
						}
						seq->nb_layers = 3;
						if (seq->regparam)
							seq->regparam = realloc(seq->regparam, seq->nb_layers * sizeof(regdata *));
						seq->needs_saving = TRUE;
					}
				}
				else if (line[1] == 'F') {
					seq->type = SEQ_FITSEQ;
#ifdef HAVE_FFMS2
					seq->ext = get_com_ext(seq->fz) + 1;
#endif
					if (seq->fitseq_file) break;
					seq->fitseq_file = malloc(sizeof(struct fits_sequence));
					fitseq_init_struct(seq->fitseq_file);
					GString *fileString = g_string_new(filename);
					g_string_append(fileString, get_com_ext(seq->fz));
					seq->fitseq_file->filename = g_string_free(fileString, FALSE);
					if (fitseq_open(seq->fitseq_file->filename, seq->fitseq_file, READONLY)) {
						g_free(seq->fitseq_file->filename);
						free(seq->fitseq_file);
						seq->fitseq_file = NULL;
						goto error;
					}
				}
#ifdef HAVE_FFMS2
				else if (line[1] == 'A') {
					seq->type = SEQ_AVI;
					if (seq->film_file) break;
					seq->film_file = malloc(sizeof(struct film_struct));
					int nb_film = get_nb_film_ext_supported();
					gchar *filmname = NULL;

					for (int ii = 0; ii < nb_film; ii++) {
						GString *filmString;
						filmString = g_string_new(seqfilename);
						filmString = g_string_truncate(filmString, strlen(seqfilename) - 3);
						filmString = g_string_append(filmString, supported_film[ii].extension);
						g_free(filmname);
						filmname = g_string_free(filmString, FALSE);

						/* test for extension in lowercase, else in uppercase */
						if (g_file_test(filmname, G_FILE_TEST_EXISTS)) {
							break;
						} else {
							int len_ext;
							gchar *upcase;

							g_free(filmname);
							filmname = NULL;

							filmString = g_string_new(seqfilename);
							g_string_truncate(filmString, strlen(seqfilename) - 3);
							len_ext = strlen(supported_film[ii].extension);
							upcase = g_ascii_strup(supported_film[ii].extension, len_ext);
							filmString = g_string_append(filmString, upcase);
							filmname = g_string_free(filmString, FALSE);
							g_free(upcase);

							if (g_file_test(filmname, G_FILE_TEST_EXISTS)) {
								break;
							}
						}
					}

					if (film_open_file(filmname, seq->film_file)) {
						free(seq->film_file);
						seq->film_file = NULL;
						g_free(filmname);
						goto error;
					}
					else {
						film_display_info(seq->film_file);
						seq->ext = get_filename_ext(seq->film_file->filename);
						g_free(filmname);
					}
				}
				else seq->ext = "fit";
#endif
				break;

			case 'U': // for versions up to 4
				siril_log_message(_("Seq file had info to upscale at stacking. This must now be passed as a stacking option\n"));
				break;
			case 'M':
				/* stats may not exist for all images and layers so we use
				 * indices for them, the line is Mx-y with x the layer number
				 * and y the image index */

				if (line[1] == '*') {
					/* these are stats for the CFA channel, the star is a
					 * way to differentiate stats belonging to CFA and
					 * those belonging to the demosaiced red channel, both
					 * would have layer number 0 otherwise */
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.pref.debayer.open_debayer) {
						siril_debug_print("- stats: using CFA stats\n");
						to_backup = 0;
					} else {
						siril_debug_print("- stats: backing up CFA stats\n");
						to_backup = 1;
					}
					current_layer = 0;
				}
				else {
					to_backup = 0;
					if (seq->type == SEQ_SER && ser_is_cfa(seq->ser_file) &&
							!com.pref.debayer.open_debayer) {
						to_backup = 1;
						siril_debug_print("- stats: backing up demosaiced stats\n");
					}
					current_layer = line[1] - '0';
				}

				if (current_layer < 0 || current_layer > 9 || line[2] != '-') {
					fprintf(stderr, "readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				/*if (current_layer >= seq->nb_layers) {
				// it may happen when opening a CFA file in monochrome
				break;
				}*/
				stats = NULL;
				allocate_stats(&stats);
				nb_tokens = sscanf(line + 3,
						"%d %ld %ld %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
						&image,
						&(stats->total),
						&(stats->ngoodpix),
						&(stats->mean),
						&(stats->median),
						&(stats->sigma),
						&(stats->avgDev),
						&(stats->mad),
						&(stats->sqrtbwmv),
						&(stats->location),
						&(stats->scale),
						&(stats->min),
						&(stats->max),
						&(stats->normValue),
						&(stats->bgnoise));
				if (nb_tokens == 15) {
					if (to_backup)
						add_stats_to_seq_backup(seq, image, current_layer, stats);
					else add_stats_to_seq(seq, image, current_layer, stats);
					free_stats(stats);	// we unreference it here
				} else {
					free_stats(stats);
					fprintf(stderr, "readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				break;
			case 'O':
				current_layer = line[1] - '0';
				if (!seq->ostats) {
					seq->ostats = alloc_ostats(seq->nb_layers, seq->number);
					if (!seq->ostats) {
						PRINT_ALLOC_ERR;
						goto error;
					}
				}
				overlap_stats_t ostat = { 0 };
				nb_tokens = sscanf(line + 3,
					"%d %d %d %d %d %d %d %d %zu %g %g %g %g %g %g %g %g",
					&(ostat.i),
					&(ostat.j),
					&(ostat.areai.x),
					&(ostat.areai.y),
					&(ostat.areaj.x),
					&(ostat.areaj.y),
					&(ostat.areai.w),
					&(ostat.areai.h),
					&(ostat.Nij),
					&(ostat.medij),
					&(ostat.medji),
					&(ostat.madij),
					&(ostat.madji),
					&(ostat.locij),
					&(ostat.locji),
					&(ostat.scaij),
					&(ostat.scaji));
				if (nb_tokens != 17) {
					fprintf(stderr, "readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				ostat.areaj.w = ostat.areai.w;
				ostat.areaj.h = ostat.areai.h;
				int ijth = get_ijth_pair_index(seq->number, ostat.i, ostat.j);
				seq->ostats[current_layer][ijth] = ostat;
				break;
		}
	}
	if (!allocated) {
		siril_log_message(_("The sequence file %s seems to be corrupted\n"), seqfilename);
		goto error;
	}
	seq->needs_saving = FALSE;	// loading stats sets it to true
	fclose(seqfile);
	seq->end = seq->imgparam[seq->number-1].filenum;
	seq->current = -1;
	fix_selnum(seq, TRUE);

	// copy some regparam_bkp to regparam if it applies
	if (ser_is_cfa(seq->ser_file) && com.pref.debayer.open_debayer &&
			seq->regparam_bkp && seq->regparam_bkp[0] &&
			seq->regparam && seq->nb_layers == 3 && !seq->regparam[1]) {
		siril_log_color_message(_("%s: Copying registration data from non-demosaiced layer to green layer\n"), "salmon", seqfilename);
		seq->regparam[1] = calloc(seq->number, sizeof(regdata));
		for (image = 0; image < seq->number; image++) {
			memcpy(&seq->regparam[1][image], &seq->regparam_bkp[0][image], sizeof(regdata));
		}
	}


	free(seqfilename);
	return seq;
error:
	fclose(seqfile);
	if (seq->seqname)
		free(seq->seqname);
	free(seq);
	siril_log_message(_("Could not load sequence %s\n"), name);

	free(seqfilename);
	return NULL;
}

/* Saves the sequence in the seqname.seq file. */
int writeseqfile(sequence *seq){
	char *filename;
	FILE *seqfile;
	int i, layer;

	if (!seq->seqname || seq->seqname[0] == '\0') return 1;
	if (!seq->imgparam) return 1;
	filename = malloc(strlen(seq->seqname)+5);
	sprintf(filename, "%s.seq", seq->seqname);
	seqfile = g_fopen(filename, "w+t");
	if (seqfile == NULL) {
		perror("writeseqfile, fopen");
		fprintf(stderr, "Writing sequence file: cannot open %s for writing\n", filename);
		free(filename);
		return 1;
	}
	fprintf(stdout, "Writing sequence file %s\n", filename);
	free(filename);

	fprintf(seqfile,"#Siril sequence file. Contains list of images, selection, registration data and statistics\n");
	fprintf(seqfile,"#S 'sequence_name' start_index nb_images nb_selected fixed_len reference_image version variable_size fz_flag\n");
	fprintf(seqfile,"S '%s' %d %d %d %d %d %d %d %d\n",
			seq->seqname, seq->beg, seq->number, seq->selnum, seq->fixed,
			seq->reference_image, CURRENT_SEQFILE_VERSION, seq->is_variable, seq->fz);
	if (seq->type != SEQ_REGULAR) {
		char type;
		switch (seq->type) {
			default: // cannot happen
			case SEQ_SER: type = 'S'; break;
#ifdef HAVE_FFMS2
			case SEQ_AVI: type = 'A'; break;
#endif
			case SEQ_FITSEQ: type = 'F'; break;
		}
		/* sequence type, not needed for regular, S for ser, A for avi */
		fprintf(seqfile, "T%c\n", type);
	}

	fprintf(seqfile, "L %d\n", seq->nb_layers);

	for(i = 0; i < seq->number; ++i){
		if (seq->is_variable) {
			fprintf(seqfile,"I %d %d %d,%d\n",
					seq->imgparam[i].filenum,
					seq->imgparam[i].incl,
					seq->imgparam[i].rx,
					seq->imgparam[i].ry);
		} else {
			fprintf(seqfile,"I %d %d\n",
					seq->imgparam[i].filenum,
					seq->imgparam[i].incl);
		}
	}

	for (layer = 0; layer < seq->nb_layers; layer++) {
		if (seq->regparam && seq->regparam[layer]) {
			if (layer_has_distortion(seq, layer)) {
				if (seq->distoparam[layer].index == DISTO_FILE)
					fprintf(seqfile, "D%c %d %s\n",
					seq->cfa_opened_monochrome ? '*' : '0' + layer,
					DISTO_FILE,
					seq->distoparam[layer].filename);
				else if (seq->distoparam[layer].index == DISTO_FILES)
					fprintf(seqfile, "D%c %d\n",
					seq->cfa_opened_monochrome ? '*' : '0' + layer,
					DISTO_FILES);
				else if (seq->distoparam[layer].index == DISTO_MASTER) {
					fprintf(seqfile, "D%c %d %s\n",
					seq->cfa_opened_monochrome ? '*' : '0' + layer,
					DISTO_MASTER,
					seq->distoparam[layer].filename);
				}
				else if (seq->distoparam[layer].index == DISTO_FILE_COMET) {
					fprintf(seqfile, "D%c %d %.3f %.3f %s\n",
					seq->cfa_opened_monochrome ? '*' : '0' + layer,
					DISTO_FILE_COMET,
					seq->distoparam[layer].velocity.x,
					seq->distoparam[layer].velocity.y,
					seq->distoparam[layer].filename ? seq->distoparam[layer].filename : "");
				}
			}
			for (i = 0; i < seq->number; ++i) {
				fprintf(seqfile, "R%c %g %g %g %g %g %d H %g %g %g %g %g %g %g %g %g\n",
						seq->cfa_opened_monochrome ? '*' : '0' + layer,
						seq->regparam[layer][i].fwhm,
						seq->regparam[layer][i].weighted_fwhm,
						seq->regparam[layer][i].roundness,
						seq->regparam[layer][i].quality,
						seq->regparam[layer][i].background_lvl,
						seq->regparam[layer][i].number_of_stars,
						seq->regparam[layer][i].H.h00,
						seq->regparam[layer][i].H.h01,
						seq->regparam[layer][i].H.h02,
						seq->regparam[layer][i].H.h10,
						seq->regparam[layer][i].H.h11,
						seq->regparam[layer][i].H.h12,
						seq->regparam[layer][i].H.h20,
						seq->regparam[layer][i].H.h21,
						seq->regparam[layer][i].H.h22
					);
			}
		}
		if (seq->stats && seq->stats[layer]) {
			for (i = 0; i < seq->number; ++i) {
				if (!seq->stats[layer][i]) continue;

				fprintf(seqfile, "M%c-%d %ld %ld %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
						seq->cfa_opened_monochrome ? '*' : '0' + layer, i,
						seq->stats[layer][i]->total,
						seq->stats[layer][i]->ngoodpix,
						seq->stats[layer][i]->mean,
						seq->stats[layer][i]->median,
						seq->stats[layer][i]->sigma,
						seq->stats[layer][i]->avgDev,
						seq->stats[layer][i]->mad,
						seq->stats[layer][i]->sqrtbwmv,
						seq->stats[layer][i]->location,
						seq->stats[layer][i]->scale,
						seq->stats[layer][i]->min,
						seq->stats[layer][i]->max,
						seq->stats[layer][i]->normValue,
						seq->stats[layer][i]->bgnoise);

			}
		}
	}
	for (layer = 0; layer < 3; layer++) {
		if (seq->regparam_bkp && seq->regparam_bkp[layer]) {
			for (i = 0; i < seq->number; ++i) {
				fprintf(seqfile, "R%c %g %g %g %g %g %d H %g %g %g %g %g %g %g %g %g\n",
						seq->cfa_opened_monochrome ? '0' + layer : '*',
						seq->regparam_bkp[layer][i].fwhm,
						seq->regparam_bkp[layer][i].weighted_fwhm,
						seq->regparam_bkp[layer][i].roundness,
						seq->regparam_bkp[layer][i].quality,
						seq->regparam_bkp[layer][i].background_lvl,
						seq->regparam_bkp[layer][i].number_of_stars,
						seq->regparam_bkp[layer][i].H.h00,
						seq->regparam_bkp[layer][i].H.h01,
						seq->regparam_bkp[layer][i].H.h02,
						seq->regparam_bkp[layer][i].H.h10,
						seq->regparam_bkp[layer][i].H.h11,
						seq->regparam_bkp[layer][i].H.h12,
						seq->regparam_bkp[layer][i].H.h20,
						seq->regparam_bkp[layer][i].H.h21,
						seq->regparam_bkp[layer][i].H.h22
					);
			}
		}
		if (seq->stats_bkp && seq->stats_bkp[layer]) {
			for (i = 0; i < seq->number; ++i) {
				if (!seq->stats_bkp[layer][i]) continue;

				fprintf(seqfile, "M%c-%d %ld %ld %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
						seq->cfa_opened_monochrome ? '0' + layer : '*', i,
						seq->stats_bkp[layer][i]->total,
						seq->stats_bkp[layer][i]->ngoodpix,
						seq->stats_bkp[layer][i]->mean,
						seq->stats_bkp[layer][i]->median,
						seq->stats_bkp[layer][i]->sigma,
						seq->stats_bkp[layer][i]->avgDev,
						seq->stats_bkp[layer][i]->mad,
						seq->stats_bkp[layer][i]->sqrtbwmv,
						seq->stats_bkp[layer][i]->location,
						seq->stats_bkp[layer][i]->scale,
						seq->stats_bkp[layer][i]->min,
						seq->stats_bkp[layer][i]->max,
						seq->stats_bkp[layer][i]->normValue,
						seq->stats_bkp[layer][i]->bgnoise);
			}
		}
	}
	if (seq->ostats) {
		for (layer = 0; layer < seq->nb_layers; layer++) {
			int Npairs = seq->number * (seq->number - 1) / 2;
			for (i = 0; i < Npairs; i++) {
				if (seq->ostats[layer][i].i == -1)
					continue;
				fprintf(seqfile, "O%d %d %d %d %d %d %d %d %d %zu %g %g %g %g %g %g %g %g\n",
					layer,
					seq->ostats[layer][i].i,
					seq->ostats[layer][i].j,
					seq->ostats[layer][i].areai.x,
					seq->ostats[layer][i].areai.y,
					seq->ostats[layer][i].areaj.x,
					seq->ostats[layer][i].areaj.y,
					seq->ostats[layer][i].areai.w,
					seq->ostats[layer][i].areai.h,
					seq->ostats[layer][i].Nij,
					seq->ostats[layer][i].medij,
					seq->ostats[layer][i].medji,
					seq->ostats[layer][i].madij,
					seq->ostats[layer][i].madji,
					seq->ostats[layer][i].locij,
					seq->ostats[layer][i].locji,
					seq->ostats[layer][i].scaij,
					seq->ostats[layer][i].scaji);
			}
		}
	}

	fclose(seqfile);
	seq->needs_saving = FALSE;
	return 0;
}

gboolean existseq(const char *name){
	char *filename;
	GStatBuf sts;
	if (!name || name[0] == '\0') return FALSE;
	filename = malloc(strlen(name)+5);
	sprintf(filename, "%s.seq", name);
	if(g_stat(filename, &sts)==0){
		free(filename);
		return TRUE;
	}
	free(filename);
	return FALSE;
}

/* try to create the sequence file for the newly found sequence */
int buildseqfile(sequence *seq, int force_recompute) {
	image_type imagetype;
	int i;
	char *filename;
	imgdata *oldparam;

	if (seq->end <= 0 || !seq->seqname || seq->seqname[0] == '\0') return 1;
	if (existseq(seq->seqname) && !force_recompute) {
		fprintf(stdout, "seqfile '%s.seq' already exists, not overwriting\n", seq->seqname);
		return 0;
	}

	if (force_recompute) {
		for (i = 0; i < seq->nb_layers; i++)
			clear_stats(seq, i);
	}

	filename = malloc(strlen(seq->seqname) + 20);
	if (filename == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	if (seq->type == SEQ_REGULAR) {
		get_possible_image_filename(seq, seq->beg, filename);
		// check if the sequence begins at first_index

		if (stat_file(filename, &imagetype, NULL) || imagetype != TYPEFITS) {
			siril_log_message(_("The sequence %s doesn't start at the frame number %d"
					" with the specified fixed size index (%d). Cannot load.\n"),
					seq->seqname, seq->beg, seq->fixed);
			free(filename);
			return 1;
		}
	}

	int alloc_size = 30;
	//seq->number = 0;
	// fill in one pass: realloc needed
	if (seq->type == SEQ_REGULAR) {
		if (seq->end - seq->beg < 111)
			alloc_size = seq->end - seq->beg + 1;	// last index IS included
	} else alloc_size = seq->end - seq->beg + 1;		// always continuous
	oldparam = seq->imgparam;
	if ((seq->imgparam = realloc(seq->imgparam, alloc_size*sizeof(imgdata))) == NULL) {
		fprintf(stderr, "Could not reallocate image parameters structure in sequence\n");
		if (oldparam) free(oldparam);
		free(filename);
		return 2;
	}
	for (i = seq->beg; i <= seq->end; i++) {
		if (seq->type == SEQ_REGULAR) {
			get_possible_image_filename(seq, i, filename);
			if (!stat_file(filename, &imagetype, NULL) && imagetype == TYPEFITS) {
				if (seq->number + 1 > alloc_size - 1) {
					alloc_size += 25;
					oldparam = seq->imgparam;
					if (!(seq->imgparam = realloc(seq->imgparam, alloc_size * sizeof(imgdata)))) {
						PRINT_ALLOC_ERR;
						if (oldparam) free(oldparam);
						free(filename);
						return 2;
					}
				}
				seq->imgparam[seq->number].filenum = i;
				seq->imgparam[seq->number].incl = SEQUENCE_DEFAULT_INCLUDE;
				seq->imgparam[seq->number].date_obs = NULL;
				seq->number++;
			}
		} else {
			seq->imgparam[i].filenum = i;
			seq->imgparam[i].incl = SEQUENCE_DEFAULT_INCLUDE;
			seq->imgparam[i].date_obs = NULL;
		}
	}
#if SEQUENCE_DEFAULT_INCLUDE == TRUE
	seq->selnum = seq->number;
#else
	seq->selnum = 0;
#endif
	writeseqfile(seq);

	fprintf(stdout, "Sequence found: %s %d->%d\n", seq->seqname, seq->beg, seq->end);
	free(filename);
	return 0;
}

