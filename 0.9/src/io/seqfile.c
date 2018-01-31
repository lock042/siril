/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include <unistd.h>
#include <assert.h>

#include "core/siril.h"
#include "io/ser.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
#include "io/films.h"
#endif

/* name is sequence filename, with or without .seq extension
 * It should always be used with seq_check_basic_data() because on first loading
 * of a .seq that was created from scan of the filesystem, number of layers and
 * image size are unknown and some properties of the sequence are null or unset.
 * Returns NULL if the sequence could not be loaded.
 */
sequence * readseqfile(const char *name){
	char line[512], *scanformat;
	char filename[512], *seqfilename;
	int i, nbsel;
	FILE *seqfile;
	int allocated=0;
	int current_layer = -1;
	sequence *seq;

	if (!name) return NULL;
	fprintf(stdout, "Reading sequence file `%s'.\n", name);

	if(!ends_with(name, ".seq")){
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
	i=0;
	while (fgets(line, 511, seqfile)) {
		switch (line[0]) {
			case '#':
				continue;
			case 'S':
				/* The double quote as sequence name is a sequence with no name.
				 * Such sequences don't exist anymore. */
				assert(line[2] != '"');
				if (line[2] == '\'')	/* new format, quoted string */
					scanformat = "'%511[^']' %d %d %d %d %d";
				else scanformat = "%511s %d %d %d %d %d";

				if(sscanf(line+2, scanformat,
							filename, &seq->beg, &seq->number,
							&seq->selnum, &seq->fixed,
							&seq->reference_image) != 6 ||
						allocated != 0){
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				if (seq->number == 0) {
					fprintf(stderr, "readseqfile: sequence is empty?\n");
					goto error;
				}
				seq->seqname = strdup(filename);
				seq->imgparam = calloc(seq->number, sizeof(imgdata));
				allocated = 1;
				break;

			case 'L':
				/* for now, the L line stores the number of layers for each image. */
				if (line[1] == ' ') {
					if (sscanf(line+2, "%d", &seq->nb_layers) != 1) {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
						goto error;
					}
					/* seq->nb_layers can be -1 when the sequence has not
					 * been opened for the first time */
					if (seq->nb_layers >= 1) {
						seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
						seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
					}
				} else if (line[1] >= '0' && line[1] <= '9') {
					/* in the future, wavelength and name of each layer will be added here */
				}
				break;
			case 'I':
				seq->imgparam[i].stats = malloc(sizeof(imstats));
				/* new format: with stats, if already computed, else it's old format */
				int nb_tokens = sscanf(line + 2,
						"%d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
						&(seq->imgparam[i].filenum),
						&(seq->imgparam[i].incl),
						&(seq->imgparam[i].stats->mean),
						&(seq->imgparam[i].stats->median),
						&(seq->imgparam[i].stats->sigma),
						&(seq->imgparam[i].stats->avgDev),
						&(seq->imgparam[i].stats->mad),
						&(seq->imgparam[i].stats->sqrtbwmv),
						&(seq->imgparam[i].stats->location),
						&(seq->imgparam[i].stats->scale),
						&(seq->imgparam[i].stats->min),
						&(seq->imgparam[i].stats->max));
				if (nb_tokens == 12) {
					if (seq->nb_layers == 1)
						strcpy(seq->imgparam[i].stats->layername, "B&W");
					else	strcpy(seq->imgparam[i].stats->layername, "Red");
				} else {
					if (nb_tokens == 2) {
						/* previously was the simple image description with no line key */
						free(seq->imgparam[i].stats);
						seq->imgparam[i].stats = NULL;
					}
					else {
						fprintf(stderr,"readseqfile: sequence file format error: %s\n", line);
						goto error;
					}
				}
				++i;
				break;
			case 'R':
				/* registration info */
				current_layer = line[1] - '0';
				if (current_layer < 0 || current_layer > 9) {
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				if (seq->regparam[current_layer] == NULL) {
					seq->regparam[current_layer] = calloc(seq->number, sizeof(regdata));
					if (seq->regparam[current_layer] == NULL) {
						fprintf(stderr, "readseqfile: could not allocate registration data\n");
						goto error;
					}
					i = 0;
				}
				if (i >= seq->number) {
					fprintf(stderr, "\nreadseqfile ERROR: out of array bounds in reg info!\n\n");
				} else {
					int nb_tokens = sscanf(line+3, "%f %f %g %g %g %g %lg",
							&(seq->regparam[current_layer][i].shiftx),
							&(seq->regparam[current_layer][i].shifty),
							&(seq->regparam[current_layer][i].rot_centre_x),
							&(seq->regparam[current_layer][i].rot_centre_y),
							&(seq->regparam[current_layer][i].angle),
							&(seq->regparam[current_layer][i].fwhm),
							&(seq->regparam[current_layer][i].quality));
					if (nb_tokens != 7) {
						if (nb_tokens == 3) {
							// old format, with quality as third token
							seq->regparam[current_layer][i].rot_centre_x = 0.0f;
							// the rest is already zero due to the calloc
						} else {
							fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
							goto error;
						}
					}
					++i;
				}
				break;
			case 'T':
				/* sequence type (several files or a single file) */
				if (line[1] == 'S') {
					seq->type = SEQ_SER;
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
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
					else ser_display_info(seq->ser_file);
				}
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
				else if (line[1] == 'A') {
					seq->type = SEQ_AVI;
					if (seq->film_file) break;
					seq->film_file = malloc(sizeof(struct film_struct));
					int i = 0, nb_film = get_nb_film_ext_supported();

					gchar *filmname;
					while (i < nb_film) {
						GString *filmString;
						/* test for extension in lowercase */
						filmString = g_string_new(seqfilename);
						g_string_truncate(filmString, strlen(seqfilename) - 3);
						g_string_append(filmString, supported_film[i].extension);
						filmname = g_string_free(filmString, FALSE);

						if (g_file_test(filmname, G_FILE_TEST_EXISTS)) {
							break;
						} else {
							int len_ext;
							gchar *upcase;

							g_free(filmname);

							filmString = g_string_new(seqfilename);
							g_string_truncate(filmString, strlen(seqfilename) - 3);
							len_ext = strlen(supported_film[i].extension);
							upcase = g_ascii_strup(supported_film[i].extension, len_ext);
							g_string_append(filmString, upcase);
							filmname = g_string_free(filmString, FALSE);
							g_free(upcase);

							if (g_file_test(filmname, G_FILE_TEST_EXISTS)) {
								break;
							}
						}
						i++;
					}

					if (film_open_file(filmname, seq->film_file)) {
						free(seq->film_file);
						seq->film_file = NULL;
						goto error;
					}
					else {
						film_display_info(seq->film_file);
						seq->ext = strdup(get_filename_ext(seq->film_file->filename));
					}
				}
				else seq->ext = "fit";
#endif
				break;
			case 'U':
				/* up-scale factor for stacking. Used in simplified stacking for
				 * shift-only registrated sequences, up-scale will be done at
				 * stack-time. */
				if (line[1] == ' ' &&
						sscanf(line+2, "%lg", &seq->upscale_at_stacking) != 1) {
					fprintf(stderr,"readseqfile: sequence file format error: %s\n",line);
					goto error;
				}
				break;
		}
	}
	if (!allocated) {
		siril_log_message(_("The file seems to be corrupted\n"));
		goto error;
	}
	fclose(seqfile);
	seq->end = seq->imgparam[seq->number-1].filenum;
	seq->current = -1;
	for (i=0, nbsel=0; i<seq->number; i++)
		if (seq->imgparam[i].incl)
			nbsel++;
	if (nbsel != seq->selnum) {
		siril_log_message(_("Fixing the selection number in the .seq file (%d) to the actual value (%d) (not saved)\n"), seq->selnum, nbsel);
		seq->selnum = nbsel;
	}
	update_used_memory();
	free(seqfilename);
	return seq;
error:
	fclose(seqfile);
	if (seq->seqname)
		free(seq->seqname);
	free(seq);
	siril_log_message(_("Could not load sequence %s\n"), name);
	update_used_memory();
	free(seqfilename);
	return NULL;
}

/* Saves the sequence in the seqname.seq file. */
int writeseqfile(sequence *seq){
	char *filename;
	FILE *seqfile;
	int i,j;

	if (!seq->seqname || seq->seqname[0] == '\0') return 1;
	filename = malloc(strlen(seq->seqname)+5);
	sprintf(filename, "%s.seq", seq->seqname);
	seqfile = fopen(filename, "w+");	// g_fopen won't work (on WINDOWS).
	if (seqfile == NULL) {
		perror("writeseqfile, fopen");
		fprintf(stderr, "Writing sequence file: cannot open %s for writing\n", filename);
		free(filename);
		return 1;
	}
	fprintf(stdout, "Writing sequence file %s\n", filename);
	free(filename);

	fprintf(seqfile,"#Siril sequence file. Contains list of files (images), selection, and registration data\n");
	fprintf(seqfile,"#S 'sequence_name' start_index nb_images nb_selected fixed_len reference_image\n");
	fprintf(stderr,"S '%s' %d %d %d %d %d\n", 
			seq->seqname, seq->beg, seq->number, seq->selnum, seq->fixed, seq->reference_image);
	fprintf(seqfile,"S '%s' %d %d %d %d %d\n", 
			seq->seqname, seq->beg, seq->number, seq->selnum, seq->fixed, seq->reference_image);
	if (seq->type != SEQ_REGULAR) {
		/* sequence type, not needed for regular, S for ser, A for avi */
		fprintf(stderr, "T%c\n", seq->type == SEQ_SER ? 'S' : 'A');
		fprintf(seqfile, "T%c\n", seq->type == SEQ_SER ? 'S' : 'A');
	}

	if (seq->upscale_at_stacking != 1.0) {
		// until we have a real drizzle
		fprintf(stderr, "U %g\n", seq->upscale_at_stacking);
		fprintf(seqfile, "U %g\n", seq->upscale_at_stacking);
	}

	fprintf(stderr, "L %d\n", seq->nb_layers);
	fprintf(seqfile, "L %d\n", seq->nb_layers);

	for(i=0; i < seq->number; ++i){
		if (seq->imgparam[i].stats) {
			fprintf(seqfile,"I %d %d %g %g %g %g %g %g %g %g %g %g\n",
					seq->imgparam[i].filenum, 
					seq->imgparam[i].incl,
					seq->imgparam[i].stats->mean,
					seq->imgparam[i].stats->median,
					seq->imgparam[i].stats->sigma,
					seq->imgparam[i].stats->avgDev,
					seq->imgparam[i].stats->mad,
					seq->imgparam[i].stats->sqrtbwmv,
					seq->imgparam[i].stats->location,
					seq->imgparam[i].stats->scale,
					seq->imgparam[i].stats->min,
					seq->imgparam[i].stats->max);
		} else {
			fprintf(seqfile,"I %d %d\n", seq->imgparam[i].filenum, 
					seq->imgparam[i].incl);
		}
	}

	for(j=0; j < seq->nb_layers; j++) {
		if (seq->regparam[j]) {
			for (i=0; i < seq->number; ++i) {
				/*fprintf(stderr, "R%d %f %f %g %g %g %g %g\n", j,
						seq->regparam[j][i].shiftx,
						seq->regparam[j][i].shifty,
						seq->regparam[j][i].rot_centre_x,
						seq->regparam[j][i].rot_centre_y,
						seq->regparam[j][i].angle,
						seq->regparam[j][i].fwhm,
						seq->regparam[j][i].quality
						);*/
				fprintf(seqfile, "R%d %f %f %g %g %g %g %g\n", j,
						seq->regparam[j][i].shiftx,
						seq->regparam[j][i].shifty,
						seq->regparam[j][i].rot_centre_x,
						seq->regparam[j][i].rot_centre_y,
						seq->regparam[j][i].angle,
						seq->regparam[j][i].fwhm,
						seq->regparam[j][i].quality
				       );
			}
		}
	}
	fclose(seqfile);
	seq->needs_saving = FALSE;
	return 0;
}

gboolean existseq(const char *name){
	char *filename;
	struct stat sts;
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
		fprintf(stderr,"seqfile '%s.seq' already exists, not recomputing\n", seq->seqname);
		return 0;
	}

	filename = malloc(strlen(seq->seqname) + 20);
	if (filename == NULL) {
		printf("alloc error: buildseqfile\n");
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
				if (seq->number+1 > alloc_size-1) {
					alloc_size += 25;
					oldparam = seq->imgparam;
					if (!(seq->imgparam = realloc(seq->imgparam, alloc_size*sizeof(imgdata)))) {
						fprintf(stderr, "Could not reallocate image parameters structure in sequence\n");
						if (oldparam) free(oldparam);
						free(filename);
						return 2;
					}
				}
				seq->imgparam[seq->number].filenum = i;
				seq->imgparam[seq->number].incl = SEQUENCE_DEFAULT_INCLUDE;
				seq->imgparam[seq->number].stats = NULL;
				seq->imgparam[seq->number].date_obs = NULL;
				seq->number++;
			}
		} else {
			seq->imgparam[i].filenum = i;
			seq->imgparam[i].incl = SEQUENCE_DEFAULT_INCLUDE;
			seq->imgparam[i].stats = NULL;
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

