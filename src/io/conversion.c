/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#ifdef _WIN32
#include <windows.h>
#include "core/OS_utils.h" // for siril_real_path()
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/demosaicing.h"
#include "io/films.h"
#include "io/fits_sequence.h"
#include "io/image_format_fits.h"
#include "io/ser.h"
#include "io/seqwriter.h"
#include "io/sequence.h"
#include "io/FITS_symlink.h"
#include "gui/progress_and_log.h"
#include "conversion.h"

#ifdef HAVE_LIBRAW
#include <libraw/libraw_version.h>
#endif

static unsigned int supported_filetypes = 0;	// initialized by initialize_converters()

// NULL-terminated array, initialized by initialize_converters(), used only by stat_file
char *supported_extensions[MAX_EXTENSIONS];

supported_raw_list supported_raw[] = {
	{"dng",	"Adobe"},
	{"mos",	"Aptus"},
	{"cr2",	"Canon"},
	{"crw",	"Canon"},
#ifdef HAVE_LIBRAW
#if LIBRAW_VERSION > LIBRAW_MAKE_VERSION(0, 19, 5)
	{"cr3",	"Canon"},
#endif
#endif
	{"bay",	"Casio"},		// Not tested
	{"erf",	"Epson"},
	{"raf",	"Fuji"},
	{"3fr",	"Hasselblad"},
	{"kdc",	"Kodak"},
	{"dcr",	"Kodak"},
	{"mef",	"Mamiya"},
	{"mrw",	"Minolta"},
	{"nef",	"Nikon"},
	{"nrw",	"Nikon"},
	{"orf",	"Olympus"},
	{"raw",	"Leica"},
	{"rw2",	"Panasonic"},
	{"pef",	"Pentax"},
	{"ptx",	"Pentax"},		// Not tested
	{"x3f",	"Sigma"},		// Not supported yet
	{"srw",	"Samsung"},
	{"arw",	"Sony"}
};

int get_nb_raw_supported() {
	return G_N_ELEMENTS(supported_raw);
}

/* This function is used with command line only */
void list_format_available() {
	puts("=======================================================");
	puts("[            Supported image file formats             ]");
	puts("=======================================================");
	puts("FITS\t(*.fit, *.fits, *.fts)");
	puts("BMP\t(*.bmp)");
	puts("NetPBM\t(*.ppm, *.pgm, *.pnm)");
	puts("PIC\t(*.pic)");
#ifdef HAVE_LIBRAW
	printf("RAW\t(");
	int i, nb_raw;

	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		printf("*.%s", supported_raw[i].extension);
		if (i != nb_raw - 1)
			printf(", ");
	}
	printf(")\n");
#endif

#ifdef HAVE_LIBTIFF
	puts("TIFF\t(*.tif, *.tiff)");
#endif
#ifdef HAVE_LIBXISF
	puts("XISF\t(*.xisf)");
#endif
#ifdef HAVE_LIBJPEG
	puts("JPEG\t(*.jpg, *.jpeg)");
#endif
#ifdef HAVE_LIBJXL
	puts("JPEG XL\t(*.jxl)");
#endif
#ifdef HAVE_LIBPNG
	puts("PNG\t(*.png)");
#endif
#ifdef HAVE_LIBHEIF
	puts("HEIF\t(*.heic, *.heif)");
	puts("AVIF\t(*.avif)");
#endif
}


/**************************Public functions***********************************************************/

/* initialize converters (utilities used for different image types importing) *
 * updates the label listing the supported input file formats, and modifies the
 * list of file types used in convflags */
gchar *initialize_converters() {
	GString *string;
	gchar *text;
	int count_ext = 0;

	/* internal converters */
	supported_filetypes |= TYPEBMP;
	string = g_string_new(_("BMP images, "));
	supported_filetypes |= TYPEPIC;
	string = g_string_append(string, _("PIC images (IRIS), "));
	supported_filetypes |= TYPEPNM;
	string = g_string_append(string, _("PGM and PPM binary images"));

	supported_extensions[count_ext++] = ".fit";
	supported_extensions[count_ext++] = ".fits";
	supported_extensions[count_ext++] = ".fts";
	supported_extensions[count_ext++] = ".fit.fz";
	supported_extensions[count_ext++] = ".fits.fz";
	supported_extensions[count_ext++] = ".fts.fz";
	supported_extensions[count_ext++] = ".bmp";
	supported_extensions[count_ext++] = ".ppm";
	supported_extensions[count_ext++] = ".pgm";
	supported_extensions[count_ext++] = ".pnm";
	supported_extensions[count_ext++] = ".pic";


#ifdef HAVE_LIBRAW
	int i, nb_raw;

	supported_filetypes |= TYPERAW;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("RAW images"));

	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		// add the '.' to the extension
		char *ext = malloc(strlen(supported_raw[i].extension) + 2 * sizeof(char));
		sprintf(ext, ".%s", supported_raw[i].extension);
		supported_extensions[count_ext+i] = ext;
	}
	count_ext += nb_raw;
#endif
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("FITS-CFA images"));

#ifdef HAVE_FFMS2
	supported_filetypes |= TYPEAVI;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("Films"));
#endif

	supported_filetypes |= TYPESER;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("SER sequences"));

	/* library converters (detected by configure) */
#ifdef HAVE_LIBTIFF
	supported_filetypes |= TYPETIFF;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("TIFF images"));
	supported_extensions[count_ext++] = ".tif";
	supported_extensions[count_ext++] = ".tiff";
#endif

#ifdef HAVE_LIBXISF
	supported_filetypes |= TYPEXISF;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("XISF images"));
	supported_extensions[count_ext++] = ".xisf";
#endif

#ifdef HAVE_LIBJPEG
	supported_filetypes |= TYPEJPG;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("JPG images"));
	supported_extensions[count_ext++] = ".jpg";
	supported_extensions[count_ext++] = ".jpeg";
#endif

#ifdef HAVE_LIBJXL
	supported_filetypes |= TYPEJXL;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("JPEG XL images"));
	supported_extensions[count_ext++] = ".jxl";
#endif

#ifdef HAVE_LIBPNG
	supported_filetypes |= TYPEPNG;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("PNG images"));
	supported_extensions[count_ext++] = ".png";
#endif

#ifdef HAVE_LIBHEIF
	supported_filetypes |= TYPEHEIF;
	supported_filetypes |= TYPEAVIF;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("HEIF images, AVIF images"));
	supported_extensions[count_ext++] = ".heic";
	supported_extensions[count_ext++] = ".heif";
	supported_extensions[count_ext++] = ".avif";
#endif
	supported_extensions[count_ext++] = NULL;		// NULL-terminated array

	string = g_string_append(string, ".");
	text = g_string_free(string, FALSE);

	siril_log_message(_("Supported file types: %s\n"), text);
	return text;
}

/******************************************************************************
 *                                                                            *
 *            ALL CODE BELOW IS RELATED TO THE CONVERSION PROCESS             *
 *                                                                            *
 * The conversion takes as input a list of files and some conversion options. *
 * The process is based on a reader, that opens and reads frames from input   *
 * files, and a writer, that writes the frames to the output file format.     *
 * All reads and write are programmed in the main thread in                   *
 * convert_thread_worker(), then they are executed by a thread pool in        *
 * pool_worker().                                                             *
 *****************************************************************************/

static int check_for_raw_extensions(const char *extension) {
	int i, nb_raw;
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		if (!g_ascii_strcasecmp(extension, supported_raw[i].extension))
			return 0;
	}
	return 1;
}

/* returns the image_type for the extension without the dot, only if it is supported by
 * the current instance of Siril. */
image_type get_type_for_extension(const char *extension) {
	if (!extension) return TYPEUNDEF;
	if ((supported_filetypes & TYPEBMP) && !g_ascii_strcasecmp(extension, "bmp")) {
		return TYPEBMP;
	} else if ((supported_filetypes & TYPEJPG) &&
			(!g_ascii_strcasecmp(extension, "jpg") || !g_ascii_strcasecmp(extension, "jpeg"))) {
		return TYPEJPG;
	} else if ((supported_filetypes & TYPEJXL) &&
			(!g_ascii_strcasecmp(extension, "jxl"))) {
		return TYPEJXL;
	} else if ((supported_filetypes & TYPEHEIF) &&
			(!g_ascii_strcasecmp(extension, "heic") || !g_ascii_strcasecmp(extension, "heif"))) {
		return TYPEHEIF;
	} else if ((supported_filetypes & TYPEAVIF) &&
			(!g_ascii_strcasecmp(extension, "avif"))) {
		return TYPEAVIF;
	} else if ((supported_filetypes & TYPETIFF) &&
			(!g_ascii_strcasecmp(extension, "tif") || !g_ascii_strcasecmp(extension, "tiff"))) {
		return TYPETIFF;
	} else if ((supported_filetypes & TYPEPNG) && !g_ascii_strcasecmp(extension, "png")) {
		return TYPEPNG;
	} else if ((supported_filetypes & TYPEPNM) &&
			(!g_ascii_strcasecmp(extension, "pnm") || !g_ascii_strcasecmp(extension, "ppm") ||
			 !g_ascii_strcasecmp(extension, "pgm"))) {
		return TYPEPNM;
	} else if ((supported_filetypes & TYPEPIC) && !g_ascii_strcasecmp(extension, "pic")){
		return TYPEPIC;
	} else if ((supported_filetypes & TYPERAW) && !check_for_raw_extensions(extension)) {
		return TYPERAW;
#ifdef HAVE_FFMS2
	} else if ((supported_filetypes & TYPEAVI) && !check_for_film_extensions(extension)) {
		return TYPEAVI;
#endif
	} else if ((supported_filetypes & TYPESER) && !g_ascii_strcasecmp(extension, "ser")) {
		return TYPESER;
	} else if (!g_ascii_strcasecmp(extension, "fit") || !g_ascii_strcasecmp(extension, "fits") ||
			!g_ascii_strcasecmp(extension, "fts")) {
		return TYPEFITS;
	} else if (!g_ascii_strcasecmp(extension, "fit.fz") || !g_ascii_strcasecmp(extension, "fits.fz") ||
			!g_ascii_strcasecmp(extension, "fts.fz")) {
		return TYPEFITS;
	} else if (!g_ascii_strcasecmp(extension, "xisf")) {
		return TYPEXISF;
	}
	return TYPEUNDEF; // not recognized or not supported
}

/* open the file with path source from any image type and load it into the given FITS object */
int any_to_fits(image_type imagetype, const char *source, fits *dest,
		gboolean interactive, gboolean force_float) {
	int retval = 0;

	switch (imagetype) {
		case TYPEFITS:
			retval = (readfits(source, dest, NULL, force_float) != 0);
			break;
		case TYPEBMP:
			retval = (readbmp(source, dest) < 0);
			break;
		case TYPEPIC:
			retval = (readpic(source, dest) < 0);
			break;
#ifdef HAVE_LIBTIFF
		case TYPETIFF:
			retval = (readtif(source, dest, force_float, TRUE) < 0);
			break;
#endif
#ifdef HAVE_LIBXISF
		case TYPEXISF:
			retval = (readxisf(source, dest, force_float) < 0);
			break;
#endif
		case TYPEPNM:
			retval = (import_pnm_to_fits(source, dest) < 0);
			break;
#ifdef HAVE_LIBJPEG
		case TYPEJPG:
			retval = (readjpg(source, dest) < 0);
			break;
#endif
#ifdef HAVE_LIBJXL
		case TYPEJXL:
			retval = (readjxl(source, dest) < 0);
			break;
#endif
#ifdef HAVE_LIBHEIF
		case TYPEHEIF:
		case TYPEAVIF:
			/* need to retrieve all return values */
			retval = readheif(source, dest, interactive);
			break;
#endif
#ifdef HAVE_LIBPNG
		case TYPEPNG:
			retval = (readpng(source, dest) < 0);
			break;
#endif
#ifdef HAVE_LIBRAW
		case TYPERAW:
			{
				const char *src = source;
#ifdef _WIN32
				char *rsrc = siril_real_path(source);
				if (rsrc != NULL) {
					src  = rsrc;
				}
#endif
				retval = (open_raw_files(src , dest) < 0);
#ifdef _WIN32
				if (rsrc != NULL) {
					g_free(rsrc);
				}
#endif
			}
			break;
#endif
		case TYPESER:
		case TYPEAVI:
			siril_log_message(_("Requested converting a sequence file to single FITS image, should not happen\n"));
			retval = 1;
			break;
		case TYPEUNDEF:
		default:	// when the ifdefs are not compiled, default happens!
			siril_log_message(_("Error opening %s: file type not supported.\n"), source);
			retval = 1;
	}

	return retval;
}

#ifdef HAVE_FFMS2
int convert_single_film_to_ser(sequence *seq) {
	char **files_to_convert = malloc(1 * sizeof(char *));

	if (!files_to_convert) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	files_to_convert[0] = g_strdup(seq->film_file->filename);

	struct _convert_data *args = calloc(1, sizeof(struct _convert_data));
	args->start = 1;
	args->list = files_to_convert;
	args->total = 1;
	gchar *basename = remove_ext_from_filename(seq->film_file->filename);
	args->destroot = g_strdup_printf("%s_converted.ser", basename);
	g_free(basename);
	args->input_has_a_seq = TRUE;
	args->input_has_a_film = TRUE;
	args->debayer = FALSE;
	args->output_type = SEQ_SER;
	args->multiple_output = FALSE;
	args->make_link = FALSE;
	gettimeofday(&(args->t_start), NULL);
	if (!start_in_new_thread(convert_thread_worker, args)) {
		g_strfreev(args->list);
		g_free(args->destroot);
		free(args);
		return 1;
	}
	return 0;
}
#endif

typedef enum {
	/* for the reader data provider */
	GOT_READ_ERROR,
	GOT_OK_SEQ,
	GOT_OK_LAST_IN_SEQ,
	GOT_OK_FILE,
	GOT_OK_LAST_FILE,
	GOT_OK_LAST_IN_SEQ_LAST_FILE,

	/* for next file opening */
	OPEN_ERROR,
	OPEN_ERROR_AND_STOP,
	OPEN_OK,
	OPEN_SEQ,
	OPEN_NOT_A_SEQ,

	/* for the actual reader */
	NOT_READ,
	READ_FAILED,
	CAN_BE_LINKED,
	READ_OK
} seqread_status;

typedef enum {
	GOT_WRITE_ERROR,
	GOT_OK_WRITE,
	WRITE_FAILED,
	WRITE_OK
} seqwrite_status;

struct reader_seq_counter {
	gint count;
	gint close_sequence_after_read; // used as boolean
};

struct writer_seq_counter {
	gint count;
	gint close_sequence_after_write; // used as boolean
};

/* single image reading information */
struct reader_data {
	/* read from an opened sequence, whichever is not NULL */
	struct ser_struct *ser;
	fitseq *fitseq;
	void **threads; // for fitseq read
	int nb_threads;	// size of threads
#ifdef HAVE_FFMS2
	struct film_struct *film;
#endif
	int index;
	struct reader_seq_counter *seq_count;

	/* or read from a file to open */
	char *filename;
	// allow_link means that the input and ouputs are a FITS file, no transformation needed
	gboolean allow_link;
	gboolean allow_32bits;

	gboolean debayer;
};

/* single image writing information */
struct writer_data {
	/* write to an opened sequence, whichever is not NULL */
	struct ser_struct *ser;
	fitseq *fitseq;
	int index;
	struct writer_seq_counter *seq_count;
	gboolean have_seqwriter;

	/* or write to a FITS files sequence */
	gchar *filename;

	gint *converted_files; // points to convert_status->converted_files
};

struct readwrite_data {
	struct reader_data *reader;
	struct writer_data *writer;
};

/* conversion internal state, current considered image */
typedef struct {
	/* input */
	struct _convert_data *args;
	struct ser_struct *current_ser;
	fitseq *current_fitseq;
#ifdef HAVE_FFMS2
	struct film_struct *current_film;
#endif
	int next_file;		// index in the list of input files
	int next_image_in_file;	// index in an input sequence
	struct reader_seq_counter *readseq_count;

	/* output */
	struct ser_struct *output_ser;
	fitseq *output_fitseq;
	int output_file_number;	// also serves as number of converted images
	int next_image_in_output;
	struct writer_seq_counter *writeseq_count;
	gboolean allow_link;
	gboolean allow_32bits;

	int number_of_threads;	// size of threads, size of pool
	void **threads; // for fitseq read

	/* counters, atomic access */
	gint nb_input_images;	// for reporting
	gint failed_images, converted_images, converted_files;
	gint fatal_error;	// used as boolean

	gint first;		// to count the reads for memory concerns
} convert_status;

static void pool_worker(gpointer data, gpointer user_data);
static void open_next_input_seq(convert_status *conv);
static seqread_status open_next_input_sequence(const char *src_filename, convert_status *convert, gboolean test_only);
static seqwrite_status open_next_output_seq(const struct _convert_data *args, convert_status *conv);
static seqread_status get_next_read_details(convert_status *conv, struct reader_data *reader);
static seqwrite_status get_next_write_details(struct _convert_data *args, convert_status *conv,
		struct writer_data *writer, gboolean end_of_input_seq, gboolean last_file_and_image);
static gchar *create_sequence_filename(sequence_type output_type, const char *destroot, int index);
static seqwrite_status write_image(fits *fit, struct writer_data *writer);
static void compute_nb_images_fit_mem(fits *fit, gboolean debayer, int *nb_threads, int *nb_images);
static void print_reader(struct reader_data *reader);
static void print_writer(struct writer_data *writer);

static void init_report(struct _convert_data *args);
static void report_file_conversion(struct _convert_data *args, struct readwrite_data *rwarg);
static void write_conversion_report(struct _convert_data *args);

// if _convert_data->update_GUI is TRUE,it also updates the sequences list
// i.e. the only GTK function in this file
static gboolean end_convert(gpointer p) {
	struct _convert_data *args = (struct _convert_data *) p;
	struct timeval t_end;

	if (!args->retval && get_thread_run() && args->nb_converted_files > 0) {
		// load the sequence unless it's in another directory
		char *converted_seqname = NULL;
		if (!string_is_a_path(args->destroot)) {
			if (args->output_type != SEQ_REGULAR) {
				int extidx = get_extension_index(args->destroot);
				if (extidx > 0 && extidx < INT_MAX - 5) {
					converted_seqname = malloc(extidx + 5);
					strncpy(converted_seqname, args->destroot, extidx);
					strcpy(converted_seqname+extidx, ".seq");
				} else {
					converted_seqname = strdup(args->destroot);
				}
			} else {
				converted_seqname = strdup(args->destroot);
			}
		}
		if (converted_seqname) {
			gboolean seqfilecreated = create_one_seq(args->destroot, args->output_type); // this forces creating the .seq file (#1519)
			if (!seqfilecreated) { // just a fallback
				check_seq();
			}
			if (args->update_GUI) {
				update_sequences_list(converted_seqname);
			}
		}
		free(converted_seqname);
	}

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	stop_processing_thread();
	free(args->destroot);
	free(args->report);
	free(args);
	return FALSE;
}

/* memory and threads management:
 * until we read the first image, we don't know how many of them can fit into memory.
 * conversion uses a thread pool. It's initialized to 1 thread only, after the first image is read, it's
 * resized to how many images can fit in memory, also depending on the image writer allocation.
 * The writer is also allocated a size of 1 at first, resized after the first read.
 * Reading or writing to film formats doesn't support multi-threading.
 */
gpointer convert_thread_worker(gpointer p) {
	struct _convert_data *args = (struct _convert_data *) p;
	struct writer_data *writer = NULL;
	args->nb_converted_files = 0;
	args->retval = 0;
	gboolean allow_symlink = args->output_type == SEQ_REGULAR && test_if_symlink_is_ok(TRUE);
	args->make_link &= allow_symlink;
	if (args->multiple_output && args->output_type != SEQ_SER && args->output_type != SEQ_FITSEQ) {
		siril_log_message(_("disabling incompatible multiple output option in conversion\n"));
		args->multiple_output = FALSE;
	}
	if (args->output_type == SEQ_SER || args->output_type == SEQ_FITSEQ) {
		seqwriter_set_max_active_blocks(1);
	}
	set_progress_bar_data(_("Converting files"), PROGRESS_RESET);
	init_report(args);

	char *newdestroot = normalize_seqname(args->destroot, args->output_type == SEQ_REGULAR);
	free(args->destroot);
	args->destroot = strdup(newdestroot);
	gchar *seqname = g_strdup_printf("%s%s", args->destroot, ".seq");
	if (g_unlink(seqname))
		siril_debug_print("Error in g_unlink()\n");
	if (args->output_type == SEQ_REGULAR) {
		// to make sure destroot has an extension (will be removed when creating the filenames)
		free(args->destroot);
		args->destroot = strdup(seqname);
	}
	g_free(seqname);
	free(newdestroot);

	convert_status convert = { 0 };
	convert.args = args;
	convert.output_file_number = args->start ? args->start : args->start + 1;
	convert.number_of_threads = com.max_thread;
	convert.threads = calloc(com.max_thread, sizeof(void *));
	convert.allow_link = args->make_link;
	convert.allow_32bits = args->output_type != SEQ_SER && !com.pref.force_16bit;

	g_mutex_init(&args->pool_mutex);
	g_cond_init(&args->pool_cond);
	args->pool = g_thread_pool_new(pool_worker, &convert, 1, FALSE, NULL);
	open_next_input_seq(&convert);
	open_next_output_seq(args, &convert);
	do {
		struct reader_data *reader = calloc(1, sizeof(struct reader_data));
		seqread_status rstatus = get_next_read_details(&convert, reader);
		if (rstatus == GOT_READ_ERROR) {
			siril_debug_print("got reader error\n");
			free(reader);
			break;
		}
		print_reader(reader);
		g_atomic_int_inc(&convert.nb_input_images);

		writer = calloc(1, sizeof(struct writer_data));
		seqwrite_status wstatus = get_next_write_details(args, &convert, writer,
				rstatus == GOT_OK_LAST_IN_SEQ || rstatus == GOT_OK_LAST_IN_SEQ_LAST_FILE,
				rstatus == GOT_OK_LAST_IN_SEQ_LAST_FILE || rstatus == GOT_OK_LAST_FILE);
		if (wstatus == GOT_WRITE_ERROR) {
			siril_debug_print("got writer error\n");
			free(reader);
			free(writer);
			writer = NULL;
			break;
		}
		print_writer(writer);
		struct readwrite_data *rwarg = malloc(sizeof(struct readwrite_data));
		rwarg->reader = reader;
		rwarg->writer = writer;
		report_file_conversion(args, rwarg);

		if (!g_thread_pool_push(args->pool, rwarg, NULL)) {
			siril_log_message(_("Failed to queue image conversion task, aborting"));
			free(reader);
			free(writer);
			writer = NULL;
			free(rwarg);
			break;
		}
		if (rstatus == GOT_OK_LAST_FILE || rstatus == GOT_OK_LAST_IN_SEQ_LAST_FILE)
			break;
		if (rstatus == GOT_OK_LAST_IN_SEQ || rstatus == GOT_OK_FILE) {
			siril_debug_print("last image of the sequence reached, opening next sequence\n");
			open_next_input_seq(&convert);
		}
		// reader is freed elsewhere
	} while (com.run_thread);
	siril_debug_print("conversion scheduling loop finished, waiting for first read task to signal\n");
	g_mutex_lock(&args->pool_mutex);
	while (convert.first <= 0) {
		g_cond_wait(&args->pool_cond, &args->pool_mutex);
	}
	g_mutex_unlock(&args->pool_mutex);
	siril_debug_print("waiting for conversion tasks to finish\n");
	// this cannot be called before the g_thread_pool_set_max_threads() call, hence the pool_cond
	g_thread_pool_free(args->pool, FALSE, TRUE);
	siril_debug_print("conversion tasks finished\n");
	g_mutex_clear(&args->pool_mutex);
	g_cond_clear(&args->pool_cond);

	/* clean-up and reporting */
	g_strfreev(args->list);
	args->nb_converted_files = convert.converted_files;
	if (args->output_type == SEQ_REGULAR) {
		if (convert.fatal_error)
			siril_log_message(_("Conversion ended with error, %d/%d input files converted\n"), args->nb_converted_files, args->total);
		else {
			if (convert.nb_input_images == convert.converted_images)
				siril_log_message(_("Conversion succeeded, %d file(s) created for %d input file(s) (%d image(s) converted, %d failed)\n"), args->nb_converted_files, args->total, convert.converted_images, convert.failed_images);
			else siril_log_message(_("Conversion aborted, %d file(s) created for %d input file(s) (%d image(s) converted, %d failed)\n"), args->nb_converted_files, args->total, convert.converted_images, convert.failed_images);
			write_conversion_report(args);
		}
	} else {
		if (convert.fatal_error)
			siril_log_message(_("Conversion ended with error, %d output files created\n"), args->nb_converted_files);
		else {
			gboolean success = TRUE;
			if (!args->multiple_output && args->nb_converted_files == 1)
				siril_log_message(_("Conversion succeeded, %d file(s) created for %d input file(s) (%d image(s) converted, %d failed)\n"), args->nb_converted_files, args->total, convert.converted_images, convert.failed_images);
			else if (args->multiple_output && convert.nb_input_images == args->nb_converted_files)
				siril_log_message(_("Conversion succeeded, %d file(s) created for %d input file(s)\n"), args->nb_converted_files, args->total);
			else {
				siril_log_message(_("Conversion aborted, %d file(s) created for %d input file(s)\n"), args->nb_converted_files, args->total);
				success = FALSE;
			}
			if (success)
				write_conversion_report(args);
		}
	}
	// TODO still need to understand why, in case of error while writing a frame,
	// sometimes we do create the file and sometimes the error is caught and we get convert.fatal_error to 1...
	free(convert.output_fitseq);
	free(convert.output_ser);
	args->update_GUI = TRUE;
	if (!siril_add_idle(end_convert, args)) {
		args->update_GUI = FALSE;
		end_convert(args);
	}

	return NULL;
}

static struct reader_seq_counter *get_new_read_counter() {
	struct reader_seq_counter *counter = malloc(sizeof(struct reader_seq_counter));
	counter->count = 0;
	counter->close_sequence_after_read = 0;
	return counter;
}

static struct writer_seq_counter *get_new_write_counter() {
	struct writer_seq_counter *counter = malloc(sizeof(struct writer_seq_counter));
	counter->count = 0;
	counter->close_sequence_after_write = 0;
	return counter;
}

static void count_reader(struct reader_seq_counter *counter, gboolean close_after_read) {
	g_atomic_int_inc(&counter->count);
	g_atomic_int_set(&counter->close_sequence_after_read, (int)close_after_read);
}

static void count_writer(struct writer_seq_counter *counter, gboolean close_after_write) {
	g_atomic_int_inc(&counter->count);
	g_atomic_int_set(&counter->close_sequence_after_write, (int)close_after_write);
}

static void finish_read_seq(struct reader_data *reader) {
	if (!reader->seq_count) return;
	gboolean last = g_atomic_int_dec_and_test(&reader->seq_count->count);
	if (last && g_atomic_int_get(&reader->seq_count->close_sequence_after_read)) {
		if (reader->ser) {
			siril_debug_print("Closing input SER sequence %s\n", reader->ser->filename);
			ser_close_file(reader->ser);
		}
		else if (reader->fitseq) {
			siril_debug_print("Closing input FITS sequence file %s\n", reader->fitseq->filename);
			fitseq_close_file(reader->fitseq);
		}
#ifdef HAVE_FFMS2
		else if (reader->film) {
			siril_debug_print("Closing input film %s\n", reader->film->filename);
			film_close_file(reader->film);
		}
#endif
		free(reader->seq_count);
	}
}

static void finish_write_seq(struct writer_data *writer, gboolean success) {
	if (!writer->seq_count) return;
	gboolean last = g_atomic_int_dec_and_test(&writer->seq_count->count);
	if (last && g_atomic_int_get(&writer->seq_count->close_sequence_after_write)) {
		if (writer->ser) {
			siril_debug_print("closing write SER sequence%s\n", success ? "" : " and deleting the file");
			if (success) {
				if(!ser_write_and_close(writer->ser))
					g_atomic_int_inc(writer->converted_files);
			}
			else ser_close_and_delete_file(writer->ser);
		}
		else if (writer->fitseq) {
			siril_debug_print("closing write FITS sequence%s\n", success ? "" : " and deleting the file");
			if (success) {
				if (!fitseq_close_file(writer->fitseq))
					g_atomic_int_inc(writer->converted_files);
			}
			else fitseq_close_and_delete_file(writer->fitseq);
		}
		free(writer->seq_count);
	}
}

/* the reader part, reader arg must be zeroed */
static seqread_status get_next_read_details(convert_status *conv, struct reader_data *reader) {
	seqread_status retval = GOT_READ_ERROR;
	reader->debayer = conv->args->debayer;

	// first, check for already opened sequence files
	if (conv->current_ser) {
		reader->ser = conv->current_ser;
		reader->index = conv->next_image_in_file++;
		reader->seq_count = conv->readseq_count;
		count_reader(reader->seq_count, conv->next_image_in_file == conv->current_ser->frame_count);
		if (conv->next_image_in_file == conv->current_ser->frame_count) {
			retval = conv->next_file == conv->args->total ?
				GOT_OK_LAST_IN_SEQ_LAST_FILE : GOT_OK_LAST_IN_SEQ;
			conv->current_ser = NULL;
			conv->next_image_in_file = 0;
		}
		else retval = GOT_OK_SEQ;
	} else if (conv->current_fitseq) {
		reader->fitseq = conv->current_fitseq;
		reader->index = conv->next_image_in_file++;
		reader->threads = conv->threads;
		reader->nb_threads = conv->number_of_threads;
		reader->seq_count = conv->readseq_count;
		count_reader(reader->seq_count, conv->next_image_in_file == conv->current_fitseq->frame_count);
		if (conv->next_image_in_file == conv->current_fitseq->frame_count) {
			retval = conv->next_file == conv->args->total ?
				GOT_OK_LAST_IN_SEQ_LAST_FILE : GOT_OK_LAST_IN_SEQ;
			conv->current_fitseq = NULL;
			conv->next_image_in_file = 0;
		}
		else retval = GOT_OK_SEQ;
	}
#ifdef HAVE_FFMS2
	else if (conv->current_film) {
		reader->film = conv->current_film;
		reader->index = conv->next_image_in_file++;
		reader->seq_count = conv->readseq_count;
		count_reader(reader->seq_count, conv->next_image_in_file == conv->current_film->frame_count);
		if (conv->next_image_in_file == conv->current_film->frame_count) {
			retval = conv->next_file == conv->args->total ?
				GOT_OK_LAST_IN_SEQ_LAST_FILE : GOT_OK_LAST_IN_SEQ;
			conv->current_film = NULL;
			conv->next_image_in_file = 0;
		}
		else retval = GOT_OK_SEQ;
	}
#endif
	else {
		// else, read next file
		seqread_status next_status;
		do {
			char *filename = conv->args->list[conv->next_file++];
			// it should not be a sequence here
			next_status = open_next_input_sequence(filename, conv, TRUE);
			if (next_status == OPEN_NOT_A_SEQ) {
				if (conv->args->multiple_output) {
					siril_log_message(_("Ignoring an image file '%s' in the inputs as multiple outputs is enabled\n"), filename);
					retval = OPEN_ERROR;
					break;
				}
				reader->filename = filename;
				reader->allow_link = conv->allow_link;
				reader->allow_32bits = conv->allow_32bits;
				if (conv->next_file == conv->args->total)
					retval = GOT_OK_LAST_FILE;
				else retval = GOT_OK_FILE;
			}
			else siril_log_message(_("Skipping input file %s (failed to be opened)\n"), filename);
		} while (next_status == OPEN_ERROR && conv->next_file < conv->args->total);
	}
	return retval;
}

static void open_next_input_seq(convert_status *conv) {
	seqread_status status;
	do {
		const char *filename = conv->args->list[conv->next_file];
		status = open_next_input_sequence(filename, conv, FALSE);
		if (status == OPEN_ERROR || status == OPEN_ERROR_AND_STOP) {
			siril_log_color_message(_("File %s was not recognised as readable by Siril, skipping\n"), "salmon", filename);
			g_atomic_int_inc(&conv->failed_images);
			g_atomic_int_set(&conv->fatal_error, 1);
			if (status == OPEN_ERROR_AND_STOP) break;
		}
		else if (status == OPEN_OK) {
			conv->next_file++;
		}
		// else, it's OPEN_NOT_A_SEQ and we don't do anything here
	} while (status == OPEN_ERROR);
}

/* open the file with path source from any image type and load it into a new FITS object */
static fits *any_to_new_fits(image_type imagetype, const char *source, gboolean force_debayer, gboolean allow_32bits) {
	int retval = 0;
	fits *tmpfit = calloc(1, sizeof(fits));
	retval = any_to_fits(imagetype, source, tmpfit, FALSE, FALSE);

	if (!retval) {
		if (!allow_32bits && tmpfit->type == DATA_FLOAT) {
			siril_log_color_message(_("Converting 32 bits images (from %s) to 16 bits is not supported, ignoring file.\n"), "salmon", source);
			retval = 1;
		}
		else retval = debayer_if_needed(imagetype, tmpfit, force_debayer);
	}

	if (retval) {
		clearfits(tmpfit);
		free(tmpfit);
		return NULL;
	}

	return tmpfit;
}

static int get_thread_id(struct reader_data *reader) {
	void *self = g_thread_self();
	for (int i = 0; i < reader->nb_threads; i++) {
		if (!reader->threads[i])
			reader->threads[i] = self;
		if (reader->threads[i] == self)
			return i;
	}
	siril_debug_print("ERROR: could not find thread for fitseq reading\n");
	return -1;
}

static fits *read_fit(struct reader_data *reader, seqread_status *retval) {	// reentrant
	fits *fit = NULL;
	if (reader->ser) {
		fit = calloc(1, sizeof(fits));
		if (ser_read_frame(reader->ser, reader->index, fit, FALSE, com.pref.debayer.open_debayer))
			*retval = READ_FAILED;
		else *retval = READ_OK;
		finish_read_seq(reader);
	} else if (reader->fitseq) {
		fit = calloc(1, sizeof(fits));
		if (fitseq_read_frame(reader->fitseq, reader->index, fit, FALSE, get_thread_id(reader)))
			*retval = READ_FAILED;
		else {
			debayer_if_needed(TYPEFITS, fit, FALSE);
			*retval = READ_OK;
		}

		finish_read_seq(reader);
	}
#ifdef HAVE_FFMS2
	else if (reader->film) {
		fit = calloc(1, sizeof(fits));
		if (film_read_frame(reader->film, reader->index, fit) != FILM_SUCCESS)
			*retval = READ_FAILED;
		else *retval = READ_OK;
		finish_read_seq(reader);
	}
#endif
	else if (reader->filename) {
		const char *src_ext = get_filename_ext(reader->filename);
		image_type imagetype = get_type_for_extension(src_ext);
		if (imagetype == TYPEFITS && reader->allow_link) {
			*retval = CAN_BE_LINKED;
			return NULL;	// do not free reader, we need it for links
		} else {
			fit = any_to_new_fits(imagetype, reader->filename,
					reader->debayer, reader->allow_32bits);
			*retval = fit ? READ_OK : READ_FAILED;
			if (*retval == READ_OK) {
				/* Copy original filename in the header. Truncate if needed */
				/** 65 is the maximum length taking into account total length of card with, keyword and spaces and without comment */
				copy_filename(g_filename_display_basename(reader->filename), fit->keywords.filename, 65);
			}
		}
	}
	else *retval = NOT_READ;
	free(reader);
	return fit;
}

static int make_link(struct readwrite_data *rwdata) {
	int retval = 1;
	siril_debug_print("making link: %s -> %s\n", rwdata->reader->filename, rwdata->writer->filename);
	if (rwdata->writer->filename) {
		if (!symlink_uniq_file(rwdata->reader->filename, rwdata->writer->filename, TRUE))
			retval = 0;
		free(rwdata->reader);
	}
	return retval;
}

// to be called to unblock the threads if readjust_memory_limits cannot be called
static void signal_memory_limit(convert_status *conv) {
	if (g_atomic_int_add(&conv->first, 1))
		return;
	siril_debug_print("unblocking the main thread after conversion error\n");
	g_mutex_lock(&conv->args->pool_mutex);
	g_cond_signal(&conv->args->pool_cond);
	g_mutex_unlock(&conv->args->pool_mutex);
}

// see the memory and threads notes higher in the file
static void readjust_memory_limits(convert_status *conv, fits *fit) {
	if (g_atomic_int_add(&conv->first, 1))
		return;
	if (conv->args->input_has_a_film)
		goto unlock_end;
	int nb_threads, nb_images;
	compute_nb_images_fit_mem(fit, conv->args->debayer, &nb_threads, &nb_images);
	if (nb_threads <= 0)
		goto unlock_end;
	siril_log_message("%d image(s) can be processed in parallel\n", nb_threads);
	g_thread_pool_set_max_threads(conv->args->pool, nb_threads, NULL);
	if (conv->args->output_type == SEQ_SER || conv->args->output_type == SEQ_FITSEQ)
		seqwriter_set_max_active_blocks(nb_images);
unlock_end:
	g_mutex_lock(&conv->args->pool_mutex);
	g_cond_signal(&conv->args->pool_cond);
	g_mutex_unlock(&conv->args->pool_mutex);
}

static void handle_error(struct readwrite_data *rwdata) {
	siril_debug_print("conversion aborted or failed, cancelling this thread\n");
	seqwriter_release_memory();
	if (rwdata->reader) {
		finish_read_seq(rwdata->reader);
		free(rwdata->reader);
	}
	finish_write_seq(rwdata->writer, FALSE);
	free(rwdata->writer);
	free(rwdata);
}

static void pool_worker(gpointer data, gpointer user_data) {
	struct readwrite_data *rwdata = (struct readwrite_data *)data;
	convert_status *conv = (convert_status *)user_data;
	seqread_status read_status;
	if (rwdata->writer->have_seqwriter)
		seqwriter_wait_for_memory();

	if (!get_thread_run() || g_atomic_int_get(&conv->fatal_error)) {
		handle_error(rwdata);
		signal_memory_limit(conv);
		return;
	}

	fits *fit = read_fit(rwdata->reader, &read_status);
	if (read_status == CAN_BE_LINKED) {
		if (make_link(rwdata)) {
			g_atomic_int_inc(&conv->failed_images);
			g_atomic_int_set(&conv->fatal_error, 1);
		}
		else {
			g_atomic_int_inc(&conv->converted_images);
			g_atomic_int_inc(&conv->converted_files);
		}
		free(rwdata);	// reader and writer are freed in their function
		if (fit) {
			clearfits(fit);
			free(fit);
		}
		signal_memory_limit(conv);
		return;
	}
	else if (!fit || read_status == NOT_READ || read_status == READ_FAILED) {
		siril_debug_print("read error, ignoring image\n");
		g_atomic_int_inc(&conv->failed_images);
		finish_write_seq(rwdata->writer, FALSE);
		if (rwdata->writer->have_seqwriter)
			seqwriter_release_memory();
		free(rwdata);
		if (fit) {
			clearfits(fit);
			free(fit);
		}
		signal_memory_limit(conv);
		return;
	}
	readjust_memory_limits(conv, fit);

	if (!get_thread_run() || g_atomic_int_get(&conv->fatal_error)) {
		rwdata->reader = NULL;
		clearfits(fit);
		free(fit);
		if (rwdata->writer->have_seqwriter)
			handle_error(rwdata);
		return;
	}

	seqwrite_status write_status = write_image(fit, rwdata->writer);
	// clearfits is managed in write_image or sequence writer
	free(rwdata);	// reader and writer are freed in their function
	if (write_status == WRITE_FAILED) {
		g_atomic_int_inc(&conv->failed_images);
		g_atomic_int_set(&conv->fatal_error, 1);
	}
	else g_atomic_int_inc(&conv->converted_images);

	double percent = (double)g_atomic_int_get(&conv->converted_images) /
		(double)g_atomic_int_get(&conv->nb_input_images);
	set_progress_bar_data(NULL, percent);
	return;
}

static seqwrite_status get_next_write_details(struct _convert_data *args, convert_status *conv,
		struct writer_data *writer, gboolean end_of_input_seq, gboolean last_file_and_image) {
	writer->converted_files = &conv->converted_files;
	if (!args->multiple_output) {
		if (args->output_type == SEQ_SER) {
			if (!conv->output_ser) {
				conv->output_ser = calloc(1, sizeof(struct ser_struct));
				ser_init_struct(conv->output_ser);
				gchar *dest = g_str_has_suffix(args->destroot, ".ser") ? g_strdup(args->destroot) : g_strdup_printf("%s.ser", args->destroot);
				if (ser_create_file(dest, conv->output_ser, TRUE, NULL)) {
					siril_log_message(_("Creating the SER file `%s' failed, aborting.\n"), args->destroot);
					g_free(dest);
					return GOT_WRITE_ERROR;
				}
				g_free(dest);
				conv->next_image_in_output = 0;
				conv->writeseq_count = get_new_write_counter();
			}
			writer->ser = conv->output_ser;
			writer->index = conv->next_image_in_output++;
			writer->have_seqwriter = TRUE;
			writer->seq_count = conv->writeseq_count;
			count_writer(writer->seq_count, last_file_and_image);
			return GOT_OK_WRITE;
		}
		else if (args->output_type == SEQ_FITSEQ) {
			if (!conv->output_fitseq) {
				conv->output_fitseq = calloc(1, sizeof(struct fits_sequence));
				char *dest = g_str_has_suffix(args->destroot, com.pref.ext) ? args->destroot : g_strdup_printf("%s%s", args->destroot, com.pref.ext);
				if (fitseq_create_file(dest, conv->output_fitseq,
							args->input_has_a_seq ? -1 : args->total)) {
					siril_log_message(_("Creating the FITS sequence file `%s' failed, aborting.\n"), args->destroot);
					return GOT_WRITE_ERROR;
				}
				conv->next_image_in_output = 0;
				conv->writeseq_count = get_new_write_counter();
			}
			writer->fitseq = conv->output_fitseq;
			writer->index = conv->next_image_in_output++;
			writer->have_seqwriter = TRUE;
			writer->seq_count = conv->writeseq_count;
			count_writer(writer->seq_count, last_file_and_image);
			return GOT_OK_WRITE;
		}
		else {
			g_assert(args->output_type == SEQ_REGULAR);
			writer->filename = create_sequence_filename(SEQ_REGULAR, args->destroot, conv->output_file_number++);
			return GOT_OK_WRITE;
		}
	} else {
		if (args->output_type == SEQ_SER) {
			writer->ser = conv->output_ser;
			writer->index = conv->next_image_in_output++;
			writer->seq_count = conv->writeseq_count;
			count_writer(writer->seq_count, end_of_input_seq);
			if (end_of_input_seq && open_next_output_seq(args, conv) == GOT_WRITE_ERROR)
				return GOT_WRITE_ERROR;
			return GOT_OK_WRITE;
		}
		else if (args->output_type == SEQ_FITSEQ) {
			writer->fitseq = conv->output_fitseq;
			writer->index = conv->next_image_in_output++;
			writer->seq_count = conv->writeseq_count;
			count_writer(writer->seq_count, end_of_input_seq);
			if (end_of_input_seq && open_next_output_seq(args, conv) == GOT_WRITE_ERROR)
				return GOT_WRITE_ERROR;
			return GOT_OK_WRITE;
		}
	}
	return GOT_WRITE_ERROR;
}

static seqwrite_status open_next_output_seq(const struct _convert_data *args, convert_status *conv) {
	if (args->multiple_output) {
		if (args->output_type == SEQ_SER) {
			if (conv->next_file != conv->args->total) {
				gchar *dest_filename = create_sequence_filename(SEQ_SER, args->destroot, conv->output_file_number++);
				conv->output_ser = calloc(1, sizeof(struct ser_struct));
				ser_init_struct(conv->output_ser);
				if (ser_create_file(dest_filename, conv->output_ser, TRUE, NULL)) {
					siril_log_message(_("Creating the SER file `%s' failed, aborting.\n"), dest_filename);
					g_free(dest_filename);
					return GOT_WRITE_ERROR;
				}
				g_free(dest_filename);
				conv->next_image_in_output = 0;
				conv->writeseq_count = get_new_write_counter();
			}
		}
		else if (args->output_type == SEQ_FITSEQ) {
			if (conv->next_file != conv->args->total) {
				gchar *dest_filename = create_sequence_filename(SEQ_FITSEQ, args->destroot, conv->output_file_number++);
				conv->output_fitseq = calloc(1, sizeof(struct fits_sequence));
				if (fitseq_create_file(dest_filename, conv->output_fitseq, -1)) {
					siril_log_message(_("Creating the FITS sequence file `%s' failed, aborting.\n"), dest_filename);
					g_free(dest_filename);
					return GOT_WRITE_ERROR;
				}
				g_free(dest_filename);
				conv->next_image_in_output = 0;
				conv->writeseq_count = get_new_write_counter();
			}
		}
	}
	return GOT_OK_WRITE;
}

static seqread_status open_next_input_sequence(const char *src_filename, convert_status *convert, gboolean test_only) {
	const char *src_ext = get_filename_ext(src_filename);
	gchar *name = g_path_get_basename(src_filename);
	image_type imagetype = get_type_for_extension(src_ext);
	if (imagetype == TYPEUNDEF) {
		g_free(name);
		return OPEN_ERROR;
	}
#ifdef HAVE_FFMS2
	if (imagetype == TYPEAVI) {
		if (test_only) {
			g_free(name);
			return OPEN_SEQ;
		}
		if (convert->current_film) {
			siril_debug_print("error: opening a film while the previous was still here\n");
			g_free(name);
			return OPEN_ERROR;
		}
		convert->current_film = calloc(1, sizeof(struct film_struct));
		siril_log_message(_("Reading %s\n"), src_filename);
		if (film_open_file(src_filename, convert->current_film) != FILM_SUCCESS) {
			siril_log_message(_("Error while opening film %s, aborting.\n"), src_filename);
			free(convert->current_film);
			convert->current_film = NULL;
			g_free(name);
			return OPEN_ERROR;
		}
		convert->readseq_count = get_new_read_counter();
		g_free(name);
		return OPEN_OK;
	}
#endif
	else if (imagetype == TYPESER) {
		if (test_only) {
			g_free(name);
			return OPEN_SEQ;
		}
		if (convert->current_ser) {
			siril_debug_print("error: opening a SER while the previous was still here\n");
			g_free(name);
			return OPEN_ERROR;
		}
		convert->current_ser = malloc(sizeof(struct ser_struct)); // no need for calloc as init follows
		ser_init_struct(convert->current_ser);
		siril_log_message(_("Reading %s\n"), src_filename);
		if (ser_open_file(src_filename, convert->current_ser)) {
			siril_log_message(_("Error while opening ser file %s, aborting.\n"), src_filename);
			free(convert->current_ser);
			convert->current_ser = NULL;
			g_free(name);
			return OPEN_ERROR_AND_STOP;
		}
		convert->readseq_count = get_new_read_counter();
		g_free(name);
		return OPEN_OK;
	}
	else if (imagetype == TYPEFITS && fitseq_is_fitseq(name, NULL)) {
		if (test_only) return OPEN_SEQ;
		if (convert->current_fitseq) {
			siril_debug_print("error: opening a FITSEQ while the previous was still here\n");
			g_free(name);
			return OPEN_ERROR;
		}
		convert->current_fitseq = malloc(sizeof(fitseq)); // no need for calloc as init follows
		fitseq_init_struct(convert->current_fitseq);
		siril_log_message(_("Reading %s\n"), src_filename);
		if (fitseq_open(name, convert->current_fitseq, READONLY)) {
			siril_log_message(_("Error while opening ser file %s, ignoring file.\n"), src_filename);
			free(convert->current_fitseq);
			convert->current_fitseq = NULL;
			g_free(name);
			return OPEN_ERROR;
		}
		g_free(name);

		if (!convert->allow_32bits && get_data_type(convert->current_fitseq->bitpix) == DATA_FLOAT) {
			siril_log_color_message(_("Converting 32 bits images (from %s) to 16 bits is not supported, ignoring file.\n"), "salmon", src_filename);
			fitseq_close_file(convert->current_fitseq);
			free(convert->current_fitseq);
			convert->current_fitseq = NULL;
			return OPEN_ERROR_AND_STOP;
		}
		convert->readseq_count = get_new_read_counter();
		return OPEN_OK;
	}
	g_free(name);
	return OPEN_NOT_A_SEQ;
}

static gchar *create_sequence_filename(sequence_type output_type, const char *destroot, int index) {
	char *destroot_noext = remove_ext_from_filename(destroot);
	char dest_end = destroot_noext[strlen(destroot_noext)-1];
	gchar *output = NULL;
	gboolean append_underscore = dest_end != '_' && dest_end != '-' && (dest_end < '0' || dest_end > '9');
	switch (output_type) {
		case SEQ_REGULAR:
			output = g_strdup_printf("%s%s%05d%s", destroot_noext, append_underscore ? "_" : "", index, com.pref.ext);
			break;
		case SEQ_SER:
			output = g_strdup_printf("%s%s%03d.ser", destroot_noext, append_underscore ? "_" : "", index);
			break;
		case SEQ_FITSEQ:
			output = g_strdup_printf("%s%s%03d%s", destroot_noext, append_underscore ? "_" : "", index, com.pref.ext);
			break;
		default:
			siril_log_color_message(_("output sequence type unknown, aborting\n"), "red");
	}
	free(destroot_noext);
	return output;
}


/* the writer part */
static seqwrite_status write_image(fits *fit, struct writer_data *writer) {
	seqwrite_status retval = WRITE_FAILED;
	if (writer->ser) {
		if (ser_write_frame_from_fit(writer->ser, fit, writer->index)) {
			siril_log_color_message(_("Error while converting to SER (no space left?)\n"), "red");
		}
		else retval = WRITE_OK;
		finish_write_seq(writer, retval == WRITE_OK);
	}
	else if (writer->fitseq) {
		if (fitseq_write_image(writer->fitseq, fit, writer->index)) {
			siril_log_color_message(_("Error while converting to FITSEQ (no space left?)\n"), "red");
		}
		else retval = WRITE_OK;
		finish_write_seq(writer, retval == WRITE_OK);
	}
	else if (writer->filename) {
		if (savefits(writer->filename, fit)) {
			siril_log_color_message(_("Error while converting to FITS\n"), "red");
			/* We test if there's any space left. If not, we remove the image */
			if (!is_space_disk_available(com.wd)) {
				siril_log_color_message(_("No space left!\n"), "red");
				g_unlink(writer->filename);
			}
		}
		else {
			retval = WRITE_OK;
			g_atomic_int_inc(writer->converted_files);
		}
		clearfits(fit);
		free(fit);
		g_free(writer->filename);
	}
	else {
		siril_log_color_message(_("Error while converting, unknown output\n"), "red");
	}


	free(writer);
	return retval;
}

/* similar to compute_nb_images_fit_memory from sequence.c, but uses a FITS as input (same size as the
 * output files), not a sequence. It computes how many threads can be created and how many images fit in
 * memory with those threads, for the writer */
static void compute_nb_images_fit_mem(fits *fit, gboolean debayer, int *nb_threads, int *nb_images) {
	int max_memory_MB = get_max_memory_in_MB();
	/* image size only changes in case of debayer and in this case, it also needs
	 * more memory to do the debayer:
	 * 	O(3n) for ushort (the output image ushort and 1 float for processing)
	 * 	O(2n) for float (the output image and 1 for processing)
	 */
	uint64_t memory_per_output_image = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	uint64_t memory_per_processed_image;
	if (fit->type == DATA_FLOAT) {
		memory_per_output_image *= sizeof(float);
		memory_per_processed_image = memory_per_output_image;
		if (debayer)
			memory_per_processed_image *= 2;
	} else {
		memory_per_output_image *= sizeof(WORD);
		memory_per_processed_image = memory_per_output_image;
		if (debayer)
			memory_per_processed_image *= 3;
	}
	unsigned int memory_per_image_MB = memory_per_output_image / BYTES_IN_A_MB;
	if (memory_per_image_MB == 0)
		memory_per_image_MB = 1;
	max_memory_MB += memory_per_image_MB;	// one is already loaded
	unsigned int memory_per_thread_MB = memory_per_processed_image / BYTES_IN_A_MB;
	if (memory_per_thread_MB == 0)
		memory_per_thread_MB= 1;
	siril_debug_print("Memory required per image: %u MB. Max memory: %d MB\n", memory_per_thread_MB, max_memory_MB);
	*nb_threads = max_memory_MB / memory_per_thread_MB;
	if (*nb_threads >= com.max_thread)
		*nb_threads = com.max_thread;
	*nb_images = *nb_threads + (max_memory_MB - *nb_threads * memory_per_thread_MB) / memory_per_image_MB;
	if (*nb_images >= com.max_thread * 3)
		*nb_images = com.max_thread * 3;
}

static void print_reader(struct reader_data *reader) {
	if (reader->ser)
		siril_debug_print("I> reader: %s image %d%s%s\n", reader->ser->filename,
				reader->index,
				reader->seq_count->close_sequence_after_read ? " (close after read)" : "",
				reader->debayer ? " (debayer)" : "");
	else if (reader->fitseq)
		siril_debug_print("I> reader: %s image %d%s%s\n", reader->fitseq->filename,
				reader->index,
				reader->seq_count->close_sequence_after_read ? " (close after read)" : "",
				reader->debayer ? " (debayer)" : "");
#ifdef HAVE_FFMS2
	else if (reader->film)
		siril_debug_print("I> reader: %s image %d%s%s\n", reader->film->filename,
				reader->index,
				reader->seq_count->close_sequence_after_read ? " (close after read)" : "",
				reader->debayer ? " (debayer)" : "");
#endif
	else siril_debug_print("I> reader: %s%s\n", reader->filename,
				reader->debayer ? " (debayer)" : "");
}

static void print_writer(struct writer_data *writer) {
	if (writer->ser)
		siril_debug_print("O> writer: %s image %d%s\n", writer->ser->filename,
				writer->index,
				writer->seq_count->close_sequence_after_write ? " (close after write)" : "");
	else if (writer->fitseq)
		siril_debug_print("O> writer: %s image %d%s\n", writer->fitseq->filename,
				writer->index,
				writer->seq_count->close_sequence_after_write ? " (close after write)" : "");
	else siril_debug_print("O> writer: %s\n", writer->filename);
}

static void init_report(struct _convert_data *args) {
	args->report = malloc(args->total * sizeof(char *));
	args->report_length = 0;
}

static void report_file_conversion(struct _convert_data *args, struct readwrite_data *rwarg) {
	gchar *str = NULL;
	if (rwarg->reader->filename) {
		if (rwarg->writer->filename) {
			str = g_strdup_printf("'%s' -> '%s'\n", rwarg->reader->filename, rwarg->writer->filename);
		}
		else if (rwarg->writer->fitseq) {
			str = g_strdup_printf("'%s' -> '%s' image %d\n", rwarg->reader->filename, rwarg->writer->fitseq->filename, rwarg->writer->index);
		}
		else if (rwarg->writer->ser) {
			str = g_strdup_printf("'%s' -> '%s' image %d\n", rwarg->reader->filename, rwarg->writer->ser->filename, rwarg->writer->index);
		}
	}
	if (str) {
		args->report[args->report_length++] = str;
	}
}

static void write_conversion_report(struct _convert_data *args) {
	if (args->report_length <= 0)
		return;
	gchar *filename;
	char *filename_noext = remove_ext_from_filename(args->destroot);
	if (g_str_has_suffix(filename_noext, "_"))
		filename = g_strdup_printf("%sconversion.txt", filename_noext);
	else filename = g_strdup_printf("%s_conversion.txt", filename_noext);
	FILE *fd = g_fopen(filename, "w+");
	g_free(filename);
	free(filename_noext);
	if (!fd)
		return;

	for (int i = 0; i < args->report_length; i++)
		if (fputs(args->report[i], fd) == EOF)
			break;

	for (int i = 0; i < args->report_length; i++)
		g_free(args->report[i]);
	fclose(fd);
}

