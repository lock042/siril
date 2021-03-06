#ifndef CONVERSION_H
#define CONVERSION_H

#include <glib.h>
#include "core/siril.h" // for image_type

#define XTRANS_1 4
#define XTRANS_2 5

typedef struct {
	char *extension;			// name of the extension of raw
	char *manufacturer;			// name of the manufacturer
} supported_raw_list;

struct _convert_data {
	struct timeval t_start;
	gchar **list;
	int total;	// length of list
	int start;	// offset in output number
	gboolean input_has_a_seq;
	gboolean input_has_a_film;
	gboolean make_link;
	gchar *destroot;
	gboolean debayer;
	sequence_type output_type;
	gboolean multiple_output;	// multiple SER output

	int retval;
	int nb_converted_files;
};

#define MAX_EXTENSIONS 50	// actual size of supported_extensions

extern supported_raw_list supported_raw[];	//supported raw extensions
extern char *supported_extensions[MAX_EXTENSIONS];

int get_nb_raw_supported();
void list_format_available();
image_type get_type_for_extension(const char *extension);
gchar *initialize_converters();
gpointer convert_thread_worker(gpointer p);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean force_debayer);
int any_to_fits(image_type imagetype, const char *source, fits *dest, gboolean interactive, gboolean force_float, gboolean debayer);
#ifdef HAVE_FFMS2
int convert_single_film_to_ser(sequence *seq);
#endif

#endif
