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
/* functions and variables defined only for linking */

#include "../core/siril.h"
#include "../core/pipe.h"

/* the global variables of the whole project (replacing main.c) */
cominfo com;	// the main data struct
guiinfo gui;	// the gui data struct
fits gfit;	// currently loaded image
char **supported_extensions;

gboolean sequence_is_loaded() {
        fprintf(stderr, "ERROR: calling undefined function sequence_is_loaded\n");
	return FALSE;
}

void free_sequence(sequence *seq, gboolean free_seq_too) {
        fprintf(stderr, "ERROR: calling undefined function free_sequence\n");
}

void close_sequence(int loading_another) {
        fprintf(stderr, "ERROR: calling undefined function close_sequence\n");
}

sequence * readseqfile(const char *name){
        fprintf(stderr, "ERROR: calling undefined function readseqfile\n");
	return NULL;
}

imstats* free_stats(imstats *stat) {
	fprintf(stderr, "ERROR: calling undefined function free_stats\n");
	return NULL;
}

void add_stats_to_fit(fits *fit, int layer, imstats *stat) {
	fprintf(stderr, "ERROR: calling undefined function add_stats_to_fit\n");
}

int pipe_send_message(pipe_message msgtype, pipe_verb verb, const char *arg) {
        fprintf(stderr, "ERROR: calling undefined function pipe_send_message\n");
	return 0;
}

double fit_get_min(fits *fit, int layer) {
        fprintf(stderr, "ERROR: calling undefined function fit_get_min\n");
	return 0.0;
}

double fit_get_max(fits *fit, int layer) {
        fprintf(stderr, "ERROR: calling undefined function fit_get_max\n");
	return 1.0;
}

int image_find_minmax(fits *fit) {
        fprintf(stderr, "ERROR: calling undefined function image_find_minmax\n");
	return 0;
}

image_type get_type_for_extension(const char *extension) {
        fprintf(stderr, "ERROR: calling undefined function get_type_for_extension\n");
	return TYPEFITS;
}

int single_image_is_loaded() {
        fprintf(stderr, "ERROR: calling undefined function single_image_is_loaded\n");
	return 0;
}

void set_GUI_MEM(guint64 used, const gchar *label) {
        fprintf(stderr, "ERROR: calling undefined function set_GUI_MEM\n");
}

void set_GUI_DiskSpace(gint64 space, const gchar *label) {
        fprintf(stderr, "ERROR: calling undefined function set_GUI_DiskSpace\n");
}

void gui_function(set_GUI_CWD, NULL) {
        fprintf(stderr, "ERROR: calling undefined function set_GUI_CWD\n");
}

void update_used_memory() {
        fprintf(stderr, "ERROR: calling undefined function update_used_memory\n");
}

void *lookup_widget(char *w) {
        fprintf(stderr, "ERROR: calling undefined function lookup_widget\n");
	return NULL;
}

void invalidate_stats_from_fit(fits *fit) {
        fprintf(stderr, "ERROR: calling undefined function invalidate_stats_from_fit\n");
}

const gchar* siril_get_system_data_dir() {
        fprintf(stderr, "ERROR: calling undefined function siril_get_system_data_dir\n");
	return NULL;
}

int fitseq_is_fitseq(const char *filename, int *frames) {
        fprintf(stderr, "ERROR: calling undefined function fitseq_is_fitseq\n");
	return 0;
}

void full_stats_invalidation_from_fit(fits *fit) {
        fprintf(stderr, "ERROR: calling undefined function full_stats_invalidation_from_fit\n");
}

gboolean load_WCS_from_fits(fits* fit) {
        fprintf(stderr, "ERROR: calling undefined function load_WCS_from_fits\n");
	return FALSE;
}


void free_wcs(fits *fit) {
        fprintf(stderr, "ERROR: calling undefined function free_wcs\n");
}

SirilWidget *siril_file_chooser_save(GtkWindow *parent, GtkFileChooserAction action) {
        fprintf(stderr, "ERROR: calling undefined function siril_file_chooser_save\n");
	return NULL;
}

gint siril_dialog_run(SirilWidget *widgetdialog) {
        fprintf(stderr, "ERROR: calling undefined function siril_dialog_run\n");
	return 0;
}

gboolean siril_confirm_dialog(gchar *title, gchar *msg) {
        fprintf(stderr, "ERROR: calling undefined function siril_confirm_dialog\n");
	return FALSE;
}

int process_clear(int nb) {
        fprintf(stderr, "ERROR: calling undefined function process_clear\n");
	return 0;
}

void siril_widget_destroy(SirilWidget *widgetdialog) {
        fprintf(stderr, "ERROR: calling undefined function siril_widget_destroy\n");
}

void gui_log_message(const char* msg, const char* color) {
        fprintf(stderr, "ERROR: calling undefined function gui_log_message\n");
}


void wcs_pc_to_cd(double pc[2][2], double cdelt[2], double cd[2][2]) {
        fprintf(stderr, "ERROR: calling undefined function wcs_pc_to_cd\n");
}

gboolean has_wcsdata(fits *fit) {
        fprintf(stderr, "ERROR: calling undefined function has_wcsdata\n");
	return FALSE;
}

double parse_hms(const char *objctra) {
        fprintf(stderr, "ERROR: calling undefined function parse_hms\n");
	return 0.0;
}

double parse_dms(const char *objctdec) {
        fprintf(stderr, "ERROR: calling undefined function parse_dms\n");
	return 0.0;
}

int check_threading(const threading_type *t) {
	if (*t != SINGLE_THREADED)
		return 2;
	return 1;
}

int limit_threading(threading_type *t, int set, int max) {
	return check_threading(t);
}

#ifdef DUMMY_LOG
char* siril_log_message(const char* format, ...) {
	va_list args;
	va_start(args, format);
	printf("\t");
	vprintf(format, args);
	va_end(args);
	return "";
}

char* siril_log_color_message(const char* format, const char* color, ...) {
	va_list args;
	va_start(args, color);
	printf("\t");
	vprintf(format, args);
	va_end(args);
	return "";
}
#endif
