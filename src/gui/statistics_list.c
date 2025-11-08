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


#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "io/sequence.h"
#include "io/single_image.h"

static GtkListStore *list_store = NULL;

static const char *first_colour[] = { "WhiteSmoke", "#1B1B1B" };
static const char *second_colour[] = { "Powder Blue", "#39394A" };

enum {
	COLUMN_NAME,		// string
	COLUMN_RVALUE,		// converted to string, not pure double
	COLUMN_GVALUE,		// converted to string, not pure double
	COLUMN_BVALUE,		// converted to string, not pure double
	COLUMN_COLOR,		// string
	N_COLUMNS
};

char *statName[] = {
		N_("mean"),
		N_("median"),
		N_("sigma"),
		N_("avgDev"),
		N_("MAD"),
		N_("sqrt(BWMV)"),
		N_("min"),
		N_("max"),
};

static void get_statlist_store() {
	if (list_store == NULL)
		list_store = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "liststoreStat"));
}


static void display_stat(const double *value, const double *normalization, char *format, int nblayer, int i, data_type type) {
	char rvalue[20], gvalue[20], bvalue[20];
	GtkTreeIter iter;
	int color = (com.pref.gui.combo_theme == 0) ? 1 : 0;

	sprintf(rvalue, format, value[RLAYER] / normalization[RLAYER]);

	if (nblayer > 1) {
		sprintf(gvalue, format, value[GLAYER] / normalization[GLAYER]);
		sprintf(bvalue, format, value[BLAYER] / normalization[BLAYER]);
	} else {
		sprintf(gvalue, "--");
		sprintf(bvalue, "--");
	}

	gtk_list_store_append(list_store, &iter);
	gtk_list_store_set(list_store, &iter, COLUMN_NAME, _(statName[i]),
			COLUMN_RVALUE, rvalue,
			COLUMN_GVALUE, gvalue,
			COLUMN_BVALUE, bvalue,
			COLUMN_COLOR, i % 2 ? second_colour[color] : first_colour[color],
			-1);
}

static double get_value_from_stat(imstats *stat, int index) {
	if (!stat) return 0.0;
	switch (index) {
	case 0:
		return stat->mean;
		break;
	case 1:
		return stat->median;
		break;
	case 2:
		return stat->sigma;
		break;
	case 3:
		return stat->avgDev;
		break;
	case 4:
		return stat->mad;
		break;
	case 5:
		return stat->sqrtbwmv;
		break;
	case 6:
		return stat->min;
		break;
	case 7:
		return stat->max;
		break;
	}
	return 0.0;
}

static void init_dialog() {
	static GtkTreeSelection *selection = NULL;
	static GtkWidget *cfacheck = NULL;
	get_statlist_store();
	if (!selection) {
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "treeview-selection9"));
		cfacheck = lookup_widget("cfastatsCheckButton");
	}

	gtk_list_store_clear(list_store);
	gtk_widget_set_sensitive(cfacheck, gfit.rx > 0 && gfit.naxes[2] == 1 && gfit.keywords.bayer_pattern[0] != '\0');
}

static void add_chan_stats_to_list(imstats **stat, int nblayer, data_type type, gboolean normalized) {
	double normalization[] = { 1.0, 1.0, 1.0 };
	char format[6];

	for (int i = 0; i < G_N_ELEMENTS(statName); i++) {
		double value[3] = { 0.0 };

		value[RLAYER] = get_value_from_stat(stat[RLAYER], i);
		if (nblayer > 1) {
			value[GLAYER] = get_value_from_stat(stat[GLAYER], i);
			value[BLAYER] = get_value_from_stat(stat[BLAYER], i);
		}

		if (normalized) {
			normalization[RLAYER] = stat[RLAYER] ? stat[RLAYER]->normValue : 0.0;
			if (nblayer > 1) {
				normalization[GLAYER] = stat[GLAYER] ? stat[GLAYER]->normValue : 0.0;
				normalization[BLAYER] = stat[BLAYER] ? stat[BLAYER]->normValue : 0.0;
			}
			if (value[RLAYER] < 1E-5 && value[RLAYER] > 0.0) {
				sprintf(format, "%%.5e");
			} else {
				sprintf(format, "%%.7lf");
			}
		} else {
			if (type == DATA_FLOAT) {
				/* by default it is shown in ushort mode */
				normalization[RLAYER] = 1.0 / USHRT_MAX_DOUBLE;
				if (nblayer > 1) {
					normalization[GLAYER] = 1.0 / USHRT_MAX_DOUBLE;
					normalization[BLAYER] = 1.0 / USHRT_MAX_DOUBLE;
				}
			}
			sprintf(format, "%%.1lf");
		}
		display_stat(value, normalization, format, nblayer, i, type);
	}
}

void on_statButtonClose_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("StatWindow");
}

void computeStat() {
	GtkToggleButton *checkButton, *cfaCheckButton;
	GtkLabel *statNameLabel, *statSelecLabel;
	gboolean normalized, use_cfa;
	int channel;
	gchar *name, *selection;

	checkButton = GTK_TOGGLE_BUTTON(lookup_widget("statCheckButton"));
	cfaCheckButton = GTK_TOGGLE_BUTTON(lookup_widget("cfastatsCheckButton"));
	statNameLabel = GTK_LABEL(lookup_widget("statNameLabel"));
	statSelecLabel = GTK_LABEL(lookup_widget("statSelecLabel"));
	normalized = gtk_toggle_button_get_active(checkButton);
	use_cfa = gtk_toggle_button_get_active(cfaCheckButton);

	if (single_image_is_loaded())
		name = g_strdup_printf("%s", com.uniq->filename);
	else if (sequence_is_loaded())
		name = g_strdup_printf(_("Image %d/%d from the sequence %s"),
				com.seq.current, com.seq.number + 1, com.seq.seqname);
	else
		name = g_strdup_printf(_("unknown image"));

	gtk_label_set_text(statNameLabel, name);
	g_free(name);

	if (com.selection.h && com.selection.w) {
		selection = g_strdup_printf(_("Size of selection in pixel: (%d,%d)"),
				com.selection.w, com.selection.h);
	} else {
		selection = g_strdup_printf(_("No selection"));
	}

	gtk_label_set_text(statSelecLabel, selection);
	g_free(selection);

	init_dialog();

	int nb_channels = (int)gfit.naxes[2];
	if (use_cfa) {
		if (nb_channels == 1 && gfit.keywords.bayer_pattern[0] != '\0')
			if ((com.selection.h || com.selection.w) && (com.selection.w < 2 || com.selection.h < 2))
				use_cfa = FALSE;
			else
				nb_channels = 3;
		else {
			use_cfa = FALSE;
		}
	}
	imstats *stat[3] = { NULL, NULL, NULL };
	for (channel = 0; channel < nb_channels; channel++) {
		int super_chan = channel;
		if (use_cfa)
			super_chan = -channel - 1;
		stat[channel] = statistics(NULL, -1, &gfit, super_chan, &com.selection, STATS_MAIN, MULTI_THREADED);
		if (!stat[channel]) {
			siril_log_message(_("Error: statistics computation failed.\n"));
		}
	}
	add_chan_stats_to_list(stat, nb_channels, gfit.type, normalized);

	for (channel = 0; channel < nb_channels; channel++) {
		free_stats(stat[channel]);
	}
}

void on_statCheckButton_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	set_cursor_waiting(FALSE);
}

void on_statButtonRun_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	set_cursor_waiting(FALSE);
}
