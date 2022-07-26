#include <gtk/gtk.h>
#include "algos/photometry.h"
#include "algos/siril_wcs.h"

gchar *open_nina_lc_dialog() {
	GtkWidget *dialog = gtk_dialog_new_with_buttons("Light curve with NINA star list", NULL,
			0, "_OK", GTK_RESPONSE_ACCEPT, "_Cancel", GTK_RESPONSE_REJECT, NULL);
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

	GtkWidget *button = gtk_file_chooser_button_new (_("Select the NINA star list file"),
			GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(button), com.wd);
	g_object_set(G_OBJECT(button), "margin", 15, NULL);
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("CSV file (*.csv)"));
	gtk_file_filter_add_pattern(f, "*.csv");
	gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(button), f);

	GtkWidget *label = gtk_label_new(_("Process a sequence to get a light curve on a star using the list of reference stars created by the NINA exoplanet plugin"));
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	g_object_set(G_OBJECT(label), "margin", 15, NULL);

	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_set_spacing(GTK_BOX(content_area), 20);
	gtk_container_add(GTK_CONTAINER(content_area), label);
	gtk_container_add(GTK_CONTAINER(content_area), button);
	gtk_widget_show_all(dialog);

	int retval = gtk_dialog_run(GTK_DIALOG(dialog));
	gtk_widget_hide(dialog);
	if (retval == GTK_RESPONSE_ACCEPT) {
		return gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(button));
	}
	return NULL;
}

void on_nina_button_clicked(GtkButton *button, gpointer user_data) {
	gchar *file_path = open_nina_lc_dialog();
	printf("returned %s\n", file_path);
	return;

#if 0
	if (!has_wcs(&gfit)) {
		// find reference image and load it
		return;
	}
	struct light_curve_args *args = malloc(sizeof(struct light_curve_args));
	if (parse_nina_stars_file_using_WCS(args, nina_file, &gfit)) {
		// fail
	}
	args->seq = &com.seq;
	args->layer = isrgb(&gfit) ? 1 : 0;

	start_in_new_thread(light_curve_worker, args);
#endif
}
