#include "siril_file_chooser_button.h"

struct _SirilFileChooserButton {
    GtkButton parent_instance;
    GtkWidget *label;
    GtkWidget *mimetype_icon;
    GtkWidget *file_icon;
};

// Signal handler for file selection
static void on_file_selected(GtkFileChooserNative *native, gpointer user_data) {
    SirilFileChooserButton *file_chooser_button = SIRIL_FILE_CHOOSER_BUTTON(user_data);
    GtkFileChooser *chooser = GTK_FILE_CHOOSER(native);
    gchar *filename = gtk_file_chooser_get_filename(chooser);

    // Handle the selected file here
    g_print("Selected file: %s\n", filename);

    // Update the label, mimetype icon, and file icon
    const gchar *display_text = (filename != NULL) ? filename : "(None)";
    gtk_label_set_text(GTK_LABEL(file_chooser_button->label), display_text);

    GdkPixbuf *file_icon = (filename != NULL) ? gtk_icon_theme_load_icon(gtk_icon_theme_get_default(), "document-open", 16, GTK_ICON_LOOKUP_FORCE_SIZE, NULL) : NULL;
    gtk_image_set_from_pixbuf(GTK_IMAGE(file_chooser_button->file_icon), file_icon);

    gchar *mime_type = (filename != NULL) ? g_content_type_guess(filename, NULL, 0, NULL) : NULL;
    GdkPixbuf *mimetype_icon = (mime_type != NULL) ? gdk_pixbuf_new_from_file(g_content_type_get_icon(mime_type), NULL) : NULL;
    gtk_image_set_from_pixbuf(GTK_IMAGE(file_chooser_button->mimetype_icon), mimetype_icon);

    if (file_icon != NULL) {
        g_object_unref(file_icon);
    }
    if (mimetype_icon != NULL) {
        g_object_unref(mimetype_icon);
    }
    g_free(filename);
    g_free(mime_type);
}

// Callback function to open the file chooser
static void open_file_chooser(GtkButton *button, gpointer user_data) {
    SirilFileChooserButton *file_chooser_button = SIRIL_FILE_CHOOSER_BUTTON(user_data);

    // Retrieve parameters for gtk_file_chooser_native_new
    const gchar *title = "Open File";
    GtkWindow *parent = GTK_WINDOW(gtk_widget_get_toplevel(GTK_WIDGET(file_chooser_button)));
    GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
    const gchar *accept_label = "_Open";
    const gchar *cancel_label = "_Cancel";

    GtkFileChooserNative *native = gtk_file_chooser_native_new(title, parent, action, accept_label, cancel_label);

    // Set the current folder before running the dialog
    const gchar *current_folder = gtk_file_chooser_get_current_folder(GTK_FILE_CHOOSER(native));
    siril_file_chooser_set_current_folder(file_chooser_button, current_folder);

    g_signal_connect(native, "response", G_CALLBACK(on_file_selected), file_chooser_button);
    g_signal_connect(native, "destroy", G_CALLBACK(g_object_unref), NULL);

    gtk_native_dialog_run(GTK_NATIVE_DIALOG(native));
}

// Class initialization function
static void siril_file_chooser_button_class_init(SirilFileChooserButtonClass *klass) {
    GTK_BUTTON_CLASS(klass)->clicked = open_file_chooser;
}

// Instance initialization function
static void siril_file_chooser_button_init(SirilFileChooserButton *button) {
    button->label = gtk_label_new("(None)");
    button->mimetype_icon = gtk_image_new();
    button->file_icon = gtk_image_new();

    GtkWidget *box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_pack_start(GTK_BOX(box), button->mimetype_icon, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(box), button->label, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(box), button->file_icon, FALSE, FALSE, 0);

    gtk_container_add(GTK_CONTAINER(button), box);
}

// Function to set the maximum permissible width of the label
void siril_file_chooser_button_set_width_chars(SirilFileChooserButton *button, gint width_chars) {
    g_return_if_fail(SIRIL_IS_FILE_CHOOSER_BUTTON(button));
    GtkWidget *label = button->label;
    gtk_label_set_width_chars(GTK_LABEL(label), width_chars);
}

// Function to set the current folder in the file chooser
void siril_file_chooser_set_current_folder(SirilFileChooserButton *button, const gchar *folder_path) {
    g_return_if_fail(SIRIL_IS_FILE_CHOOSER_BUTTON(button));
    g_return_if_fail(folder_path != NULL);

    GtkFileChooser *file_chooser = GTK_FILE_CHOOSER(gtk_file_chooser_native_get_file_chooser(GTK_FILE_CHOOSER_NATIVE(button)));
    gtk_file_chooser_set_current_folder(file_chooser, folder_path);
}

// Function to create a new SirilFileChooserButton with parameters
SirilFileChooserButton *siril_file_chooser_button_new(const gchar *title, GtkFileChooserAction action) {
    return g_object_new(SIRIL_TYPE_FILE_CHOOSER_BUTTON,
                        "title", title,
                        "parent", parent,
                        "action", action,
                        "accept-label", accept_label,
                        "cancel-label", cancel_label,
                        NULL);
}

// Macro to define the custom widget type
G_DEFINE_TYPE(SirilFileChooserButton, siril_file_chooser_button, GTK_TYPE_BUTTON);
