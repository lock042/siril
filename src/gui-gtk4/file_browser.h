#ifndef SRC_GUI_GTK4_FILE_BROWSER_H_
#define SRC_GUI_GTK4_FILE_BROWSER_H_

#include <gtk/gtk.h>

/* SirilFileBrowser
 *
 * A custom modal file picker dialog built on the GTK4 list/grid APIs
 * (GtkDirectoryList → SortListModel → FilterListModel → SingleSelection
 * → GtkListView).  Designed for use cases where GtkFileDialog's lack of
 * a preview slot is a blocker — e.g. picking a FITS file to remix
 * against another, or the main Open File dialog.
 *
 * Synchronous API: siril_file_browser_run() spins a local GMainLoop and
 * returns GTK_RESPONSE_ACCEPT or GTK_RESPONSE_CANCEL.
 *
 * Preview rendering is provided either by the built-in default handler
 * (GdkTexture-loads any format GTK knows, including raster formats; FITS
 * via Siril's loader) or by a user-supplied callback.
 */

typedef struct _SirilFileBrowser SirilFileBrowser;

/* Preview callback.  Invoked on the main thread whenever the selection
 * changes (path is NULL when nothing is selected, or selected is a
 * directory).  The callback should set the preview image on `picture`
 * (e.g. via gtk_picture_set_paintable / _set_filename) and may also
 * update `metadata_label`.  Use siril_file_browser_default_preview()
 * for the built-in handler. */
typedef void (*SirilFileBrowserPreview)(const gchar *path,
                                        GtkPicture  *picture,
                                        GtkLabel    *metadata_label,
                                        gpointer     user_data);

SirilFileBrowser *siril_file_browser_new(GtkWindow *parent, const gchar *title);

void siril_file_browser_set_initial_folder (SirilFileBrowser *fb, const gchar *path);
void siril_file_browser_set_initial_file   (SirilFileBrowser *fb, const gchar *path);
void siril_file_browser_add_filter_pattern (SirilFileBrowser *fb,
                                            const gchar *title,
                                            const gchar *pattern,
                                            gboolean     set_default);

/* Enable multi-file selection.  Must be called before _run().  Ctrl/Shift
 * clicks then extend the selection; double-click still activates a single
 * item (directory navigation or single-file accept).  Use
 * siril_file_browser_get_paths() to retrieve every selected path. */
void siril_file_browser_set_select_multiple(SirilFileBrowser *fb, gboolean multi);

/* Show a "Debayer" check button in the action row.  When toggled, the
 * browser updates com.pref.debayer.open_debayer and mirrors the change
 * onto the shared `demosaicingButton` (Convert tab), and vice-versa —
 * the toggle reflects the pref's current state on first show.  Use this
 * only on Open / Convert image pickers where the user might want to
 * flip the debayer-on-open behaviour without leaving the dialog. */
void siril_file_browser_set_show_debayer_toggle(SirilFileBrowser *fb, gboolean show);

/* Provide a custom preview callback.  When unset (default), the browser
 * uses siril_file_browser_default_preview() which handles common image
 * formats via GdkTexture and FITS via Siril's reader. */
void siril_file_browser_set_preview_callback(SirilFileBrowser       *fb,
                                             SirilFileBrowserPreview cb,
                                             gpointer                user_data);

/* Built-in preview handler — exported so callers can compose it. */
void siril_file_browser_default_preview(const gchar *path,
                                        GtkPicture  *picture,
                                        GtkLabel    *metadata_label,
                                        gpointer     user_data);

/* Synchronous run.  Returns GTK_RESPONSE_ACCEPT or GTK_RESPONSE_CANCEL. */
gint     siril_file_browser_run      (SirilFileBrowser *fb);
gchar   *siril_file_browser_get_path (SirilFileBrowser *fb);
/* Multi-select result.  Caller g_free's each entry and frees the list. */
GSList  *siril_file_browser_get_paths(SirilFileBrowser *fb);

/* No siril_file_browser_destroy: the browser is a process-wide singleton,
 * hidden on close and reset on next _new().  See file_browser.c for the
 * background — this design sidesteps the macOS AppKit teardown bug. */

/* Inline image-picker button helper.
 *
 * Replaces GTK3's GtkFileChooserButton (which GTK4 removed) for the
 * specific case where the .ui has a plain GtkButton that should open a
 * SirilFileBrowser (with image preview) when clicked.  On accept:
 *   - The chosen absolute path is stashed on the button via
 *     g_object_set_data("siril-path", ...) so the existing
 *     siril_file_chooser_get_filename(GtkFileChooser*) shim continues to
 *     read it back transparently.
 *   - The button label is updated to the file's basename.
 *   - `on_picked` (if non-NULL) is invoked on the main thread with the
 *     newly chosen path so the consumer can load / cache as needed.
 */
typedef void (*SirilImageButtonCallback)(GtkWidget   *button,
                                         const gchar *path,
                                         gpointer     user_data);

void siril_image_button_init(GtkWidget               *button,
                             const gchar             *title,
                             const gchar             *filter_title,
                             const gchar             *filter_pattern,
                             SirilImageButtonCallback on_picked,
                             gpointer                 user_data);

#endif /* SRC_GUI_GTK4_FILE_BROWSER_H_ */
