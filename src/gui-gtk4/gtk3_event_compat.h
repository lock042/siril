/*
 * GTK3 → GTK4 transition shim.
 *
 * In GTK3 each input event had its own type (GdkEventButton, GdkEventScroll,
 * GdkEventMotion, GdkEventKey, GdkEventCrossing).  GTK4 folds them all into
 * a single opaque GdkEvent.  The full migration is tracked in Phases 8/8B
 * of gtk4-migration-plan.txt; this header lets headers and signatures that
 * still reference the legacy type names compile by aliasing them onto
 * GdkEvent.  Field access on these types (event->x, event->button, …) is
 * not provided here — it must be rewritten via GtkEventController APIs.
 */
#ifndef _GTK3_EVENT_COMPAT_H_
#define _GTK3_EVENT_COMPAT_H_

#include <gtk/gtk.h>

typedef GdkEvent GdkEventButton;
typedef GdkEvent GdkEventScroll;
typedef GdkEvent GdkEventMotion;
typedef GdkEvent GdkEventKey;
typedef GdkEvent GdkEventCrossing;

/* Phase 13: GtkClipboard typedef removed; all gui-gtk4/ sites now use
 * GdkClipboard directly via gdk_display_get_clipboard()/
 * gtk_widget_get_clipboard(). */

/* Phase 14B: GTK4 changed gtk_file_chooser_set_current_folder to take
 * (chooser, GFile*, GError**) instead of (chooser, const char *path).
 * The string-based wrapper siril_file_chooser_set_current_folder_path
 * lives in utils.c next to its sibling helpers; declared in utils.h. */

/* GTK4 removed GDK_DOUBLE_BUTTON_PRESS / GDK_TRIPLE_BUTTON_PRESS from the
 * GdkEventType enum.  The replacement is the n_press counter on the
 * GtkGestureClick controller.  Until Phase 8B reroutes those reads to the
 * controller API, several files still use the old enum names as in-memory
 * sentinels (for storing/comparing "double-click" vs "single-click" mouse
 * configuration).  Define them as out-of-range values that cannot collide
 * with real GdkEventType members. */
#ifndef GDK_DOUBLE_BUTTON_PRESS
#define GDK_DOUBLE_BUTTON_PRESS  ((GdkEventType)0x1000)
#endif
#ifndef GDK_TRIPLE_BUTTON_PRESS
#define GDK_TRIPLE_BUTTON_PRESS  ((GdkEventType)0x1001)
#endif

#endif /* _GTK3_EVENT_COMPAT_H_ */
