#ifndef _IMAGE_DISPLAY_H_
#define _IMAGE_DISPLAY_H_

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/gui_iface.h"   /* SirilRedrawType, REDRAW_OVERLAY/REDRAW_IMAGE/REDRAW_ALL */

/* Backward-compatibility alias: existing callers use remap_type; new code
 * should use SirilRedrawType from core/gui_iface.h directly. */
typedef SirilRedrawType remap_type;

void check_gfit_profile_identical_to_monitor();

void allocate_hd_remap_indices();
void hd_remap_indices_cleanup();

void initialize_image_display();

void copy_roi_into_gfit();

void remap_all();
/* Invalidate gfit stats/histogram and, if on the GTK main thread, also remap
 * the Cairo buffers.  Safe to call from any thread; non-main-thread callers
 * have the remap deferred to end_gfit_operation(). */
void redraw(remap_type doremap);	// redraw the image
void queue_redraw(remap_type doremap); // call redraw from another thread
void queue_redraw_and_wait_for_it(remap_type doremap); // call redraw from another thread and wait for it
gboolean redraw_mask_idle(gpointer p);
void queue_redraw_mask(); // queue a redraw of the mask only
void block_drawarea_handlers(void);   // block viewport draw signal handlers (GTK main thread only)
void unblock_drawarea_handlers(void); // unblock viewport draw signal handlers (GTK main thread only)
void install_drawarea_draw_funcs(void); // wire viewports to redraw_drawingarea (GTK4)

/* Custom GtkWidget that renders an image viewport via GtkSnapshot + GdkTexture
 * (Phase 1 of the image-display GPU migration).  Replaces GtkDrawingArea in
 * siril.ui so we can append a GdkTexture for the image and a Cairo subnode
 * for overlays in a single snapshot pass.  The type must be registered with
 * GType before gtk_builder reads the XML, see siril_image_view_register(). */
#define SIRIL_TYPE_IMAGE_VIEW (siril_image_view_get_type())
G_DECLARE_FINAL_TYPE(SirilImageView, siril_image_view, SIRIL, IMAGE_VIEW, GtkWidget)

/* Ensure the SirilImageView GType is registered.  Called once from
 * application startup before any GtkBuilder XML is parsed. */
void siril_image_view_register(void);

double get_zoom_val();	// for image_interactions

point get_center_of_vport();
void add_image_and_label_to_cairo(cairo_t *cr, int vport);

gboolean get_context_rotation_matrix(double rotation, cairo_matrix_t *transform, gboolean invert); //computes rotation matrix about center of com.selection

#endif

