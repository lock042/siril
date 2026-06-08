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

#include "core/siril.h"
#include "core/proto.h"
#include "io/sequence.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/sequence_list.h"
#include "registration/registration.h"
#include "opencv/opencv.h"

/* Dimensions of the persisted reference-frame pixel buffer
 * (gui.refimage_regbuffer).  The buffer holds the whole reference frame's
 * display pixels (a plain malloc — no Cairo size limit); only small windows
 * of it are ever turned into Cairo surfaces, in make_ref_window().  Valid
 * only while gui.refimage_regbuffer != NULL. */
static int ref_w = 0, ref_h = 0, ref_stride = 0;

static GtkCheckButton *regprev_check_display_ref = NULL;
static GtkWidget *regprev_label_regref = NULL;
static GtkDropDown *regprev_layer_combo = NULL;
static GtkToggleButton *regprev_toggle1 = NULL;
static GtkToggleButton *regprev_toggle2 = NULL;
static GtkSpinButton *regprev_spin_shiftx = NULL;
static GtkSpinButton *regprev_spin_shifty = NULL;
static GtkDropDown *regprev_seq_combo = NULL;

static void registration_preview_init_statics(void) {
	if (regprev_check_display_ref) return;
	regprev_check_display_ref = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_displayref"));
	regprev_label_regref = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelRegRef"));
	regprev_layer_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboboxreglayer"));
	regprev_toggle1 = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "toggle_reg_manual1"));
	regprev_toggle2 = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "toggle_reg_manual2"));
	regprev_spin_shiftx = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbut_shiftx"));
	regprev_spin_shifty = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbut_shifty"));
	regprev_seq_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "seqlist_dialog_combo"));
}

/* Build a w×h RGB24 surface holding the image-space window (ox, oy) of
 * viewport `vport`'s live display buffer, materialising only the tiles it
 * covers.  Pixels outside the image stay black.  Caller owns the surface. */
static cairo_surface_t *make_view_window(int vport, int ox, int oy, int w, int h) {
	if (w <= 0 || h <= 0)
		return NULL;
	cairo_surface_t *surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, w, h);
	if (cairo_surface_status(surf) != CAIRO_STATUS_SUCCESS) {
		cairo_surface_destroy(surf);
		return NULL;
	}
	cairo_surface_flush(surf);
	guchar *data = cairo_image_surface_get_data(surf);
	const int stride = cairo_image_surface_get_stride(surf);
	memset(data, 0, (size_t) stride * h);
	siril_image_view_copy_region(vport, ox, oy, w, h, data, stride);
	cairo_surface_mark_dirty(surf);
	return surf;
}

/* As make_view_window, but the source is the persisted reference buffer
 * (gui.refimage_regbuffer) rather than a live viewport. */
static cairo_surface_t *make_ref_window(int ox, int oy, int w, int h) {
	if (!gui.refimage_regbuffer || ref_w <= 0 || ref_h <= 0 || w <= 0 || h <= 0)
		return NULL;
	cairo_surface_t *surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, w, h);
	if (cairo_surface_status(surf) != CAIRO_STATUS_SUCCESS) {
		cairo_surface_destroy(surf);
		return NULL;
	}
	cairo_surface_flush(surf);
	guchar *data = cairo_image_surface_get_data(surf);
	const int stride = cairo_image_surface_get_stride(surf);
	memset(data, 0, (size_t) stride * h);
	const int ix0 = MAX(0, ox);
	const int ix1 = MIN(ref_w, ox + w);
	for (int row = 0; row < h; row++) {
		const int iy = oy + row;
		if (iy < 0 || iy >= ref_h || ix1 <= ix0)
			continue;
		memcpy(data + (size_t) row * stride + (size_t)(ix0 - ox) * 4,
		       gui.refimage_regbuffer + (size_t) iy * ref_stride + (size_t) ix0 * 4,
		       (size_t)(ix1 - ix0) * 4);
	}
	cairo_surface_mark_dirty(surf);
	return surf;
}

/* GtkDrawingAreaDrawFunc — was the GTK3 "draw" signal handler bound
 * from siril.ui on drawingarea_reg_manual_preview1/2.  GTK4 dropped
 * that signal; the new way is gtk_drawing_area_set_draw_func() with
 * this 5-arg signature.  Wired from callbacks.c after the
 * gui.preview_area[] pointers are populated. */
void redraw_preview(GtkDrawingArea *area, cairo_t *cr,
                    int width, int height, gpointer data) {
	GtkWidget *widget = GTK_WIDGET(area);
	int current_preview, shiftx = 0, shifty = 0;
	const int area_width = width;
	const int area_height = height;
	gboolean display_ref_image;

	registration_preview_init_statics();

	if (widget == gui.preview_area[0]) current_preview = 0;
	else if (widget == gui.preview_area[1]) current_preview = 1;
	else {
		fprintf(stderr, "Uninitialized com.preview_area or unknown drawing area!\n");
		return;
	}
	// update previewW and previewH in case the drawing area was resized
	com.seq.previewW[current_preview] = area_width;
	com.seq.previewH[current_preview] = area_height;

	/* fill preview background */
	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	/* No preview centre selected yet: show the placeholder label. */
	if (com.seq.previewX[current_preview] < 0) {
		(void) gtk_widget_get_state_flags(widget);
		PangoLayout *layout;
		gchar *msg;
		gdouble scale;
		GdkRGBA color;
		gint w, h;
		const int alloc_w = gtk_widget_get_width(widget);
		const int alloc_h = gtk_widget_get_height(widget);

		layout = gtk_widget_create_pango_layout(widget, NULL);

		if (sequence_is_loaded()) {
			msg = g_strdup_printf(_("Preview %d"), current_preview + 1);
		} else {
			msg = g_strdup(_("Load\nsequences"));
		}
		pango_layout_set_markup(layout, msg, -1);
		g_free(msg);
		pango_layout_set_alignment(layout, PANGO_ALIGN_CENTER);

		pango_layout_get_pixel_size(layout, &w, &h);
		scale = MIN(((gdouble) alloc_w / 3.0) / (gdouble) w,
				((gdouble) alloc_h / 3.0) / (gdouble) h);

		gtk_widget_get_color(widget, &color);
		gdk_cairo_set_source_rgba(cr, &color);

		cairo_move_to(cr, (alloc_w - (w * scale)) / 2,
				(alloc_h - (h * scale)) / 2);

		cairo_scale(cr, scale, scale);

		pango_cairo_show_layout(cr, layout);

		g_object_unref(layout);
		return;
	}

	display_ref_image = gui.refimage_regbuffer
			&& siril_toggle_get_active(GTK_WIDGET(regprev_check_display_ref))
			&& !gtk_widget_get_visible(regprev_label_regref);

	int cvport = select_vport(gui.cvport);
	if (com.seq.regparam && com.seq.regparam[cvport]) {
		double dx, dy;
		translation_from_H(com.seq.regparam[cvport][com.seq.current].H, &dx, &dy);
		shiftx = round_to_int(dx);
		shifty = round_to_int(dy);
		if (shiftx == INT_MIN) { // mainly to avoid static checker warning
			siril_log_debug("Error: image #%d has a wrong shift x value\n", com.seq.current + 1);
			shiftx += 1;
		}
		if (shifty == INT_MIN) { // mainly to avoid static checker warning
			siril_log_debug("Error: image #%d has a wrong shift y value\n", com.seq.current + 1);
			shifty += 1;
		}
	}

	const int pX = com.seq.previewX[current_preview];
	const int pY = com.seq.previewY[current_preview];

	/* Current image: render exactly the visible window from the live
	 * display buffer.  The window origin folds in the registration shift so
	 * the geometry matches the former full-image-surface + translate path. */
	cairo_surface_t *cur_surf = make_view_window(gui.cvport,
			pX - shiftx - area_width / 2,
			pY + shifty - area_height / 2,
			area_width, area_height);
	if (cur_surf) {
		cairo_set_source_surface(cr, cur_surf, 0, 0);
		cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
		cairo_paint(cr);
		cairo_surface_destroy(cur_surf);
	}

	/* Reference overlay (unshifted), from the persisted reference buffer. */
	if (display_ref_image) {
		cairo_surface_t *ref_surf = make_ref_window(
				pX - area_width / 2,
				pY - area_height / 2,
				area_width, area_height);
		if (ref_surf) {
			cairo_set_source_surface(cr, ref_surf, 0, 0);
			cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
			cairo_paint_with_alpha(cr, 0.5);
			cairo_surface_destroy(ref_surf);
		}
	}
}

/* vport can be -1 if the correct viewport should be tested */
void test_and_allocate_reference_image(int vport) {
	registration_preview_init_statics();
	if (vport == -1)
		vport = gtk_drop_down_get_selected(regprev_layer_combo);

	if (sequence_is_loaded() && com.seq.current == com.seq.reference_image
			&& gtk_drop_down_get_selected(regprev_layer_combo) == vport && vport < gfit->naxes[2]) {
		/* This is the registration layer and the reference frame: persist
		 * its display pixels for the alignment preview.  We keep a full
		 * pixel copy (a plain malloc, no Cairo size limit) and assemble it
		 * with siril_image_view_copy_region so it works in lazy mode too
		 * (where there is no contiguous view->buf); only small windows of
		 * it are later turned into Cairo surfaces in make_ref_window(). */
		const int rx = (int) gfit->rx;
		const int ry = (int) gfit->ry;
		const int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, rx);
		guchar *newbuf = realloc(gui.refimage_regbuffer, (size_t) stride * ry);
		if (newbuf == NULL) {
			PRINT_ALLOC_ERR;
			return;
		}
		gui.refimage_regbuffer = newbuf;
		ref_w = rx;
		ref_h = ry;
		ref_stride = stride;
		memset(gui.refimage_regbuffer, 0, (size_t) stride * ry);
		if (!siril_image_view_copy_region(vport, 0, 0, rx, ry,
				gui.refimage_regbuffer, stride)) {
			fprintf(stderr, "Error capturing the reference image for alignment preview.\n");
			free(gui.refimage_regbuffer);
			gui.refimage_regbuffer = NULL;
			ref_w = ref_h = ref_stride = 0;
			return;
		}
		fprintf(stdout, "Saved the reference frame buffer for alignment preview.\n");
		enable_view_reference_checkbox(TRUE);
	}
}

gboolean redraw_previews(gpointer user_data) {
	int i;
	if (com.script) return FALSE;
	for (i = 0; i < PREVIEW_NB; i++)
		gtk_widget_queue_draw(gui.preview_area[i]);
	return FALSE;
}

void clear_previews() {
	registration_preview_init_statics();
	for (int i = 0; i < PREVIEW_NB; i++) {
		com.seq.previewX[i] = -1;
		com.seq.previewY[i] = -1;
		com.seq.previewW[i] = 0;
		com.seq.previewH[i] = 0;
		if (gui.preview_area[i])
			gtk_widget_queue_draw(gui.preview_area[i]);
	}
	siril_toggle_set_active(GTK_WIDGET(regprev_toggle1), FALSE);
	siril_toggle_set_active(GTK_WIDGET(regprev_toggle2), FALSE);
}

void set_preview_area(int preview_area, int centerX, int centerY) {
	/* Records the preview centre; the visible window is rendered on demand
	 * in redraw_preview() from the live display buffer (and, for the
	 * overlay, the persisted reference buffer).  Called from the mouse
	 * click callback and on image change. */
	com.seq.previewX[preview_area] = centerX;
	com.seq.previewY[preview_area] = centerY;
	com.seq.previewW[preview_area] = gtk_widget_get_width(gui.preview_area[preview_area]);
	com.seq.previewH[preview_area] = gtk_widget_get_height(gui.preview_area[preview_area]);

	gtk_widget_queue_draw(gui.preview_area[preview_area]);
}

void on_toggle_preview_toggled(GtkToggleButton *toggle, gpointer user_data) {

	if (sequence_is_loaded()) {

		registration_preview_init_statics();
		if (siril_toggle_get_active(GTK_WIDGET(toggle))) {
			if (toggle == regprev_toggle1)
				mouse_status = MOUSE_ACTION_SELECT_PREVIEW1;
			else
				mouse_status = MOUSE_ACTION_SELECT_PREVIEW2;
		} else {
			/* deactivate preview */
			int preview_area;
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			if (toggle == regprev_toggle1)
				preview_area = 0;
			else
				preview_area = 1;
			com.seq.previewX[preview_area] = -1;
			com.seq.previewY[preview_area] = -1;
			com.seq.previewH[preview_area] = 0;
			com.seq.previewW[preview_area] = 0;
			// queue for redraw
			gtk_widget_queue_draw(gui.preview_area[preview_area]);
		}
	}
}

void on_checkbutton_displayref_toggled(GtkCheckButton *togglebutton, gpointer user_data) {
	gui_function(redraw_previews, NULL);
}

/* display registration data (shift{x|y} for now) in the manual adjustments */
void adjust_reginfo() {
	gboolean set_sensitive;
	gint cvport;

	registration_preview_init_statics();
	GtkSpinButton *spin_shiftx = regprev_spin_shiftx;
	GtkSpinButton *spin_shifty = regprev_spin_shifty;
	GtkDropDown *seqcombo = regprev_seq_combo;

	cvport = gtk_drop_down_get_selected(GTK_DROP_DOWN(seqcombo));
	if (cvport < 0) return;

	g_signal_handlers_block_by_func(spin_shiftx, on_spinbut_shift_value_change, NULL);
	g_signal_handlers_block_by_func(spin_shifty, on_spinbut_shift_value_change, NULL);
	if (com.seq.regparam == NULL || com.seq.regparam[cvport] == NULL) {
		gtk_spin_button_set_value(spin_shiftx, 0.);
		gtk_spin_button_set_value(spin_shifty, 0.);
	} else {
		double dx, dy;
		translation_from_H(com.seq.regparam[cvport][com.seq.current].H, &dx, &dy);
		gtk_spin_button_set_value(spin_shiftx, round_to_int(dx));
		gtk_spin_button_set_value(spin_shifty, round_to_int(dy));
	}
	g_signal_handlers_unblock_by_func(spin_shiftx, on_spinbut_shift_value_change, NULL);
	g_signal_handlers_unblock_by_func(spin_shifty, on_spinbut_shift_value_change, NULL);
	set_sensitive = (com.seq.current != com.seq.reference_image);
	gtk_widget_set_sensitive(GTK_WIDGET(spin_shiftx), set_sensitive);
	gtk_widget_set_sensitive(GTK_WIDGET(spin_shifty), set_sensitive);
}

void on_spinbut_shift_value_change(GtkSpinButton *spinbutton, gpointer user_data) {
	static GtkSpinButton *spin_shiftx = NULL;
	static GtkDropDown *cbbt_layers = NULL;
	int current_layer, new_value;
	if (spin_shiftx == NULL) {
		spin_shiftx = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbut_shiftx"));
		cbbt_layers = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboboxreglayer"));
	}
	if (com.seq.regparam == NULL) {
		/* allocated when the number of layers is loaded from the sequence,
		 * = shouldn't happen */
		fprintf(stderr, "regparam not allocated, sequence not loaded, never displayed or malformed\n");
		return;
	}

	current_layer = gtk_drop_down_get_selected(cbbt_layers);
	activate_tab(current_layer);

	if (com.seq.regparam[current_layer] == NULL) {
		printf("Allocating registration data for this layer\n");
		com.seq.regparam[current_layer] = calloc(com.seq.number, sizeof(regdata));
		if (com.seq.regparam[current_layer] == NULL) {
			PRINT_ALLOC_ERR;
			return;
		}
		cvGetEye(&com.seq.regparam[current_layer][com.seq.reference_image].H);
	}

	new_value = gtk_spin_button_get_value_as_int(spinbutton);
	double shiftx, shifty;
	translation_from_H(com.seq.regparam[current_layer][com.seq.current].H, &shiftx, &shifty);
	if (spinbutton == spin_shiftx)
		shiftx = new_value;
	else shifty = new_value;
	com.seq.regparam[current_layer][com.seq.current].H = H_from_translation(shiftx, shifty);
	writeseqfile(&com.seq);
	update_seqlist(current_layer);
	fill_sequence_list(&com.seq, current_layer, FALSE);	// update list with new regparam
	gui_function(redraw_previews, NULL);
}

/* enables or disables the "display reference" checkbox in registration preview */
void enable_view_reference_checkbox(gboolean status) {
	registration_preview_init_statics();
	GtkWidget *widget = GTK_WIDGET(regprev_check_display_ref);
	if (status && gtk_widget_get_sensitive(widget))
		return;	// may be already enabled but deactivated by user, don't force it again
	gtk_widget_set_sensitive(widget, status);
	gtk_widget_set_visible(regprev_label_regref, !status);
	siril_toggle_set_active(GTK_WIDGET(regprev_check_display_ref), status);
}
