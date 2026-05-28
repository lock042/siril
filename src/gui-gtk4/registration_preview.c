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
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/sequence_list.h"
#include "registration/registration.h"
#include "opencv/opencv.h"

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

/* GtkDrawingAreaDrawFunc — was the GTK3 "draw" signal handler bound
 * from siril.ui on drawingarea_reg_manual_preview1/2.  GTK4 dropped
 * that signal; the new way is gtk_drawing_area_set_draw_func() with
 * this 5-arg signature.  Wired from callbacks.c after the
 * gui.preview_area[] pointers are populated. */
void redraw_preview(GtkDrawingArea *area, cairo_t *cr,
                    int width, int height, gpointer data) {
	GtkWidget *widget = GTK_WIDGET(area);
	int current_preview, shiftx = 0, shifty = 0;
	guint area_width = (guint) width;
	guint area_height = (guint) height;
	gboolean display_ref_image;

	registration_preview_init_statics();
	display_ref_image = gui.refimage_regbuffer && gui.refimage_surface
			&& siril_toggle_get_active(GTK_WIDGET(regprev_check_display_ref))
							&& !gtk_widget_get_visible(regprev_label_regref);

	if (widget == gui.preview_area[0]) current_preview = 0;
	else if (widget == gui.preview_area[1]) current_preview = 1;
	else {
		fprintf(stderr, "Uninitialized com.preview_area or unknown drawing area!\n");
		return;
	}
	// update previewW and previewH in case the drawing area was resized
	com.seq.previewW[current_preview] = area_width;
	com.seq.previewH[current_preview] = area_height;

	/* fill preview bacground */
	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	/* display current image with shifts */
	if (!gui.preview_surface[current_preview]) {
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
	cairo_translate(cr, area_width / 2.0 - com.seq.previewX[current_preview],
			area_height/2.0-com.seq.previewY[current_preview]);
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
	if (shiftx || shifty)
		cairo_translate(cr, shiftx, -shifty);
	cairo_set_source_surface(cr, gui.preview_surface[current_preview], 0, 0);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);

	/* display reference image */
	if (display_ref_image) {
		/* already translated above but undo the shift translate */
		if (shiftx || shifty)
			cairo_translate(cr, -shiftx, shifty);
		cairo_set_source_surface(cr, gui.refimage_surface, 0, 0);
		cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
		cairo_paint_with_alpha(cr, 0.5);
	}
}

/* vport can be -1 if the correct viewport should be tested */
void test_and_allocate_reference_image(int vport) {
	registration_preview_init_statics();
	if (vport == -1)
		vport = gtk_drop_down_get_selected(regprev_layer_combo);

	if (sequence_is_loaded() && com.seq.current == com.seq.reference_image
			&& gtk_drop_down_get_selected(regprev_layer_combo) == vport && vport < gfit->naxes[2]) {
		/* this is the registration layer and the reference frame,
		 * save the buffer for alignment preview */
		struct image_view *view = &gui.view[vport];
		/* Use actual surface dimensions — may be smaller than gfit when the
		 * image exceeds Cairo's 32767-px limit (gui.surface_scale < 1.0). */
		const int surf_h = view->buf_height;
		const int surf_w = view->buf_stride / 4; /* stride = sw*4 for RGB24 */
		if (!gui.refimage_regbuffer || !gui.refimage_surface) {
			guchar *oldbuf = gui.refimage_regbuffer;
			gui.refimage_regbuffer = realloc(gui.refimage_regbuffer,
					view->buf_stride * surf_h * sizeof(guchar));
			if (gui.refimage_regbuffer == NULL) {
				PRINT_ALLOC_ERR;
				if (oldbuf)
					free(oldbuf);
				return;
			}

			if (gui.refimage_surface)
				cairo_surface_destroy(gui.refimage_surface);
			gui.refimage_surface = cairo_image_surface_create_for_data(
					gui.refimage_regbuffer, CAIRO_FORMAT_RGB24, surf_w,
					surf_h, view->buf_stride);
			if (cairo_surface_status(gui.refimage_surface)
					!= CAIRO_STATUS_SUCCESS) {
				fprintf(stderr,
						"Error creating the Cairo image surface for the reference image.\n");
				cairo_surface_destroy(gui.refimage_surface);
				gui.refimage_surface = NULL;
			} else {
				fprintf(stdout,
						"Saved the reference frame buffer for alignment preview.\n");
				enable_view_reference_checkbox(TRUE);
			}
		}
		memcpy(gui.refimage_regbuffer, view->buf,
				view->buf_stride * surf_h * sizeof(guchar));
		cairo_surface_flush(gui.refimage_surface);
		cairo_surface_mark_dirty(gui.refimage_surface);
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
	/* free alignment preview data */
	for (int i = 0; i < PREVIEW_NB; i++) {
		if (gui.preview_surface[i]) {
			cairo_surface_destroy(gui.preview_surface[i]);
			gui.preview_surface[i] = NULL;
		}
	}

	registration_preview_init_statics();
	siril_toggle_set_active(GTK_WIDGET(regprev_toggle1), FALSE);
	siril_toggle_set_active(GTK_WIDGET(regprev_toggle2), FALSE);
}

void set_preview_area(int preview_area, int centerX, int centerY) {
	/* equivalent to remap() for full image visualization, called in the
	 * mouse click callback and image change.
	 * load the preview area from reference image, from current image,
	 * sets the cairo_surface_t */
	guint area_width = gtk_widget_get_width (gui.preview_area[preview_area]);
	guint area_height = gtk_widget_get_height (gui.preview_area[preview_area]);

	com.seq.previewX[preview_area] = centerX;
	com.seq.previewY[preview_area] = centerY;
	com.seq.previewW[preview_area] = area_width;
	com.seq.previewH[preview_area] = area_height;

	struct image_view *view = &gui.view[gui.cvport];
	/* Surface dimensions may be smaller than gfit when gui.surface_scale < 1. */
	const int surf_h = view->buf_height;
	const int surf_w = view->buf_stride / 4;
	if (!gui.preview_surface[preview_area] ||
			cairo_image_surface_get_width(gui.preview_surface[preview_area]) != surf_w ||
			cairo_image_surface_get_height(gui.preview_surface[preview_area]) != surf_h) {
		if (gui.preview_surface[preview_area])
			cairo_surface_destroy(gui.preview_surface[preview_area]);
		gui.preview_surface[preview_area] =
			cairo_image_surface_create_for_data(view->buf,
					CAIRO_FORMAT_RGB24,
					surf_w, surf_h,
					view->buf_stride);
		if (cairo_surface_status(gui.preview_surface[preview_area]) != CAIRO_STATUS_SUCCESS) {
			fprintf(stderr, "Error creating the Cairo image surface for preview %d\n", preview_area);
			cairo_surface_destroy(gui.preview_surface[preview_area]);
			gui.preview_surface[preview_area] = NULL;
		}
	}
	// queue for redraw
	gtk_widget_queue_draw(gui.preview_area[preview_area]);
	// flush to ensure all writing to the image was done and redraw the surface
	//cairo_surface_flush (com.preview_surface[preview_area]);
	//cairo_surface_mark_dirty (com.preview_surface[preview_area]);
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
			cairo_surface_destroy(gui.preview_surface[preview_area]);
			gui.preview_surface[preview_area] = NULL;
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
