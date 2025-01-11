/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/image_interactions.h"
#include "gui/registration_preview.h"
#include "gui/sequence_list.h"
#include "registration/registration.h"

gboolean redraw_preview(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int current_preview, shiftx = 0, shifty = 0;
	static GtkToggleButton *check_display_ref = NULL;
	static GtkWidget *labelRegRef = NULL;
	guint area_width = gtk_widget_get_allocated_width (widget);
	guint area_height = gtk_widget_get_allocated_height (widget);
	gboolean display_ref_image;

	if (check_display_ref == NULL) {
		check_display_ref = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_displayref"));
		labelRegRef = lookup_widget("labelRegRef");
	}
	display_ref_image = gui.refimage_regbuffer && gui.refimage_surface
			&& gtk_toggle_button_get_active(check_display_ref)
							&& !gtk_widget_get_visible(labelRegRef);

	if (widget == gui.preview_area[0]) current_preview = 0;
	else if (widget == gui.preview_area[1]) current_preview = 1;
	else {
		fprintf(stderr, "Uninitialized com.preview_area or unknown drawing area!\n");
		return TRUE;
	}
	// update previewW and previewH in case the drawing area was resized
	com.seq.previewW[current_preview] = area_width;
	com.seq.previewH[current_preview] = area_height;

	/* fill preview bacground */
	cairo_rectangle(cr, 0, 0, area_width, area_height);
	cairo_fill(cr);

	/* display current image with shifts */
	if (!gui.preview_surface[current_preview]) {
		GtkStyleContext *context = gtk_widget_get_style_context(widget);
		GtkStateFlags state = gtk_widget_get_state_flags(widget);
		PangoLayout *layout;
		gchar *msg;
		GtkAllocation allocation;
		gdouble scale;
		GdkRGBA color;
		gint w, h;

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
		gtk_widget_get_allocation(widget, &allocation);

		scale = MIN(((gdouble) allocation.width / 3.0) / (gdouble ) w,
				((gdouble) allocation.height / 3.0) / (gdouble ) h);

		gtk_style_context_get_color(context, state, &color);
		gdk_cairo_set_source_rgba(cr, &color);

		cairo_move_to(cr, (allocation.width - (w * scale)) / 2,
				(allocation.height - (h * scale)) / 2);

		cairo_scale(cr, scale, scale);

		pango_cairo_show_layout(cr, layout);

		g_object_unref(layout);
		return TRUE;
	}
	cairo_translate(cr, area_width / 2.0 - com.seq.previewX[current_preview],
			area_height/2.0-com.seq.previewY[current_preview]);
	if (com.seq.regparam) {
		double dx, dy;
		translation_from_H(com.seq.regparam[com.seq.current].H, &dx, &dy);
		shiftx = round_to_int(dx);
		shifty = round_to_int(dy);
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

	return FALSE;
}

/* vport can be -1 if the correct viewport should be tested */
void test_and_allocate_reference_image(int vport) {
	static GtkComboBox *cbbt_layers = NULL;
	if (cbbt_layers == NULL) {
		cbbt_layers = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
	}
	if (vport == -1)
		vport = gtk_combo_box_get_active(cbbt_layers);

	if (sequence_is_loaded() && com.seq.current == com.seq.reference_image
			&& gtk_combo_box_get_active(cbbt_layers) == vport && vport < gfit.naxes[2]) {
		/* this is the registration layer and the reference frame,
		 * save the buffer for alignment preview */
		struct image_view *view = &gui.view[vport];
		if (!gui.refimage_regbuffer || !gui.refimage_surface) {
			guchar *oldbuf = gui.refimage_regbuffer;
			gui.refimage_regbuffer = realloc(gui.refimage_regbuffer,
					view->full_surface_stride * gfit.ry * sizeof(guchar));
			if (gui.refimage_regbuffer == NULL) {
				PRINT_ALLOC_ERR;
				if (oldbuf)
					free(oldbuf);
				return;
			}

			if (gui.refimage_surface)
				cairo_surface_destroy(gui.refimage_surface);
			gui.refimage_surface = cairo_image_surface_create_for_data(
					gui.refimage_regbuffer, CAIRO_FORMAT_RGB24, gfit.rx,
					gfit.ry, view->full_surface_stride);
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
				view->full_surface_stride * gfit.ry * sizeof(guchar));
		cairo_surface_flush(gui.refimage_surface);
		cairo_surface_mark_dirty(gui.refimage_surface);
	}
}


void redraw_previews() {
	int i;
	if (com.script) return;
	for (i = 0; i < PREVIEW_NB; i++)
		gtk_widget_queue_draw(gui.preview_area[i]);
}

void clear_previews() {
	/* free alignment preview data */
	for (int i = 0; i < PREVIEW_NB; i++) {
		if (gui.preview_surface[i]) {
			cairo_surface_destroy(gui.preview_surface[i]);
			gui.preview_surface[i] = NULL;
		}
	}

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_reg_manual1")), FALSE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_reg_manual2")), FALSE);
}

void set_preview_area(int preview_area, int centerX, int centerY) {
	/* equivalent to remap() for full image visualization, called in the
	 * mouse click callback and image change.
	 * load the preview area from reference image, from current image,
	 * sets the cairo_surface_t */
	guint area_width = gtk_widget_get_allocated_width (gui.preview_area[preview_area]);
	guint area_height = gtk_widget_get_allocated_height (gui.preview_area[preview_area]);

	com.seq.previewX[preview_area] = centerX;
	com.seq.previewY[preview_area] = centerY;
	com.seq.previewW[preview_area] = area_width;
	com.seq.previewH[preview_area] = area_height;

	struct image_view *view = &gui.view[gui.cvport];
	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx) !=
			view->full_surface_stride ||
			gfit.ry != view->full_surface_height ||
			!gui.preview_surface[preview_area]) {
		if (gui.preview_surface[preview_area])
			cairo_surface_destroy(gui.preview_surface[preview_area]);
		gui.preview_surface[preview_area] =
			cairo_image_surface_create_for_data(view->buf,
					CAIRO_FORMAT_RGB24,
					gfit.rx, gfit.ry,
					view->full_surface_stride);
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

		static GtkToggleButton *preview1 = NULL;

		if (preview1 == NULL)
			preview1 = GTK_TOGGLE_BUTTON(lookup_widget("toggle_reg_manual1"));
		if (gtk_toggle_button_get_active(toggle)) {
			if (toggle == preview1)
				mouse_status = MOUSE_ACTION_SELECT_PREVIEW1;
			else
				mouse_status = MOUSE_ACTION_SELECT_PREVIEW2;
		} else {
			/* deactivate preview */
			int preview_area;
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			if (toggle == preview1)
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

void on_checkbutton_displayref_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	redraw_previews();
}

/* display registration data (shift{x|y} for now) in the manual adjustments */
void adjust_reginfo() {
	GtkSpinButton *spin_shiftx, *spin_shifty;
	gboolean set_sensitive;

	spin_shiftx = GTK_SPIN_BUTTON(lookup_widget("spinbut_shiftx"));
	spin_shifty = GTK_SPIN_BUTTON(lookup_widget("spinbut_shifty"));

	g_signal_handlers_block_by_func(spin_shiftx, on_spinbut_shift_value_change, NULL);
	g_signal_handlers_block_by_func(spin_shifty, on_spinbut_shift_value_change, NULL);
	if (com.seq.regparam == NULL) {
		gtk_spin_button_set_value(spin_shiftx, 0.);
		gtk_spin_button_set_value(spin_shifty, 0.);
	} else {
		double dx, dy;
		translation_from_H(com.seq.regparam[com.seq.current].H, &dx, &dy);
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
	static GtkComboBox *cbbt_layers = NULL;
	int current_layer, new_value;
	if (spin_shiftx == NULL) {
		spin_shiftx = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbut_shiftx"));
		cbbt_layers = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "comboboxreglayer"));
	}
	if (com.seq.regparam == NULL) {
		/* allocated when the number of layers is loaded from the sequence,
		 * = shouldn't happen */
		fprintf(stderr, "regparam not allocated, sequence not loaded, never displayed or malformed\n");
		return;
	}

	current_layer = gtk_combo_box_get_active(cbbt_layers);
	activate_tab(current_layer);

	if (com.seq.regparam == NULL) {
		printf("Allocating registration data\n");
		com.seq.regparam = calloc(com.seq.number, sizeof(regdata));
		if (com.seq.regparam == NULL) {
			PRINT_ALLOC_ERR;
			return;
		}
	}

	new_value = gtk_spin_button_get_value_as_int(spinbutton);
	double shiftx, shifty;
	translation_from_H(com.seq.regparam[com.seq.current].H, &shiftx, &shifty);
	if (spinbutton == spin_shiftx)
		shiftx = new_value;
	else shifty = new_value;
	com.seq.regparam[com.seq.current].H = H_from_translation(shiftx, shifty);
	writeseqfile(&com.seq);
	update_seqlist();
	fill_sequence_list(&com.seq, FALSE);	// update list with new regparam
	redraw_previews();
}

/* enables or disables the "display reference" checkbox in registration preview */
void enable_view_reference_checkbox(gboolean status) {
	static GtkToggleButton *check_display_ref = NULL;
	static GtkWidget *widget = NULL, *labelRegRef = NULL;
	if (check_display_ref == NULL) {
		check_display_ref = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_displayref"));
		widget = GTK_WIDGET(check_display_ref);
		labelRegRef = lookup_widget("labelRegRef");
	}
	if (status && gtk_widget_get_sensitive(widget))
		return;	// may be already enabled but deactivated by user, don't force it again
	gtk_widget_set_sensitive(widget, status);
	gtk_widget_set_visible(labelRegRef, !status);
	gtk_toggle_button_set_active(check_display_ref, status);
}

