#ifndef _IMAGE_INTERACTIONS_H_
#define _IMAGE_INTERACTIONS_H_

#include "core/siril.h"

typedef void (*selection_update_callback)();
typedef void (*star_selection_callback)(pointi);
gboolean update_zoom(gdouble x, gdouble y, double scale);
void update_zoom_fit_button();

void register_selection_update_callback(selection_update_callback f);
void unregister_selection_update_callback(const selection_update_callback f);

void new_selection_zone();
void delete_selected_area();
void reset_display_offset();
void reset_zoom_default();
void update_zoom_label();

void enforce_ratio_and_clamp();

gboolean display_quick_photo();

/* mouse behaviour */
typedef enum {
	MOUSE_ACTION_NONE,
	MOUSE_ACTION_SELECT_REG_AREA,
	MOUSE_ACTION_SELECT_PREVIEW1,
	MOUSE_ACTION_SELECT_PREVIEW2,
	MOUSE_ACTION_DRAW_SAMPLES,
	MOUSE_ACTION_PHOTOMETRY,
	MOUSE_ACTION_GET_COMP_CENTER_COORDINATE,
	MOUSE_ACTION_CUT_SELECT,
	MOUSE_ACTION_CUT_WN1,
	MOUSE_ACTION_CUT_WN2,
} mouse_status_enum;

typedef enum {
	CUT_NOT_CUTTING,
	CUT_UNCONSTRAINED,
	CUT_VERT_OR_HORIZ
} cut_method;

void init_mouse();

extern mouse_status_enum mouse_status;	// defined in registration_preview.c

extern gpointer cut_profile(); // graph a pixel intensity profile between two point,
										 // like the IRIS cut function. Defined in cut.c

#endif
