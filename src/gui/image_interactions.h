#ifndef _IMAGE_INTERACTIONS_H_
#define _IMAGE_INTERACTIONS_H_

#include "core/siril.h"

typedef void (*selection_update_callback)();
typedef void (*star_selection_callback)(pointi);

gboolean update_zoom(gdouble x, gdouble y, double scale);
void update_zoom_fit_button();

void register_selection_update_callback(selection_update_callback f);
void unregister_selection_update_callback(const selection_update_callback f);

void set_mouse_event_mask();
gboolean new_selection_zone(gpointer user_data);
void delete_selected_area();
void reset_display_offset();
void reset_menu_toggle_button();
void reset_zoom_default();
void update_zoom_label();
void enforce_ratio_and_clamp();
gboolean display_quick_photo();
GdkModifierType get_primary();
gboolean is_over_the_left_side_of_sel(pointi zoomed, double zoom);
gboolean is_over_the_right_side_of_sel(pointi zoomed, double zoom);
gboolean is_over_the_bottom_of_sel(pointi zoomed, double zoom);
gboolean is_over_the_top_of_sel(pointi zoomed, double zoom);
gboolean is_inside_of_sel(pointi zoomed, double zoom);

/* mouse behaviour */
/* Do not assign a value to x_REF_MAX, it will be set automatically.
 * NEVER change the assigned number of these enums, as they are used to
 * store config. Changing them will lead to unexpected behaviour (which is
 * why they are specifically defined with values, so if we ever need to delete
 * a function and its enum (unlikely, but...), the rest will retain their values.
 */
typedef enum _mouse_function_ref {
	MOUSE_REF_NULL = 0,
	MOUSE_REF_MAIN_ACTION = 1,
	MOUSE_REF_DRAG = 2,
	MOUSE_REF_MEASURE = 3,
	MOUSE_REF_OPEN_UNLOADED = 4,
	MOUSE_REF_ZOOM = 5,
	MOUSE_REF_PHOTOMETRY_BOX = 6,
	MOUSE_REF_SECOND_ACTION = 7,
	MOUSE_REF_SAVE = 8,
	MOUSE_REF_SAVE_AS = 9,
	MOUSE_REF_SNAPSHOT = 10,
	MOUSE_REF_UNDO = 11,
	MOUSE_REF_REDO = 12,
	MOUSE_REF_FINDSTAR = 13,
	MOUSE_REF_PLATESOLVE = 14,
	MOUSE_REF_SPCC = 15,
// Insert additional mouse function references here
	MOUSE_REF_MAX
} mouse_function_ref;

typedef enum _scroll_function_ref {
	SCROLL_REF_NULL = 0,
	SCROLL_REF_ZOOM = 1,
	SCROLL_REF_SEQ = 2,
	SCROLL_REF_TABS = 3,
// Insert additional scroll function references here
	SCROLL_REF_MAX
} scroll_function_ref;

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

typedef enum _SirilScrollDirection {
	MOUSE_HORIZ_SCROLL,
	MOUSE_VERTICAL_SCROLL,
	MOUSE_SMOOTH_SCROLL
} SirilScrollDirection;

typedef struct _mouse_data {
	GtkWidget *widget;
	GdkEventButton *event;
	double zoom;
	point evpos;
	pointi zoomed;
	gboolean inside;
	mouse_status_enum *mouse_status;
	cut_method *cutting;
} mouse_data;

typedef struct _scroll_data {
	GdkEventScroll *event;
	SirilScrollDirection direction;
} scroll_data;

typedef gboolean (*mouse_function)(mouse_data *);
typedef gboolean (*scroll_function)(scroll_data *);

/* Contains the immutable related elements of a mouse action. A
 * mouse_function_metadata can be added to a mouse_action along with button
 * number, event type and modifiers to form a completely configured mouse action
 */
typedef struct _mouse_function_metadata {
	mouse_function function;
	const gchar* name;
	const gchar* tooltip;
	gboolean doubleclick_compatible;
	mouse_function_ref reference;
} mouse_function_metadata;

/* Contains the immutable related elements of a scroll action. A
 * scroll_function_metadata can be added to a scroll_action along with
 * direction and modifiers to form a completely configured scroll action
 */
typedef struct _scroll_function_metadata {
	scroll_function function;
	const gchar* name;
	const gchar* tooltip;
	scroll_function_ref reference;
} scroll_function_metadata;

/* Contains everything required to define a complete mouse action configured
 * item, and once populated is added to the GList iterated in the drawingarea
 * button press event handler.
 */
typedef struct _mouse_action {
	guint button;
	GdkEventType type;
	GdkModifierType state;
	const mouse_function_metadata *data;
} mouse_action;

/* Contains everything required to define a complete scroll action configured
 * item, and once populated is added to the GList iterated in the drawingarea
 * scroll event handler.
 */
typedef struct _scroll_action {
	SirilScrollDirection direction;
	GdkModifierType state;
	const scroll_function_metadata *data;
} scroll_action;

/* Contains everything required to define a release action. In this case no
 * metadata is required as the action is automatically added from the button
 * press function, so we only need to configure the button and function pointer.
 */
typedef struct _release_action {
	guint button;
	mouse_function function;
} release_action;

void init_mouse();
void initialize_mouse_actions();
void initialize_scroll_actions();
void load_or_initialize_mouse_actions();
void load_or_initialize_scroll_actions();
mouse_action* create_mouse_action(guint button, GdkEventType type, GdkModifierType state, const mouse_function_metadata *metadata);
scroll_action* create_scroll_action(GdkModifierType state, const scroll_function_metadata *metadata, SirilScrollDirection direction);
extern mouse_status_enum mouse_status;	// defined in registration_preview.c
void register_release_callback(mouse_function function, guint button);
void clear_release_callback();
GSList *mouse_actions_config_to_list(GSList *config);
GSList *scroll_actions_config_to_list(GSList *config);

GSList *mouse_actions_list_to_config(GSList *actions);
GSList *scroll_actions_list_to_config(GSList *actions);

#define CLICK_TEXT N_("Single click")
#define DOUBLE_CLICK_TEXT N_("Double click")
#define SHIFT_TEXT N_("Shift")
#ifdef __APPLE__
#define CTRL_TEXT N_("Compose")
#define CTRL_SHIFT_TEXT N_("Compose Shift")
#else
#define CTRL_TEXT N_("Ctrl")
#define CTRL_SHIFT_TEXT N_("Ctrl Shift")
#endif
#define NO_MODIFIER_TEXT N_("None")
#define SCROLL_HORIZ_TEXT N_("Horizontal")
#define SCROLL_VERTICAL_TEXT N_("Vertical")
#define NEW_ENTRY_TEXT N_("New entry")

#endif
