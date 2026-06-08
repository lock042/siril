#ifndef _IMAGE_INTERACTIONS_H_
#include <gtk/gtk.h>
#define _IMAGE_INTERACTIONS_H_

#include "core/siril.h"
#include "gui-gtk4/gtk3_event_compat.h"

typedef void (*selection_update_callback)();
typedef void (*star_selection_callback)(pointi);

gboolean update_zoom(gdouble x, gdouble y, double scale);
void update_zoom_fit_button();

void register_selection_update_callback(selection_update_callback f);
void unregister_selection_update_callback(const selection_update_callback f);

void attach_drawingarea_event_controllers(GtkWidget *area);
gboolean new_selection_zone(gpointer user_data);
void delete_selected_area();
void reset_display_offset();
void reset_menu_toggle_button();
void reset_zoom_default();
void update_zoom_label();
gboolean update_zoom_label_idle(gpointer user_data);
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
	MOUSE_ACTION_NONE, // 0
	MOUSE_ACTION_SELECT_REG_AREA, // 1
	MOUSE_ACTION_SELECT_PREVIEW1, // 2
	MOUSE_ACTION_SELECT_PREVIEW2, // 3
	MOUSE_ACTION_DRAW_SAMPLES, // 4
	MOUSE_ACTION_PHOTOMETRY, // 5
	MOUSE_ACTION_GET_COMP_CENTER_COORDINATE, // 6
	MOUSE_ACTION_CUT_SELECT, // 7
	MOUSE_ACTION_CUT_WN1, // 8
	MOUSE_ACTION_CUT_WN2, // 9
	MOUSE_ACTION_DRAW_POLY, // 10
	MOUSE_ACTION_ADD_POLY_TO_MASK, // 11
	MOUSE_ACTION_CLEAR_POLY_FROM_MASK, // 12
	MOUSE_ACTION_SAMPLE_MASK_COLOR // 13
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

/* Siril-owned event-type and modifier-mask enums used in saved configs
 * and in the mouse_action / scroll_action match path.  GdkEventType
 * was renumbered between GTK3 and GTK4 (BUTTON_PRESS went 4→2) and the
 * multi-press types were removed from the public API; binding the
 * on-disk format to either GDK vintage is fragile.  Instead, registered
 * actions store SIRIL_* values that are translated from live GDK input
 * at event-handler entry.
 *
 * Values are pinned to the GTK3 GdkEventType / GdkModifierType numbers
 * (BUTTON_PRESS=4, 2BUTTON_PRESS=5, …; SHIFT=1, CONTROL=4, MOD1=8).
 * That choice means a Siril 1.4 config.1.4.ini — which serialised raw
 * GdkEventType / GdkModifierType integers under GTK3 — loads directly
 * with no value translation, on Linux and Windows.  A 1.4 macOS install
 * needs one tiny fix-up because the primary modifier was saved there as
 * GDK_META_MASK (0x10000000) rather than CONTROL=4; see the OS_OSX
 * branch in mouse_actions_config_to_list / scroll_actions_config_to_list.
 *
 * NEVER renumber these — they are written to disk. */
typedef enum _SirilEventType {
	SIRIL_EVENT_BUTTON_PRESS   = 4,
	SIRIL_EVENT_2BUTTON_PRESS  = 5,
	SIRIL_EVENT_3BUTTON_PRESS  = 6,
	SIRIL_EVENT_BUTTON_RELEASE = 7,
} SirilEventType;

typedef enum _SirilModifierMask {
	SIRIL_MOD_SHIFT   = 1u << 0,  /* matches GDK_SHIFT_MASK */
	SIRIL_MOD_PRIMARY = 1u << 2,  /* matches GDK_CONTROL_MASK (Linux/Windows primary) */
	SIRIL_MOD_ALT     = 1u << 3,  /* matches GDK_MOD1_MASK / GDK_ALT_MASK */
} SirilModifierMask;

static inline SirilEventType siril_event_type_from_n_press(int n_press) {
	return (n_press >= 3) ? SIRIL_EVENT_3BUTTON_PRESS
	     : (n_press == 2) ? SIRIL_EVENT_2BUTTON_PRESS
	     :                  SIRIL_EVENT_BUTTON_PRESS;
}

/* Translate a live GdkModifierType bitmask (from the event controller)
 * into the SIRIL_MOD_* bits the action store uses.  Filters out modifiers
 * we don't care about (NumLock, CapsLock, button-state bits, …) — replaces
 * the historical static `eventmask` AND filter at the comparison site. */
static inline SirilModifierMask siril_mods_from_gdk(GdkModifierType gdk_state) {
	SirilModifierMask m = 0;
	if (gdk_state & GDK_SHIFT_MASK)   m |= SIRIL_MOD_SHIFT;
#ifdef OS_OSX
	if (gdk_state & GDK_META_MASK)    m |= SIRIL_MOD_PRIMARY;
#else
	if (gdk_state & GDK_CONTROL_MASK) m |= SIRIL_MOD_PRIMARY;
#endif
	if (gdk_state & GDK_ALT_MASK)     m |= SIRIL_MOD_ALT;
	return m;
}

typedef struct _mouse_data {
	GtkWidget *widget;
	GdkEventButton *event;       /* legacy; opaque under GTK4 */
	/* Phase 8 transition: the GTK3 callback used to read button/state/x/y
	 * from data->event->X.  Under GTK4 GdkEvent is opaque, so the gesture
	 * dispatcher fills these directly. */
	guint button;
	GdkModifierType state;
	GdkEventType type;
	double x;
	double y;
	double zoom;
	point evpos;
	pointi zoomed;
	gboolean inside;
	mouse_status_enum *mouse_status;
	cut_method *cutting;
} mouse_data;

typedef struct _scroll_data {
	GdkEventScroll *event;            /* legacy; opaque under GTK4 */
	GdkScrollDirection direction_raw; /* GDK4 enum; was event->direction */
	GdkModifierType state;            /* was event->state */
	double x;                         /* was event->x */
	double y;                         /* was event->y */
	double dx;                        /* smooth-scroll dx (GTK4 controller value) */
	double dy;                        /* smooth-scroll dy */
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
	SirilEventType type;
	SirilModifierMask state;
	const mouse_function_metadata *data;
} mouse_action;

/* Contains everything required to define a complete scroll action configured
 * item, and once populated is added to the GList iterated in the drawingarea
 * scroll event handler.
 */
typedef struct _scroll_action {
	SirilScrollDirection direction;
	SirilModifierMask state;
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
mouse_status_enum get_mouse_status();
void initialize_mouse_actions();
void init_draw_poly();
void init_add_poly_to_mask();
void init_clear_poly_from_mask();
void initialize_scroll_actions();
void load_or_initialize_mouse_actions();
void load_or_initialize_scroll_actions();
mouse_action* create_mouse_action(guint button, SirilEventType type, SirilModifierMask state, const mouse_function_metadata *metadata);
scroll_action* create_scroll_action(SirilModifierMask state, const scroll_function_metadata *metadata, SirilScrollDirection direction);
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
