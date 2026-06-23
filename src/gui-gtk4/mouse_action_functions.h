#ifndef _MOUSE_ACTION_FUNCTIONS_H_
#define _MOUSE_ACTION_FUNCTIONS_H_

void cache_widgets();

/* CHECKLIST FOR ADDING NEW MOUSE FUNCTIONS OR SCROLL FUNCTIONS
 *
 * Note: the required typedef structs are defined in image_interactions.h
 *
 * 1. In this header define the action name, the mouse_function_metadata and the
 *    function that the action will call, following the existing entries as templates.
 *
 * 2. Add an entry to the mouse_function_ref or scroll_function_ref enum in
 *    image_interactions.h
 *
 * 3. Initialize the metadata and the function in mouse_action_functions.c
 *
 * 4. Add the pointer to the mouse_function_metadata or scroll_function_metadata
 *    to the array in mouse_action_preferences.c
 */

// Function metadata
#define NAME_MAIN_ACTION N_("Main multifunction action")
extern mouse_function_metadata main_action;
gboolean main_action_click(mouse_data *data);

#define NAME_DRAG_ACTION N_("Drag")
extern mouse_function_metadata drag_action;
gboolean drag_click(mouse_data *data);

#define NAME_MEASURE_ACTION N_("Measure")
extern mouse_function_metadata measure_action;
gboolean measure_click (mouse_data *data);

#define NAME_OPEN_IF_UNLOADED_ACTION N_("Open if nothing is loaded")
extern mouse_function_metadata open_if_unloaded_action;
gboolean open_if_nothing_loaded(mouse_data *data);

#define NAME_ZOOM_ACTION N_("Configurable zoom toggle")
extern mouse_function_metadata zoom_action;
gboolean zoom_function (mouse_data *data);

#define NAME_PHOTOMETRY_BOX_ACTION N_("Photometry box")
extern mouse_function_metadata photometry_box_action;
gboolean photometry_box_set(mouse_data *data);

#define NAME_SECOND_ACTION N_("Mouse Menu / Remove")
extern mouse_function_metadata second_action;
gboolean second_action_click(mouse_data *data);

#define NAME_NULL_ACTION ""
extern mouse_function_metadata null_action;
gboolean mouse_nullfunction(mouse_data *data);

#define NAME_SAVE_ON_CLICK_ACTION N_("Save")
extern mouse_function_metadata save_on_click_action;
gboolean save_on_click(mouse_data *data);

#define NAME_SAVE_AS_ON_CLICK_ACTION N_("Save As")
extern mouse_function_metadata save_as_on_click_action;
gboolean save_as_on_click(mouse_data *data);

#define NAME_SNAPSHOT_ON_CLICK_ACTION N_("Snapshot")
extern mouse_function_metadata snapshot_on_click_action;
gboolean snapshot_on_click(mouse_data *data);

#define NAME_UNDO_ON_CLICK_ACTION N_("Undo")
extern mouse_function_metadata undo_on_click_action;
gboolean undo_on_click(mouse_data *data);

#define NAME_REDO_ON_CLICK_ACTION N_("Redo")
extern mouse_function_metadata redo_on_click_action;
gboolean redo_on_click(mouse_data *data);

#define NAME_FINDSTAR_ON_CLICK_ACTION N_("Find stars")
extern mouse_function_metadata findstar_on_click_action;
gboolean findstar_on_click(mouse_data *data);

#define NAME_PLATESOLVE_ON_CLICK_ACTION N_("Plate solve")
extern mouse_function_metadata platesolve_on_click_action;
gboolean platesolve_on_click(mouse_data *data);

#define NAME_SPCC_ON_CLICK_ACTION N_("SPCC")
extern mouse_function_metadata spcc_on_click_action;
gboolean spcc_on_click(mouse_data *data);

// Scroll functions
#define NAME_SCROLL_ZOOM_ACTION N_("Scroll to zoom")
extern scroll_function_metadata scroll_zooms_action;
gboolean scroll_zooms(scroll_data *data);

#define NAME_SCROLL_SEQ_ACTION N_("Scroll through sequence")
extern scroll_function_metadata scroll_traverses_seq_action;
gboolean scroll_traverses_sequence(scroll_data *data);

#define NAME_SCROLL_TAB_ACTION N_("Scroll through GUI tabs")
extern scroll_function_metadata scroll_changes_tab_action;
gboolean scroll_changes_tab(scroll_data *data);

extern scroll_function_metadata scroll_null_action;
gboolean scroll_nullfunction(scroll_data *data);

#endif
