#ifndef SRC_GUI_MPP_AP_EDITOR_H_
#define SRC_GUI_MPP_AP_EDITOR_H_

#include <gtk/gtk.h>

/* True if the AP editor dialog is currently visible. Used by the mouse
 * handlers in mouse_action_functions.c to gate edit-mode behaviour. */
gboolean mpp_ap_editor_is_open(void);

/* Index of the AP currently being dragged; -1 when not dragging.
 * Mouse handlers set this on click, motion reads, release callback clears. */
int  mpp_ap_editor_get_drag_idx(void);
void mpp_ap_editor_set_drag_idx(int idx);

/* Index of the AP under the cursor for hover-highlighting; -1 when not
 * hovering. Set by the motion handler when in MOUSE_ACTION_EDIT_APS mode,
 * cleared on dialog hide. The overlay draw uses this to colour the
 * targeted AP differently so the user knows which one a click would
 * affect. */
int  mpp_ap_editor_get_hover_idx(void);
void mpp_ap_editor_set_hover_idx(int idx);

/* Recompute the "Current APs: N" label from the cached run. Safe to call
 * before the dialog has ever been shown (no-op until widget cached). */
void mpp_ap_editor_refresh_count_label(void);

#endif /* SRC_GUI_MPP_AP_EDITOR_H_ */
