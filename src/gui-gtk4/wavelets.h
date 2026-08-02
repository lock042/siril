#ifndef SRC_GUI_WAVELETS_H_
#define SRC_GUI_WAVELETS_H_

void apply_wavelets_cancel(void);

/* Invalidate the tool's computed decomposition when the active FLIS layer
 * changes; called from the active-layer reconciler (gui_iface_impl.c). */
void wavelets_active_layer_changed(void);

#endif /* SRC_GUI_WAVELETS_H_ */
