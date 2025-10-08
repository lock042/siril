#ifndef _COMPOSITION_H_
#define _COMPOSITION_H_

void open_compositing_window();
void reset_compositing_module();
gboolean valid_rgbcomp_seq();
int crop_rgbcomp_seq();
#ifndef _PROCESSING_H_
struct generic_seq_args;
#endif
#ifndef SRC_ALGOS_GEOMETRY_H_
gint64 crop_compute_size_hook(struct generic_seq_args *args, int nb_frames);
int crop_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads);
int crop_finalize_hook(struct generic_seq_args *args);
#endif

#endif
