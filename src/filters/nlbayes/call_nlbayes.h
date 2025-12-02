#ifndef SRC_CALL_NLBAYES_H_
#define SRC_CALL_NLBAYES_H_

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

typedef struct denoise_args {
	destructor destroy_fn;  // Must be first member for generic_img_args compatibility
	fits *fit;
	float modulation;
	unsigned int sos;
	int da3d;
	float rho;
	gboolean do_anscombe;
	gboolean do_cosme;
	gboolean suppress_artefacts;
	gboolean previewing;
} denoise_args;

/* Allocator and destructor functions */
EXTERNC struct denoise_args *new_denoise_args();
EXTERNC void free_denoise_args(void *args);

/* Image processing hook */
EXTERNC int denoise_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

/* Idle functions */
EXTERNC gboolean denoise_preview_idle(gpointer p);
EXTERNC gboolean denoise_apply_idle(gpointer p);

/* Core denoising function */
EXTERNC int do_nlbayes(fits *fit, const float modulation, unsigned sos, int da3d,
                       const float rho, const gboolean do_anscombe);

/* Dialog close function */
EXTERNC void close_denoise();

#ifdef __cplusplus
}
#endif

#endif /* SRC_CALL_NLBAYES_H_ */
