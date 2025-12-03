#ifndef SRC_FILTERS_BILAT_H_
#define SRC_FILTERS_BILAT_H_

typedef enum {
	EP_BILATERAL,
	EP_GUIDED
} ep_filter_t;

struct epfargs {
	destructor destroy_fn;  // Must be first member
	fits *fit; // just a reference, not freed in free_epf_args
	fits *guidefit;
	gboolean guide_needs_freeing;
	double d;
	double sigma_col;
	double sigma_space;
	double mod;
	ep_filter_t filter;
	gboolean verbose;
	gboolean applying;
};

gchar *epf_log_hook(gpointer p, log_hook_detail detail);

/* Allocator and destructor functions */
struct epfargs *new_epf_args();
void free_epf_args(void *args);

/* Image processing hook */
int epf_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

/* Idle functions */
gboolean epf_preview_idle(gpointer p);
gboolean epf_apply_idle(gpointer p);

void epf_change_between_roi_and_image();
void apply_epf_cancel();

#endif /* SRC_FILTERS_BILAT_H_ */
