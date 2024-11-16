#ifndef SRC_FILTERS_BILAT_H_
#define SRC_FILTERS_BILAT_H_

typedef enum {
	EP_BILATERAL,
	EP_GUIDED
} ep_filter_t;

struct epfargs {
	fits *fit;
	fits *guidefit;
	gboolean guide_needs_freeing;
	double d;
	double sigma_col;
	double sigma_space;
	double mod;
	ep_filter_t filter;
	gboolean verbose;
};

gpointer epfhandler(gpointer args);
gpointer epf_filter (gpointer args);
int edge_preserving_filter(struct epfargs *args);
void epf_change_between_roi_and_image();
void apply_epf_cancel();

#endif /* SRC_GUI_ASINH_H_ */
