#ifndef _HIST_H_
#define _HIST_H_

#include "filters/mtf.h"
#include "filters/ght.h"

#define NO_STRETCH_SET_YET 0
#define HISTO_STRETCH 1
#define GHT_STRETCH 2

struct mtf_data {
	void (*destroy_fn)(void *args);  // First member - destructor
	fits *fit;
	sequence *seq;
	gboolean linked;
	struct mtf_params params;
	struct mtf_params uparams[3]; // for unlinked stretch
	char *seqEntry;
	gboolean auto_display_compensation;
	gboolean is_preview;
};

struct ght_data {
	void (*destroy_fn)(void *args);  // First member - destructor
	fits *fit;
	sequence *seq;
	struct ght_params *params_ght;
	char *seqEntry;
	gboolean auto_display_compensation;
	gboolean is_preview;
};

typedef enum {
	SCALE_LOW,
	SCALE_MID,
	SCALE_HI
} ScaleType;

struct mtf_data* create_mtf_data();
struct ght_data* create_ght_data();
void destroy_mtf_data(void *args);
void destroy_ght_data(void *args);

gchar *invmtf_log_hook(gpointer p, log_hook_detail detail);
gchar *mtf_log_hook(gpointer p, log_hook_detail detail);
gchar *ght_log_hook(gpointer p, log_hook_detail detail);

int invmtf_single_image_hook(struct generic_img_args *args, fits *fit, int threads);
int mtf_single_image_hook(struct generic_img_args *args, fits *fit, int threads);
int ght_single_image_hook(struct generic_img_args *args, fits *fit, int threads);
gboolean mtf_single_image_idle(gpointer p);
gboolean ght_single_image_idle(gpointer p);

gsl_histogram* computeHisto(fits*, int);
gsl_histogram* computeHisto_Selection(fits*, int, rectangle *);
void histo_change_between_roi_and_image();
void compute_histo_for_fit(fits *thefit);
void invalidate_gfit_histogram();
void update_gfit_histogram_if_needed();
void apply_histo_cancel();
void toggle_histogram_window_visibility(int _invocation);

gboolean on_button_histo_close_clicked(GtkButton *button, gpointer user_data); // callback needed

void on_histoMidEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoShadEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoHighEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histo_toggled(GtkToggleButton *togglebutton, gpointer user_data);

void apply_mtf_to_sequence(struct mtf_data *mtf_args);
void apply_ght_to_sequence(struct ght_data *ght_args);

void erase_histo_display(cairo_t *cr, int width, int height);
void draw_grid(cairo_t *cr, int width, int height);
void display_scale(cairo_t *cr, int width, int height);
void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width, int height, double zoomH, double zoomV, gboolean isOrig, gboolean is_log);
#endif
