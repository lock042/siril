#ifndef SRC_ALGOS_ASTROMETRY_SOLVER_H_
#define SRC_ALGOS_ASTROMETRY_SOLVER_H_

#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "io/remote_catalogues.h"
#include "algos/search_objects.h"
#include "registration/matching/degtorad.h"

#define BRIGHTEST_STARS 2500
#define AT_MATCH_CATALOG_NBRIGHT 60
//#define CROP_ALLOWANCE 1.20


typedef enum {
	ERROR_PLATESOLVE = 1,
	ERROR_PHOTOMETRY = 10,
} platesolve_error_type;

typedef enum {
	LIMIT_MAG_AUTO,
	LIMIT_MAG_AUTO_WITH_OFFSET,
	LIMIT_MAG_ABSOLUTE
} limit_mag_mode;

typedef enum {
	SOLVER_SIRIL,
	SOLVER_LOCALASNET
} platesolve_solver;

struct astrometry_data {
	/* user input */
	fits *fit;		// the image
	double pixel_size;	// pixel size in µm
	double focal_length;	// focal length in mm
	platesolve_solver solver;	// the solver being used (siril or localasnet for now)
	siril_catalogue *ref_stars; // siril_catalogue containing query parameters and results
	gboolean coords_forced; // target coords are forced (used for seqplatesolve)
	gboolean blind; // if this flag is false, ref_stars conesearch is done once (mainly when catalogs are not local)
	SirilWorldCS *cat_center;	// starting point for the search
	gboolean downsample;	// downsample image before solving
	gboolean autocrop;	// crop image if fov is larger than 5 degrees
	gboolean flip_image;	// Flip at the end if detected mirrored
	gboolean manual;	// use stars already detected by user, in com.stars
	psf_star **stars;	// alternate source for manual stars (vincent's special)
	limit_mag_mode mag_mode;// automatically limit magnitude of the catalog
	double magnitude_arg;	// if not automatic, use this limit magnitude
	gboolean verbose;	// display all information
	gboolean for_sequence;	// sequence operation, don't free everything
	gchar *filename;	// the name of the file being processed
	int rx_solver;		// width of the image being solved (accounting for downscale if any)
	int ry_solver;		// height of the image being solved (accounting for downscale if any)
	double scalefactor;	// scale factor accounting for downscale if any
	int trans_order; // order of the polynomial fit (if > 1, it includes distortions)
	double blindradius; // radius of the cone if blind

	/* program-processed input, by process_plate_solver_input() */
	double scale;		// scale (resolution) in arcsec per pixel
	double used_fov;	// field of view for the solved image region (arcmin)
	rectangle solvearea;	// area in case of manual selection or autocrop
	gboolean uncentered;	// solvearea is not centered with image
	gboolean asnet_checked;	// local asnet availability already checked

	/* results */
	int ret;		// return value
	gchar *message;		// error message
	gboolean image_flipped;	// image has been flipped
};

void open_astrometry_dialog();
void close_astrometry_dialog();
void process_plate_solver_input(struct astrometry_data *args);
int fill_plate_solver_structure_from_GUI(struct astrometry_data *args);
void wcs_pc_to_cd(double pc[][2], const double cdelt[2], double cd[][2]);
gpointer plate_solver(gpointer p);
double compute_mag_limit_from_fov(double fov_degrees);
gboolean confirm_delete_wcs_keywords(fits *fit);
void reframe_astrometry_data(fits *fit, Homography H);

void set_focal_and_pixel_pitch();

void start_sequence_astrometry(sequence *seq, struct astrometry_data *args);

gboolean asnet_is_available();
void reset_asnet_version();

/* for the GUI */
double get_resolution(double focal, double pixel);
double get_radius_deg(double resolution, int rx, int ry);

gboolean has_any_keywords();
SirilWorldCS *get_eqs_from_header(fits *fit);
double get_fov_arcmin(double resolution, int rx, int ry);

const char *catalog_to_str(siril_cat_index cat);

/* from the GUI */
gboolean end_process_catsearch(gpointer p);
void update_coords();
gboolean end_plate_solver(gpointer p);

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data);

#endif /* SRC_ALGOS_ASTROMETRY_SOLVER_H_ */
