#ifndef SRC_ALGOS_ASTROMETRY_SOLVER_H_
#define SRC_ALGOS_ASTROMETRY_SOLVER_H_

#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "io/remote_catalogues.h"
#include "algos/search_objects.h"
#include "registration/matching/degtorad.h"

#define BRIGHTEST_STARS 2000
#define AT_MATCH_CATALOG_NBRIGHT 60

typedef enum {
	LIMIT_MAG_AUTO,
	LIMIT_MAG_AUTO_WITH_OFFSET,
	LIMIT_MAG_ABSOLUTE
} limit_mag_mode;

typedef enum {
	SOLVER_SIRIL,
	SOLVER_LOCALASNET
} platesolve_solver;

typedef enum {
	FORCED_CENTER,
	FORCED_PIXEL,
	FORCED_FOCAL
} forced_metadata_t;

typedef enum {
	SOLVE_LINONLY = -2, // siril or asnet solvers returned a linear solution only instead of the required SIP order
	SOLVE_LASTSOLVE = -1, // siril internal solver found a solution but the last centered solve did not succeed
	// generic
	SOLVE_OK,
	SOLVE_NO_MATCH, //solution->message = g_strdup(_("Could not match stars from the catalogue"));
	SOLVE_CANCELLED, //solution->message = g_strdup(_("Cancelled"));
	// platesolve common errors
	SOLVE_DOWNSAMPLE, //args->message = g_strdup(_("Not enough memory"));
	SOLVE_NOTENOUGHSTARS, //		args->message = g_strdup_printf(_("There are not enough stars picked in the image. least %d are needed."), AT_MATCH_STARTN_LINEAR);
	// siril solver
	SOLVE_INVALID_TRANS, 		//solution->message = g_strdup(_("Transformation matrix is invalid, solve failed"));
	SOLVE_PROJ, //solution->message = g_strdup(_("Reprojecting catalog failed."));
	SOLVE_TRANS_UPDATE, //g_strdup(_("Updating trans failed."));
	SOLVE_NEAR_NO_MATCH, // siril internal solver did not find a solution using near solver
	SOLVE_NEAR_TIMEOUT, // siril internal near solver timed out
	SOLVE_NEAR_THREADS, // siril internal near solver problem allocating threads

	// asnet solver
	SOLVE_ASNET_PROC, // asnet errors prior to solve

} siril_solve_errcode;

struct astrometry_data {
	/* user input */
	fits *fit;		// the image
	double pixel_size;	// pixel size in µm
	double focal_length;	// focal length in mm
	platesolve_solver solver;	// the solver being used (siril or localasnet for now)
	siril_catalogue *ref_stars; // siril_catalogue containing query parameters and results
	SirilWorldCS *cat_center;	// starting point for the search
	gboolean downsample;	// downsample image before solving
	gboolean autocrop;	// crop image if fov is larger than 5 degrees
	gboolean flip_image;	// Flip at the end if detected mirrored
	gboolean manual;	// use stars already detected by user, in com.stars
	gboolean nocache;	// solve each image with its metadata for a sequence (do not use a common catalogue)
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
	double searchradius; // radius of the cone if nearsearch in degrees
	gboolean forced_metadata[3]; // flags using for seq, to indicate if center, pixel and focal where forced
	gboolean force; // flag to force solving again already solved images from sequence
	gchar *distofilename; // flag to save distortion file

	/* program-processed input, by process_plate_solver_input() */
	double scale;		// scale (resolution) in arcsec per pixel
	double used_fov;	// field of view for the solved image region (arcmin)
	rectangle solvearea;	// area in case of manual selection or autocrop
	gboolean uncentered;	// solvearea is not centered with image
	gboolean asnet_checked;	// local asnet availability already checked
	gboolean near_solve; // if this flag is false, ref_stars conesearch is done once (mainly when catalogs are not local)
	gboolean asnet_blind_pos; // if this flag is true, no position is passed to asnet, the solve is blind in position
	gboolean asnet_blind_res; // if this flag is true, no resolution is passed to asnet, the solve is blind in scale
	int numthreads; //nb of threads that can be used

	/* results */
	int ret;		// return value
	gboolean image_flipped;	// image has been flipped
	int seqprogress; // used to log number of solved images when processing a sequence
	int seqskipped; // used to log number of skipped images (already solved) when processing a sequence
	gboolean update_reg; // used to compute regdata when processing a sequence
	int layer; // the layer which will contain the regdata
	struct starfinder_data *sfargs;		// star finder configuration for peaker
	struct wcsprm *WCSDATA;
};

typedef struct {
	psf_star **stars;
	siril_catalogue *search_cat;
	int nb_stars;
	double scale;
	point *centers;
	double radius;
	int N;
	gboolean verbose;
	// modified by the pool workers
	gint progress;
	gint solved;
	point *center;
} near_solve_data;

void open_astrometry_dialog();
void close_astrometry_dialog();
void process_plate_solver_input(struct astrometry_data *args);
int fill_plate_solver_structure_from_GUI(struct astrometry_data *args);
void wcs_pc_to_cd(double pc[][2], const double cdelt[2], double cd[][2]);
gpointer plate_solver(gpointer p);
double compute_mag_limit_from_position_and_fov(double ra, double dec, double fov_degrees, int Nstars);
gboolean confirm_delete_wcs_keywords(fits *fit);
void reframe_wcs(struct wcsprm *wcs, Homography *H);
void reframe_astrometry_data(fits *fit, Homography *H);
void flip_bottom_up_astrometry_data(fits *fit);
void update_wcsdata_from_wcs(fits *fit);

void init_astrometry();
void reset_astrometry_checks();
void initialize_ips_dialog();
void free_astrometry_data(struct astrometry_data *args);

void start_sequence_astrometry(sequence *seq, struct astrometry_data *args);

gboolean asnet_is_available();
void reset_asnet_version();

/* for the GUI */
double get_resolution(double focal, double pixel);
double get_radius_deg(double resolution, int rx, int ry);

gboolean has_any_keywords();
SirilWorldCS *get_eqs_from_header(fits *fit);
double get_fov_arcmin(double resolution, int rx, int ry);
gchar *platesolve_msg(struct astrometry_data *args);

/* from the GUI */
gboolean end_process_catsearch(gpointer p);
void update_coords();
gboolean end_plate_solver(gpointer p);
gboolean end_platesolve_sequence(gpointer p);

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data);

#endif /* SRC_ALGOS_ASTROMETRY_SOLVER_H_ */
