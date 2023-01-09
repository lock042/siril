#ifndef SRC_ALGOS_ASTROMETRY_SOLVER_H_
#define SRC_ALGOS_ASTROMETRY_SOLVER_H_

#include "core/siril.h"
#include "core/siril_world_cs.h"

#include "registration/matching/degtorad.h"

#define BRIGHTEST_STARS 2500
#define AT_MATCH_CATALOG_NBRIGHT   60
//#define CROP_ALLOWANCE 1.20

#define RADtoASEC (3600.0 * 180.0 / M_PI)

#define CDSSESAME "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame"
#define VIZIERSESAME "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame"
#define SIMBADSESAME "http://simbad.u-strasbg.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=TSV&query=SELECT basic.OID, ra, dec, main_id FROM basic JOIN ident ON ident.oidref = oid WHERE id ='"
#define EPHEMCC "http://ssp.imcce.fr/webservices/miriade/api/ephemcc.php?"
#define SKYBOT "http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?"


typedef enum {
	TYCHO2,
	NOMAD,
	GAIADR3,
	PPMXL,
	BRIGHT_STARS,
	APASS,
	LOCAL = 99
} online_catalog;

typedef enum {
	RESOLVER_UNSET = -1,
	RESOLVER_NED = 0,
	RESOLVER_SIMBAD,
	RESOLVER_VIZIER,
	RESOLVER_NUMBER,
} resolver_t;


typedef enum {
	QUERY_SERVER_CDS,
	QUERY_SERVER_VIZIER,
	QUERY_SERVER_SIMBAD,
	QUERY_SERVER_EPHEMCC,
	QUERY_SERVER_SKYBOT,
} query_server;

struct astrometry_data {
	/* user input */
	fits *fit;		// the image
	double pixel_size;	// pixel size in µm
	double focal_length;	// focal length in mm
	online_catalog onlineCatalog;	// choice of catalog for the plate solve
	SirilWorldCS *cat_center;	// starting point for the search
	gboolean downsample;	// downsample mage before solving
	gboolean autocrop;	// crop image if fov is larger than 5 degrees
	gboolean flip_image;	// Flip at the end if detected mirrored
	gboolean manual;	// use stars already detected by user, in com.stars
	gboolean auto_magnitude;// automatically limit magnitude of the catalog
	double forced_magnitude;// if not automatic, use this limit magnitude
	gboolean for_photometry_cc;	// proceeed to PCC after a successful plate solve
	struct photometric_cc_data *pcc;// PCC configuration

	/* program-processed input, by process_plate_solver_input() */
	double limit_mag;	// limit magnitude to sear for in the catalog
	double scale;		// scale (resolution) in arcsec per pixel
	double used_fov;	// field of view for the solved image region (arcmin)
	gboolean use_local_cat;	// use local catalogues if installed
	GFile *catalog_file;	// downloaded file containing raw catalog data
	gchar *catalogStars;	// file name of the projected catalog
	rectangle solvearea;	// area in case of manual selection or autocrop
	gboolean uncentered;	// solvearea is not centered with image

	/* results */
	int ret;		// return value
	gchar *message;		// error message
	gboolean image_flipped;	// image has been flipped
	SirilWorldCS *new_center; // the image center found by the solve
};

struct sky_object {
	gchar *name;
	double radius;
	int maxRecords;
	SirilWorldCS *world_cs;
	point imageCenter;
	gboolean south;
};

void open_astrometry_dialog();
gchar *search_in_online_catalogs(const gchar *object, query_server server);
void process_plate_solver_input(struct astrometry_data *args);
int fill_plate_solver_structure_from_GUI(struct astrometry_data *args);
void wcs_cd_to_pc(double cd[][2], double pc[][2], double cdelt[2]);
void wcs_pc_to_cd(double pc[][2], const double cdelt[2], double cd[][2]);
gpointer match_catalog(gpointer p);
double compute_mag_limit_from_fov(double fov_degrees);

gboolean confirm_delete_wcs_keywords(fits *fit);
void flip_bottom_up_astrometry_data(fits *fit);
void reframe_astrometry_data(fits *fit, Homography H);

void set_focal_and_pixel_pitch();

/* for the GUI */
double get_resolution(double focal, double pixel);
double get_radius_deg(double resolution, int rx, int ry);
void free_Platedobject();
int parse_content_buffer(char *buffer, struct sky_object *obj);
gchar *search_in_online_conesearch(struct astrometry_data *args);
gboolean has_nonzero_coords();
gboolean has_any_keywords();
SirilWorldCS *get_eqs_from_header(fits *fit);
GFile *download_catalog(online_catalog onlineCatalog, SirilWorldCS *catalog_center, double radius, double mag);
gchar *get_catalog_url(SirilWorldCS *center, double mag_limit, double dfov, int type);

/* from the GUI */
void update_coords();
gboolean end_plate_solver(gpointer p);

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data);

#endif /* SRC_ALGOS_ASTROMETRY_SOLVER_H_ */
