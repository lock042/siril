#if !defined(MISC_H)
#define MISC_H

#include "core/siril.h"

/*
 *
 * FILE: misc.h
 *
 * DESCRIPTION:
 * Definitions for structures and functions used in support
 * of the matching code.  Much of this is material that appears
 * in SHIVA, the SDSS environment.
 *
 */


#define SH_SUCCESS        0          /* indicates that all went well */
#define SH_GENERIC_ERROR  1          /* indicates that error occurred */


	/* a buffer used for parsing command-line arguments */
#define CMDBUFLEN       500

   /* ignore any lines in input files that start with this */
#define COMMENT_CHAR   '#'

   /* data files can have this many data columns, at most */
#define MAX_DATA_COL    30

   /* each column in the data file can have at most this many characters */
#define MAX_COL_LENGTH 50


   /* possibilities for the "order" field of TRANS structures */
#define AT_TRANS_UNDEFINED   0           /* when not specified */
#define AT_TRANS_LINEAR      1           /* linear terms only */
#define AT_TRANS_QUADRATIC   2           /* linear plus quadratic */
#define AT_TRANS_CUBIC       3           /* linear plus quadratic plus cubic */
#define AT_TRANS_QUARTIC     4           /* linear plus quadratic plus cubic plus quartic*/
#define AT_TRANS_QUINTIC     5           /* linear plus quadratic plus cubic plus quartic plus quintic*/

   /*
    * little wrappers around 'siril_malloc' and 'siril_free' functions.
    */

void *
shMalloc(int nbytes);

void
shFree(void *vptr);

void
shError(char *format, ...);

void
shFatal(char *format, ...);

void
shDebug(int level, char *format, ...);

   /*
    * This is a preprocessor macro, which acts like a function.  The user
    * calls it with one argument, a condition which should evaluate to
    * 0 or 1.  If it evaluates to 1, then nothing happens; but if it
    * evaluates to 0, then the program prints an error message, giving
    * location of the error, and halts execution with an error code.
    *
    * Thus, one typically uses it to make a 'sanity check', such as
    * making sure that a pointer (which really, really should have
    * a valid value) isn't NULL:
    *
    *     fp = fopen("file", "r");
    *        ...
    *     shAssert(fp != NULL);
    *     fgets(line, LINELEN, fp);
    *        ...
    */
#define shAssert(x) if((x)!=1){fprintf(stderr,"assertion fails in file %s, line %d\n",__FILE__,__LINE__);}


   /*
    * This is a generic transformation from one coordinate system to another.
    * Given the measured (x, y), the transformed coords (x', y') are
    * calculated like this:
    *
    *   if linear terms only:
    *
    *       x' = X00 + X10*x + X01*y
    *       y' = Y00 + Y10*x + Y01*y
    *
    *   if linear plus quadratic terms,
    *
    *      x' =  X00 + X10*x + X01*y + X20*xx + X11*xy + X02*yy
    *      y' =  Y00 + Y10*x + Y01*y + Y20*xx + Y11*xy + Y02*yy
    *
    *   if linear plus quadratic plus cubic,
    *
    *      x' =  X00 + X10*x + X01*y + X20*xx + X11*xy + X02*yy + X30*xxx + X21*xxy + X12*xyy + X03*yyy
    *      y' =  Y00 + Y10*x + Y01*y + Y20*xx + Y11*xy + Y02*yy + Y30*xxx + Y21*xxy + Y12*xyy + Y03*yyy
    * 
    *   if linear plus quadratic plus cubic plus quartic,
    *
    *      x' =  X00 + X10*x + X01*y + X20*xx + X11*xy + X02*yy + X30*xxx + X21*xxy + X12*xyy + X03*yyy + X40*xxxx + X31*xxxy + X22*xxyy + X13*xyyy + X04*yyyy
    *      y' =  Y00 + Y10*x + Y01*y + Y20*xx + Y11*xy + Y02*yy + Y30*xxx + Y21*xxy + Y12*xyy + Y03*yyy + Y40*xxxx + Y31*xxxy + Y22*xxyy + Y13*xyyy + Y04*yyyy
    *
    *
    *  The 'order' field of the TRANS structure signals which
    *  of the above cases is to be used.
	 *
	 *  The 'nr' field contains the number of pairs ultimately used to
	 *  determine the transform between the coordinate system.
	 *
	 *  The 'nm' field contains the number of pairs which -- after
	 *  the transform has been determined and applied -- match up,
	 *  and can be used to determine the quality of the fit.
	 *
	 *  The 'sig' field holds the standard deviation of the separation
	 *  between matched stars in the two sets, after they have been
	 *  transformed into the coordinate system of the second list.
	 *  Only items which were used to derive the TRANS are included
	 *  in this calculation.
	 *
	 *  The 'sx' and 'sy' fields hold the standard deviation of the
	 *  1-D separations between corresponding items in the two sets,
	 *  after they have been transformed into the coordinate system
	 *  of the second list.  A single iteration of 3-sigma clipping
	 *  is used in the calculation.
	 *  If the user specifies the "recalc" option, then the "sx" and "sy"
	 *  values will be based on all items which can be matched in
	 *  the two frames -- not just those few which were used to derive
	 *  the TRANS.
    */

typedef struct Trans {
  int id;
  int order;
  double x00, x10, x01, x20, x11, x02, x30, x21, x12, x03, x40, x31, x22, x13, x04, x50, x41, x32, x23, x14, x05;
  double y00, y10, y01, y20, y11, y02, y30, y21, y12, y03, y40, y31, y22, y13, y04, y50, y41, y32, y23, y14, y05;
  int nr;
  int nm;
  double sig;
  double sx;
  double sy;
} TRANS;

void atTransOrderSet(int order);
int atTransOrderGet(void);
TRANS *atTransNew(void);
Homography *atHNew(void);
void print_trans(TRANS *trans);
void atTransDel(TRANS *trans);
void atHDel(Homography *H);
void print_H(Homography *H);

   /*
    * create a new s_star structure
    */

struct s_star *
atStarNew(double x, double y, double mag, double BV);

   /*
    * converts a psf_star list to a s_star list usable by the match functions
    */

int
get_stars(psf_star **s, int n, int *num_stars, struct s_star **list);

   /*
    * updates stars positions using an updated psf_star catalog
    */
void
update_stars_positions(struct s_star **old_list, int n_old, psf_star **s);

   /*
    * creates a list of stars on a mesh
    */

struct s_star *
create_grid_list(int rx, int ry, int nbpoints);

void
free_stars(struct s_star **head);

int
is_blank(char *line);

int
get_value(char *str, double *val);


#endif    /* MISC_H */
