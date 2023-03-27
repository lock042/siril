

/*-------------------------------------------------------------------------*/
/**
  @file     gnuplot_i.c
  @author   N. Devillard
  @date Sep 1998
  @version  $Revision: 2.10 $
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.

*/
/*--------------------------------------------------------------------------*/

/*
    $Id: gnuplot_i.c,v 2.10 2003/01/27 08:58:04 ndevilla Exp $
    $Author: ndevilla $
    $Date: 2003/01/27 08:58:04 $
    $Revision: 2.10 $
 */

/*
    $New version modified for Siril $
    $Author: Cyril Richard $
    $Date: 2017
    $Revision: 2.12 $
 */

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/

#include "gnuplot_i.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <unistd.h>

#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#include <gio/gwin32inputstream.h>
#else
#include <sys/types.h> // for waitpid(2)
#include <sys/wait.h> // for waitpid(2)
#include <gio/gunixinputstream.h>
#endif

#include <glib.h> // g_get_tmp_dir
#include <glib/gstdio.h>

#include "gui/plot.h"
#include "core/siril_log.h"

#ifdef _WIN32
#define GNUPLOT_BIN "gnuplot.exe"
#else
#define GNUPLOT_BIN "gnuplot"
#endif

static gboolean gnuplot_is_in_path = FALSE;

/*********************** finding gnuplot first **********************/
static gchar *siril_get_gnuplot_bin() {
    if (gnuplot_is_in_path)
        return g_strdup(GNUPLOT_BIN);
    return g_build_filename(com.pref.gnuplot_dir, GNUPLOT_BIN, NULL);
}

#if defined (_WIN32) || defined(OS_OSX)
gboolean gnuplot_is_available() {
    gchar *bin = siril_get_gnuplot_bin();
    if (!bin) return FALSE;

    gboolean is_available = g_file_test(bin, G_FILE_TEST_EXISTS);
    g_free(bin);

    return is_available;
}

#else
/* returns true if the command gnuplot is available */
gboolean gnuplot_is_available() {
    gchar *str = g_strdup_printf("%s -e > /dev/null 2>&1", GNUPLOT_BIN);

    int retval = system(str);
    g_free(str);
    if (WIFEXITED(retval)) {
        gnuplot_is_in_path = TRUE;
        return 0 == WEXITSTATUS(retval);
    }

    gchar *bin = siril_get_gnuplot_bin();
    gboolean is_available = g_file_test(bin, G_FILE_TEST_EXISTS);
    g_free(bin);

    return is_available;
}
#endif

/*---------------------------------------------------------------------------
                                Defines
 ---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------
                          Prototype Functions
 ---------------------------------------------------------------------------*/

/**
 * Creates a temporary file name for writing
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 *
 * @return char const * Pointer to file name string.
 */
char const * gnuplot_tmpfile(gnuplot_ctrl * handle);

/**
 * Plot a temporary file.
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 * @param tmp_filename
 * @param title
 */
void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename, char const* title, int x_offset);

/*-------------------------------------------------------------------------*/
/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   void

  Closes gnuplot by calling an exit command and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.
  This is meant to be called when plot are not displayed

 */
/*--------------------------------------------------------------------------*/

void gnuplot_exit(gnuplot_ctrl * handle)
{
    gnuplot_cmd(handle, "exit");
    free(handle);
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    gnuplot tmpfile watcher. Monitors for tmp files that are finished
			with and closes them when required.
  @param    gpointer user_data. Pointer to gnuplot_ctrl handle.
  @return   GINT_TO_POINTER(1)

  NOTE:
  TODO:		This function currently triggers an Address Sanitizer error with
			debug builds. It runs correctly without Address Sanitizer active.
			Need to see if there is a way to fix its interaction with AS.
*/
/*--------------------------------------------------------------------------*/

gpointer tmpwatcher (gpointer user_data) {
	printf("tmpwatcher started\n");
	gnuplot_ctrl* handle = (gnuplot_ctrl*) user_data;
	GInputStream *stream = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(handle->child_fd), FALSE);
#else
	stream = g_unix_input_stream_new(handle->child_fd, FALSE);
#endif
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL))) {
		printf("No. of tmp files: %d\n", handle->ntmp);
		printf("Buffer: %s\n", buffer);
		if (!handle->ntmp)
			continue;
		gchar *arg = buffer;
		if (g_str_has_prefix(buffer, "Done ")) {
			printf("Received Done message ntmp = %d\n", handle->ntmp);
			arg += 5;
			for (int i = 0 ; i < handle->ntmp ; i++) {
				printf("%s / %s\n", arg, handle->tmp_filename_tbl[i]);
				if (!g_strcmp0(arg, handle->tmp_filename_tbl[i])) {
					g_unlink(handle->tmp_filename_tbl[i]);
					printf("Reaped file: i = %d, filename = %s\n", i, arg);
					g_free(handle->tmp_filename_tbl[i]);
					handle->tmp_filename_tbl[i] = NULL;
					for (int j = i ; j < handle->ntmp - 1 ; j++) {
						g_free(handle->tmp_filename_tbl[j]);
						handle->tmp_filename_tbl[j] = g_strdup(handle->tmp_filename_tbl[j+1]);
					}
					g_free(handle->tmp_filename_tbl[handle->ntmp - 1]);
					handle->tmp_filename_tbl[handle->ntmp - 1] = NULL;
					handle->ntmp = handle->ntmp - 1;
					break;
				}
			}
			if (handle->ntmp == 0) {
				// Close the gnuplot handle when the temp file
				// counter drops to 0
				g_free(buffer);
				g_object_unref(data_input);
				g_object_unref(stream);
				gnuplot_exit(handle);
				return GINT_TO_POINTER(1);
			}
		} else if (g_str_has_prefix(buffer, "Terminate")) {
			if (handle->ntmp) {
				for (int i = 0 ; i < handle->ntmp ; i++) {
					g_unlink(handle->tmp_filename_tbl[i]);
				}
				handle->ntmp = 0;
			}
			g_free(buffer);
			g_object_unref(data_input);
			g_object_unref(stream);
			gnuplot_exit(handle);
			return GINT_TO_POINTER(1);
		}
		g_free(buffer);
		buffer = NULL;
	}
	printf("exiting tmpwatcher\n");
	g_object_unref(data_input);
	g_object_unref(stream);
	return GINT_TO_POINTER(1);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Opens up a gnuplot session, ready to receive commands.
  @param    keep_plot_alive Flag to keep plots opened after gnuplot process is closed
  @return   Newly allocated gnuplot control structure.

  This opens up a new gnuplot session, ready for input. The struct
  controlling a gnuplot session should remain opaque and only be
  accessed through the provided functions.

  The session must be closed using gnuplot_close().
 */
/*--------------------------------------------------------------------------*/

gnuplot_ctrl * gnuplot_init(gboolean keep_plot_alive)
{
    gnuplot_ctrl *  handle ;
    int i;

    /*
     * Structure initialization:
     */
    handle = (gnuplot_ctrl*)malloc(sizeof(gnuplot_ctrl)) ;
    handle->nplots = 0 ;
    gnuplot_setstyle(handle, "points") ;
    handle->ntmp = 0 ;
	handle->thread = NULL;

    gchar *bin = siril_get_gnuplot_bin();
    gchar* bin2[3];
    bin2[0] = bin;
    bin2[2] = NULL;
    // passing the option --persist keeps the plot opened even after gnuplot process has been closed
	bin2[1] = (keep_plot_alive) ? "--persist" : NULL;
    printf("%s\n", bin2[0]);
    /* call gnuplot */
    gint child_stdin, child_stdout, child_stderr;
    GPid child_pid;
    g_autoptr(GError) error = NULL;

    g_spawn_async_with_pipes(NULL, bin2, NULL,
            G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_SEARCH_PATH,
            NULL, NULL, &child_pid, &child_stdin, &child_stdout,
            &child_stderr, &error);
    if (error != NULL) {
        siril_log_color_message(_("Spawning gnuplot failed: %s\n"), "red", error->message);
        g_free(bin);
        return NULL;
    }

    handle->gnucmd = fdopen(child_stdin, "w");
	handle->gnumon = fdopen(child_stderr, "r");
	handle->child_fd = child_stderr;
	handle->thread = g_thread_new("gplotwatcher", tmpwatcher, handle);

    g_free(bin);
    if (handle->gnucmd == NULL) {
        fprintf(stderr, "error starting gnuplot, is gnuplot or gnuplot.exe in your path?\n") ;
        free(handle);
        return NULL;
    }

    for (i=0;i<GP_MAX_TMP_FILES; i++)
    {
        handle->tmp_filename_tbl[i] = NULL;
    }
    return handle;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   void

  Closes gnuplot by calling an exit command and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.
  This is meant to be called when plot are not displayed

 */
/*--------------------------------------------------------------------------*/

void gnuplot_close(gnuplot_ctrl * handle)
{
	gnuplot_cmd(handle, "print \"Terminate\"");
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Notifies the tmpwatcher that a temporary file can be reaped
  @param    handle Gnuplot session control handle.
  @param    filename Filename to reap.
  @return   void

  Notifed the tmpwatcher thread to reap a temporary file, remove it from the
  list of GNUplot temporary files in the index and decrement the count

 */
/*--------------------------------------------------------------------------*/


void gnuplot_rmtmpfile(gnuplot_ctrl * handle, const char *filename)
{
	gchar *cmd = g_strdup_printf("print \"Done %s\"", filename);
	printf("Calling gnuplot_cmd\n");
	gnuplot_cmd(handle, cmd);
	g_free(cmd);
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   gboolean

  Closes gnuplot by calling an exit command and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.
  This is meant to be called with g_idle_add, when plot are displayed and need to survive

 */
/*--------------------------------------------------------------------------*/

gboolean gnuplot_close_idle(gpointer p) {
    siril_debug_print("closing gnuplot in idle mode\n");
    gnuplot_ctrl *handle = (gnuplot_ctrl *) p;
    gnuplot_close(handle);
    return FALSE;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Sends a command to an active gnuplot session.
  @param    handle Gnuplot session control handle
  @param    cmd    Command to send, same as a printf statement.

  This sends a string to an active gnuplot session, to be executed.
  There is strictly no way to know if the command has been
  successfully executed or not.
  The command syntax is the same as printf.

  Examples:

  @code
  gnuplot_cmd(g, "plot %d*x", 23.0);
  gnuplot_cmd(g, "plot %.18e * cos(%.18e * x)", 32.0, -3.0);
  @endcode

  Since the communication to the gnuplot process is run through
  a standard Unix pipe, it is only unidirectional. This means that
  it is not possible for this interface to query an error status
  back from gnuplot.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_cmd(gnuplot_ctrl *  handle, char const *  cmd, ...)
{
    va_list ap ;

    va_start(ap, cmd);
    vfprintf(handle->gnucmd, cmd, ap);
    va_end(ap);

    fputs("\n", handle->gnucmd) ;
    fflush(handle->gnucmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Change the plotting style of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_style Plotting-style to use (character string)
  @return   void

  The provided plotting style is a character string. It must be one of
  the following:

  - lines
  - points
  - linespoints
  - impulses
  - dots
  - steps
  - errorbars
  - boxes
  - boxeserrorbars
 */
/*--------------------------------------------------------------------------*/

void gnuplot_setstyle(gnuplot_ctrl * h, char * plot_style)
{
    if (strcmp(plot_style, "lines") &&
        strcmp(plot_style, "points") &&
        strcmp(plot_style, "linespoints") &&
        strcmp(plot_style, "impulses") &&
        strcmp(plot_style, "dots") &&
        strcmp(plot_style, "steps") &&
        strcmp(plot_style, "errorbars") &&
        strcmp(plot_style, "boxes") &&
        strcmp(plot_style, "boxerrorbars")) {
        fprintf(stderr, "warning: unknown requested style: using points\n") ;
        strcpy(h->pstyle, "points") ;
    } else {
        strncpy(h->pstyle, plot_style, 31); // 32 char fixed buffer, 1 char for the NULL
    }
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the title of a gnuplot session.
  @author Cyril Richard (03/03/2017)

  @param    h Gnuplot session control handle.
  @param    title Character string to use for title.
  @return   void

  Sets the title for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_title(gnuplot_ctrl * h, char * title)
{
    gnuplot_cmd(h, "set title \"%s\"", title) ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the x label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for X label.
  @return   void

  Sets the x label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_xlabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set xlabel \"%s\"", label) ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Reverse x axis of a gnuplot session.
  @param    h Gnuplot session control handle.
  @return   void

  Reverse x axis of a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_reverse_xaxis(gnuplot_ctrl * h)
{
    gnuplot_cmd(h, "set xrange [] reverse") ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the y label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Y label.
  @return   void

  Sets the y label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_ylabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set ylabel \"%s\"", label) ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Reverse y axis of a gnuplot session.
  @param    h Gnuplot session control handle.
  @return   void

  Reverse y axis of a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_reverse_yaxis(gnuplot_ctrl * h)
{
    gnuplot_cmd(h, "set yrange [] reverse") ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Resets a gnuplot session (next plot will erase previous ones).
  @param    h Gnuplot session control handle.
  @return   void

  Resets a gnuplot session, i.e. the next plot will erase all previous
  ones.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_resetplot(gnuplot_ctrl * h)
{
    if (h->ntmp) {
        for (int i = 0; i < h->ntmp; i++) {
            if (g_remove(h->tmp_filename_tbl[i]) == -1)
                siril_debug_print("g_remove() failed\n");
            free(h->tmp_filename_tbl[i]);
            h->tmp_filename_tbl[i] = NULL;

        }
    }
    h->ntmp = 0 ;
    h->nplots = 0 ;
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Plots a 2d graph from a list of doubles.
  @param    handle  Gnuplot session control handle.
  @param    d       Array of doubles.
  @param    n       Number of values in the passed array.
  @param    title   Title of the plot.
  @return   void

  Plots out a 2d graph from a list of doubles. The x-coordinate is the
  index of the double in the list, the y coordinate is the double in
  the list.

  Example:

  @code
    gnuplot_ctrl    *h ;
    double          d[50] ;
    int             i ;

    h = gnuplot_init(TRUE) ;
    for (i=0 ; i<50 ; i++) {
        d[i] = (double)(i*i) ;
    }
    gnuplot_plot_x(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_x(
    gnuplot_ctrl    *   handle,
    double          *   d,
    int                 n,
    char            *   title
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || d==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = g_fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n ; i++) {
      fprintf(tmpfd, "%.18e\n", d[i]);
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title,0);
	gnuplot_rmtmpfile(handle,tmpfname);
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Plot a 2d graph from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    n           Number of doubles in x (assumed the same as in y).
  @param    title       Title of the plot.
  @return   void

  Plots out a 2d graph from a list of points. Provide points through a list
  of x and a list of y coordinates. Both provided arrays are assumed to
  contain the same number of values.

  @code
    gnuplot_ctrl    *h ;
    double          x[50] ;
    double          y[50] ;
    int             i ;

    h = gnuplot_init(TRUE) ;
    for (i=0 ; i<50 ; i++) {
        x[i] = (double)(i)/10.0 ;
        y[i] = x[i] * x[i] ;
    }
    gnuplot_plot_xy(h, x, y, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_xy(
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    int                 n,
    char            *   title
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = g_fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    // Write Title
    if (title != NULL)
    {
        fprintf(tmpfd, "%s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e\n", x[i], y[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_xy_from_datfile(handle,tmpfname);
	gnuplot_rmtmpfile(handle,tmpfname);
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Plot a 2d graph from a list of points with errors.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    yerr        Pointer to a list of y coordinate errors.
  @param    n           Number of doubles in x (assumed the same as in y).
  @param    title       Title of the plot.
  @return   void

  Plots out a 2d graph from a list of points. Provide points through a list
  of x and a list of y coordinates. Both provided arrays are assumed to
  contain the same number of values.

 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_xyyerr(
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    double          *   yerr,
    int                 n,
    char            *   title,
    int             x_offset /* the entire part of julian date, useful for plotting */
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || yerr==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = g_fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%14.6f %8.6f %8.6f\n", x[i], y[i], yerr[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title, x_offset);
	gnuplot_rmtmpfile(handle,tmpfname);
    return ;
}


void gnuplot_plot_xyyerr_from_datfile(
    gnuplot_ctrl * handle,
    const char   * datfile,
    char         * title,
    int            x_offset
)
{
    gnuplot_plot_atmpfile(handle, datfile, title, x_offset);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Open a new session, plot a signal, close the session.
  @param    title   Plot title
  @param    style   Plot style
  @param    label_x Label for X
  @param    label_y Label for Y
  @param    x       Array of X coordinates
  @param    y       Array of Y coordinates (can be NULL)
  @param    n       Number of values in x and y.
  @return

  This function opens a new gnuplot session, plots the provided
  signal as an X or XY signal depending on a provided y, waits for
  a carriage return on stdin and closes the session.

  It is Ok to provide an empty title, empty style, or empty labels for
  X and Y. Defaults are provided in this case.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_once(
  char    *   title,
  char    *   style,
  char    *   label_x,
  char    *   label_y,
  double  *   x,
  double  *   y,
  int         n
)
{
  gnuplot_ctrl    *   handle ;

  if (x==NULL || n<1) return ;

  if ((handle = gnuplot_init(TRUE)) == NULL) return ;
  if (style!=NULL) {
      gnuplot_setstyle(handle, style);
  } else {
      gnuplot_setstyle(handle, "lines");
  }
  if (label_x!=NULL) {
      gnuplot_set_xlabel(handle, label_x);
  } else {
      gnuplot_set_xlabel(handle, "X");
  }
  if (label_y!=NULL) {
      gnuplot_set_ylabel(handle, label_y);
  } else {
      gnuplot_set_ylabel(handle, "Y");
  }
  if (y==NULL) {
      gnuplot_plot_x(handle, x, n, title);
  } else {
      gnuplot_plot_xy(handle, x, y, n, title);
  }
  printf("press ENTER to continue\n");
  while (getchar()!='\n') {}
  gnuplot_close(handle);
  return ;
}

void gnuplot_plot_slope(
    gnuplot_ctrl    *   handle,
    double              a,
    double              b,
    char            *   title
)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    gnuplot_cmd(handle, "%s %.18e * x + %.18e title \"%s\" with %s",
                  cmd, a, b, title, handle->pstyle) ;

    handle->nplots++ ;
    return ;
}


void gnuplot_plot_equation(
    gnuplot_ctrl    *   h,
    char            *   equation,
    char            *   title
)
{
    char const *    cmd    = (h->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    gnuplot_cmd(h, "%s %s title \"%s\" with %s",
                  cmd, equation, title, h->pstyle) ;
    h->nplots++ ;
    return ;
}


int gnuplot_write_x_csv(
    char const * fileName,
    double const * d,
    int n,
    char const * title)
{
    int     i;
    FILE*   fileHandle;

    if (fileName==NULL || d==NULL || (n<1))
    {
        return -1;
    }

    fileHandle = g_fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%d, %.18e\n", i, d[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_xy_csv(
    char const *        fileName,
    double const    *   x,
    double const    *   y,
    int                 n,
    char const      *   title)
{
    int     i ;
    FILE*   fileHandle;

    if (fileName==NULL || x==NULL || y==NULL || (n<1))
    {
        return -1;
    }

    fileHandle = g_fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%.18e, %.18e\n", x[i], y[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_xy_dat(
    char const *        fileName,
    double const    *   x,
    double const    *   y,
    int                 n,
    char const      *   title)
{
    int     i ;
    FILE*   fileHandle;

    if (fileName==NULL || x==NULL || y==NULL || (n<1))
    {
        return -1;
    }

    fileHandle = g_fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Title
    if (title != NULL)
    {
        fprintf(fileHandle, "%s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%8.6f %8.6f\n", x[i], y[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_xrgb_dat(
    char const *        fileName,
    double const    *   x,
    double const    *   r,
    double const    *   g,
    double const    *   b,
    int                 n,
    char const      *   title)
{
    int     i ;
    FILE*   fileHandle;

    if (fileName==NULL || x==NULL || r==NULL || g == NULL || b == NULL || (n<1))
    {
        return -1;
    }

    fileHandle = g_fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "%s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%8.6f %8.6f %8.6f %8.6f\n", x[i], r[i], g[i], b[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_xyyerr_dat(char const *fileName, double const *x,
        double const *y, double const *yerr, int n, char const *title) {
    if (!fileName || !x || !y || !yerr || n < 1) {
        return -1;
    }

    FILE *fileHandle = g_fopen(fileName, "w");
    if (!fileHandle) {
        perror("creating data file");
        return -1;
    }

    // Write Comment.
    if (title)
        fprintf(fileHandle, "# %s\n", title);

    /* Write data to this file  */
    for (int i=0; i<n; i++) {
        fprintf(fileHandle, "%14.6f %8.6f %8.6f\n", x[i], y[i], yerr[i]);
    }

    fclose(fileHandle);
    return 0;
}

int gnuplot_write_multi_csv(
    char const *        fileName,
    double const    **  xListPtr,
    int                 n,
    int                 numColumns,
    char const      *   title)
{
    int     i;
    int     j;
    FILE*   fileHandle;

    if (fileName==NULL || xListPtr==NULL || (n<1) || numColumns <1)
    {
        return -1;
    }

    for (j=0;j<numColumns;j++)
    {
        if (xListPtr[j] == NULL)
        {
            return -1;
        }
    }

    fileHandle = g_fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%d, %.18e", i, xListPtr[0][i]) ;
        for (j=1;j<numColumns;j++)
        {
            fprintf(fileHandle, ", %.18e", xListPtr[j][i]) ;
        }
        fprintf(fileHandle, "\n");
    }

    fclose(fileHandle) ;

    return 0;
}

char const * gnuplot_tmpfile(gnuplot_ctrl * handle)
{
    static char const * tmp_filename_template = "gnuplot_tmpdatafile_XXXXXX";
    char *              tmp_filename = NULL;
    char const        * tmp_dir = g_get_tmp_dir();
#ifndef _WIN32
    int                 unx_fd;
#endif // #ifndef _WIN32

    assert(handle->tmp_filename_tbl[handle->ntmp] == NULL);

    /* Open one more temporary file? */
    if (handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return NULL;
    }

/* Due to a Windows behavior and Mingw temp file name,
 * we escapes the special characters by inserting a '\' before them */
#ifdef _WIN32
    gchar *tmp = g_build_filename(tmp_dir, tmp_filename_template, NULL);
    tmp_filename = g_strescape(tmp, NULL);
    g_free(tmp);
#else
    tmp_filename = g_build_filename(tmp_dir, tmp_filename_template, NULL);
#endif

#ifdef _WIN32
    if (_mktemp(tmp_filename) == NULL)
    {
        return NULL;
    }
#else // #ifdef _WIN32
    unx_fd = mkstemp(tmp_filename);
    if (unx_fd == -1)
    {
        free(tmp_filename);
        return NULL;
    }
    close(unx_fd);

#endif // #ifdef _WIN32

    handle->tmp_filename_tbl[handle->ntmp] = tmp_filename;
    handle->ntmp ++;
    return tmp_filename;
}

void gnuplot_plot_xy_from_datfile(gnuplot_ctrl * handle, char const* tmp_filename)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    gnuplot_cmd(handle, "%s \"%s\" using ($1):($2) with %s title columnheader",
		   cmd, tmp_filename, handle->pstyle);
    handle->nplots++ ;
    return ;
}

void gnuplot_plot_xrgb_from_datfile(gnuplot_ctrl * handle, char const* tmp_filename)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    gnuplot_cmd(handle, "%s for [col=2:4] \"%s\" using ($1):col with %s title columnheader",
		   cmd, tmp_filename, handle->pstyle);
    handle->nplots++ ;
    return ;
}
void gnuplot_plot_xy_datfile_to_png(gnuplot_ctrl * handle, char const* dat_filename,
		char const *curve_title, char const* png_filename)
{
    gnuplot_cmd(handle, "set term png size 800,600");
    gnuplot_cmd(handle, "set output \"%s\"", png_filename);

    if (curve_title && curve_title[0] != '\0')
	    gnuplot_cmd(handle, "plot \"%s\" using ($1):($2) title \"%s\" with %s", dat_filename,
			    curve_title, handle->pstyle);
    else
	    gnuplot_cmd(handle, "plot \"%s\" with %s", dat_filename,
			    handle->pstyle);
}

void gnuplot_plot_xrgb_datfile_to_png(gnuplot_ctrl * handle, char const* dat_filename,
		char const *curve_title, char const* png_filename)
{
    gnuplot_cmd(handle, "set term png size 800,600");
    gnuplot_cmd(handle, "set output \"%s\"", png_filename);

    if (curve_title && curve_title[0] != '\0')
	    gnuplot_cmd(handle, "plot \"%s\" using ($1):($2):($3):($4) title \"%s\" with %s", dat_filename,
			    curve_title, handle->pstyle);
    else
	    gnuplot_cmd(handle, "plot \"%s\" with %s", dat_filename,
			    handle->pstyle);
}

void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename, char const* title, int x_offset)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;
    gnuplot_cmd(handle, "%s \"%s\" using ($1 - %d):($2):($3) title \"%s\" with %s",
           cmd, tmp_filename, x_offset, title, handle->pstyle);
    handle->nplots++ ;
    return ;
}

void gnuplot_plot_datfile_to_png(gnuplot_ctrl * handle, char const* dat_filename,
        char const *curve_title, int offset, char const* png_filename)
{
    gnuplot_cmd(handle, "set term png size 800,600");
    gnuplot_cmd(handle, "set output \"%s\"", png_filename);

    if (curve_title && curve_title[0] != '\0')
        gnuplot_cmd(handle, "plot \"%s\" using ($1 - %d):($2):($3) title \"%s\" with %s", dat_filename,
                offset, curve_title, handle->pstyle);
    else
        gnuplot_cmd(handle, "plot \"%s\" with %s", dat_filename,
                handle->pstyle);
}

void gnuplot_multiplot_3xy(gnuplot_ctrl * handle, double *x, double *y1, double *y2, double *y3, int n)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y1==NULL || y2 == NULL || y3 == NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = g_fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e %.18e %.18e\n", x[i], y1[i], y2[i], y3[i]) ;
    }
    fclose(tmpfd) ;

	char *curve_title = strdup("Title");
	gnuplot_cmd(handle, "set multiplot layout 3,1 rowsfirst");
	gnuplot_cmd(handle, "plot \"%s\" using ($1):($2) title \"%s\" with %s", tmpfname,
		curve_title, handle->pstyle);
	gnuplot_cmd(handle, "plot \"%s\" using ($1):($3) title \"%s\" with %s", tmpfname,
		curve_title, handle->pstyle);
	gnuplot_cmd(handle, "plot \"%s\" using ($1):($4) title \"%s\" with %s", tmpfname,
		curve_title, handle->pstyle);
	gnuplot_cmd(handle, "unset multiplot");
	free(curve_title);
}

/* vim: set ts=4 et sw=4 tw=75 */
