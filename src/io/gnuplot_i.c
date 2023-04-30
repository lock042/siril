

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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>	// Required for definition of NAN

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

#include "gnuplot_i.h"
#include "gui/plot.h"
#include "core/siril_log.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"

#ifdef _WIN32
#define GNUPLOT_BIN "gnuplot.exe"
#else
#define GNUPLOT_BIN "gnuplot"
#endif

// Uncomment the following line for lots of debug messages
#define GPLOT_DEBUG

static gboolean gnuplot_is_in_path = FALSE;

/*********************** finding gnuplot first **********************/
static gchar *siril_get_gnuplot_bin() {
    if (gnuplot_is_in_path)
        return g_strdup(GNUPLOT_BIN);
    return g_build_filename(com.pref.gnuplot_dir, GNUPLOT_BIN, NULL);
}

gchar* gnuplot_version_is_bad() {
	gchar* retval = NULL;
	gchar* bin2[3];
	gchar* child_stdout = NULL;
	bin2[0] = siril_get_gnuplot_bin();
	bin2[1] = "--version";
	bin2[2] = NULL;
	g_autoptr(GError) error = NULL;
	g_spawn_sync(NULL, bin2, NULL, G_SPAWN_SEARCH_PATH | G_SPAWN_STDERR_TO_DEV_NULL, NULL, NULL, &child_stdout, NULL, NULL, &error);
	if (error) {
		retval = g_strdup(_("Error: failed to execute GNUplot to check version\n"));
	} else {
		if (g_strstr_len(child_stdout, -1, "gnuplot")) {
			g_strchomp(child_stdout);
			g_strchug(child_stdout);
			gchar** chunks = g_strsplit(child_stdout, " ", 5);
			if (g_strv_length(chunks) != 4) {
				retval = g_strdup_printf(_("Could not determine version from version string %s\n"), child_stdout);
			}
			double ver = g_ascii_strtod(chunks[1], NULL);
			double rev = g_ascii_strtod(chunks[3], NULL);
			g_strfreev(chunks);
			if (ver < 5.4 || (ver == 5.4 && rev < 6)) {
#ifdef _WIN32
				retval = g_strdup(_("Error: Windows requires GNUplot >= 5.4.6 which fixes a critical bug that prevents its use with Siril. Please update your GNUplot installation"));
#endif
				siril_debug_print("Detected GNUplot version that would cause an error on Windows\n");
			}
		}
	}
	g_free(child_stdout);
	return retval;
}

#if defined (_WIN32) || defined(OS_OSX)
gboolean gnuplot_is_available() {
	gchar *bin = siril_get_gnuplot_bin();
    if (!bin) return FALSE;

    gboolean is_available = g_file_test(bin, G_FILE_TEST_EXISTS);
    g_free(bin);
	if (is_available) {
		gchar *msg = gnuplot_version_is_bad();
		if (msg) {
			if (!com.script) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Bad GNUplot version"), msg);
				control_window_switch_to_tab(OUTPUT_LOGS);
			} else {
				siril_log_color_message("%s\n", "red", msg);
			}
			is_available = FALSE;
			g_free(msg);
		}
	}
    return is_available;
}

#else
/* returns true if the command gnuplot is available */
gboolean gnuplot_is_available() {
	gboolean is_available;
    gchar *str = g_strdup_printf("%s -e > /dev/null 2>&1", GNUPLOT_BIN);

    int retval = system(str);
    g_free(str);
    if (WIFEXITED(retval)) {
        gnuplot_is_in_path = TRUE;
        is_available = (0 == WEXITSTATUS(retval));
    } else {
		gchar *bin = siril_get_gnuplot_bin();
		is_available = g_file_test(bin, G_FILE_TEST_EXISTS);
		g_free(bin);
	}
	if (is_available) {
		gchar *msg = gnuplot_version_is_bad();
		if (msg) {
			if (!com.script) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Bad GNUplot version"), msg);
				control_window_switch_to_tab(OUTPUT_LOGS);
			} else {
				siril_log_color_message("%s\n", "red", msg);
			}
			is_available = FALSE;
			g_free(msg);
		}
	}
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
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    gnuplot tmpfile watcher. Monitors for tmp files that are finished
			with and closes them when required.
  @param    gpointer user_data. Pointer to gnuplot_ctrl handle.
  @return   GINT_TO_POINTER(1)

*/
/*--------------------------------------------------------------------------*/

gpointer tmpwatcher (gpointer user_data) {
#ifdef GPLOT_DEBUG
	siril_debug_print("tmpwatcher started\n");
#endif
	gnuplot_ctrl* handle = (gnuplot_ctrl*) user_data;
	GInputStream *stream = NULL;
	g_autoptr(GError) error = NULL;
	g_autoptr(GError) error2 = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(handle->child_fd_stderr), FALSE);
#else
	stream = g_unix_input_stream_new(handle->child_fd_stderr, FALSE);
#endif
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL))) {
#ifdef GPLOT_DEBUG
		siril_debug_print("No. of tmp files: %d\n", handle->ntmp);
		siril_debug_print("Buffer: %s\n", buffer);
#endif
		gchar *arg = buffer;
		if (g_str_has_prefix(buffer, "Reap ")) {
#ifdef GPLOT_DEBUG
			siril_debug_print("Received Reap message ntmp = %d\n", handle->ntmp);
#endif
			if (!handle->ntmp)
				continue;
			arg += 5;
			for (int i = 0 ; i < handle->ntmp ; i++) {
#ifdef GPLOT_DEBUG
				siril_debug_print("%s / %s\n", arg, handle->tmp_filename_tbl[i]);
#endif
				if (!g_strcmp0(arg, handle->tmp_filename_tbl[i])) {
					if (g_unlink(handle->tmp_filename_tbl[i]))
						siril_debug_print("Error in g_unlink()\n");
#ifdef GPLOT_DEBUG
					siril_debug_print("Reaped file: i = %d, filename = %s\n", i, arg);
#endif
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
		}
		else if (g_str_has_prefix(buffer, "Terminate")) {
			gnuplot_cmd(handle, "exit gnuplot\n");
			if (handle->ntmp) {
				for (int i = 0 ; i < handle->ntmp ; i++) {
					if (g_unlink(handle->tmp_filename_tbl[i]))
						siril_debug_print("Error in g_unlink()\n");
					free(handle->tmp_filename_tbl[i]);
					handle->tmp_filename_tbl[i] = NULL;
				}
			}
			free(handle->tmp_filename_tbl);
			handle->tmp_filename_tbl = NULL;
			handle->ntmp = 0;
			g_free(buffer);
			g_object_unref(data_input);
			g_object_unref(stream);

			if (!g_close(handle->child_fd_stdin, &error))
				siril_debug_print("%s\n", error->message);
			if (!g_close(handle->child_fd_stderr, &error2))
				siril_debug_print("%s\n", error->message);
			handle->running = FALSE; // Don't free the handle here, it will be freed in gnuplot_close()
			return GINT_TO_POINTER(1);
		}
		g_free(buffer);
		buffer = NULL;
	}
	g_object_unref(data_input);
	g_object_unref(stream);
	return GINT_TO_POINTER(1);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Opens up a gnuplot session, ready to receive commands.
  @param    None
  @return   Newly allocated gnuplot control structure.

  This opens up a new gnuplot session, ready for input. The struct
  controlling a gnuplot session should remain opaque and only be
  accessed through the provided functions.

  The session must be closed using gnuplot_close().
 */
/*--------------------------------------------------------------------------*/

static void child_watch_cb(GPid pid, gint status, gpointer user_data) {
	// This handles cleanup if the GNUplot program dies
	// e.g. if the user closes it using "s" or closing the window
	gnuplot_ctrl* handle = (gnuplot_ctrl*) user_data;
	if (!handle) {
		g_spawn_close_pid(pid);
		return;
	}
#ifdef GPLOT_DEBUG
	siril_debug_print("Closing handle %lu via callback\n", (size_t) handle->thread);
#endif
	if (handle->ntmp) {
		for (int i = 0 ; i < handle->ntmp ; i++) {
			if (g_unlink(handle->tmp_filename_tbl[i]))
				siril_debug_print("Error in g_unlink()\n");
			free(handle->tmp_filename_tbl[i]);
			handle->tmp_filename_tbl[i] = NULL;
		}
	}
	free(handle->tmp_filename_tbl);
	handle->tmp_filename_tbl = NULL;
	handle->ntmp = 0;
	handle->running = FALSE;
	g_autoptr(GError) error = NULL;
	g_autoptr(GError) error2 = NULL;
	if (!g_close(handle->child_fd_stdin, &error))
		siril_debug_print("%s\n", error->message);
	if (!g_close(handle->child_fd_stderr, &error2))
		siril_debug_print("%s\n", error2->message);
	null_handle_in_com_gnuplot_handles(handle);
	free(handle);
	handle = NULL;
	g_spawn_close_pid(pid);
	return;
}

gnuplot_ctrl * gnuplot_init()
{
    gnuplot_ctrl *  handle ;

    /*
     * Structure initialization:
     */
    handle = (gnuplot_ctrl*)malloc(sizeof(gnuplot_ctrl)) ;
	handle->tmp_filename_tbl = calloc(1, sizeof(char*));
	handle->tmp_filename_tbl[0] = NULL;
	handle->ntmp = 0;
	handle->nplots = 0;
	handle->replot = FALSE;
    gnuplot_setstyle(handle, "points") ;
    handle->ntmp = 0 ;
	handle->thread = NULL;

    gchar *bin = siril_get_gnuplot_bin();
    gchar* bin2[3];
    bin2[0] = bin;
    bin2[1] = "--slow";
	bin2[2] = NULL;
#ifdef GPLOT_DEBUG
	siril_debug_print("GNUplot executable: %s\n", bin2[0]);
#endif
	/* call gnuplot */
    gint child_stdin, child_stderr;
    GPid child_pid;
    g_autoptr(GError) error = NULL;

    g_spawn_async_with_pipes(NULL, bin2, NULL,
            G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD,
            NULL, NULL, &child_pid, &child_stdin, NULL,
            &child_stderr, &error);
    if (error != NULL) {
        siril_log_color_message(_("Spawning gnuplot failed: %s\n"), "red", error->message);
        g_free(bin);
		free(handle);
        return NULL;
    }
	g_child_watch_add(child_pid, child_watch_cb, handle);
	handle->running = TRUE;
    handle->gnucmd = fdopen(child_stdin, "w");
	handle->child_fd_stderr = child_stderr;
    handle->child_fd_stdin = child_stdin;
    handle->child_pid = child_pid;
	handle->thread = g_thread_new("gplotwatcher", tmpwatcher, handle);

    g_free(bin);
    if (handle->gnucmd == NULL) {
        fprintf(stderr, "error starting gnuplot, is gnuplot or gnuplot.exe in your path?\n") ;
        free(handle);
        return NULL;
    }
    // Add handle to the list of gnuplot handles
    com.gnuplot_handles = realloc(com.gnuplot_handles, (com.num_gnuplot_handles + 1) * sizeof(gnuplot_ctrl*));
	com.gnuplot_handles[com.num_gnuplot_handles] = handle;
	com.num_gnuplot_handles++;

	gchar *cmd = g_strdup("bind \"Close\" \"unset output ; exit gnuplot\"\n");
	gnuplot_cmd(handle, cmd);
	g_free(cmd);

    return handle;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Exits all GNUplot handles in com.gnuplot_handles
  @param    NULL
  @return   void
 */
/*--------------------------------------------------------------------------*/

void exit_com_gnuplot_handles() {
	for (int i = com.num_gnuplot_handles - 1; i >= 0 ; i--) {
		if (com.gnuplot_handles[i]) {
			gnuplot_close(com.gnuplot_handles[i]);
		}
	}
}

void null_handle_in_com_gnuplot_handles(gnuplot_ctrl* handle) {
	if (!com.gnuplot_handles)
		return;
	for (int i = 0; i < com.num_gnuplot_handles ; i++) {
		if (com.gnuplot_handles && com.gnuplot_handles[i] && com.gnuplot_handles[i] == handle) {
			com.gnuplot_handles[i] = NULL;
			for (int j = i ; j < com.num_gnuplot_handles - 1 ; j++) {
				com.gnuplot_handles[i] = com.gnuplot_handles[i+1];
			}
			break;
		}
	}
	com.num_gnuplot_handles--;
	com.gnuplot_handles = realloc(com.gnuplot_handles, com.num_gnuplot_handles * sizeof(gnuplot_ctrl*));
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

	while (TRUE) {
		g_usleep(1000);
	if (!handle->running)
			break;
	}
	null_handle_in_com_gnuplot_handles(handle);
	free(handle);
	handle = NULL;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Declares a filename as a GNUplot temporary file
			This allows controlled reaping using gnuplot_rmtmpfile
  @param    handle Gnuplot session control handle.
  @param    filename Filename to reap.
  @return   void

 */
/*--------------------------------------------------------------------------*/

void gnuplot_declaretmpfile(gnuplot_ctrl *handle, char *filename) {
    assert(handle->tmp_filename_tbl[handle->ntmp] == NULL);
	handle->tmp_filename_tbl = realloc(handle->tmp_filename_tbl, (handle->ntmp + 2) * sizeof(char*));
	handle->tmp_filename_tbl[handle->ntmp] = strdup(filename);
	handle->tmp_filename_tbl[handle->ntmp + 1] = NULL;
	handle->ntmp++;
#ifdef GPLOT_DEBUG
	siril_debug_print("GNUplot tmpfile %s declared, new ntmp %d\n", filename, handle->ntmp);
#endif
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
	gchar *cmd = g_strdup_printf("print \"Reap %s\"", filename);
#ifdef GPLOT_DEBUG
	siril_debug_print("Calling gnuplot_cmd\n");
#endif
	gnuplot_cmd(handle, cmd);
	g_free(cmd);
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

int gnuplot_cmd(gnuplot_ctrl *  handle, char const *  cmd, ...)
{
    int retval;
	va_list ap ;

    va_start(ap, cmd);
    vfprintf(handle->gnucmd, cmd, ap);
    va_end(ap);

    retval = fputs("\n", handle->gnucmd) ;
    fflush(handle->gnucmd) ;
    return retval;
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

void gnuplot_resetplot(gnuplot_ctrl * handle)
{
    if (handle->ntmp) {
        for (int i = 0; i < handle->ntmp; i++) {
            if (g_remove(handle->tmp_filename_tbl[i]) == -1)
#ifdef GPLOT_DEBUG
                siril_debug_print("g_remove() failed\n");
#endif
			free(handle->tmp_filename_tbl[i]);
            handle->tmp_filename_tbl[i] = NULL;
        }
    }
    handle->ntmp = 0 ;
    handle->nplots = 0 ;
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

    h = gnuplot_init() ;
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
	handle->nplots++;
//	gnuplot_rmtmpfile(handle,tmpfname);
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

    h = gnuplot_init() ;
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
	handle->nplots++;
//	gnuplot_rmtmpfile(handle,tmpfname);
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
	handle->nplots++;
//	gnuplot_rmtmpfile(handle,tmpfname);
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
// This function should not be used in Siril as it uses printf and getchar for
// control and indication
/*
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

  if ((handle = gnuplot_init()) == NULL) return ;
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
*/

void gnuplot_plot_slope(
    gnuplot_ctrl    *   handle,
    double              a,
    double              b,
    char            *   title
)
{
    char const *    cmd    = (handle->replot && handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    gnuplot_cmd(handle, "%s %.18e * x + %.18e title \"%s\" with %s",
                  cmd, a, b, title, handle->pstyle) ;

    handle->nplots++ ;
    return ;
}


void gnuplot_plot_equation(
    gnuplot_ctrl    *   handle,
    char            *   equation,
    char            *   title
)
{
    char const *    cmd    = (handle->replot && handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    gnuplot_cmd(handle, "%s %s title \"%s\" with %s",
                  cmd, equation, title, handle->pstyle) ;
    handle->nplots++ ;
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

int gnuplot_write_xcfa_dat(
    char const *        fileName,
    double const    *   x,
    double const    *   cfa0,
    double const    *   cfa1,
    double const    *   cfa2,
    double const    *   cfa3,
    int                 n,
    char const      *   title)
{
    int     i ;
    FILE*   fileHandle;

    if (fileName==NULL || x==NULL || cfa0==NULL || cfa1 == NULL || cfa2 == NULL || cfa3 == NULL || (n<1))
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
        fprintf(fileHandle, "%8.6f %8.6f %8.6f %8.6f %8.6f\n", x[i], cfa0[i], cfa1[i], cfa2[i], cfa3[i]) ;
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
    char const *    cmd    = (handle->replot && handle->nplots > 0) ? "replot" : "plot";
//    gnuplot_cmd(handle, "set term wxt raise persist");
    gnuplot_cmd(handle, "%s \"%s\" using ($1):($2) with %s title columnheader",
		   cmd, tmp_filename, handle->pstyle);
    handle->nplots++ ;
    return ;
}

void gnuplot_plot_xrgb_from_datfile(gnuplot_ctrl * handle, char const* tmp_filename)
{
    char const *    cmd    = (handle->replot && handle->nplots > 0) ? "replot" : "plot";
//    gnuplot_cmd(handle, "set term wxt raise persist");
    gnuplot_cmd(handle, "%s for [col=2:4] \"%s\" using ($1):col with %s title columnheader",
		   cmd, tmp_filename, handle->pstyle);
    handle->nplots++ ;
    return ;
}

void gnuplot_plot_xcfa_from_datfile(gnuplot_ctrl * handle, char const* tmp_filename)
{
    char const *    cmd    = (handle->replot && handle->nplots > 0) ? "replot" : "plot";
    gnuplot_cmd(handle, "%s for [col=2:5] \"%s\" using ($1):col with %s title columnheader",
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
	    gnuplot_cmd(handle, "plot \"%s\" using ($1):($2) with %s title \"%s\"", dat_filename,
			    handle->pstyle, curve_title);
    else
	    gnuplot_cmd(handle, "plot \"%s\" with %s", dat_filename,
			    handle->pstyle);
}

void gnuplot_plot_xy_datfile_colheader_to_png(gnuplot_ctrl * handle, char const* dat_filename,
		char const *curve_title, char const* png_filename)
{
    gnuplot_cmd(handle, "set term png size 800,600");
    gnuplot_cmd(handle, "set output \"%s\"", png_filename);

    if (curve_title && curve_title[0] != '\0')
	    gnuplot_cmd(handle, "plot \"%s\" using ($1):($2) with %s title columnheader", dat_filename,
			    handle->pstyle);
    else
	    gnuplot_cmd(handle, "plot \"%s\" with %s", dat_filename,
			    handle->pstyle);
}

void gnuplot_plot_xrgb_datfile_to_png(gnuplot_ctrl * handle, char const* dat_filename,
		char const* png_filename)
{
    gnuplot_cmd(handle, "set term png size 800,600");
    gnuplot_cmd(handle, "set output \"%s\"", png_filename);

	gnuplot_cmd(handle, "plot for [col=2:4] \"%s\" using ($1):col with %s title columnheader",
				dat_filename, handle->pstyle);
}

void gnuplot_plot_xcfa_datfile_to_png(gnuplot_ctrl * handle, char const* dat_filename,
		char const* png_filename)
{
    gnuplot_cmd(handle, "set term png size 800,600");
    gnuplot_cmd(handle, "set output \"%s\"", png_filename);

	gnuplot_cmd(handle, "plot for [col=2:5] \"%s\" using ($1):col with %s title columnheader",
				dat_filename, handle->pstyle);
}

void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename, char const* title, int x_offset)
{
    char const *    cmd    = (handle->replot && handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;
    gnuplot_cmd(handle, "%s \"%s\" using ($1 - %d):($2):($3) title \"%s\" with %s",
           cmd, tmp_filename, x_offset, title, handle->pstyle);
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
/*
 * Not used in Siril yet
void gnuplot_multiplot_3xy(gnuplot_ctrl * handle, double *x, double *y1, double *y2, double *y3, int n)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y1==NULL || y2 == NULL || y3 == NULL || (n<1)) return ;

    // Open temporary file for output
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = g_fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    // Write data to this file
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
*/
/* vim: set ts=4 et sw=4 tw=75 */
