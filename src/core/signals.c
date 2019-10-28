/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <signal.h>
#include <gio/gio.h>

#include "siril.h"
#include "signals.h"

#ifdef HAVE_EXCHNDL
#include <windows.h>
#include <time.h>
#include <exchndl.h>

static LPTOP_LEVEL_EXCEPTION_FILTER g_prevExceptionFilter = NULL;

static LONG WINAPI  siril_sigfatal_handler (PEXCEPTION_POINTERS pExceptionInfo);

#else

static void siril_sigfatal_handler(gint sig_num) G_GNUC_NORETURN;

#endif

void siril_init_signal_handlers(gchar **backtrace_file) {
	time_t t;
	gchar *filename;
	gchar *dir;

#ifdef _WIN32
  /* This has to be the non-roaming directory (i.e., the local
     directory) as backtraces correspond to the binaries on this
     system. */
  dir = g_build_filename (g_get_user_data_dir (),
                          SIRILDIR, SIRIL_USER_VERSION, "CrashLog",
                          NULL);
#else
	dir = g_build_filename("/tmp/", "CrashLog", NULL);
#endif

	time(&t);
	filename = g_strdup_printf("%s-crash-%" G_GUINT64_FORMAT ".txt",
	PACKAGE_NAME, t);
	*backtrace_file = g_build_filename(dir, filename, NULL);
	g_free(filename);
	g_free(dir);

#ifdef _WIN32
  /* Use Dr. Mingw (dumps backtrace on crash) if it is available. Do
   * nothing otherwise on Win32.
   * The user won't get any stack trace from glib anyhow.
   * Without Dr. MinGW, It's better to let Windows inform about the
   * program error, and offer debugging (if the user has installed MSVC
   * or some other compiler that knows how to install itself as a
   * handler for program errors).
   */

#ifdef HAVE_EXCHNDL
  /* Order is very important here. We need to add our signal handler
   * first, then run ExcHndlInit() which will add its own handler, so
   * that ExcHnl's handler runs first since that's in FILO order.
   */
  if (! g_prevExceptionFilter)
    g_prevExceptionFilter = SetUnhandledExceptionFilter (siril_sigfatal_handler);

  ExcHndlInit ();
  ExcHndlSetLogFileNameA (*backtrace_file);

#endif /* HAVE_EXCHNDL */

#else

	/* Handle fatal signals */

	/* these are handled by siril_terminate() */
	signal(SIGHUP, siril_sigfatal_handler);
	signal(SIGINT, siril_sigfatal_handler);
	signal(SIGQUIT, siril_sigfatal_handler);
	signal(SIGTERM, siril_sigfatal_handler);

	/* these are handled by siril_fatal_error() */
	signal(SIGABRT, siril_sigfatal_handler);
	signal(SIGBUS, siril_sigfatal_handler);
	signal(SIGSEGV, siril_sigfatal_handler);
	signal(SIGFPE, siril_sigfatal_handler);

	/* Ignore SIGPIPE because plug_in.c handles broken pipes */
	signal(SIGPIPE, SIG_IGN);

#endif /* _WIN32 */
}

#ifdef _WIN32

#ifdef HAVE_EXCHNDL
static LONG WINAPI
siril_sigfatal_handler (PEXCEPTION_POINTERS pExceptionInfo)
{
  EXCEPTION_RECORD *er;
  int               fatal;

  if (pExceptionInfo == NULL ||
      pExceptionInfo->ExceptionRecord == NULL)
    return EXCEPTION_CONTINUE_SEARCH;

  er = pExceptionInfo->ExceptionRecord;
  fatal = I_RpcExceptionFilter (er->ExceptionCode);

  /* IREF() returns EXCEPTION_CONTINUE_SEARCH for fatal exceptions */
  if (fatal == EXCEPTION_CONTINUE_SEARCH)
    {
      /* Just in case, so that we don't loop or anything similar, just
       * re-establish previous handler.
       */
      SetUnhandledExceptionFilter (g_prevExceptionFilter);

      /* Now process the exception. */
      siril_fatal_error ("unhandled exception");
    }

  if (g_prevExceptionFilter && g_prevExceptionFilter != siril_sigfatal_handler)
    return g_prevExceptionFilter (pExceptionInfo);
  else
    return EXCEPTION_CONTINUE_SEARCH;
}
#endif

#else

/* siril core signal handler for fatal signals */

static void siril_sigfatal_handler(gint sig_num) {
	switch (sig_num) {
	case SIGHUP:
	case SIGINT:
	case SIGQUIT:
	case SIGTERM:
		siril_terminate(g_strsignal(sig_num));
		break;

	case SIGABRT:
	case SIGBUS:
	case SIGSEGV:
	case SIGFPE:
	default:
		siril_fatal_error(g_strsignal(sig_num));
		break;
	}
}
