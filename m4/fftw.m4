# fftw has several ways to be compiled in an optimized way and the way we link
# against it changes with the architecture
# In this file we detect which works and set the correct CFLAGS and LDFLAGS
#
# From the documentation:
# programs using the parallel complex transforms should be linked with
# -lfftw3_threads -lfftw3 -lm on Unix, or -lfftw3_omp -lfftw3 -lm if you
# compiled with OpenMP. You will also need to link with whatever library is
# responsible for threads on your system (e.g. -lpthread on GNU/Linux) or
# include whatever compiler flag enables OpenMP (e.g. -fopenmp with gcc). 

siril_save_CFLAGS="$CFLAGS"
siril_save_LIBS="$LIBS"

dnl we check for the availability of the lib first
PKG_CHECK_MODULES(FFTW, [fftw3f])

dnl then we check for its threading capabilities, first OpenMP
AC_MSG_CHECKING([for fftw (fftw3f openmp threads)])
CFLAGS="$CFLAGS $FFTW_CFLAGS"
LIBS="$LIBS -lfftw3f_omp"
AC_LANG([C])
AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <fftw3.h>],
			       [fftwf_init_threads();])],
              [have_fftw_openmp="yes"],
              [have_fftw_openmp="no"],
              [have_fftw_openmp="unknown (cross-compiling)"])
CFLAGS="$siril_save_CFLAGS"
LIBS="$siril_save_LIBS"
AC_MSG_RESULT($have_fftw_openmp)

if test "x$have_fftw_openmp" = "xyes"; then
	AC_MSG_NOTICE([fftw3 will use OpenMP optimizations])
	AC_DEFINE([HAVE_FFTW3F_OMP], [1], [Defined if FFTW3 threading is done with OpenMP])
	FFTW_LIBS="$FFTW_LIBS -lfftw3f_omp"
else
	dnl then with pthreads
	AC_MSG_CHECKING([for fftw (fftw3f threads arg)])
	CFLAGS="$CFLAGS $FFTW_CFLAGS"
	LIBS="$LIBS -lfftw3f_threads"
	AC_LANG([C])
	AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <fftw3.h>],
				       [fftwf_init_threads();])],
	              [have_fftw_threads="yes"],
	              [have_fftw_threads="no"],
	              [have_fftw_threads="unknown (cross-compiling)"])
	CFLAGS="$siril_save_CFLAGS"
	LIBS="$siril_save_LIBS"
	AC_MSG_RESULT($have_fftw_threads)

	if test "x$have_fftw_threads" = "xyes"; then
		AC_MSG_NOTICE([fftw3 will use pthreads threading])
		AC_DEFINE([HAVE_FFTW3F_THREADS], [1], [Defined if FFTW3 threading is done with pthread])
		FFTW_LIBS="$FFTW_LIBS -lfftw3f_threads"
	else
		dnl then with no linker flags
		AC_MSG_CHECKING([for fftw (fftw3f threads no arg)])
		CFLAGS="$CFLAGS $FFTW_CFLAGS"
		LIBS="$LIBS"
		AC_LANG([C])
		AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <fftw3.h>],
					       [fftwf_init_threads();])],
		              [have_fftw_nothreads="yes"],
		              [have_fftw_nothreads="no"],
		              [have_fftw_nothreads="unknown (cross-compiling)"])
		CFLAGS="$siril_save_CFLAGS"
		LIBS="$siril_save_LIBS"
		AC_MSG_RESULT($have_fftw_nothreads)
		AC_MSG_WARN([libfftw does not seem to support threading])
	fi
fi

