dnl Process this file with autoconf to produce a configure script.
m4_define([siril_major_version], [1])
m4_define([siril_minor_version], [4])
m4_define([siril_micro_version], [0])
m4_define([siril_prerelease_version], [beta3])
m4_define([siril_version],
          [siril_major_version.siril_minor_version.siril_micro_version])
m4_define([siril_version_full],
          [siril_major_version().siril_minor_version().siril_micro_version()m4_bpatsubst(siril_prerelease_version(), [^], [-])])

m4_define([siril_unstable],
          m4_if(m4_eval(siril_minor_version % 2), [1], [yes], [no]))
m4_define([siril_stable],
          m4_if(m4_eval(siril_minor_version % 2), [0], [yes], [no]))

dnl change siril_version_full by siril_version for stable release
AC_INIT([siril], [siril_version_full],
dnl AC_INIT([siril], [siril_version],
        [https://gitlab.com/free-astro/siril/issues],
        [siril],
        [https://www.siril.org])

AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_MACRO_DIR([m4])

SIRIL_UNSTABLE=siril_unstable
SIRIL_MAJOR_VERSION=siril_major_version
SIRIL_MINOR_VERSION=siril_minor_version
SIRIL_MICRO_VERSION=siril_micro_version

dnl ================================================================
dnl Gettext stuff.
dnl ================================================================
IT_PROG_INTLTOOL([0.50])
AC_SUBST([GETTEXT_PACKAGE], [siril])
AC_SUBST(SIRIL_UNSTABLE)
AC_SUBST(SIRIL_MAJOR_VERSION)
AC_SUBST(SIRIL_MINOR_VERSION)
AC_SUBST(SIRIL_MICRO_VERSION)

dnl Define these values in config.h
AC_DEFINE([SIRIL_MAJOR_VERSION], [siril_major_version], [Major version of Siril])
AC_DEFINE([SIRIL_MINOR_VERSION], [siril_minor_version], [Minor version of Siril])
AC_DEFINE([SIRIL_MICRO_VERSION], [siril_micro_version], [Micro version of Siril])

remember_set_CFLAGS="$CFLAGS"

AC_PROG_CC

m4_include(m4/openmp.m4)
AC_OPENMP

AC_LANG([C])
AC_LANG([C++])

# remove -g flag from default CFLAGS
if test "x$remember_set_CFLAGS" = "x"; then
	if test "$CFLAGS" = "-g -O2"; then
		CFLAGS="-O2"
	elif test "$CFLAGS" = "-g"; then
		CFLAGS=""
	fi
fi

CFLAGS="-D_FILE_OFFSET_BITS=64 -Werror=unknown-pragmas -Werror=unused-label -fno-common -ftree-vectorize -fno-math-errno $CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CFLAGS -std=c++17"

if test "x$SIRIL_UNSTABLE" = "xyes"; then
  AC_DEFINE(SIRIL_UNSTABLE, 1,
            [Define to 1 if this is an unstable version of Siril])
fi
AM_CONDITIONAL(SIRIL_UNSTABLE, test "x$SIRIL_UNSTABLE" = "xyes")

dnl Set PACKAGE_SOURCE_DIR in config.h.
packagesrcdir=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED(PACKAGE_SOURCE_DIR, "${packagesrcdir}", [Package source dir])

# Checks for header files.
AC_CHECK_HEADERS([execinfo.h sys/vfs.h sys/vmount.h sys/statfs.h sys/statvfs.h sys/filio.h])

# Checks for reallocarray function (used in kplot)
AC_CHECK_FUNCS(reallocarray, AC_DEFINE([HAVE_REALLOCARRAY], [1], [reallocarray is available]), )
# Checks for library functions.
AC_CHECK_FUNCS(timegm gmtime_r)

AC_CHECK_FUNCS(backtrace, , AC_CHECK_LIB(execinfo, backtrace))

dnl Set PACKAGE_DOC_DIR in config.h.
dnl AC_DEFINE_UNQUOTED(PACKAGE_DOC_DIR, "${ac_default_prefix}/share/doc/${PACKAGE}", "package doc dir")

dnl On some platforms, std::atomic needs a helper library
AC_MSG_CHECKING(whether -latomic is needed)
AC_LINK_IFELSE([AC_LANG_SOURCE([[
  #include <atomic>
  #include <cstdint>
  std::atomic<std::int64_t> v;
  int main() {
    return v;
  }
]])], STD_ATOMIC_NEED_LIBATOMIC=no, STD_ATOMIC_NEED_LIBATOMIC=yes)
AC_MSG_RESULT($STD_ATOMIC_NEED_LIBATOMIC)
if test "x$STD_ATOMIC_NEED_LIBATOMIC" = xyes; then
  LIBATOMIC_LIBS="-latomic"
fi
AC_SUBST([LIBATOMIC_LIBS])

dnl ================ PLATFORM CHECKS ================

dnl Compute the canonical host-system type variable (need config.guess and config.sub)
AC_CANONICAL_HOST

AC_MSG_CHECKING([if compiling for Win32])
case "$host_os" in
  mingw* | cygwin*)
    platform_win32=yes
    ;;
  *)
    platform_win32=no
    ;;
esac
AC_MSG_RESULT([$platform_win32])

dnl ================ MANDATORY DEPENDENCIES ================
dnl check math lib
AC_CHECK_LIB(m, log, LIBS="$LIBS -lm")

dnl check GLIB version
PKG_CHECK_MODULES(GLIB, glib-2.0 >= 2.56.0)
case "$host_os" in
  mingw*)
    PKG_CHECK_MODULES([GIO],[gio-2.0 >= 2.56.0 gio-windows-2.0 >= 2.56.0])
    ;;
  *)
    PKG_CHECK_MODULES([GIO],[gio-2.0 >= 2.56.0 gio-unix-2.0 >= 2.56.0])
    ;;
esac

dnl check GTK
PKG_CHECK_MODULES(GTK, gtk+-3.0 >= 3.20.0)

dnl check lcms2
PKG_CHECK_MODULES(LCMS2, lcms2 >= 2.14)

dnl check cairo
PKG_CHECK_MODULES(CAIRO, cairo)

dnl check wcslib
PKG_CHECK_MODULES(WCSLIB, wcslib >= 7.12)

dnl check fftw3 and its threading capabilities
m4_include([m4/fftw.m4])

dnl check GNU Scientific Library
PKG_CHECK_MODULES(GSL, [gsl < 2],
	[AC_DEFINE([HAVE_GSL_1], [1], [Use GSL 1])],
	[PKG_CHECK_MODULES(GSL, [gsl >= 2],
		[AC_DEFINE([HAVE_GSL_2], [1], [Use GSL 2])])
	])

dnl check cfitsio
PKG_CHECK_MODULES(CFITSIO, [cfitsio])

dnl Check for gtksourceview-4
PKG_CHECK_MODULES(GTKSOURCEVIEW, [gtksourceview-4])

dnl check opencv
AC_PROG_CXX
if test -n "$CXX" && ( test "X$CXX" != "Xno" &&
	( (test "X$CXX" = "Xg++" && `g++ -v >/dev/null 2>&1` ) ||
	(test "X$CXX" != "Xg++"))) ; then
	AC_PROG_CXXCPP

	dnl check for C++14 support
	m4_include([m4/ax_cxx_compile_stdcxx.m4])
	AX_CXX_COMPILE_STDCXX(14, noext, mandatory)

	PKG_CHECK_MODULES([OPENCV], [opencv4 >= 4.4.0],
			  [AC_DEFINE([HAVE_CV44], [1], [Using OpenCV4.4])]
			  [OPENCV_LIBS="-lopencv_core -lopencv_imgproc -lopencv_calib3d -lopencv_stitching"],
			  [PKG_CHECK_MODULES([OPENCV], [opencv4 >= 4.2.0],
					     [AC_MSG_WARN([OpenCV 4.2 found, translation-only registration will be disabled - requires 4.4])]
					     [OPENCV_LIBS="-lopencv_core -lopencv_imgproc -lopencv_calib3d -lopencv_stitching"],
					     AC_MSG_ERROR([opencv not found.]))])
else
	AC_MSG_ERROR([no C++ compiler found])
fi

PKG_CHECK_MODULES(LIBRTPROCESS, [rtprocess], have_rtprocess="yes", have_rtprocess="no")
if test "x$have_rtprocess" = "xno"; then
  AC_MSG_NOTICE([librtprocess not found])
  if test "$cross_compiling" != "yes"; then
    dnl check CMake used for librtprocess
    AC_PATH_PROG(CMAKE, cmake, no)
    if test "x$CMAKE" = xno; then
      AC_MSG_ERROR([CMAKE not found (needed to built internal librtprocess)])
    fi

    CMAKE_GENERATOR=
    if test "$platform_win32" = "yes"; then
      CMAKE_GENERATOR="-G \"MSYS Makefiles\""
    fi
    AC_SUBST([CMAKE_GENERATOR])

    AC_MSG_NOTICE([librtprocess will be built statically])
  else
    # When cross-compiling, you should build rtprocess yourself and
    # install it yourself.
    AC_MSG_ERROR([librtprocess not found])
  fi
fi
AM_CONDITIONAL(HAVE_LIBRTPROCESS, test "x$have_rtprocess" = "xyes")

if test "$cross_compiling" != "yes"; then
	dnl check CMake used for htmesh
	AC_PATH_PROG(CMAKE, cmake, no)
	if test "x$CMAKE" = xno; then
		AC_MSG_ERROR([CMAKE not found (needed to built internal htmesh)])
	fi

	CMAKE_GENERATOR=
	if test "$platform_win32" = "yes"; then
		CMAKE_GENERATOR="-G \"MSYS Makefiles\""
	fi
	AC_SUBST([CMAKE_GENERATOR])

	AC_MSG_NOTICE([htmesh will be built statically])
else
	# When cross-compiling, you should build rtprocess yourself and
	# install it yourself.
	AC_MSG_ERROR([htmesh cannot be built])
fi
AM_CONDITIONAL(HAVE_HTMESH, test "x$have_htmesh" = "xyes")

if test "$cross_compiling" != "yes"; then
	dnl check CMake used for healpix_cxx
	AC_PATH_PROG(CMAKE, cmake, no)
	if test "x$CMAKE" = xno; then
		AC_MSG_ERROR([CMAKE not found (needed to built internal healpix_cxx)])
	fi

	CMAKE_GENERATOR=
	if test "$platform_win32" = "yes"; then
		CMAKE_GENERATOR="-G \"MSYS Makefiles\""
	fi
	AC_SUBST([CMAKE_GENERATOR])

	AC_MSG_NOTICE([healpix_cxx will be built statically])
else
	# When cross-compiling, you should build healpix_cxx yourself and
	# install it yourself.
	AC_MSG_ERROR([healpix_cxx cannot be built])
fi
AM_CONDITIONAL(HAVE_HEALPIX_CXX, test "x$have_healpix_cxx" = "xyes")

PKG_CHECK_MODULES(LIBYYJSON, [yyjson >= 0.10.0], have_yyjson="yes", have_yyjson="no")
if test "x$have_yyjson" = "xno"; then
  AC_MSG_NOTICE([libyyjson not found])
  if test "$cross_compiling" != "yes"; then
    dnl check CMake used for yyjson
    AC_PATH_PROG(CMAKE, cmake, no)
    if test "x$CMAKE" = xno; then
      AC_MSG_ERROR([CMAKE not found (needed to built internal yyjson)])
    fi

    CMAKE_GENERATOR=
    if test "$platform_win32" = "yes"; then
      CMAKE_GENERATOR="-G \"MSYS Makefiles\""
    fi
    AC_SUBST([CMAKE_GENERATOR])

    AC_MSG_NOTICE([yyjson will be built statically])
  else
    # When cross-compiling, you should build rtprocess yourself and
    # install it yourself.
    AC_MSG_ERROR([libyyjson not found])
  fi
fi
AM_CONDITIONAL([HAVE_YYJSON], [test "x$have_yyjson" = "xyes"])

dnl ================ OPTIONAL DEPENDENCIES ================

dnl check raw lib
PKG_CHECK_MODULES(LIBRAW, [libraw >= 0.17],
		  [AC_DEFINE([HAVE_LIBRAW], [1], [Using RAW images])],
		  AC_MSG_WARN([libraw not found. Not using RAW importer.]))

dnl check libexiv2
PKG_CHECK_MODULES(EXIV2, [exiv2 >= 0.25],
		  [AC_DEFINE([HAVE_EXIV2], [1], [Using exif previews])],
		  AC_MSG_WARN([exiv2 not found. Not extracting image previews.]))

dnl check tiff lib
PKG_CHECK_MODULES(LIBTIFF, [libtiff-4],
		  [AC_DEFINE([HAVE_LIBTIFF], [1], [Using TIFF images])],
		  AC_MSG_WARN([libtiff not found. Not using TIFF importer and exporter.]))

dnl check jpeg lib
dnl PKG_CHECK_MODULES(JPEG, [libjpeg])
AC_CHECK_LIB(jpeg, jpeg_mem_src, [],
	     AC_MSG_WARN([libjpeg not found. Not using JPEG importer and exporter.]))

dnl check png lib
PKG_CHECK_MODULES(LIBPNG, [libpng >= 1.6],
		  [AC_DEFINE([HAVE_LIBPNG], [1], [Using PNG images])],
		  AC_MSG_WARN([libpng not found. Not using PNG importer and exporter.]))

dnl check heif lib
PKG_CHECK_MODULES(LIBHEIF, [libheif],
		  [AC_DEFINE([HAVE_LIBHEIF], [1], [Using HEIF images])],
		  AC_MSG_WARN([libheif not found. Not using HEIF importer and exporter.]))

dnl check heif lib
PKG_CHECK_MODULES(LIBXISF, [libxisf],
		  [AC_DEFINE([HAVE_LIBXISF], [1], [Using XISF images >= 0.2.7])],
		  AC_MSG_WARN([libxisf not found. Not using XISF importer and exporter.]))

dnl check ffms2
PKG_CHECK_MODULES(FFMS2, [ffms2],
		  [AC_DEFINE([HAVE_FFMS2], [1], [Using FFMS2])],
		  AC_MSG_WARN([libffms2 not found. Not using AVI importer.]))

dnl checking several ffmpeg libraries
PKG_CHECK_MODULES(FFMPEG, [libavformat libavutil >= 55.20 libavcodec libswscale libswresample],
		  [AC_DEFINE([HAVE_FFMPEG], [1], [film export supported])],
		  AC_MSG_WARN([ffmpeg libav libraries not found or too old. Not using film export formats.]))

dnl checking libcurl
curl_avail=no
m4_include([m4/libcurl.m4])
LIBCURL_CHECK_CONFIG([yes], [7.10.0],
		     [AC_DEFINE([HAVE_LIBCURL], [1], [Using Libcurl])] curl_avail=yes, [])

dnl Check for libgit
PKG_CHECK_MODULES([LIBGIT], [libgit2],
		  [AC_DEFINE([HAVE_LIBGIT2], [1], [Using online scripts])],
		  AC_MSG_WARN([libgit2 was not found. Not using online scripts.]))

dnl ================================================================
dnl Collect some arguments.
dnl ================================================================
AC_ARG_ENABLE(debug,
	      [AS_HELP_STRING([--enable-debug],
			      [enable debug build [default=no]])],,
			      [enable_debug="no"])

AM_CONDITIONAL([SIRIL_DEBUG],[test "x$enable_debug" = "xyes"])

if test "x$enable_debug" = "xyes"; then
	AC_DEFINE([SIRIL_OUTPUT_DEBUG], [1], [Using debug output])
	CFLAGS="-g $CFLAGS"
	CXXFLAGS="$CFLAGS"
fi

#####################################
# Check for Platform and architecture
#####################################
platform_family="other"
PATHSEP=':'
AC_MSG_CHECKING([platform and architecture])
case "$host_os" in
	linux-**)
		platform=unix
		platform_family="linux"
		;;
	freebsd-*)
		platform=unix
		platform_family="bsd"
		;;
	cygwin*)
		platform=windows
		platform_family="windows"
		;;
	mingw*)
		platform=windows
		platform_family="windows"
		PATHSEP=';'
		CFLAGS=" $CFLAGS -D_WIN32_WINNT=0x0600 -DWINVER=0x0600 "
		CXXFLAGS=" $CXXFLAGS -D_WIN32_WINNT=0x0600 -DWINVER=0x0600 "
		case "$host_cpu" in
			x86_64)
				;;
			*)
				WIN32_LARGE_ADDRESS_AWARE='-Wl,--large-address-aware'
				;;
		esac
		;;
	darwin*)
		platform=osx
		platform_family="macos"
		AC_SUBST(EXTRA_LIBS,'-framework AppKit')
		CFLAGS="$CFLAGS -ObjC"
		;;
esac
AC_MSG_RESULT([$platform])
AM_CONDITIONAL([PLATFORM_UNIX], [test "$platform" = "unix"])
AM_CONDITIONAL([PLATFORM_WIN32], [test "$platform" = "windows"])
AM_CONDITIONAL([PLATFORM_OSX], [test "$platform" = "osx"])

AC_DEFINE_UNQUOTED(CPU_ARCH, "$host_cpu", [The CPU family this was built for])
AC_DEFINE_UNQUOTED(SIRIL_BUILD_PLATFORM_FAMILY, "$platform_family", [The OS family this was built for])

AC_SUBST(WIN32_LARGE_ADDRESS_AWARE)
AC_SUBST(PATHSEP)

if test "$platform" = "windows"; then
	AC_CHECK_PROG(ms_librarian, lib.exe, yes, no)
	AC_CHECK_TOOL(WINDRES, windres, :)
else
	WINDRES=":"
	AC_PROG_LN_S
fi

AM_CONDITIONAL(MS_LIB_AVAILABLE, test "x$ms_librarian" = xyes)
AM_CONDITIONAL(HAVE_WINDRES, test "x$WINDRES" != "x:")

AC_SUBST(WINDRES)

#####################################
# Check for gdk windowing style
#####################################
AC_MSG_CHECKING([for native gdk windows style])
gdk_targets=`$PKG_CONFIG --variable=targets gdk-3.0`
gdk_windowing="x11"
for target in $gdk_targets;
do
	case "$target" in
		quartz)
			gdk_windowing=quartz
			;;
		win32)
			gdk_windowing=win32
			;;
	esac
done

AC_MSG_RESULT([$gdk_windowing])
AM_CONDITIONAL(GDK_WINDOWING_QUARTZ, test "$gdk_windowing" = "quartz")
AM_CONDITIONAL(GDK_WINDOWING_WIN32, test "$gdk_windowing" = "win32")
AM_CONDITIONAL(GDK_WINDOWING_X11, test "$gdk_windowing" = "x11")

if test "$platform" = "osx"; then
	AC_DEFINE([OS_OSX],[1],[Defined if os is Mac OSX])

	PKG_CHECK_MODULES(GTK_MAC, gtk-mac-integration-gtk3 >= 2.0.8)

	AC_SUBST(GTK_MAC_CFLAGS)
	AC_SUBST(GTK_MAC_LIBS)
fi

if test `uname` != "Darwin"; then
	if (test "$platform" = "windows") ; then
		AC_SUBST(EXTRA_LIBS,'-lpsapi')
	fi
fi

#################
# Bundle packages
#################

AC_MSG_CHECKING([whether we build a relocatable package])
AC_ARG_ENABLE(relocatable-bundle,
              [  --enable-relocatable-bundle
                          build with resources considered bundled under the same prefix
                          (default=auto)],,
              enable_relocatable_bundle=auto)

if test "x$enable_relocatable_bundle" != xno &&
   test "x$enable_relocatable_bundle" != xyes; then
  # By default, assume building for Windows or macOS everything to be on
  # the same prefix and can be relocated.
  # On other platforms, build-time paths are meaningful.
  if test "$platform" = "osx" ||
     test "$platform" = "windows"; then
     enable_relocatable_bundle=yes
  else
     enable_relocatable_bundle=no
  fi
fi
if test "x$enable_relocatable_bundle" = xyes; then
  AC_DEFINE(ENABLE_RELOCATABLE_RESOURCES, 1,
            [Define to 1 if resources are considered bundled under the same prefix])
fi
AM_CONDITIONAL(ENABLE_RELOCATABLE_RESOURCES, test x$enable_relocatable_bundle = xyes)

AC_MSG_RESULT([$enable_relocatable_bundle])

AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile src/Makefile subprojects/Makefile data/Makefile pixmaps/Makefile pixmaps/icons/Makefile platform-specific/Makefile platform-specific/linux/Makefile po/Makefile.in])
AC_OUTPUT
