# This file is part of Siril, an astronomy image processor.
# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
#
# Siril is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Siril is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siril. If not, see <http://www.gnu.org/licenses/>.


################################################################################
# Meson main file
################################################################################


################################################################################
# Define project

project('siril',
  ['c','cpp'],
  version : '1.3.5',
  license : ['GPL3'],
  meson_version: '>=0.56.0',
  default_options : [
    'c_std=gnu11',
    'cpp_std=c++17',
    'buildtype=debugoptimized',
  ],
)

rootInclude = include_directories('.')
srcInclude  = include_directories('src')

warnings = []

have_fftw3fomp = false
have_fftw3fthreads = false
have_fftw3fthreadsnolinkerarg = false

enable_openmp    = get_option('openmp')
enable_libcurl   = get_option('libcurl')
enable_exiv      = get_option('exiv2')
enable_libraw    = get_option('libraw')
enable_libtiff   = get_option('libtiff')
enable_libjpeg   = get_option('libjpeg')
enable_libpng    = get_option('libpng')
enable_libheif   = get_option('libheif')
enable_libXISF   = get_option('libXISF')
enable_ffms2     = get_option('ffms2')
enable_ffmpeg    = get_option('ffmpeg')
enable_criterion = get_option('criterion')
enable_libjxl    = get_option('libjxl')
enable_libgit2    = get_option('libgit2')

# General version
siril_version = meson.project_version()

version_array = siril_version.split('.')
siril_version_major = version_array[0].to_int()
siril_version_minor = version_array[1].to_int()
siril_version_micro = version_array[2]

stable = (siril_version_minor % 2 == 0)

## Get version from git
git = find_program('git', required: false)

if not git.found()
  message('Git not installed or git directory not found!')
  git_version = 'Unknown'
  git_last_commit_year = 'Unknown'
else
  r = run_command('git', 'show', '-s', '--format="%ci %h"','-n1', check: false)
  if r.returncode() == 0
    git_info = r.stdout().strip().split('"')[1]
    git_info_array = git_info.split(' ')
    git_version = git_info_array[3]
    git_last_commit_year = git_info_array[0].split('-')[0]
    message('Siril version:' + git_version)
    errortxt = r.stderr().strip()
  else
    message('Cannot retrieve git version')
    message(r.stderr().strip())
    git_version = 'Unknown'
    git_last_commit_year = 'Unknown'
  endif
endif

################################################################################
# Host system detection

host_cpu_family = host_machine.cpu_family()
host_os = host_machine.system().to_lower()
message('Host os: ' + host_os)

platform_linux = (
  host_os.contains('linux')
)

platform_windows = (
  host_os.contains('mingw') or
  host_os.contains('cygwin') or
  host_os.contains('windows')
)

platform_osx = (
  host_os.contains('machten') or
  host_os.contains('rhapsody') or
  host_os.contains('darwin')
)

platform_bsd = (
  host_os.contains('bsd')
)

################################################################################
# Meson modules and compilers

## Modules
cmake = import('cmake')
i18n = import('i18n')
if platform_windows
  winm = import('windows')
endif
gnome = import('gnome')

## Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

################################################################################
# Set compiler and linker flags

siril_cpp_flag = []
siril_c_flag = []
siril_link_arg = []


# -fno-math-errno is used to speedup code, this for example speeds up all sqrt() functions.
# There is no case where we check the math errno. So it’s safe to use this flag for siril.
siril_c_flag += ['-DHAVE_CONFIG_H', '-fcommon', '-fno-math-errno', '-D_FILE_OFFSET_BITS=64']
siril_cpp_flag += ['-DHAVE_CONFIG_H', '-fcommon', '-fno-math-errno', '-D_FILE_OFFSET_BITS=64']

if platform_linux
  siril_link_arg += ['-rdynamic']
elif platform_bsd
  siril_link_arg += '-Wl,--no-as-needed' # workaround for https://gitlab.com/free-astro/siril/-/issues/614 (Bug#251206)
  siril_link_arg += '-rdynamic'
  siril_link_arg += '-lexecinfo'
elif platform_windows
  #TODO: Do Meson provide more elegant way?
  siril_link_arg += '-Wl,--export-all-symbols'
  siril_link_arg += '-lpsapi'
  siril_link_arg += '-ldbghelp'
  siril_link_arg += '-lbcrypt'
  siril_c_flag += ['-D_WIN32_WINNT=0x0600', '-DWINVER=0x0600']
  siril_cpp_flag += ['-D_WIN32_WINNT=0x0600', '-DWINVER=0x0600']
elif platform_osx
  # on OSX ObjC and C sources are mixed so adding objc to the linkflags
  siril_link_arg += ['-Wl,-framework,Foundation', '-Wl,-framework,AppKit', '-Wl,-framework,Security', '-ObjC']
  add_project_arguments('-ObjC', language: [ 'c'])
  add_project_link_arguments(siril_link_arg, language: [ 'c', 'cpp' ])
else
  error('Unsupported OS.')
endif

# sse2 assembly
have_sse2 = false
if cc.has_argument('-msse2') and cpp.has_argument('-msse2')
  message('sse2 assembly available')
  if cc.compiles('asm ("punpckhwd %xmm0,%xmm1");')
    have_sse2 = true
  endif
endif

# sse4.1 assembly
have_sse4_1 = false
if cc.has_argument('-msse4.1') and cpp.has_argument('-msse4.1')
  message('sse4.1 assembly available')
  if cc.compiles('asm ("pmovzxbd %xmm0,%xmm1");')
    have_sse4_1 = true
  endif
endif

################################################################################
# Directories

siril_share = get_option('prefix') + '/' + get_option('datadir') + '/'
siril_datadir = siril_share + meson.project_name() + '/'
siril_locale = get_option('prefix') + '/' + get_option('localedir') + '/'

################################################################################
# Dependencies

if get_option('relocatable-bundle') == 'yes'
  relocatable_bundle = true
elif get_option('relocatable-bundle') == 'no'
  relocatable_bundle = false
else # == 'platform-default'
  # By default, assume building for Windows or macOS everything to be on
  # the same prefix and can be relocated.
  # On other platforms, build-time paths are meaningful.
  if platform_windows or platform_osx
      relocatable_bundle = true
  else
      relocatable_bundle = false
  endif
endif

# Dependencies
no_dep = dependency('', required: false)

## Core dependencies
dbghelp = platform_windows ? cc.find_library('dbghelp') : no_dep
glib_dep = dependency('glib-2.0', required : true, version: '>= 2.56.0')
gio_dep = dependency('gio-2.0', required : true, version: '>= 2.56.0')
if platform_windows
  gio_platform_dep = dependency('gio-windows-2.0', required : true, version: '>= 2.56.0')
else
  gio_platform_dep = dependency('gio-unix-2.0', required : true, version: '>= 2.56.0')
endif
gtk_dep = dependency('gtk+-3.0', required : true, version: '>= 3.22.0')
gsl_dep = dependency('gsl', required : true)
m_dep = cc.find_library('m', required : false)
# lcms2 dependency for color management
lcms_dep = dependency('lcms2', required : true, version: '>= 2.14')
# wcslib mandatory - uses subproject fallback to build if not found
wcslib_dep = dependency('wcslib', fallback : ['wcslib', 'wcslib_dep'], required : true, version: '>= 7.12')
json_glib_dep = dependency('json-glib-1.0', required : true, version: '>=1.2.6')

fftw3f_dep = dependency('fftw3f', required : true)
# threaded fftw3f
if cc.links('''#include <fftw3.h>
               int main () {
                 fftwf_init_threads();
                 return 0;
               }''', dependencies : [fftw3f_dep], name : 'fftw3f linked with pthreads and combined_threads')
  have_fftw3fthreadsnolinkerarg = true
elif cc.links('''#include <fftw3.h>
               int main () {
                 fftwf_init_threads();
                 return 0;
               }''', dependencies : [fftw3f_dep], args : ['-lfftw3f_omp'], name : 'fftw3f linked with openmp')
  have_fftw3fomp = true
elif cc.links('''#include <fftw3.h>
               int main () {
                 fftwf_init_threads();
                 return 0;
               }''', dependencies : [fftw3f_dep], args : ['-lfftw3f_threads'], name : 'fftw3f linked with pthreads')
  have_fftw3fthreads = true
endif

if platform_windows
  cfitsio_dep = dependency('cfitsio', required : true, version: '>= 4.4.0') # for #1324
else
  cfitsio_dep = dependency('cfitsio', required : true)
endif
cairo_dep = dependency('cairo', required : true)
kplot_dep = dependency('kplot', fallback : ['kplot', 'kplot_dep'], required : true)
opencv4_dep = dependency('opencv4', required : true, version: '>= 4.2.0')
# checking required opencv modules
opencv_symbols = ['HAVE_OPENCV_CORE', 'HAVE_OPENCV_IMGPROC', 'HAVE_OPENCV_CALIB3D', 'HAVE_OPENCV_STITCHING']
opencv_incdir = '-I'+opencv4_dep.get_variable(pkgconfig : 'includedir')
foreach s : opencv_symbols
  cpp.has_header_symbol('opencv2/opencv_modules.hpp', s, args: [opencv_incdir], required : true)
endforeach

rtprocess_dep = dependency('rtprocess', required : false)
if not rtprocess_dep.found()
  opt_var = cmake.subproject_options()

  # Call CMake with `-DBUILD_SHARED_LIBS=OFF`
  opt_var.add_cmake_defines({'BUILD_SHARED_LIBS': false})
  opt_var.add_cmake_defines({'OPTION_OMP': false})
  if get_option('buildtype').startswith('debug')
    opt_var.add_cmake_defines({'CMAKE_BUILD_TYPE': 'debug'})
  endif
  sub_proj = cmake.subproject('librtprocess', options: opt_var)

  rtprocess_dep = sub_proj.dependency('rtprocess')
  rtprocess_inc = sub_proj.include_directories('rtprocess')
else
  siril_cpp_flag += '-DHAVE_LIBRTPROCESS'
endif

htmesh_dep = dependency('htmesh', required : false)
if not htmesh_dep.found()
  opt_var = cmake.subproject_options()

  # Call CMake with `-DCMAKE_CXX_FLAGS="-I$(abs_srcdir)/htmesh"
  opt_var.add_cmake_defines({'CMAKE_CXX_FLAGS': '-I'+meson.current_source_dir()+'/subprojects/htmesh'})
if platform_windows
  opt_var.append_compile_args('cpp', '-Wno-unknown-pragmas')
endif
  sub_proj = cmake.subproject('htmesh', options: opt_var)

  htmesh_dep = sub_proj.dependency('htmesh')
  htmesh_inc = sub_proj.include_directories('htmesh')
else
  siril_cpp_flag += '-DHAVE_HTMESH'
endif

# healpix dependency to support Gaia catalogues
# healpix_dep = dependency('healpix_cxx', required: false)
# if not healpix_dep.found()
#   opt_var = cmake.subproject_options()

#   # Call CMake with `-DCMAKE_CXX_FLAGS="-I$(abs_srcdir)/htmesh"
#   opt_var.add_cmake_defines({'CMAKE_CXX_FLAGS': '-I'+meson.current_source_dir()+'/subprojects/healpix_cxx'})
#   sub_proj = cmake.subproject('healpix_cxx', options: opt_var)

#   healpix_dep = sub_proj.dependency('healpix_cxx')
#   healpix_inc = sub_proj.include_directories('healpix_cxx')

# endif

healpix_dep = dependency('healpix_cxx', fallback : ['healpix_cxx', 'healpix_cxx_dep'], required : true)

## Optional dependencies

openmp_dep = []
if enable_openmp
  if platform_osx
    message('Set OpenMP flags for macOS')
    siril_link_arg += ['-lomp']
    siril_cpp_flag += ['-Xpreprocessor', '-fopenmp']
    siril_c_flag += ['-Xpreprocessor', '-fopenmp']
  else
  openmp_dep = dependency('openmp', required : true)
  endif
else
  message('OpenMP disabled')
endif

# Now we know OpenMP availability we can choose FFTW flags
if enable_openmp
  if have_fftw3fomp
    siril_c_flag += '-DHAVE_FFTW3F_OMP'
    siril_cpp_flag += '-DHAVE_FFTW3F_OMP'
    siril_link_arg += ['-lfftw3f_omp']
    message('FFTW linked against OpenMP')
  elif have_fftw3fthreads
    siril_c_flag += '-DHAVE_FFTW3F_THREADS'
    siril_cpp_flag += '-DHAVE_FFTW3F_THREADS'
    siril_link_arg += ['-lfftw3f_threads']
    message('FFTW linked against pthreads: will use OpenMP callback')
  elif have_fftw3fthreadsnolinkerarg
# Note in this instance FFTW has been built with the combined_threads
# option which means pthreads is in use
    siril_c_flag += '-DHAVE_FFTW3F_THREADS'
    siril_cpp_flag += '-DHAVE_FFTW3F_THREADS'
    message('FFTW linked against pthreads: will use OpenMP callback')
  endif
elif have_fftw3fthreads
  siril_c_flag += '-DHAVE_FFTW3F_THREADS'
  siril_cpp_flag += '-DHAVE_FFTW3F_THREADS'
  siril_link_arg += ['-lfftw3f_threads']
  message('FFTW linked against pthreads')
elif have_fftw3fthreadsnolinkerarg
  siril_c_flag += '-DHAVE_FFTW3F_THREADS'
  siril_cpp_flag += '-DHAVE_FFTW3F_THREADS'
  message('FFTW linked against pthreads')
else
  message('FFTW does not support multithreading')
endif

exiv_dep = no_dep
if enable_exiv
  exiv_dep = dependency('exiv2', required : true, version: '>= 0.25')
endif

libraw_dep = no_dep
if enable_libraw
  libraw_dep = dependency('libraw', required : true)
endif

libtiff_dep = no_dep
if enable_libtiff
  libtiff_dep = dependency('libtiff-4', required : true)
endif

libjpeg_dep = no_dep
if enable_libjpeg
  libjpeg_dep = dependency('libjpeg', required : true)
endif

libjxl_dep = no_dep
libjxl_threads_dep = no_dep
if enable_libjxl
  libjxl_dep = dependency('libjxl', required : true)
  libjxl_threads_dep = dependency('libjxl_threads', required : false)
endif

libpng_dep = no_dep
if enable_libpng
  libpng_dep = dependency('libpng', required : true, version: '>= 1.6.0')
endif

libheif_dep = no_dep
if enable_libheif
  libheif_dep = dependency('libheif', required : true, version: '>= 1.8.0')
endif

libXISF_dep = no_dep
if enable_libXISF
  libXISF_dep = dependency('libxisf', required : true, version: '>= 0.2.7')
endif

libgit2_dep = no_dep
if enable_libgit2
  libgit2_dep = dependency('libgit2', required : true)
endif

ffms2_dep = no_dep
if enable_ffms2
  ffms2_dep = dependency('ffms2', required : true)
endif

ffmpeg = []
if enable_ffmpeg
  ffmpeg = [
    dependency('libavformat', required : true),
    dependency('libavcodec', required : true),
    dependency('libavutil', version : '>= 55.20', required : true),
    dependency('libswscale', required : true),
    dependency('libswresample', required : true)
  ]
endif

libcurl_dep = no_dep
if enable_libcurl
  libcurl_dep = dependency('libcurl', required : true)
endif

criterion_dep = no_dep
if enable_criterion
  criterion_dep = dependency('criterion', required : true)
endif

################################################################################
# Configuration

## General and package configuration
conf_data = configuration_data()

conf_data.set_quoted('CPU_ARCH', host_cpu_family)

conf_data.set('SIRIL_UNSTABLE', stable ? false : 1)
conf_data.set('SIRIL_OUTPUT_DEBUG', get_option('buildtype').startswith('debug') ? 1 : false)
conf_data.set('ENABLE_RELOCATABLE_RESOURCES', relocatable_bundle ? 1 : false)
conf_data.set('OS_OSX', platform_osx ? 1 : false)

conf_data.set('SIRIL_MAJOR_VERSION', siril_version_major)
conf_data.set('SIRIL_MINOR_VERSION', siril_version_minor)
pre_release_arr = siril_version_micro.split('-')
conf_data.set('SIRIL_MICRO_VERSION', pre_release_arr[0])

conf_data.set_quoted('PACKAGE', meson.project_name(), description : 'Name of package.')
conf_data.set_quoted('PACKAGE_NAME', meson.project_name(), description : 'Define to the full name of this package.')
conf_data.set_quoted('PACKAGE_VERSION', siril_version, description : 'Define to the home page for this package.')
conf_data.set_quoted('VERSION', siril_version, description : 'Define to the home page for this package.')
conf_data.set_quoted('PACKAGE_SOURCE_DIR', 'src', description : 'Package source dir.')
conf_data.set_quoted('PACKAGE_STRING', meson.project_name() + ' ' + siril_version, description : 'Define to the full name and version of this package.')
conf_data.set_quoted('PACKAGE_TARNAME', meson.project_name(), description : 'Define to the one symbol short name of this package.')
conf_data.set_quoted('PACKAGE_URL', 'https://www.siril.org', description : 'Define to the home page for this package.')
conf_data.set_quoted('PACKAGE_BUGREPORT', 'https://gitlab.com/free-astro/siril/issues', description : 'Define to the address where bug reports for this package should be sent.')
conf_data.set_quoted('PACKAGE_DATA_DIR', siril_datadir, description : 'Define to the data directory.')
conf_data.set_quoted('LOCALEDIR', siril_locale, description : 'Define to the locale directory.')

# Man page platform specific

if platform_windows
  conf_data.set('MAN_CONFIG_PATH', '.B %LOCALAPPDATA%\esiril\econfig.ini')
  conf_data.set('MAN_CSS_PATH', '.B C:\eProgram Files\eSiril\eshare\esiril\esiril.css')
elif platform_osx
  conf_data.set('MAN_CONFIG_PATH', '.B ~/Library/Application Support/siril/config.ini')
  conf_data.set('MAN_CSS_PATH', '.B ' + siril_share + 'siril/siril.css')
else
  conf_data.set('MAN_CONFIG_PATH', '.B ~/.config/siril/config.ini')
  conf_data.set('MAN_CSS_PATH', '.B ' + siril_share + 'siril/siril.css')
endif


# Build identifiers #

if platform_linux
  conf_data.set_quoted('SIRIL_BUILD_PLATFORM_FAMILY', 'linux')
elif platform_windows
  conf_data.set_quoted('SIRIL_BUILD_PLATFORM_FAMILY', 'windows')
elif platform_osx
  conf_data.set_quoted('SIRIL_BUILD_PLATFORM_FAMILY', 'macos')
elif platform_bsd
  conf_data.set_quoted('SIRIL_BUILD_PLATFORM_FAMILY', 'bsd')
else
  conf_data.set_quoted('SIRIL_BUILD_PLATFORM_FAMILY', 'other')
endif

if get_option('flatpak-id') != ''
  conf_data.set_quoted('FLATPAK_ID', get_option('flatpak-id'))
endif

# Check for available headers
foreach header : [
    { 'm': 'HAVE_EXECINFO_H',     'v': 'execinfo.h' },
    { 'm': 'HAVE_SYS_STATVFS_H',  'v': 'sys/statvfs.h' },
    { 'm': 'HAVE_SYS_VFS_H',      'v': 'sys/vfs.h' },
    { 'm': 'HAVE_SYS_PARAM_H',    'v': 'sys/param.h' },
  ]
  conf_data.set(header['m'], cc.has_header(header['v']) ? 1 : false)
endforeach

## Dependencies configuration
if opencv4_dep.version().version_compare('>=4.4.0')
  conf_data.set('HAVE_CV44', 1, description : 'Using OpenCV 4.4 and all registration features.')
else
  warning('OpenCV 4.4 was not found, using a lower version that does not support translation-only registration instead.')
endif

if gsl_dep.found()
  if gsl_dep.version().version_compare('>=2.0')
    conf_data.set('HAVE_GSL_1', 0, description : 'Use GSL 1.')
    conf_data.set('HAVE_GSL_2', 1, description : 'Use GSL 2.')
  else
    conf_data.set('HAVE_GSL_1', 1, description : 'Use GSL 1.')
    conf_data.set('HAVE_GSL_2', 0, description : 'Use GSL 2.')
  endif
else
  conf_data.set('HAVE_GSL_1', 0, description : 'Use GSL 1.')
  conf_data.set('HAVE_GSL_2', 0, description : 'Use GSL 2.')
endif

if exiv_dep.found()
  conf_data.set('HAVE_EXIV2', true, description : 'Using exif previews.')
else
  conf_data.set('HAVE_EXIV2', false, description : 'Using exif previews.')
endif

if libraw_dep.found()
  conf_data.set('HAVE_LIBRAW', true, description : 'Using RAW images.')
else
  conf_data.set('HAVE_LIBRAW', false, description : 'Using RAW images.')
endif

if libtiff_dep.found()
  conf_data.set('HAVE_LIBTIFF', true, description : 'Using TIFF images.')
else
  conf_data.set('HAVE_LIBTIFF', false, description : 'Using TIFF images.')
endif

if libjpeg_dep.found()
  conf_data.set('HAVE_LIBJPEG', true, description : 'Using JPEG images.')
else
  conf_data.set('HAVE_LIBJPEG', false, description : 'Using JPEG images.')
endif

if libjxl_dep.found()
  conf_data.set('HAVE_LIBJXL', true, description : 'Using JPEG XL images.')
else
  conf_data.set('HAVE_LIBJXL', false, description : 'Using JPEG XL images.')
endif

if libjxl_threads_dep.found()
  conf_data.set('HAVE_LIBJXL_THREADS', true, description : 'Using JPEG XL threading.')
else
  conf_data.set('HAVE_LIBJXL_THREADS', false, description : 'Using JPEG XL threading.')
endif

if libpng_dep.found()
  conf_data.set('HAVE_LIBPNG', true, description : 'Using PNG images.')
else
  conf_data.set('HAVE_LIBPNG', false, description : 'Using PNG images.')
endif

if libheif_dep.found()
  conf_data.set('HAVE_LIBHEIF', true, description : 'Using HEIF.')
else
  conf_data.set('HAVE_LIBHEIF', false, description : 'Using HEIF.')
endif

if libXISF_dep.found()
  conf_data.set('HAVE_LIBXISF', true, description : 'Using XISF.')
else
  conf_data.set('HAVE_LIBXISF', false, description : 'Using XISF.')
endif

if libgit2_dep.found()
  conf_data.set('HAVE_LIBGIT2', true, description : 'Using libgit2.')
else
  conf_data.set('HAVE_LIBGIT2', false, description : 'Using libgit2.')
endif

ffmpeg_found = false
foreach dep : ffmpeg
  if not dep.found()
    ffmpeg_found = false
    break
  else
    ffmpeg_found = true
  endif
endforeach
conf_data.set('HAVE_FFMPEG', ffmpeg_found, description : 'Using FFMS.')

if ffms2_dep.found()
  conf_data.set('HAVE_FFMS2', true, description : 'Using FFMS.')
else
  conf_data.set('HAVE_FFMS2', false, description : 'Using FFMS.')
endif

if libcurl_dep.found()
  conf_data.set('HAVE_LIBCURL', true, description : 'Using Curl.')
else
  conf_data.set('HAVE_LIBCURL', false, description : 'Using Curl.')
endif

conf_data.set('HAVE_CXX11', true, description : 'Define to 1 if cpp compiler is C++11 compliant.')

## Git configuration
git_data = configuration_data()
git_data.set_quoted('SIRIL_GIT_VERSION_ABBREV', git_version)
git_data.set_quoted('SIRIL_GIT_LAST_COMMIT_YEAR', git_last_commit_year)

################################################################################
# Externals dependencies

externals_dep = declare_dependency(
  dependencies: [
    dbghelp,
    m_dep,
    glib_dep,
    gio_dep,
    gio_platform_dep,
    gtk_dep,
    lcms_dep,
    libcurl_dep,
    libgit2_dep,
    gsl_dep,
    opencv4_dep,
    rtprocess_dep,
    htmesh_dep,
    healpix_dep,
    exiv_dep,
    wcslib_dep,
    fftw3f_dep,
    cfitsio_dep,
    json_glib_dep,
    libraw_dep,
    libtiff_dep,
    libjpeg_dep,
    libpng_dep,
    libheif_dep,
    libXISF_dep,
    libjxl_dep,
    libjxl_threads_dep,
    ffms2_dep,
    openmp_dep,
    kplot_dep,
    ffmpeg
  ]
)

################################################################################
# Sources directories

subdir('platform-specific')
subdir('build')
subdir('src')
subdir('pixmaps')
subdir('data')
subdir('po')
subdir('scripts')
subdir('3rdparty_licenses')

################################################################################
# Install misc

man_files = [
  configure_file(
    input : meson.project_name() + '.1',
    output: meson.project_name() + '.1',
    configuration: conf_data,
  ),
  configure_file(
    input : meson.project_name() + '.1',
    output: meson.project_name() + '-cli.1',
    configuration: conf_data,
  ),
]

install_man(man_files)
install_data(['AUTHORS', 'LICENSE.md'],
             rename : [
	       'doc/siril/AUTHORS',
	       'doc/siril/LICENSE.md',
	     ],
             install_dir : siril_share)


################################################################################
# Build summary

summary(
  {
    'Building Siril with prefix': get_option('prefix'),
  }, section: 'Directories'
)

summary(
  {
    'sse2'      : have_sse2,
    'sse4.1'    : have_sse4_1,
  }, section: 'Processor extensions'
)

summary(
  {
    'openmp'    : get_option('openmp'),
  }, section: 'Optional features'
)

summary(
  {
    'exiv2'             : exiv_dep.found(),
    'libraw'            : libraw_dep.found(),
    'TIFF'              : libtiff_dep.found(),
    'XISF'              : libXISF_dep.found(),
    'JPEG'              : libjpeg_dep.found(),
    'JPEG XL'           : libjxl_dep.found(),
    'PNG'               : libpng_dep.found(),
    'HEIF'              : libheif_dep.found(),
    'FFMS2'             : ffms2_dep.found(),
    'libgit2'           : libgit2_dep.found(),
    'ffmpeg'            : ffmpeg_found,
    'libcurl'           : libcurl_dep.found(),
    'criterion'         : criterion_dep.found(),
  }, section: 'Optional dependencies'
)

if warnings.length() > 0
  warning('Warnings occurred during configuration')
  foreach warning : warnings
    warning(warning)
  endforeach
endif
