SIRIL
=====

> Copyright &copy; 2012-2020, Team free-astro
> <<https://free-astro.org/index.php/Siril>>
> <<https://www.siril.org>>

Summary
-------
SIRIL is an astronomical image processing tool.

It is specially tailored for noise reduction and improving the signal/noise
ratio of an image from multiple captures, as required in astronomy.
SIRIL can align automatically or manually, stack and enhance pictures from various file formats,
even image sequence files (films and SER files).
It works well with limited system resources, like in embedded platforms, but is
also very fast when run on more powerful computers.

Contributors are welcome. Programming language is C, with parts in C++.
Main development is done with most recent versions of libraries.

Requirements
------------
For compilation, these tools are needed in addition to the base development packages:
 * **meson**
 * **ninja**
 * **cmake**

Then, mandatory build dependencies:
 * **GTK+ 3**, (>= 3.20) as GUI toolkit
 * **cfitsio** for FITS image read and write
 * **fftw3** for Fourier transforms
 * **GSL** (The GNU Scientific Library) for PSF implementation, histograms and background extraction
 * **A C++ compiler** for opencv code and avi exporter
 * **libopencv** for various image transformation algorithms (>= 4.4, 4.2 is possbile without some shift-only registration)
 * **exiv2** to manage image metadata

SIRIL works internally with FITS files, but other file formats can be used as
input and converted using the conversion tab of the control window. Some file
formats are handled internally, like BMP, PPM and SER, some require external
libraries or programs listed below. Libraries need to be present at compilation
time, or their support won't be included.

 * **json-glib-1.0**, (>= 1.2.6) for Siril update check
 * **libraw** for DSLR RAW files import
 * **libffms2** for films import (any format supported by ffmpeg)
 * **libtiff** (>= 4) for TIFF format support
 * **libjpeg** or compatible libraries like libjpeg-turbo for JPEG format support
 * **libheif** for HEIF format files import
 * **libpng** (>= 1.6) for PNG format support
 * **libavformat**, **libavutil** (>= 55.20), **libavcodec**, **libswscale** and **libswresample** for avi export (usually provided by ffmpeg)
 * **libcurl** for web interaction. Useless on GNU-Linux MUST be installed on macOS and Windows platform as GIO is broken
 * **wcslib** for some astrometry utilities
 * **gnuplot** for photometry graphs output
 * **libconfig** (>= 1.4) to read old configuration files (not used since 1.1)
 * **criterion** for unit testing (development)

All these libraries and programs are available in most Linux distributions and
free systems, maybe with the exception of ffms2 that is not as popular as
others and may need to be compiled.

Scripting
---------

SIRIL accepts commands from the graphical command line, from scripts as a file
that contains a sequence of commands, or from a named pipe. The list of
supported commands is documented
[here](https://free-astro.org/index.php?title=Siril:Commands). We recommend to use
the **siril-cli** binary for that as no X-server is needed.

Some general purpose scripts have been made by the developers and some power
users, and are provided with the source code and some installers. When they are
in some default directories or in the directories defined by the user in the
settings, scripts appear in a top-menu of the main window. See [this
page](https://free-astro.org/index.php?title=Siril:scripts) for more
information about scripts and a list of scripts ready for use.

The named pipe is only enabled when using SIRIL in a non-graphical way and when
the `-p` command is passed to the program command line.

Source download
---------------

You can get SIRIL source code from the release archives on their webpage, or the latest version from git:

    git clone https://gitlab.com/free-astro/siril.git

So far, we are using submodule for the use of some algorithms. You must therefore run the following commands:

    cd siril
    git submodule sync --recursive
    git submodule update --init --recursive

Building SIRIL for GNU/Linux
----------------------------
The process now uses the meson build system that is faster and more modern than autotools.
Run with the following commands:

    meson --buildtype release _build
    ninja -C _build
    ninja -C _build install

To install into an arbitrary path - for testing, or where root privilege
is not available - a prefix may be supplied to meson:
`meson --prefix /tmp/siril-testing-something ...` Then build and install,
as normal.

To update Siril, run the following commands

    git pull --recurse-submodules
    ninja -C _build install

To uninstall Siril, run the following command:

    ninja -C _build uninstall

Note that a binary package for stable version of SIRIL is maintained for Debian.
PPA repositories for Ubuntu and Linux Mint maintained by SIRIL's authors are
available in **ppa:lock042/siril**.

See the [download](https://free-astro.org/index.php?title=Siril:releases) page
of the current version for other packages that could be available.

Building SIRIL for macOS
------------------------
We provide a dmg installer on the [website](https://www.siril.org/download/),
but you can also install SIRIL from sources using homebrew.

    brew install siril

SIRIL on Microsoft Windows
----------------
SIRIL is supported on Microsoft Windows since version 0.9.8.  We provide binary files
in an installer.

You can also build it from source yourself with msys2, it is documented
[here](https://free-astro.org/index.php?title=Siril:install#Installing_on_Windows).

Translation instructions for SIRIL
----------------------------------
We now use weblate for the translation of the application. Go to this [address](https://weblate.pixls.us/projects/siril/) and create an account.

Notes on SIRIL FITS image format
--------------------------------
[Flexible Image Transport System (FITS)](https://en.wikipedia.org/wiki/FITS) is an open
standard defining a digital file format useful for storage, transmission and processing
of scientific and other images.
FITS is the most commonly used digital file format in astronomy.

Since FITS is a container and doesn't specify the order and size of data, it's
useful to fix it at some point. Currently, SIRIL uses 32-bit floating point per
channel values (TFLOAT), and images are stored channel after channel on a
bottom-to-top, left-to-right order. The convention chosen is the same as professional
tools, like ds9 (Harvard Smithsonian Center for Astrophysics) and fv
(FITS viewer from NASA) that store images bottom-up too. More details are
described [here](https://free-astro.org/index.php?title=Siril:FITS_orientation).

All files imported and converted in SIRIL or files exported by SIRIL are in this
FITS format, except sequence files like SER and films, which are read from the
file and converted on-the-fly. Films should be converted to SER now to process
them, many parallel operations are unsupported on them.

Notes on image sequence files
-----------------------------
SIRIL makes a strong case for the use SER sequences against the generic film
containers that are not well suited for astronomy data and that may not be read
the same way by different players. SIRIL can convert any film format supported
by FFMS2 (probably all ffmpeg formats, which is a lot) to SER, and even any
image sequence to SER.

SIRIL supports SER v3. See [this page](https://free-astro.org/index.php/SER) for more details.

Useful links
------------
 * [Project Homepage](https://www.siril.org)
 * [Documentation](https://free-astro.org/siril_doc-en)
 * [Forum](https://discuss.pixls.us/siril)
 * [Releases and Downloads](https://free-astro.org/index.php?title=Siril:releases)
 * [Report a bug](https://gitlab.com/free-astro/siril/issues)
 * [Supported commands](https://free-astro.org/index.php?title=Siril:Commands)
