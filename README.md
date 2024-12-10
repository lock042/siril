# Siril

> Copyright &copy; 2012-2023, Team free-astro
> <<https://free-astro.org/index.php/Siril>>
> <<https://www.siril.org>>

## Summary

Siril is an astronomical image processing tool.

It is specially tailored for noise reduction and improving the signal/noise
ratio of an image from multiple captures, as required in astronomy.
Siril can align automatically or manually, stack and enhance pictures from various file formats,
even image sequence files (films and SER files).
It works well with limited system resources, like in embedded platforms, but is
also very fast when run on more powerful computers.

Contributors are welcome. Programming language is C, with parts in C++.
Main development is done with most recent versions of libraries.

If you use Siril in your work, please cite this software using the following information:
C. Richard et al., Journal of Open Source Software, 2024, 9(102), 7242. DOI:
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07242/status.svg)](https://doi.org/10.21105/joss.07242).

## Requirements

For compilation, these tools are needed in addition to the base development packages:
- **meson**
- **ninja**
- **cmake**

Then, mandatory build dependencies:
- **glib-2.0** (>= 2.56.0) Glib Convenience Library
- **GTK+ 3**, (>= 3.20) as GUI toolkit
- **gtksourceview4** to provide context highlighting in the script editor
- **lcms2** for color space management
- **cfitsio** for FITS image read and write
- **wcslib** (>=7.12) to handle astrometric data
- **fftw3** for Fourier transforms
- **GSL** (The GNU Scientific Library) for PSF implementation, histograms and
background extraction
- **A C++ compiler** for opencv code and avi exporter
- **libopencv** for various image transformation algorithms (>= 4.4, 4.2 is
possible without some shift-only registration)
- **json-glib-1.0**, (>= 1.2.6) for Siril update check, spectrophotometry
color calibration and metadata output

Siril works internally with FITS files, but other file formats can be used as
input and converted using the conversion tab of the control window. Some file
formats are handled internally, like BMP, PPM and SER, some require external
libraries listed below. Libraries need to be present at compilation time, or
their support won't be included.

- **libcurl** for web access
- **exiv2** to get thumbnails from files
- **libraw** for DSLR RAW files import
- **libffms2** for films import (any format supported by ffmpeg)
- **libtiff** (>= 4) for TIFF format support
- **libXISF** (>=0.2.7) and **zstd** for XISF format support
- **libjpeg** or compatible libraries like libjpeg-turbo for JPEG format support
- **libjxl** for JPEG XL format support
- **libheif** for HEIF format files import
- **libpng** (>= 1.6) for PNG format support
- **libavformat**, **libavutil** (>= 55.20), **libavcodec**, **libswscale** and **libswresample** for avi export (usually provided by ffmpeg)
- **libgit2** for git integration to sync with the official siril-scripts repository
- **criterion** for unit testing with meson (development)

All these libraries and programs are available in most Linux distributions and
free systems, maybe with the exception of ffms2 that is not as popular as the
others and may need to be compiled.

At runtime, you need a functional Python installation (>=3.9) including the
python3-venv and python3-pip modules, as well as python3-tk to support scripts with
GUIs. If you are using a prebuilt Siril package these will be included, but if you
are compiling from source you need to ensure these are available by installing the
appropriate packages for your operating system.

## Scripting

Siril accepts commands from the graphical command line, from scripts as a file
that contains a sequence of commands, or from a named pipe. The list of
supported commands is documented
[here](https://free-astro.org/index.php?title=Siril:Commands). We recommend to use
the `siril-cli` binary for that as no X-server is needed.

Some general purpose scripts have been made by the developers and some power
users, and are provided with the source code and some installers. When they are
in some default directories or in the directories defined by the user in the
settings, scripts appear in a top-menu of the main window. See [this
page](https://free-astro.org/index.php?title=Siril:scripts) for more
information about scripts and a list of scripts ready for use.

The named pipe is only enabled when using Siril in a non-graphical mode and when
the `-p` argument is passed to the program on the command line.

## Download binaries

We maintain binary packages of the latest stable version of Siril. A full list
of all releases is available on
[free-astro](https://free-astro.org/index.php?title=Siril:releases) as well as
[Siril's website](https://siril.org/download/). The available packages per
relesae might differ.

- Ubuntu and Linux Mint: [`ppa:lock042/siril`](https://launchpad.net/~lock042/+archive/ubuntu/siril)
- Windows: see [Downloads](https://siril.org/download/)
- macOS: see [Downloads](https://siril.org/download/)

## Download source

You can get Siril's source code from the release archives on the webpage or
fetch the latest version from our [repository on
GitLab](https://gitlab.com/free-astro/siril):

```bash
git clone --recurse-submodules https://gitlab.com/free-astro/siril.git
```

## Building Siril for GNU/Linux

Siril uses the Meson build system. Run the following commands:

```bash
# adjust the prefix to wherever you want to install Siril
meson setup --prefix /usr/local --buildtype release _build
ninja -C _build
ninja -C _build install
```

To update your build/installation, run the following commands:

```bash
git pull
git submodule update
ninja -C _build install
```

To uninstall Siril, run the following command:

```bash
ninja -C _build uninstall
```

Using meson to build siril requires all optional dependencies to be available or explicitly
disabled on the meson command line adding `-Dexiv2=false` for example. The autotools way
still only enables dependencies that are found and is available using autogen.sh.

## Building Siril for macOS

The official JHBuild-based build scripts are available in the
[siril_macos](https://gitlab.com/free-astro/siril_macos) repository.

Alternatively, you can install Siril via [Homebrew](https://brew.sh). Please
note that this is not maintained by Siril developers.

```bash
brew install siril
```

## Building Siril for Windows

The build process using [msys2](https://www.msys2.org) is documented
[here](https://free-astro.org/index.php?title=Siril:install#Building_on_Windows_with_msys2).

## Translation instructions for Siril

If you're interested in contributing to the translation of the application and documentation,
we encourage you to use [Weblate](https://weblate.pixls.us/projects/siril/). Weblate is a
powerful web-based translation tool that allows for easy collaboration and efficient
translation workflow.

## Notes on Siril FITS image format

[Flexible Image Transport System (FITS)](https://en.wikipedia.org/wiki/FITS) is an open
standard defining a digital file format useful for storage, transmission and processing
of scientific and other images.
FITS is the most commonly used digital file format in astronomy.

Since FITS is a container and doesn't specify the order and size of data, it's
useful to fix it at some point. Currently, Siril uses 32-bit floating point per
channel values (TFLOAT), and images are stored channel after channel on a
bottom-to-top, left-to-right order. The convention chosen is the same as professional
tools, like ds9 (Harvard Smithsonian Center for Astrophysics) and fv
(FITS viewer from NASA) that store images bottom-up too. More details are
described [here](https://free-astro.org/index.php?title=Siril:FITS_orientation).

All files imported and converted in Siril or files exported by Siril are in this
FITS format, except sequence files like SER and films, which are read from the
file and converted on-the-fly. Films should be converted to SER now to process
them, many parallel operations are unsupported on them.

## Notes on image sequence files

Siril makes a strong case for the use SER sequences against the generic film
containers that are not well suited for astronomy data and that may not be read
the same way by different players. Siril can convert any film format supported
by FFMS2 (probably all ffmpeg formats, which is a lot) to SER, and even any
image sequence to SER.

Siril supports SER v3. See [this page](https://free-astro.org/index.php/SER) for more details.

## Useful links

- [Project Homepage](https://www.siril.org)
- [Documentation](https://siril.rtfd.io)
- [Forum](https://discuss.pixls.us/siril)
- [Releases and Downloads](https://free-astro.org/index.php?title=Siril:releases)
- [Report a bug](https://gitlab.com/free-astro/siril/issues)
- [Commands reference](https://siril.readthedocs.io/en/latest/genindex.html)

## License

[GPL-3.0-or-later](LICENSE.md)
