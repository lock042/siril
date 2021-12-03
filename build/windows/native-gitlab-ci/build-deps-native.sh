#!/bin/bash

set -e

if [[ "$MSYSTEM" == "MINGW32" ]]; then
    export ARTIFACTS_SUFFIX="-w32"
    export MSYS2_ARCH="i686"
    export MSYS_PREFIX="/c/msys64/mingw32/"
else
    export ARTIFACTS_SUFFIX="-w64"
    export MSYS2_ARCH="x86_64"
    export MSYS_PREFIX="/c/msys64/mingw64/"
fi

# Update everything
pacman --noconfirm -Suy

# Install the required packages
pacman --noconfirm -S --needed \
    base-devel \
    mingw-w64-$MSYS2_ARCH-toolchain \
    mingw-w64-$MSYS2_ARCH-meson \
    \
    mingw-w64-$MSYS2_ARCH-fftw \
    mingw-w64-$MSYS2_ARCH-exiv2 \
    mingw-w64-$MSYS2_ARCH-gtk3 \
    mingw-w64-$MSYS2_ARCH-libconfig \
    mingw-w64-$MSYS2_ARCH-gsl \
    mingw-w64-$MSYS2_ARCH-opencv \
    mingw-w64-$MSYS2_ARCH-libheif \
    mingw-w64-$MSYS2_ARCH-ffms2 \
    mingw-w64-$MSYS2_ARCH-cfitsio

export GIT_DEPTH=1
export SIRIL_PREFIX="`realpath ./_install`${ARTIFACTS_SUFFIX}"
export PATH="$SIRIL_PREFIX/bin:$PATH"
export PKG_CONFIG_PATH="${SIRIL_PREFIX}/lib/pkgconfig:$PKG_CONFIG_PATH"
export PKG_CONFIG_PATH="${SIRIL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="${SIRIL_PREFIX}/lib:${LD_LIBRARY_PATH}"
export ACLOCAL_FLAGS="-I/c/msys64/mingw64/share/aclocal"
export XDG_DATA_DIRS="${SIRIL_PREFIX}/share:/mingw64/share/"


mkdir _build
cd _build
meson -Dprefix="${SIRIL_PREFIX}" \
      --buildtype=release
ninja
ninja install
cd ../..
