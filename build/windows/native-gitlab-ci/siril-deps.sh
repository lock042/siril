#!/bin/bash

MSYSTEM=${MSYSTEM:-MINGW64}

if [ "$MSYSTEM" = "UCRT64" ]; then
    PREFIX="mingw-w64-ucrt-x86_64"
elif [ "$MSYSTEM" = "MINGW64" ]; then
    PREFIX="mingw-w64-x86_64"
else
    echo "Unsupported MSYSTEM: $MSYSTEM"
    exit 1
fi

# development system, compiler toolchain and build tools
pacman -S --noconfirm  --needed \
base-devel \
${PREFIX}-toolchain \
${PREFIX}-cmake \
git \
automake \
${PREFIX}-autotools \
${PREFIX}-meson \
${PREFIX}-ninja \
${PREFIX}-ccache

# Siril required dependencies (either Siril or subprojects)
pacman -S --noconfirm  --needed \
${PREFIX}-glib2 \
${PREFIX}-gtk3 \
${PREFIX}-gtksourceview4 \
${PREFIX}-gsl \
${PREFIX}-lcms2 \
${PREFIX}-fftw \
${PREFIX}-cfitsio \
${PREFIX}-opencv \
${PREFIX}-cairo \

# Siril optional dependencies
pacman -S --noconfirm  --needed \
${PREFIX}-exiv2 \
${PREFIX}-libraw \
${PREFIX}-libtiff \
${PREFIX}-libjpeg-turbo \
${PREFIX}-libjxl \
${PREFIX}-libpng \
${PREFIX}-libheif \
${PREFIX}-libxisf \
${PREFIX}-libgit2-winhttp \
${PREFIX}-ffmpeg \
${PREFIX}-ffms2 \
${PREFIX}-curl \
${PREFIX}-sqlite3 \

