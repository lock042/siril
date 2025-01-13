# development system, compiler toolchain and build tools
pacman -S --noconfirm  --needed \
base-devel \
mingw-w64-x86_64-toolchain \
mingw-w64-x86_64-cmake \
git \
automake \
mingw-w64-x86_64-autotools \
mingw-w64-x86_64-meson \
mingw-w64-x86_64-ninja \
mingw-w64-x86_64-ccache

# Siril required dependencies (either Siril or subprojects)
pacman -S --noconfirm  --needed \
mingw-w64-x86_64-glib2 \
mingw-w64-x86_64-gtk3 \
mingw-w64-x86_64-gsl \
mingw-w64-x86_64-lcms2 \
mingw-w64-x86_64-json-glib \
mingw-w64-x86_64-fftw \
mingw-w64-x86_64-cfitsio \
mingw-w64-x86_64-opencv \
mingw-w64-x86_64-cairo \

# Siril optional dependencies
pacman -S --noconfirm  --needed \
mingw-w64-x86_64-exiv2 \
mingw-w64-x86_64-libraw \
mingw-w64-x86_64-libtiff \
mingw-w64-x86_64-libjpeg-turbo \
mingw-w64-x86_64-libjxl \
mingw-w64-x86_64-libpng \
mingw-w64-x86_64-libheif \
mingw-w64-x86_64-libxisf \
mingw-w64-x86_64-libgit2-winhttp \
mingw-w64-x86_64-ffmpeg \
mingw-w64-x86_64-ffms2 \
mingw-w64-x86_64-curl

# # Build LibRaw from github
# mkdir _deps && cd _deps
# git clone --depth 1 https://github.com/LibRaw/LibRaw.git
# cd LibRaw
# autoreconf -fi && \
# ./configure --disable-examples --disable-static && \
# make install -j$(nproc) || exit 1
# cd ..

