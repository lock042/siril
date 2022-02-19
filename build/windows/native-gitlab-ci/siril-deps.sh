pacman --noconfirm -S --needed base-devel \
mingw-w64-x86_64-toolchain \
mingw-w64-x86_64-cmake \
git \
automake \
mingw-w64-x86_64-lcms2 \
mingw-w64-x86_64-curl \
mingw-w64-x86_64-json-glib \
mingw-w64-x86_64-meson \
mingw-w64-x86_64-ninja \
mingw-w64-x86_64-fftw \
mingw-w64-x86_64-exiv2 \
mingw-w64-x86_64-gtk3 \
mingw-w64-x86_64-libconfig \
mingw-w64-x86_64-gsl \
mingw-w64-x86_64-opencv \
mingw-w64-x86_64-libheif \
mingw-w64-x86_64-ffms2 \
mingw-w64-x86_64-cfitsio \

# Build LibRaw from github
mkdir _deps && cd _deps
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
/usr/bin/autoreconf -fi && \
./configure --disable-examples --disable-static && \
make install || exit 1
cd ..

