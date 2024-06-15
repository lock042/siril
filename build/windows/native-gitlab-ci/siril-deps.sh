pacman --noconfirm -S --needed base-devel \
mingw-w64-x86_64-toolchain \
mingw-w64-x86_64-cmake \
git \
automake \
mingw-w64-x86_64-autotools \
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

mkdir _deps && cd _deps

# Fetch older version of cfitsio (4.3.1-1)
wget https://repo.msys2.org/mingw/mingw64/mingw-w64-x86_64-cfitsio-1~4.3.1-1-any.pkg.tar.zst
pacman --noconfirm -U mingw-w64-x86_64-cfitsio-1~4.3.1-1-any.pkg.tar.zst
cd ..

# Build LibRaw from github
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
autoreconf -fi && \
./configure --disable-examples --disable-static && \
make install -j$(nproc) || exit 1
cd ..

