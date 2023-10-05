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
mingw-w64-x86_64-gsl \
mingw-w64-x86_64-opencv \
mingw-w64-x86_64-libheif \
mingw-w64-x86_64-ffms2 \
mingw-w64-x86_64-cfitsio \
mingw-w64-x86_64-libgit2-winhttp \

# Build LibRaw from github
mkdir _deps && cd _deps
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
autoreconf -fi && \
./configure --disable-examples --disable-static && \
make install -j$(nproc) || exit 1
cd ..

# Build libXISF from git rep
git clone https://gitea.nouspiro.space/nou/libXISF.git
cd libXISF
mkdir -p build && cd build
cmake -G "MSYS Makefiles" -DCMAKE_INSTALL_PREFIX="$MSYSTEM_PREFIX" -DCMAKE_BUILD_TYPE="Release" ..
make install || exit 1
cd ..
