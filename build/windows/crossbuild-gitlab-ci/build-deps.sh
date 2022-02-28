crossroad source msys2

mkdir _deps && cd _deps
# Install deps from crossroad
crossroad install fftw \
                  exiv2 \
                  libconfig \
                  gsl \
                  opencv \
                  libheif \
                  ffms2 \
                  cfitsio \
                  lcms2

# Build LibRaw from github
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
autoreconf -fi && \
crossroad ./configure --disable-examples --disable-static && \
make install || exit 1
cd ..

cd ..
# Install librtprocess from here
mkdir subprojects/librtprocess/_build && cd subprojects/librtprocess/_build
crossroad cmake -G Ninja -DCMAKE_BUILD_TYPE="Release" -DBUILD_SHARED_LIBS=OFF .. && ninja && ninja install
cd ../../..

if [ $? -ne 0 ]; then
  echo "Installation of pre-built dependencies failed.";
  exit 1;
fi