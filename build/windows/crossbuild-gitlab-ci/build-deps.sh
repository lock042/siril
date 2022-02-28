crossroad source msys2

mkdir _deps && cd _deps

# Build LibRaw from github
crossroad install lcms2
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
autoreconf -fi && \
crossroad ./configure --disable-examples --disable-static && \
make install || exit 1
cd ..

# Install librtprocess from here
cd ..
mkdir subprojects/librtprocess/_build && cd subprojects/librtprocess/_build
crossroad cmake -G Ninja -DCMAKE_BUILD_TYPE="Release" -DBUILD_SHARED_LIBS=OFF .. && ninja && ninja install
cd ../../..

cd _deps
# Install deps from crossroad
crossroad install fftw \
                  exiv2 \
                  gtk3 \
                  libconfig \
                  gsl \
                  opencv \
                  libheif \
                  ffms2 \
                  cfitsio \



if [ $? -ne 0 ]; then
  echo "Installation of pre-built dependencies failed.";
  exit 1;
fi