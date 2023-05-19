crossroad source msys2

mkdir _deps && cd _deps

crossroad install lcms2 \
                  gtk3 \
                  fftw \
                  exiv2 \
                  libconfig \
                  gsl \
                  opencv\
                  libheif \
                  ffms2 \
                  cfitsio \
# need to uninstall crt-git
# probably same root cause as https://github.com/msys2/MINGW-packages/issues/10837
# otherwise, it's messing up all the subsequent builds 
crossroad uninstall crt-git

# Build LibRaw from github
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
patch --binary -p1 < ../../build/windows/crossbuild-gitlab-ci/libraw.patch
autoreconf -fi && \
crossroad ./configure --disable-examples --disable-static && \
make install || exit 1
cd ..

# Build libXISF from git rep
git clone https://gitea.nouspiro.space/nou/libXISF.git
cd libXISF
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr ..
make install || exit 1
cd ../..

if [ $? -ne 0 ]; then
  echo "Installation of pre-built dependencies failed.";
  exit 1;
fi