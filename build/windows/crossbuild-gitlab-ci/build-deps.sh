crossroad source msys2
mkdir _deps && cd _deps

crossroad install fftw \
                  exiv2 \
                  gtk3 \
                  libconfig \
                  gsl \
                  opencv \
                  libheif \
                  ffms2 \
                  cfitsio

wget ftp://ftp.atnf.csiro.au/pub/software/wcslib/wcslib.tar.bz2
tar xfv wcslib.tar.bz2 && cd wcslib-7.3.1
patch -p1 -i ../../build/windows/crossbuild-gitlab-ci/wcslib.patch
crossroad configure LIBS=\"-pthread -lcurl -lm\" --without-pgplot --disable-fortran --prefix=${INSTALL_PREFIX}
make
make install || exit 1
