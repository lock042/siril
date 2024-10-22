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
                  curl \
                  librsvg \
                  libraw \
# need to uninstall crt-git
# probably same root cause as https://github.com/msys2/MINGW-packages/issues/10837
# otherwise, it's messing up all the subsequent builds 
crossroad uninstall crt-git

