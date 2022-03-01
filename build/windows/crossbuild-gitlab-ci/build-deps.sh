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

# disable LibRaw for now
# Build LibRaw from github
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
autoreconf -fi && \
crossroad ./configure --disable-examples && \
make install || exit 1
cd ..


if [ $? -ne 0 ]; then
  echo "Installation of pre-built dependencies failed.";
  exit 1;
fi