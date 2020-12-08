crossroad source msys2

crossroad install fftw exiv2 gtk3 libconfig gsl opencv libheif ffms2 zlib

wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.48.tar.gz
tar xfv cfitsio-3.48.tar.gz
mkdir cfitsio-3.48/build 
cd cfitsio-3.48
patch -p1 -i ../build/windows/crossbuild-gitlab-ci/cfitsio-cmake.patch 
cd build 
crossroad cmake .. -DUSE_PTHREADS=ON -DCMAKE_BUILD_TYPE=Release
make
make install

cd
wget https://indilib.org/jdownloads/wcslib/wcslib-7.3.1.tar.gz
tar xfv wcslib-7.3.1.tar.gz && cd wcslib-7.3.1
mkdir build && cd build 
crossroad cmake .. -DCMAKE_BUILD_TYPE=Release
make
make install

