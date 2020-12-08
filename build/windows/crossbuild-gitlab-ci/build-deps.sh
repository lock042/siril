crossroad source msys2
mkdir _deps && cd _deps

crossroad install fftw exiv2 gtk3 libconfig gsl opencv libheif ffms2 zlib cfitsio

wget https://indilib.org/jdownloads/wcslib/wcslib-7.3.1.tar.gz
tar xfv wcslib-7.3.1.tar.gz && cd wcslib-7.3.1
mkdir build && cd build 
crossroad cmake .. -DCMAKE_BUILD_TYPE=Release
make
make install || exit 1
