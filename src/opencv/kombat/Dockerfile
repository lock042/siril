# This file allows former contributor of KOMBAT to get an isolated
# compilation/test environment for the algorithm.

# To use this:
#
#  a) set $my_dirs variable according to directories you want to be accessible from your Image,
#      my_dirs="-v /home/fred/git/siril:/src -v /media:/media -v /home/fred/git/kombat/:/kombat"
#
#  b) Build/update the image (if needed):
#      docker build --no-cache . -t kombat-dev-img
#
#  c) Run it:
#     docker run -it --rm --memory 8G -e DISPLAY=${DISPLAY} -v /tmp/.X11-unix:/tmp/.X11-unix --net host $my_dirs kombat-dev-img /bin/bash
#     cd /src                                                     # ofc, only if you set up Siril sources from host in /src,
#     make && make ofc, only install && siril # or make -j8 if you want parallel compilation.

FROM ubuntu

WORKDIR /tmp
RUN apt-get update

RUN ln -fs /usr/share/zoneinfo/Europe/Paris /etc/localtime
RUN apt-get install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get install -y g++-10
RUN apt-get install -y intltool libglib2.0-dev libjson-glib-dev libgtk-3-dev fftw3-dev libgsl-dev libcfitsio-dev libconfig-dev libtiff-dev libavutil-dev libswscale-dev libswresample-dev libavformat-dev libexiv2-dev

RUN apt-get install -y build-essential checkinstall cmake pkg-config yasm

RUN apt-get install -y wget git unzip

RUN apt-get install -y wcslib-dev

WORKDIR /tmp
RUN wget -O opencv.tar.gz https://github.com/opencv/opencv/archive/4.5.5/opencv-4.5.5.tar.gz
# RUN wget -O opencv-contrib.tar.gz https://github.com/opencv/opencv_contrib/archive/4.5.5/opencv_contrib-4.5.5.tar.gz
RUN tar xzvf opencv.tar.gz && rm opencv.tar.gz
# RUN unzip opencv-contrib.zip && rm opencv-contrib.zip

RUN cd opencv-4.5.5 && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=/usr      \
      -DCMAKE_BUILD_TYPE=Release       \
      -DENABLE_CXX11=ON                \
      -DBUILD_PERF_TESTS=OFF           \
      -DOPENCV_GENERATE_PKGCONFIG=YES  \
      -DWITH_XINE=ON                   \
      -DBUILD_TESTS=OFF                \
      -DENABLE_PRECOMPILED_HEADERS=OFF \
      -DCMAKE_SKIP_RPATH=ON            \
      -DBUILD_WITH_DEBUG_INFO=OFF      \
      -Wno-dev  ..                     && \
	make -j 4 && make install
