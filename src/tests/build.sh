#!/bin/sh
CC=gcc
LD=gcc
CFLAGS="-Wall -I.. -g -fopenmp `pkg-config --cflags gtk+-3.0` `pkg-config --cflags cfitsio` `pkg-config --cflags gsl`"
LDFLAGS="-fopenmp -rdynamic `pkg-config --libs gtk+-3.0 cfitsio gsl` -lopencv_core -lopencv_imgproc -lopencv_calib3d -lm \
 `pkg-config --libs libavformat libavutil libavcodec libswscale libswresample libraw libtiff-4 libjpeg fftw3 libconfig libcurl libpng` `pkg-config --libs ffms2`"
SLDFLAGS="../../deps/kplot/libkplot.a"

set -x
# $CC $CFLAGS -c -o compare_fits.o compare_fits.c &&
# $CC $CFLAGS -c -o dummy.o dummy.c &&
# $LD $LDFLAGS -o compare_fits compare_fits.o dummy.o ../io/image_format_fits.o ../core/utils.o ../gui/progress_and_log.o
# 
# $CC $CFLAGS -c -o sorting.o sorting.c &&
# $CC $CFLAGS -DUSE_ALL_SORTING_ALGOS -c -o ../algos/sorting.o ../algos/sorting.c &&
# $LD $LDFLAGS -o sorting sorting.o ../algos/sorting.o

$CC $CFLAGS -c -o zone_compare.o zone_compare.c &&
#$CC $CFLAGS -c -o dummy.o dummy.c &&
#$LD $LDFLAGS -o zone_compare zone_compare.o dummy.o ../io/image_format_fits.o ../core/utils.o ../gui/progress_and_log.o ../algos/statistics.o ../algos/sorting.o ../algos/quantize.o

$LD $LDFLAGS -o zone_compare zone_compare.o ../libsiril.a $SLDFLAGS
