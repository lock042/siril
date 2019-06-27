#!/bin/sh
CC=clang
LD=clang
CFLAGS="-Wall -I.. -g -fopenmp `pkg-config --cflags gtk+-3.0` `pkg-config --cflags cfitsio` `pkg-config --cflags gsl`"
LDFLAGS="-fopenmp -rdynamic `pkg-config --libs gtk+-3.0` `pkg-config --libs cfitsio` `pkg-config --libs gsl` -lm"

set -x
# $CC $CFLAGS -c -o compare_fits.o compare_fits.c &&
# $CC $CFLAGS -c -o dummy.o dummy.c &&
# $LD $LDFLAGS -o compare_fits compare_fits.o dummy.o ../io/image_format_fits.o ../core/utils.o ../gui/progress_and_log.o
# 
# $CC $CFLAGS -c -o sorting.o sorting.c &&
# $CC $CFLAGS -DUSE_ALL_SORTING_ALGOS -c -o ../algos/sorting.o ../algos/sorting.c &&
# $LD $LDFLAGS -o sorting sorting.o ../algos/sorting.o

$CC $CFLAGS -c -o zone_compare.o zone_compare.c &&
$CC $CFLAGS -c -o dummy.o dummy.c &&
$LD $LDFLAGS -o zone_compare zone_compare.o dummy.o ../io/image_format_fits.o ../core/utils.o ../gui/progress_and_log.o ../algos/statistics.o ../algos/sorting.o ../algos/quantize.o
