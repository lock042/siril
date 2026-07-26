CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -I../../src/io/healpix

SRC = ../../src/io/healpix/xp_continuous.cpp ../../src/io/healpix/gaia_xp_design.c
HDR = ../../src/io/healpix/xp_continuous.h ../../src/io/healpix/gaia_xp_design.h

xpcts_validate: xpcts_validate.cpp $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -o $@ xpcts_validate.cpp $(SRC)

clean:
	rm -f xpcts_validate
.PHONY: clean
