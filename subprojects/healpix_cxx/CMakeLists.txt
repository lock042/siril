## This project has been copied from
## https://sourceforge.net/p/healpix/code/HEAD/tree/branches/branch_v383r1259/src/cxx/
## and was modified to compile on its own, providing the necessary functions
## mainly to query healpixels

project(healpix_cxx CXX)
cmake_minimum_required(VERSION 3.10)

#LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")


# add sources
file(GLOB SRC_FILES *.cc)
#build a shared library
ADD_LIBRARY(healpix_cxx STATIC ${SRC_FILES})
set_target_properties(healpix_cxx PROPERTIES VERSION 8.3.0 SOVERSION 4)

# Glob headers
file(GLOB SRC_HEADERS *.h)
#install headers
INSTALL(FILES ${SRC_HEADERS} DESTINATION include/healpix_cxx)

#install library
INSTALL(TARGETS healpix_cxx LIBRARY DESTINATION lib
                            ARCHIVE DESTINATION lib)