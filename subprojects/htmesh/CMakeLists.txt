# This project has been copied from
# https://invent.kde.org/education/kstars/-/tree/master/kstars/htmesh
# and was slightly modified to compile on its own
project(htmesh)
cmake_minimum_required(VERSION 3.5.0)

INCLUDE(CheckFunctionExists)
check_function_exists(drand48 HAVE_DRAND48)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config-htmesh.h.cmake ${CMAKE_CURRENT_SOURCE_DIR}/config-htmesh.h )

### HTMesh library
SET(HTMesh_LIB_SRC
    MeshIterator.cpp
    HtmRange.cpp
    HtmRangeIterator.cpp
    RangeConvex.cpp
    SkipList.cpp
    SkipListElement.cpp
    SpatialConstraint.cpp
#   SpatialDomain.cpp
    SpatialEdge.cpp
    SpatialException.cpp
    SpatialIndex.cpp
    SpatialVector.cpp
    HTMesh.cpp
    MeshBuffer.cpp)

IF (UNITY_BUILD)
    ENABLE_UNITY_BUILD(HTMesh HTMesh_LIB_SRC 10 cpp)
ENDIF ()

SET(HTMesh_LIB_SRC ${HTMesh_LIB_SRC}
   cc_aux.c)

IF (ANDROID)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF ()

add_library(htmesh STATIC ${HTMesh_LIB_SRC})

if (BUILD_PYKSTARS)
  set_target_properties(htmesh PROPERTIES POSITION_INDEPENDENT_CODE ON)
ENDIF ()

add_definitions(-D_EXPORTING)

set_target_properties(htmesh PROPERTIES
    VERSION 1.0.0
    SOVERSION 1)

# It is not necessary to install the static library as we just link against it.
#if (NOT ANDROID)
#    install(TARGETS htmesh ${KDE_INSTALL_TARGETS_DEFAULT_ARGS} )
#endif ()

# If you wish to compile the HTMesh perl wrapper, uncomment this, rebuild and copy the library into /usr/lib/, because we will use it as a shared object. See README in the perl wrapper directory (kstars/kstars/data/tools/HTMesh-*) for more details.
#set_property(TARGET htmesh PROPERTY POSITION_INDEPENDENT_CODE YES)

