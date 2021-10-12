#!/usr/bin/env bash
#
# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# Build a working version of Inkscape and install it into our tree (install
# prefix is VER_DIR). Also patch Inkscape's main binary and library as we do
# not have a valid @rpath at this point.

### includes ###################################################################

# shellcheck disable=SC1090 # can't point to a single source here
for script in "$(dirname "${BASH_SOURCE[0]}")"/0??-*.sh; do
  source "$script";
done

### settings ###################################################################

# shellcheck disable=SC2034 # this is from ansi_.sh
ANSI_TERM_ONLY=false   # use ANSI control characters even if not in terminal

error_trace_enable

### main #######################################################################

#-------------------------------------------------------- (re-) configure ccache

# This is required when using the precompiled toolset as ccache will not have
# been setup before (it happens in 110-sysprep.sh).

ccache_configure

#------------------------------------------------------- (re-) configure JHBuild

# This allows compiling Inkscape with a different setup than what the toolset
# was built with, most importantly a different SDK.

jhbuild_configure

#---------------------------------------------------------------- build Inkscape

if ! $CI_GITLAB; then     # not running GitLab CI

  if [ -d "$SIRIL_DIR"/.git ]; then   # Already cloned?
    # Special treatment 1/2 for local builds: Leave it up to the
    # user if they need a fresh clone. This way we enable making changes
    # to the code and running the build again.
    echo_w "using existing repository in $SIRIL_DIR"
    echo_w "Do 'rm -rf $SIRIL_DIR' if you want a fresh clone."
    sleep 5
  else
    git clone \
      --branch "$SIRIL_BRANCH" \
      --depth 10 \
      --recurse-submodules \
      --single-branch \
      "$SIRIL_URL" "$SIRIL_DIR"
  fi

  if ! $CI && [ -d "$SIRIL_BLD_DIR" ]; then   # Running locally and build exists?
    # Special treatment 2/2 for local builds: remove the build directory
    # to ensure clean builds.
    rm -rf "$SIRIL_BLD_DIR"
  fi
fi

mkdir -p "$SIRIL_BLD_DIR"

ninja -C "$SIRIL_BLD_DIR"

ninja_run
ninja_run install
ninja_run tests

#--------------------------------------------- make library link paths canonical

# Most libraries are linked to with their fully qualified paths. Some of them
# have been linked to using '@rpath', which we are not setting, therefore
# breaking the binary. Since we want Inkscape to work in unpackaged form as
# well, we adjust all offending paths to use qualified paths.
# (They will be re-adjusted later to use relative location inside the
# application bundle.)
#
# example 'ldd /Users/Shared/work/0.47/lib/inkscape/libinkscape_base.dylib':
#
#   /Users/Shared/work/0.47/lib/inkscape/libinkscape_base.dylib:
#     @rpath/libinkscape_base.dylib (compatibility version 0.0.0, ...     <- fix
#     @rpath/libboost_filesystem.dylib (compatibility version 0.0....     <- fix
#     @rpath/lib2geom.1.1.0.dylib (compatibility version 1.1.0, cu...     <- fix
#     /Users/Shared/work/0.47/lib/libharfbuzz.0.dylib (compatibili...     <- ok
#     /Users/Shared/work/0.47/lib/libpangocairo-1.0.0.dylib (compa...     <- ok
#     /Users/Shared/work/0.47/lib/libcairo.2.dylib (compatibility ...     <- ok
#     ...

#for binary in $BIN_DIR/inkscape \
#              $LIB_DIR/inkscape/libinkscape_base.dylib; do
#  for lib in $(otool -L "$binary" |
#               grep "@rpath/" |
#               awk '{ print $1 }'); do
#    # This here is the reason we require GNU's 'find', as the macOS one doesn't
#    # pick up the files from the bottom layer of our union mount.
#    lib_canonical=$(find "$LIB_DIR" -maxdepth 2 -name "$(basename "$lib")")
#    lib_change_path "$lib_canonical" "$binary"
#  done
#done
