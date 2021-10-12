# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This file contains everything related to setup Ninja build system.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced

### variables ##################################################################

# https://github.com/ninja-build/ninja
NINJA_VER=1.8.2
NINJA_URL=https://github.com/ninja-build/ninja/releases/download/v$NINJA_VER/\
ninja-mac.zip

if [ -z "$NINJA_OPTS" ]; then
  NINJA_OPTS=   # yes, empty by default
fi

### functions ##################################################################

function ninja_install
{
  if [ -f "$BIN_DIR"/ninja ]; then
    echo_i "already installed"
  else
    curl -o "$PKG_DIR"/"$(basename $NINJA_URL)" -L "$NINJA_URL"
    unzip -d "$BIN_DIR" "$PKG_DIR/$(basename "$NINJA_URL")"
  fi
}

function ninja_run
{
  local args="$*"

  # shellcheck disable=SC2086 # we want word splitting
  ninja $NINJA_OPTS $args
}