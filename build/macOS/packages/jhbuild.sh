# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This file contains everything related to setup JHBuild.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced

### variables ##################################################################

#----------------------------------------------------------------------- JHBuild

export JHBUILDRC=$ETC_DIR/jhbuildrc
export JHBUILDRC_CUSTOM=$JHBUILDRC-custom

# JHBuild build system (3.38.0+ from master branch because of specific patch)
# https://gitlab.gnome.org/GNOME/jhbuild
# https://wiki.gnome.org/Projects/Jhbuild/Introduction
JHBUILD_VER=a896cbf404461cab979fa3cd1c83ddf158efe83b
JHBUILD_URL=https://gitlab.gnome.org/GNOME/jhbuild/-/archive/$JHBUILD_VER/\
jhbuild-$JHBUILD_VER.tar.bz2

# This Python runtime powers JHBuild on system that do not provide Python 3.
JHBUILD_PYTHON_VER_MAJOR=3
JHBUILD_PYTHON_VER_MINOR=8
JHBUILD_PYTHON_VER=$JHBUILD_PYTHON_VER_MAJOR.$JHBUILD_PYTHON_VER_MINOR
JHBUILD_PYTHON_URL="https://gitlab.com/dehesselle/python_macos/-/jobs/\
artifacts/master/raw/python_${JHBUILD_PYTHON_VER//.}_$(uname -p).tar.xz?\
job=python${JHBUILD_PYTHON_VER//.}:$(uname -p)"

### functions ##################################################################

function jhbuild_install_python
{
  mkdir -p "$OPT_DIR"
  curl -L "$JHBUILD_PYTHON_URL" | tar -C "$OPT_DIR" -x

  local python_bin_dir
  python_bin_dir=$OPT_DIR/Python.framework/Versions/Current/bin
  ln -s "$python_bin_dir"/python$JHBUILD_PYTHON_VER_MAJOR "$BIN_DIR"
  ln -s "$python_bin_dir"/python$JHBUILD_PYTHON_VER "$BIN_DIR"
  ln -s "$python_bin_dir"/pip$JHBUILD_PYTHON_VER_MAJOR "$BIN_DIR"
}

function jhbuild_install
{
  # The safest way is to use our custom Python, even if the system provides one.
  jhbuild_install_python

  # Without this, JHBuild won't be able to access https links later because
  # Apple's Python won't be able to validate certificates.
  certifi_install

  # Download JHBuild.
  local archive
  archive=$PKG_DIR/$(basename $JHBUILD_URL)
  curl -o "$archive" -L "$JHBUILD_URL"
  tar -C "$SRC_DIR" -xf "$archive"
  JHBUILD_DIR=$SRC_DIR/jhbuild-$JHBUILD_VER

  # Determine the real path to the Python interpreter so we can hardcode it into
  # the file. This way we can pick up the desired interpreter via environment
  # lookup once - i.e. right now - and stick to it.
  local python_interpreter
  python_interpreter=$(python3 -c \
    "import os, sys; print(os.path.realpath(sys.executable))")

  # Create 'jhbuild' executable. This code has been adapted from
  # https://gitlab.gnome.org/GNOME/gtk-osx/-/blob/master/gtk-osx-setup.sh
  cat <<EOF > "$BIN_DIR/jhbuild"
#!$python_interpreter
# -*- coding: utf-8 -*-

import sys
import os
import builtins

sys.path.insert(0, '$JHBUILD_DIR')
pkgdatadir = None
datadir = None
import jhbuild
srcdir = os.path.abspath(os.path.join(os.path.dirname(jhbuild.__file__), '..'))

builtins.__dict__['PKGDATADIR'] = pkgdatadir
builtins.__dict__['DATADIR'] = datadir
builtins.__dict__['SRCDIR'] = srcdir

import jhbuild.main
jhbuild.main.main(sys.argv[1:])
EOF

  chmod 755 "$BIN_DIR"/jhbuild

  # Install JHBuild's external dependencies.
  meson_install
  ninja_install
}

function jhbuild_configure
{
  # Copy JHBuild configuration files. We can't use 'cp' because that doesn't
  # with the ramdisk overlay.
  mkdir -p "$(dirname "$JHBUILDRC")"
  # shellcheck disable=SC2094 # not the same file
  cat "$SELF_DIR"/jhbuild/"$(basename "$JHBUILDRC")"        > "$JHBUILDRC"
  # shellcheck disable=SC2094 # not the same file
  cat "$SELF_DIR"/jhbuild/"$(basename "$JHBUILDRC_CUSTOM")" > "$JHBUILDRC_CUSTOM"

  {
    # set moduleset directory
    echo "modulesets_dir = '$SELF_DIR/jhbuild'"

    # basic directory layout
    echo "buildroot = '$BLD_DIR'"
    echo "checkoutroot = '$SRC_DIR'"
    echo "prefix = '$VER_DIR'"
    echo "tarballdir = '$PKG_DIR'"
    echo "top_builddir = '$VAR_DIR/jhbuild'"

    # set macOS SDK
    echo "setup_sdk(sdkdir=\"$SDKROOT\")"

    # set release build
    echo "setup_release()"

    # enable ccache
    echo "os.environ[\"CC\"] = \"$BIN_DIR/gcc\""
    echo "os.environ[\"OBJC\"] = \"$BIN_DIR/gcc\""
    echo "os.environ[\"CXX\"] = \"$BIN_DIR/g++\""

    # certificates for https
    echo "os.environ[\"SSL_CERT_FILE\"] = \
      \"$LIB_DIR/python$JHBUILD_PYTHON_VER/site-packages/certifi/cacert.pem\""
    echo "os.environ[\"REQUESTS_CA_BUNDLE\"] = \
      \"$LIB_DIR/python$JHBUILD_PYTHON_VER/site-packages/certifi/cacert.pem\""

    # user home directory
    echo "os.environ[\"HOME\"] = \"$HOME\""

    # less noise on the terminal if not CI
    if ! $CI; then
      echo "quiet_mode = True"
      echo "progress_bar = True"
    fi

  } >> "$JHBUILDRC_CUSTOM"
}
