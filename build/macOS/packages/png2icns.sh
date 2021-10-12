# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# Convert png to icns.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced

### variables ##################################################################

# https://github.com/bitboss-ca/png2icns
PNG2ICNS_VER=0.1
PNG2ICNS_URL=https://github.com/bitboss-ca/png2icns/archive/\
v$PNG2ICNS_VER.tar.gz

### functions ##################################################################

function png2icns_install
{
  local archive
  archive=$PKG_DIR/$(basename $PNG2ICNS_URL)
  curl -o "$archive" -L "$PNG2ICNS_URL"
  tar -C "$SRC_DIR" -xf "$archive"
  ln -s "$SRC_DIR"/png2icns-$PNG2ICNS_VER/png2icns.sh "$BIN_DIR"
}
