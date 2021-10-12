# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# Convert svg to icns.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced

### variables ##################################################################

# Nothing here.

### functions ##################################################################

function svg2icns_install
{
  cairosvg_install
  png2icns_install
}

function svg2icns
{
  local svg_file=$1
  local icns_file=$2

  local png_file
  png_file=$TMP_DIR/$(basename -s .svg "$svg_file").png

  # svg to png
  jhbuild run cairosvg -f png -s 1 -o "$png_file" "$svg_file"

  # png to icns
  cd "$TMP_DIR" || exit 1   # png2icns.sh outputs to current directory
  png2icns.sh "$png_file"

  mv "$(basename -s .png "$png_file")".icns "$icns_file"
}
