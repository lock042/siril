# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This file contains everything related to setup Meson build system.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced

### variables ##################################################################

# https://mesonbuild.com
MESON_PIP=meson==0.57.1

### functions ##################################################################

function meson_install
{
  pip3 install --prefix "$VER_DIR" "$MESON_PIP"
}
