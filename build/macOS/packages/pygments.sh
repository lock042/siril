# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This file contains everything related to setup Pygments.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced

### variables ##################################################################

# https://pygments.org
PYGMENTS_PIP=pygments==2.8.1

### functions ##################################################################

function pygments_install
{
  pip3 install --prefix "$VER_DIR" "$PYGMENTS_PIP"
}
