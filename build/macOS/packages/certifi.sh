# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# Install Mozilla root certificates to facilitate SSL certificate checks.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced

### variables ##################################################################

# https://pypi.org/project/certifi/
CERTIFI_PIP=certifi   # unversioned on purpose

### functions ##################################################################

function certifi_install
{
  pip3 install --prefix "$VER_DIR" "$CERTIFI_PIP"
}
