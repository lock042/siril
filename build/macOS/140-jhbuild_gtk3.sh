#!/usr/bin/env bash
#
# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# Install GTK3 libraries (and their dependencies).

### includes ###################################################################

# shellcheck disable=SC1090 # can't point to a single source here
for script in "$(dirname "${BASH_SOURCE[0]}")"/0??-*.sh; do
  source "$script";
done

### settings ###################################################################

# Nothing here.

### main #######################################################################

jhbuild build python3 # avoid https://gitlab.gnome.org/GNOME/gtk-osx/-/issues/32

pygments_install  # required by gtk-doc

jhbuild build \
  meta-gtk-osx-bootstrap \
  meta-gtk-osx-gtk3 \
  gtkmm3
