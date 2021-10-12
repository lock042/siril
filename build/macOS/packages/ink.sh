# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This file contains everything related to Inkscape.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced
# shellcheck disable=SC2034 # no exports desired

### variables ##################################################################

#----------------------------------------------- source directory and git branch

# If we're running inside Inkscape's official CI, the repository is already
# there and we adjust INK_DIR accordingly.
# If not, check if a custom repository location and/or branch has been
# specified in the environment.

if $CI_GITLAB; then   # running GitLab CI
  INK_DIR=$SELF_DIR/../..
else                  # not running GitLab CI
  INK_DIR=$SRC_DIR/inkscape

  # Allow using a custom Inkscape repository and branch.
  if [ -z "$INK_URL" ]; then
    INK_URL=https://gitlab.com/inkscape/inkscape
  fi

  # Allow using a custom branch.
  if [ -z "$INK_BRANCH" ]; then
    INK_BRANCH=master
  fi
fi

INK_BLD_DIR=$BLD_DIR/$(basename "$INK_DIR")  # we build out-of-tree

#------------------------------------ Python runtime to be bundled with Inkscape

# Inkscape will be bundled with its own (customized) Python 3 runtime to make
# the core extensions work out-of-the-box. This is independent from whatever
# Python is running JHBuild or getting built as a dependency.
#
# We are only pinning major and minor versions here, not the patch level.
# Patch level is determined by whatever is current in the python_macos
# project.

INK_PYTHON_VER_MAJOR=3
INK_PYTHON_VER_MINOR=8
INK_PYTHON_VER=$INK_PYTHON_VER_MAJOR.$INK_PYTHON_VER_MINOR
INK_PYTHON_URL="https://gitlab.com/dehesselle/python_macos/-/jobs/\
artifacts/master/raw/python_${INK_PYTHON_VER//.}_$(uname -p).tar.xz?\
job=python${INK_PYTHON_VER//.}:inkscape:$(uname -p)"

# Python packages are also built externally (on a system running the oldest
# supported OS for better backward compatiblity) and included here.

INK_PYTHON_WHEELS_VER=0.51
INK_PYTHON_WHEELS_URL=https://github.com/dehesselle/mibap_wheels/releases/\
download/v$INK_PYTHON_WHEELS_VER/wheels.tar.xz

#----------------------------------- Python packages to be bundled with Inkscape

# https://pypi.org/project/cssselect/
INK_PYTHON_CSSSELECT=cssselect==1.1.0

# https://pypi.org/project/lxml/
INK_PYTHON_LXML=lxml==4.6.3

# https://pypi.org/project/numpy/
INK_PYTHON_NUMPY=numpy==1.20.3

# https://pypi.org/project/PyGObject/
INK_PYTHON_PYGOBJECT="\
  PyGObject==3.40.1\
  pycairo==1.20.0\
"

# https://pypi.org/project/pyserial/
INK_PYTHON_PYSERIAL=pyserial==3.5

# https://pypi.org/project/scour/
INK_PYTHON_SCOUR="\
  scour==0.38.2\
  six==1.16.0\
"

# https://pypi.org/project/urllib3
INK_PYTHON_URLLIB3=urllib3==1.26.5

#------------------------------------------- application bundle directory layout

INK_APP_DIR=$ARTIFACT_DIR/Inkscape.app

INK_APP_CON_DIR=$INK_APP_DIR/Contents
INK_APP_RES_DIR=$INK_APP_CON_DIR/Resources
INK_APP_FRA_DIR=$INK_APP_CON_DIR/Frameworks
INK_APP_BIN_DIR=$INK_APP_RES_DIR/bin
INK_APP_ETC_DIR=$INK_APP_RES_DIR/etc
INK_APP_EXE_DIR=$INK_APP_CON_DIR/MacOS
INK_APP_LIB_DIR=$INK_APP_RES_DIR/lib

INK_APP_SITEPKG_DIR=$INK_APP_LIB_DIR/python$INK_PYTHON_VER/site-packages

### functions ##################################################################

function ink_get_version
{
  local file=$INK_DIR/CMakeLists.txt
  local ver_major
  ver_major=$(grep INKSCAPE_VERSION_MAJOR "$file" | head -n 1 | awk '{ print $2+0 }')
  local ver_minor
  ver_minor=$(grep INKSCAPE_VERSION_MINOR "$file" | head -n 1 | awk '{ print $2+0 }')
  local ver_patch
  ver_patch=$(grep INKSCAPE_VERSION_PATCH "$file" | head -n 1 | awk '{ print $2+0 }')
  local ver_suffix
  ver_suffix=$(grep INKSCAPE_VERSION_SUFFIX "$file" | head -n 1 | awk '{ print $2 }')

  ver_suffix=${ver_suffix%\"*}   # remove "double quote and everything after" from end
  ver_suffix=${ver_suffix#\"}   # remove "double quote" from beginning

  # shellcheck disable=SC2086 # they are integers
  echo $ver_major.$ver_minor.$ver_patch"$ver_suffix"
}

function ink_get_repo_shorthash
{
  # do it the same way as in CMakeScripts/inkscape-version.cmake
  git -C "$INK_DIR" rev-parse --short HEAD
}

function ink_pipinstall
{
  local packages=$1
  local wheels_dir=$2   # optional
  local options=$3      # optional

  if [ -z "$wheels_dir" ]; then
    wheels_dir=$PKG_DIR
  fi

  # turn package names into filenames of our wheels
  local wheels
  for package in $packages; do
    wheels="$wheels $(eval echo "$wheels_dir"/"${package%==*}"*.whl)"
  done

  local PATH_ORIGINAL=$PATH
  export PATH=$INK_APP_FRA_DIR/Python.framework/Versions/Current/bin:$PATH

  # shellcheck disable=SC2086 # we need word splitting here
  pip$INK_PYTHON_VER_MAJOR install \
    --prefix "$INK_APP_RES_DIR" \
    --ignore-installed \
    $options \
    $wheels

  export PATH=$PATH_ORIGINAL
}

function ink_pipinstall_cssselect
{
  local wheels_dir=$1
  local options=$2

  ink_pipinstall "$INK_PYTHON_CSSSELECT" "$wheels_dir" "$options"
}

function ink_pipinstall_lxml
{
  local wheels_dir=$1
  local options=$2

  ink_pipinstall "$INK_PYTHON_LXML" "$wheels_dir" "$options"

  lib_change_paths \
    @loader_path/../../.. \
    "$INK_APP_LIB_DIR" \
    "$INK_APP_SITEPKG_DIR"/lxml/etree.cpython-"${INK_PYTHON_VER/./}"-darwin.so \
    "$INK_APP_SITEPKG_DIR"/lxml/objectify.cpython-"${INK_PYTHON_VER/./}"-darwin.so
}

function ink_pipinstall_numpy
{
  local wheels_dir=$1
  local options=$2

  ink_pipinstall "$INK_PYTHON_NUMPY" "$wheels_dir" "$options"

  sed -i '' '1s/.*/#!\/usr\/bin\/env python'"$INK_PYTHON_VER_MAJOR"'/' \
    "$INK_APP_BIN_DIR"/f2py
  sed -i '' '1s/.*/#!\/usr\/bin\/env python'"$INK_PYTHON_VER_MAJOR"'/' \
    "$INK_APP_BIN_DIR"/f2py$INK_PYTHON_VER_MAJOR
  sed -i '' '1s/.*/#!\/usr\/bin\/env python'"$INK_PYTHON_VER_MAJOR"'/' \
    "$INK_APP_BIN_DIR"/f2py$INK_PYTHON_VER
}

function ink_pipinstall_pygobject
{
  local wheels_dir=$1
  local options=$2

  ink_pipinstall "$INK_PYTHON_PYGOBJECT" "$wheels_dir" "$options"

  lib_change_paths \
    @loader_path/../../.. \
    "$INK_APP_LIB_DIR" \
    "$INK_APP_SITEPKG_DIR"/gi/_gi.cpython-"${INK_PYTHON_VER/./}"-darwin.so \
    "$INK_APP_SITEPKG_DIR"/gi/_gi_cairo.cpython-"${INK_PYTHON_VER/./}"-darwin.so \
    "$INK_APP_SITEPKG_DIR"/cairo/_cairo.cpython-"${INK_PYTHON_VER/./}"-darwin.so
}

function ink_pipinstall_pyserial
{
  local wheels_dir=$1
  local options=$2

  ink_pipinstall "$INK_PYTHON_PYSERIAL" "$wheels_dir" "$options"

  find "$INK_APP_SITEPKG_DIR"/serial -type f -name "*.pyc" -exec rm {} \;
  sed -i '' '1s/.*/#!\/usr\/bin\/env python3/' "$INK_APP_BIN_DIR"/pyserial-miniterm
  sed -i '' '1s/.*/#!\/usr\/bin\/env python3/' "$INK_APP_BIN_DIR"/pyserial-ports
}

function ink_pipinstall_scour
{
  local wheels_dir=$1
  local options=$2

  ink_pipinstall "$INK_PYTHON_SCOUR" "$wheels_dir" "$options"

  sed -i '' '1s/.*/#!\/usr\/bin\/env python3/' "$INK_APP_BIN_DIR"/scour
}

function ink_pipinstall_urllib3
{
  local wheels_dir=$1
  local options=$2

  ink_pipinstall "$INK_PYTHON_URLLIB3" "$wheels_dir" "$options"
}

function ink_download_python
{
  curl -o "$PKG_DIR"/"$(basename "${INK_PYTHON_URL%\?*}")" -L "$INK_PYTHON_URL"
}

function ink_install_python
{
  mkdir "$INK_APP_FRA_DIR"
  tar -C "$INK_APP_FRA_DIR" -xf "$PKG_DIR"/"$(basename "${INK_PYTHON_URL%\?*}")"

  # link it to INK_APP_BIN_DIR so it'll be in PATH for the app
  mkdir -p "$INK_APP_BIN_DIR"
  # shellcheck disable=SC2086 # it's an integer
  ln -sf ../../Frameworks/Python.framework/Versions/Current/bin/\
python$INK_PYTHON_VER_MAJOR "$INK_APP_BIN_DIR"

  # create '.pth' file inside Framework to include our site-packages directory
  # shellcheck disable=SC2086 # it's an integer
  echo "../../../../../../../Resources/lib/python$INK_PYTHON_VER/site-packages"\
    > "$INK_APP_FRA_DIR"/Python.framework/Versions/Current/lib/\
python$INK_PYTHON_VER/site-packages/inkscape.pth
}

# shellcheck disable=SC2086 # we need word splitting here
function ink_build_wheels
{
  jhbuild run pip3 install wheel
  jhbuild run pip3 wheel $INK_PYTHON_CSSSELECT -w "$PKG_DIR"
  jhbuild run pip3 wheel --no-binary :all: $INK_PYTHON_LXML -w "$PKG_DIR"
  jhbuild run pip3 wheel $INK_PYTHON_NUMPY     -w "$PKG_DIR"
  jhbuild run pip3 wheel $INK_PYTHON_PYGOBJECT -w "$PKG_DIR"
  jhbuild run pip3 wheel $INK_PYTHON_PYSERIAL  -w "$PKG_DIR"
  jhbuild run pip3 wheel $INK_PYTHON_SCOUR     -w "$PKG_DIR"
  jhbuild run pip3 wheel $INK_PYTHON_URLLIB3   -w "$PKG_DIR"
}

function ink_download_wheels
{
  curl \
    -o "$PKG_DIR"/"$(basename "${INK_PYTHON_WHEELS_URL%\?*}")" \
    -L "$INK_PYTHON_WHEELS_URL"
}
