# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This file contains everything related to Siril.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced
# shellcheck disable=SC2034 # no exports desired

### variables ##################################################################

#----------------------------------------------- source directory and git branch

# If we're running inside Siril's official CI, the repository is already
# there and we adjust SIRIL_DIR accordingly.
# If not, check if a custom repository location and/or branch has been
# specified in the environment.

if $CI_GITLAB; then   # running GitLab CI
  SIRIL_DIR=$SELF_DIR/../..
else                  # not running GitLab CI
  SIRIL_DIR=$SRC_DIR/siril

  # Allow using a custom Siril repository and branch.
  if [ -z "$SIRIL_URL" ]; then
    SIRIL_URL=https://gitlab.com/free-astro/siril
  fi

  # Allow using a custom branch.
  if [ -z "$SIRIL_BRANCH" ]; then
    SIRIL_BRANCH=master
  fi
fi

SIRIL_BLD_DIR=$BLD_DIR/$(basename "$SIRIL_DIR")  # we build out-of-tree

#------------------------------------ Python runtime to be bundled with Siril

# Siril will be bundled with its own (customized) Python 3 runtime to make
# the core extensions work out-of-the-box. This is independent from whatever
# Python is running JHBuild or getting built as a dependency.
#
# We are only pinning major and minor versions here, not the patch level.
# Patch level is determined by whatever is current in the python_macos
# project.

SIRIL_PYTHON_VER_MAJOR=3
SIRIL_PYTHON_VER_MINOR=8
SIRIL_PYTHON_VER=$SIRIL_PYTHON_VER_MAJOR.$SIRIL_PYTHON_VER_MINOR
SIRIL_PYTHON_URL="https://gitlab.com/dehesselle/python_macos/-/jobs/\
artifacts/master/raw/python_${SIRIL_PYTHON_VER//.}_$(uname -p).tar.xz?\
job=python${SIRIL_PYTHON_VER//.}:inkscape:$(uname -p)"

# Python packages are also built externally (on a system running the oldest
# supported OS for better backward compatiblity) and included here.

SIRIL_PYTHON_WHEELS_VER=0.51
SIRIL_PYTHON_WHEELS_URL=https://github.com/dehesselle/mibap_wheels/releases/\
download/v$SIRIL_PYTHON_WHEELS_VER/wheels.tar.xz

#----------------------------------- Python packages to be bundled with Siril

# https://pypi.org/project/cssselect/
SIRIL_PYTHON_CSSSELECT=cssselect==1.1.0

# https://pypi.org/project/lxml/
SIRIL_PYTHON_LXML=lxml==4.6.3

# https://pypi.org/project/numpy/
SIRIL_PYTHON_NUMPY=numpy==1.20.3

# https://pypi.org/project/PyGObject/
SIRIL_PYTHON_PYGOBJECT="\
  PyGObject==3.40.1\
  pycairo==1.20.0\
"

# https://pypi.org/project/pyserial/
SIRIL_PYTHON_PYSERIAL=pyserial==3.5

# https://pypi.org/project/scour/
SIRIL_PYTHON_SCOUR="\
  scour==0.38.2\
  six==1.16.0\
"

# https://pypi.org/project/urllib3
SIRIL_PYTHON_URLLIB3=urllib3==1.26.5

#------------------------------------------- application bundle directory layout

SIRIL_APP_DIR=$ARTIFACT_DIR/Siril.app

SIRIL_APP_CON_DIR=$SIRIL_APP_DIR/Contents
SIRIL_APP_RES_DIR=$SIRIL_APP_CON_DIR/Resources
SIRIL_APP_FRA_DIR=$SIRIL_APP_CON_DIR/Frameworks
SIRIL_APP_BIN_DIR=$SIRIL_APP_RES_DIR/bin
SIRIL_APP_ETC_DIR=$SIRIL_APP_RES_DIR/etc
SIRIL_APP_EXE_DIR=$SIRIL_APP_CON_DIR/MacOS
SIRIL_APP_LIB_DIR=$SIRIL_APP_RES_DIR/lib

SIRIL_APP_SITEPKG_DIR=$SIRIL_APP_LIB_DIR/python$SIRIL_PYTHON_VER/site-packages

### functions ##################################################################

function siril_get_version
{
  local file=$SIRIL_DIR/meson.build
  
  local ver_major
  ver_major=$(grep "version : " meson.build | head -n 1  | awk '{ print $3 }' | tr -d "'," | cut -d . -f 1)
  local ver_minor
  ver_minor=$(grep "version : " meson.build | head -n 1  | awk '{ print $3 }' | tr -d "'," | cut -d . -f 2)
  local ver_patch
  ver_patch=$(grep "version : " meson.build | head -n 1  | awk '{ print $3 }' | tr -d "'," | cut -d . -f 3)

  # shellcheck disable=SC2086 # they are integers
  echo $ver_major.$ver_minor.$ver_patch
}

function siril_get_repo_shorthash
{
  git -C "$SIRIL_DIR" rev-parse --short HEAD
}

function siril_pipinstall
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
  export PATH=$SIRIL_APP_FRA_DIR/Python.framework/Versions/Current/bin:$PATH

  # shellcheck disable=SC2086 # we need word splitting here
  pip$SIRIL_PYTHON_VER_MAJOR install \
    --prefix "$SIRIL_APP_RES_DIR" \
    --ignore-installed \
    $options \
    $wheels

  export PATH=$PATH_ORIGINAL
}

function siril_pipinstall_cssselect
{
  local wheels_dir=$1
  local options=$2

  siril_pipinstall "$SIRIL_PYTHON_CSSSELECT" "$wheels_dir" "$options"
}

function siril_pipinstall_lxml
{
  local wheels_dir=$1
  local options=$2

  siril_pipinstall "$SIRIL_PYTHON_LXML" "$wheels_dir" "$options"

  lib_change_paths \
    @loader_path/../../.. \
    "$SIRIL_APP_LIB_DIR" \
    "$SIRIL_APP_SITEPKG_DIR"/lxml/etree.cpython-"${SIRIL_PYTHON_VER/./}"-darwin.so \
    "$SIRIL_APP_SITEPKG_DIR"/lxml/objectify.cpython-"${SIRIL_PYTHON_VER/./}"-darwin.so
}

function siril_pipinstall_numpy
{
  local wheels_dir=$1
  local options=$2

  siril_pipinstall "$SIRIL_PYTHON_NUMPY" "$wheels_dir" "$options"

  sed -i '' '1s/.*/#!\/usr\/bin\/env python'"$SIRIL_PYTHON_VER_MAJOR"'/' \
    "$SIRIL_APP_BIN_DIR"/f2py
  sed -i '' '1s/.*/#!\/usr\/bin\/env python'"$SIRIL_PYTHON_VER_MAJOR"'/' \
    "$SIRIL_APP_BIN_DIR"/f2py$SIRIL_PYTHON_VER_MAJOR
  sed -i '' '1s/.*/#!\/usr\/bin\/env python'"$SIRIL_PYTHON_VER_MAJOR"'/' \
    "$SIRIL_APP_BIN_DIR"/f2py$SIRIL_PYTHON_VER
}

function siril_pipinstall_pygobject
{
  local wheels_dir=$1
  local options=$2

  siril_pipinstall "$SIRIL_PYTHON_PYGOBJECT" "$wheels_dir" "$options"

  lib_change_paths \
    @loader_path/../../.. \
    "$SIRIL_APP_LIB_DIR" \
    "$SIRIL_APP_SITEPKG_DIR"/gi/_gi.cpython-"${SIRIL_PYTHON_VER/./}"-darwin.so \
    "$SIRIL_APP_SITEPKG_DIR"/gi/_gi_cairo.cpython-"${SIRIL_PYTHON_VER/./}"-darwin.so \
    "$SIRIL_APP_SITEPKG_DIR"/cairo/_cairo.cpython-"${SIRIL_PYTHON_VER/./}"-darwin.so
}

function siril_pipinstall_pyserial
{
  local wheels_dir=$1
  local options=$2

  siril_pipinstall "$SIRIL_PYTHON_PYSERIAL" "$wheels_dir" "$options"

  find "$SIRIL_APP_SITEPKG_DIR"/serial -type f -name "*.pyc" -exec rm {} \;
  sed -i '' '1s/.*/#!\/usr\/bin\/env python3/' "$SIRIL_APP_BIN_DIR"/pyserial-miniterm
  sed -i '' '1s/.*/#!\/usr\/bin\/env python3/' "$SIRIL_APP_BIN_DIR"/pyserial-ports
}

function siril_pipinstall_scour
{
  local wheels_dir=$1
  local options=$2

  siril_pipinstall "$SIRIL_PYTHON_SCOUR" "$wheels_dir" "$options"

  sed -i '' '1s/.*/#!\/usr\/bin\/env python3/' "$SIRIL_APP_BIN_DIR"/scour
}

function siril_pipinstall_urllib3
{
  local wheels_dir=$1
  local options=$2

  siril_pipinstall "$SIRIL_PYTHON_URLLIB3" "$wheels_dir" "$options"
}

function siril_download_python
{
  curl -o "$PKG_DIR"/"$(basename "${SIRIL_PYTHON_URL%\?*}")" -L "$SIRIL_PYTHON_URL"
}

function siril_install_python
{
  mkdir "$SIRIL_APP_FRA_DIR"
  tar -C "$SIRIL_APP_FRA_DIR" -xf "$PKG_DIR"/"$(basename "${SIRIL_PYTHON_URL%\?*}")"

  # link it to SIRIL_APP_BIN_DIR so it'll be in PATH for the app
  mkdir -p "$SIRIL_APP_BIN_DIR"
  # shellcheck disable=SC2086 # it's an integer
  ln -sf ../../Frameworks/Python.framework/Versions/Current/bin/\
python$SIRIL_PYTHON_VER_MAJOR "$SIRIL_APP_BIN_DIR"

  # create '.pth' file inside Framework to include our site-packages directory
  # shellcheck disable=SC2086 # it's an integer
  echo "../../../../../../../Resources/lib/python$SIRIL_PYTHON_VER/site-packages"\
    > "$SIRIL_APP_FRA_DIR"/Python.framework/Versions/Current/lib/\
python$SIRIL_PYTHON_VER/site-packages/siril.pth
}

# shellcheck disable=SC2086 # we need word splitting here
function siril_build_wheels
{
  jhbuild run pip3 install wheel
  jhbuild run pip3 wheel $SIRIL_PYTHON_CSSSELECT -w "$PKG_DIR"
  jhbuild run pip3 wheel --no-binary :all: $SIRIL_PYTHON_LXML -w "$PKG_DIR"
  jhbuild run pip3 wheel $SIRIL_PYTHON_NUMPY     -w "$PKG_DIR"
  jhbuild run pip3 wheel $SIRIL_PYTHON_PYGOBJECT -w "$PKG_DIR"
  jhbuild run pip3 wheel $SIRIL_PYTHON_PYSERIAL  -w "$PKG_DIR"
  jhbuild run pip3 wheel $SIRIL_PYTHON_SCOUR     -w "$PKG_DIR"
  jhbuild run pip3 wheel $SIRIL_PYTHON_URLLIB3   -w "$PKG_DIR"
}

function siril_download_wheels
{
  curl \
    -o "$PKG_DIR"/"$(basename "${SIRIL_PYTHON_WHEELS_URL%\?*}")" \
    -L "$SIRIL_PYTHON_WHEELS_URL"
}
