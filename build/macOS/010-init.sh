# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This is the main initialization file to setup the environment. Its purpose is
#   - to provide some basic configuration variables which all other files
#     depend upon
#   - source all other scripts
#   - run a few essential checks to see if we're good
#
# It's meant to be sourced by all other scripts and supposed to be a "passive"
# file, i.e. it defines variables and functions but does not do anything on its
# own. However, this is only 99% true at the moment as the above mentioned
# checks are capable of calling it quits (search for 'exit') if a few very
# fundamental things appear to be broken.

### includes ###################################################################

# Nothing here.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced
# shellcheck disable=SC2034 # we only use export if we really need it

### main #######################################################################

#--------------------------------------------------------------- toolset version

VERSION=0.51

#-------------------------------------------------------------- target OS by SDK

# If SDKROOT is set, use that. If it is not set, try to select the 10.13 SDK
# (which is our minimum system requirement/target) and fallback to whatever
# SDK is available as the default one.

if [ -z "$SDKROOT" ]; then
  SDKROOT=$(xcodebuild -version -sdk macosx10.13 Path 2>/dev/null ||
            xcodebuild -version -sdk macosx Path)
fi
export SDKROOT

#--------------------------------------------------------------------- detect CI

if [ -z "$CI" ]; then   # Both GitHub and GitLab set this.
  CI=false
  CI_GITHUB=false
  CI_GITLAB=false
else
  CI=true

  if [ -z "$CI_PROJECT_NAME" ]; then  # This is a GitLab variable.
    CI_GITHUB=true
    CI_GITLAB=false
  else
    CI_GITHUB=false
    CI_GITLAB=true
  fi
fi

#------------------------------------------------------------- directories: work

# This is the main directory where all the action takes place below. The
# default, being directly below /Users/Shared, is guaranteed user-writable
# and present on every macOS system.

if [ -z "$WRK_DIR" ]; then
  WRK_DIR=/Users/Shared/work
fi

#---------------------------- directories: FSH-like layout for the build toolset

VER_DIR=$WRK_DIR/$VERSION
BIN_DIR=$VER_DIR/bin
ETC_DIR=$VER_DIR/etc
INC_DIR=$VER_DIR/include
LIB_DIR=$VER_DIR/lib
OPT_DIR=$VER_DIR/opt
VAR_DIR=$VER_DIR/var
BLD_DIR=$VAR_DIR/build
PKG_DIR=$VAR_DIR/cache/pkgs
SRC_DIR=$VER_DIR/usr/src
TMP_DIR=$VER_DIR/tmp

export HOME=$VER_DIR/home   # Yes, we redirect the user's home.

#---------------------------------------------- directories: temporary locations

export TMP=$TMP_DIR
export TEMP=$TMP_DIR
export TMPDIR=$TMP_DIR   # TMPDIR is the common macOS default.

#-------------------------------------------------------------- directories: XDG

export XDG_CACHE_HOME=$VAR_DIR/cache  # instead ~/.cache
export XDG_CONFIG_HOME=$ETC_DIR       # instead ~/.config

#-------------------------------------------------------------- directories: pip

export PIP_CACHE_DIR=$XDG_CACHE_HOME/pip       # instead ~/Library/Caches/pip
export PIPENV_CACHE_DIR=$XDG_CACHE_HOME/pipenv # instead ~/Library/Caches/pipenv

#--------------------------------------------------------- directories: artifact

# In CI mode, the artifacts are placed into the respective project repositories
# so they can be picked up from there. In non-CI mode the artifacts are
# placed in VER_DIR.

if   $CI_GITHUB; then
  ARTIFACT_DIR=$GITHUB_WORKSPACE
elif $CI_GITLAB; then
  ARTIFACT_DIR=$CI_PROJECT_DIR
else
  ARTIFACT_DIR=$VER_DIR
fi

#------------------------------------------------------ directories: set up path

export PATH=$BIN_DIR:/usr/bin:/bin:/usr/sbin:/sbin

#------------------------------------------------------------- directories: self

# We want a fully qualified path to our directory in canonicalized form. Since
# we neither have 'readlink -f' nor 'realpath' on macOS, we use Python.
# There is a fallback that solely relies on BASH_SOURCE for systems that
# do not provide python3 (we will provide our own via 110-sysprep.sh).

SELF_DIR=$(dirname \
  "$(python3 -c "import os; print(os.path.realpath('${BASH_SOURCE[0]}'))" \
    2>/dev/null || echo "${BASH_SOURCE[0]}")")

#------------------------------------------ source functions from bash_d library

# Things I do not want to re-invent and/or am sharing between other projects
# of mine come from bash_d.
# https://github.com/dehesselle/bash_d

INCLUDE_DIR=$SELF_DIR/bash_d
# shellcheck source=bash_d/1_include.sh
source "$INCLUDE_DIR"/1_include.sh
include_file echo.sh
include_file error.sh
include_file lib.sh
include_file sed.sh

#----------------------------------------------------------- source our packages

for package in "$SELF_DIR"/packages/*.sh; do
  # shellcheck disable=SC1090 # can't point to a single source here
  source "$package"
done

# shellcheck disable=SC2034 # this is from ansi_.sh
ANSI_TERM_ONLY=false   # use ANSI control characters even if not in terminal

#---------------------------------------------------- check if WRK_DIR is usable

# shellcheck disable=SC2046 # result is integer
if  [ $(mkdir -p "$WRK_DIR" 2>/dev/null; echo $?) -eq 0 ] &&
    [ -w "$WRK_DIR" ] ; then
  : # WRK_DIR has been created or was already there and is writable
else
  echo_e "WRK_DIR not usable: $WRK_DIR"
  exit 1
fi

#----------------------------------------------------- check for presence of SDK

if [ ! -d "$SDKROOT" ]; then
  echo_e "SDK not found: $SDKROOT"
  exit 1
fi

#---------------------------------------------------- check recommended versions

sys_check_versions
