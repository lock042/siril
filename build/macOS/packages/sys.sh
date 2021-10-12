# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# System related things like checking versions of macOS and Xcode.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced
# shellcheck disable=SC2034 # no exports desired

### variables ##################################################################

SYS_MACOS_VER=$(sw_vers -productVersion)
SYS_MACOS_VER_RECOMMENDED=10.15.7

SYS_SDK_VER=$(/usr/libexec/PlistBuddy -c "Print \
:DefaultProperties:MACOSX_DEPLOYMENT_TARGET" "$SDKROOT"/SDKSettings.plist)
SYS_SDK_VER_RECOMMENDED=10.13

SYS_XCODE_VER=$( (xcodebuild -version 2>/dev/null || echo "Xcode n/a") | grep Xcode | awk '{ print $2 }')
SYS_XCODE_VER_RECOMMENDED=12.4

### functions ##################################################################

function sys_check_versions
{
  # Check version recommendations.

  if [ "$SYS_SDK_VER" != "$SYS_SDK_VER_RECOMMENDED" ]; then
    echo_w "recommended      SDK version: $SYS_SDK_VER_RECOMMENDED"
    echo_w "       your      SDK version: $SYS_SDK_VER"
  fi

  if [ "$SYS_XCODE_VER" != "$SYS_XCODE_VER_RECOMMENDED" ]; then
    echo_w "recommended    Xcode version: $SYS_XCODE_VER_RECOMMENDED"
    echo_w "       your    Xcode version: $SYS_XCODE_VER"
  fi

  if [ "$SYS_MACOS_VER" != "$SYS_MACOS_VER_RECOMMENDED" ]; then
    echo_w "recommended    macOS version: $SYS_MACOS_VER_RECOMMENDED"
    echo_w "       your    macOS version: $SYS_MACOS_VER"
  fi
}

function sys_create_log
{
  # Create release.log file.

  mkdir -p "$VAR_DIR"/log

  for var in SYS_MACOS_VER SYS_SDK_VER SYS_XCODE_VER VERSION WRK_DIR; do
    echo "$var = $(eval echo \$$var)" >> "$VAR_DIR"/log/release.log
  done
}
