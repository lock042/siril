# SPDX-FileCopyrightText: 2021 René de Hesselle <dehesselle@web.de>
#
# SPDX-License-Identifier: GPL-2.0-or-later

### description ################################################################

# This file contains everything related to the toolset.

### settings ###################################################################

# shellcheck shell=bash # no shebang as this file is intended to be sourced
# shellcheck disable=SC2034 # no exports desired

### variables ##################################################################

TOOLSET_VER=$VERSION

TOOLSET_VOLNAME=mibap_v$TOOLSET_VER

# A disk image containing a built version of the whole toolset.
# https://github.com/dehesselle/mibap
TOOLSET_URL=https://github.com/dehesselle/mibap/releases/download/\
v$TOOLSET_VER/$TOOLSET_VOLNAME.dmg

TOOLSET_REPO_DIR=$WRK_DIR/repo   # persistent storage for downloaded dmg

if [ -z "$TOOLSET_OVERLAY_FILE" ]; then
  TOOLSET_OVERLAY_FILE=ram
fi

TOOLSET_OVERLAY_SIZE=3   # writable overlay, unit in GiB

### functions ##################################################################

function toolset_install
{
  local toolset_dmg
  toolset_dmg=$TOOLSET_REPO_DIR/$(basename "$TOOLSET_URL")

  if [ -f "$toolset_dmg" ]; then
    echo_i "toolset found: $toolset_dmg"
  else
    # File not present on disk, we need to download.
    echo_i "downloading: $TOOLSET_URL"
    toolset_download
  fi

  toolset_mount "$toolset_dmg" "$VER_DIR"

  # Sadly, there are some limitations involved with union-mounting:
  #   - Files are not visible to macOS' versions 'ls' or 'find'.
  #     (The GNU versions do work though.)
  #   - You cannot write in a location without having written to its
  #     parent location first. That's why we need to pre-create all directories
  #     below.
  #
  # Shadow-mounting a dmg is not a feasible alternative due to its
  # bad write-performance.
  #
  # Prepare a script for mass-creating directories. We have to do this before
  # untion-mounting as macOS' 'find' won't see the directories anymore after.
  # (GNU's 'find' does)
  find "$VER_DIR" -type d ! -path "$BLD_DIR/*" ! -path "$SRC_DIR/*" \
      -exec echo "mkdir {}" \; > "$WRK_DIR"/create_dirs.sh
    sed -i "" "1d" "$WRK_DIR"/create_dirs.sh   # remove first line ("file exists")
  chmod 755 "$WRK_DIR"/create_dirs.sh

  # create writable overlay
  if [ "$TOOLSET_OVERLAY_FILE" = "ram" ]; then    # overlay in memory
    device=$(hdiutil attach -nomount \
      ram://$((TOOLSET_OVERLAY_SIZE * 1024 * 2048)) | xargs)
    newfs_hfs -v "overlay" "$device" >/dev/null
    echo_i "$TOOLSET_OVERLAY_SIZE GiB ram attached to $device"
  else                                            # overlay on disk
    hdiutil create -size 3g -fs HFS+ -nospotlight \
      -volname overlay "$TOOLSET_OVERLAY_FILE"
    echo_i "$TOOLSET_OVERLAY_SIZE GiB sparseimage attached to $device"
    device=$(hdiutil attach -nomount "$TOOLSET_OVERLAY_FILE" |
      grep "^/dev/disk" | grep "Apple_HFS" | awk '{ print $1 }')
  fi

  mount -o nobrowse,rw,union -t hfs "$device" "$VER_DIR"
  echo_i "$device mounted at $VER_DIR"

  # create all directories inside overlay
  "$WRK_DIR"/create_dirs.sh
  rm "$WRK_DIR"/create_dirs.sh
}

function toolset_uninstall
{
  toolset_unmount "$VER_DIR"

  if [ -f "$TOOLSET_OVERLAY_FILE" ]; then
    rm "$TOOLSET_OVERLAY_FILE"
  fi
}

function toolset_download
{
  if [ ! -d "$TOOLSET_REPO_DIR" ]; then
    mkdir -p "$TOOLSET_REPO_DIR"
  fi

  curl -o "$TOOLSET_REPO_DIR"/"$(basename "$TOOLSET_URL")" -L "$TOOLSET_URL"
}

function toolset_mount
{
  local toolset_dmg=$1
  local mountpoint=$2

  echo_i "mounting compressed disk image may take some time..."

  if [ ! -d "$VER_DIR" ]; then
    mkdir -p "$VER_DIR"
  fi

  if [ -z "$mountpoint" ]; then
    hdiutil attach "$toolset_dmg"
  else
    local device
    device=$(hdiutil attach -nomount "$toolset_dmg" | grep "^/dev/disk" |
      grep "Apple_HFS" | awk '{ print $1 }')
    echo_i "toolset attached to $device"
    mount -o nobrowse,noowners,ro -t hfs "$device" "$mountpoint"
    echo_i "$device mounted at $VER_DIR"
  fi
}

function toolset_unmount
{
  local mountpoint=$1

  while : ; do   # unmount everything (in reverse order)
    local disk
    disk=$(mount | grep "$mountpoint" | tail -n1 | awk '{ print $1 }')
    disk=${disk%s[0-9]}  # cut off slice specification

    if [ ${#disk} -eq 0 ]; then
      break   # nothing to do here
    else
      # We loop over the 'eject' since it occasionally timeouts.
      while ! diskutil eject "$disk" > /dev/null; do
        echo_w "retrying to eject $disk in 5 seconds"
        sleep 5
      done

      echo_i "ejected $disk"
    fi
  done
}

function toolset_create_dmg
{
  local target_dir=$1

  if [ -z "$target_dir" ]; then
    target_dir=$ARTIFACT_DIR
  fi

  if [ "$target_dir" = "$VER_DIR" ]; then
    echo_e "not allowed: target_dir = VER_DIR"
    exit 1
  fi

  # remove files to reduce size
  rm -rf "${BLD_DIR:?}"/*
  find "$SRC_DIR" -mindepth 1 -maxdepth 1 -type d \
    ! -name 'gtk-mac-bundler*' -a \
    ! -name 'jhbuild*' -a \
    ! -name 'png2icns*' \
    -exec rm -rf {} \;

  # create dmg and sha256, print sha256
  cd "$WRK_DIR" || exit 1

  hdiutil create -fs HFS+ -ov -format UDBZ \
    -srcfolder "$VERSION" \
    -volname "$TOOLSET_VOLNAME" \
    "$target_dir/$TOOLSET_VOLNAME".dmg
  shasum -a 256 "$target_dir/$TOOLSET_VOLNAME.dmg" > \
                "$target_dir/$TOOLSET_VOLNAME.dmg".sha256
  cat "$target_dir/$TOOLSET_VOLNAME.dmg.sha256"
}

function toolset_copy
{
  local target=$1

  local toolset_dmg
  toolset_dmg=$TOOLSET_REPO_DIR/$(basename "$TOOLSET_URL")

  if [ -f "$toolset_dmg" ]; then
    echo_i "toolset found: $toolset_dmg"
  else
    # File not present on disk, we need to download.
    echo_i "downloading: $TOOLSET_URL"
    toolset_download
  fi

  toolset_mount "$toolset_dmg"

  echo_i "copying files..."
  rsync -a /Volumes/"$TOOLSET_VOLNAME"/ "$target"

  toolset_unmount /Volumes/"$TOOLSET_VOLNAME"
}