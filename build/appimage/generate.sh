#!/bin/bash

# Section 1: Environment Setup
########################################################################
set -e  # Exit on error
set -x  # Print commands for debugging
export DEBIAN_FRONTEND=noninteractive

# Store base directory
BASE_DIR="$PWD"
BUILDDIR="$BASE_DIR/build/appimage/build"
APPDIR="$BUILDDIR/appdir"
PREFIX="/usr"

[Previous sections remain the same until the end of AppImage creation]

# Section 6: AppImage Generation
########################################################################
cd "$BUILDDIR"

# Download linuxdeployqt
wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage

# Extract AppImage tool
./linuxdeployqt-continuous-x86_64.AppImage --appimage-extract
LINUXDEPLOY="$PWD/squashfs-root/AppRun"

# Verify desktop file exists and is valid
DESKTOP_FILE="$APPDIR/usr/share/applications/org.siril.Siril.desktop"
if [ ! -f "$DESKTOP_FILE" ]; then
    echo "ERROR: Desktop file not found: $DESKTOP_FILE"
    exit 1
fi

# Create AppImage with minimal deployment options (since we know this works)
"$LINUXDEPLOY" \
    "$DESKTOP_FILE" \
    -verbose=3 \
    -appimage \
    -bundle-non-qt-libs || {
        echo "AppImage creation failed"
        exit 1
    }

# Section 7: AppImage Finalization
########################################################################
# Create all necessary parent directories
mkdir -p "$BASE_DIR/dist"

# Find the generated AppImage files
APPIMAGE_FILES=(Siril*.AppImage*)
if [ ${#APPIMAGE_FILES[@]} -eq 0 ]; then
    echo "ERROR: No AppImage files found"
    ls -la
    exit 1
fi

# Move each AppImage file individually
for file in "${APPIMAGE_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "Moving $file to $BASE_DIR/dist/"
        cp "$file" "$BASE_DIR/dist/" || {
            echo "Failed to copy $file to dist directory"
            exit 1
        }
        # Verify the copy succeeded before removing the original
        if [ -f "$BASE_DIR/dist/$(basename "$file")" ]; then
            rm "$file"
        else
            echo "ERROR: Failed to verify copied file: $BASE_DIR/dist/$(basename "$file")"
            exit 1
        }
    else
        echo "WARNING: AppImage file not found: $file"
    fi
done

# Verify the files were moved successfully
if [ "$(ls -A "$BASE_DIR/dist/")" ]; then
    echo "AppImage files successfully moved to dist directory:"
    ls -la "$BASE_DIR/dist/"
    echo "AppImage creation and deployment completed successfully"
else
    echo "ERROR: No files found in dist directory after move operation"
    exit 1
fi
