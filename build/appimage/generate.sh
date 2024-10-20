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

# Python version and download URL
PYTHON_VERSION="3.9.18"
PYTHON_SRC_URL="https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz"
PYTHON_SHORT_VERSION=$(echo $PYTHON_VERSION | cut -d. -f1-2)

# Ensure clean build directory
rm -rf "$APPDIR"
mkdir -p "$APPDIR"

[Previous sections remain the same until Section 6]

# Section 6: AppImage Generation
########################################################################
cd "$BUILDDIR"

# Download linuxdeployqt
wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage

# Verify desktop file exists and is valid
DESKTOP_FILE="$APPDIR/usr/share/applications/org.siril.Siril.desktop"
if [ ! -f "$DESKTOP_FILE" ]; then
    echo "ERROR: Desktop file not found: $DESKTOP_FILE"
    exit 1
fi

# Ensure all required executables exist
SIRIL_BIN="$APPDIR/usr/bin/siril"
PYTHON_BIN="$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}"
PYTHON_WRAPPER="$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}.sh"

for bin in "$SIRIL_BIN" "$PYTHON_BIN" "$PYTHON_WRAPPER"; do
    if [ ! -x "$bin" ]; then
        echo "ERROR: Required executable not found or not executable: $bin"
        exit 1
    fi
done

# Prepare GDK pixbuf loader arguments
linuxdeployqtargs=()
while IFS= read -r -d '' so; do
    if [ -f "$so" ]; then
        linuxdeployqtargs+=("-executable=$(readlink -f "$so")")
    fi
done < <(find "$APPDIR/usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders" -name "*.so" -print0 2>/dev/null || true)

# Add library directories to LD_LIBRARY_PATH
export LD_LIBRARY_PATH="$APPDIR/usr/lib:$APPDIR/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH"

# Create destination directory for AppImage
mkdir -p "$BASE_DIR/dist"

# Generate AppImage with verbose output
./linuxdeployqt-continuous-x86_64.AppImage --appimage-extract-and-run \
    "$DESKTOP_FILE" \
    -verbose=2 \
    -appimage \
    -unsupported-bundle-everything \
    "${linuxdeployqtargs[@]}" \
    -executable="$(readlink -f "$PYTHON_BIN")" \
    -executable="$(readlink -f "$PYTHON_WRAPPER")" \
    -bundle-non-qt-libs \
    -updateinformation="guess" || {
        echo "linuxdeployqt failed with exit code $?"
        # Print contents of AppDir for debugging
        echo "AppDir contents:"
        ls -la "$APPDIR"
        exit 1
    }

# Move AppImage to dist directory
mv Siril*.AppImage* "$BASE_DIR/dist/" || {
    echo "Failed to move AppImage to dist directory"
    exit 1
}

echo "AppImage creation completed successfully"
