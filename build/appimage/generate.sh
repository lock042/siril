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

# Section 2: Build Configuration
########################################################################
# Configure build with meson
meson ${BUILDDIR} \
    --prefix=${PREFIX} \
    --buildtype=release \
    -Drelocatable-bundle=yes

# Build using ninja
ninja -C ${BUILDDIR} -j$(nproc)

# Install to AppDir
DESTDIR="$APPDIR" ninja -C ${BUILDDIR} -j$(nproc) install

# Create and setup AppRun script
cat > "$APPDIR/AppRun" << 'EOF'
#!/bin/bash

# Find the AppDir
HERE="$(dirname "$(readlink -f "${0}")")"
APPDIR="${HERE}"

# Export path
export PATH="${APPDIR}/usr/bin:${PATH}"
export LD_LIBRARY_PATH="${APPDIR}/usr/lib:${APPDIR}/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}"

# Add any additional environment variables here
export XDG_DATA_DIRS="${APPDIR}/usr/share:${XDG_DATA_DIRS}"
export GDK_PIXBUF_MODULE_FILE="${APPDIR}/usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0/2.10.0/loaders.cache"
export GDK_PIXBUF_MODULEDIR="${APPDIR}/usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0/2.10.0/loaders"

# Python environment will be added here by sed commands later
# PYTHON ENV MARKER - DO NOT REMOVE THIS LINE

# Execute the application
exec "${APPDIR}/usr/bin/siril" "$@"
EOF

chmod +x "$APPDIR/AppRun"

# Copy icon
cp "$APPDIR/usr/share/icons/hicolor/scalable/apps/org.siril.Siril.svg" "$APPDIR/org.siril.Siril.svg"

# Section 3: Python Installation
########################################################################
# Create Python build directory
PYTHON_BUILD_DIR="$BUILDDIR/python_build"
mkdir -p "$PYTHON_BUILD_DIR"
cd "$PYTHON_BUILD_DIR"

# Download and extract Python
wget "$PYTHON_SRC_URL"
tar xzf "Python-${PYTHON_VERSION}.tgz"
cd "Python-${PYTHON_VERSION}"

# Configure Python
./configure \
    --prefix=/usr \
    --enable-optimizations \
    --with-ensurepip=install \
    --enable-shared \
    LDFLAGS="-Wl,-rpath='\$\$ORIGIN/../lib'"

# Build Python
make -j$(nproc)

# Install Python directly to AppDir
make install DESTDIR="$APPDIR"

# Return to base directory and clean up Python build
cd "$BASE_DIR"
rm -rf "$PYTHON_BUILD_DIR"

# Create Python configuration
mkdir -p "$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}"
cat > "$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}/_sysconfigdata.py" << EOF
build_time_vars = {
    'PYTHONPATH': '\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}',
    'prefix': '\$APPDIR/usr',
    'exec_prefix': '\$APPDIR/usr',
    'LIBDIR': '\$APPDIR/usr/lib',
}
EOF

# Create Python wrapper script
mkdir -p "$APPDIR/usr/bin"
cat > "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}.sh" << EOF
#!/bin/bash
export PYTHONHOME="\$APPDIR/usr"
export PYTHONPATH="\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}:\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}/site-packages"
export LD_LIBRARY_PATH="\$APPDIR/usr/lib:\$APPDIR/usr/lib/x86_64-linux-gnu:\$LD_LIBRARY_PATH"
exec "\$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}" "\$@"
EOF
chmod +x "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}.sh"

# Add Python environment to AppRun
sed -i "/# PYTHON ENV MARKER/a export PYTHONPATH=\"\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}:\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}/site-packages\"" "$APPDIR/AppRun"
sed -i "/# PYTHON ENV MARKER/a export PYTHONHOME=\"\$APPDIR/usr\"" "$APPDIR/AppRun"

# Section 4: System Configuration Files
########################################################################
# Bundle font configuration
mkdir -p "$APPDIR/etc/fonts/"
cp /etc/fonts/fonts.conf "$APPDIR/etc/fonts/"

# Bundle SSL certificates
mkdir -p "$APPDIR/etc/ssl/"
cp -rf /etc/ssl/certs "$APPDIR/etc/ssl/"
cp -rf /etc/ssl/openssl.cnf "$APPDIR/etc/ssl/"
mkdir -p "$APPDIR/etc/ssl/private"

# Section 5: GLib Schema Configuration
########################################################################
# Install schema files
mkdir -p "$APPDIR/usr/share/glib-2.0/schemas/"
cp -f /usr/share/glib-2.0/schemas/*.xml "$APPDIR/usr/share/glib-2.0/schemas/" || true

# Compile schemas (with error handling)
if [ -d "$APPDIR/usr/share/glib-2.0/schemas/" ]; then
    # First, ensure all required enum files are present
    for schema in "$APPDIR/usr/share/glib-2.0/schemas/"*.xml; do
        if [ -f "$schema" ]; then
            # Check if schema references missing enums and remove if so
            if grep -q "enum id='org.gnome.desktop.GDesktop" "$schema"; then
                rm "$schema"
            fi
        fi
    done

    # Now compile the remaining schemas
    (cd "$APPDIR/usr/share/glib-2.0/schemas/" && glib-compile-schemas .)
fi


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
        fi
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
