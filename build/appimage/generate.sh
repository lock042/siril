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
