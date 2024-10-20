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

# Function to verify AppImage tool
verify_appimage_tool() {
    local tool="$1"
    echo "Verifying AppImage tool: $tool"
    if ! "$tool" --version; then
        echo "Failed to run AppImage tool"
        # Try extracting if direct execution fails
        "$tool" --appimage-extract
        LINUXDEPLOY_DIR="squashfs-root"
        if [ -d "$LINUXDEPLOY_DIR" ]; then
            echo "Using extracted AppImage tool"
            LINUXDEPLOY="$PWD/$LINUXDEPLOY_DIR/AppRun"
            return 0
        fi
        return 1
    fi
}

# Download linuxdeployqt
wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage

# Verify and potentially extract linuxdeployqt
LINUXDEPLOY="./linuxdeployqt-continuous-x86_64.AppImage"
verify_appimage_tool "$LINUXDEPLOY"

# Verify desktop file exists and is valid
DESKTOP_FILE="$APPDIR/usr/share/applications/org.siril.Siril.desktop"
if [ ! -f "$DESKTOP_FILE" ]; then
    echo "ERROR: Desktop file not found: $DESKTOP_FILE"
    exit 1
fi

# Ensure all executables exist and are actually executable
for bin in "$APPDIR/usr/bin/siril" "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}" "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}.sh"; do
    if [ ! -x "$bin" ]; then
        echo "ERROR: $bin is not executable"
        ls -l "$bin"
        chmod +x "$bin"
    fi
done

# Verify AppDir structure
echo "AppDir structure before linuxdeployqt:"
find "$APPDIR" -type f -exec file {} \;

# Create log directory
mkdir -p "$BUILDDIR/logs"

# Try AppImage creation with different options
(
    # First attempt: Basic deployment
    "$LINUXDEPLOY" \
        "$DESKTOP_FILE" \
        -verbose=3 \
        -appimage \
        -unsupported-bundle-everything \
        -executable="$(readlink -f "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}")" \
        -executable="$(readlink -f "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}.sh")" \
        -bundle-non-qt-libs \
        2>"$BUILDDIR/logs/linuxdeployqt.log" || {

        # Second attempt: With additional options
        echo "First attempt failed, trying with additional options..."
        "$LINUXDEPLOY" \
            "$DESKTOP_FILE" \
            -verbose=3 \
            -appimage \
            -unsupported-bundle-everything \
            -executable="$(readlink -f "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}")" \
            -executable="$(readlink -f "$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}.sh")" \
            -bundle-non-qt-libs \
            -no-translations \
            -no-plugins \
            2>"$BUILDDIR/logs/linuxdeployqt_attempt2.log" || {

            # Third attempt: Minimal deployment
            echo "Second attempt failed, trying minimal deployment..."
            "$LINUXDEPLOY" \
                "$DESKTOP_FILE" \
                -verbose=3 \
                -appimage \
                -bundle-non-qt-libs \
                2>"$BUILDDIR/logs/linuxdeployqt_attempt3.log" || {
                    echo "All attempts failed. Collecting diagnostic information..."
                    echo "=== Last 50 lines of first attempt log ==="
                    tail -n 50 "$BUILDDIR/logs/linuxdeployqt.log"
                    echo "=== Last 50 lines of second attempt log ==="
                    tail -n 50 "$BUILDDIR/logs/linuxdeployqt_attempt2.log"
                    echo "=== Last 50 lines of third attempt log ==="
                    tail -n 50 "$BUILDDIR/logs/linuxdeployqt_attempt3.log"
                    echo "=== AppDir final structure ==="
                    find "$APPDIR" -type f -exec file {} \;
                    exit 1
                }
            }
        }
)

# If we get here, one of the attempts succeeded
echo "AppImage creation completed successfully"
mv Siril*.AppImage* "$BASE_DIR/dist/" || {
    echo "Failed to move AppImage to dist directory"
    exit 1
}
