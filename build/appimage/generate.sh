#!/bin/bash

# Section 1: Environment Setup
########################################################################
set -e  # Exit on error
export DEBIAN_FRONTEND=noninteractive

# Python version and download URL
PYTHON_VERSION="3.9.18"
PYTHON_SRC_URL="https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz"
PYTHON_SHORT_VERSION=$(echo $PYTHON_VERSION | cut -d. -f1-2)

# Section 2: Build Configuration
########################################################################
BUILDDIR="build/appimage/build"
APPDIR="$BUILDDIR/appdir"
PREFIX="/usr"

# Ensure clean build directory
rm -rf "$APPDIR"
mkdir -p "$APPDIR"

# Configure build with meson
meson ${BUILDDIR} \
    --prefix=${PREFIX} \
    --buildtype=release \
    -Drelocatable-bundle=yes

# Build using ninja
ninja -C ${BUILDDIR} -j$(nproc)

# Install to AppDir
DESTDIR="$PWD/$APPDIR" ninja -C ${BUILDDIR} -j$(nproc) install

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
PYTHON_BUILD_DIR="$APPDIR/python_build"
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

# Clean up Python build
cd "$APPDIR"
rm -rf python_build

# Create Python configuration
mkdir -p "usr/lib/python${PYTHON_SHORT_VERSION}"
cat > "usr/lib/python${PYTHON_SHORT_VERSION}/_sysconfigdata.py" << EOF
build_time_vars = {
    'PYTHONPATH': '\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}',
    'prefix': '\$APPDIR/usr',
    'exec_prefix': '\$APPDIR/usr',
    'LIBDIR': '\$APPDIR/usr/lib',
}
EOF

# Create Python wrapper script
mkdir -p usr/bin
cat > "usr/bin/python${PYTHON_SHORT_VERSION}.sh" << EOF
#!/bin/bash
export PYTHONHOME="\$APPDIR/usr"
export PYTHONPATH="\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}:\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}/site-packages"
export LD_LIBRARY_PATH="\$APPDIR/usr/lib:\$APPDIR/usr/lib/x86_64-linux-gnu:\$LD_LIBRARY_PATH"
exec "\$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}" "\$@"
EOF
chmod +x "usr/bin/python${PYTHON_SHORT_VERSION}.sh"

# Configure AppRun with Python environment
sed -i "2i export PYTHONHOME=\"\$APPDIR/usr\"" AppRun
sed -i "3i export PYTHONPATH=\"\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}:\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}/site-packages\"" AppRun

# Section 4: System Configuration Files
########################################################################
# Bundle font configuration
mkdir -p etc/fonts/
cp /etc/fonts/fonts.conf etc/fonts/

# Bundle SSL certificates
mkdir -p etc/ssl/
cp -rf /etc/ssl/certs etc/ssl/
cp -rf /etc/ssl/openssl.cnf etc/ssl/
mkdir -p etc/ssl/private

# Section 5: GLib Schema Configuration
########################################################################
# Install schema files
mkdir -p usr/share/glib-2.0/schemas/
cp -f /usr/share/glib-2.0/schemas/*.xml usr/share/glib-2.0/schemas/ || true

# Compile schemas (with error handling)
if [ -d usr/share/glib-2.0/schemas/ ]; then
    # First, ensure all required enum files are present
    for schema in usr/share/glib-2.0/schemas/*.xml; do
        if [ -f "$schema" ]; then
            # Check if schema references missing enums and remove if so
            if grep -q "enum id='org.gnome.desktop.GDesktop" "$schema"; then
                rm "$schema"
            fi
        fi
    done

    # Now compile the remaining schemas
    (cd usr/share/glib-2.0/schemas/ && glib-compile-schemas .)
fi

# Section 6: AppImage Generation
########################################################################
cd "$BUILDDIR"

# Download linuxdeployqt
wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage

# Prepare GDK pixbuf loader arguments
linuxdeployqtargs=()
while IFS= read -r -d '' so; do
    linuxdeployqtargs+=("-executable=$(readlink -f "$so")")
done < <(find appdir/usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders -name "*.so" -print0)

# Generate AppImage
./linuxdeployqt-continuous-x86_64.AppImage --appimage-extract-and-run appdir/usr/share/applications/org.siril.Siril.desktop \
    -appimage -unsupported-bundle-everything \
    "${linuxdeployqtargs[@]}" \
    -executable="$(readlink -f "appdir/usr/bin/python${PYTHON_SHORT_VERSION}")" \
    -executable="$(readlink -f "appdir/usr/bin/python${PYTHON_SHORT_VERSION}.sh")"

# Move AppImage to parent directory
mv Siril*.AppImage* ../
