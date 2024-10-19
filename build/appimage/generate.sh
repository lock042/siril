#!/bin/bash

# Section 1: Environment Setup
########################################################################
# Set non-interactive mode for apt to prevent prompts during installation
export DEBIAN_FRONTEND=noninteractive

# Python version and download URL
PYTHON_VERSION="3.9.18"  # Specify the exact version
PYTHON_SRC_URL="https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz"

# Section 2: Build Configuration and Installation
########################################################################
BUILDDIR=build/appimage/build
PREFIX=/usr

# Configure build with meson
meson ${BUILDDIR} \
    --prefix=${PREFIX} \
    --buildtype=release \
    -Drelocatable-bundle=yes   # Enable relocatable binaries for AppImage

# Build using ninja with parallel jobs
ninja -C ${BUILDDIR} -j$(nproc)

# Install to temporary AppDir
DESTDIR=$PWD/${BUILDDIR}/appdir ninja -C ${BUILDDIR} -j$(nproc) install
find ${BUILDDIR}/appdir/   # List all installed files

# Change to build directory and setup AppImage structure
cd ${BUILDDIR}
cp ../AppRun appdir/AppRun
chmod +x appdir/AppRun    # Make AppRun executable
cp ./appdir/usr/share/icons/hicolor/scalable/apps/org.siril.Siril.svg ./appdir/org.siril.Siril.svg

# Section 3: Dependencies Bundling
########################################################################
cd appdir/

# Helper function to download and extract debian packages
apt_bundle() {
    apt-get download "$@"              # Download packages
    find *.deb -exec dpkg-deb -x {} . \;   # Extract contents
    find *.deb -delete                # Clean up .deb files
}

# Update package lists
apt update

# Bundle core system libraries
apt_bundle libc6
# Modify library path to prevent loading from system
sed -i -e 's|/usr|/xxx|g' lib/x86_64-linux-gnu/ld-linux-x86-64.so.2

# Bundle GTK dependencies
apt_bundle librsvg2-common libgdk-pixbuf2.0-0 heif-gdk-pixbuf heif-thumbnailer

# Copy and configure GDK pixbuf loaders
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/* usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/
sed -i -e 's|/usr/lib/x86_64-linux-gnu/gdk-pixbuf-.*/.*/loaders/||g' usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache

# Install build dependencies for Python
apt_bundle build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev curl \
    libncursesw5-dev xz-utils tk-dev libxml2-dev \
    libxmlsec1-dev libffi-dev liblzma-dev

# Download and compile Python from source
mkdir -p python_build
cd python_build
wget $PYTHON_SRC_URL
tar xzf Python-${PYTHON_VERSION}.tgz
cd Python-${PYTHON_VERSION}

# Configure Python build with optimization and pip support
./configure --prefix=/usr \
    --enable-optimizations \
    --with-ensurepip=install \
    --enable-shared \
    LDFLAGS="-Wl,-rpath='\$\$ORIGIN/../lib'"

# Compile Python (using all available cores)
make -j$(nproc)

# Install Python to the AppDir
make install DESTDIR=$PWD/../..

cd ../..
rm -rf python_build

# Create Python configuration files
PYTHON_SHORT_VERSION=$(echo $PYTHON_VERSION | cut -d. -f1-2)
cat > usr/lib/python${PYTHON_SHORT_VERSION}/_sysconfigdata.py << EOF
build_time_vars = {
    'PYTHONPATH': '$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}',
    'prefix': '$APPDIR/usr',
    'exec_prefix': '$APPDIR/usr',
    'LIBDIR': '$APPDIR/usr/lib',
}
EOF

# Create a wrapper script to set Python environment variables
cat > usr/bin/python${PYTHON_SHORT_VERSION}.sh << EOF
#!/bin/bash
export PYTHONHOME="\$APPDIR/usr"
export PYTHONPATH="\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}:\$APPDIR/usr/lib/python${PYTHON_SHORT_VERSION}/site-packages"
export LD_LIBRARY_PATH="\$APPDIR/usr/lib:\$APPDIR/usr/lib/x86_64-linux-gnu:\$LD_LIBRARY_PATH"
exec "\$APPDIR/usr/bin/python${PYTHON_SHORT_VERSION}" "\$@"
EOF
chmod +x usr/bin/python${PYTHON_SHORT_VERSION}.sh

# Modify AppRun to set Python environment
sed -i '2i export PYTHONHOME="$APPDIR/usr"' AppRun
sed -i '3i export PYTHONPATH="$APPDIR/usr/lib/python'${PYTHON_SHORT_VERSION}':$APPDIR/usr/lib/python'${PYTHON_SHORT_VERSION}'/site-packages"' AppRun

# Bundle font configuration
mkdir -p etc/fonts/
cp /etc/fonts/fonts.conf etc/fonts/

# Bundle SSL certificates
mkdir -p etc/ssl/
cp -rf /etc/ssl/* etc/ssl/

# Setup and compile GLib schemas
apt_bundle gnome-settings-daemon-common
mkdir -p usr/share/glib-2.0/schemas/
cp /usr/share/glib-2.0/schemas/*.gschema.xml usr/share/glib-2.0/schemas/
if [ -d usr/share/glib-2.0/schemas/ ] ; then
  ( cd usr/share/glib-2.0/schemas/ ; glib-compile-schemas . )
fi
cd -

# Section 4: AppImage Generation
########################################################################
# Download linuxdeployqt tool
wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage

# Prepare arguments for bundling GDK pixbuf loaders
linuxdeployqtargs=()
for so in $(find \
    appdir/usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders \
    -name \*.so); do
    linuxdeployqtargs+=("-executable=$(readlink -f "$so")")
done

# Generate the final AppImage
./linuxdeployqt-continuous-x86_64.AppImage --appimage-extract-and-run appdir/usr/share/applications/org.siril.Siril.desktop \
  -appimage -unsupported-bundle-everything \
  "${linuxdeployqtargs[@]}" \
  -executable=$(readlink -f "appdir/usr/bin/python${PYTHON_SHORT_VERSION}") \
  -executable=$(readlink -f "appdir/usr/bin/python${PYTHON_SHORT_VERSION}.sh")

# Move the generated AppImage to parent directory
mv Siril*.AppImage* ../
