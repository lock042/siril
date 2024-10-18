#!/bin/bash
set -e

# Configuration
APP="siril"
ARCH="x86_64"
PYTHON_VERSION="3.10"

# Create application root directory
APP_DIR="${APP}.AppDir"
mkdir -p "$APP_DIR"
cd "$APP_DIR"

# Create directory structure
mkdir -p usr/bin
mkdir -p usr/lib
mkdir -p usr/lib/python${PYTHON_VERSION}
mkdir -p usr/lib/python${PYTHON_VERSION}/site-packages
mkdir -p usr/lib/${ARCH}-linux-gnu
mkdir -p usr/share
mkdir -p etc/fonts

# Copy main application binary and its dependencies
cp /usr/bin/${APP} usr/bin/
copy_dependencies() {
    local binary="$1"
    local dest_dir="$2"
    ldd "$binary" | grep "=> /" | awk '{print $3}' | while read -r lib; do
        if [ ! -f "$dest_dir/$(basename "$lib")" ]; then
            cp -L "$lib" "$dest_dir/"
        fi
    done
}

# Copy main application dependencies
copy_dependencies "/usr/bin/${APP}" "usr/lib"

# Copy Python and its dependencies
cp /usr/bin/python${PYTHON_VERSION} usr/bin/
ln -sf python${PYTHON_VERSION} usr/bin/python3
ln -sf python3 usr/bin/python

# Copy Python standard library
cp -r /usr/lib/python${PYTHON_VERSION}/* usr/lib/python${PYTHON_VERSION}/

# Copy Python dependencies
copy_dependencies "usr/bin/python${PYTHON_VERSION}" "usr/lib"

# Copy Python site-packages
cp -r /usr/lib/python3/dist-packages usr/lib/python${PYTHON_VERSION}/ || true
cp -r /usr/lib/python${PYTHON_VERSION}/site-packages/* usr/lib/python${PYTHON_VERSION}/site-packages/ || true

# Copy GTK and GLib dependencies
cp -r /usr/lib/${ARCH}-linux-gnu/gdk-pixbuf-* usr/lib/${ARCH}-linux-gnu/
cp -r /usr/lib/${ARCH}-linux-gnu/girepository-1.0 usr/lib/${ARCH}-linux-gnu/
cp -r /usr/lib/${ARCH}-linux-gnu/gconv usr/lib/${ARCH}-linux-gnu/

# Copy necessary shared data
cp -r /usr/share/glib-2.0 usr/share/
cp -r /usr/share/icons usr/share/
cp -r /usr/share/${APP} usr/share/
cp -r /etc/fonts/fonts.conf etc/fonts/

# Generate GDK pixbuf loaders cache
GDK_PIXBUF_MODULEDIR=$(ls -d usr/lib/${ARCH}-linux-gnu/gdk-pixbuf-*/*/loaders/)
GDK_PIXBUF_MODULE_FILE="${GDK_PIXBUF_MODULEDIR}/loaders.cache"
gdk-pixbuf-query-loaders > "$GDK_PIXBUF_MODULE_FILE"
sed -i "s|/usr/lib/${ARCH}-linux-gnu/|./usr/lib/${ARCH}-linux-gnu/|g" "$GDK_PIXBUF_MODULE_FILE"

# Create AppRun script
cat > AppRun << 'EOL'
#!/bin/bash
HERE="$(dirname "$(readlink -f "${0}")")"

# Allow AppRun or the AppImage to be symlinked to
if [ ! -z $APPIMAGE ] ; then
    BINARY_NAME=$(basename "$ARGV0")
else
    BINARY_NAME=$(basename "$0")
fi

if [ ! -z "$1" ] && [ -e "$HERE/usr/bin/$1" ] ; then
    MAIN="$HERE/usr/bin/$1"
    shift
elif [ -e "$HERE/usr/bin/$BINARY_NAME" ] ; then
    MAIN="$HERE/usr/bin/$BINARY_NAME"
else
    MAIN="$HERE/usr/bin/siril"
fi

# Python environment settings
export PYTHONHOME="${HERE}/usr"
export PYTHONPATH="${HERE}/usr/lib/python3.10:${HERE}/usr/lib/python3.10/site-packages:${PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
export LD_LIBRARY_PATH="${HERE}/usr/lib:${LD_LIBRARY_PATH}"

# Platform settings
PLATFORM=x86_64-linux-gnu

# Glib/Gtk environment
export GCONV_PATH="$HERE/usr/lib/${PLATFORM}/gconv"
export FONTCONFIG_FILE="$HERE/etc/fonts/fonts.conf"
export GTK_EXE_PREFIX="$HERE/usr"
export GDK_PIXBUF_MODULEDIR=$(readlink -f "$HERE"/usr/lib/${PLATFORM}/gdk-pixbuf-*/*/loaders/)
export GDK_PIXBUF_MODULE_FILE=$(readlink -f "$HERE"/usr/lib/${PLATFORM}/gdk-pixbuf-*/*/loaders.cache)
export GI_TYPELIB_PATH="${HERE}/usr/lib/${PLATFORM}/girepository-1.0"
export XDG_DATA_DIRS="$HERE/usr/share:${XDG_DATA_DIRS:-/usr/local/share:/usr/share}"

# Library paths
LIBRARY_PATH="$HERE/lib/${PLATFORM}"
LIBRARY_PATH+=":$HERE/usr/lib/${PLATFORM}"
LIBRARY_PATH+=":$HERE/usr/lib"
LIBRARY_PATH+=":$GDK_PIXBUF_MODULEDIR"

export LD_LIBRARY_PATH="$GDK_PIXBUF_MODULEDIR:$LD_LIBRARY_PATH"
export PATH="${HERE}/usr/bin:${PATH}"

exec "$HERE/lib/${PLATFORM}/ld-linux-x86-64.so.2" --inhibit-cache --library-path "$LIBRARY_PATH" "$MAIN" "$@"
EOL

chmod +x AppRun

# Copy desktop file and icon
cp /usr/share/applications/${APP}.desktop .
cp /usr/share/icons/hicolor/256x256/apps/${APP}.png .

# Copy ld-linux
mkdir -p lib/${ARCH}-linux-gnu/
cp /lib/${ARCH}-linux-gnu/ld-linux-x86-64.so.2 lib/${ARCH}-linux-gnu/

# Fix permissions
chmod +x usr/bin/*

echo "AppDir structure created. Creating AppImage..."
cd ..

# Download appimagetool if not present
if [ ! -f appimagetool-x86_64.AppImage ]; then
    wget "https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage"
    chmod +x appimagetool-x86_64.AppImage
fi

# Create AppImage
VERSION=$(date +%Y%m%d)
./appimagetool-x86_64.AppImage --no-appstream ${APP}.AppDir ${APP}-${VERSION}-${ARCH}.AppImage

echo "AppImage created: ${APP}-${VERSION}-${ARCH}.AppImage"
