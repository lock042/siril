#!/bin/bash

########################################################################
# Install build-time and run-time dependencies
########################################################################

export DEBIAN_FRONTEND=noninteractive

########################################################################
# Build Siril and install to appdir/
########################################################################

BUILDDIR=build/appimage/build
PREFIX=/usr

# Capture the project root before any cd, so steps below that operate
# inside appdir/ can still reference build/UV_VERSION etc.
PROJECT_ROOT="$PWD"

meson ${BUILDDIR} \
    --prefix=${PREFIX} \
    --buildtype=release \
    -Drelocatable-bundle=yes


ninja -C ${BUILDDIR} -j$(nproc)
DESTDIR=$PWD/${BUILDDIR}/appdir ninja -C ${BUILDDIR} -j$(nproc) install; find ${BUILDDIR}/appdir/
cd ${BUILDDIR}
cp ../AppRun appdir/AppRun ; chmod +x appdir/AppRun
cp ./appdir/usr/share/icons/hicolor/scalable/apps/org.siril.Siril.svg ./appdir/org.siril.Siril.svg

cd appdir/

########################################################################
# Bundle everyhting
# to allow the AppImage to run on older systems as well
########################################################################

apt_bundle() {
    apt-get download "$@"
    find *.deb -exec dpkg-deb -x {} . \;
    find *.deb -delete
}

# Bundle all of glibc; this should eventually be done by linuxdeployqt
apt update
apt_bundle libc6

# Make absolutely sure it will not load stuff from /lib or /usr
sed -i -e 's|/usr|/xxx|g' lib/x86_64-linux-gnu/ld-linux-x86-64.so.2

# Bundle Gdk pixbuf loaders without which the bundled Gtk does not work;
# this should eventually be done by linuxdeployqt
apt_bundle librsvg2-common libgdk-pixbuf2.0-0 heif-gdk-pixbuf heif-thumbnailer
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/* usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/
sed -i -e 's|/usr/lib/x86_64-linux-gnu/gdk-pixbuf-.*/.*/loaders/||g' usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache

# Bundle fontconfig settings
mkdir -p etc/fonts/
cp /etc/fonts/fonts.conf etc/fonts/

# Bundle ssl certificates
mkdir -p etc/ssl/
cp -rf /etc/ssl/* etc/ssl/

# Compile GLib schemas if the subdirectory is present in the AppImage
# AppRun has to export GSETTINGS_SCHEMA_DIR for this to work
apt_bundle gnome-settings-daemon-common
mkdir -p usr/share/glib-2.0/schemas/
cp /usr/share/glib-2.0/schemas/*.gschema.xml usr/share/glib-2.0/schemas/
if [ -d usr/share/glib-2.0/schemas/ ] ; then
  ( cd usr/share/glib-2.0/schemas/ ; glib-compile-schemas . )
fi

########################################################################
# Bundle uv (Astral's standalone static binary; no shared-library deps).
#
# Required so the per-script venv path (src/io/siril_pythonmodule.c)
# can find uv at runtime. AppRun:28 already prepends ${HERE}/usr/bin to
# PATH, so g_find_program_in_path("uv") will pick it up automatically.
#
# Version pinned in build/UV_VERSION (single source of truth). Bump
# both the version and the sha256 below in lockstep.
########################################################################
UV_VERSION=$(cat "${PROJECT_ROOT}/build/UV_VERSION")
UV_SHA256_x86_64="74947fe2c03315cf07e82ab3acc703eddef01aba4d5232a98e4c6825ec116131"
UV_URL="https://github.com/astral-sh/uv/releases/download/${UV_VERSION}/uv-x86_64-unknown-linux-gnu.tar.gz"
wget -c -nv "${UV_URL}" -O uv.tar.gz
echo "${UV_SHA256_x86_64}  uv.tar.gz" | sha256sum -c -
mkdir -p tmp-uv
tar -xzf uv.tar.gz -C tmp-uv --strip-components=1
install -m 0755 tmp-uv/uv usr/bin/uv
install -m 0755 tmp-uv/uvx usr/bin/uvx
rm -rf tmp-uv uv.tar.gz

cd -

########################################################################
# Generate AppImage
########################################################################

wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage

linuxdeployqtargs=()

for so in $(find \
    appdir/usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders \
    -name \*.so); do
    linuxdeployqtargs+=("-executable=$(readlink -f "$so")")
done

./linuxdeployqt-continuous-x86_64.AppImage --appimage-extract-and-run appdir/usr/share/applications/org.siril.Siril.desktop \
  -appimage -unsupported-bundle-everything \
  "${linuxdeployqtargs[@]}"

mv Siril*.AppImage* ../

