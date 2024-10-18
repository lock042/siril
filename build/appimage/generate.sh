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
# Bundle everything
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

# Add Python and essential packages
apt_bundle python3-minimal python3-pip python3-wheel python3-setuptools
# apt_bundle python3-numpy python3-scipy python3-astropy

# Make absolutely sure it will not load stuff from /lib or /usr
sed -i -e 's|/usr|/xxx|g' lib/x86_64-linux-gnu/ld-linux-x86-64.so.2

# Bundle Gdk pixbuf loaders without which the bundled Gtk does not work;
# this should eventually be done by linuxdeployqt
apt_bundle librsvg2-common libgdk-pixbuf2.0-0 heif-gdk-pixbuf heif-thumbnailer
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/* usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/
sed -i -e 's|/usr/lib/x86_64-linux-gnu/gdk-pixbuf-.*/.*/loaders/||g' usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache

# Bundle Python standard library and dist-packages
mkdir -p usr/lib/python3/
cp -r /usr/lib/python3/* usr/lib/python3/
mkdir -p usr/lib/python3/dist-packages/
cp -r /usr/lib/python3/dist-packages/* usr/lib/python3/dist-packages/

# Fix Python paths in scripts
find usr/bin -type f -exec sed -i -e "s|/usr/bin/python3|/xxx/bin/python3|g" {} \;

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

# Create Python symlinks
ln -sf python3 usr/bin/python
for f in usr/bin/python3.*; do
    if [ -f "$f" ]; then
        ln -sf python3 "$f"
    fi
done

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

# Add Python shared libraries to deployment
for so in $(find \
    appdir/usr/lib/python3* \
    -name \*.so); do
    linuxdeployqtargs+=("-executable=$(readlink -f "$so")")
done

./linuxdeployqt-continuous-x86_64.AppImage --appimage-extract-and-run appdir/usr/share/applications/org.siril.Siril.desktop \
  -appimage -unsupported-bundle-everything \
  "${linuxdeployqtargs[@]}"
mv Siril*.AppImage* ../
