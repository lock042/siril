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
    -Drelocatable-bundle=yes \
    -Denable-libcurl=yes


ninja -C ${BUILDDIR} -j$(nproc)
DESTDIR=$PWD/${BUILDDIR}/appdir ninja -C ${BUILDDIR} -j$(nproc) install; find ${BUILDDIR}/appdir/
cd ${BUILDDIR}
cp ../AppRun appdir/AppRun ; chmod +x appdir/AppRun
cp ./appdir/usr/share/icons/hicolor/scalable/apps/org.free_astro.siril.svg ./appdir/org.free_astro.siril.svg

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
apt_bundle librsvg2-common libgdk-pixbuf2.0-0
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/* usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders/
cp /usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/
sed -i -e 's|/usr/lib/x86_64-linux-gnu/gdk-pixbuf-.*/.*/loaders/||g' usr/lib/x86_64-linux-gnu/gdk-pixbuf-*/*/loaders.cache

# Bundle fontconfig settings
mkdir -p etc/fonts/
cp /etc/fonts/fonts.conf etc/fonts/

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

./linuxdeployqt-continuous-x86_64.AppImage --appimage-extract-and-run appdir/usr/share/applications/org.free_astro.siril.desktop \
  -appimage -unsupported-bundle-everything \
  "${linuxdeployqtargs[@]}"

mv Siril*.AppImage* ../

