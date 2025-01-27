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
    --buildtype=debug \
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

##############################################################################

# Bundle GTK and core dependencies
apt_bundle \
    libgtk-3-0 \
    libglib2.0-0 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libcairo2 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libepoxy0 \
    libxkbcommon0 \
    libwayland-client0 \
    libwayland-cursor0 \
    libwayland-egl1 \
    libx11-6 \
    libxcomposite1 \
    libxdamage1 \
    libxext6 \
    libxfixes3 \
    libxinerama1 \
    libxrandr2 \
    libgtk-3-common \
    gvfs-common \
    gvfs-libs \
    gvfs-daemons \
    gvfs

# Copy gvfs libs
mkdir -p appdir/usr/lib/x86_64-linux-gnu/
cp /usr/lib/x86_64-linux-gnu/libgvfscommon.so* appdir/usr/lib/x86_64-linux-gnu/

# Copy GTK modules
#mkdir -p usr/lib/x86_64-linux-gnu/gtk-3.0/
#cp -r /usr/lib/x86_64-linux-gnu/gtk-3.0/modules usr/lib/x86_64-linux-gnu/gtk-3.0/
#cp -r /usr/lib/x86_64-linux-gnu/gtk-3.0/immodules usr/lib/x86_64-linux-gnu/gtk-3.0/

# Copy GTK settings schemas
mkdir -p usr/share/glib-2.0/schemas
cp /usr/share/glib-2.0/schemas/org.gtk.Settings.* usr/share/glib-2.0/schemas/

# Compile GTK immodules cache
#gtk-query-immodules-3.0 > usr/lib/x86_64-linux-gnu/gtk-3.0/3.0.0/immodules.cache
#sed -i -e 's|/usr/lib/x86_64-linux-gnu/||g' usr/lib/x86_64-linux-gnu/gtk-3.0/3.0.0/immodules.cache

###############################################################################

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

