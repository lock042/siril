# Package shared ressources
mkdir -p ${INSTALL_PREFIX}/etc/ssl
cp -fr ${W64_OUT}/etc/ssl ${INSTALL_PREFIX}/etc
cp -fr ${W64_OUT}/share/glib-2.0/ ${INSTALL_PREFIX}/share
cp -fr ${W64_OUT}/share/icons/ ${INSTALL_PREFIX}/share
cp -fr ${W64_OUT}/share/locale/ ${INSTALL_PREFIX}/share
# Package executable
cp -fr ${W64_OUT}/bin/gdbus.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/glib-compile-schemas.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gdk-pixbuf-query-loaders.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gdk-pixbuf-pixdata.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gdk-pixbuf-thumbnailer.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gspawn-win64-helper.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gspawn-win64-helper-console.exe ${INSTALL_PREFIX}/bin/
cp build/windows/crossbuild-gitlab-ci/siril-wrapper.cmd ${INSTALL_PREFIX}/bin/
# Package lib
cp -fr ${W64_OUT}/lib/gdk-pixbuf-2.0 ${INSTALL_PREFIX}/lib/
# Package dll with a Python script using objdump
python3 build/windows/dll_link.py ${INSTALL_PREFIX}/bin/siril.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${INSTALL_PREFIX}/bin/gdbus.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${W64_OUT}/bin/gdk-pixbuf-query-loaders.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${W64_OUT}/bin/gdk-pixbuf-pixdata.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${W64_OUT}/bin/gdk-pixbuf-thumbnailer.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${W64_OUT}/bin/gspawn-win64-helper.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${W64_OUT}/bin/gspawn-win64-helper-console.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${W64_OUT}/bin/glib-compile-schemas.exe ${W64_OUT}/ ${INSTALL_PREFIX}
python3 build/windows/dll_link.py ${W64_OUT}/lib/gdk-pixbuf-2.0/2.10.0/loaders/pixbufloader_svg.dll ${W64_OUT}/ ${INSTALL_PREFIX}
cd "${INSTALL_PREFIX}/bin"
./gdk-pixbuf-query-loaders.exe --update-cache
./glib-compile-schemas.exe --targetdir="../share/glib-2.0/schemas" "../share/glib-2.0/schemas"