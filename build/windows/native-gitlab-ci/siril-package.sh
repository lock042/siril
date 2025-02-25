# Package shared ressources
mkdir -p ${INSTALL_PREFIX}/etc/ssl
mkdir -p ${INSTALL_PREFIX}/python
cp -fr ${W64_OUT}/etc/ssl ${INSTALL_PREFIX}/etc
cp -fr ${W64_OUT}/share/glib-2.0/ ${INSTALL_PREFIX}/share
cp -fr ${W64_OUT}/share/icons/ ${INSTALL_PREFIX}/share
cp -fr ${W64_OUT}/share/locale/ ${INSTALL_PREFIX}/share
cp -fr ${W64_OUT}/share/gtksourceview-4/ ${INSTALL_PREFIX}/share
cp -fr ${PY_OUT}/* ${INSTALL_PREFIX}/python
# Package executable
cp -fr ${W64_OUT}/bin/gdbus.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/glib-compile-schemas.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gdk-pixbuf-query-loaders.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gdk-pixbuf-pixdata.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gdk-pixbuf-thumbnailer.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gspawn-win64-helper.exe ${INSTALL_PREFIX}/bin/
cp -fr ${W64_OUT}/bin/gspawn-win64-helper-console.exe ${INSTALL_PREFIX}/bin/
# Package lib
cp -fr ${W64_OUT}/lib/gdk-pixbuf-2.0 ${INSTALL_PREFIX}/lib/
# Package dll with a Python script using objdump
python3 build/windows/dll_link.py ${INSTALL_PREFIX}/bin/siril.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_siril.txt
python3 build/windows/dll_link.py ${INSTALL_PREFIX}/bin/gdbus.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_gdbus.txt
python3 build/windows/dll_link.py ${W64_OUT}/bin/gdk-pixbuf-query-loaders.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_gdk-pixbuf-query-loaders.txt
python3 build/windows/dll_link.py ${W64_OUT}/bin/gdk-pixbuf-pixdata.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_gdk-pixbuf-pixdata.txt
python3 build/windows/dll_link.py ${W64_OUT}/bin/gdk-pixbuf-thumbnailer.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_gdk-pixbuf-thumbnailer.txt
python3 build/windows/dll_link.py ${W64_OUT}/bin/gspawn-win64-helper.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_gspawn-win64-helper.txt
python3 build/windows/dll_link.py ${W64_OUT}/bin/gspawn-win64-helper-console.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_gspawn-win64-helper-console.txt
python3 build/windows/dll_link.py ${W64_OUT}/bin/glib-compile-schemas.exe ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_glib-compile-schemas.txt
python3 build/windows/dll_link.py ${W64_OUT}/lib/gdk-pixbuf-2.0/2.10.0/loaders/pixbufloader_svg.dll ${W64_OUT}/ ${INSTALL_PREFIX} -l ${CI_PROJECT_DIR}/DLLcache/DLLlist_${SDW_NAME}_pixbufloader_svg.txt