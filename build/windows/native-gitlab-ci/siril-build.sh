git submodule update --init

cd ./subprojects/yyjson
dos2unix ../../build/yyjson_wchar.patch
patch --binary -p1 < ../../build/yyjson_wchar.patch
cd ../..

mkdir _build && cd _build
meson setup --buildtype=release --prefix="${INSTALL_PREFIX}" || { echo "Error: failed at meson step"; exit 1; }
ninja install || { echo "Error: failed at ninja step"; exit 1; }
cd ..
