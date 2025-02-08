git submodule update --init

mkdir _build && cd _build
meson setup --buildtype=release --prefix="${INSTALL_PREFIX}" || { echo "Error: failed at meson step"; exit 1; }
ninja install || { echo "Error: failed at ninja step"; exit 1; }
cd ..
