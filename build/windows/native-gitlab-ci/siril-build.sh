git submodule update --init

mkdir _build && cd _build
meson setup --buildtype=release --prefix="${INSTALL_PREFIX}"
ninja install
cd ..
