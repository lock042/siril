git submodule update --init

meson _build --buildtype=release --prefix="${INSTALL_PREFIX}"
ninja -C _build install
cp -fr /c/msys64/mingw64 ../${W64_OUT}
