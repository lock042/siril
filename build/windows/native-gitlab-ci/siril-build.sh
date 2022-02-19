git submodule update --init

mkdir _build && cd _build
meson --buildtype=release --prefix="${INSTALL_PREFIX}"
ninja install
cp -fr /c/msys64/mingw64 ../${W64_OUT}
