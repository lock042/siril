git clone https://gitlab.com/free-astro/siril.git
cd siril
git submodule update --init

meson _build
ninja -C _build install