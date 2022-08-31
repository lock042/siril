#!/bin/sh

version=1.0.4
siril_dir="siril-$version"

echo 'cloning'
git clone -q -c advice.detachedHead=false https://gitlab.com/free-astro/siril.git --branch $version --single-branch $siril_dir
cd $siril_dir
git submodule update --init
rm -rf .git .gitlab subprojects/librtprocess/.git build/flatpak/shared-modules/.git

# add commands reference from the website
echo 'getting command reference'
wget -nv 'https://free-astro.org/index.php?title=Siril:Commands&oldid=8045' -O - | perl -0777 -pe 's/<script.*?script>//gs' | perl -pe 's/<link .*?\/>//gs' > commands_reference.html

cd ..

tar jcvf $siril_dir.tar.bz2 $siril_dir
sha256sum $siril_dir.tar.bz2

