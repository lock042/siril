#!/bin/sh

version="1.4.0-beta3"
siril_dir="siril-$version"

echo 'cloning'
git clone -q -c advice.detachedHead=false https://gitlab.com/free-astro/siril.git --branch "$version" --single-branch $siril_dir
cd $siril_dir
git submodule update --init
rm -rf .git .gitlab subprojects/librtprocess/.git build/flatpak/shared-modules/.git

# add commands reference from the website
echo 'getting command reference'
mkdir -p commands_reference/_images
wget -nv 'https://siril.readthedocs.io/en/stable/Commands.html' -O - | perl -0777 -pe 's/<script.*?script>//gs' | perl -pe 's/<link .*?\/>//gs' > commands_reference/commands_reference.html

grep -o '<img [^>]*src="[^"]*\.svg"[^>]*>' commands_reference/commands_reference.html | \
sed -E 's/<img [^>]*src="([^"]*\.svg)"[^>]*>/\1/g' | \
awk -v baseurl="https://siril.readthedocs.io/en/stable/" '{print baseurl $0}' | \
sort | uniq > svgs.txt

wget -nv -P commands_reference/_images -i svgs.txt

rm svgs.txt

cd ..

tar jcvf $siril_dir.tar.bz2 $siril_dir
sha256sum $siril_dir.tar.bz2
