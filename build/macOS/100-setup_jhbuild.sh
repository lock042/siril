#!/bin/bash

cd $HOME
mkdir -p ~/project && mkdir -p ~/.config
curl -L 'https://gitlab.com/free-astro/siril-macos-build/-/raw/mainline/jhbuildrc-gtk-osx-gimp-2.99' > ~/.config/jhbuildrc-custom
curl 'https://gitlab.com/free-astro/siril-macos-build/-/raw/mainline/gtk-osx-setup.sh' > ~/project/gtk-osx-setup.sh
chmod +x gtk-osx-setup.sh
echo 'export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH:$HOME/.new_local/bin"' >> ~/.profile
echo 'export ARCHFLAGS="-arch x86_64"' >> ~/.profile
# PYTHON variable seems to be incorrectly set
echo 'export PYTHON=/Library/Frameworks/Python.framework/Versions/3.9/bin/python3' >> ~/.profile
# Unclear why this is needed, but if missing, jhbuild fails in some circumstances
echo 'export PYENV_VERSION="3.9.7"' >> ~/.profile
source ~/.profile
PIPENV_YES=1 ~/project/gtk-osx-setup.sh
$HOME/.new_local/bin/jhbuild bootstrap-gtk-osx-gimp
cat ~/.profile

if [ $? -ne 0 ]; then
  echo "Installation of jhbuild failed.";
  exit 1;
fi
