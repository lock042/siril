#!/bin/bash

source ~/.profile && jhbuild build icu meta-gtk-osx-freetype meta-gtk-osx-bootstrap meta-gtk-osx-gtk3

if [ $? -ne 0 ]; then
  echo "Installation of pre-built dependencies failed.";
  exit 1;
fi
