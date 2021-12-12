# Build SiriL/macOS inside CircleCI

This repository contains files related to SiriL/macOS build using CircleCI.

## Build process description

To build SiriL/macOS we are using [fork](https://gitlab.gnome.org/samm-git/gtk-osx/tree/gimp)
of the [gtk-osx](https://gitlab.gnome.org/GNOME/gtk-osx) project (`gimp` branch).
Fork adds modules related to SiriL and some specific patches to GTK.
Currently build is done using CircleCI.

## Before you starting

I found that GTK build process on macOS is very fragile. If you have any other build system (brew, MacPorts) installed - try to remove it first or at least isolate from the JHBuild environment as much as you can.

I was able to get working builds in the VirtualBox VM, it works stable enough for me.

## Steps in the CircleCI [config.yml](https://gitlab.gnome.org/Infrastructure/gimp-macos-build/blob/master/.circleci/config.yml) are:

- Setting up macOS 10.9 SDK. This is needed to ensure that SiriL can run on macOS 10.9+. See [this article](https://smallhacks.wordpress.com/2018/11/11/how-to-support-old-osx-version-with-a-recent-xcode/) for the details.
- Setting up JHBuild with a custom `~/.config/jhbuildrc-custom` file (see https://github.com/GNOME/gimp-macos-build/blob/master/jhbuildrc-gtk-osx-gimp). As part of the setup, it is running `bootstrap-gtk-osx-gimp` JHBuild command to compile required modules to run JHBuild. JHBuild is using Python3 venv to run.
- Installs [fork of the gtk-mac-bundler](https://github.com/samm-git/gtk-mac-bundler/tree/fix-otool) - the tool which helps to create macOS application bundles for the GTK apps. The only difference with official one is [this PR](https://github.com/jralls/gtk-mac-bundler/pull/10)
- Installing all gtk-osx and SiriL dependencies using JHBuild
- Building SiriL (from the git).
- Importing signing certificate/key from the environment variables
- Launching `build.sh` which:
  - Building package using `gtk-mac-bundler`
  - Signing all binaries *(not yet enabled)*
  - Creating a DMG package using [create-dmg](https://github.com/andreyvit/create-dmg) tool and signing it
- Notarizing package using Apple `altool` utility *(not yet enabled)*
- Uploading a DMG to the CircleCI build artifacts

## Other related links

 - [Gtk-OSX](https://gitlab.gnome.org/GNOME/gtk-osx/) project to simplify building MacOS application bundles for Gtk+-based applications
 - [gimp-plugins-collection](https://github.com/aferrero2707/gimp-plugins-collection)
 - CircleCI [siril-macos-build project](https://circleci.com/gh/samm-git/siril-macos-build)
