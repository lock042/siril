#!/bin/bash

cd /Library/Developer/CommandLineTools/SDKs
sudo curl -L 'https://github.com/phracker/MacOSX-SDKs/releases/download/10.15/MacOSX10.12.sdk.tar.xz' | sudo tar -xzf -
echo 'export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX10.12.sdk' > ~/.profile
echo 'export MACOSX_DEPLOYMENT_TARGET=10.12' >> ~/.profile
