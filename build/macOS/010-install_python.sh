#!/bin/bash

cd ~/
curl -L 'https://www.python.org/ftp/python/3.9.7/python-3.9.7-macos11.pkg' > python-3.9.7-macosx11.pkg
sudo installer -pkg python-3.9.7-macosx11.pkg -target /
cd /Applications/Python\ 3.9/
./Install\ Certificates.command
