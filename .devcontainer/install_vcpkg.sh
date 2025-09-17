#!/bin/bash

cd /usr/local
sudo git clone https://github.com/microsoft/vcpkg.git
sudo chown -R $USER /usr/local/vcpkg
cd vcpkg
./bootstrap-vcpkg.sh
