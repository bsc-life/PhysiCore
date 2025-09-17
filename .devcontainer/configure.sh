#!/bin/bash

# Install LLVM
sudo bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

# Install vcpkg
cd /usr/local
sudo git clone https://github.com/microsoft/vcpkg.git
sudo chown -R $USER /usr/local/vcpkg
cd vcpkg
./bootstrap-vcpkg.sh
