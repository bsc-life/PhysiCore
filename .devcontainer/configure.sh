#!/bin/bash

# Install LLVM
wget https://apt.llvm.org/llvm.sh
chmod +x ./llvm.sh
sudo ./llvm.sh 20 all
rm llvm.sh

# Install vcpkg
cd /usr/local
sudo git clone https://github.com/microsoft/vcpkg.git
sudo chown -R $USER /usr/local/vcpkg
cd vcpkg
./bootstrap-vcpkg.sh
