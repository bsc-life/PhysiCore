#!/bin/bash

# Install LLVM
wget https://apt.llvm.org/llvm.sh
chmod +x ./llvm.sh
sudo ./llvm.sh 20 all
rm llvm.sh

# Bootstrapping vcpkg
cd vcpkg
./bootstrap-vcpkg.sh
