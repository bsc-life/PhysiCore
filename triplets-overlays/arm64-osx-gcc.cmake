set(VCPKG_TARGET_ARCHITECTURE arm64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)

set(VCPKG_CMAKE_SYSTEM_NAME Darwin)
set(VCPKG_OSX_ARCHITECTURES arm64)

set(VCPKG_C_COMPILER "/opt/homebrew/bin/gcc-15")
set(VCPKG_CXX_COMPILER "/opt/homebrew/bin/g++-15")

set(VCPKG_CMAKE_CONFIGURE_OPTIONS "-DCMAKE_C_COMPILER=${VCPKG_C_COMPILER}"
                                  "-DCMAKE_CXX_COMPILER=${VCPKG_CXX_COMPILER}")

set(VCPKG_C_FLAGS "-fcf-protection=none")
set(VCPKG_CXX_FLAGS "-fcf-protection=none")
