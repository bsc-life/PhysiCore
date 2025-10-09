vcpkg_from_git(
  OUT_SOURCE_PATH SOURCE_PATH URL
  https://github.com/ParaCoToUl/noarr-structures.git REF
  5f6b1bdba7720a8067a86264fdee35e0f7b1d5c2)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})

vcpkg_cmake_install()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
