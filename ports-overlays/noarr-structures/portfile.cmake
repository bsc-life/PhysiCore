vcpkg_from_git(
  OUT_SOURCE_PATH SOURCE_PATH URL
  https://github.com/ParaCoToUl/noarr-structures.git REF
  4bdf239b6e407e9c939b80e7718da2009962bac4)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})

vcpkg_cmake_install()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
