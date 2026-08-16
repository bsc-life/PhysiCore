vcpkg_from_git(
  OUT_SOURCE_PATH SOURCE_PATH URL
  https://github.com/ParaCoToUl/noarr-structures.git REF
  914846c383f90e64c1c9ccf0a327170e07c3fcac)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH})

vcpkg_cmake_install()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
