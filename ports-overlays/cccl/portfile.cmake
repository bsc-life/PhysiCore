vcpkg_from_git(
  OUT_SOURCE_PATH SOURCE_PATH URL https://github.com/NVIDIA/cccl.git REF
  0fc280aa5ba7b8c780294b1c9906a25c6753387b)

vcpkg_cmake_configure(
  SOURCE_PATH
  ${SOURCE_PATH}
  OPTIONS
  -DCCCL_ENABLE_UNSTABLE=OFF
  -DCCCL_ENABLE_LIBCUDACXX=OFF
  -DCCCL_ENABLE_CUB=OFF
  -DCCCL_ENABLE_THRUST=OFF
  -DCCCL_ENABLE_CUDAX=OFF
  -DCCCL_ENABLE_TESTING=OFF
  -DCCCL_ENABLE_EXAMPLES=OFF
  -Dlibcudacxx_ENABLE_INSTALL_RULES=ON
  -DCUB_ENABLE_INSTALL_RULES=ON
  -DThrust_ENABLE_INSTALL_RULES=ON)

vcpkg_cmake_install()

file(GLOB_RECURSE HEADER_SEARCH_FILES
     "${CURRENT_PACKAGES_DIR}/lib/cmake/*/*-header-search.cmake")

# CCCL expects that the cmake files will stay under lib/cmake/PORT, but we will
# be moving them to share/PORT. So we need to change the relative paths set in
# the header search files.
foreach(_file IN LISTS HEADER_SEARCH_FILES)
  file(READ "${_file}" _contents)
  string(REPLACE "set(from_install_prefix \"../../../\")"
                 "set(from_install_prefix \"../../\")" _contents "${_contents}")
  file(WRITE "${_file}" "${_contents}")
endforeach()

vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/cccl PACKAGE_NAME cccl
                         DO_NOT_DELETE_PARENT_CONFIG_PATH)
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/thrust PACKAGE_NAME thrust
                         DO_NOT_DELETE_PARENT_CONFIG_PATH)
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/cub PACKAGE_NAME cub
                         DO_NOT_DELETE_PARENT_CONFIG_PATH)
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/libcudacxx PACKAGE_NAME
                         libcudacxx DO_NOT_DELETE_PARENT_CONFIG_PATH)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/lib")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
