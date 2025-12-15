vcpkg_from_git(
  OUT_SOURCE_PATH SOURCE_PATH URL https://github.com/NVIDIA/cccl.git REF
  db78d421a0ce92ed159400d800268574aded45e3)

vcpkg_cmake_configure(
  SOURCE_PATH
  ${SOURCE_PATH}
  OPTIONS
  -DCMAKE_CUDA_ARCHITECTURES=all-major
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

# Relocate headers under a single cccl/ prefix: include/cccl/{thrust,cub,cuda}
file(MAKE_DIRECTORY "${CURRENT_PACKAGES_DIR}/include/cccl")
set(_cccl_header_dirs cub thrust cuda)
foreach(_dir IN LISTS _cccl_header_dirs)
  if(EXISTS "${CURRENT_PACKAGES_DIR}/include/${_dir}")
    file(COPY "${CURRENT_PACKAGES_DIR}/include/${_dir}" DESTINATION "${CURRENT_PACKAGES_DIR}/include/cccl")
    file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/include/${_dir}")
  endif()
endforeach()

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
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/include/cccl/cub/detail/ptx-json")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
