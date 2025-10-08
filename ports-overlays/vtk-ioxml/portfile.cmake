set(VTK_SHORT_VERSION 9.5)
vcpkg_from_git(
  OUT_SOURCE_PATH SOURCE_PATH URL https://github.com/Kitware/VTK.git REF
  4043626a4c42c536cf4563dc50321579b99d424a)

vcpkg_cmake_configure(
  SOURCE_PATH
  "${SOURCE_PATH}"
  OPTIONS
  -DBUILD_SHARED_LIBS=OFF
  -DVTK_BUILD_ALL_MODULES=OFF
  -DVTK_ENABLE_REMOTE_MODULES=OFF
  -DVTK_MODULE_ENABLE_VTK_IOXML=YES
  -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT
  -DVTK_GROUP_ENABLE_Rendering=NO
  -DVTK_GROUP_ENABLE_Web=NO
  -DVTK_GROUP_ENABLE_MPI=NO
  -DVTK_GROUP_ENABLE_Imaging=NO
  -DVTK_GROUP_ENABLE_Views=NO
  -DVTK_GROUP_ENABLE_Qt=NO
  -DVTK_BUILD_TESTING=OFF
  -DVTK_BUILD_EXAMPLES=OFF
  -DVTK_ENABLE_LOGGING=OFF
  -DVTK_ENABLE_WRAPPING=OFF)

vcpkg_cmake_build()

vcpkg_cmake_install()

vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/vtk-${VTK_SHORT_VERSION})

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")

vcpkg_copy_tool_dependencies("${CURRENT_PACKAGES_DIR}/tools/vtk")

file(RENAME "${CURRENT_PACKAGES_DIR}/share/licenses"
     "${CURRENT_PACKAGES_DIR}/share/${PORT}/licenses")
vcpkg_install_copyright(
  FILE_LIST
  "${SOURCE_PATH}/Copyright.txt"
  COMMENT
  [[
This file presents the top-level Copyright.txt.
Additional licenses and notes are located in the licenses directory.
]])
