vcpkg_from_github(
  OUT_SOURCE_PATH
  SOURCE_PATH
  REPO
  Kitware/VTK
  REF
  v9.5.2
  SHA512
  a65aeac33770d2249c7646b05f65c7b7dd19b7e7075fccfcc882d7fdd230f7a0d5a4a0ec36addc2b84c2b2b4841bb6d52716525c3434c69ff28dd67f3565bf0c
)

vcpkg_cmake_configure(
  SOURCE_PATH
  "${SOURCE_PATH}"
  OPTIONS
  -DBUILD_SHARED_LIBS=OFF
  -DVTK_BUILD_ALL_MODULES=OFF
  -DVTK_MODULE_ENABLE_VTK_IOXMLParser=YES
  -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT
  -DVTK_GROUP_ENABLE_Rendering=NO
  -DVTK_GROUP_ENABLE_Web=NO
  -DVTK_GROUP_ENABLE_MPI=NO
  -DVTK_GROUP_ENABLE_Imaging=NO
  -DVTK_GROUP_ENABLE_Views=NO
  -DVTK_GROUP_ENABLE_Qt=NO
  -DVTK_BUILD_TESTING=OFF
  -DVTK_BUILD_EXAMPLES=OFF
  -DVTK_ENABLE_LOGGING=OFF)

vcpkg_cmake_build()

vcpkg_cmake_install()

# Fix CMake package config
vcpkg_cmake_config_fixup(CONFIG_PATH share/vtk)

# Clean up unnecessary files file(REMOVE_RECURSE
# "${CURRENT_PACKAGES_DIR}/debug/include") file(REMOVE_RECURSE
# "${CURRENT_PACKAGES_DIR}/debug/share")

# License
file(
  INSTALL "${SOURCE_PATH}/LICENSE"
  DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}"
  RENAME copyright)
