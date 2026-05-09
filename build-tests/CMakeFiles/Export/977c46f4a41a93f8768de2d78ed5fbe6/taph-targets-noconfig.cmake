#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "taph::taph" for configuration ""
set_property(TARGET taph::taph APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(taph::taph PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib64/libtaph.a"
  )

list(APPEND _cmake_import_check_targets taph::taph )
list(APPEND _cmake_import_check_files_for_taph::taph "${_IMPORT_PREFIX}/lib64/libtaph.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
