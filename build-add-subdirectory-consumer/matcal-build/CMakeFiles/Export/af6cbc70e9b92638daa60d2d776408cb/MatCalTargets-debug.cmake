#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MatCal::Linalg" for configuration "Debug"
set_property(TARGET MatCal::Linalg APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MatCal::Linalg PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmatcal_linalg.a"
  )

list(APPEND _cmake_import_check_targets MatCal::Linalg )
list(APPEND _cmake_import_check_files_for_MatCal::Linalg "${_IMPORT_PREFIX}/lib/libmatcal_linalg.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
