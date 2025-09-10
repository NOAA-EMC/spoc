# Include the list of binary files in the SPOC_AUX
include(spoc_aux_binary_files.cmake)

# Decide if SPOC_AUX Binary files need to be downloaded or simply copied.
if(DEFINED ENV{SPOC_AUX_BINARY_SOURCE_DIR})
  set(DOWNLOAD_TARBALL OFF)
  set(SPOC_AUX_BIN_FIX_DIR $ENV{SPOC_AUX_BINARY_SOURCE_DIR})
else()
  set(DOWNLOAD_TARBALL ON)
  set(SPOC_AUX_BIN_FIX_DIR ${CMAKE_CURRENT_BINARY_DIR})
endif()

# Set URL, sha256sum and tar name of the file to download.
set(URL "https://ftp.emc.ncep.noaa.gov/static_files/public/spoc_aux")
set(SHA "8c40deefdff14396909010bbe4ca3e96b39b493de60a43e48e665c1af3ccc459")
set(TAR "spoc_aux-8c40de.tgz")

# Check if the tarball is already downloaded.
if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${TAR})
  file(SHA256 ${CMAKE_CURRENT_BINARY_DIR}/${TAR} SHA256SUM)
  if("${SHA256SUM}" STREQUAL "${SHA}")
    message(STATUS "${TAR} found, skipping download ...")
    set(DOWNLOAD_TARBALL OFF)
  else()
    message(STATUS "Incomplete or invalid ${TAR}, re-downloading ...")
    set(DOWNLOAD_TARBALL ON)
  endif()
endif()

if(DOWNLOAD_TARBALL)
  message(STATUS "Downloading ${TAR} from ${URL}")

  # Download binary fix files and check SHA256
  file(DOWNLOAD
    ${URL}/${TAR}
    ${CMAKE_CURRENT_BINARY_DIR}/${TAR}
    INACTIVITY_TIMEOUT 30
    TIMEOUT 600
    SHOW_PROGRESS
    STATUS _DOWNLOAD_STATUS
    EXPECTED_HASH SHA256=${SHA}
    )

  list(GET _DOWNLOAD_STATUS 0 _STATUS)
  list(GET _DOWNLOAD_STATUS 1 _MESSAGE)
  if(_STATUS)
    message(FATAL_ERROR "Error downloading ${TAR}: ${_MESSAGE}")
  endif()

  # Extract downloaded tarball.
  file(ARCHIVE_EXTRACT INPUT ${TAR})
  message(DEBUG "Download ${TAR} succeeded and extracted")
endif()

# Copy from the downloaded or central storage location
message(STATUS "Copying SPOC_AUX binary fix files from ${SPOC_AUX_BIN_FIX_DIR}")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Big_Endian)
foreach(file IN LISTS spoc_aux_binaries)
  message(DEBUG "copying ... ${file}")
  configure_file(${SPOC_AUX_BIN_FIX_DIR}/${file} ${CMAKE_CURRENT_SOURCE_DIR}/${file} COPYONLY)
  # TODO: Replace with file(COPY_FILE) once minimum CMake is >= 3.21.
  # file(COPY_FILE ${SPOC_AUX_BIN_FIX_DIR}/${file} ${CMAKE_CURRENT_SOURCE_DIR}/${file} ONLY_IF_DIFFERENT)
endforeach()
