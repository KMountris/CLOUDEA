#--------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.8.0 FATAL_ERROR)
PROJECT(
    CLOUDEA VERSION 0.1.0
    LANGUAGES CXX
)

# Set project's CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Library compilation options.
option(BUILD_SHARED_LIBS "Build ${PROJECT_NAME} as shared library" ON)
option(BUILD_APPS "Build ${PROJECT_NAME} applications" ON)
option(BUILD_DOC "Build ${PROJECT_NAME} documentation" ON)

option(CME_PAPER_PROGS "Scripts for CME paper" OFF)
option(${PROJECT_NAME}_USE_CGAL "Build ${PROJECT_NAME} with CGAL libraries dependency" ON)

# General set up.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default to Release build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "${PROJECT_NAME} Build type. Options: Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
message(STATUS "${PROJECT_NAME} build type set to ${CMAKE_BUILD_TYPE}")

# Set build and installation directories
include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# Set the install prefix for the system's architecture.
if (WIN32)
  if (CMAKE_GENERATOR_PLATFORM EQUAL x64)
    set(CMAKE_INSTALL_PREFIX "C:/Program Files/${PROJECT_NAME}" CACHE PATH "Installation directory prefix for ${PROJECT_NAME}")
  else()
    set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/${PROJECT_NAME}" CACHE PATH "Installation directory prefix for ${PROJECT_NAME}")
  endif(CMAKE_GENERATOR_PLATFORM EQUAL x64)
else()
  set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory prefix for ${PROJECT_NAME}")
endif (WIN32)
message(STATUS "${PROJECT_NAME} will be installed in ${CMAKE_INSTALL_PREFIX}")

# Allow overriding the installation directories
set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME} CACHE PATH "Installation directory for ${PROJECT_NAME} libraries")
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME} CACHE PATH "Installation directory for ${PROJECT_NAME} applications")
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Installation directory for ${PROJECT_NAME} header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKEDIR CMake/${PROJECT_NAME})
else()
  set(DEF_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for ${PROJECT_NAME} CMake files")

# Make relative paths, if any, absolute.
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()
#------------------------------------------------------#

# Set correct RPATH.
if(APPLE)
  set(base_point @loader_path)
else()
  set(base_point $ORIGIN)
endif()
file(RELATIVE_PATH rel_dir 
    ${CMAKE_CURRENT_BINARY_DIR}/${INSTALL_BINDIR}
    ${CMAKE_CURRENT_BINARY_DIR}/${INSTALL_LIBDIR}
)
set(CMAKE_INSTALL_RPATH ${base_point} ${base_point}/${rel_dir})

#--- Report Installation directories to user
foreach(p LIB BIN INCLUDE CMAKE)
  file(TO_NATIVE_PATH ${INSTALL_${p}DIR} _path)
  message(STATUS "Installing ${p} components in ${_path}")
  unset(_path)
endforeach()

# To choose system
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(STATUS "Configuring on/for Linux")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  message(STATUS "Configuring on/for macOS")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  message(STATUS "Configuring on/for Windows")
else()
  message(STATUS "Configuring on/for ${CMAKE_SYSTEM_NAME}")
endif()

# Set compile options for supported compilers
if(MSVC)
    set(${PROJECT_NAME}_COMPILE_OPTIONS "/W4" "$<$<CONFIG:RELEASE>:/O2>")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(${PROJECT_NAME}_COMPILE_OPTIONS "-Wall" "-Wextra" "-Wdouble-promotion" "-ansi" "-pedantic"
    "-march=native" "-frounding-math" "-fsignaling-nans" "-fPIC" "-stdlib=libc++" "$<$<CONFIG:RELEASE>:-O3>")
else()
    set(${PROJECT_NAME}_COMPILE_OPTIONS "-Wall" "-Wextra" "-Wdouble-promotion" "-ansi" "-pedantic"
      "-march=native" "-frounding-math" "-fsignaling-nans" "$<$<CONFIG:RELEASE>:-O3>" "-fPIC" "-pthread")
endif()
add_compile_options(${${PROJECT_NAME}_COMPILE_OPTIONS})


# Compile the library components.
add_subdirectory(src)

# Uninstall custom command.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cloudea_uninstall.cmake.in
               ${CMAKE_BINARY_DIR}/cloudea_uninstall.cmake IMMEDIATE  @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cloudea_uninstall.cmake)

# Compile the applications.
if (BUILD_APPS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/apps)
endif()

# Compile the documentation.
if (BUILD_DOC)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/doc)
endif ()
