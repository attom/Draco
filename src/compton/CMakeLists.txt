#-----------------------------*-cmake-*----------------------------------------#
# file   src/compton/CMakeLists.txt
# author Kendra Keady
# date   2017 February 28
# brief  Instructions for building compton Makefiles.
# note   Copyright (C) 2017-2018 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( compton CXX )

# This packages requires libcompton.  See draco/config/vendor_libraries.cmake
# and draco/config/FindCOMPTON.cmake to manage how the build system discovers
# this TPL.  If not found, only a stub header file will be installed.

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

configure_file( config.h.in ${PROJECT_BINARY_DIR}/compton/config.h )

set( sources "${PROJECT_SOURCE_DIR}/Compton.cc" )
set( headers
  ${PROJECT_BINARY_DIR}/compton/config.h
  ${PROJECT_SOURCE_DIR}/Compton.hh )

if( TARGET COMPTON::compton )

  # Make the header files available in the IDE.
  if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
    list( APPEND sources ${headers} )
  endif()

  # -------------------------------------------------------------------------- #
  # Build package library
  # -------------------------------------------------------------------------- #

  add_component_library(
    TARGET       Lib_compton
    TARGET_DEPS  "Lib_c4;Lib_dsxx;COMPTON::compton"
    LIBRARY_NAME compton
    SOURCES      "${sources}" )
  # generated include directive files (config.h)
  target_include_directories( Lib_compton
    PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

  # -------------------------------------------------------------------------- #
  # Installation instructions
  # -------------------------------------------------------------------------- #

  install( TARGETS Lib_compton EXPORT draco-targets DESTINATION
    ${DBSCFGDIR}lib )
  install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/compton )

  # -------------------------------------------------------------------------- #
  # Unit tests
  # -------------------------------------------------------------------------- #

  if( BUILD_TESTING )
    add_subdirectory( test )
  endif()

  # -------------------------------------------------------------------------- #
  # Autodoc
  # -------------------------------------------------------------------------- #

  process_autodoc_pages()

else()

  # If no compton library, only install headers (stubs).
  install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/compton )

endif()

#------------------------------------------------------------------------------#
# End compton/CMakeLists.txt
#------------------------------------------------------------------------------#
