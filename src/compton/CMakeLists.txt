#-----------------------------*-cmake-*----------------------------------------#
# file   src/compton/CMakeLists.txt
# author Kendra Keady
# date   2017 February 28
# brief  Instructions for building compton Makefiles.
# note   Copyright (C) 2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( compton CXX )

##---------------------------------------------------------------------------##
# Requires NWA
##---------------------------------------------------------------------------##
find_package( NWA QUIET )
set_package_properties( NWA PROPERTIES
   DESCRIPTION "Access multigroup Compton scattering data."
   TYPE OPTIONAL
   PURPOSE "Required for bulding the compton component."
   )

if( NWA_FOUND )
# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # component headers
   ${PROJECT_BINARY_DIR}      # compton/config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ and compton header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
   ${ode_BINARY_DIR}          # ode/config.h
   ${NWA_INCLUDE_DIR}
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_compton
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME compton
   SOURCES      "${sources}"
   VENDOR_LIST "NWA"
   VENDOR_LIBS "${NWA_LIBRARY}"
   VENDOR_INCLUDE_DIRS "${NWA_INCLUDE_DIR}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_compton  EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/compton )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
  add_subdirectory( test )
endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

endif() # NWA_FOUND
