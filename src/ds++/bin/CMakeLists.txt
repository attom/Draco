#-----------------------------*-cmake-*----------------------------------------#
# file   ds++/bin/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 Nov 6
# brief  Provide.
# note   Copyright (C) 2016-2018, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( dsxx_bin CXX )

if( WIN32 )

# ---------------------------------------------------------------------------- #
# Build binaries
# ---------------------------------------------------------------------------- #

add_component_executable(
  TARGET      Exe_testCrashHandler
  TARGET_DEPS Lib_dsxx
  SOURCES     ${PROJECT_SOURCE_DIR}/testCrashHandler.cc
  PREFIX       Draco
  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Exe_testCrashHandler EXPORT draco-targets
  DESTINATION ${DBSCFGDIR}bin )

endif() # WIN32

# ---------------------------------------------------------------------------- #
# end ds++/bin/CMakeLists.txt
# ---------------------------------------------------------------------------- #
