#-----------------------------*-cmake-*----------------------------------------#
# file   mesh/test/CMakeLists.txt
# author Ryan Wollaeger <wollaeger@lanl.gov>
# date   Thursday, Jun 07, 2018, 15:44 pm
# brief  Generate build project files for mesh/test.
# note   Copyright (C) 2018, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( mesh_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_parallel_tests(
    SOURCES "${test_sources}"
    PE_LIST "1"
    DEPS    "Lib_c4;Lib_mesh" )

#------------------------------------------------------------------------------#
# End mesh/test/CMakeLists.txt
#------------------------------------------------------------------------------#
