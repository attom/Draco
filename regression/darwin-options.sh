#!/bin/bash
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/sripts/ccscs-options.sh
## Date  : Tuesday, May 16, 2017, 21:42 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
##
## Summary: Capture ccscs-specific features that plug into regression-master.sh.
##---------------------------------------------------------------------------##

# Main Options
export machine_class="darwin"
export machine_name_short="darwin"
export machine_name_long="Darwin Hetrogeneous Cluster"
#platform_extra_params="clang coverage fulldiagnostics gcc630 nr perfbench scalar static valgrind vtest"
platform_extra_params="arm gpu-volta knl fulldiagnostics perfbench power9 vtest"
pem_match=`echo $platform_extra_params | sed -e 's/[ ]/|/g'`

##---------------------------------------------------------------------------##
## End ccscs-options.sh
##---------------------------------------------------------------------------##
