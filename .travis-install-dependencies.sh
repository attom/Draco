#!/bin/bash
##---------------------------------------------------------------------------##
## File  : ./travis-install-dependencies.sh
## Date  : Tuesday, Sep 20, 2016, 11:50 am
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2018, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# Install tools and libraries that the Draco build system requires that are
# not provided by travis or by apt-get.

# preliminaries and environment
set -e
source regression/scripts/common.sh

topdir=`pwd` # /home/travis/build/lanl/Draco

# Return integer > 0 if 'develop' branch is found.
function find_dev_branch
{
  set -f
  git branch -a | grep -c develop
  set +f
}

# printenv

if [[ ${STYLE} ]]; then

  # Ensure the 'develop' branch is available.  In some cases (merge a branch
  # that lives at github.com/lanl), the develop branch is missing in the
  # travis checkout. Since we only test files that are modified when comapred to
  # the 'develop' branch, the develop branch must be available locally.
  num_dev_branches_found=`find_dev_branch`
  if [[ $num_dev_branches_found == 0 ]]; then
    echo "no develop branches found."
    # Register the develop branch in draco/.git/config
    run "git config --local remote.origin.fetch +refs/heads/develop:refs/remotes/origin/develop"
    # Download the meta-data for the 'develop' branch
    run "git fetch"
    # Create a local tracking branch
    run "git branch -t develop origin/develop"
  fi

  # clang-format and git-clang-format
  # https://blog.kowalczyk.info/article/k/how-to-install-latest-clang-6.0-on-ubuntu-16.04-xenial-wsl.html
  echo " "
  echo "Clang-format"
  run "wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | sudo apt-key add -"
  run "sudo add-apt-repository 'deb http://apt.llvm.org/${DISTRO}/ llvm-toolchain-${DISTRO}-${CLANG_FORMAT_VER} main'"
  run "sudo apt-get update"
  run "sudo apt-get install -y clang-format-${CLANG_FORMAT_VER}"
  run "cd ${VENDOR_DIR}/bin"
  if [[ -x /usr/bin/clang-format-${CLANG_FORMAT_VER} ]]; then
    run "sudo ln -s /usr/bin/clang-format-${CLANG_FORMAT_VER} clang-format"
    run "sudo ln -s /usr/bin/git-clang-format-${CLANG_FORMAT_VER} git-clang-format"
  else
    die "Didn't find /usr/bin/clang-format-${CLANG_FORMAT_VER}"
  fi
  run "cd $topdir"

else

  # Random123
  echo " "
  echo "Random123"
  run "cd $VENDOR_DIR"
  run "wget -q http://www.deshawresearch.com/downloads/download_random123.cgi/Random123-${RANDOM123_VER}.tar.gz"
  run "tar -xvf Random123-${RANDOM123_VER}.tar.gz &> build-r123.log"
  echo "Please set RANDOM123_INC_DIR=$VENDOR_DIR/Random123-${RANDOM123_VER}/include"
  run "ls $VENDOR_DIR/Random123-${RANDOM123_VER}/include"

  # CMake
  echo " "
  echo "CMake"
  run "cd $VENDOR_DIR"
  run "wget -q --no-check-certificate http://www.cmake.org/files/v${CMAKE_VER:0:3}/cmake-${CMAKE_VER}.tar.gz"
  run "tar -xzf cmake-${CMAKE_VER}.tar.gz &> build-cmake.log"

  # Numdiff
  echo " "
  echo "Numdiff"
  run "cd $VENDOR_DIR"
  run "wget -q http://mirror.lihnidos.org/GNU/savannah/numdiff/numdiff-${NUMDIFF_VER}.tar.gz"
  run "tar -xvf numdiff-${NUMDIFF_VER}.tar.gz >& build-numdiff.log"
  run "mkdir numdiff-build"
  run "cd numdiff-build"
  run "../numdiff-${NUMDIFF_VER}/configure --prefix=${VENDOR_DIR}/numdiff-${NUMDIFF_VER} && make >> build-numdiff.log 2>&1"
  run "make -j 4 install"
  run "cd $topdir"

  # GSL
  echo " "
  echo "GSL"
  run "cd $VENDOR_DIR"
  run "wget -q http://mirror.switch.ch/ftp/mirror/gnu/gsl/gsl-${GSL_VER}.tar.gz"
  run "tar -xvf gsl-${GSL_VER}.tar.gz &> build-gsl.log"
  run "mkdir -p gsl-build"
  run "cd gsl-build"
  run "CC=gcc-${GCCVER} CFLAGS=-fPIC ../gsl-${GSL_VER}/configure --with-pic --enable-static --disable-shared --prefix=${VENDOR_DIR}/gsl-${GSL_VER} >> build-gsl.log 2>&1"
  run "make -j 4"
  run "make install"
  run "cd $topdir"

  # # OpenMPI
  # echo " "
  # echo "OpenMPI"
  # run "wget -q --no-check-certificate https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-${OPENMPI_VER}.tar.gz"
  # run "tar -zxf openmpi-${OPENMPI_VER}.tar.gz > build-openmpi.log"
  # run "cd openmpi-${OPENMPI_VER}"
  # run "./configure --enable-mpi-thread-multiple --quiet >> build-openmpi.log 2>&1"
  # # run "travis_wait 20 make"
  # run "make >> build-openmpi.log 2>&1"
  # run "sudo make install"
  # run "sudo sh -c 'echo \"/usr/local/lib\n/usr/local/lib/openmpi\" > /etc/ld.so.conf.d/openmpi.conf'"
  # run "sudo ldconfig"
  # run "cd $topdir"

  echo "Done building and/or installing TPLs"

fi

# Finish up and report
