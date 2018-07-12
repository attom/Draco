//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstofpstream.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::determinate_swap and c4::indeterminate_swap functions
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/ofpstream.hh"
#include "ds++/Release.hh"
#include <cmath>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstofpstream(UnitTest &ut) {

  unsigned const pid = rtt_c4::node();

  ofpstream out("tstofpstream.txt");

  out << "MPI rank " << pid << " reporting ..." << endl;

  out.send();

  out.shrink_to_fit();

  out << "MPI rank " << pid << " reporting a second time ..." << endl;

  out.shrink_to_fit();
  out.send();

  ut.passes("completed serialized write without hanging or segfaulting");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstofpstream(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstofpstream.cc
//---------------------------------------------------------------------------//
