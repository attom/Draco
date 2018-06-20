//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/bin/xthi.cc
 * \author Mike Berry <mrberry@lanl.gov>, Kelly Thompson <kgt@lanl.gov>,
 *         Tim Kelley <tkelley@lanl.gov.
 * \date   Tuesday, Jun 05, 2018, 17:12 pm
 * \brief  Print MPI rank, thread number and core affinity bindings.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * Rewritten by Tim Kelley to run C++11 std::threads You may override
 * \c NUM_WORKERS on the compile command line.  For example to run with 4 worker
 * threads:
 *
 * \code
 * $ ./ythi 4
 * \endcode
 *
 * The default is 1 worker thread (over and above the host thread)
 */
//---------------------------------------------------------------------------//

#include "c4/C4_Functions.hh"
#include "ds++/SystemCall.hh"
#include <atomic>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

namespace rtt_c4 {

//----------------------------------------------------------------------------//
/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
std::string cpuset_to_string(void) {

  // return value;
  std::ostringstream cpuset;

  // local storage
  cpu_set_t coremask;
  (void)sched_getaffinity(0, sizeof(coremask), &coremask);

  size_t entry_made = 0;
  for (int i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, &coremask)) {
      int run = 0;
      entry_made = 1;
      for (int j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, &coremask))
          run++;
        else
          break;
      }
      if (run == 0)
        cpuset << i << ",";
      else
        cpuset << i << "-" << i + run << ",";
      i += run;
    }
  }
  return cpuset.str().substr(0, cpuset.str().length() - entry_made);
}

} // end namespace rtt_c4

/**\brief After atomic bool changes to true, print out some thread info. */
void run_thread(std::atomic<bool> &signal, std::string const &hostname,
                int const rank, int const simple_thread_id) {
  while (!signal) {
  }
  std::string cpuset = rtt_c4::cpuset_to_string();
  std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5)
            << rank << ", Thread " << std::setfill('0') << std::setw(3)
            << simple_thread_id << ", core affinity = " << cpuset << std::endl;
  return;
}

//----------------------------------------------------------------------------//
int main(int argc, char **argv) {
  size_t YTHI_NUM_WORKERS = 1;

  if (argc > 1) {
    YTHI_NUM_WORKERS = std::stoi(argv[1]);
  }

  rtt_c4::initialize(argc, argv);
  int const rank = rtt_c4::node();
  std::string const hostname = rtt_dsxx::draco_gethostname();

  std::vector<std::atomic<bool>> signals(YTHI_NUM_WORKERS);
  std::vector<std::thread> threads(YTHI_NUM_WORKERS);

  for (size_t i = 0; i < YTHI_NUM_WORKERS; ++i) {
    signals[i].store(false);
    threads[i] = std::thread(run_thread, std::ref(signals[i]),
                             std::ref(hostname), rank, i + 1);
  }
  std::string cpuset = rtt_c4::cpuset_to_string();
  int host_thread(0);
  std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5)
            << rank << ", Thread " << std::setfill('0') << std::setw(3)
            << host_thread << ", core affinity = " << cpuset << std::endl;
  for (size_t i = 0; i < YTHI_NUM_WORKERS; ++i) {
    signals[i].store(true);
    threads[i].join();
  }

  rtt_c4::finalize();
  return (0);
}

//----------------------------------------------------------------------------//
// End c4/bin/ythi.cc
//----------------------------------------------------------------------------//
