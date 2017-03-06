//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/Compton_NWA.hh
 * \author Kendra Keady
 * \date   Mon Feb 27 2017
 * \brief  Header file for compton NWA interface
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __compton_Compton_NWA_hh__
#define __compton_Compton_NWA_hh__

// C++ standard library dependencies
#include <iostream>
#include <memory>
#include <vector>
// headers provided in compton NWA include directory
#include "etemp_interp.hh"
#include "multigroup_compton_data.hh"

namespace rtt_compton {
//===========================================================================//
/*!
 * \class Compton_NWA
 *
 * \brief Provides access to relativistic Compton scattering angle and
 * multigroup frequency distributions from the NWA project.
 *
 * This interface class allows the client to:
 * 1) access (interpolate) data from existing multigroup NWA libraries
 * 2) build new multigroup libraries from existing NWA pointwise libraries
 * 3) obtain auxiliary information for existing multigroup libraries
 *    (electron temperature bounds, frequency group structures, etc)
 *
 * This class is designed to be used with the NWA library and headers. If this
 * are not found at the CMake configure step, the lib_compton portion of draco
 * will not be build
 *
 * <b>User's environment</b>
 *
 * Cmake searches for the NWA library/include headers during the configuration
 * step. The script that does this is located at:
 *
 * /draco/config/findNWA.cmake
 *
 * TODO: The NWA libs can be built with or without openMP. We need to decide
 * which will be released (or if both will exist on a system), and ensure
 * draco links to the correct one based on the preproc macro OPENMP_FOUND.
 */

/*!
 * \example compton/test/tCompton_NWA.cc
 *
 * This unit test demonstrates the two methods for constructing a Compton_NWA
 * object, and exercises all routines for interpolation and data access.
*/

class Compton_NWA {

private:
  //! Shared pointer to an electron interpolation object:
  std::shared_ptr<etemp_interp> ei;

public:
  //! Constructor for an existing multigroup library
  explicit Compton_NWA(const std::string &filehandle);

  //! Constructor to build a multigroup library from an existing pointwise file
  Compton_NWA(const std::string &file, const std::vector<double> &group_bounds,
              const size_t n_xi);

  //! Interpolation of all data to a certain electron temperature:
  std::vector<std::vector<std::vector<double>>> interpolate(const double etemp);

  //! Retrieve group structure for the given library:
  std::vector<double> get_group_bounds() {
    return ei->get_Cdata()->get_group_bds();
  }

  //! Retrieve min electron temperature for the given library:
  double get_min_etemp() { return ei->get_min_etemp(); }

  //! Retrieve max electron temperature for the given library:
  double get_max_etemp() { return ei->get_max_etemp(); }

  //! Retrieve number of groups in the given multigroup structure:
  size_t get_num_groups() { return ei->get_Cdata()->get_n_grps(); }

  //! Retrieve number of angular moments/evaluation points in the lib data:
  size_t get_num_xi() { return ei->get_Cdata()->get_n_xi_pts(); }

  //! Retrieve electron temperature eval points (diagnostic use)
  std::vector<double> get_etemp_pts() { return ei->get_etemp_pts(); }
};
}

#endif
