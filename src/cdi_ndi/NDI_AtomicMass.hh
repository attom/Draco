//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_AtomicMass.hh
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AtomicMass class definition.
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef cdi_ndi_NDI_AtomicMass_hh
#define cdi_ndi_NDI_AtomicMass_hh

#include "ds++/Assert.hh"
#include "ds++/path.hh"
#include "units/PhysicalConstexprs.hh"
#include <string>

namespace rtt_cdi_ndi {

//============================================================================//
/*!
 * \class NDI_AtomicMass
 *
 * \brief Class for getting atomic mass weights by ZAID from NDI data using a
 *        method due to T. Saller. For more details on NDI, see
 *        https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *        Currently only multigroup data is supported, continuous energy data
 *        is probably best added through a refactor.
 * \example cdi_ndi/test/tstNDI_AtomicMass.cc
 */
//============================================================================//
class NDI_AtomicMass {

public:
  //! Constructor (default gendir path)
  NDI_AtomicMass();

  /*!
   * \brief Constructor for NDI atomic mass weight reader, using custom path to
   *        NDI gendir file.
   * \param[in] gendir_path_in path to gendir file
   */
  NDI_AtomicMass(const std::string &gendir_path_in)
      : gendir_path(gendir_path_in) {
    Require(rtt_dsxx::fileExists(gendir_path));
  }

  //! Retrieve atomic mass weight for isotope with given ZAID
  double get_amw(const int zaid) const;

private:
  //! Path to gendir file
  std::string gendir_path;

  //! Unit system
  const rtt_units::PhysicalConstexprs<rtt_units::CGS> pc;
};

} // namespace rtt_cdi_ndi

#endif // cdi_ndi_NDI_AtomicMass_hh

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_AtomicMass.hh
//----------------------------------------------------------------------------//
