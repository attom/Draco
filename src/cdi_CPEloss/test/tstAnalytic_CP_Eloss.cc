//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_CPEloss/test/tstAnalytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss test.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "cdi_CPEloss/Analytic_CP_Eloss.hh"
#include "cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_Spitzer_Eloss_Model.hh"
#include "cdi/CDI.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"
#include <sstream>

using namespace std;

using rtt_cdi::CDI;
using rtt_cdi_cpeloss::Analytic_CP_Eloss;
using rtt_cdi_cpeloss::Analytic_Eloss_Model;
//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

void KP_alpha_test(rtt_dsxx::UnitTest &ut) {

  // Deuterium:
  int32_t deuterium_zaid = 1002;
  double deuterium_mass = 3.34358e-24;
  rtt_cdi::CParticle target_in(deuterium_zaid, deuterium_mass);
  // Alpha particle:
  int32_t alpha_zaid = 2004;
  double alpha_mass = 6.64424e-24;
  rtt_cdi::CParticle projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> model_in(
      new rtt_cdi_cpeloss::Analytic_KP_Alpha_Eloss_Model(target_in,
                                                         projectile_in));

  Analytic_CP_Eloss eloss_mod(model_in, target_in, projectile_in,
                              rtt_cdi::CPModelAngleCutoff::NONE);

  // Check that basic accessors return correct result:
  // Analytic model should match that passed to constructor
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(
      eloss_mod.getEloss(1., 10., 1.),
      eloss_mod.get_Analytic_Model()->calculate_eloss(1., 10., 1.), 1.0e-3));

  // Model type better be analytic:
  FAIL_IF_NOT(eloss_mod.getModelType() == rtt_cdi::CPModelType::ANALYTIC_ETYPE);

  // NOT tabular data
  FAIL_IF(eloss_mod.is_data_in_tabular_form());

  // All sizes should be 0 (again, not tabular data)
  FAIL_IF_NOT(eloss_mod.getTemperatureGrid().size() ==
              eloss_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_mod.getDensityGrid().size() == eloss_mod.getNumDensities());
  FAIL_IF_NOT(eloss_mod.getEnergyGrid().size() == eloss_mod.getNumEnergies());

  // Data file name should be an empty string:
  FAIL_IF_NOT(eloss_mod.getDataFilename().empty());

  // Check that accessors return the correct target and projectile:
  FAIL_IF_NOT(target_in.get_zaid() == eloss_mod.getTarget().get_zaid());
  FAIL_IF_NOT(projectile_in.get_zaid() == eloss_mod.getProjectile().get_zaid());

  // Check that accessor returns the correct model angle cutoff
  FAIL_IF_NOT(eloss_mod.getModelAngleCutoff() ==
              rtt_cdi::CPModelAngleCutoff::NONE);

  // Get eloss values for some sample data:
  {
    double T = 1.0;    // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 212.287, 1.0e-3));
  }

  // Point near the maximum on rho in [0, 20] g/cc, T in [0, 100] keV:
  {
    double T = 15.0;   // keV
    double rho = 20.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 78.5845, 1.0e-4));
  }

  // Point where the range fit goes negative -- eloss should return 1.0e25:
  // (equivalent to the particle range limiting to a very small number)
  {
    double T = 90.0;   // keV
    double rho = 17.5; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.0e25));
  }

  if (ut.numFails == 0)
    PASSMSG("KP_Alpha CPEloss test passes.");
  else
    FAILMSG("KP_Alpha CPEloss test fails.");

  return;
}

void Spitzer_test(rtt_dsxx::UnitTest &ut) {

  // Deuterium:
  int32_t deuterium_zaid = 1002;
  double deuterium_mass = 3.34358e-24;
  rtt_cdi::CParticle target_in(deuterium_zaid, deuterium_mass);
  // Alpha particle:
  int32_t alpha_zaid = 2004;
  double alpha_mass = 6.64424e-24;
  rtt_cdi::CParticle projectile_in(alpha_zaid, alpha_mass);

  std::shared_ptr<Analytic_Eloss_Model> model_in(
      new rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model(target_in,
                                                        projectile_in));

  Analytic_CP_Eloss eloss_mod(model_in, target_in, projectile_in,
                              rtt_cdi::CPModelAngleCutoff::NONE);

  // Check that basic accessors return correct result:
  // Analytic model should match that passed to constructor
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(
      eloss_mod.getEloss(1., 10., 1.),
      eloss_mod.get_Analytic_Model()->calculate_eloss(1., 10., 1.), 1.0e-3));

  // Model type better be analytic:
  FAIL_IF_NOT(eloss_mod.getModelType() == rtt_cdi::CPModelType::ANALYTIC_ETYPE);

  // NOT tabular data
  FAIL_IF(eloss_mod.is_data_in_tabular_form());

  // All sizes should be 0 (again, not tabular data)
  FAIL_IF_NOT(eloss_mod.getTemperatureGrid().size() ==
              eloss_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_mod.getDensityGrid().size() == eloss_mod.getNumDensities());
  FAIL_IF_NOT(eloss_mod.getEnergyGrid().size() == eloss_mod.getNumEnergies());

  // Data file name should be an empty string:
  FAIL_IF_NOT(eloss_mod.getDataFilename().empty());

  // Check that accessors return the correct target and projectile:
  FAIL_IF_NOT(target_in.get_zaid() == eloss_mod.getTarget().get_zaid());
  FAIL_IF_NOT(projectile_in.get_zaid() == eloss_mod.getProjectile().get_zaid());

  // Check that accessor returns the correct model angle cutoff
  FAIL_IF_NOT(eloss_mod.getModelAngleCutoff() ==
              rtt_cdi::CPModelAngleCutoff::NONE);

  // Get eloss values for some sample data:
  {
    double T = 1.0;    // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(eloss_coeff, 9.823342835413329303e+05, 1.0e-8));
  }

  // Check for large value when 2kT > E
  {
    double T = 1.e2;   // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(eloss_coeff, 6.241509343260178341e+23, 1.0e-8));
  }

  if (ut.numFails == 0)
    PASSMSG("Spitzer CPEloss test passes.");
  else
    FAILMSG("Spitzer CPEloss test fails.");

  return;
}

//----------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    KP_alpha_test(ut);
    Spitzer_test(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// end of tstAnalytic_CP_Eloss.cc
//----------------------------------------------------------------------------//
