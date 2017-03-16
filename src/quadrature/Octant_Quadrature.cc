//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Octant_Quadrature.cc
 * \author Kent Budge
 * \date   Friday, Nov 30, 2012, 08:27 am
 * \brief  Implementation for Octant_Quadrature
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id: Octant_Quadrature.cc 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include "Octant_Quadrature.hh"

#include "ds++/Soft_Equivalence.hh"
#include "ds++/to_string.hh"
#include "units/PhysicalConstants.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace rtt_quadrature {
using namespace rtt_dsxx;

//---------------------------------------------------------------------------------------//

bool Octant_Quadrature::has_axis_assignments() const {
  return has_axis_assignments_;
}


//---------------------------------------------------------------------------------------//
unsigned
Octant_Quadrature::create_aligned_ordinates_(vector<double>& mu,
                                             vector<double>& eta,
                                             vector<double>& wt) const {
  // Find axis-aligned ordinates, allowing for degenerates and replicates
  // Create and reorder mu,eta,wt starting with these ordinates in 3-D

  unsigned const octantOrdinates = mu.size();
  Check(octantOrdinates > 0);
  Check(eta.size() == octantOrdinates);
  Check(wt.size() == octantOrdinates);

  // First, we find the weight of these ordinates, combining equivalents
  double wtPolarAligned = 0.0;
  double wtAzimuthalAligned = 0.0;
  unsigned numUnaligned = 0;
  for (size_t n=0; n < octantOrdinates; ++n) {
    // TODO: change abs to std::fabs, just in case
    if ((abs(mu[n]) == 1.0) || (abs(eta[n]) == 1.0)) {
      // Weights for directions (1,0,0) and (0,1,0) will be combined
      wtAzimuthalAligned += wt[n];
    }
    else if ((abs(mu[n]) == 0.0) && (abs(eta[n]) == 0.0)) {
      wtPolarAligned += wt[n];
    }
    else {
      numUnaligned++;
    }
  }

  // If any ordinates are aligned, we reorder and add the 3-D aligned ordinates
  unsigned const numAligned = ((wtPolarAligned > 0.0) ? 2 : 0) + ((wtAzimuthalAligned > 0.0) ? 4 : 0);
  unsigned const octantAndAlignedOrdinates = numUnaligned + numAligned;
  Check(numUnaligned <= octantOrdinates);
  if (numUnaligned != octantOrdinates) {
    vector<double> muCopy(octantAndAlignedOrdinates, 0.0);
    vector<double> etaCopy(octantAndAlignedOrdinates, 0.0);
    vector<double> wtCopy(octantAndAlignedOrdinates, 0.0);
    size_t ordinate = 0;
    // Add all 3-D polar aligned ordinates
    if (wtPolarAligned > 0.0) {
      // direction = (0,0,+1)
      muCopy[ordinate] = 0.0;
      etaCopy[ordinate] = 0.0;
      wtCopy[ordinate] = wtPolarAligned;
      ++ordinate;
      // direction = (0,0,-1)
      muCopy[ordinate] = 0.0;
      etaCopy[ordinate] = 0.0;
      wtCopy[ordinate] = wtPolarAligned;
      ++ordinate;
    }
    // Add all 3-D azimuthal aligned ordinates
    if (wtAzimuthalAligned > 0.0) {
      // direction = (+1,0,0)
      muCopy[ordinate] = 1.0;
      etaCopy[ordinate] = 0.0;
      wtCopy[ordinate] = wtAzimuthalAligned;
      ++ordinate;
      // direction = (0,1,0)
      muCopy[ordinate] = 0.0;
      etaCopy[ordinate] = 1.0;
      wtCopy[ordinate] = wtAzimuthalAligned;
      ++ordinate;
      // direction = (-1,0,0)
      muCopy[ordinate] = -1.0;
      etaCopy[ordinate] = 0.0;
      wtCopy[ordinate] = wtAzimuthalAligned;
      ++ordinate;
      // direction = (0,-1,0)
      muCopy[ordinate] = 0.0;
      etaCopy[ordinate] = -1.0;
      wtCopy[ordinate] = wtAzimuthalAligned;
      ++ordinate;
    }
    // Add unaligned ordinates in the octant
    for (size_t n=0; n < octantOrdinates; ++n) {
      if (!((abs(mu[n]) == 1.0) || (abs(eta[n]) == 1.0) ||
            (abs(mu[n]) == 0.0) || (abs(eta[n]) == 0.0))) {
        Check(n < mu.size() && n < eta.size() && n < wt.size());
        Check(ordinate < muCopy.size() && ordinate < etaCopy.size() &&
              ordinate < wtCopy.size());
        muCopy[ordinate] = mu[n];
        etaCopy[ordinate] = eta[n];
        wtCopy[ordinate] = wt[n];
        ++ordinate;
      }
    }
    mu.swap(muCopy);
    eta.swap(etaCopy);
    wt.swap(wtCopy);
  }

  return numAligned;
}

//---------------------------------------------------------------------------------------//
vector<Ordinate> Octant_Quadrature::create_ordinates_(
    unsigned const dimension, Geometry const geometry, double const norm,
    unsigned const mu_axis, unsigned const eta_axis,
    bool const include_starting_directions,
    bool const include_extra_directions) const {
  using rtt_dsxx::soft_equiv;

  // We build the 3-D first, then edit as appropriate.

  vector<double> mu, eta, wt;

  create_octant_ordinates_(mu, eta, wt);

  unsigned const numAligned = create_aligned_ordinates_(mu, eta, wt);
  bool const hasAlignedPolar = (numAligned == 2 || numAligned == 6);
  bool const hasAlignedAzimuthal = (numAligned == 4 || numAligned == 6);

  unsigned const octantOrdinates = mu.size();
  unsigned const numUnaligned = octantOrdinates - numAligned;
  Check(numAligned <= octantOrdinates);
  Check(octantOrdinates > 0);
  Check(eta.size() == octantOrdinates);
  Check(wt.size() == octantOrdinates);

  unsigned const numOrdinates = numUnaligned * 8 + numAligned;
  mu.resize(numOrdinates);
  eta.resize(numOrdinates);
  wt.resize(numOrdinates);

  // Evaluate unaligned mu and eta for octants 2-4
  for (size_t octant = 2; octant <= 4; ++octant)
    for (size_t i = 0; i < numUnaligned; ++i) {
      size_t const n = i + numAligned;
      size_t const m = (octant - 1) * numUnaligned + n;
      Check(n < mu.size() && n < eta.size() && n < wt.size());
      Check(m < mu.size() && m < eta.size() && m < wt.size());
      switch (octant) {
      case 2:
        mu[m] = -mu[n];
        eta[m] = eta[n];
        wt[m] = wt[n];
        break;

      case 3:
        mu[m] = -mu[n];
        eta[m] = -eta[n];
        wt[m] = wt[n];
        break;

      case 4:
        mu[m] = mu[n];
        eta[m] = -eta[n];
        wt[m] = wt[n];
        break;
      default:
        Insist(false, "Octant value should only be 2, 3 or 4 in this loop.");
        break;
      }
    }

  // Evaluate unaligned mu and eta for octants 5-8
  for (size_t i = 0; i < 4 * numUnaligned; ++i) {
    size_t const n = i + numAligned;
    size_t const m = n + 4 * numUnaligned;
    Check(n < mu.size() && n < eta.size() && n < wt.size());
    Check(m < mu.size() && m < eta.size() && m < wt.size());
    mu[m] = mu[n];
    eta[m] = eta[n];
    wt[m] = wt[n];
  }

  // Evaluate xi for all octants
  vector<double> xi(numOrdinates);
  {
      size_t n = 0;
      if (hasAlignedPolar) {
          xi[n++] = 1.0;
          xi[n++] = -1.0;
      }
      if (hasAlignedAzimuthal) {
          xi[n++] = 0.0;
          xi[n++] = 0.0;
          xi[n++] = 0.0;
          xi[n++] = 0.0;
      }
  }
  // Octants 2-4, unaligned ordinates
  for (size_t i = 0; i < 4 * numUnaligned; ++i) {
    size_t const n = i + numAligned;
    Check(n < xi.size());
    xi[n] = std::sqrt(1.0 - (mu[n] * mu[n] + eta[n] * eta[n]));
  }
  // Octants 5-8, unaligned ordinates
  for (size_t i = 0; i < 4 * numUnaligned; ++i) {
    size_t const n = i + numAligned;
    size_t const m = n + 4 * numUnaligned;
    Check(n < xi.size());
    Check(m < xi.size());
    xi[m] = -xi[n];
  }

  vector<Ordinate> Result;

  if (dimension == 3) {
    map_axes_(mu_axis, eta_axis, mu, eta, xi);

    // Copy all ordinates into the result vector.
    Result.resize(numOrdinates);
    std::cout << "___ Full sphere\n";
    for (size_t i = 0; i < numOrdinates; ++i) {
      Result[i] = Ordinate(mu[i], eta[i], xi[i], wt[i]);
    }
  } else if (dimension == 2 || geometry == rtt_mesh_element::AXISYMMETRIC) {
    map_axes_(mu_axis, eta_axis, mu, eta, xi);

    // Copy the half-sphere
    unsigned const upperBoundNumOrdinates = 4 * numUnaligned + numAligned;
    Result.resize(upperBoundNumOrdinates);
    unsigned m = 0;
    if (dimension == 1) std::cout << "___ Quarter sphere (axisymmetric)\n";
    else std::cout << "___ Half sphere\n";
    for (size_t i = 0; i < numOrdinates; ++i) {
      if (xi[i] >= 0.0 && (dimension > 1 || eta[i] >= 0.0)) {
        // degeneracy: reduction in number of ordinates from 3-D to half-sphere
        // for a type of ordinates
        int degeneracy;
        if (abs(mu[i]) == 1.0) degeneracy = 1;
        else if (abs(xi[i]) == 1.0) degeneracy = 2;
        else if (dimension == 1 && abs(eta[i]) == 1.0) degeneracy = 2;
        else if (abs(eta[i]) == 1.0) degeneracy = 1;
        else if (dimension == 1) degeneracy = 4; // dimension == 1 and unaligned
        else degeneracy = 2; // dimension == 2 and unaligned
        Check(m < Result.size());
        Result[m++] = Ordinate(mu[i], eta[i], xi[i], wt[i] * degeneracy);
      }
    }
    Result.resize(m);

    // Add starting directions if appropriate
    if (true || !hasAlignedAzimuthal) {
      add_2D_starting_directions_(geometry, include_starting_directions,
                                  include_extra_directions, Result);
    }
  } else {
    Check(dimension == 1 && geometry != rtt_mesh_element::AXISYMMETRIC);
    Check(mu_axis == 2);

    // Swap eta and xi
    map_axes_(0, 2, mu, eta, xi);

    // Only need the quarter sphere
    unsigned const upperBoundNumOrdinates = 2 * numUnaligned + numAligned;
    Result.resize(upperBoundNumOrdinates);
    unsigned m = 0;
    std::cout << "___ Quarter sphere (spherical)\n";
    for (size_t i = 0; i < numOrdinates; ++i) {
      if (mu[i] >= 0.0 && xi[i] >= 0.0) {
        int degeneracy;
        if (abs(mu[i])==1.0) degeneracy = 2;
        else if (abs(eta[i])==1.0) degeneracy = 1;
        else if (abs(xi[i])==1.0) degeneracy = 2;
        else degeneracy = 4; // unaligned
        Check(m < Result.size());
        Result[m++] = Ordinate(mu[i], eta[i], xi[i], wt[i] * degeneracy);
        //std::cout << mu[i] << " " << eta[i] << " " << xi[i] << " " << wt[i] << " " << degeneracy << "\n";
      }
    }

    // Sort
    std::sort(Result.begin(), Result.end(), Ordinate_Set::level_compare);

    // Now sum around the axis.
    m = 0;
    double eta = Result[0].eta();
    double sum = Result[0].wt();
    for (unsigned i = 1; i < upperBoundNumOrdinates; ++i) {
      double old_eta = eta;
      eta = Result[i].eta();
      if (!soft_equiv(eta, old_eta)) {
        // New level
        Result[m++] = Ordinate(old_eta, sum);
        sum = Result[i].wt();
      } else {
        // Still on old level
        sum += Result[i].wt();
      }
    }
    // Final level
    Result[m++] = Ordinate(eta, sum);
    Result.resize(m);

    // Add starting directions if appropriate
    if (!hasAlignedPolar) {
      add_1D_starting_directions_(geometry, include_starting_directions,
                                  include_extra_directions, Result);
    }
  }

  unsigned const numOrdinatesFinal = Result.size();

  // Normalize the quadrature set
  double wsum = 0.0;
  for (size_t n = 0; n <= numOrdinatesFinal - 1; ++n)
    wsum = wsum + Result[n].wt();

  if (dimension == 1 && geometry != rtt_mesh_element::AXISYMMETRIC) {
    for (size_t n = 0; n <= numOrdinatesFinal - 1; ++n)
      Result[n] = Ordinate(Result[n].mu(), Result[n].wt() * (norm / wsum));
  } else {
    for (size_t n = 0; n <= numOrdinatesFinal - 1; ++n)
      Result[n] = Ordinate(Result[n].mu(), Result[n].eta(), Result[n].xi(),
                           Result[n].wt() * (norm / wsum));
  }
  return Result;
}

//---------------------------------------------------------------------------------------//
vector<Ordinate> Octant_Quadrature::create_ordinates_(
    unsigned dimension, Geometry geometry, double norm,
    bool include_starting_directions, bool include_extra_directions) const {
  unsigned mu_axis(0), eta_axis(0);
  if (has_axis_assignments_) {
    mu_axis = mu_axis_;
    eta_axis = eta_axis_;
  } else {
    switch (dimension) {
    case 1:
      switch (geometry) {
      case rtt_mesh_element::AXISYMMETRIC:
        mu_axis = 0;
        eta_axis = 2;
        break;

      default:
        mu_axis = 2;
        eta_axis = 1;
        break;
      }
      break;

    case 2:
      switch (geometry) {
      case rtt_mesh_element::AXISYMMETRIC:
        mu_axis = 0;
        eta_axis = 2;
        break;

      default:
        mu_axis = 0;
        eta_axis = 1;
        break;
      }
      break;

    case 3:
      mu_axis = 0;
      eta_axis = 1;
      break;

    default:
      Insist(false, "bad case");
    }
  }
  return create_ordinates_(dimension, geometry, norm, mu_axis, eta_axis,
                           include_starting_directions,
                           include_extra_directions);
}

//---------------------------------------------------------------------------------------//
/*!
 * Pure virtual used in conjuction with child implementations, for common
 * features.
 */

string Octant_Quadrature::as_text(string const &indent) const {
  string Result;

  if (has_axis_assignments_) {
    Result += indent + "  axis assignments, mu = " + to_string(mu_axis_) +
              " eta = " + to_string(eta_axis_);
  }

  Result += indent + "end";

  return Result;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                 end of Octant_Quadrature.cc
//---------------------------------------------------------------------------------------//
