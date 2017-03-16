//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Product_Chebyshev_Lobatto.cc
 * \author Andrew T. Till
 * \date   Mon Mar  6 10:19:52 2017
 * \brief  A class for Product Chebyshev-Gauss-Lobatto quadrature sets.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Product_Chebyshev_Lobatto.hh"
#include "Lobatto.hh"

#include "ds++/to_string.hh"
#include "units/MathConstants.hh"

// TODO: RM
#include <iostream>

namespace rtt_quadrature {
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
string Product_Chebyshev_Lobatto::name() const {
  return "Product Chebyshev Lobatto";
}

//---------------------------------------------------------------------------//
string Product_Chebyshev_Lobatto::parse_name() const { return "product clo"; }

//---------------------------------------------------------------------------//
Quadrature_Class Product_Chebyshev_Lobatto::quadrature_class() const {
  return OCTANT_QUADRATURE;
}

//---------------------------------------------------------------------------//
unsigned Product_Chebyshev_Lobatto::number_of_levels() const {
  return sn_order_;
}

//---------------------------------------------------------------------------//
string Product_Chebyshev_Lobatto::as_text(string const &indent) const {
  string Result = indent + "type = " + parse_name() + indent + "  order = " +
                  to_string(sn_order()) + " " + to_string(azimuthal_order_) +
                  Octant_Quadrature::as_text(indent);

  return Result;
}

//---------------------------------------------------------------------------------------//
bool Product_Chebyshev_Lobatto::is_open_interval() const {
  // Lobatto is one of our few closed interval quadratures.
  return false;
}

//---------------------------------------------------------------------------//
void Product_Chebyshev_Lobatto::create_octant_ordinates_(
    vector<double> &mu, vector<double> &eta, vector<double> &wt) const {
  using std::fabs;
  using std::sqrt;
  using std::cos;
  using rtt_dsxx::soft_equiv;

  // The number of quadrature levels is equal to the requested SN order.
  size_t levels = sn_order();

  // We build one octant of the quadrature
  size_t numOrdinates = (levels - 2) * azimuthal_order_ / 4 + 1;

  // Force the direction vectors to be the correct length.
  mu.resize(numOrdinates);
  eta.resize(numOrdinates);
  wt.resize(numOrdinates);

  std::shared_ptr<Lobatto> GLo(new Lobatto(sn_order_));

  // NOTE: this aligns the gauss points with the x-axis (r-axis in cylindrical
  // coords)

  // The first ordinate points in the polar direction
  mu[0] = 0.0;
  eta[0] = 0.0;
  wt[0] = GLo->wt(0);
  //RM
  {
      const unsigned ij = 0;
      const double xi = sqrt(1.0 - mu[ij] * mu[ij] - eta[ij] * eta[ij]);
      std::cout << "ATT DBG " << ij << " " << mu[ij] << " " << eta[ij] << " " << xi << " " << wt[ij] << "\n";
  }

  unsigned icount = 1;

  for (unsigned i = 1; i < levels / 2; ++i) {
    double xmu = GLo->mu(i);
    double xwt = GLo->wt(i);
    double xsr = sqrt(1.0 - xmu * xmu);

    // RM
    //std::cout << "ATT DBG " << i << " " << xmu << " " << xwt << " " << xsr << "\n";

    for (unsigned j = 0; j < azimuthal_order_ / 2; ++j) {
      unsigned ordinate = icount;

      mu[ordinate] =
          xsr * cos(rtt_units::PI * (2.0 * j + 1.0) / azimuthal_order_ / 2.0);
      eta[ordinate] =
          xsr * sin(rtt_units::PI * (2.0 * j + 1.0) / azimuthal_order_ / 2.0);
      wt[ordinate] = xwt / (2.0 * azimuthal_order_);

      //RM
      const unsigned ij = ordinate;
      const double xi = sqrt(1.0 - mu[ij] * mu[ij] - eta[ij] * eta[ij]);
      std::cout << "ATT DBG " << ij << " " << mu[ij] << " " << eta[ij] << " " << xi << " " << wt[ij] << "\n";

      ++icount;
    }
  }
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
// end of Product_Chebyshev_Lobatto.cc
//---------------------------------------------------------------------------//
