//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Product_Chebyshev_Lobatto__Parser.cc
 * \author Andrew T. Till
 * \date   Mon Mar 6 10:21:50 2017
 * \brief  A class representing an product Chebyshev-Lobatto quadrature set.
 * \note   Copyright 2016-2017 Los Alamos National Security, LLC. All rights
 *         reserved.  */
//---------------------------------------------------------------------------//

#include "Product_Chebyshev_Lobatto.hh"
#include "parser/utilities.hh"
#include <iostream>

namespace rtt_quadrature {
using namespace rtt_parser;

//---------------------------------------------------------------------------//
std::shared_ptr<Quadrature>
Product_Chebyshev_Lobatto::parse(Token_Stream &tokens) {
  // Takes two numbers, first the number of Gauss-Lobatto points, which is the
  // SN order
  Token token = tokens.shift();
  tokens.check_syntax(token.text() == "order", "expected an order");

  unsigned sn_order = parse_positive_integer(tokens);
  tokens.check_semantics(sn_order > 2,
                         "polar order must be greater than two");
  tokens.check_semantics(sn_order % 2 == 0, "polar order must be even");

  // The second number is the number of azimuthal points on each level
  // corresponding to the Gauss-Lobatto points

  unsigned azimuthal_order = parse_positive_integer(tokens);
  tokens.check_semantics(azimuthal_order > 0,
                         "azimuthal order must be greater than zero");
  tokens.check_semantics(azimuthal_order % 2 == 0,
                         "azimuthal order must be even");

  bool has_axis_assignments;
  unsigned mu_axis, eta_axis;
  Octant_Quadrature::parse(tokens, has_axis_assignments, mu_axis, eta_axis);

  if (has_axis_assignments)
    return std::shared_ptr<Quadrature>(new Product_Chebyshev_Lobatto(
        sn_order, azimuthal_order, mu_axis, eta_axis));
  else
    return std::shared_ptr<Quadrature>(
        new Product_Chebyshev_Lobatto(sn_order, azimuthal_order));
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
// end of quadrature/Product_Chebyshev_Lobatto.cc
//---------------------------------------------------------------------------//
