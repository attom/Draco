//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   compton2/Compton2.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Implementation file for compton interface
 * \note   Copyright (C) 2017-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

// headers provided in draco:
#include "compton2/Compton2.hh"
#include "c4/C4_Functions.hh" //node()
#include "c4/global.hh"
#include "ds++/Assert.hh"

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <iostream>

using UINT = size_t;
using FP = double;
using vec = std::vector<FP>;

namespace rtt_compton2 {

Compton2::Compton2(std::string filename)
    : num_temperatures_(0U), num_groups_(0U), num_leg_moments_(0U),
      num_evals_(0U), Ts_(0U), Egs_(0U), first_groups_(0U), indexes_(0U),
      data_(0U), derivs_(0U) {
  Require(filename.length() > 0U);
  int rank = rtt_c4::node();
  constexpr int bcast_rank = 0;
  int errcode = 0;

  if (rank == bcast_rank) {
    errcode = read_binary(filename);
  }
  broadcast_MPI(errcode);

  Ensure(check_class_invariants());
}

void Compton2::broadcast_MPI(int errcode) {
  rtt_c4::global_max(errcode);
  Insist(errcode == 0, "Non-zero errorcode. Exiting.");

  int rank = rtt_c4::node();
  constexpr int bcast_rank = 0;

  // Broadcast sizes
  UINT data_size = data_.size();
  std::array<UINT, 5> pack = {num_temperatures_, num_groups_, num_leg_moments_,
                              num_evals_, data_size};
  rtt_c4::broadcast(&pack[0], pack.size(), bcast_rank);
  num_temperatures_ = pack[0];
  num_groups_ = pack[1];
  num_leg_moments_ = pack[2];
  num_evals_ = pack[3];
  data_size = pack[4];

  // Derived sizes
  const UINT tsz = num_temperatures_;
  const UINT egsz = num_groups_ + 1U;
  const UINT fgsz = num_temperatures_ * num_groups_;
  const UINT isz = fgsz + 1U;
  const UINT dsz = data_size;

  // Broadcast grids
  if (rank != bcast_rank)
    Ts_.resize(tsz);
  rtt_c4::broadcast(&Ts_[0], tsz, 0);
  if (rank != bcast_rank)
    Egs_.resize(egsz);
  rtt_c4::broadcast(&Egs_[0], egsz, 0);

  // Broadcast sparse data structures
  if (rank != bcast_rank)
    first_groups_.resize(fgsz);
  rtt_c4::broadcast(&first_groups_[0], fgsz, 0);
  if (rank != bcast_rank)
    indexes_.resize(isz);
  rtt_c4::broadcast(&indexes_[0], isz, 0);

  // Broadcast data itself
  if (rank != bcast_rank)
    data_.resize(dsz);
  rtt_c4::broadcast(&data_[0], dsz, 0);
  if (rank != bcast_rank)
    derivs_.resize(dsz);
  rtt_c4::broadcast(&derivs_[0], dsz, 0);
}

int Compton2::read_binary(std::string filename) {
  // Read
  auto fin = std::ifstream(filename, std::ios::in | std::ios::binary);

  // Ensure valid type
  char expected[] = " csk ";
  char file_type[sizeof(expected)];
  fin.read(file_type, sizeof(file_type));
  if (std::strcmp(file_type, expected) != 0) {
    std::cerr << "Expecting binary file " << filename << " to start with '"
              << expected << "' but got '" << file_type << "'";
    std::cerr << std::endl;
    return 1;
  }

  UINT binary_ordering, version_major, version_minor;
  fin.read(reinterpret_cast<char *>(&version_major), sizeof(UINT));
  fin.read(reinterpret_cast<char *>(&version_minor), sizeof(UINT));
  fin.read(reinterpret_cast<char *>(&binary_ordering), sizeof(UINT));
  if (version_major != 1 || binary_ordering > 1) {
    std::cerr << "Expecting a CSK binary file (version 1) with ordering 0 or 1 "
                 "but got "
              << version_major << " with ordering " << binary_ordering;
    std::cerr << std::endl;
    return 2;
  }

  constexpr UINT n = 7;
  std::array<UINT, n> szs;
  for (UINT i = 0; i < n; ++i)
    fin.read(reinterpret_cast<char *>(&szs[i]), sizeof(szs[i]));
  UINT j = 0;
  UINT tsz = szs[j++];
  UINT gsz = szs[j++];
  UINT lsz = szs[j++];
  UINT esz = szs[j++]; // number of evals (points, really)
  UINT fgsz = szs[j++];
  UINT isz = szs[j++];
  UINT dsz = szs[j++];

  num_temperatures_ = tsz;
  num_groups_ = gsz;
  num_leg_moments_ = lsz;
  num_evals_ = esz;
  UINT egsz = gsz + 1;

  std::cout << "DBG num_temperatures_ " << num_temperatures_ << '\n';
  std::cout << "DBG num_groups_ " << num_groups_ << '\n';
  std::cout << "DBG num_leg_moments_ " << num_leg_moments_ << '\n';
  std::cout << "DBG num_evals_ " << num_evals_ << '\n';
  std::cout << "DBG len(first_groups_) " << fgsz << '\n';
  std::cout << "DBG len(indexes_) " << isz << '\n';
  std::cout << "DBG len(data_/derivs_) " << dsz << '\n';

  Ts_.resize(tsz);
  for (UINT i = 0; i < tsz; ++i)
    fin.read(reinterpret_cast<char *>(&Ts_[i]), sizeof(Ts_[i]));

  Egs_.resize(egsz);
  for (UINT i = 0; i < egsz; ++i)
    fin.read(reinterpret_cast<char *>(&Egs_[i]), sizeof(Egs_[i]));

  first_groups_.resize(fgsz);
  for (UINT i = 0; i < fgsz; ++i)
    fin.read(reinterpret_cast<char *>(&first_groups_[i]),
             sizeof(first_groups_[i]));

  indexes_.resize(isz);
  for (UINT i = 0; i < isz; ++i)
    fin.read(reinterpret_cast<char *>(&indexes_[i]), sizeof(indexes_[i]));

  data_.resize(dsz);
  for (UINT i = 0; i < dsz; ++i)
    fin.read(reinterpret_cast<char *>(&data_[i]), sizeof(data_[i]));

  derivs_.resize(dsz);
  for (UINT i = 0; i < dsz; ++i)
    fin.read(reinterpret_cast<char *>(&derivs_[i]), sizeof(derivs_[i]));

  fin.close();
  return 0;
}

void Compton2::interp_matvec(vec &x, const vec &leftscale,
                             const vec &rightscale, double Te_keV,
                             bool zeroth_moment_only) const {
  // TODO: Redo interface? Not in-place??

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

void Compton2::interp_matvec_transpose(vec &xT, const vec &leftscale,
                                       const vec &rightscale, double Te_keV,
                                       bool zeroth_moment_only) const {
  // TODO: Redo interface? Not in-place??

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

void Compton2::interp_sparse_inscat(Sparse_Compton_Matrix &inscat,
                                    const vec &leftscale, const vec &rightscale,
                                    double Te_keV,
                                    bool zeroth_moment_only) const {
  // TODO: Allow offset to fill in data?

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;
  // TODO: implement
}

void Compton2::interp_dense_inscat(vec &inscat, const vec &leftscale,
                                   const vec &rightscale, double Te_keV,
                                   bool zeroth_moment_only) const {
  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  UINT const end_leg = zeroth_moment_only ? 1U : num_leg_moments_;

  inscat.resize(num_groups_ * num_groups_ * end_leg);
  std::fill(inscat.begin(), inscat.end(), 0.0);

  // TODO: implement
}

void Compton2::interp_linear_outscat(vec &outscat, const vec &leftscale,
                                     const vec &rightscale,
                                     double Te_keV) const {
  outscat.resize(num_groups_);
  std::fill(outscat.begin(), outscat.end(), 0.0);

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  // TODO: implement
}

void Compton2::interp_nonlinear_diff(vec &nldiff, const vec &leftscale,
                                     const vec &rightscale, const vec &flux,
                                     double flux_scale, double Te_keV) const {
  // Need "const" that changes based on DBC level?
  // Require(leftscale.len() == num_groups_);
  // Require(rightscale.len() == num_groups_);
  // Require(flux.len() == num_groups_);
  // Require(flux_scale > 0.0);
  // Require(Te_keV >= 0.0);

  nldiff.resize(num_groups_);
  std::fill(nldiff.begin(), nldiff.end(), 0.0);

  const vec &L = leftscale;
  const vec &R = rightscale;
  FP const fscale = flux_scale;
  FP const T = std::min(Ts_[num_temperatures_ - 1], std::max(Ts_[0], Te_keV));
  // TODO: implement

  // Ensure(nldiff.len() == num_groups_);
}

} // namespace rtt_compton2

//----------------------------------------------------------------------------//
// End compton2/Compton2.cc
//----------------------------------------------------------------------------//
