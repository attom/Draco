//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Builder_pt.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Tuesday, Jul 03, 2018, 11:52 am
 * \brief  Draco_Mesh_Builder class header file.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Draco_Mesh_Builder.t.hh"
#include "RTT_Format_Reader/RTT_Format_Reader.hh"

namespace rtt_mesh {

template class Draco_Mesh_Builder<rtt_RTT_Format_Reader::RTT_Format_Reader>;

} // end namespace rtt_mesh

//---------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Builder_pt.cc
//---------------------------------------------------------------------------//
