//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/X3D_Draco_Mesh_Reader.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>, Kendra Keady
 * \date   Wednesday, Jul 11, 2018, 14:24 pm
 * \brief  X3D_Draco_Mesh_Reader header file.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_mesh_X3D_Draco_Mesh_Reader_hh
#define rtt_mesh_X3D_Draco_Mesh_Reader_hh

#include "Draco_Mesh_Reader.hh"
#include "ds++/Assert.hh"
#include "ds++/config.h"
#include <map>
#include <string>
#include <vector>

namespace rtt_mesh {

//===========================================================================//
/*!
 * \class X3D_Draco_Mesh_Reader
 *
 * \brief Read x3d file format and supply data to Draco_Mesh_Builder.
 *
 * This class parses mesh data from an x3d file format and manipulates it into a
 * form that is compatible to the Draco_Mesh_Builder class.  The parsing of the
 * file stream to the map follows work by Kendra Keady on a parsing framework.
 *
 * For a description of the X3D mesh file layout, see:
 * https://xcp-confluence.lanl.gov/display/SIMC/Ingen->Flag+Data+Transfer
 *
 * \todo: Boundary condition data is evidently parsed separately in X3D, so
 * meshes generated from this reader will not have side flag data (which ids
 * boundary conditions.
 *
 * \todo: Consider using the Class_Parse_Table formalism developed by Kent Budge
 * as an alternative.
 */
//===========================================================================//

class X3D_Draco_Mesh_Reader : public Draco_Mesh_Reader {
public:
  // >>> TYPEDEFS
  typedef std::pair<std::string, std::vector<std::string>> Parsed_Element;
  typedef std::vector<Parsed_Element> Parsed_Elements;

private:
  // >>> DATA

  //! File name
  const std::string filename;

  //! Boundary file names (optional data)
  const std::vector<std::string> bdy_filenames;

  //! Vector of all parsed key-value data pairs (includes valueless delimiters)
  Parsed_Elements parsed_pairs;

  //! Header data map (header in x3d file)
  std::map<std::string, std::vector<int>> x3d_header_map;

  //! Node coordinates
  std::map<int, std::vector<double>> x3d_coord_map;

  //! Face-to-node map
  std::map<int, std::vector<int>> x3d_facenode_map;

  //! Cell-to-face map
  std::map<int, std::vector<int>> x3d_cellface_map;

  //! Side-to-node map (0-based indices, unlike other maps)
  std::map<int, std::vector<int>> x3d_sidenode_map;

public:
  //! Constructor
  DLL_PUBLIC_mesh
  X3D_Draco_Mesh_Reader(const std::string &filename_,
                        const std::vector<std::string> &bdy_filenames_ = {});

  // >>> SERVICES

  DLL_PUBLIC_mesh void read_mesh();

  // >>> ACCESSORS

  // header data
  unsigned get_process() const { return x3d_header_map.at("process")[0] - 1; }
  unsigned get_numdim() const { return x3d_header_map.at("numdim")[0]; }
  unsigned get_numcells() const { return x3d_header_map.at("elements")[0]; }
  unsigned get_numnodes() const { return x3d_header_map.at("nodes")[0]; }

  // coord data
  std::vector<double> get_nodecoord(size_t node) const {
    return x3d_coord_map.at(node + 1);
  }

  // accessors with deferred implementations
  unsigned get_celltype(size_t cell) const;
  std::vector<int> get_cellnodes(size_t cell) const;

  // data needed from x3d boundary file (?)
  // \todo: parse x3d boundary file
  unsigned get_numsides() const { return x3d_sidenode_map.size(); }
  unsigned get_sidetype(size_t side) const {
    return x3d_sidenode_map.at(side).size();
  }
  unsigned get_sideflag(size_t /*side*/) const { return 0; }
  std::vector<int> get_sidenodes(size_t side) const {
    return x3d_sidenode_map.at(side);
  }

private:
  // >>> SUPPORT FUNCTIONS

  Parsed_Elements::const_iterator find_iter_of_key(const Parsed_Elements &pairs,
                                                   std::string key,
                                                   size_t start = 0);

  template <typename KT, typename VT>
  std::map<KT, std::vector<VT>> map_x3d_block(const std::string &block_name,
                                              int &dist);

  template <typename KT> KT convert_key(const std::string &skey);

  std::vector<int> get_facenodes(int face) const;

  void read_bdy_files();
};

//---------------------------------------------------------------------------//
// EXPLICIT SPECIALIZATIONS
//---------------------------------------------------------------------------//
template <>
DLL_PUBLIC_mesh std::string
X3D_Draco_Mesh_Reader::convert_key<std::string>(const std::string &skey);

} // end namespace rtt_mesh

// implementation header file
#include "X3D_Draco_Mesh_Reader.i.hh"

#endif // rtt_mesh_X3D_Draco_Mesh_Reader_hh

//---------------------------------------------------------------------------//
// end of mesh/X3D_Draco_Mesh_Reader.hh
//---------------------------------------------------------------------------//
