//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Draco_Mesh.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:38 pm
 * \brief  Draco_Mesh class header file.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_mesh_Draco_Mesh_hh
#define rtt_mesh_Draco_Mesh_hh

#include "ds++/config.h"
#include "mesh_element/Geometry.hh"
#include <map>
#include <vector>

namespace rtt_mesh {

//===========================================================================//
/*!
 * \class Draco_Mesh
 *
 * \brief General unstructured mesh class.
 *
 * The Draco_Mesh class takes cell-node (or cell-vertex) data, and generates a
 * mesh with layout (cell adjacency) information.  This class also provides
 * basic services, including access to cell information.  This mesh is based on
 * an unstructured mesh implementation by Kent Budge.
 *
 * Two important features for a fully realized Draco_Mesh are the following:
 * 1) Layout, which stores cell connectivity and hence the mesh topology.
 *    a) It has an internal layout containing local cell-to-cell linkage,
 *    b) and a boundary layout with side and off-process linkage.
 * 2) Geometry, which implies a metric for distance between points.
 *
 * Possibly temporary features:
 * 1) The cell_type_ vector (argument to the constructor) is currently taken to
 *    be the number of nodes per cell.
 * 2) The layout data structure(s) will proabably be moved to a separate class,
 *    where accessors might be used on a flattened version.
 */
//===========================================================================//

class Draco_Mesh {
public:
  // >>> TYPEDEFS
  typedef rtt_mesh_element::Geometry Geometry;
  typedef std::map<unsigned,
                   std::vector<std::pair<unsigned, std::vector<unsigned>>>>
      Layout;
  typedef std::map<unsigned,
                   std::vector<std::pair<unsigned, std::vector<unsigned>>>>
      Boundary_Layout;

private:
  // >>> DATA

  // Dimension
  const unsigned dimension;

  // Geometry enumeration
  const Geometry geometry;

  // Number of cells
  const unsigned num_cells;

  // Number of nodes
  const unsigned num_nodes;

  // Side set flag (can be used for mapping BCs to sides)
  const std::vector<unsigned> side_set_flag;

  // vector subscripted with node index with coordinate vector
  const std::vector<std::vector<double>> node_coord_vec;

  // Layout of mesh: vector index is cell index, vector element is
  // description of cell's adjacency to other cells in the mesh.
  Layout cell_to_cell_linkage;
  // TODO: update this to a full layout class.

  // Boundary layout of mesh
  Boundary_Layout cell_to_side_linkage;

public:
  //! Constructor.
  DLL_PUBLIC_mesh Draco_Mesh(unsigned dimension_, Geometry geometry_,
                             const std::vector<unsigned> &cell_type_,
                             const std::vector<unsigned> &cell_to_node_linkage_,
                             const std::vector<unsigned> &side_set_flag_,
                             const std::vector<unsigned> &side_node_count_,
                             const std::vector<unsigned> &side_to_node_linkage_,
                             const std::vector<double> &coordinates_,
                             const std::vector<unsigned> &global_node_number_);

  // >>> ACCESSORS

  unsigned get_dimension() const { return dimension; }
  Geometry get_geometry() const { return geometry; }
  unsigned get_num_cells() const { return num_cells; }
  unsigned get_num_nodes() const { return num_nodes; }
  const Layout &get_cc_linkage() const { return cell_to_cell_linkage; }
  const Boundary_Layout &get_cs_linkage() const { return cell_to_side_linkage; }

  // >>> SERVICES

private:
  // >>> SUPPORT FUNCTIONS

  //! Calculate (merely de-serialize) the vector of node coordinates
  std::vector<std::vector<double>>
  compute_node_coord_vec(const std::vector<double> &coordinates) const;

  //! Calculate the cell-to-cell linkage (layout)
  // TODO: add layout class and complete temporary version of this function.
  void compute_cell_to_cell_linkage(
      const std::vector<unsigned> &cell_type,
      const std::vector<unsigned> &cell_to_node_linkage,
      const std::vector<unsigned> &side_node_count,
      const std::vector<unsigned> &side_to_node_linkage);
};

} // end namespace rtt_mesh

#endif // rtt_mesh_Draco_Mesh_hh

//---------------------------------------------------------------------------//
// end of mesh/Draco_Mesh.hh
//---------------------------------------------------------------------------//
