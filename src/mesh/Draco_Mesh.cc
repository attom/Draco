//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Draco_Mesh.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:38 pm
 * \brief  Draco_Mesh class implementation file.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Draco_Mesh.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <iostream>
#include <numeric>

namespace rtt_mesh {

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Draco_Mesh constructor.
 *
 * \param[in] dimension_ dimension of mesh
 * \param[in] geometry_ enumerator of possible coordinate system geometries
 * \param[in] cell_type_ number of vertices for each cell
 * \param[in] cell_to_node_linkage_ serialized map of cell indices to node
 * indices.
 * \param[in] side_set_flag_ map of side indices (per cell) to side flag (global
 * index for a side).
 * \param[in] side_node_count_ number of nodes per each cell on a side of
 * the mesh.
 * \param[in] side_to_node_linkage_ serialized map of side indices (per side
 * cell) to node indices.
 * \param[in] coordinates_ serialized map of node index to coordinate values.
 * \param[in] global_node_number_ map of local to global node index (vector
 * subscript is local node index and value is global node index; for one
 * process, this is the identity map).
 */
Draco_Mesh::Draco_Mesh(unsigned dimension_, Geometry geometry_,
                       const std::vector<unsigned> &cell_type_,
                       const std::vector<unsigned> &cell_to_node_linkage_,
                       const std::vector<unsigned> &side_set_flag_,
                       const std::vector<unsigned> &side_node_count_,
                       const std::vector<unsigned> &side_to_node_linkage_,
                       const std::vector<double> &coordinates_,
                       const std::vector<unsigned> &global_node_number_)
    : dimension(dimension_), geometry(geometry_), num_cells(cell_type_.size()),
      num_nodes(global_node_number_.size()), side_set_flag(side_set_flag_),
      node_coord_vec(compute_node_coord_vec(coordinates_)) {

  // Require(dimension_ <= 3);
  // TODO: generalize mesh generation to 1D,3D (and uncomment requirment above)
  Insist(dimension_ == 2, "dimension_ != 2");

  // require some constraints on vector sizes
  Require(cell_to_node_linkage_.size() ==
          std::accumulate(cell_type_.begin(), cell_type_.end(), 0u));
  Require(
      side_to_node_linkage_.size() ==
      std::accumulate(side_node_count_.begin(), side_node_count_.end(), 0u));
  Require(coordinates_.size() == dimension_ * global_node_number_.size());

  // build the layout (or linkage) of the mesh
  compute_cell_to_cell_linkage(cell_type_, cell_to_node_linkage_,
                               side_node_count_, side_to_node_linkage_);
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 *
 * \param[in] coordinates serialized map of node index to coordinate values
 * (passed from constructor).
 * \return a vector of vectors of size=dimension of coordinates.
 */
std::vector<std::vector<double>> Draco_Mesh::compute_node_coord_vec(
    const std::vector<double> &coordinates) const {

  Require(coordinates.size() == dimension * num_nodes);

  // resize this class's coordinate data member
  std::vector<std::vector<double>> ret_node_coord_vec(
      num_nodes, std::vector<double>(dimension));

  // de-serialize the vector of node coordinates
  std::vector<double>::const_iterator ncv_first = coordinates.begin();
  for (unsigned node = 0; node < num_nodes; ++node) {

    // use the cell_type to create a vector of node indices for this cell
    std::vector<double> coord_vec(ncv_first, ncv_first + dimension);

    // resize each entry to the number of dimensions
    ret_node_coord_vec[node] = coord_vec;

    // increment pointer
    ncv_first += dimension;
  }

  Ensure(ncv_first == coordinates.end());

  return ret_node_coord_vec;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 *
 * \param[in] cell_type number of vertices per cell.
 * \param[in] cell_to_node_linkage serial map of cell index to node indices.
 * \param[in] side_node_count number of verices per side.
 * \param[in] side_to_node_linkage serial map of side index to node indices.
 */
void Draco_Mesh::compute_cell_to_cell_linkage(
    const std::vector<unsigned> &cell_type,
    const std::vector<unsigned> &cell_to_node_linkage,
    const std::vector<unsigned> &side_node_count,
    const std::vector<unsigned> &side_to_node_linkage) {

  Require(cell_type.size() == num_cells);
  Require(cell_to_node_linkage.size() ==
          std::accumulate(cell_type.begin(), cell_type.end(), 0u));
  Require(side_to_node_linkage.size() ==
          std::accumulate(side_node_count.begin(), side_node_count.end(), 0u));

  // STEP 1: create de-serialized map of cell index to node indices

  std::map<unsigned, std::vector<unsigned>> cell_to_node_map;

  // initialize pointers into cell_to_node_linkage vector
  std::vector<unsigned>::const_iterator cn_first = cell_to_node_linkage.begin();
  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // use the cell_type to create a vector of node indices for this cell
    std::vector<unsigned> node_vec(cn_first, cn_first + cell_type[cell]);

    // increment the map with this cell and vector-of-nodes entry
    cell_to_node_map[cell] = node_vec;

    // update the vector pointer
    cn_first += cell_type[cell];
  }

  Check(cell_to_node_map.size() == num_cells);

  // STEP 2: create a node-to-cell map (inverse of step 1)

  std::map<unsigned, std::vector<unsigned>> node_to_cell_map;

  // push cell index to vector for each node
  for (unsigned cell = 0; cell < num_cells; ++cell) {
    for (auto node : cell_to_node_map[cell]) {
      node_to_cell_map[node].push_back(cell);
    }
  }

  // STEP 3: create a node-to-side map

  std::map<unsigned, std::vector<unsigned>> node_to_side_map;

  // get the number of sides
  const unsigned num_sides = side_node_count.size();
  Check(dimension == 2 ? side_to_node_linkage.size() == 2 * num_sides : true);

  // TODO: extend to 1D, 3D (merely loop over side_node_count)
  unsigned node_offset = 0;
  for (unsigned side = 0; side < num_sides; ++side) {

    // at the relevant node indices increment the side vectors
    node_to_side_map[side_to_node_linkage[node_offset]].push_back(side);
    node_to_side_map[side_to_node_linkage[node_offset + 1]].push_back(side);

    // increment offset
    node_offset += side_node_count[side];
  }

  Check(node_to_side_map.size() == num_sides);

  // STEP 4: identify faces and create cell to cell map

  // TODO: amend to include ghost faces
  // TODO: global face index?
  // TODO: extend to 1D, 3D

  // identify faces per cell and create cell-to-cell linkage
  // in 2D, faces will always have 2 nodes
  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // create a vector of all possible node pairs
    unsigned nper_cell = cell_to_node_map[cell].size();
    unsigned num_pairs = nper_cell * (nper_cell - 1) / 2;
    std::vector<std::vector<unsigned>> vec_node_vec(num_pairs,
                                                    std::vector<unsigned>(2));
    // TODO: change type of vec_node_vec?
    unsigned k = 0;
    for (unsigned i = 0; i < nper_cell; ++i) {
      for (unsigned j = i + 1; j < nper_cell; ++j) {

        // set kth pair entry to the size 2 vector of nodes
        vec_node_vec[k] = {cell_to_node_map[cell][i],
                           cell_to_node_map[cell][j]};

        // increment k
        k++;
      }
    }

    Check(k == num_pairs);

    // check if each pair constitutes a face
    for (unsigned l = 0; l < num_pairs; ++l) {

      // get adjacent cells from node-to-cell map
      // TODO: add DbC to ensure these are sorted from step 2
      const std::vector<unsigned> &vert0_cells =
          node_to_cell_map[vec_node_vec[l][0]];
      const std::vector<unsigned> &vert1_cells =
          node_to_cell_map[vec_node_vec[l][1]];

      // find common cells (set_intersection is low-complexity)
      // TODO: reserve size for cells_in_common
      std::vector<unsigned> cells_in_common;
      std::set_intersection(vert0_cells.begin(), vert0_cells.end(),
                            vert1_cells.begin(), vert1_cells.end(),
                            std::back_inserter(cells_in_common));

      // these nodes should have at least cell index "cell" in common
      Check(cells_in_common.size() >= 1);

      // populate cell-to-cell linkage
      // TODO: populate face-to-cell (and inverse) map here?
      if (cells_in_common.size() > 1) {
        for (auto oth_cell : cells_in_common) {
          if (oth_cell != cell)
            cell_to_cell_linkage[cell].push_back(
                std::make_pair(oth_cell, vec_node_vec[l]));
        }
      }

      // check if this cell pair has a side flag
      // TODO: add DbC to ensure these are sorted from step 3
      const std::vector<unsigned> &vert0_sides =
          node_to_side_map[vec_node_vec[l][0]];
      const std::vector<unsigned> &vert1_sides =
          node_to_side_map[vec_node_vec[l][1]];

      // find common sides
      std::vector<unsigned> sides_in_common;
      std::set_intersection(vert0_sides.begin(), vert0_sides.end(),
                            vert1_sides.begin(), vert1_sides.end(),
                            std::back_inserter(sides_in_common));

      Check(sides_in_common.size() <= 1);
      if (sides_in_common.size() > 0) {
        // populate cell-to-side linkage
        // TODO: replace with typedef Draco_Layout::Boundary_Layout
        cell_to_side_linkage[cell].push_back(
            std::make_pair(sides_in_common[0], vec_node_vec[l]));
      }

      // TODO: add check for ghost face
    }
  }

  // STEP 5: instantiate the full layout
  // TODO: finish Draco_Layout class

  Ensure(cell_to_cell_linkage.size() == num_cells);
  Ensure(cell_to_side_linkage.size() == num_cells);
}
} // end namespace rtt_mesh

//---------------------------------------------------------------------------//
// end of mesh/Draco_Mesh.cc
//---------------------------------------------------------------------------//
