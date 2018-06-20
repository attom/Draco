//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstDraco_Mesh.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:43 pm
 * \brief  Draco_Mesh class unit test.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "mesh/Draco_Mesh.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <bitset>

using rtt_mesh::Draco_Mesh;

//---------------------------------------------------------------------------//
// AUXILIARY FUNCTIONS AND CLASSES
//---------------------------------------------------------------------------//

std::vector<unsigned>
flatten_cn_linkage(const unsigned num_cells, const unsigned num_nodes_per_cell,
                   const Draco_Mesh::Layout &layout,
                   const Draco_Mesh::Boundary_Layout &bd_layout) {
  // this function assumes a fixed number of nodes per cell
  std::vector<unsigned> cn_linkage;
  cn_linkage.reserve(num_cells * num_nodes_per_cell);
  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // node vector
    std::vector<unsigned> node_vec;
    // insert node indices from interior cell-to-cell linkage
    for (auto lpair : layout.at(cell)) {
      node_vec.insert(node_vec.end(), lpair.second.begin(), lpair.second.end());
    }

    // insert node indices from boundary cell-to-side linkage
    for (auto lpair : bd_layout.at(cell)) {
      node_vec.insert(node_vec.end(), lpair.second.begin(), lpair.second.end());
    }

    // sort and unique
    std::sort(node_vec.begin(), node_vec.end());
    auto last = std::unique(node_vec.begin(), node_vec.end());
    node_vec.erase(last, node_vec.end());

    // insert unique nodes into linkage array
    cn_linkage.insert(cn_linkage.end(), node_vec.begin(), node_vec.end());
  }

  // return the flattened cell-to-node linkage
  return cn_linkage;
}

std::vector<unsigned>
flatten_sn_linkage(const unsigned num_cells, const unsigned num_sides,
                   const Draco_Mesh::Boundary_Layout &bd_layout) {
  // this function assumes 2 nodes per side
  std::vector<unsigned> sn_linkage(num_sides * 2);

  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // insert node indices from boundary cell-to-side linkage
    for (auto lpair : bd_layout.at(cell)) {

      // get the side
      unsigned side = lpair.first;

      // set the vertices
      sn_linkage[2 * side] = lpair.second[0];
      sn_linkage[2 * side + 1] = lpair.second[1];
    }
  }

  // return the flattened side-to-node linage
  return sn_linkage;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// 2D Cartesian mesh construction test
void cartesian_mesh_2d(rtt_c4::ParallelUnitTest &ut) {

  // instantiate a mesh, to be defined below
  std::shared_ptr<Draco_Mesh> mesh;

  //>>> SET UP CELL AND NODE DATA

  // set the number of cells and nodes
  const size_t num_xdir = 2;
  const size_t num_ydir = 1;
  const size_t num_cells = num_xdir * num_ydir;
  const size_t num_nodes = (num_xdir + 1) * (num_ydir + 1);

  // set the number of sides and side flags
  const size_t poff = num_xdir + num_ydir; // parallel side offset
  const size_t num_sides = 2 * poff;

  // TODO: parse mesh data instead of hard-coding
  // TODO: make mesh data setup a builder and move to separate file.

  // use two dimensions and Cartesian geometry
  const unsigned dim = 2;
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::CARTESIAN;

  // set cell type (quadrilateral cells)
  const unsigned num_nodes_per_cell = 4;
  std::vector<unsigned> cell_type(num_cells, num_nodes_per_cell);

  // set the cell-to-node linkage counterclockwise about each cell
  std::vector<unsigned> cell_to_node_linkage(num_cells * num_nodes_per_cell);
  for (size_t j = 0; j < num_ydir; ++j) {
    for (size_t i = 0; i < num_xdir; ++i) {

      // get a serial cell index
      size_t cell = i + num_xdir * j;

      // set each node entry per cell
      cell_to_node_linkage[4 * cell] = cell + j;
      cell_to_node_linkage[4 * cell + 1] = cell + j + 1;
      cell_to_node_linkage[4 * cell + 2] = cell + num_xdir + 1 + j + 1;
      cell_to_node_linkage[4 * cell + 3] = cell + num_xdir + 1 + j;
    }
  }

  // set two nodes per each side cell (only possibility in 2D)
  std::vector<unsigned> side_node_count(num_sides, 2);

  // calculate arrays storing data about the sides of the mesh
  std::vector<unsigned> side_set_flag(num_sides);
  std::vector<unsigned> side_to_node_linkage(2 * num_sides);
  // ... over top and bottom faces
  for (size_t i = 0; i < num_xdir; ++i) {
    // bottom face
    side_set_flag[i] = 1;
    side_to_node_linkage[2 * i] = i;
    side_to_node_linkage[2 * i + 1] = i + 1;
    // top face
    side_set_flag[i + poff] = 3;
    side_to_node_linkage[2 * (i + poff)] = num_nodes - 1 - i;
    side_to_node_linkage[2 * (i + poff) + 1] = num_nodes - 1 - i - 1;
  }
  // ... over left and right faces
  for (size_t j = 0; j < num_ydir; ++j) {
    // right face
    side_set_flag[j + num_xdir] = 2;
    side_to_node_linkage[2 * (j + num_xdir)] = num_xdir * (j + 1) + j;
    side_to_node_linkage[2 * (j + num_xdir) + 1] = num_xdir * (j + 2) + j + 1;
    // left face
    side_set_flag[j + poff + num_xdir] = 4;
    side_to_node_linkage[2 * (j + poff + num_xdir)] =
        num_nodes - 1 - num_xdir * (j + 1) - j;
    side_to_node_linkage[2 * (j + poff + num_xdir) + 1] =
        num_nodes - 1 - num_xdir * (j + 2) - (j + 1);
  }

  // set some coordinates and global node indices
  std::vector<double> coordinates(dim * num_nodes);
  std::vector<unsigned> global_node_number(num_nodes);
  for (unsigned i = 0; i < num_nodes; ++i) {

    // TODO: generalize coordinate generation
    std::bitset<8> b2(i);
    coordinates[dim * i] = static_cast<double>(b2[0]);
    coordinates[1 + dim * i] = static_cast<double>(b2[1]);

    // TODO: generalize for multi-processing
    global_node_number[i] = i;
  }

  // build the mesh
  mesh.reset(new Draco_Mesh(
      dim, geometry, cell_type, cell_to_node_linkage, side_set_flag,
      side_node_count, side_to_node_linkage, coordinates, global_node_number));

  // check that the scalar data is correct
  if (mesh->get_dimension() != 2)
    ITFAILS;
  if (mesh->get_geometry() != Draco_Mesh::Geometry::CARTESIAN)
    ITFAILS;
  if (mesh->get_num_cells() != num_cells)
    ITFAILS;
  if (mesh->get_num_nodes() != num_nodes)
    ITFAILS;

  // get the layout generated by the mesh
  const Draco_Mesh::Layout &layout = mesh->get_cc_linkage();

  // check that the layout has been generated
  if (layout.size() != num_cells)
    ITFAILS;

  // check that each cell has the correct neighbors
  {
    std::map<unsigned, std::vector<unsigned>> test_cell_map;
    for (unsigned j = 0; j < num_ydir; ++j) {
      for (unsigned i = 0; i < num_xdir; ++i) {

        // calculate the cell index
        unsigned cell = i + j * num_xdir;

        // calculate neighbor cell indices
        if (i > 0)
          test_cell_map[cell].push_back(cell - 1);
        if (i < num_xdir - 1)
          test_cell_map[cell].push_back(cell + 1);
        if (j > 0)
          test_cell_map[cell].push_back(cell - num_xdir);
        if (j < num_ydir - 1)
          test_cell_map[cell].push_back(cell + num_xdir);
      }
    }

    for (unsigned cell = 0; cell < num_cells; ++cell) {

      // get number of faces per cell in layout
      const unsigned num_faces = layout.at(cell).size();

      // check that the number of faces per cell is correct
      if (num_faces != test_cell_map[cell].size())
        ITFAILS;

      // check that cell neighbors are correct
      for (unsigned face = 0; face < num_faces; ++face) {
        if (layout.at(cell)[face].first != test_cell_map[cell][face])
          ITFAILS;
      }
    }
  }

  // get the boundary layout generated by the mesh
  const Draco_Mesh::Boundary_Layout &bd_layout = mesh->get_cs_linkage();

  // check that the boundary (or side) layout has been generated
  if (bd_layout.size() != num_cells)
    ITFAILS;

  // check that cell-to-node linkage data is correct
  {
    std::vector<unsigned> test_cn_linkage =
        flatten_cn_linkage(num_cells, num_nodes_per_cell, layout, bd_layout);

    // check that cn_linkage is a permutation of the original cell-node linkage
    std::vector<unsigned>::const_iterator cn_first =
        cell_to_node_linkage.begin();
    std::vector<unsigned>::const_iterator test_cn_first =
        test_cn_linkage.begin();
    for (unsigned cell = 0; cell < num_cells; ++cell) {

      // nodes must only be permuted at the cell level
      if (!std::is_permutation(test_cn_first, test_cn_first + cell_type[cell],
                               cn_first, cn_first + cell_type[cell]))
        ITFAILS;

      // update the iterators
      cn_first += cell_type[cell];
      test_cn_first += cell_type[cell];
    }
  }

  // check that each cell has the correct sides
  {
    std::vector<unsigned> test_sn_linkage =
        flatten_sn_linkage(num_cells, num_sides, bd_layout);

    // check that cn_linkage is a permutation of the original cell-node linkage
    std::vector<unsigned>::const_iterator sn_first =
        side_to_node_linkage.begin();
    std::vector<unsigned>::const_iterator test_sn_first =
        test_sn_linkage.begin();
    for (unsigned side = 0; side < num_sides; ++side) {

      // check that sn_linkage is a permutation of the original side-node linkage
      if (!std::is_permutation(test_sn_first,
                               test_sn_first + side_node_count[side], sn_first,
                               sn_first + side_node_count[side]))
        ITFAILS;

      // update the iterators
      sn_first += side_node_count[side];
      test_sn_first += side_node_count[side];
    }
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("2D Draco_Mesh tests ok.");
  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Insist(rtt_c4::nodes() == 1, "This test only uses 1 PE.");
    cartesian_mesh_2d(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of mesh/test/tstDraco_Mesh.cc
//---------------------------------------------------------------------------//
