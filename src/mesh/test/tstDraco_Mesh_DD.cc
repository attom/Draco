//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstDraco_Mesh_DD.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Sunday, Jun 24, 2018, 14:38 pm
 * \brief  Draco_Mesh class unit test.
 * \note   Copyright (C) 2018-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Test_Mesh_Interface.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"

using rtt_mesh::Draco_Mesh;
using rtt_mesh_test::Test_Mesh_Interface;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// 2D Cartesian domain-decomposed mesh construction test
void cartesian_mesh_2d_dd(rtt_c4::ParallelUnitTest &ut) {

  // use Cartesian geometry
  const Draco_Mesh::Geometry geometry = Draco_Mesh::Geometry::CARTESIAN;

  //>>> SET UP CELL AND NODE DATA

  // set the number of cells and nodes
  const size_t num_xdir = 1;
  const size_t num_ydir = 1;

  // generate a constainer for data needed in mesh construction
  std::shared_ptr<Test_Mesh_Interface> mesh_iface;
  if (rtt_c4::node() == 0) {
    mesh_iface.reset(new Test_Mesh_Interface(num_xdir, num_ydir, {0, 1, 3, 4}));
  } else {
    mesh_iface.reset(
        new Test_Mesh_Interface(num_xdir, num_ydir, {1, 2, 4, 5}, 1.0, 0.0));
  }

  // set ghost data
  std::vector<unsigned> ghost_cell_type = {2};
  std::vector<int> ghost_cell_number = {0};
  std::vector<unsigned> ghost_cell_to_node_linkage(2);
  std::vector<int> ghost_cell_rank(1);
  if (rtt_c4::node() == 0) {
    ghost_cell_to_node_linkage = {1, 3};
    ghost_cell_rank = {1};
  } else {
    ghost_cell_to_node_linkage = {2, 0};
    ghost_cell_rank = {0};
  }

  // short-cut to some arrays
  const std::vector<unsigned> &cell_type = mesh_iface->cell_type;
  const std::vector<unsigned> &cell_to_node_linkage =
      mesh_iface->cell_to_node_linkage;
  const std::vector<unsigned> &side_node_count = mesh_iface->side_node_count;
  const std::vector<unsigned> &side_to_node_linkage =
      mesh_iface->side_to_node_linkage;

  // instantiate the mesh
  std::shared_ptr<Draco_Mesh> mesh(new Draco_Mesh(
      mesh_iface->dim, geometry, cell_type, cell_to_node_linkage,
      mesh_iface->side_set_flag, side_node_count, side_to_node_linkage,
      mesh_iface->coordinates, mesh_iface->global_node_number,
      mesh_iface->face_type, ghost_cell_type, ghost_cell_to_node_linkage,
      ghost_cell_number, ghost_cell_rank));

  // check that the scalar data is correct
  FAIL_IF_NOT(mesh->get_dimension() == 2);
  FAIL_IF_NOT(mesh->get_geometry() == Draco_Mesh::Geometry::CARTESIAN);
  FAIL_IF_NOT(mesh->get_num_cells() == mesh_iface->num_cells);
  FAIL_IF_NOT(mesh->get_num_nodes() == mesh_iface->num_nodes);

  // check that the vector data is correct
  FAIL_IF_NOT(mesh->get_ghost_cell_numbers() == ghost_cell_number);
  FAIL_IF_NOT(mesh->get_ghost_cell_ranks() == ghost_cell_rank);

  // get the layout generated by the mesh
  const Draco_Mesh::Layout &layout = mesh->get_cc_linkage();

  // cell-to-cell linkage on each node (1 cell on each node, 0 on-node nhbrs)
  FAIL_IF_NOT(layout.size() == 0);

  // get the boundary layout generated by the mesh
  const Draco_Mesh::Layout &bd_layout = mesh->get_cs_linkage();

  // check that the boundary (or side) layout has been generated
  FAIL_IF_NOT(bd_layout.size() == mesh_iface->num_cells);

  // get the ghost layout generated by the mesh
  const Draco_Mesh::Layout &go_layout = mesh->get_cg_linkage();

  // check that the ghost cell layout has been generated
  FAIL_IF_NOT(go_layout.size() == mesh_iface->num_cells);

  // check that cell-to-node linkage data is correct
  {
    std::vector<unsigned> test_cn_linkage =
        mesh_iface->flatten_cn_linkage(layout, bd_layout, go_layout);

    // check that cn_linkage is a permutation of the original cell-node linkage
    std::vector<unsigned>::const_iterator cn_first =
        cell_to_node_linkage.begin();
    std::vector<unsigned>::const_iterator test_cn_first =
        test_cn_linkage.begin();
    for (unsigned cell = 0; cell < mesh_iface->num_cells; ++cell) {

      // nodes must only be permuted at the cell level
      FAIL_IF_NOT(std::is_permutation(test_cn_first,
                                      test_cn_first + cell_type[cell], cn_first,
                                      cn_first + cell_type[cell]));

      // update the iterators
      cn_first += cell_type[cell];
      test_cn_first += cell_type[cell];
    }
  }

  // check that ghost-cell-to-node linkage data is correct
  {
    std::vector<unsigned> test_gn_linkage =
        mesh_iface->flatten_sn_linkage(go_layout);

    // check that cn_linkage is a permutation of the original cell-node linkage
    std::vector<unsigned>::const_iterator gn_first =
        ghost_cell_to_node_linkage.begin();
    std::vector<unsigned>::const_iterator test_gn_first =
        test_gn_linkage.begin();
    size_t const num_ghost_cells = ghost_cell_type.size();
    for (unsigned ghost = 0; ghost < num_ghost_cells; ++ghost) {

      // check that sn_linkage is a permutation of original ghost-node linkage
      FAIL_IF_NOT(std::is_permutation(
          test_gn_first, test_gn_first + ghost_cell_type[ghost], gn_first,
          gn_first + ghost_cell_type[ghost]));

      // update the iterators
      gn_first += ghost_cell_type[ghost];
      test_gn_first += ghost_cell_type[ghost];
    }
  }

  // successful test output
  if (ut.numFails == 0)
    PASSMSG("2D domain-decomposed Draco_Mesh tests ok.");
  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    Insist(rtt_c4::nodes() == 2, "This test only uses 2 PE.");
    cartesian_mesh_2d_dd(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of mesh/test/tstDraco_Mesh_DD.cc
//---------------------------------------------------------------------------//
