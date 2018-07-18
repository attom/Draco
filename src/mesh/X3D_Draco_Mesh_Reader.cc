//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/X3D_Draco_Mesh_Reader.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>, Kendra Keady
 * \date   Thursday, Jul 12, 2018, 08:46 am
 * \brief  X3D_Draco_Mesh_Reader class implementation file.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "X3D_Draco_Mesh_Reader.hh"
#include "ds++/DracoStrings.hh"
#include <fstream>
#include <iostream>

namespace rtt_mesh {

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief X3D_Draco_Mesh_Reader constructor.
 *
 * \param[in] filename_ name of file to be parsed
 */
X3D_Draco_Mesh_Reader::X3D_Draco_Mesh_Reader(const std::string filename_)
    : filename(filename_) {
  // check for valid file name
  Insist(filename_.size() > 0, "No file name supplied.");
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 */
void X3D_Draco_Mesh_Reader::read_mesh() {

  Require(parsed_pairs.size() == 0);

  // STEP 1: open the file stream

  std::ifstream x3dfile(filename.c_str());

  // file must exist and be readable
  Insist(x3dfile.is_open(), "Failed to find or open specified X3D mesh file.");

  // STEP 2: parse file token stream into an initial vector of string pairs

  std::vector<std::pair<std::string, std::string>> raw_pairs;

  while (!x3dfile.eof()) {
    // get a line from the file
    std::string data_line;
    std::getline(x3dfile, data_line);

    // trim outermost white space
    data_line = rtt_dsxx::trim(data_line);

    // ignore empty lines
    if (data_line.size() == 0)
      continue;

    // find last char index of key
    size_t key_offset = data_line.find(' ');

    // first string without whitespace will be key
    std::string key = data_line.substr(0, key_offset);

    // check if anything is after key, to set as value
    std::string value = "";
    if (key_offset != std::string::npos)
      value = data_line.substr(key_offset);

    // add pairing even if value string is size 0 (gets headers and footers)
    raw_pairs.push_back(std::pair<std::string, std::string>(key, value));
  }

  // STEP 3: close the file stream

  x3dfile.close();

  // STEP 4: split strings in value at whitespace

  for (auto rpair : raw_pairs) {

    // split the string value into a vector
    std::vector<std::string> tmp_vec = rtt_dsxx::tokenize(rpair.second);

    // add string-vector pair to intermediate map
    parsed_pairs.push_back(Parsed_Element(rpair.first, tmp_vec));
  }

  // STEP 5: derive mesh data maps from parsed_pairs

  // keep track of point parsed_pairs data
  int dist = 0;
  if (parsed_pairs[0].first != "header")
    dist++;

  // parse x3d header block and generate x3d_header_map
  Remember(int dist_old = dist);
  x3d_header_map = map_x3d_block<std::string, int>("header", dist);
  Check(dist > 0);
  Check(dist > dist_old);

  // parse x3d node coordinate block and generate x3d_coord_map
  Remember(dist_old = dist);
  x3d_coord_map = map_x3d_block<int, double>("nodes", dist);
  Check(dist > dist_old);

  // parse x3d face-to-node index map
  Remember(dist_old = dist);
  x3d_facenode_map = map_x3d_block<int, int>("faces", dist);
  Check(dist > dist_old);

  // parse x3d cell-to-face index map
  Remember(dist_old = dist);
  x3d_cellface_map = map_x3d_block<int, int>("cells", dist);
  Check(dist > dist_old);

  Ensure(parsed_pairs.size() > 0);
  Ensure(x3d_header_map.size() > 0);
  Ensure(x3d_coord_map.size() > 0);
  Ensure(x3d_facenode_map.size() > 0);
  Ensure(x3d_cellface_map.size() > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return number of nodes for a given cell.
 *
 * \param[in] cell index of cell
 *
 * \return number of nodes for cell
 */
unsigned X3D_Draco_Mesh_Reader::get_celltype(size_t cell) const {

  // get the list of cell nodes
  const std::vector<int> node_indexes = get_cellnodes(cell);

  // merely the size of the vector of unique nodes
  unsigned num_nodes_pc = node_indexes.size();

  Ensure(num_nodes_pc > 0);
  return num_nodes_pc;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the vector of node indices for a given cell.
 *
 * \param[in] cell index of cell
 *
 * \return vector of int node indices
 */
std::vector<int> X3D_Draco_Mesh_Reader::get_cellnodes(size_t cell) const {

  Require(cell < static_cast<size_t>(x3d_header_map.at("elements")[0]));

  // x3d file's node, face, and cell indexes start from 1
  const std::vector<int> &cell_data = x3d_cellface_map.at(cell + 1);
  const size_t num_faces = cell_data[0];

  // calculate number of nodes for this cell
  std::vector<int> node_indexes;

  for (size_t i = 1; i <= num_faces; ++i) {

    // get the face index, which will by key for face-to-node map
    int face = cell_data[i];

    // number of nodes is first value after face index in x3d file
    const std::vector<int> &face_data = x3d_facenode_map.at(face);
    const size_t num_nodes = face_data[0];

    // push each node instance onto node vector (subtract 1 to get 0-based node)
    for (size_t j = 1; j <= num_nodes; ++j)
      node_indexes.push_back(face_data[j] - 1);
  }

  // reduce return vector to unique node entries
  std::sort(node_indexes.begin(), node_indexes.end());
  node_indexes.erase(std::unique(node_indexes.begin(), node_indexes.end()),
                     node_indexes.end());

  Ensure(node_indexes.size() > 0);
  return node_indexes;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Find iterator of a vector of pairs at a key value.
 *
 * \param[in] pairs const reference to the vector of pairs
 * \param[in] key string key to search
 * \param[in] start index from which to start search for key
 *
 * \return iterator to pair with key
 */
X3D_Draco_Mesh_Reader::Parsed_Elements::const_iterator
X3D_Draco_Mesh_Reader::find_iter_of_key(const Parsed_Elements &pairs,
                                        std::string key, size_t start) {
  auto start_it = pairs.begin() + start;
  auto it =
      std::find_if(start_it, pairs.end(),
                   [&key](const Parsed_Element &p) { return p.first == key; });
  return it;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Simply return key as specialization to numerical conversion.
 *
 * \param[in] skey string key
 *
 * \return numerical key of type "KT"
 */
template <>
std::string
X3D_Draco_Mesh_Reader::convert_key<std::string>(const std::string &skey) {
  std::string ret_key = skey;
  return ret_key;
}

} // end namespace rtt_mesh

//---------------------------------------------------------------------------//
// end of mesh/X3D_Draco_Mesh_Reader.cc
//---------------------------------------------------------------------------//
