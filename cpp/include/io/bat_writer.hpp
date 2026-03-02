#pragma once

#include "bat/bat_coords.hpp"
#include "structure/bfs_tree.hpp"

#include <string>

namespace geometry {

/// Write BAT coordinates to a text file.
///
/// The file stores the anchor Cartesian positions and, for every other atom,
/// its BFS reference indices (parent / grandparent / great-grandparent) and
/// (bond, angle, torsion) values.  All information needed by
/// reconstruct_cartesian() is present, so no original Structure or BFSTree is
/// required at read time.
///
/// \param filename  Path to write (created or overwritten).
/// \param bat       BAT coordinates produced by compute_bat().
/// \param tree      The same BFS spanning forest used to produce bat.
/// \throws std::runtime_error on I/O error.
void write_bat(const std::string& filename,
               const BATCoords&   bat,
               const BFSTree&     tree);

/// Holds everything returned by read_bat() — sufficient for reconstruct_cartesian().
struct BATFile {
    BATCoords bat;   ///< Anchor Cartesian positions + per-atom (bond, angle, torsion).
    BFSTree   tree;  ///< Only n and refs[] are populated (parent/depth/etc. are empty).
};

/// Read a BAT file produced by write_bat().
///
/// \param filename  Path to read.
/// \returns         BATFile with bat and tree populated.
/// \throws std::runtime_error on I/O error or malformed input.
BATFile read_bat(const std::string& filename);

} // namespace geometry
