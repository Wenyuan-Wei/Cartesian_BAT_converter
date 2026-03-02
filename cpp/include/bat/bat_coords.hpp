#pragma once

#include "cartesian/vec3.hpp"
#include "structure/bfs_tree.hpp"
#include "structure/protein_root.hpp"
#include "structure/structure.hpp"

#include <cstddef>
#include <vector>

namespace geometry {

/// BAT (Bond-Angle-Torsion) coordinate representation of a protein structure.
///
/// The three backbone atoms N, CA, C of the first residue are stored as
/// Cartesian anchors; every other atom is represented by (bond, angle, torsion)
/// relative to its BFS parent/grandparent/great-grandparent.
///
/// Virtual-reference convention (applied when BFS refs are npos):
///   missing grandparent g       -> anchor.N_idx
///   missing great-grandparent h -> anchor.C_idx, unless anchor.C_idx == parent,
///                                  in which case anchor.N_idx
///
/// This convention is symmetric: compute_bat and reconstruct_cartesian both
/// apply it, guaranteeing exact round-trip recovery of all atom positions.
struct BATCoords {

    /// Cartesian anchor: N/CA/C stored as absolute positions.
    /// These define the global position and orientation of the structure.
    /// Indices refer to positions in the original Structure::atoms vector.
    struct Anchor {
        std::size_t N_idx  = BFSTree::npos;
        std::size_t CA_idx = BFSTree::npos;
        std::size_t C_idx  = BFSTree::npos;
        Vec3 N_pos{};
        Vec3 CA_pos{};
        Vec3 C_pos{};
    };
    Anchor anchor{};

    /// One entry per non-anchor atom, in BFS traversal order.
    struct Entry {
        std::size_t atom_idx;  ///< index into Structure::atoms
        double bond;           ///< distance to BFS parent (Angstrom)
        double angle;          ///< bond angle at parent (radians, in [0, pi])
        double torsion;        ///< dihedral angle (radians, in (-pi, pi])
    };
    std::vector<Entry> entries;
};

/// Compute BAT coordinates from a structure and its BFS spanning tree.
///
/// \param st     Molecular structure (atoms + bonds, already bond-inferred).
/// \param tree   BFS spanning forest built with first_root = frame.CA.
/// \param frame  Backbone frame providing the N/CA/C Cartesian anchor.
///
/// \throws std::runtime_error if frame is invalid, tree/structure sizes
///         mismatch, or a non-anchor disconnected-component root is found.
BATCoords compute_bat(
    const Structure&    st,
    const BFSTree&      tree,
    const ProteinFrame& frame
);

} // namespace geometry
