#pragma once

#include "bat/bat_coords.hpp"
#include "cartesian/vec3.hpp"
#include "structure/bfs_tree.hpp"

#include <cstddef>
#include <vector>

namespace geometry {

/// Place one atom using the NeRF (Natural Extension Reference Frame) formula.
///
/// Given three reference atoms h (great-grandparent), g (grandparent),
/// p (parent), and the BAT parameters, returns the Cartesian position of v
/// such that:
///   distance(v, p)       = bond
///   angle(v, p, g) at p  = angle   (radians, [0, pi])
///   dihedral(v, p, g, h) = torsion (radians, (-pi, pi])
///
/// Precondition: g != p (non-zero G->P bond vector).
/// If H is collinear with the G-P axis, an arbitrary perpendicular is used
/// for the torsion reference (torsion is ill-defined in that case anyway).
Vec3 place_atom(const Vec3& h, const Vec3& g, const Vec3& p,
                double bond, double angle, double torsion);

/// Reconstruct Cartesian positions from BAT coordinates.
///
/// Applies the identical virtual-reference convention as compute_bat, so the
/// round-trip  Structure -> compute_bat -> reconstruct_cartesian  is exact.
///
/// \param bat   BAT coordinates produced by compute_bat().
/// \param tree  The same BFS spanning forest that was used to produce bat.
/// \returns     Vector of positions, indexed identically to Structure::atoms
///              (i.e., positions[i] is the position of atom i).
///
/// \throws std::runtime_error if the anchor indices are out of range.
std::vector<Vec3> reconstruct_cartesian(
    const BATCoords& bat,
    const BFSTree&   tree
);

} // namespace geometry
