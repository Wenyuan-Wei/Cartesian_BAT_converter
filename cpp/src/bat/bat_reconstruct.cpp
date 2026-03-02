#include "bat/bat_reconstruct.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace geometry {

// ---------------------------------------------------------------------------
// Local Vec3 math helpers (keep this file self-contained; vec3_ops is private)
// ---------------------------------------------------------------------------
namespace {

inline double dot3(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return { a.y*b.z - a.z*b.y,
             a.z*b.x - a.x*b.z,
             a.x*b.y - a.y*b.x };
}

inline double norm3(const Vec3& v) {
    return std::sqrt(dot3(v, v));
}

inline Vec3 normalize3(const Vec3& v) {
    double n = norm3(v);
    return { v.x/n, v.y/n, v.z/n };
}

} // anonymous namespace

// ---------------------------------------------------------------------------
Vec3 place_atom(const Vec3& h, const Vec3& g, const Vec3& p,
                double bond, double angle_rad, double torsion_rad) {
    // Build a right-handed local frame at P:
    //   e_x : direction from G to P  (parent bond axis)
    //   e_y : component of G->H perpendicular to e_x  (torsion reference)
    //   e_z : cross(e_x, e_y)

    Vec3 gp { p.x - g.x, p.y - g.y, p.z - g.z };
    Vec3 e_x = normalize3(gp);

    Vec3 gh { h.x - g.x, h.y - g.y, h.z - g.z };
    double proj = dot3(gh, e_x);
    Vec3 v { gh.x - proj*e_x.x,
             gh.y - proj*e_x.y,
             gh.z - proj*e_x.z };
    double v_norm = norm3(v);

    Vec3 e_y;
    if (v_norm > 1e-9) {
        e_y = { v.x/v_norm, v.y/v_norm, v.z/v_norm };
    } else {
        // H is collinear with G-P; pick an arbitrary perpendicular.
        Vec3 arb = (std::abs(e_x.x) < 0.9) ? Vec3{1.0, 0.0, 0.0}
                                             : Vec3{0.0, 1.0, 0.0};
        double arb_proj = dot3(arb, e_x);
        Vec3 vv { arb.x - arb_proj*e_x.x,
                  arb.y - arb_proj*e_x.y,
                  arb.z - arb_proj*e_x.z };
        e_y = normalize3(vv);
    }

    Vec3 e_z = cross3(e_x, e_y);

    // NeRF displacement in local frame, then rotate to global.
    double dx = -bond * std::cos(angle_rad);
    double dy =  bond * std::sin(angle_rad) * std::cos(torsion_rad);
    double dz =  bond * std::sin(angle_rad) * std::sin(torsion_rad);

    return {
        p.x + e_x.x*dx + e_y.x*dy + e_z.x*dz,
        p.y + e_x.y*dx + e_y.y*dy + e_z.y*dz,
        p.z + e_x.z*dx + e_y.z*dy + e_z.z*dz
    };
}

// ---------------------------------------------------------------------------
std::vector<Vec3> reconstruct_cartesian(
    const BATCoords& bat,
    const BFSTree&   tree
) {
    if (bat.anchor.N_idx  >= tree.n ||
        bat.anchor.CA_idx >= tree.n ||
        bat.anchor.C_idx  >= tree.n) {
        throw std::runtime_error(
            "reconstruct_cartesian: anchor indices out of range for tree size " +
            std::to_string(tree.n));
    }

    const std::size_t aN  = bat.anchor.N_idx;
    const std::size_t aCA = bat.anchor.CA_idx;
    const std::size_t aC  = bat.anchor.C_idx;

    std::vector<Vec3> pos(tree.n);
    pos[aN]  = bat.anchor.N_pos;
    pos[aCA] = bat.anchor.CA_pos;
    pos[aC]  = bat.anchor.C_pos;

    // Entries are in BFS order, so every parent is placed before its children.
    for (const auto& e : bat.entries) {
        const std::size_t v = e.atom_idx;
        const BFSTree::Ref3& r = tree.refs[v];
        const std::size_t p = r.p;   // always valid (compute_bat already checked)

        // Same virtual-reference convention as compute_bat.
        std::size_t g = r.g;
        if (g == BFSTree::npos) g = aN;

        std::size_t h = r.h;
        if (h == BFSTree::npos)
            h = (aC != p && aC != g) ? aC : aN;

        pos[v] = place_atom(pos[h], pos[g], pos[p],
                            e.bond, e.angle, e.torsion);
    }

    return pos;
}

} // namespace geometry
