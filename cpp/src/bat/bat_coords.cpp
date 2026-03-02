#include "bat/bat_coords.hpp"
#include "cartesian/geometry.hpp"

#include <stdexcept>
#include <string>

namespace geometry {

BATCoords compute_bat(
    const Structure&    st,
    const BFSTree&      tree,
    const ProteinFrame& frame
) {
    if (!frame.valid()) {
        throw std::runtime_error("compute_bat: invalid protein frame");
    }
    if (tree.n != st.atoms.size()) {
        throw std::runtime_error(
            "compute_bat: BFS tree size (" + std::to_string(tree.n) +
            ") does not match structure atom count (" +
            std::to_string(st.atoms.size()) + ")"
        );
    }

    BATCoords bat;

    // --- Fill anchor ---
    bat.anchor.N_idx  = frame.N;
    bat.anchor.CA_idx = frame.CA;
    bat.anchor.C_idx  = frame.C;
    bat.anchor.N_pos  = st.atoms[frame.N].position;
    bat.anchor.CA_pos = st.atoms[frame.CA].position;
    bat.anchor.C_pos  = st.atoms[frame.C].position;

    const std::size_t aN  = frame.N;
    const std::size_t aCA = frame.CA;
    const std::size_t aC  = frame.C;

    bat.entries.reserve(tree.n >= 3 ? tree.n - 3 : 0);

    // --- Walk atoms in BFS order ---
    for (const std::size_t v : tree.order) {

        // Skip the three Cartesian anchor atoms.
        if (v == aN || v == aCA || v == aC) continue;

        const BFSTree::Ref3& r = tree.refs[v];
        const std::size_t p = r.p;

        if (p == BFSTree::npos) {
            // v is a BFS root that is not the CA anchor —
            // a disconnected component that we cannot handle.
            throw std::runtime_error(
                "compute_bat: atom " + std::to_string(v) +
                " is a BFS component root but is not the CA anchor; "
                "the structure must be a single connected molecule"
            );
        }

        // Resolve effective grandparent g.
        // Missing when v is at depth 1 (parent is the root CA).
        std::size_t g = r.g;
        if (g == BFSTree::npos) {
            g = aN;  // virtual grandparent: anchor N
        }

        // Resolve effective great-grandparent h.
        // Missing when v is at depth 1 or 2.
        std::size_t h = r.h;
        if (h == BFSTree::npos) {
            // Prefer anchor C; fall back to anchor N if C would collide with p or g.
            h = (aC != p && aC != g) ? aC : aN;
        }

        const Vec3& pos_v = st.atoms[v].position;
        const Vec3& pos_p = st.atoms[p].position;
        const Vec3& pos_g = st.atoms[g].position;
        const Vec3& pos_h = st.atoms[h].position;

        BATCoords::Entry e;
        e.atom_idx = v;
        e.bond     = geometry::distance(pos_v, pos_p);
        e.angle    = geometry::angle(pos_v, pos_p, pos_g);
        e.torsion  = geometry::dihedral(pos_v, pos_p, pos_g, pos_h);

        bat.entries.push_back(e);
    }

    return bat;
}

} // namespace geometry
