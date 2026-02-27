#include "structure/bond_inference.hpp"
#include "chem/covalent_radii.hpp"
#include "chem/valence.hpp"
#include "structure/structure.hpp"
#include "vec3_ops.hpp"

#include <unordered_map>
#include <vector>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <string>
#include <sstream>
#include <algorithm>

namespace geometry {

static inline std::string atom_label(const Atom& a) {
    std::ostringstream oss;
    oss << a.residue_name << " " << a.chain_id << a.residue_id << " " << a.name
        << " (idx " << a.index << ", el " << a.element << ")";
    return oss.str();
}

struct CellKey {
    int x, y, z;
    bool operator==(const CellKey& o) const { return x==o.x && y==o.y && z==o.z; }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& k) const noexcept {
        // Simple hash combiner
        std::size_t h = 1469598103934665603ull;
        auto mix = [&](int v) {
            h ^= static_cast<std::size_t>(v) + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        };
        mix(k.x); mix(k.y); mix(k.z);
        return h;
    }
};

static inline CellKey cell_of(const Vec3& p, double cell_size) {
    return CellKey{
        static_cast<int>(std::floor(p.x / cell_size)),
        static_cast<int>(std::floor(p.y / cell_size)),
        static_cast<int>(std::floor(p.z / cell_size))
    };
}

std::vector<Bond> infer_bonds(const Structure& st, const BondInferOptions& opt) {
    const auto& atoms = st.atoms;
    const std::size_t n = atoms.size();
    if (n == 0) return {};

    // Cell size: use the bond-search radius (neighbor cutoff) for conservative neighbor search
    const double cell_size = opt.max_bond_search_distance;

    // Build grid: cell -> atom indices
    std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> grid;
    grid.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        const auto key = cell_of(atoms[i].position, cell_size);
        grid[key].push_back(i);
    }

    std::vector<Bond> bonds;
    bonds.reserve(n * 2); // rough guess for protein

    // Track degree for valence enforcement
    std::vector<int> degree(n, 0);

    auto maybe_add_bond = [&](std::size_t i, std::size_t j) {
        if (i == j) return;
        if (j < i) std::swap(i, j);

        const Atom& ai = atoms[i];
        const Atom& aj = atoms[j];

        if (!opt.include_hydrogen && (ai.element == "H" || aj.element == "H")) return;

        const double d = vec3_ops::distance(ai.position, aj.position);
        if (d < opt.min_distance) {
            throw std::runtime_error(
                "Atom overlap/invalid geometry: distance " + std::to_string(d) +
                " Å between " + atom_label(ai) + " and " + atom_label(aj)
            );
        }

        // Prefilter: skip pairs outside the bond-search radius
        if (d > opt.max_bond_search_distance) return;

        const double ri = chem::covalent_radius_angstrom(ai.element);
        const double rj = chem::covalent_radius_angstrom(aj.element);
        const double cutoff = ri + rj + opt.tolerance_angstrom;

        if (d <= cutoff) {
            // Tentatively accept; enforce valence caps immediately if enabled
            if (opt.enforce_valence) {
                const int max_i = chem::max_valence(ai.element);
                const int max_j = chem::max_valence(aj.element);
                if (degree[i] + 1 > max_i) {
                    throw std::runtime_error(
                        "Valence exceeded for " + atom_label(ai) +
                        ": degree would become " + std::to_string(degree[i] + 1) +
                        " > max " + std::to_string(max_i)
                    );
                }
                if (degree[j] + 1 > max_j) {
                    throw std::runtime_error(
                        "Valence exceeded for " + atom_label(aj) +
                        ": degree would become " + std::to_string(degree[j] + 1) +
                        " > max " + std::to_string(max_j)
                    );
                }
            }

            bonds.push_back(Bond{i, j});
            degree[i] += 1;
            degree[j] += 1;
        }
    };

    // Iterate atoms and neighbor cells
    for (std::size_t i = 0; i < n; ++i) {
        const auto ci = cell_of(atoms[i].position, cell_size);

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    const CellKey ck{ci.x + dx, ci.y + dy, ci.z + dz};
                    auto it = grid.find(ck);
                    if (it == grid.end()) continue;

                    for (std::size_t j : it->second) {
                        if (j <= i) continue; // avoid duplicates
                        maybe_add_bond(i, j);
                    }
                }
            }
        }
    }

    return bonds;
}

} // namespace geometry