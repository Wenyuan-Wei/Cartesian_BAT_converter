#pragma once
#include <vector>
#include <cstddef>
#include "structure/structure.hpp"

namespace geometry {

struct BondInferOptions {
    double tolerance_angstrom = 0.35;     // added to sum of covalent radii
    double min_distance = 0.60;           // below this => error (overlap)
    double max_bond_search_distance = 2.5;      // above this → skip entirely (neighbor cutoff)
    bool include_hydrogen = true;         // if false, treat H as non-bonding or ignore
    bool enforce_valence = true;          // hard error if exceeded
};

std::vector<Bond> infer_bonds(const Structure& st, const BondInferOptions& opt);

} // namespace geometry