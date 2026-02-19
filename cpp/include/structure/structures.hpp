#pragma once
#include <string>
#include <vector>
#include "geometry/vec3.hpp"  // for geometry::Vec3

namespace geometry {

struct Atom {
    int index;
    std::string name;
    std::string residue_name;
    int residue_id;
    char chain_id;
    std::string element;
    Vec3 position;
};

struct Structure {
    std::vector<Atom> atoms;

    Structure select_chain(char chain) const;
};

} // namespace geometry
