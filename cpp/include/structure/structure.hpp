#pragma once
#include <string>
#include <vector>
#include <cstddef>
#include "cartesian/vec3.hpp"  // for geometry::Vec3

namespace geometry {

struct Atom {
    int index;
    std::string name;
    std::string residue_name;
    int residue_id;
    char chain_id;
    std::string element;
    Vec3 position;
    char alt_loc = ' ';
    double occupancy = 1.0;
};

struct Bond { std::size_t i, j; };

struct Structure {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;

    Structure select_chain(char chain) const;
};

} // namespace geometry
