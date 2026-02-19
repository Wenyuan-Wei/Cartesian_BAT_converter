#pragma once
#include <string>
#include "structure/structures.hpp"

namespace geometry {

struct PDBReader {
    static Structure read(const std::string& filename, char chain = '\0');
};

} // namespace geometry
