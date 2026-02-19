#pragma once
#include <string>
#include "structure.hpp"

namespace geometry {

struct PDBReader {
    static Structure read(const std::string& filename);
};

} // namespace geometry
