#pragma once
#include <vector>
#include <cstddef>
#include "structure/structure.hpp"

namespace geometry {

std::vector<std::vector<std::size_t>> build_adjacency(const Structure& st);

// (optional next)
std::vector<int> connected_component_labels(const Structure& st);

} // namespace geometry