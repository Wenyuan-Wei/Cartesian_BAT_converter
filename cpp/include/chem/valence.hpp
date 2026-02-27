#pragma once
#include <string_view>
#include <stdexcept>

namespace chem {

// Max valence for strict protein inference.
// These are “hard caps” to catch nonsense.
// You can tune N and S if needed.
inline int max_valence(std::string_view el) {
    if (el == "H")  return 1;
    if (el == "C")  return 4;
    if (el == "N")  return 4; // including NH4+/protonated cases
    if (el == "O")  return 2;
    if (el == "S")  return 2; // strict for cysteine/disulfide
    throw std::runtime_error("Unsupported element for valence: '" + std::string(el) + "'");
}

} // namespace chem