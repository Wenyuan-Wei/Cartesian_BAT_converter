#pragma once
#include <string_view>
#include <stdexcept>

namespace chem {

// Covalent radii (Å). Small table for protein elements.
// Values are “single-bond typical”; good enough for inference thresholds.
inline double covalent_radius_angstrom(std::string_view el) {
    if (el == "H")  return 0.31;
    if (el == "C")  return 0.76;
    if (el == "N")  return 0.71;
    if (el == "O")  return 0.66;
    if (el == "S")  return 1.05;
    throw std::runtime_error("Unsupported element for covalent radius: '" + std::string(el) + "'");
}

} // namespace chem