#pragma once

#include "structure/structure.hpp"   // for Structure, Atom
#include <array>
#include <cstddef>
#include <optional>

namespace geometry {

struct ProteinFrame {
    std::size_t N  = static_cast<std::size_t>(-1);
    std::size_t CA = static_cast<std::size_t>(-1);
    std::size_t C  = static_cast<std::size_t>(-1);

    bool valid() const {
        const auto npos = static_cast<std::size_t>(-1);
        return N != npos && CA != npos && C != npos;
    }
};

// If chain is std::nullopt, choose the first chain (lexicographically) that yields a valid frame.
ProteinFrame select_first_backbone_frame(
    const Structure& st,
    std::optional<char> chain = std::nullopt
);

// Convenience: choose BFS root (usually CA)
inline std::size_t default_bfs_root_from_frame(const ProteinFrame& f) {
    return f.CA;
}

} // namespace geometry