#include "structure/structure_graph.hpp"
#include <stdexcept>

namespace geometry {

std::vector<std::vector<std::size_t>> build_adjacency(const Structure& st) {
    const std::size_t n = st.atoms.size();
    std::vector<std::vector<std::size_t>> adj(n);

    for (const auto& b : st.bonds) {
        if (b.i >= n || b.j >= n) {
            throw std::runtime_error("Bond index out of range in build_adjacency()");
        }
        adj[b.i].push_back(b.j);
        adj[b.j].push_back(b.i);
    }
    return adj;
}

} // namespace geometry