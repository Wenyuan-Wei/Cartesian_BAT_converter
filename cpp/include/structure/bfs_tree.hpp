#pragma once

#include <cstddef>
#include <vector>

namespace geometry {

struct BFSTree {
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

    // Graph size
    std::size_t n = 0;

    // Spanning forest structure
    std::vector<std::size_t> parent;     // parent[v] or npos if root of a component
    std::vector<int> depth;              // depth[v], -1 if unreachable (shouldn’t happen)
    std::vector<std::size_t> root_of;    // component root for each v (handy for debugging / grouping)
    std::vector<std::size_t> order;      // BFS discovery order over the whole forest

    // Forward traversal convenience
    std::vector<std::vector<std::size_t>> children; // children[u] in tree (not original graph)

    // BAT helper: for each node v, refs[v] = {p, g, h}
    // If not available (v is root/near-root), entries will be npos.
    struct Ref3 { std::size_t p, g, h; };
    std::vector<Ref3> refs;
};

/**
 * Build a deterministic BFS spanning forest.
 *
 * - adj: undirected adjacency list over nodes [0..n-1]
 * - first_root: BFS starts here, then continues with next unvisited node(s) as new roots
 * - sort_neighbors: if true, sort adjacency lists for determinism (recommended for BAT)
 */
BFSTree build_bfs_forest(
    const std::vector<std::vector<std::size_t>>& adj,
    std::size_t first_root,
    bool sort_neighbors = true
);

/**
 * Convenience wrapper: spanning tree only (assumes graph is connected).
 * If graph is disconnected, it still returns a forest; only the first component is rooted at root.
 */
inline BFSTree build_bfs_tree(
    const std::vector<std::vector<std::size_t>>& adj,
    std::size_t root,
    bool sort_neighbors = true
) {
    return build_bfs_forest(adj, root, sort_neighbors);
}

} // namespace geometry