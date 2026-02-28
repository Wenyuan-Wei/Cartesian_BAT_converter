#include "structure/bfs_tree.hpp"

#include <algorithm>
#include <queue>
#include <stdexcept>

namespace geometry {

static void validate_adj(const std::vector<std::vector<std::size_t>>& adj) {
    const std::size_t n = adj.size();
    for (std::size_t u = 0; u < n; ++u) {
        for (std::size_t v : adj[u]) {
            if (v >= n) throw std::runtime_error("Adjacency contains out-of-range node index");
        }
    }
}

BFSTree build_bfs_forest(
    const std::vector<std::vector<std::size_t>>& adj_in,
    std::size_t first_root,
    bool sort_neighbors
) {
    validate_adj(adj_in);

    const std::size_t n = adj_in.size();
    if (n == 0) return BFSTree{};
    if (first_root >= n) throw std::runtime_error("first_root out of range in build_bfs_forest()");

    // Copy adjacency only if we want to sort/dedup for determinism.
    std::vector<std::vector<std::size_t>> adj = adj_in;
    if (sort_neighbors) {
        for (auto& nbrs : adj) {
            std::sort(nbrs.begin(), nbrs.end());
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        }
    }

    BFSTree t;
    t.n = n;
    t.parent.assign(n, BFSTree::npos);
    t.depth.assign(n, -1);
    t.root_of.assign(n, BFSTree::npos);
    t.children.assign(n, {});
    t.order.reserve(n);

    // We’ll fill refs after parents are known.
    t.refs.assign(n, BFSTree::Ref3{BFSTree::npos, BFSTree::npos, BFSTree::npos});

    auto bfs_from_root = [&](std::size_t root) {
        std::queue<std::size_t> q;
        t.parent[root] = BFSTree::npos;
        t.depth[root] = 0;
        t.root_of[root] = root;
        q.push(root);

        while (!q.empty()) {
            const std::size_t u = q.front();
            q.pop();
            t.order.push_back(u);

            for (std::size_t v : adj[u]) {
                if (t.depth[v] != -1) continue; // already discovered
                t.parent[v] = u;
                t.depth[v] = t.depth[u] + 1;
                t.root_of[v] = root;
                t.children[u].push_back(v);
                q.push(v);
            }
        }
    };

    // First component rooted at first_root
    bfs_from_root(first_root);

    // Remaining components (if any)
    for (std::size_t v = 0; v < n; ++v) {
        if (t.depth[v] == -1) bfs_from_root(v);
    }

    // Build BAT-style reference triples (p, g, h) = (parent, grandparent, great-grandparent)
    for (std::size_t v = 0; v < n; ++v) {
        const std::size_t p = t.parent[v];
        if (p == BFSTree::npos) continue;
        const std::size_t g = t.parent[p];
        if (g == BFSTree::npos) { t.refs[v] = {p, BFSTree::npos, BFSTree::npos}; continue; }
        const std::size_t h = t.parent[g];
        t.refs[v] = {p, g, h}; // h may be npos, that’s OK
    }

    return t;
}

} // namespace geometry