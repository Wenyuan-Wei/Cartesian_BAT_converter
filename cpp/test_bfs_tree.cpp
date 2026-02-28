#include "structure/bfs_tree.hpp"

#include <cassert>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using geometry::BFSTree;

// Helper: count roots in forest
static std::size_t count_roots(const BFSTree& t) {
    std::size_t roots = 0;
    for (std::size_t v = 0; v < t.n; ++v) {
        if (t.parent[v] == BFSTree::npos) ++roots;
    }
    return roots;
}

// Helper: count tree edges (non-root nodes)
static std::size_t count_tree_edges(const BFSTree& t) {
    std::size_t e = 0;
    for (std::size_t v = 0; v < t.n; ++v) {
        if (t.parent[v] != BFSTree::npos) ++e;
    }
    return e;
}

// Helper: verify each child points back to parent and depths are consistent
static void verify_parent_child_consistency(const BFSTree& t) {
    for (std::size_t u = 0; u < t.n; ++u) {
        for (std::size_t v : t.children[u]) {
            assert(v < t.n);
            assert(t.parent[v] == u);
            assert(t.depth[v] == t.depth[u] + 1);
        }
    }
}

// Helper: verify no cycles in parent pointers (each node climbs to a root)
static void verify_no_cycles(const BFSTree& t) {
    for (std::size_t start = 0; start < t.n; ++start) {
        std::size_t slow = start;
        std::size_t fast = start;

        auto step = [&](std::size_t x) -> std::size_t {
            if (x == BFSTree::npos) return BFSTree::npos;
            return t.parent[x];
        };

        // Floyd cycle detection on parent links
        while (true) {
            slow = step(slow);
            fast = step(step(fast));
            if (slow == BFSTree::npos || fast == BFSTree::npos) break;
            assert(slow != fast && "Cycle detected in parent pointers!");
        }
    }
}

// Helper: verify BFS depth minimality (optional, but strong check)
// For each node v, depth[v] should equal shortest path distance from its component root.
static void verify_depths_are_shortest_paths(
    const std::vector<std::vector<std::size_t>>& adj,
    const BFSTree& t
) {
    const std::size_t n = adj.size();
    assert(t.n == n);

    // group nodes by component root
    for (std::size_t v = 0; v < n; ++v) {
        const std::size_t root = t.root_of[v];
        assert(root != BFSTree::npos);

        // BFS from root to compute true shortest distances (within full graph adjacency)
        std::vector<int> dist(n, -1);
        std::queue<std::size_t> q;
        dist[root] = 0;
        q.push(root);

        while (!q.empty()) {
            auto u = q.front(); q.pop();
            for (auto w : adj[u]) {
                if (dist[w] == -1) {
                    dist[w] = dist[u] + 1;
                    q.push(w);
                }
            }
        }

        assert(dist[v] == t.depth[v]);
    }
}

static void test_cycle_graph_connected() {
    // Graph: 0-1-2-3-0 (cycle) plus chord 1-3 (extra cycle)
    // BFS root = 0
    std::vector<std::vector<std::size_t>> adj(4);
    adj[0] = {1, 3};
    adj[1] = {0, 2, 3};
    adj[2] = {1, 3};
    adj[3] = {2, 0, 1};

    auto t = geometry::build_bfs_forest(adj, /*first_root=*/0, /*sort_neighbors=*/true);

    assert(t.n == 4);
    assert(count_roots(t) == 1);
    assert(count_tree_edges(t) == 3); // n-1 for connected component

    // Root properties
    assert(t.parent[0] == BFSTree::npos);
    assert(t.depth[0] == 0);

    verify_parent_child_consistency(t);
    verify_no_cycles(t);
    verify_depths_are_shortest_paths(adj, t);

    // refs sanity for a small graph: root has npos refs; children have p set
    assert(t.refs[0].p == BFSTree::npos);
    for (std::size_t v = 1; v < t.n; ++v) {
        assert(t.refs[v].p == t.parent[v]);
    }
}

static void test_disconnected_forest() {
    // Component A: 0-1-2
    // Component B: 3-4
    // Component C: 5 alone
    std::vector<std::vector<std::size_t>> adj(6);
    adj[0] = {1};
    adj[1] = {0, 2};
    adj[2] = {1};
    adj[3] = {4};
    adj[4] = {3};
    adj[5] = {};

    auto t = geometry::build_bfs_forest(adj, /*first_root=*/0, /*sort_neighbors=*/true);

    assert(t.n == 6);
    assert(count_roots(t) == 3);          // components: {0,1,2}, {3,4}, {5}
    assert(count_tree_edges(t) == 3);     // n - components = 6 - 3 = 3

    verify_parent_child_consistency(t);
    verify_no_cycles(t);
    verify_depths_are_shortest_paths(adj, t);

    // root_of grouping checks
    assert(t.root_of[0] == t.root_of[1] && t.root_of[1] == t.root_of[2]);
    assert(t.root_of[3] == t.root_of[4]);
    assert(t.root_of[5] == 5); // singleton root is itself in our builder
}

static void test_determinism_with_scrambled_neighbors() {
    // Same connected graph as first test, but adjacency lists are scrambled.
    std::vector<std::vector<std::size_t>> adj(4);
    adj[0] = {3, 1};
    adj[1] = {3, 2, 0};
    adj[2] = {3, 1};
    adj[3] = {1, 0, 2};

    auto t1 = geometry::build_bfs_forest(adj, 0, true);
    auto t2 = geometry::build_bfs_forest(adj, 0, true);

    // With sort_neighbors=true, parents and order should match exactly run-to-run.
    assert(t1.parent == t2.parent);
    assert(t1.order  == t2.order);
}

int main() {
    test_cycle_graph_connected();
    test_disconnected_forest();
    test_determinism_with_scrambled_neighbors();

    std::cout << "PASS: bfs_tree spanning forest tests\n";
    return 0;
}