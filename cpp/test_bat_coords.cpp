#include "bat/bat_coords.hpp"
#include "structure/bfs_tree.hpp"
#include "structure/protein_root.hpp"
#include "structure/structure.hpp"
#include "structure/structure_graph.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static constexpr double PI = 3.141592653589793238462643383;

static void assert_close(double got, double expected, double eps = 1e-10,
                          const char* label = "") {
    if (std::abs(got - expected) > eps) {
        std::cerr << "FAIL [" << label << "]: expected " << expected
                  << "  got " << got
                  << "  diff " << (got - expected) << "\n";
        std::abort();
    }
}

// Build a minimal Atom with the fields select_first_backbone_frame needs.
static geometry::Atom make_atom(int serial, const std::string& name,
                                const std::string& elem, int resid, char chain,
                                double x, double y, double z) {
    geometry::Atom a;
    a.index        = serial;
    a.name         = name;
    a.element      = elem;
    a.residue_name = "ALA";
    a.residue_id   = resid;
    a.chain_id     = chain;
    a.alt_loc      = ' ';
    a.occupancy    = 1.0;
    a.position     = geometry::Vec3{x, y, z};
    return a;
}

// ---------------------------------------------------------------------------
// Test: basic BAT computation with virtual-reference substitution
//
// Structure (6 atoms, chain A, residue 1):
//
//   idx  name  position
//    0    N    (0, 0, 0)   <-- anchor N
//    1    CA   (1, 0, 0)   <-- anchor CA (BFS root)
//    2    C    (1, 1, 0)   <-- anchor C
//    3    CB   (1, 0, 1)   <-- depth-1 child of CA
//    4    O    (2, 1, 0)   <-- depth-1 child of C (which IS anchor C)
//    5    CG   (2, 0, 1)   <-- depth-2 child of CB
//
// Bonds: N-CA, CA-C, CA-CB, C-O, CB-CG
//
// BFS from CA (sorted neighbors):
//   order: CA(1), N(0), C(2), CB(3), O(4), CG(5)
//
// refs:
//   [0] N  = {CA=1, npos, npos}
//   [1] CA = {npos, npos, npos}  root
//   [2] C  = {CA=1, npos, npos}
//   [3] CB = {CA=1, npos, npos}
//   [4] O  = {C=2,  CA=1, npos}
//   [5] CG = {CB=3, CA=1, npos}
//
// Entries (anchor atoms 0,1,2 skipped): CB, O, CG
//
// Expected values (hand-verified):
//
//   CB (3): p=CA(1), g→N(0), h→C(2)
//     bond  = dist((1,0,1),(1,0,0))         = 1.0
//     angle = angle(CB,CA,N) at CA          = pi/2
//     torsion = dihedral(CB,CA,N,C)         = +pi/2
//       [n1=(0,1,0) n2=(0,0,-1) b2hat=(-1,0,0)
//        cross(n1,n2)=(-1,0,0)  y=dot((-1,0,0),(-1,0,0))=+1  x=0]
//
//   O  (4): p=C(2), g=CA(1), h→C(2)? C==p → h=N(0)
//     bond  = dist((2,1,0),(1,1,0))         = 1.0
//     angle = angle(O,C,CA) at C            = pi/2
//     torsion = dihedral(O,C,CA,N)          = pi
//       [n1=(0,0,1) n2=(0,0,-1) b2hat=(0,-1,0)
//        cross(n1,n2)=(0,0,0)   y=0  x=dot(n1,n2)=-1]
//
//   CG (5): p=CB(3), g=CA(1), h→C(2) [C≠p=CB, C≠g=CA]
//     bond  = dist((2,0,1),(1,0,1))         = 1.0
//     angle = angle(CG,CB,CA) at CB         = pi/2
//     torsion = dihedral(CG,CB,CA,C)        = -pi/2
//       [n1=(0,-1,0) n2=(1,0,0) b2hat=(0,0,-1)
//        cross(n1,n2)=(0,0,1)   y=dot((0,0,1),(0,0,-1))=-1  x=0]
// ---------------------------------------------------------------------------
static void test_basic_bat_computation() {
    geometry::Structure st;
    st.atoms = {
        make_atom(1, "N",  "N", 1, 'A', 0, 0, 0),  // idx 0
        make_atom(2, "CA", "C", 1, 'A', 1, 0, 0),  // idx 1
        make_atom(3, "C",  "C", 1, 'A', 1, 1, 0),  // idx 2
        make_atom(4, "CB", "C", 1, 'A', 1, 0, 1),  // idx 3
        make_atom(5, "O",  "O", 1, 'A', 2, 1, 0),  // idx 4
        make_atom(6, "CG", "C", 1, 'A', 2, 0, 1),  // idx 5
    };
    st.bonds = {
        {0, 1},  // N - CA
        {1, 2},  // CA - C
        {1, 3},  // CA - CB
        {2, 4},  // C  - O
        {3, 5},  // CB - CG
    };

    auto adj   = geometry::build_adjacency(st);
    auto frame = geometry::select_first_backbone_frame(st, 'A');
    auto tree  = geometry::build_bfs_forest(adj, frame.CA, /*sort=*/true);
    auto bat   = geometry::compute_bat(st, tree, frame);

    // --- Anchor ---
    assert(bat.anchor.N_idx  == 0);
    assert(bat.anchor.CA_idx == 1);
    assert(bat.anchor.C_idx  == 2);
    assert_close(bat.anchor.N_pos.x,  0.0, 1e-15, "N_pos.x");
    assert_close(bat.anchor.CA_pos.x, 1.0, 1e-15, "CA_pos.x");
    assert_close(bat.anchor.C_pos.y,  1.0, 1e-15, "C_pos.y");

    // --- Number of entries (6 atoms - 3 anchors = 3) ---
    assert(bat.entries.size() == 3);

    // --- Entry 0: CB (atom_idx 3) ---
    {
        const auto& e = bat.entries[0];
        assert(e.atom_idx == 3);
        assert_close(e.bond,    1.0,     1e-10, "CB bond");
        assert_close(e.angle,   PI/2.0,  1e-10, "CB angle");
        assert_close(e.torsion, PI/2.0,  1e-10, "CB torsion");
    }

    // --- Entry 1: O (atom_idx 4) ---
    // Virtual h: anchor C == p(=C), so h falls back to anchor N.
    {
        const auto& e = bat.entries[1];
        assert(e.atom_idx == 4);
        assert_close(e.bond,    1.0,    1e-10, "O bond");
        assert_close(e.angle,   PI/2.0, 1e-10, "O angle");
        assert_close(e.torsion, PI,     1e-10, "O torsion");
    }

    // --- Entry 2: CG (atom_idx 5) ---
    // Virtual h: anchor C (idx 2), which is not p(=CB=3) nor g(=CA=1). OK.
    {
        const auto& e = bat.entries[2];
        assert(e.atom_idx == 5);
        assert_close(e.bond,    1.0,     1e-10, "CG bond");
        assert_close(e.angle,   PI/2.0,  1e-10, "CG angle");
        assert_close(e.torsion, -PI/2.0, 1e-10, "CG torsion");
    }
}

// ---------------------------------------------------------------------------
// Test: throws on invalid frame
// ---------------------------------------------------------------------------
static void test_throws_on_invalid_frame() {
    geometry::Structure st;
    st.atoms = { make_atom(1, "N", "N", 1, 'A', 0, 0, 0) };
    st.bonds = {};

    auto adj  = geometry::build_adjacency(st);
    auto tree = geometry::build_bfs_forest(adj, 0, true);

    geometry::ProteinFrame bad_frame;  // all fields remain npos -> invalid()
    bool threw = false;
    try {
        (void)geometry::compute_bat(st, tree, bad_frame);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw && "Expected runtime_error for invalid frame");
}

// ---------------------------------------------------------------------------
// Test: throws on size mismatch between tree and structure
// ---------------------------------------------------------------------------
static void test_throws_on_size_mismatch() {
    geometry::Structure st;
    st.atoms = { make_atom(1, "N", "N", 1, 'A', 0, 0, 0),
                 make_atom(2, "CA","C", 1, 'A', 1, 0, 0),
                 make_atom(3, "C", "C", 1, 'A', 1, 1, 0) };
    st.bonds = { {0,1}, {1,2} };

    // Build tree for a 2-atom structure, pass 3-atom structure
    geometry::Structure small;
    small.atoms = { st.atoms[0], st.atoms[1] };
    small.bonds = { {0, 1} };
    auto adj_small  = geometry::build_adjacency(small);

    auto frame = geometry::select_first_backbone_frame(st, 'A');
    // tree built from 2-node adjacency; structure has 3 atoms -> mismatch
    auto tree_small = geometry::build_bfs_forest(adj_small, 1, true);

    bool threw = false;
    try {
        (void)geometry::compute_bat(st, tree_small, frame);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw && "Expected runtime_error for size mismatch");
}

// ---------------------------------------------------------------------------
int main() {
    test_basic_bat_computation();
    test_throws_on_invalid_frame();
    test_throws_on_size_mismatch();

    std::cout << "PASS: bat_coords computation tests\n";
    return 0;
}
