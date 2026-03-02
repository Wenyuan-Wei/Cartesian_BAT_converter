#include "bat/bat_coords.hpp"
#include "bat/bat_reconstruct.hpp"
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

static void assert_close(double got, double expected, double eps,
                          const char* label) {
    if (std::abs(got - expected) > eps) {
        std::cerr << "FAIL [" << label << "]: expected " << expected
                  << "  got " << got
                  << "  diff " << (got - expected) << "\n";
        std::abort();
    }
}

static void assert_vec3_close(const geometry::Vec3& got,
                               const geometry::Vec3& expected,
                               double eps, const char* label) {
    double dx = got.x - expected.x;
    double dy = got.y - expected.y;
    double dz = got.z - expected.z;
    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (dist > eps) {
        std::cerr << "FAIL [" << label << "]: expected ("
                  << expected.x << ", " << expected.y << ", " << expected.z
                  << ")  got ("
                  << got.x << ", " << got.y << ", " << got.z
                  << ")  dist " << dist << "\n";
        std::abort();
    }
}

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
// Test: place_atom matches hand-verified positions
//
// Using the same 6-atom structure from test_bat_coords.cpp:
//   N=(0,0,0), CA=(1,0,0), C=(1,1,0), CB=(1,0,1), O=(2,1,0), CG=(2,0,1)
// After BFS from CA, the three non-anchor atoms and their refs:
//
//   CB: p=CA(1,0,0), g=N(0,0,0), h=C(1,1,0)   bond=1, angle=pi/2, torsion=+pi/2
//   O:  p=C(1,1,0),  g=CA(1,0,0), h=N(0,0,0)  bond=1, angle=pi/2, torsion=pi
//   CG: p=CB(1,0,1), g=CA(1,0,0), h=C(1,1,0)  bond=1, angle=pi/2, torsion=-pi/2
// ---------------------------------------------------------------------------
static void test_place_atom_known_positions() {
    static constexpr double PI  = 3.141592653589793238462643383;
    static constexpr double EPS = 1e-10;

    geometry::Vec3 N  {0, 0, 0};
    geometry::Vec3 CA {1, 0, 0};
    geometry::Vec3 C  {1, 1, 0};
    geometry::Vec3 CB {1, 0, 1};

    // CB: h=C, g=N, p=CA
    geometry::Vec3 cb_got = geometry::place_atom(C, N, CA, 1.0, PI/2.0, PI/2.0);
    assert_vec3_close(cb_got, CB, EPS, "place_atom CB");

    // O: h=N(virtual), g=CA, p=C   (anchor.C == p, so virtual h = anchor.N)
    geometry::Vec3 o_expected {2, 1, 0};
    geometry::Vec3 o_got = geometry::place_atom(N, CA, C, 1.0, PI/2.0, PI);
    assert_vec3_close(o_got, o_expected, EPS, "place_atom O");

    // CG: h=C(virtual anchor.C, not equal to p=CB or g=CA), g=CA, p=CB
    geometry::Vec3 cg_expected {2, 0, 1};
    geometry::Vec3 cg_got = geometry::place_atom(C, CA, CB, 1.0, PI/2.0, -PI/2.0);
    assert_vec3_close(cg_got, cg_expected, EPS, "place_atom CG");
}

// ---------------------------------------------------------------------------
// Test: full round-trip Cartesian -> BAT -> Cartesian
//
// Uses the same 6-atom structure. All reconstructed positions must match
// original positions to within 1e-9 Angstrom.
// ---------------------------------------------------------------------------
static void test_roundtrip_6atom() {
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

    auto pos = geometry::reconstruct_cartesian(bat, tree);

    assert(pos.size() == st.atoms.size());

    for (std::size_t i = 0; i < st.atoms.size(); ++i) {
        const geometry::Vec3& orig = st.atoms[i].position;
        std::string label = "atom[" + std::to_string(i) + "]";
        assert_vec3_close(pos[i], orig, 1e-9, label.c_str());
    }
}

// ---------------------------------------------------------------------------
// Test: throws when anchor indices are out of range
// ---------------------------------------------------------------------------
static void test_throws_on_bad_anchor() {
    // Build a valid 3-atom structure just so we have a tree of size 3.
    geometry::Structure st;
    st.atoms = {
        make_atom(1, "N",  "N", 1, 'A', 0, 0, 0),
        make_atom(2, "CA", "C", 1, 'A', 1, 0, 0),
        make_atom(3, "C",  "C", 1, 'A', 1, 1, 0),
    };
    st.bonds = { {0,1}, {1,2} };
    auto adj  = geometry::build_adjacency(st);
    auto tree = geometry::build_bfs_forest(adj, 1, true);

    // Manually craft a BATCoords with an out-of-range anchor index.
    geometry::BATCoords bad_bat;
    bad_bat.anchor.N_idx  = 99;   // out of range
    bad_bat.anchor.CA_idx = 1;
    bad_bat.anchor.C_idx  = 2;

    bool threw = false;
    try {
        (void)geometry::reconstruct_cartesian(bad_bat, tree);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw && "Expected runtime_error for out-of-range anchor");
}

// ---------------------------------------------------------------------------
int main() {
    test_place_atom_known_positions();
    test_roundtrip_6atom();
    test_throws_on_bad_anchor();

    std::cout << "PASS: bat_reconstruct round-trip tests\n";
    return 0;
}
