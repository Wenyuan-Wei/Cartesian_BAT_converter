#include "bat/bat_coords.hpp"
#include "bat/bat_reconstruct.hpp"
#include "io/bat_writer.hpp"
#include "structure/bfs_tree.hpp"
#include "structure/protein_root.hpp"
#include "structure/structure.hpp"
#include "structure/structure_graph.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Build the canonical 6-atom test structure:
//   N=(0,0,0), CA=(1,0,0), C=(1,1,0), CB=(1,0,1), O=(2,1,0), CG=(2,0,1)
// ---------------------------------------------------------------------------
static geometry::Structure make_6atom_structure() {
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
        {0, 1},  // N  - CA
        {1, 2},  // CA - C
        {1, 3},  // CA - CB
        {2, 4},  // C  - O
        {3, 5},  // CB - CG
    };
    return st;
}

// ---------------------------------------------------------------------------
// Test: write → read → reconstruct gives positions within 1e-9 Å of originals.
// ---------------------------------------------------------------------------
static void test_write_read_roundtrip() {
    auto st    = make_6atom_structure();
    auto adj   = geometry::build_adjacency(st);
    auto frame = geometry::select_first_backbone_frame(st, 'A');
    auto tree  = geometry::build_bfs_forest(adj, frame.CA, /*sort=*/true);
    auto bat   = geometry::compute_bat(st, tree, frame);

    const std::string tmp = "/tmp/test_bat_writer_out.bat";
    geometry::write_bat(tmp, bat, tree);

    geometry::BATFile loaded = geometry::read_bat(tmp);
    std::remove(tmp.c_str());

    auto pos = geometry::reconstruct_cartesian(loaded.bat, loaded.tree);

    assert(pos.size() == st.atoms.size());
    for (std::size_t i = 0; i < st.atoms.size(); ++i) {
        std::string label = "atom[" + std::to_string(i) + "]";
        assert_vec3_close(pos[i], st.atoms[i].position, 1e-9, label.c_str());
    }
}

// ---------------------------------------------------------------------------
// Test: reading a nonexistent file throws std::runtime_error.
// ---------------------------------------------------------------------------
static void test_throws_on_missing_file() {
    bool threw = false;
    try {
        geometry::read_bat("/tmp/this_file_does_not_exist_bat_writer_test.bat");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw && "Expected runtime_error for missing file");
}

// ---------------------------------------------------------------------------
// Test: writing to an unwritable path throws std::runtime_error.
// ---------------------------------------------------------------------------
static void test_throws_on_bad_write_path() {
    auto st    = make_6atom_structure();
    auto adj   = geometry::build_adjacency(st);
    auto frame = geometry::select_first_backbone_frame(st, 'A');
    auto tree  = geometry::build_bfs_forest(adj, frame.CA, true);
    auto bat   = geometry::compute_bat(st, tree, frame);

    bool threw = false;
    try {
        geometry::write_bat("/no/such/directory/out.bat", bat, tree);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw && "Expected runtime_error for bad write path");
}

// ---------------------------------------------------------------------------
int main() {
    test_write_read_roundtrip();
    test_throws_on_missing_file();
    test_throws_on_bad_write_path();

    std::cout << "PASS: bat_writer tests\n";
    return 0;
}
