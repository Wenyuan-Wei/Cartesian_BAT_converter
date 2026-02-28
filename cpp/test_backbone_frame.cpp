#include "structure/protein_root.hpp"   // or "structure/backbone_frame.hpp"
#include "structure/structure.hpp"

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

using geometry::Atom;
using geometry::Structure;
using geometry::Vec3;

static Atom make_atom(
    int index,
    const std::string& name,
    const std::string& resname,
    int resid,
    char chain,
    char alt_loc = ' ',
    double occ = 1.0
) {
    Atom a;
    a.index = index;
    a.name = name;
    a.residue_name = resname;
    a.residue_id = resid;
    a.chain_id = chain;
    a.element = (name == "N") ? "N" : "C";
    a.position = Vec3{0.0, 0.0, 0.0};
    a.alt_loc = alt_loc;
    a.occupancy = occ;
    return a;
}

static void test_picks_first_valid_residue_in_chain() {
    Structure st;

    // Residue 1 in chain A is missing C -> should skip
    st.atoms.push_back(make_atom(10, "N",  "GLY", 1, 'A'));
    st.atoms.push_back(make_atom(11, "CA", "GLY", 1, 'A'));

    // Residue 2 in chain A has N/CA/C -> should pick this
    st.atoms.push_back(make_atom(20, "N",  "ALA", 2, 'A'));
    st.atoms.push_back(make_atom(21, "CA", "ALA", 2, 'A'));
    st.atoms.push_back(make_atom(22, "C",  "ALA", 2, 'A'));

    // Another chain B also has a valid residue, but we asked for A
    st.atoms.push_back(make_atom(30, "N",  "SER", 1, 'B'));
    st.atoms.push_back(make_atom(31, "CA", "SER", 1, 'B'));
    st.atoms.push_back(make_atom(32, "C",  "SER", 1, 'B'));

    auto fr = geometry::select_first_backbone_frame(st, 'A');
    assert(fr.valid());

    // These are indices into st.atoms (vector positions), not Atom.index.
    // We inserted residue2's N/CA/C at positions 2,3,4
    assert(fr.N  == 2);
    assert(fr.CA == 3);
    assert(fr.C  == 4);
}

static void test_altloc_preference_space_over_A_over_others() {
    Structure st;

    // Residue 1 chain A: multiple CA candidates
    // Add backbone N and C normally
    st.atoms.push_back(make_atom(100, "N", "ALA", 1, 'A'));
    st.atoms.push_back(make_atom(101, "C", "ALA", 1, 'A'));

    // CA with altloc B (should NOT be chosen if ' ' or 'A' exists)
    st.atoms.push_back(make_atom(102, "CA", "ALA", 1, 'A', 'B', 1.0));

    // CA with altloc A
    st.atoms.push_back(make_atom(103, "CA", "ALA", 1, 'A', 'A', 1.0));

    // CA with altloc ' ' (best)
    st.atoms.push_back(make_atom(104, "CA", "ALA", 1, 'A', ' ', 1.0));

    auto fr = geometry::select_first_backbone_frame(st, 'A');
    assert(fr.valid());

    // N is at pos 0, C at pos 1, best CA at pos 4
    assert(fr.N  == 0);
    assert(fr.C  == 1);
    assert(fr.CA == 4);
}

static void test_occupancy_tie_breaking_within_preferred_altloc() {
    Structure st;

    st.atoms.push_back(make_atom(200, "N", "ALA", 1, 'A'));
    st.atoms.push_back(make_atom(201, "C", "ALA", 1, 'A'));

    // Both altloc ' ' so compare occupancy: higher should win
    st.atoms.push_back(make_atom(202, "CA", "ALA", 1, 'A', ' ', 0.40));
    st.atoms.push_back(make_atom(203, "CA", "ALA", 1, 'A', ' ', 0.90));

    auto fr = geometry::select_first_backbone_frame(st, 'A');
    assert(fr.valid());
    assert(fr.CA == 3); // the 0.90 occ candidate is at vector position 3
}

static void test_atom_index_tie_breaking_as_last_resort() {
    Structure st;

    st.atoms.push_back(make_atom(300, "N", "ALA", 1, 'A'));
    st.atoms.push_back(make_atom(301, "C", "ALA", 1, 'A'));

    // Same altloc + same occupancy => choose lower Atom.index
    // Insert higher Atom.index first (should lose)
    st.atoms.push_back(make_atom(999, "CA", "ALA", 1, 'A', ' ', 1.0));
    st.atoms.push_back(make_atom(100, "CA", "ALA", 1, 'A', ' ', 1.0));

    auto fr = geometry::select_first_backbone_frame(st, 'A');
    assert(fr.valid());

    // Lower Atom.index is 100, which we inserted at vector position 3
    assert(fr.CA == 3);
}

static void test_throws_if_chain_has_no_valid_frame() {
    Structure st;

    // Chain A: has N/CA but no C anywhere
    st.atoms.push_back(make_atom(10, "N",  "GLY", 1, 'A'));
    st.atoms.push_back(make_atom(11, "CA", "GLY", 1, 'A'));
    st.atoms.push_back(make_atom(20, "N",  "ALA", 2, 'A'));
    st.atoms.push_back(make_atom(21, "CA", "ALA", 2, 'A'));

    bool threw = false;
    try {
        (void)geometry::select_first_backbone_frame(st, 'A');
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw && "Expected runtime_error when no N/CA/C found in chain");
}

int main() {
    test_picks_first_valid_residue_in_chain();
    test_altloc_preference_space_over_A_over_others();
    test_occupancy_tie_breaking_within_preferred_altloc();
    test_atom_index_tie_breaking_as_last_resort();
    test_throws_if_chain_has_no_valid_frame();

    std::cout << "PASS: backbone_frame/protein_root selection tests\n";
    return 0;
}