// test_pdb_roundtrip.cpp — Integration test for the full BAT pipeline on a real PDB.
//
// Usage:
//   test_pdb_roundtrip <input.pdb> [chain]
//
// Exit codes:
//   0  all checks passed
//   1  any check failed (error written to stderr)
//
// Checks performed:
//   1. Full pipeline round-trip: PDB → compute_bat → reconstruct_cartesian
//      Every atom position must be recovered within ROUNDTRIP_TOL Å of original.
//   2. I/O round-trip: write_bat → read_bat → reconstruct_cartesian
//      Every atom position must still be within ROUNDTRIP_TOL Å of original.
//   3. Physical sanity: bond lengths, angles, and torsion values in valid ranges.
//   4. Count invariant: entries.size() == atoms.size() - 3.

#include "bat/bat_coords.hpp"
#include "bat/bat_reconstruct.hpp"
#include "io/bat_writer.hpp"
#include "io/pdb_reader.hpp"
#include "structure/bond_inference.hpp"
#include "structure/bfs_tree.hpp"
#include "structure/protein_root.hpp"
#include "structure/structure_graph.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

static constexpr double PI = 3.141592653589793238462643383;

// Round-trip tolerance: positions must recover to within this many Angstrom.
// Chosen to be tight (sub-nm), yet well above floating-point noise (~1e-13).
static constexpr double ROUNDTRIP_TOL = 1e-6;

// Physical sanity bounds on BAT values.
static constexpr double BOND_MIN  = 0.5;   // Å  (H-bonds ~0.95, covalent ~1.0–2.0)
static constexpr double BOND_MAX  = 3.0;   // Å  (disulfide S-S ~2.05 is the longest common)
static constexpr double ANGLE_MIN = 1e-3;  // rad  (never degenerate)
static constexpr double ANGLE_MAX = PI - 1e-3;

// ---------------------------------------------------------------------------
// Tiny accumulator for running statistics
// ---------------------------------------------------------------------------

struct Stats {
    double min_val =  1e300;
    double max_val = -1e300;
    double sum     = 0.0;
    double sum_sq  = 0.0;
    long long n    = 0;

    void update(double v) {
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
        sum    += v;
        sum_sq += v * v;
        ++n;
    }

    double mean() const { return n ? sum / n : 0.0; }
    double rmsd() const { return n ? std::sqrt(sum_sq / n) : 0.0; }
};

// ---------------------------------------------------------------------------
// Helper: compute per-atom position errors against original structure
// ---------------------------------------------------------------------------

static Stats position_errors(
    const std::vector<geometry::Vec3>& reconstructed,
    const geometry::Structure&         original,
    std::size_t*                       worst_idx = nullptr)
{
    Stats s;
    for (std::size_t i = 0; i < original.atoms.size(); ++i) {
        const auto& orig = original.atoms[i].position;
        double dx = reconstructed[i].x - orig.x;
        double dy = reconstructed[i].y - orig.y;
        double dz = reconstructed[i].z - orig.z;
        double err = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (worst_idx && err > s.max_val) *worst_idx = i;
        s.update(err);
    }
    return s;
}

// ---------------------------------------------------------------------------
// Main test logic
// ---------------------------------------------------------------------------

static int run(const std::string& pdb_path, std::optional<char> chain_opt) {

    const char chain_char = chain_opt.value_or('\0');
    int fail_count = 0;

    // -----------------------------------------------------------------------
    // Step 1: Read PDB
    // -----------------------------------------------------------------------
    geometry::Structure st;
    try {
        st = geometry::PDBReader::read(pdb_path, chain_char);
    } catch (const std::exception& e) {
        std::cerr << "FAIL  PDB read: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Atoms loaded      : " << st.atoms.size() << "\n";
    if (st.atoms.empty()) {
        std::cerr << "FAIL  No atoms loaded (wrong chain filter or empty file?)\n";
        return 1;
    }
    if (st.atoms.size() < 3) {
        std::cerr << "FAIL  Fewer than 3 atoms loaded; need at least N, CA, C\n";
        return 1;
    }

    // -----------------------------------------------------------------------
    // Step 2: Infer bonds
    // -----------------------------------------------------------------------
    geometry::BondInferOptions bond_opts;
    try {
        st.bonds = geometry::infer_bonds(st, bond_opts);
    } catch (const std::exception& e) {
        std::cerr << "FAIL  Bond inference: " << e.what() << "\n";
        return 1;
    }
    std::cout << "Bonds inferred    : " << st.bonds.size() << "\n";

    // -----------------------------------------------------------------------
    // Step 3: Build adjacency list and select backbone frame
    // -----------------------------------------------------------------------
    auto adj = geometry::build_adjacency(st);

    geometry::ProteinFrame frame = geometry::select_first_backbone_frame(st, chain_opt);
    if (!frame.valid()) {
        std::cerr << "FAIL  No backbone N/CA/C frame found "
                  << "(structure may lack standard backbone residues)\n";
        return 1;
    }
    std::cout << "Anchor indices    : N=" << frame.N
              << "  CA=" << frame.CA
              << "  C="  << frame.C  << "\n";
    std::cout << "Anchor atoms      : "
              << st.atoms[frame.N].name  << "(res " << st.atoms[frame.N].residue_id  << ")  "
              << st.atoms[frame.CA].name << "(res " << st.atoms[frame.CA].residue_id << ")  "
              << st.atoms[frame.C].name  << "(res " << st.atoms[frame.C].residue_id  << ")\n";

    // -----------------------------------------------------------------------
    // Step 4: BFS spanning forest + compute BAT
    // -----------------------------------------------------------------------
    auto tree = geometry::build_bfs_forest(adj, frame.CA, /*sort_neighbors=*/true);

    geometry::BATCoords bat;
    try {
        bat = geometry::compute_bat(st, tree, frame);
    } catch (const std::exception& e) {
        std::cerr << "FAIL  compute_bat: " << e.what() << "\n";
        return 1;
    }

    // -----------------------------------------------------------------------
    // CHECK 1: Entry count == atoms - 3
    // -----------------------------------------------------------------------
    {
        const std::size_t expected = st.atoms.size() - 3;
        if (bat.entries.size() != expected) {
            std::cerr << "FAIL  Entry count: expected " << expected
                      << "  got " << bat.entries.size() << "\n";
            ++fail_count;
        } else {
            std::cout << "CHECK entry_count : " << bat.entries.size()
                      << " == atoms(" << st.atoms.size() << ") - 3  PASS\n";
        }
    }

    // -----------------------------------------------------------------------
    // CHECK 2: In-memory round-trip (BAT → Cartesian)
    // -----------------------------------------------------------------------
    {
        std::vector<geometry::Vec3> pos;
        try {
            pos = geometry::reconstruct_cartesian(bat, tree);
        } catch (const std::exception& e) {
            std::cerr << "FAIL  reconstruct_cartesian: " << e.what() << "\n";
            return 1;
        }

        std::size_t worst = 0;
        Stats s = position_errors(pos, st, &worst);

        std::cout << "CHECK mem_roundtrip: "
                  << "max=" << s.max_val << " Å  "
                  << "mean=" << s.mean() << " Å  "
                  << "RMSD=" << s.rmsd() << " Å";

        if (s.max_val > ROUNDTRIP_TOL) {
            std::cout << "\n";
            std::cerr << "FAIL  In-memory round-trip: max error " << s.max_val
                      << " Å at atom " << worst
                      << " (" << st.atoms[worst].name
                      << " res " << st.atoms[worst].residue_id << ")\n";
            ++fail_count;
        } else {
            std::cout << "  PASS\n";
        }
    }

    // -----------------------------------------------------------------------
    // CHECK 3: I/O round-trip (write_bat → read_bat → reconstruct)
    // -----------------------------------------------------------------------
    {
        const std::string tmp = "/tmp/test_pdb_roundtrip_tmp.bat";
        bool io_ok = true;
        Stats s;

        try {
            geometry::write_bat(tmp, bat, tree);
            geometry::BATFile loaded = geometry::read_bat(tmp);
            std::remove(tmp.c_str());

            auto pos2 = geometry::reconstruct_cartesian(loaded.bat, loaded.tree);
            s = position_errors(pos2, st);
        } catch (const std::exception& e) {
            std::cerr << "FAIL  I/O round-trip: " << e.what() << "\n";
            io_ok = false;
            ++fail_count;
        }

        if (io_ok) {
            std::cout << "CHECK  io_roundtrip: "
                      << "max=" << s.max_val << " Å  "
                      << "mean=" << s.mean() << " Å  "
                      << "RMSD=" << s.rmsd() << " Å";

            if (s.max_val > ROUNDTRIP_TOL) {
                std::cout << "\n";
                std::cerr << "FAIL  I/O round-trip: max error " << s.max_val << " Å\n";
                ++fail_count;
            } else {
                std::cout << "  PASS\n";
            }
        }
    }

    // -----------------------------------------------------------------------
    // CHECK 4: Physical sanity of BAT values
    // -----------------------------------------------------------------------
    {
        Stats bond_s, angle_s, torsion_s;
        std::vector<std::size_t> bad_bonds, bad_angles, bad_torsions;

        for (const auto& e : bat.entries) {
            bond_s.update(e.bond);
            angle_s.update(e.angle);
            torsion_s.update(e.torsion);

            if (e.bond < BOND_MIN || e.bond > BOND_MAX)
                bad_bonds.push_back(e.atom_idx);
            if (e.angle <= ANGLE_MIN || e.angle >= ANGLE_MAX)
                bad_angles.push_back(e.atom_idx);
            if (e.torsion < -PI - 1e-9 || e.torsion > PI + 1e-9)
                bad_torsions.push_back(e.atom_idx);
        }

        auto deg = [](double r) { return r * (180.0 / PI); };

        std::cout << "Bond lengths (Å)   : "
                  << "min=" << bond_s.min_val
                  << "  max=" << bond_s.max_val
                  << "  mean=" << bond_s.mean() << "\n";
        std::cout << "Angles (deg)       : "
                  << "min=" << deg(angle_s.min_val)
                  << "  max=" << deg(angle_s.max_val)
                  << "  mean=" << deg(angle_s.mean()) << "\n";
        std::cout << "Torsions (deg)     : "
                  << "min=" << deg(torsion_s.min_val)
                  << "  max=" << deg(torsion_s.max_val) << "\n";

        bool phys_ok = true;

        if (!bad_bonds.empty()) {
            std::cerr << "FAIL  " << bad_bonds.size()
                      << " bond(s) outside (" << BOND_MIN << ", " << BOND_MAX << ") Å:";
            for (std::size_t idx : bad_bonds)
                std::cerr << "  atom" << idx << "(" << st.atoms[idx].name << ")";
            std::cerr << "\n";
            phys_ok = false;
        }
        if (!bad_angles.empty()) {
            std::cerr << "FAIL  " << bad_angles.size() << " degenerate angle(s).\n";
            phys_ok = false;
        }
        if (!bad_torsions.empty()) {
            std::cerr << "FAIL  " << bad_torsions.size() << " torsion(s) out of (-π, π].\n";
            phys_ok = false;
        }

        if (phys_ok)
            std::cout << "CHECK phys_sanity  : PASS\n";
        else
            ++fail_count;
    }

    return fail_count ? 1 : 0;
}

// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: test_pdb_roundtrip <input.pdb> [chain]\n";
        return 1;
    }

    const std::string pdb_path = argv[1];
    std::optional<char> chain_opt;
    if (argc >= 3 && argv[2][0] != '\0') {
        chain_opt = argv[2][0];
    }

    std::cout << "=== test_pdb_roundtrip";
    if (chain_opt) std::cout << "  chain=" << *chain_opt;
    std::cout << "  file=" << pdb_path << " ===\n";

    const int rc = run(pdb_path, chain_opt);
    std::cout << (rc == 0 ? "\nALL CHECKS PASSED\n" : "\nTEST FAILED\n");
    return rc;
}
