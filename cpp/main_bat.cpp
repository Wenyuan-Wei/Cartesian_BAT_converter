// main_bat.cpp — CLI driver for the Cartesian → BAT converter.
//
// Usage:
//   bat_convert  <input.pdb>  <output.bat>  [--chain <C>]
//
// Runs the full pipeline:
//   PDBReader::read → infer_bonds → build_adjacency
//   → select_first_backbone_frame → build_bfs_forest
//   → compute_bat → write_bat
//
// Exits 0 on success, 1 on any error (message written to stderr).

#include "bat/bat_coords.hpp"
#include "io/bat_writer.hpp"
#include "io/pdb_reader.hpp"
#include "structure/bond_inference.hpp"
#include "structure/bfs_tree.hpp"
#include "structure/protein_root.hpp"
#include "structure/structure_graph.hpp"

#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// Argument parsing
// ---------------------------------------------------------------------------

struct Args {
    std::string input_pdb;
    std::string output_bat;
    std::optional<char> chain;
};

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <input.pdb> <output.bat> [--chain <C>]\n";
}

static Args parse_args(int argc, char* argv[]) {
    if (argc < 3) {
        usage(argv[0]);
        std::exit(1);
    }

    Args a;
    a.input_pdb  = argv[1];
    a.output_bat = argv[2];

    for (int i = 3; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--chain") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --chain requires an argument.\n";
                std::exit(1);
            }
            std::string val = argv[++i];
            if (val.size() != 1) {
                std::cerr << "Error: --chain argument must be a single character.\n";
                std::exit(1);
            }
            a.chain = val[0];
        } else {
            std::cerr << "Error: unknown option '" << flag << "'.\n";
            usage(argv[0]);
            std::exit(1);
        }
    }
    return a;
}

// ---------------------------------------------------------------------------
// Pipeline
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    const Args args = parse_args(argc, argv);

    try {
        // Step 1: read PDB
        const char chain_hint = args.chain.value_or('\0');
        geometry::Structure st = geometry::PDBReader::read(args.input_pdb, chain_hint);
        std::cout << "Read " << st.atoms.size() << " atoms from "
                  << args.input_pdb << "\n";

        // Step 2: infer bonds
        geometry::BondInferOptions bond_opts;
        st.bonds = geometry::infer_bonds(st, bond_opts);
        std::cout << "Inferred " << st.bonds.size() << " bonds.\n";

        // Step 3: build adjacency list
        auto adj = geometry::build_adjacency(st);

        // Step 4: select backbone frame (N/CA/C anchor)
        geometry::ProteinFrame frame =
            geometry::select_first_backbone_frame(st, args.chain);
        if (!frame.valid()) {
            throw std::runtime_error(
                "Could not find a backbone N/CA/C frame in the structure.");
        }
        std::cout << "Backbone anchor: atom indices N=" << frame.N
                  << "  CA=" << frame.CA
                  << "  C="  << frame.C << "\n";

        // Step 5: BFS spanning forest rooted at CA
        auto tree = geometry::build_bfs_forest(adj, frame.CA, /*sort_neighbors=*/true);

        // Step 6: compute BAT coordinates
        geometry::BATCoords bat = geometry::compute_bat(st, tree, frame);
        std::cout << "Computed BAT coordinates for " << bat.entries.size()
                  << " non-anchor atoms.\n";

        // Step 7: write output
        geometry::write_bat(args.output_bat, bat, tree);
        std::cout << "Wrote BAT file to " << args.output_bat << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
