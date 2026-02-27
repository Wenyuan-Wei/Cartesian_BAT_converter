#include "io/pdb_reader.hpp"
#include "structure/bond_inference.hpp"
#include "structure/structure.hpp"
#include "structure/structure_graph.hpp"
#include "chem/valence.hpp"

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <stdexcept>

static inline std::uint64_t bond_key(std::size_t i, std::size_t j) {
    if (j < i) std::swap(i, j);
    return (static_cast<std::uint64_t>(i) << 32) | static_cast<std::uint64_t>(j);
}

struct ResKey {
    char chain;
    int resid;
    bool operator<(const ResKey& o) const {
        return (chain < o.chain) || (chain == o.chain && resid < o.resid);
    }
};

static void require(bool cond, const std::string& msg) {
    if (!cond) throw std::runtime_error(msg);
}

int main(int argc, char** argv) {
    try {
        if (argc < 3) {
            std::cerr << "Usage: test_bond_inference <pdb_file> <chain>\n";
            return 2;
        }

        const std::string pdb = argv[1];
        const char chain = argv[2][0];

        geometry::Structure st = geometry::PDBReader::read(pdb, chain);

        geometry::BondInferOptions opt;
        opt.include_hydrogen = true;     // ok even if no H
        opt.enforce_valence = true;

        st.bonds = geometry::infer_bonds(st, opt);
        auto adj = geometry::build_adjacency(st);

        std::cout << "Atoms: " << st.atoms.size() << "\n";
        std::cout << "Bonds: " << st.bonds.size() << "\n";

        // --- Sanity prints: neighbor counts for a few atoms ---
        auto print_atom = [&](std::size_t i) {
            const auto& a = st.atoms[i];
            std::cout << "Atom[" << i << "] "
                      << a.residue_name << " " << a.chain_id << a.residue_id << " " << a.name
                      << " el=" << a.element
                      << " deg=" << adj[i].size()
                      << "\n";
        };

        if (!st.atoms.empty()) {
            print_atom(0);
            print_atom(st.atoms.size() / 2);
            print_atom(st.atoms.size() - 1);
        }

        // --- Max degree vs valence check ---
        std::size_t max_deg = 0;
        std::size_t max_deg_idx = 0;

        for (std::size_t i = 0; i < st.atoms.size(); ++i) {
            const auto& a = st.atoms[i];
            require(!a.element.empty(), "Empty element on atom index " + std::to_string(i));

            const int vmax = chem::max_valence(a.element);
            const std::size_t deg = adj[i].size();

            if (deg > max_deg) { max_deg = deg; max_deg_idx = i; }

            if (deg > static_cast<std::size_t>(vmax)) {
                const auto& aa = st.atoms[i];
                throw std::runtime_error(
                    "Valence violation: " + aa.residue_name + " " +
                    std::string(1, aa.chain_id) + std::to_string(aa.residue_id) + " " + aa.name +
                    " el=" + aa.element +
                    " degree=" + std::to_string(deg) +
                    " > max_valence=" + std::to_string(vmax)
                );
            }
        }

        std::cout << "Max degree: " << max_deg << " at atom index " << max_deg_idx << "\n";
        if (!st.atoms.empty()) print_atom(max_deg_idx);

        // --- Build bondset for fast membership tests ---
        std::unordered_set<std::uint64_t> bondset;
        bondset.reserve(st.bonds.size() * 2);
        for (const auto& b : st.bonds) bondset.insert(bond_key(b.i, b.j));

        auto has_bond = [&](std::size_t i, std::size_t j) {
            return bondset.count(bond_key(i, j)) != 0;
        };

        // --- Residue -> atom-name -> atom-index map ---
        std::map<ResKey, std::unordered_map<std::string, std::size_t>> residues;
        residues.clear();
        for (std::size_t i = 0; i < st.atoms.size(); ++i) {
            const auto& a = st.atoms[i];
            residues[{a.chain_id, a.residue_id}][a.name] = i;
        }

        // --- Backbone invariant checks ---
        int backbone_failures = 0;
        int backbone_checked = 0;

        for (const auto& [rk, amap] : residues) {
            if (amap.count("N") && amap.count("CA") && amap.count("C") && amap.count("O")) {
                backbone_checked++;
                const std::size_t iN  = amap.at("N");
                const std::size_t iCA = amap.at("CA");
                const std::size_t iC  = amap.at("C");
                const std::size_t iO  = amap.at("O");

                if (!has_bond(iN, iCA)) backbone_failures++;
                if (!has_bond(iCA, iC)) backbone_failures++;
                if (!has_bond(iC, iO)) backbone_failures++;
            }
        }

        std::cout << "Backbone residues checked: " << backbone_checked << "\n";
        std::cout << "Backbone bond failures: " << backbone_failures << "\n";
        require(backbone_failures == 0, "Backbone invariant failed (missing N-CA/CA-C/C-O bonds).");

        // --- Peptide bond check: C(i) -- N(i+1) ---
        int peptide_failures = 0;
        int peptide_checked = 0;

        for (auto it = residues.begin(); it != residues.end(); ++it) {
            auto next = std::next(it);
            if (next == residues.end()) break;

            if (it->first.chain != next->first.chain) continue;
            if (next->first.resid != it->first.resid + 1) continue;

            const auto& amap_i = it->second;
            const auto& amap_j = next->second;

            if (amap_i.count("C") && amap_j.count("N")) {
                peptide_checked++;
                const std::size_t iC = amap_i.at("C");
                const std::size_t jN = amap_j.at("N");
                if (!has_bond(iC, jN)) peptide_failures++;
            }
        }

        std::cout << "Peptide bonds checked: " << peptide_checked << "\n";
        std::cout << "Peptide bond failures: " << peptide_failures << "\n";
        require(peptide_failures == 0, "Peptide invariant failed (missing C(i)-N(i+1) bonds).");

        std::cout << "PASS: bond inference + graph sanity checks\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}