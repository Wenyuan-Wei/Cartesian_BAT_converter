#include "structure/protein_root.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace geometry {

static bool altloc_preferred(char alt) {
    return (alt == ' ' || alt == 'A');
}

static bool is_backbone_name(const std::string& name) {
    return name == "N" || name == "CA" || name == "C";
}

static int altloc_rank(char alt) {
    if (alt == ' ') return 0;  // best
    if (alt == 'A') return 1;  // next best
    return 2;                  // all others
}

// Ranking: prefer altloc rank, then higher occupancy, then lower Atom.index.
static auto atom_candidate_key(const Atom& a) {
    return std::make_tuple(altloc_rank(a.alt_loc), -a.occupancy, a.index);
}

ProteinFrame select_first_backbone_frame(const Structure& st, std::optional<char> chain_opt) {
    if (st.atoms.empty()) {
        throw std::runtime_error("Structure has no atoms");
    }

    // Build set of chains present (deterministic order)
    std::vector<char> chains;
    chains.reserve(8);
    {
        std::unordered_map<char, bool> seen;
        for (const auto& a : st.atoms) {
            if (!seen[a.chain_id]) {
                seen[a.chain_id] = true;
                chains.push_back(a.chain_id);
            }
        }
        std::sort(chains.begin(), chains.end());
    }

    std::vector<char> chain_choices;
    if (chain_opt.has_value()) {
        chain_choices.push_back(*chain_opt);
    } else {
        chain_choices = chains;
    }

    // For each chain candidate, scan residues in increasing residue_id order.
    for (char chain_id : chain_choices) {

        // Map: residue_id -> best candidate atom for each backbone name
        // We'll store indices into st.atoms (atom vector position), not Atom.index.
        struct Picks {
            std::optional<std::size_t> N, CA, C;
        };

        // residue_id -> picks
        std::map<int, Picks> residue_picks; // map keeps residue_id sorted

        for (std::size_t i = 0; i < st.atoms.size(); ++i) {
            const Atom& a = st.atoms[i];
            if (a.chain_id != chain_id) continue;
            if (!is_backbone_name(a.name)) continue;

            auto& picks = residue_picks[a.residue_id];

            auto update_pick = [&](std::optional<std::size_t>& slot) {
                if (!slot.has_value()) {
                    slot = i;
                    return;
                }
                const Atom& cur = st.atoms[*slot];
                if (atom_candidate_key(a) < atom_candidate_key(cur)) {
                    slot = i;
                }
            };

            if (a.name == "N")  update_pick(picks.N);
            if (a.name == "CA") update_pick(picks.CA);
            if (a.name == "C")  update_pick(picks.C);
        }

        // Choose first residue_id in this chain that has all 3.
        for (const auto& [resid, picks] : residue_picks) {
            if (picks.N && picks.CA && picks.C) {
                ProteinFrame f;
                f.N  = *picks.N;
                f.CA = *picks.CA;
                f.C  = *picks.C;

                // Return indices into st.atoms; they match adjacency/bond indices in your code.
                return f;
            }
        }
    }

    // If we get here, no valid frame found.
    if (chain_opt.has_value()) {
        throw std::runtime_error("Could not find N/CA/C backbone frame in requested chain");
    }
    throw std::runtime_error("Could not find N/CA/C backbone frame in any chain");
}

} // namespace geometry