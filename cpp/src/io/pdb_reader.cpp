#include "io/pdb_reader.hpp"
#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <cctype>

namespace geometry {

static std::string trim(const std::string& s) {
    const auto start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    const auto end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

struct AtomKey {
    char chain;
    int resid;
    std::string name;

    bool operator==(const AtomKey& other) const {
        return chain == other.chain &&
               resid == other.resid &&
               name  == other.name;
    }
};

struct AtomKeyHash {
    std::size_t operator()(const AtomKey& k) const noexcept {
        std::size_t h = 1469598103934665603ull;
        auto mix = [&](std::size_t v) {
            h ^= v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        };
        mix(static_cast<std::size_t>(k.chain));
        mix(static_cast<std::size_t>(k.resid));
        mix(std::hash<std::string>{}(k.name));
        return h;
    }
};

static int altloc_rank(char alt) {
    // Lower is better
    if (alt == ' ' || alt == '\0') return 0;
    if (alt == 'A') return 1;
    return 2;
}

static bool nearly_same_pos(const geometry::Vec3& p, const geometry::Vec3& q, double eps = 1e-3) {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    const double dz = p.z - q.z;
    return (dx*dx + dy*dy + dz*dz) <= eps*eps;
}

static std::string infer_element_from_atom_field(const std::string& atom_field_4) {
    // atom_field_4 is the raw 4-char field (columns 13-16), NOT trimmed
    if (atom_field_4.size() != 4) return "";

    char c0 = atom_field_4[0];
    char c1 = atom_field_4[1];

    auto is_alpha = [](unsigned char c){ return std::isalpha(c) != 0; };

    char e = '\0';
    if (is_alpha((unsigned char)c0)) {
        e = c0;
    } else if (is_alpha((unsigned char)c1)) {
        e = c1;
    } else {
        return "";
    }

    e = (char)std::toupper((unsigned char)e);
    return std::string(1, e);
}

Structure PDBReader::read(const std::string& filename, char chain) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open PDB file: " + filename);
    }

    Structure st;
    std::string line;
    std::unordered_map<AtomKey, std::size_t, AtomKeyHash> seen;

    while (std::getline(file, line)) {
        if (line.size() < 54) continue; // need at least through Z coord

        // PDB record name is columns 1-6
        const std::string rec = line.substr(0, 6);
        // Only load ATOM records (Downstream pipelines do not support HETATM)
        if (rec != "ATOM  ") continue;

        Atom a{};

        // Fixed-width PDB columns
        // serial: 7-11
        a.index = std::stoi(line.substr(6, 5));

        // atom name: 13-16
        a.name = trim(line.substr(12, 4));

        // alt loc: 17
        a.alt_loc = (line.size() > 16) ? line[16] : ' ';

        // residue name: 18-20
        a.residue_name = trim(line.substr(17, 3));

        // chain id: 22
        a.chain_id = (line.size() > 21) ? line[21] : ' ';

        // residue seq: 23-26
        a.residue_id = std::stoi(line.substr(22, 4));

        // x,y,z: 31-38, 39-46, 47-54
        const double x = std::stod(line.substr(30, 8));
        const double y = std::stod(line.substr(38, 8));
        const double z = std::stod(line.substr(46, 8));
        a.position = Vec3{x, y, z};

        // atom occupancy: 55-60 (default 1.0)
        a.occupancy = (line.size() >= 60) ? std::stod(line.substr(54, 6)) : 1.0;

        // element: columns 77-78 (optional in PDB)
        std::string element = "";
        if (line.size() >= 78) {
            element = trim(line.substr(76, 2));
        }

        // If element column is empty, infer from atom name
        if (element.empty()) {
        
            // PDB atom name is exactly 4 characters (columns 13-16)
            const std::string raw_name = line.substr(12, 4);
        
            if (raw_name.size() != 4) {
                throw std::runtime_error(
                    "Malformed atom name field while inferring element: '" + raw_name + "'"
                );
            }
        
            // Rule: if first character is space → element is second character
            // else → element is first character
            char c0 = raw_name[0];
            char inferred = '\0';
        
            if (std::isalpha(c0)) {
                inferred = c0;
            } else {
                if (std::isalpha(raw_name[1])) {
                    inferred = raw_name[1];
                } else {
                    throw std::runtime_error(
                        "Cannot infer element from atom name: '" + raw_name +
                        "' (residue " + a.residue_name +
                        " " + std::to_string(a.residue_id) + ")"
                    );
                }
            }
        
            inferred = std::toupper(inferred);
            element = std::string(1, inferred);
        }

        // Normalize element to uppercase
        for (char& c : element) {
            c = std::toupper(c);
        }

        // Strict protein-only element validation
        static const std::unordered_set<std::string> allowed = {
            "H", "C", "N", "O", "S"
        };

        if (!allowed.count(element)) {
            throw std::runtime_error(
                "Unsupported element '" + element +
                "' in residue " + a.residue_name +
                " " + std::to_string(a.residue_id) +
                " (atom " + a.name + ")"
            );
        }

        a.element = element;

        AtomKey key{a.chain_id, a.residue_id, a.name};

        auto it = seen.find(key);
        if (it == seen.end()) {
            seen.emplace(std::move(key), st.atoms.size());
            st.atoms.push_back(std::move(a));
        } else {
            Atom& existing = st.atoms[it->second];
        
            if (a.occupancy > existing.occupancy) {
                existing = std::move(a);
            } else if (a.occupancy < existing.occupancy) {
                // keep existing
            } else {
                // equal occupancy: tie-break
                const int r_new = altloc_rank(a.alt_loc);
                const int r_old = altloc_rank(existing.alt_loc);
            
                if (r_new < r_old) {
                    existing = std::move(a);
                } else if (r_new > r_old) {
                    // keep existing
                } else if (nearly_same_pos(a.position, existing.position)) {
                    // exact duplicate (same altLoc rank + same coords) -> keep existing
                } else {
                    throw std::runtime_error(
                        "Duplicate atom with equal occupancy and ambiguous altLoc/coords for " +
                        a.residue_name + " " + std::string(1, a.chain_id) + std::to_string(a.residue_id) +
                        " " + a.name +
                        " (occ=" + std::to_string(a.occupancy) +
                        ", altLoc='" + std::string(1, a.alt_loc) + "')"
                    );
                }
            }
        }
    }

    return st.select_chain(chain);
}

} // namespace geometry
