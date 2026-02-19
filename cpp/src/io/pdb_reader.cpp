#include "geometry/pdb_reader.hpp"
#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <cctype>

namespace geometry {

static std::string trim(const std::string& s) {
    const auto start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    const auto end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

Structure PDBReader::read(const std::string& filename, char chain) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open PDB file: " + filename);
    }

    Structure st;
    std::string line;

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


        st.atoms.push_back(std::move(a));
    }

    return st.select_chain(chain);
}

} // namespace geometry
