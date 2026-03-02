#include "io/bat_writer.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace geometry {

// ---------------------------------------------------------------------------
// File format
//
// BAT 1
// NATOMS <n>
// ANCHOR_N  <N_idx>  <x> <y> <z>
// ANCHOR_CA <CA_idx> <x> <y> <z>
// ANCHOR_C  <C_idx>  <x> <y> <z>
// ENTRY <atom_idx> <p> <g> <h>  <bond> <angle> <torsion>
// ...
//
// Integers p/g/h use -1 to represent BFSTree::npos.
// Doubles are written at precision 17 for lossless IEEE 754 round-trip.
// ---------------------------------------------------------------------------

namespace {

inline long long to_signed(std::size_t v) {
    return (v == BFSTree::npos) ? -1LL : static_cast<long long>(v);
}

inline std::size_t from_signed(long long v) {
    return (v == -1LL) ? BFSTree::npos : static_cast<std::size_t>(v);
}

} // anonymous namespace

// ---------------------------------------------------------------------------
void write_bat(const std::string& filename,
               const BATCoords&   bat,
               const BFSTree&     tree) {
    std::ofstream f(filename);
    if (!f)
        throw std::runtime_error("write_bat: cannot open \"" + filename + "\"");

    f << std::setprecision(17);

    f << "BAT 1\n";
    f << "NATOMS " << tree.n << "\n";

    const auto& a = bat.anchor;
    f << "ANCHOR_N  " << a.N_idx  << "  "
      << a.N_pos.x  << "  " << a.N_pos.y  << "  " << a.N_pos.z  << "\n";
    f << "ANCHOR_CA " << a.CA_idx << "  "
      << a.CA_pos.x << "  " << a.CA_pos.y << "  " << a.CA_pos.z << "\n";
    f << "ANCHOR_C  " << a.C_idx  << "  "
      << a.C_pos.x  << "  " << a.C_pos.y  << "  " << a.C_pos.z  << "\n";

    for (const auto& e : bat.entries) {
        const BFSTree::Ref3& r = tree.refs[e.atom_idx];
        f << "ENTRY "
          << e.atom_idx        << " "
          << to_signed(r.p)    << " "
          << to_signed(r.g)    << " "
          << to_signed(r.h)    << "  "
          << e.bond            << "  "
          << e.angle           << "  "
          << e.torsion         << "\n";
    }
}

// ---------------------------------------------------------------------------
BATFile read_bat(const std::string& filename) {
    std::ifstream f(filename);
    if (!f)
        throw std::runtime_error("read_bat: cannot open \"" + filename + "\"");

    BATFile result;
    BATCoords& bat  = result.bat;
    BFSTree&   tree = result.tree;

    bool got_version = false;
    bool got_natoms  = false;
    bool got_N = false, got_CA = false, got_C = false;

    auto parse_error = [&](const std::string& msg) {
        throw std::runtime_error("read_bat(\"" + filename + "\"): " + msg);
    };

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        std::string tag;
        if (!(ss >> tag)) continue;  // whitespace-only line

        if (tag == "BAT") {
            int ver;
            ss >> ver;
            if (!ss || ver != 1)
                parse_error("unsupported version (expected BAT 1)");
            got_version = true;
        }
        else if (tag == "NATOMS") {
            ss >> tree.n;
            if (!ss) parse_error("malformed NATOMS line");
            tree.refs.assign(tree.n,
                BFSTree::Ref3{BFSTree::npos, BFSTree::npos, BFSTree::npos});
            got_natoms = true;
        }
        else if (tag == "ANCHOR_N" || tag == "ANCHOR_CA" || tag == "ANCHOR_C") {
            if (!got_natoms) parse_error("ANCHOR line before NATOMS");
            long long idx;
            double x, y, z;
            ss >> idx >> x >> y >> z;
            if (!ss) parse_error("malformed " + tag + " line");
            if (tag == "ANCHOR_N") {
                bat.anchor.N_idx  = from_signed(idx);
                bat.anchor.N_pos  = {x, y, z};
                got_N  = true;
            } else if (tag == "ANCHOR_CA") {
                bat.anchor.CA_idx = from_signed(idx);
                bat.anchor.CA_pos = {x, y, z};
                got_CA = true;
            } else {
                bat.anchor.C_idx  = from_signed(idx);
                bat.anchor.C_pos  = {x, y, z};
                got_C  = true;
            }
        }
        else if (tag == "ENTRY") {
            if (!got_N || !got_CA || !got_C)
                parse_error("ENTRY before all ANCHOR lines");
            long long atom_idx, p, g, h;
            double bond, angle, torsion;
            ss >> atom_idx >> p >> g >> h >> bond >> angle >> torsion;
            if (!ss) parse_error("malformed ENTRY line");
            std::size_t v = static_cast<std::size_t>(atom_idx);
            if (v >= tree.n) parse_error("ENTRY atom_idx out of range");
            bat.entries.push_back({v, bond, angle, torsion});
            tree.refs[v] = {from_signed(p), from_signed(g), from_signed(h)};
        }
        else {
            parse_error("unknown tag \"" + tag + "\"");
        }
    }

    if (!got_version)         parse_error("missing BAT header line");
    if (!got_natoms)          parse_error("missing NATOMS line");
    if (!got_N || !got_CA || !got_C) parse_error("missing one or more ANCHOR lines");

    return result;
}

} // namespace geometry
