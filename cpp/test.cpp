#include <iostream>
#include "cartesian/geometry.hpp"
#include "io/pdb_reader.hpp"

constexpr double PI = 3.14159265358979323846;

// Test cases for geometry functions
/*
int main() {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 0.0, 0.0};
    Vec3 c{0.0, 1.0, 0.0};

    double theta = angle(a, b, c);

    std::cout << theta * 180.0 / PI << std::endl;
}
*/

/*
int main() {
    geometry::Vec3 a{1.0, 0.0, 0.0};
    geometry::Vec3 b{0.0, 0.0, 0.0};
    geometry::Vec3 c{0.0, 1.0, 0.0};
    geometry::Vec3 d{0.0, 1.0, 1.0};

    double phi = geometry::dihedral(a, b, c, d);
    std::cout << phi * 180.0 / PI << std::endl;
}
*/

/*
int main() {
    geometry::Vec3 a{1,0,0}, b{0,0,0}, c{0,1,0}, d{0,1,1};
    std::cout << geometry::dihedral(a,b,c,d) * 180.0 / PI << "\n";
}
*/


#include <iostream>

int main() {
    auto all = geometry::PDBReader::read("6crk.pdb");          // all chains
    auto A   = geometry::PDBReader::read("6crk.pdb", 'A');     // chain A only

    std::cout << "All atoms: " << all.atoms.size() << "\n";
    std::cout << "Chain A atoms: " << A.atoms.size() << "\n";

    if (!A.atoms.empty()) {
        const auto& a = A.atoms.front();
        std::cout << "First A atom: " << a.index << " " << a.name
                  << " " << a.residue_name << " " << a.residue_id
                  << " chain " << a.chain_id << "\n";
    }
}