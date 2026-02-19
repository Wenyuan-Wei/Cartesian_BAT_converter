#include "structure/structure.hpp"

namespace geometry {
    Structure Structure::select_chain(char chain) const {
        if (chain == '\0') return *this;

        Structure out;
        out.atoms.reserve(atoms.size());
        for (const auto& a : atoms) {
            if (a.chain_id == chain) out.atoms.push_back(a);
        }
        out.atoms.shrink_to_fit();
        return out;
    }
}