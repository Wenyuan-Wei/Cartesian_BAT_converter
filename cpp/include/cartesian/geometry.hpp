#include "cartesian/vec3.hpp"
#pragma once

namespace geometry {

    double distance(const Vec3& a, const Vec3& b);

    // calculate angle given three points a, b, c with b as vertex
    double angle(const Vec3& a, const Vec3& b, const Vec3& c);

    // calculate dihedral angle given four points a, b, c, d
    double dihedral(const Vec3& a, const Vec3& b,
                    const Vec3& c, const Vec3& d);
}