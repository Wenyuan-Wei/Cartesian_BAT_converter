#pragma once
#include "cartesian/vec3.hpp"
#include <cmath>

namespace geometry::vec3_ops {
    double dot(const geometry::Vec3& a, const geometry::Vec3& b);
    geometry::Vec3 cross(const geometry::Vec3& a, const geometry::Vec3& b);
    double norm(const geometry::Vec3& v);
    double distance(const geometry::Vec3& a, const geometry::Vec3& b);
} // namespace geometry::vec3_ops