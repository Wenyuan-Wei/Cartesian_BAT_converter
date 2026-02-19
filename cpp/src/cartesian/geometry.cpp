#include "geometry/vec3.hpp"
#include "vec3_ops.hpp"
#include <cmath>
#include <algorithm>

namespace geometry {

double distance(const Vec3& a, const Vec3& b) {
    return vec3_ops::distance(a, b);
}

double angle(const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 ba{a.x - b.x, a.y - b.y, a.z - b.z};
    Vec3 bc{c.x - b.x, c.y - b.y, c.z - b.z};

    double nba = vec3_ops::norm(ba);
    double nbc = vec3_ops::norm(bc);
    if (nba == 0.0 || nbc == 0.0) return 0.0;

    double cos_theta = vec3_ops::dot(ba, bc) / (nba * nbc);
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);
    return std::acos(cos_theta);
}

double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    Vec3 b1{b.x - a.x, b.y - a.y, b.z - a.z};
    Vec3 b2{c.x - b.x, c.y - b.y, c.z - b.z};
    Vec3 b3{d.x - c.x, d.y - c.y, d.z - c.z};

    Vec3 n1 = vec3_ops::cross(b1, b2);
    Vec3 n2 = vec3_ops::cross(b2, b3);

    double b2_norm = vec3_ops::norm(b2);
    if (b2_norm == 0.0) return 0.0;

    Vec3 b2_hat{b2.x / b2_norm, b2.y / b2_norm, b2.z / b2_norm};

    double x = vec3_ops::dot(n1, n2);
    double y = vec3_ops::dot(vec3_ops::cross(n1, n2), b2_hat);

    return std::atan2(y, x);
}

} // namespace geometry
