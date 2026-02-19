#include "vec3_ops.hpp"
#include "cartesian/vec3.hpp"

namespace geometry::vec3_ops {

    double dot(const geometry::Vec3& a, const geometry::Vec3& b) {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    geometry::Vec3 cross(const geometry::Vec3& a, const geometry::Vec3& b) {
        return {
            a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x
        };
    }

    double norm(const geometry::Vec3& v) {
        return std::sqrt(dot(v, v));
    }

    double distance(const geometry::Vec3& a, const geometry::Vec3& b) {
        return norm({a.x - b.x, a.y - b.y, a.z - b.z});
    }

} // namespace geometry::vec3_ops