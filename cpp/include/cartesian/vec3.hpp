#pragma once

namespace geometry {
    struct Vec3 {
        double x, y, z;

        Vec3() = default;
        Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    };
}