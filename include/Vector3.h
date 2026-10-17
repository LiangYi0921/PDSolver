#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>

namespace PD {

struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vector3() = default;
    constexpr Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Vector arithmetic
    Vector3 operator+(const Vector3& rhs) const {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    Vector3 operator-(const Vector3& rhs) const {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }

    Vector3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    friend Vector3 operator*(double scalar, const Vector3& v) {
        return {v.x * scalar, v.y * scalar, v.z * scalar};
    }

    Vector3 operator/(double scalar) const {
        double inv = 1.0 / scalar;
        return {x * inv, y * inv, z * inv};
    }

    Vector3& operator+=(const Vector3& rhs) {
        x += rhs.x; y += rhs.y; z += rhs.z;
        return *this;
    }

    Vector3& operator-=(const Vector3& rhs) {
        x -= rhs.x; y -= rhs.y; z -= rhs.z;
        return *this;
    }

    Vector3& operator*=(double scalar) {
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    Vector3& operator/=(double scalar) {
        double inv = 1.0 / scalar;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    Vector3 operator-() const {
        return {-x, -y, -z};
    }

    // Vector operations
    double dot(const Vector3& rhs) const {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    Vector3 cross(const Vector3& rhs) const {
        return {
            y * rhs.z - z * rhs.y,
            z * rhs.x - x * rhs.z,
            x * rhs.y - y * rhs.x
        };
    }

    double normSq() const {
        return x * x + y * y + z * z;
    }

    double norm() const {
        return std::sqrt(normSq());
    }

    Vector3 normalized() const {
        double n = norm();
        if (n > 1e-12) {
            return *this / n;
        }
        return {0.0, 0.0, 0.0};
    }

    void setZero() {
        x = 0.0; y = 0.0; z = 0.0;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector3& v) {
        os << "[" << std::setw(10) << v.x << ", " << std::setw(10) << v.y << ", " << std::setw(10) << v.z << "]";
        return os;
    }
};

} // namespace PD

namespace LiangYiTest {
    using Vector3 = PD::Vector3;
}
