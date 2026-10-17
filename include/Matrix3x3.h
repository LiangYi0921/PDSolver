#pragma once
#include "Vector3.h"

namespace PD {

struct Matrix3x3 {
    // Row-major storage: m[row][col]
    double m[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    Matrix3x3() = default;

    Matrix3x3(double m00, double m01, double m02,
              double m10, double m11, double m12,
              double m20, double m21, double m22) {
        m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
        m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
        m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
    }

    static Matrix3x3 identity() {
        return Matrix3x3();
    }

    static Matrix3x3 zero() {
        Matrix3x3 res;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                res.m[r][c] = 0.0;
        return res;
    }

    Vector3 operator*(const Vector3& v) const {
        return {
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
        };
    }

    Matrix3x3 operator*(const Matrix3x3& rhs) const {
        Matrix3x3 res = zero();
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                for (int k = 0; k < 3; ++k) {
                    res.m[r][c] += m[r][k] * rhs.m[k][c];
                }
            }
        }
        return res;
    }

    Matrix3x3 operator+(const Matrix3x3& rhs) const {
        Matrix3x3 res;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                res.m[r][c] = m[r][c] + rhs.m[r][c];
        return res;
    }

    Matrix3x3 operator*(double s) const {
        Matrix3x3 res;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                res.m[r][c] = m[r][c] * s;
        return res;
    }

    Matrix3x3 transpose() const {
        return {
            m[0][0], m[1][0], m[2][0],
            m[0][1], m[1][1], m[2][1],
            m[0][2], m[1][2], m[2][2]
        };
    }

    double determinant() const {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
               m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    Matrix3x3 inverse() const {
        double det = determinant();
        if (std::abs(det) < 1e-15) {
            return identity();
        }
        double invDet = 1.0 / det;
        Matrix3x3 inv;
        inv.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * invDet;
        inv.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet;
        inv.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet;

        inv.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet;
        inv.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet;
        inv.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * invDet;

        inv.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * invDet;
        inv.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * invDet;
        inv.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * invDet;
        return inv;
    }

    // Outer product: a * b^T
    static Matrix3x3 outerProduct(const Vector3& a, const Vector3& b) {
        return {
            a.x * b.x, a.x * b.y, a.x * b.z,
            a.y * b.x, a.y * b.y, a.y * b.z,
            a.z * b.x, a.z * b.y, a.z * b.z
        };
    }

    // Rodrigues rotation matrix from rotation vector (axial vector)
    static Matrix3x3 exp(const Vector3& omega) {
        double theta = omega.norm();
        if (theta < 1e-12) {
            return identity();
        }
        Vector3 k = omega / theta;
        Matrix3x3 K(
             0.0, -k.z,  k.y,
             k.z,  0.0, -k.x,
            -k.y,  k.x,  0.0
        );
        Matrix3x3 K2 = K * K;
        double s = std::sin(theta);
        double c = 1.0 - std::cos(theta);
        return identity() + (K * s) + (K2 * c);
    }

    Vector3 col(int c) const {
        return {m[0][c], m[1][c], m[2][c]};
    }

    Vector3 row(int r) const {
        return {m[r][0], m[r][1], m[r][2]};
    }
};

} // namespace PD

namespace LiangYiTest {
    using Matrix3x3 = PD::Matrix3x3;
}
