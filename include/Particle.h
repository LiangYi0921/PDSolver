#pragma once
#include "Vector3.h"
#include "Matrix3x3.h"

namespace PD {

class Particle {
public:
    int id = 0;

    // Position coordinates
    Vector3 pos0;       // Initial coordinate (t = 0)
    Vector3 pos;        // Current coordinate x^n
    Vector3 posPrev;    // Previous coordinate x^{n-1}

    // Dynamic kinematics (Linear & Rotational)
    Vector3 vel;        // Current translational velocity v^n
    Vector3 acc;        // Current translational acceleration a^n
    Vector3 rot;        // Incremental rotation vector theta^n
    Vector3 rotPrev;    // Previous rotation vector
    Vector3 angVel;     // Current angular velocity omega^n
    Vector3 angAcc;     // Current angular acceleration alpha^n

    // Base triad orientation frame
    Matrix3x3 triad = Matrix3x3::identity();

    // Virtual mass & inertia (for explicit solvers / kinetic damping)
    Vector3 mass = {1.0, 1.0, 1.0};
    Vector3 invMass = {1.0, 1.0, 1.0};
    Vector3 inertia = {1.0, 1.0, 1.0};
    Vector3 invInertia = {1.0, 1.0, 1.0};

    // Elastic boundary springs (Translational kTrans & Rotational kRot)
    Vector3 springK = {0.0, 0.0, 0.0};
    Vector3 springKRot = {0.0, 0.0, 0.0};
    bool springCompressionOnlyZ = false;

    // Boundary condition constraints
    bool fixX = false;
    bool fixY = false;
    bool fixZ = false;
    bool fixRx = false;
    bool fixRy = false;
    bool fixRz = false;

    // Dynamic load accumulators
    Vector3 fExt;       // External applied force
    Vector3 fInt;       // Internal structural restoring resistance
    Vector3 mExt;       // External applied moment
    Vector3 mInt;       // Internal structural restoring moment

    Particle() = default;
    Particle(int id_, const Vector3& pos0_)
        : id(id_), pos0(pos0_), pos(pos0_), posPrev(pos0_) {}

    Vector3 displacement() const {
        return pos - pos0;
    }

    void resetInternalForces() {
        fInt.setZero();
        mInt.setZero();
    }

    Vector3 springForce() const {
        Vector3 disp = displacement();
        Vector3 fs(disp.x * springK.x, disp.y * springK.y, disp.z * springK.z);
        if (springCompressionOnlyZ && disp.z > 0.0) {
            fs.z = 0.0;
        }
        return fs;
    }

    Vector3 springMoment() const {
        return {rot.x * springKRot.x, rot.y * springKRot.y, rot.z * springKRot.z};
    }

    Vector3 unbalancedForce() const {
        Vector3 R = fExt - fInt - springForce();
        if (fixX) R.x = 0.0;
        if (fixY) R.y = 0.0;
        if (fixZ) R.z = 0.0;
        return R;
    }

    Vector3 unbalancedMoment() const {
        Vector3 Rm = mExt - mInt - springMoment();
        if (fixRx) Rm.x = 0.0;
        if (fixRy) Rm.y = 0.0;
        if (fixRz) Rm.z = 0.0;
        return Rm;
    }

    void setFixed(bool tx = true, bool ty = true, bool tz = true) {
        fixX = tx;
        fixY = ty;
        fixZ = tz;
    }

    void setFixedRot(bool rx = true, bool ry = true, bool rz = true) {
        fixRx = rx;
        fixRy = ry;
        fixRz = rz;
    }

    void setFullyFixed() {
        fixX = fixY = fixZ = true;
        fixRx = fixRy = fixRz = true;
    }

    void setSpring(const Vector3& kTrans, const Vector3& kRot = {0.0, 0.0, 0.0}, bool compressionOnlyZ = false) {
        springK = kTrans;
        springKRot = kRot;
        springCompressionOnlyZ = compressionOnlyZ;
    }
};

} // namespace PD

namespace LiangYiTest {
    using Particle = PD::Particle;
}
