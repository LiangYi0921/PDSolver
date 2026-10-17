#pragma once
#include "Element.h"
#include <algorithm>

namespace PD {

class TrussElement : public Element {
public:
    double E = 2.06e11;     // Elastic modulus (Pa)
    double A = 1.0e-4;      // Cross-sectional area (m^2)
    double N0 = 0.0;        // Initial pre-tension / prestress force (N)
    bool isCable = false;   // Tension-only cable element

    double L0 = 0.0;        // Initial length (m)
    double currentL = 0.0;  // Current deformed length (m)
    double axialForce = 0.0;// Current axial internal force (N) (tension > 0)
    double strain = 0.0;    // Current engineering strain

    // Local orientation frame
    Vector3 e1, e2, e3;

    TrussElement(int id_, int n1_, int n2_, double E_, double A_, double N0_ = 0.0, bool isCable_ = false)
        : E(E_), A(A_), N0(N0_), isCable(isCable_) {
        id = id_;
        node1 = n1_;
        node2 = n2_;
        type = isCable ? ElementType::Cable : ElementType::Truss;
    }

    void initialize(const std::vector<Particle>& particles) override {
        const auto& p1 = particles[node1].pos0;
        const auto& p2 = particles[node2].pos0;
        Vector3 d = p2 - p1;
        L0 = d.norm();
        currentL = L0;
        axialForce = N0;
        strain = 0.0;

        if (L0 > 1e-12) {
            e1 = d / L0;
        } else {
            e1 = {1.0, 0.0, 0.0};
        }
    }

    void computeInternalForces(std::vector<Particle>& particles) override {
        const auto& p1 = particles[node1].pos;
        const auto& p2 = particles[node2].pos;
        Vector3 d = p2 - p1;
        currentL = d.norm();

        if (currentL < 1e-12) {
            return;
        }

        e1 = d / currentL;
        strain = (currentL - L0) / L0;
        axialForce = E * A * strain + N0;

        // If tension-only cable, compressive axial force cannot be sustained
        if (isCable && axialForce < 0.0) {
            axialForce = 0.0;
        }

        Vector3 f1 = -axialForce * e1;
        Vector3 f2 =  axialForce * e1;

        particles[node1].fInt += f1;
        particles[node2].fInt += f2;
    }

    double getTranslationalStiffness() const override {
        double len = (L0 > 1e-6) ? L0 : 1.0;
        return (E * A) / len;
    }

    Vector3 getTranslationalStiffness3D() const override {
        double len = (L0 > 1e-6) ? L0 : 1.0;
        double k = (E * A) / len;
        double absX = std::abs(e1.x);
        double absY = std::abs(e1.y);
        double absZ = std::abs(e1.z);
        double sumAbs = absX + absY + absZ;
        return Vector3(k * absX * sumAbs, k * absY * sumAbs, k * absZ * sumAbs);
    }

    double getRotationalStiffness() const override {
        return 0.0;
    }
};

} // namespace PD

namespace LiangYiTest {
    using TrussElement = PD::TrussElement;
}
