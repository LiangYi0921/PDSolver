#pragma once
#include <vector>
#include "Particle.h"

namespace PD {

enum class ElementType {
    Truss,
    Beam,
    Cable
};

class Element {
public:
    int id = 0;
    int node1 = 0;
    int node2 = 0;
    ElementType type = ElementType::Truss;

    virtual ~Element() = default;

    virtual void initialize(const std::vector<Particle>& particles) = 0;

    // Computes internal nodal forces and moments and accumulates into particles[node1], particles[node2]
    virtual void computeInternalForces(std::vector<Particle>& particles) = 0;

    // Estimate maximum stiffness contribution for mass/time step stability
    virtual double getTranslationalStiffness() const = 0;
    virtual double getRotationalStiffness() const = 0;
    virtual Vector3 getTranslationalStiffness3D() const {
        double kt = getTranslationalStiffness();
        return Vector3(kt, kt, kt);
    }
    virtual Vector3 getRotationalStiffness3D() const {
        double kr = getRotationalStiffness();
        return Vector3(kr, kr, kr);
    }
};

} // namespace PD

namespace LiangYiTest {
    using ElementType = PD::ElementType;
    using Element = PD::Element;
}
