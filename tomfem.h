#pragma once

#include <array>
#include <variant>
#include <vector>

namespace tomfem {

struct Node { float x, y; };

using NodeIndex = size_t;

using Element = std::array<NodeIndex, 3>;

struct Mesh {
    std::vector<Node> nodes;
    std::vector<Element> elements;
};

using BoundaryElement = std::array<NodeIndex, 2>;

// Return the boundary of the mesh.
std::vector<BoundaryElement> boundary(const Mesh& mesh);

// Hmm... some of these boundary conditions should just be for a point?

struct DirichletBoundaryCondition {
    float temperature;
};

struct NeumannBoundaryCondition {
    float flux;
};

struct ConvectionBoundaryCondition {
    float alpha;
    float ambient_temperature;
};

struct BoundaryCondition {
    BoundaryElement boundary;
    std::variant<DirichletBoundaryCondition,
                 NeumannBoundaryCondition,
                 ConvectionBoundaryCondition> condition;
};

using BoundaryConditions = std::vector<BoundaryCondition>;

// Return the nodal values.
std::vector<float>
solve_heat_equation(const Mesh& mesh,
                    const BoundaryConditions& boundary_conditions);

} // namespace tomfem
