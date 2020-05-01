#include "tomfem.h"

namespace tomfem {

std::vector<float>
solve_heat_equation(const Mesh& mesh,
                    const BoundaryConditions& /* boundary_conditions */) {

    std::vector<float> result(mesh.nodes.size(), 1.f);

    return result;
}

} // namespace tomfem
