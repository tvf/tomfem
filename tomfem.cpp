#include "tomfem.h"

#include "solve_linear_system.h"

#include <cassert>

namespace tomfem {

//#include <iostream>
//void check_solve_linear_system() {
//
//    //Matrix A { 3, 3, { 1, 2, 3,   4, 5, 6,   7, 8, 10 } };
//    Matrix A { 3, 2, { 1, 2,   3, 4,   5, 6 } };
//    //Matrix A { 3, 3, { 0, 0, 0,   0, 0, 0,   0, 0, 0 } };
//    Matrix b { 3, 1, { 3,   3,   4 } };
//
//    Matrix solution = solve_linear_system(A, b);
//
//    std::cout << solution.num_rows  << " rows of "
//              << solution.num_columns  << " each\n";
//
//    for (size_t i = 0; i < solution.num_rows; ++i) {
//        std::cout << '[' << ' ';
//        for (size_t j = 0; j < solution.num_columns; ++j) {
//            std::cout << solution.at(i, j) << ' ';
//        }
//        std::cout << ']' << '\n';
//    }
//}


// OK what's the plan here?
//
// For now assume D = I
//
// [ ] calculate stiffness matrix K from the mesh alone
// [ ] calculate Kc according to any convection boundary conditions
//
// [ ] calculate the term from the Neumann boundary conditions
// [ ] calculate the term from the convection boundary conditions
// [ ] take the Dirichlet boundary conditions into account by doctoring K and f
//
// [X] don't bother with the load vector for now but write a function for it
// [X] solve linear system with Eigen

void add_load_vector(Matrix& f) {
    assert(f.num_columns == 1);
    // for (size_t i = 0; i < f.num_rows; ++i) {
    //     f.at(i, 0) += integral over domain of (Q Ni(x, y))
    // }
}

std::vector<float>
solve_heat_equation(const Mesh& mesh,
                    const BoundaryConditions& /* boundary_conditions */) {

    std::vector<float> result(mesh.nodes.size(), 1.f);

    //check_solve_linear_system();

    return result;
}

} // namespace tomfem
