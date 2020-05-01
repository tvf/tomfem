#include "tomfem.h"

#include "solve_linear_system.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ostream>
#include <iomanip>
#include <iostream>

namespace tomfem {

void log_matrix(std::ostream& out, const Matrix& matrix) {

    out << matrix.num_rows  << " rows of "
        << matrix.num_columns  << " each\n";

    for (size_t i = 0; i < matrix.num_rows; ++i) {
        out << '[' << ' ';

        for (size_t j = 0; j < matrix.num_columns; ++j) {
            out << std::setw(5) << matrix.at(i, j) << ' ';
        }

        out << ']' << '\n';
    }
}

// void check_solve_linear_system() {
//
//     //Matrix A { 3, 3, { 1, 2, 3,   4, 5, 6,   7, 8, 10 } };
//     Matrix A { 3, 2, { 1, 2,   3, 4,   5, 6 } };
//     //Matrix A { 3, 3, { 0, 0, 0,   0, 0, 0,   0, 0, 0 } };
//     Matrix b { 3, 1, { 3,   3,   4 } };
//
//     Matrix solution = solve_linear_system(A, b);
//
//     log_matrix(std::cout, solution);
// }

// For now assume D = I
//
// [X] calculate stiffness matrix K from the mesh alone
//
// [ ] calculate the term from the Neumann boundary conditions
// [ ] take the Dirichlet boundary conditions into account by doctoring K and f
//
// [ ] calculate Kc according to any convection boundary conditions
// [ ] calculate the term from the convection boundary conditions
//
// [X] don't bother with the load vector for now but write a function for it
// [X] solve linear system with Eigen

void add_load_vector(Matrix& f) {
    assert(f.num_columns == 1);
    // for (size_t i = 0; i < f.num_rows; ++i) {
    //     f.at(i, 0) += integral over domain of (Q Ni(x, y))
    // }
}

float length(Node v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}

float cross(Node a, Node b) {
    return a.x * b.y - a.y * b.x;
}

Node grad_Ni(size_t i, const Mesh& mesh, const Element& e) {

    Node nodes[3];
    nodes[0] = mesh.nodes[i];

    size_t index = 1;
    for (size_t v : e) {
        if (v != i) {
            nodes[index] = mesh.nodes[v];
            ++index;
        }
    }

    Node bc = {nodes[2].x - nodes[1].x, nodes[2].y - nodes[1].y};

    float len_bc = length(bc);
    bc.x /= len_bc;
    bc.y /= len_bc;

    Node ba = {nodes[0].x - nodes[1].x, nodes[0].y - nodes[1].y};

    float dot = ba.x * bc.x + ba.y * bc.y;

    Node perp_point { nodes[1].x + bc.x * dot, nodes[1].y + bc.y * dot };
    Node perp_to_ni { perp_point.x - nodes[0].x, perp_point.y - nodes[0].y };

    Node gradient = {-perp_to_ni.x, -perp_to_ni.y};
    float len_gradient = length(gradient);

    gradient.x /= len_gradient;
    gradient.y /= len_gradient;
    gradient.x /= len_gradient;
    gradient.y /= len_gradient;

    return gradient;
}

float area(Node v0, Node v1, Node v2) {

    Node ab{ v1.x - v0.x, v1.y - v0.y};
    Node ac{ v2.x - v0.x, v2.y - v0.y};

    return 0.5 * std::abs(cross(ab, ac));
}

Matrix stiffness_matrix_linear_elements(const Mesh& mesh) {

    // need integrate ( B^T * B )
    // B is 2 x n where n is num nodes
    // outcome is n x n matrix
    // each entry K(i, j) is  dNi/dx * dNj/dx + dNi/dy * dNj/dy
    //
    // these guys are piecewise constant so ought to be easy
    // if Ni and Nj don't share a triangle, 0
    // go over every triangle involving Ni and Nj (at most 2?)
    // figure out the values and multiply by triangle area

    size_t n = mesh.nodes.size();
    Matrix K{n, n, std::vector<float>(n * n, 0.f)};

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {

            // determine K(i, j) as sum of contributions from elements

            for (const Element& e : mesh.elements) {

                bool has_i = false;
                bool has_j = false;

                for (size_t v : e) {
                    if (v == i) has_i = true;
                    if (v == j) has_j = true;
                }

                if (!has_i || !has_j) continue;

                // determine the contribution and add!

                Node dNi_d = grad_Ni(i, mesh, e);
                Node dNj_d = grad_Ni(j, mesh, e);

                float value_on_element = dNi_d.x * dNj_d.x + dNi_d.y * dNj_d.y;

                float element_area = area(mesh.nodes[e[0]],
                                          mesh.nodes[e[1]],
                                          mesh.nodes[e[2]]);

                K.at(i, j) += value_on_element * element_area;
            }
        }
    }

    return K;
}

void add_neumann_boundary_conditions(const Mesh& mesh,
                                     const BoundaryConditions& bcs,
                                     Matrix *force_vector) {

    assert(force_vector->num_rows == mesh.nodes.size());

    // -integral (over bits of neumann boundary) ( Ni q )

    for (const BoundaryCondition& bc : bcs) {

        const NeumannBoundaryCondition *nbc
            = std::get_if<NeumannBoundaryCondition>(&(bc.condition));

        if (!nbc) continue;

        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            if (i == bc.boundary[0] || i == bc.boundary[1]) {

                Node boundary_start = mesh.nodes[bc.boundary[0]];
                Node boundary_end = mesh.nodes[bc.boundary[1]];

                Node boundary_segment { boundary_end.x - boundary_start.x,
                                        boundary_end.y - boundary_start.y };

                float segment_length = length(boundary_segment);
                float contribution = -0.5f * segment_length * nbc->flux;

                force_vector->at(i, 0) += contribution;
            }
        }
    }
}

std::vector<float>
solve_heat_equation(const Mesh& mesh,
                    const BoundaryConditions& boundary_conditions) {

    size_t n = mesh.nodes.size();

    std::vector<float> result(n, 1.f);

    Matrix K = stiffness_matrix_linear_elements(mesh);
    log_matrix(std::cout, K);

    Matrix f{n, 1, std::vector<float>(n, 0.f)};
    add_neumann_boundary_conditions(mesh, boundary_conditions, &f);
    log_matrix(std::cout, f);

    //check_solve_linear_system();

    return result;
}

} // namespace tomfem
