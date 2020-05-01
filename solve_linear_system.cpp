#include "solve_linear_system.h"

#include <Eigen/Dense>

namespace tomfem {

// void transpose(Matrix *m);

// Return an mx1 matrix x s.t. Ax = b.
Matrix solve_linear_system(const Matrix& Anxm, const Matrix& bnx1) {

    Eigen::MatrixXf A(Anxm.num_rows, Anxm.num_columns);
    for (size_t i = 0; i < Anxm.num_rows; ++i) {
        for (size_t j = 0; j < Anxm.num_columns; ++j) {
            A(i, j) = Anxm.at(i, j);
        }
    }

    Eigen::VectorXf b(bnx1.num_rows);
    for (size_t i = 0; i < bnx1.num_rows; ++i) {
        b(i) = bnx1.at(i, 0);
    }

    Eigen::VectorXf x = A.colPivHouseholderQr().solve(b);

    size_t eigensize = x.size();

    assert(Anxm.num_columns == eigensize);

    Matrix result{eigensize, 1, std::vector<float>(eigensize)};

    for (size_t i = 0; i < eigensize; ++i) {
        result.at(i, 0) = x(i);
    }

    return result;
}

}  // namespace tomfem
