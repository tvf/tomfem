#include <vector>

namespace tomfem {

struct Matrix {

    size_t num_rows;
    size_t num_columns;

    std::vector<float> values;

    float& at(size_t i, size_t j) {
        return values[i * num_columns + j];
    }

    float at(size_t i, size_t j) const {
        return values[i * num_columns + j];
    }
};

// void transpose(Matrix *m);

// Return an mx1 matrix x s.t. Ax = b.
Matrix solve_linear_system(const Matrix& Anxm, const Matrix& bnx1);

}  // namespace tomfem
