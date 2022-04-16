#pragma once
#include "matrix.hpp"
#include "vector.hpp"

template <typename matrix_T, typename vector_T, typename result_T>
class LinearSystem {
public:
    LinearSystem(const Matrix<matrix_T>& A, const Vector<vector_T>& b): _A(A), _b(b), _n(A.size()) {
        if (A.size() != b.size()) {
            throw "A and b does not have same dimensions";
        }
    }

    Vector<result_T> solve_reflection_method(double *first_stage_elapsed_time, 
                                             double *second_stage_elapsed_time) {
        *first_stage_elapsed_time = 0;
        *second_stage_elapsed_time = 0;
        std::vector<result_T> result(_n, 0);
        return Vector(result, VectorType::Column);
    }

    double calculate_residual(const Vector<result_T>& x) {
        return 0;
    }

    double calculate_error(const Vector<result_T>& x) {
        return 0;
    }

private:
    const Matrix<matrix_T>& _A;
    const Vector<vector_T>& _b;
    const size_t _n;
};
