#pragma once
#include "omp.h"
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

    std::unique_ptr<Vector<result_T> > solve_reflection_method(double *first_stage_elapsed_time, 
                                                               double *second_stage_elapsed_time) {
        auto first_stage_elapsed_time_first = omp_get_wtime();

        auto first_stage_elapsed_time_second = omp_get_wtime();

        auto second_stage_elapsed_time_first = omp_get_wtime();
        
        auto second_stage_elapsed_time_second = omp_get_wtime();
        *first_stage_elapsed_time = first_stage_elapsed_time_second - first_stage_elapsed_time_first;
        *second_stage_elapsed_time = second_stage_elapsed_time_second - second_stage_elapsed_time_first;
        std::vector<result_T> result_array(_n, 0);
        std::unique_ptr<Vector<result_T> > result_uptr(new ArrayVector<result_T>(_n, result_array));
        return result_uptr;
    }

    double calculate_residual(const Vector<result_T>& x) {
        return x.size();
    }

    double calculate_error(const Vector<result_T>& x) {
        return x.size();
    }

private:
    const Matrix<matrix_T>& _A;
    const Vector<vector_T>& _b;
    const size_t _n;
};
