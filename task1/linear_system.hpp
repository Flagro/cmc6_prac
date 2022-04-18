#pragma once
#include <iostream>
#include <cmath>
#include "omp.h"
#include "matrix.hpp"
#include "vector.hpp"

enum { TRIANGULAR_MATRIX_PRINT_COUNT = 5};

template <typename matrix_T, typename vector_T, typename result_T>
class LinearSystem {
public:
    LinearSystem(const Matrix<matrix_T>& A, const Vector<vector_T>& b): _A(A), _b(b), _n(A.size()) {
        if (A.size() != b.size()) {
            throw "A and b does not have same dimensions";
        }
    }

    std::unique_ptr<Vector<result_T> > solve_reflection_method(int threads_num, 
                                                               double *first_stage_elapsed_time, 
                                                               double *second_stage_elapsed_time) {
        std::vector<std::vector<result_T> > result_matrix(_n, std::vector<result_T>(_n, 0));
        for (size_t i = 0; i < _n; ++i) {
            for (size_t j = 0 ; j < _n; ++j) {
                result_matrix[i][j] = _A.get(i, j);
            }
        }
        std::vector<result_T> result_vector(_n, 0);
        for (size_t i = 0; i < _n; ++i) {
            result_vector[i] = _b.get(i);
        }
        std::vector<result_T> x(_n, 0);
        auto first_stage_elapsed_time_first = omp_get_wtime();
        for (size_t i = 0; i < _n; ++i) {
            result_T s = 0;
            #pragma omp parallel for reduction(+:s) 
            for (size_t j = i + 1; j < _n; ++j) {
                s += result_matrix[j][i] * result_matrix[j][i];
            }
            result_T cur_a_norm = sqrt(result_matrix[i][i] * result_matrix[i][i] + s);
            result_T cur_x1 = result_matrix[i][i] - cur_a_norm;
            result_T cur_x_norm = sqrt(cur_x1 * cur_x1 + s);
            x[i] = cur_x1 / cur_x_norm;
            #pragma omp parallel for
            for (size_t j = i + 1; j < _n; ++j) {
                x[j] = result_matrix[j][i] / cur_x_norm;
            }
            // doing mxv (U(x)b) = (I-2xx^*)b = b - 2x(b, x)
            result_T alpha = 0;
            #pragma omp parallel for reduction(+:alpha) 
            for (size_t j = i; j < _n; ++j) {
                alpha += result_vector[j] * x[j];
            }
            #pragma omp parallel for
            for (size_t j = i; j < _n; ++j) {
                result_vector[j] -= 2 * x[j] * alpha;
            }
            // doing mxv (U(x)A_{*j}) = A_{*j} - 2x(A_{*j}, x) for each j-th column from A
            result_matrix[i][i] = cur_a_norm;
            #pragma omp parallel for
            for (size_t j = i + 1; j < _n; ++j) {
                result_matrix[j][i] = 0;
            }
            #pragma omp parallel for
            for (size_t j = i + 1; j < _n; ++j) {
                result_T alpha = 0;
                for (size_t k = i; k < _n; ++k) {
                    alpha += result_matrix[k][j] * x[j];
                }
                for (size_t k = i; k < _n; ++k) {
                    result_matrix[k][j] -= 2 * x[j] * alpha;
                }
            }
        }

        auto first_stage_elapsed_time_second = omp_get_wtime();
        std::cout << "Triangular form:" << std::endl;
        ArrayMatrix<result_T>(_n, result_matrix).print(TRIANGULAR_MATRIX_PRINT_COUNT);
        ArrayVector<result_T>(_n, result_vector).print(TRIANGULAR_MATRIX_PRINT_COUNT);
        std::cout << std::endl;

        std::vector<result_T> result_array(_n, 0);
        auto second_stage_elapsed_time_first = omp_get_wtime();

        for (int i = (int)_n - 1; i >= 0; --i) {
            result_array[i] = result_matrix[i][i] / result_vector[i];
            for (int j = i - 1; j >= 0; --j) {
                result_vector[j] -= result_matrix[j][i] * result_array[i];
            }
        }
        auto second_stage_elapsed_time_second = omp_get_wtime();

        *first_stage_elapsed_time = first_stage_elapsed_time_second - first_stage_elapsed_time_first;
        *second_stage_elapsed_time = second_stage_elapsed_time_second - second_stage_elapsed_time_first;
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
