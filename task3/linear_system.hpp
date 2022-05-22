#pragma once
#include <iostream>
#include <cmath>
#include "omp.h"
#include "matrix.hpp"
#include "vector.hpp"
#include "cramer_rule.hpp"

enum { TRIANGULAR_MATRIX_PRINT_COUNT = 5};

template <typename matrix_T, typename vector_T, typename result_T>
class LinearSystem {
public:
    LinearSystem(const Matrix<matrix_T>& A, const Vector<vector_T>& b): _A(A), _b(b), _n(A.size()) {
        if (A.size() != b.size()) {
            throw "A and b does not have same dimensions";
        }
    }

    std::unique_ptr<Vector<result_T> > solve_cg_method(int threads_num, bool mpi_machine_used, 
                                                               double *first_stage_elapsed_time, 
                                                               double *second_stage_elapsed_time) {
        if (!mpi_machine_used) {
            std::cout << "turning off omp set dynamic" << std::endl;
            omp_set_dynamic(0);
            omp_set_num_threads(threads_num);
        }
        #pragma omp parallel
        {
            #pragma omp single
            std::cout << "Solving linear system with Reflection Method using " << omp_get_num_threads() << " threads" << std::endl;
        }
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
        
        *first_stage_elapsed_time = first_stage_elapsed_time_second - first_stage_elapsed_time_first;
        *second_stage_elapsed_time = second_stage_elapsed_time_second - second_stage_elapsed_time_first;
        std::unique_ptr<Vector<result_T> > result_uptr(new ArrayVector<result_T>(_n, result_array));
        return result_uptr;
    }

    double calculate_residual(const Vector<result_T>& x) {
        std::vector<result_T> residual(_n, 0);
        #pragma omp parallel for
        for (size_t i = 0; i < _n; ++i) {
            for (size_t j = 0; j < _n; ++j) {
                residual[i] += _A.get(i, j) * x.get(j);
            }
            residual[i] -= _b.get(i);
        }
        result_T residual_norm = 0;
        #pragma omp parallel for reduction(+:residual_norm) 
        for (size_t i = 0; i < _n; ++i) {
            residual_norm += residual[i] * residual[i];
        }
        return sqrt(residual_norm);
    }

    double calculate_error(const Vector<result_T>& x) {
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
        std::vector<std::pair<double, double> > almost_exact_result;
        try {
            almost_exact_result = solveCramer(result_matrix, result_vector);
        } catch(std::runtime_error &e) {
            std::cout << "Error occured, Not calculating the error" << std::endl;
            return 0;
        }
        std::cout << "Found an answer with cramer's rule:" << std::endl;
        for (size_t i = 0; i < _n; ++i) {
            std::cout << almost_exact_result[i].first / almost_exact_result[i].second << " ";
        }
        std::cout << std::endl;
        double error = 0;
        // comparing x_i and (a_1i / a_2i) means comparing (x_i * a_2i and a_1i)
        for (size_t i = 0; i < _n; ++i) {
            double cur_difference = (x.get(i) * almost_exact_result[i].second - almost_exact_result[i].first) / almost_exact_result[i].second;
            error += cur_difference * cur_difference;
        }
        return sqrt(error);
    }

private:
    const Matrix<matrix_T>& _A;
    const Vector<vector_T>& _b;
    const size_t _n;
};
