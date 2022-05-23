#pragma once
#include <iostream>
#include <cmath>
#include "omp.h"
#include "matrix.hpp"
#include "vector.hpp"
#include "cramer_rule.hpp"
#include "operations.hpp"

enum { TRIANGULAR_MATRIX_PRINT_COUNT = 5};

template <typename matrix_T, typename vector_T, typename result_T>
class LinearSystem {
public:
    LinearSystem(const SparseMatrix<matrix_T>& A, const DenseVector<vector_T>& b): _A(A), _b(b), _n(A.size()) {
        if (A.size() != b.size()) {
            throw "A and b does not have same dimensions";
        }
    }

    DenseVector<result_T> solve_cg_method(int threads_num, bool mpi_machine_used, 
                                          double *first_stage_elapsed_time, 
                                          double *second_stage_elapsed_time,
                                          bool calculate_bandwidth,
                                          long long max_iter,
                                          double epsilon) {
        if (!mpi_machine_used) {
            std::cout << "turning off omp set dynamic" << std::endl;
            omp_set_dynamic(0);
            omp_set_num_threads(threads_num);
        }
        #pragma omp parallel
        {
            #pragma omp single
            std::cout << "Solving linear system with CG Method using " << omp_get_num_threads() << " threads" << std::endl;
        }

        double elapsed_time;
        double bandwidth;
        DenseVector<result_T> x_previous(std::vector<result_T>(_n, 0));
        auto Ax_0 = SpMV<result_T>(_A, x_previous, &elapsed_time, calculate_bandwidth, &bandwidth);
        auto r_previous = axpby<result_T>(1, _b, -1, Ax_0, &elapsed_time, calculate_bandwidth, &bandwidth);
        DenseVector<result_T> p_previous;
        result_T ro_previous;
        bool convergence = false;
        int iteration = 1;
        do {
            DenseVector<result_T> z_k = r_previous;
            auto ro_k = dot_product<result_T>(r_previous, z_k, &elapsed_time, calculate_bandwidth, &bandwidth);
            DenseVector<result_T> p_k;
            if (iteration == 1) {
                p_k = z_k;
            } else {
                auto beta_k = ro_k / ro_previous;
                p_k = axpby<result_T>(1, z_k, beta_k, p_previous, &elapsed_time, calculate_bandwidth, &bandwidth);
            }
            auto q_k = SpMV<result_T>(_A, p_k, &elapsed_time, calculate_bandwidth, &bandwidth);
            _A.print(5);
            p_k.print(5);
            q_k.print(5);
            auto alpha_k = ro_k / dot_product<result_T>(p_k, q_k, &elapsed_time, calculate_bandwidth, &bandwidth);
            auto x_k = axpby<result_T>(1, x_previous, alpha_k, p_k, &elapsed_time, calculate_bandwidth, &bandwidth);
            auto r_k = axpby<result_T>(1, r_previous, -alpha_k, q_k, &elapsed_time, calculate_bandwidth, &bandwidth);
            if (ro_k < epsilon || iteration >= max_iter) {
                convergence = true;
            } else {
                ++iteration;
                p_previous = p_k;
                ro_previous = ro_k;
                r_previous = r_k;
            }
            x_previous = x_k;
        } while (convergence);
        
        return x_previous;
    }

    double calculate_residual(const DenseVector<result_T>& x) {
        std::vector<result_T> residual(_n, 0);
        #pragma omp parallel for
        for (size_t i = 0; i < _n; ++i) {
            for (size_t j = 0; j < _A.get_ellpack_m(); ++j) {
                residual[i] += _A.get_ellpack_val(i)[j] * x.get(_A.get_ellpack_col(i)[j]);
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

    double calculate_error(const DenseVector<result_T>& x) {
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
    const SparseMatrix<matrix_T>& _A;
    const DenseVector<vector_T>& _b;
    const size_t _n;
};
