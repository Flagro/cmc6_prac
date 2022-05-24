#pragma once
#include <iostream>
#include <cmath>
#include "omp.h"
#include "matrix.hpp"
#include "vector.hpp"
#include "cramer_rule.hpp"
#include "operations.hpp"

enum { TRIANGULAR_MATRIX_PRINT_COUNT = 5};

template <typename result_T>
struct CGSolverResult {
    DenseVector<result_T> x;
    double dot_product_elapsed_time, spmv_elapsed_time, axpby_elapsed_time;
    double total_cg_time;
    double dot_product_bandwidth, spmv_bandwidth, axpby_bandwidth;
    int iterations_taken;
    double cg_residual;

    void print(double double_print_precision) const {
        std::cout << std::setprecision(double_print_precision) << "Calculated vector x:" << std::endl;
        x.print(double_print_precision);
        std::cout << std::endl;
        std::cout << std::setprecision(double_print_precision) 
                  << "Total CG solver elapsed time: " << total_cg_time << " seconds" << std::endl
                  << "CG solver iterations taken: " << iterations_taken << std::endl
                  << "CG residual: " << cg_residual << std::endl
                  << "dot product total time: " << dot_product_elapsed_time << " seconds" << std::endl
                  << "spmv total time: " << spmv_elapsed_time << " seconds" << std::endl
                  << "axpby total time: " << axpby_elapsed_time << " seconds" << std::endl
                  << "dot product bandwidth: " << dot_product_bandwidth << " Gbps" << std::endl
                  << "spmv bandwidth: " << spmv_bandwidth << " Gbps" << std::endl
                  << "axpby bandwidth: " << axpby_bandwidth << " Gbps" << std::endl;
    }
};

template <typename matrix_T, typename vector_T, typename result_T>
class LinearSystem {
public:
    LinearSystem(const SparseMatrix<matrix_T>& A, const DenseVector<vector_T>& b): _A(A), _b(b), _n(A.size()) {
        if (A.size() != b.size()) {
            throw "A and b does not have same dimensions";
        }
    }

    CGSolverResult<result_T> solve_cg_method(int threads_num, bool mpi_machine_used,
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
        double timer_start = omp_get_wtime();
        double dot_product_elapsed_time = 0, spmv_elapsed_time = 0, axpby_elapsed_time = 0;
        double dot_product_bandwidth_bytes = 0, spmv_bandwidth_bytes = 0, axpby_bandwidth_bytes = 0;
        DenseVector<result_T> x_previous(std::vector<result_T>(_n, 0));
        auto Ax_0 = SpMV<result_T>(_A, x_previous, &spmv_elapsed_time, calculate_bandwidth, &spmv_bandwidth_bytes);
        auto r_previous = axpby<result_T>(1, _b, -1, Ax_0, &axpby_elapsed_time, calculate_bandwidth, &axpby_bandwidth_bytes);
        DenseVector<result_T> p_previous;
        result_T ro_previous;
        bool convergence = false;
        int iteration = 1;
        do {
            DenseVector<result_T> z_k = r_previous;
            auto ro_k = dot_product<result_T>(r_previous, z_k, &dot_product_elapsed_time, calculate_bandwidth, &dot_product_bandwidth_bytes);
            DenseVector<result_T> p_k;
            if (iteration == 1) {
                p_k = z_k;
            } else {
                auto beta_k = ro_k / ro_previous;
                p_k = axpby<result_T>(1, z_k, beta_k, p_previous, &axpby_elapsed_time, calculate_bandwidth, &axpby_bandwidth_bytes);
            }
            auto q_k = SpMV<result_T>(_A, p_k, &spmv_elapsed_time, calculate_bandwidth, &spmv_bandwidth_bytes);
            auto alpha_k = ro_k / dot_product<result_T>(p_k, q_k, &dot_product_elapsed_time, calculate_bandwidth, &dot_product_bandwidth_bytes);
            auto x_k = axpby<result_T>(1, x_previous, alpha_k, p_k, &axpby_elapsed_time, calculate_bandwidth, &axpby_bandwidth_bytes);
            auto r_k = axpby<result_T>(1, r_previous, -alpha_k, q_k, &axpby_elapsed_time, calculate_bandwidth, &axpby_bandwidth_bytes);
            if (ro_k < epsilon || iteration >= max_iter) {
                convergence = true;
            } else {
                ++iteration;
                p_previous = p_k;
                ro_previous = ro_k;
                r_previous = r_k;
            }
            x_previous = x_k;
        } while (!convergence);
        double timer_end = omp_get_wtime();
        CGSolverResult<result_T> result;
        result.x = x_previous;
        result.dot_product_elapsed_time = dot_product_elapsed_time;
        result.spmv_elapsed_time = spmv_elapsed_time;
        result.axpby_elapsed_time = axpby_elapsed_time;
        result.dot_product_bandwidth = dot_product_bandwidth_bytes / dot_product_elapsed_time / 1048576;
        result.spmv_bandwidth = spmv_bandwidth_bytes / spmv_elapsed_time / 1048576;
        result.axpby_bandwidth = axpby_bandwidth_bytes / axpby_elapsed_time / 1048576;
        result.iterations_taken = iteration;
        result.cg_residual = ro_previous;
        result.total_cg_time = timer_end - timer_start;
        return result;
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
