#pragma once
#include "omp.h"
#include "vector.hpp"
#include "matrix.hpp"

template<typename result_T, typename T1, typename T2>
result_T dot_product(const DenseVector<T1> &first_vector, const DenseVector<T2> &second_vector, 
        double *elapsed_time, bool calculate_bandwidth, double *bandwidth) {
    double timer_start = omp_get_wtime();
    long long bandwidth_bytes;
    if (calculate_bandwidth) {
        bandwidth_bytes += sizeof(result_T);
        bandwidth_bytes += sizeof(size_t);
        #pragma omp parallel for reduction(+:bandwidth_bytes)
        for (size_t i = 0; i < first_vector.size(); ++i) {
            bandwidth_bytes += sizeof(result_T);
            bandwidth_bytes += sizeof(T1);
            bandwidth_bytes += sizeof(T2);
        }
    }
    if (first_vector.size() != second_vector.size()) {
        throw "different vector sizes for dot prodct";
    }
    result_T result = 0;
    #pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < first_vector.size(); ++i) {
        result += first_vector.get(i) * second_vector.get(i);
    }

    double timer_end = omp_get_wtime();
    *bandwidth = bandwidth_bytes / (timer_end - timer_start) / 1000;
    *elapsed_time = timer_end - timer_start;
    return result;
}

template<typename result_T, typename T1, typename T2>
DenseVector<result_T> SpMV(const SparseMatrix<T1> &sparse_matrix, const DenseVector<T2> &dense_vector, 
        double *elapsed_time, bool calculate_bandwidth, double *bandwidth) {
    double timer_start = omp_get_wtime();
    long long bandwidth_bytes;
    if (calculate_bandwidth) {
        bandwidth_bytes += sizeof(result_T) * sparse_matrix.size();
        bandwidth_bytes += sizeof(size_t);
        #pragma omp parallel for reduction(+:bandwidth_bytes)
        for (size_t i = 0; i < dense_vector.size(); ++i) {
            bandwidth_bytes += sizeof(size_t);
            bandwidth_bytes += sizeof(T1);
            bandwidth_bytes += sizeof(T2);
            for (size_t j = 0; j < sparse_matrix.get_ellpack_m(); ++j) {
                bandwidth_bytes += sizeof(T1);
                bandwidth_bytes += sizeof(size_t);
                size_t cur_col = sparse_matrix.get_ellpack_col(i)[j];
                if (cur_col > dense_vector.size()) {
                    bandwidth_bytes += sizeof(T1);
                    bandwidth_bytes += sizeof(T2);
                }
            }
        }
    }
    if (sparse_matrix.size() != dense_vector.size()) {
        throw "different sizes for SpMV";
    }
    std::vector<result_T> result(sparse_matrix.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < dense_vector.size(); ++i) {
        result_T cur_result = 0;
        for (size_t j = 0; j < sparse_matrix.get_ellpack_m(); ++j) {
            size_t cur_col = sparse_matrix.get_ellpack_col(i)[j];
            if (cur_col > dense_vector.size()) {
                cur_result += sparse_matrix.get_ellpack_val(i)[j] * dense_vector.get(cur_col);
            }
        }
        result[i] = cur_result;
    }

    double timer_end = omp_get_wtime();
    *bandwidth = bandwidth_bytes / (timer_end - timer_start) / 1000;
    *elapsed_time = timer_end - timer_start;
    return DenseVector<result_T>(result);
}

template<typename result_T, typename T1, typename T2, typename T3, typename T4>
DenseVector<result_T> axpby(const T1 &a, const DenseVector<T2> &first_vector, 
        const T3 &b, const DenseVector<T4> &second_vector, 
        double *elapsed_time, bool calculate_bandwidth, double *bandwidth) {
    double timer_start = omp_get_wtime();
    long long bandwidth_bytes;
    if (calculate_bandwidth) {
        bandwidth_bytes += sizeof(result_T) * first_vector.size();
        bandwidth_bytes += sizeof(size_t);
        #pragma omp parallel for reduction(+:bandwidth_bytes)
        for (size_t i = 0; i < first_vector.size(); ++i) {
            bandwidth_bytes += sizeof(size_t);
            bandwidth_bytes += sizeof(T1);
            bandwidth_bytes += sizeof(T2);
            bandwidth_bytes += sizeof(T3);
            bandwidth_bytes += sizeof(T4);
        }
    }
    if (first_vector.size() != second_vector.size()) {
        throw "different sizes for axpby";
    }
    
    std::vector<result_T> result(first_vector.size());

    #pragma omp parallel for
    for (size_t i = 0; i < first_vector.size(); ++i) {
        result[i] = a * first_vector.get(i) + b * second_vector.get(i);
    }

    double timer_end = omp_get_wtime();
    *bandwidth = bandwidth_bytes / (timer_end - timer_start) / 1000;
    *elapsed_time = timer_end - timer_start;
    return DenseVector<result_T>(result);
}
