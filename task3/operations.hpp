#pragma once
#include "omp.h"
#include "vector.hpp"
#include "matrix.hpp"

template<typename result_T, typename T1, typename T2>
result_T dot_product(const DenseVector<T1> &first_vector, const DenseVector<T2> &second_vector, double *elapsed_time) {
    if (first_vector.size() != second_vector.size()) {
        throw "different vector sizes for dot prodct";
    }
    result_T result = 0;
    #pragma omp parallel for
    for (size_t i = 0; i < first_vector.size(); ++i) {
        result += first_vector.get(i) * second_vector.get(i);
    }
    return result;
}

template<typename result_T, typename T1, typename T2>
result_T dot_product_sequential(const DenseVector<T1> &first_vector, const DenseVector<T2> &second_vector, double *elapsed_time) {
    if (first_vector.size() != second_vector.size()) {
        throw "different vector sizes for dot prodct";
    }
    result_T result = 0;
    for (size_t i = 0; i < first_vector.size(); ++i) {
        result += first_vector.get(i) * second_vector.get(i);
    }
    return result;
}
