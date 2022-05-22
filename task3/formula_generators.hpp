#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include "ellpack_format.hpp"

template <typename matrix_T>
struct CubeMatrixGenerator {
public:
    CubeMatrixGenerator(const std::function<EllpackMatrix<matrix_T>(size_t, size_t, size_t)>& generate_function): generate(generate_function) {}
    const std::function<EllpackMatrix<matrix_T>(size_t, size_t, size_t)> generate;
};

template <typename vector_T>
struct VectorGenerator {
public:
    VectorGenerator(const std::function<std::vector<vector_T>(size_t)>& generate_function): generate(generate_function) {}
    const std::function<std::vector<vector_T>(size_t)> generate;
};

std::vector<CubeMatrixGenerator<double> > cube_matrix_function_generators = {
    {[](size_t n_x, size_t n_y, size_t n_z) {
        static const size_t M = 7;
        static const double PI = 3.14;
        size_t n = n_x * n_y * n_z;
        EllpackMatrix<double> result;
        result.n = n;
        result.m = M;
        result.ellpack_col = std::vector<std::vector<size_t> >(n, std::vector<size_t>(M));
        result.ellpack_val = std::vector<std::vector<double> >(n, std::vector<double>(M));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < M; ++j) {
                result.ellpack_col[i][j];
                result.ellpack_val[i][j];
            }
        }
        return result;
    }}
};

std::vector<VectorGenerator<double> > vector_function_generators = {
    {[](size_t n) {
        std::vector<double> result(n);
        for (size_t i = 0; i < n; ++i) {
            result[i] = cos(i);
        }
        return result;
    }}
};
