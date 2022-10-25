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
        static const auto non_diagonal_value = [](size_t i, size_t j) { return cos(i * j + PI); };

        size_t n = n_x * n_y * n_z;

        EllpackMatrix<double> result;
        result.n = n;
        result.m = M;
        result.ellpack_col = std::vector<std::vector<size_t> >(n, std::vector<size_t>(M));
        result.ellpack_val = std::vector<std::vector<double> >(n, std::vector<double>(M, 0.0));
        for (size_t i = 0; i < n_x; ++i) {
            for (size_t j = 0; j < n_y; ++j) {
                for (size_t k = 0; k < n_z; ++k) {
                    size_t cur_row = k * n_x * n_y + j * n_x + i;
                    std::vector<std::pair<size_t, bool> > cur_cols_ids = 
                        {{cur_row - n_x * n_y, k > 0}, 
                         {cur_row - n_x, j > 0}, 
                         {cur_row - 1, i > 0}, 
                         {cur_row, true}, 
                         {cur_row + 1, i < n_x - 1}, 
                         {cur_row + n_x, j < n_y - 1}, 
                         {cur_row + n_x * n_y, k < n_z - 1}};

                    size_t cur_col_offset = 0;
                    size_t diagonal_el_offset = 0;
                    double non_diagonal_norm_sum = 0;
                    for (const auto &el : cur_cols_ids) {
                        const auto cur_col = el.first;
                        const auto to_use = el.second;
                        if (to_use) {
                            result.ellpack_col[cur_row][cur_col_offset] = cur_col;
                            if (cur_row == cur_col) {
                                diagonal_el_offset = cur_col_offset;
                            } else {
                                result.ellpack_val[cur_row][cur_col_offset] = non_diagonal_value(cur_row, cur_col);
                                non_diagonal_norm_sum += fabs(result.ellpack_val[cur_row][cur_col_offset]);
                            }
                            ++cur_col_offset;
                        }
                    }
                    result.ellpack_val[cur_row][diagonal_el_offset] = 1.5 * non_diagonal_norm_sum;
                    for (size_t zero_val_id = cur_col_offset; zero_val_id < M; ++zero_val_id) {
                        result.ellpack_col[cur_row][zero_val_id] = n + 1 + zero_val_id - cur_col_offset;
                    }
                }
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
    }}, {[](size_t n) {
        std::vector<double> result(n);
        for (size_t i = 0; i < n; ++i) {
            result[i] = cos(i * i);
        }
        return result;
    }}, {[](size_t n) {
        std::vector<double> result(n);
        for (size_t i = 0; i < n; ++i) {
            result[i] = sin(i * i);
        }
        return result;
    }}
};
