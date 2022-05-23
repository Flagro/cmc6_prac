#pragma once
#include <vector>
#include <memory>
#include <algorithm>
#include "formula_generators.hpp"
#include "ellpack_format.hpp"

template <typename T>
class SparseMatrix {
public:
    SparseMatrix() = default;
    SparseMatrix(const EllpackMatrix<T>& ellpack_matrix) : _n(ellpack_matrix.n), _ellpack_data(ellpack_matrix) {}

    const T & get(size_t i, size_t j) const {
        static const T zero_val = 0;
        auto lower_bound_it = std::lower_bound(_ellpack_data.ellpack_col[i].begin(), _ellpack_data.ellpack_col[i].end(), j);
        if (lower_bound_it != _ellpack_data.ellpack_col[i].end() && *lower_bound_it == j) {
            return _ellpack_data.ellpack_val[i][lower_bound_it - _ellpack_data.ellpack_col[i].begin()];
        }
        return zero_val;
    }

    size_t get_ellpack_m() const {
        return _ellpack_data.m;
    }

    const std::vector<size_t> & get_ellpack_col(size_t i) const {
        return _ellpack_data.ellpack_col[i];
    }

    const std::vector<T> & get_ellpack_val(size_t i) const {
        return _ellpack_data.ellpack_val[i];
    }

    size_t size() const {
        return _n;
    }

    void print(size_t print_count) const {
        if (2 * print_count >= _n) {
            for (size_t i = 0; i < _n; ++i) {
                for (size_t j = 0; j < _n; ++j) {
                    std::cout << get(i, j);
                    if (j + 1 != _n) {
                        std::cout << " ";
                    }
                }
                std::cout << std::endl;
            }
        } else {
            for (size_t i = 0; i < print_count; ++i) {
                for (size_t j = 0; j < print_count; ++j) {
                    std::cout << get(i, j) << " ";
                }
                std::cout << ". . . ";
                for (size_t j = _n - print_count - 1; j < _n; ++j) {
                    std::cout << get(i, j);
                    if (j + 1 != _n) {
                        std::cout << " ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << ". . . " << std::endl;
            for (size_t i = _n - print_count - 1; i < _n; ++i) {
                for (size_t j = 0; j < print_count; ++j) {
                    std::cout << get(i, j) << " ";
                }
                std::cout << ". . . ";
                for (size_t j = _n - print_count - 1; j < _n; ++j) {
                    std::cout << get(i, j);
                    if (j + 1 != _n) {
                        std::cout << " ";
                    }
                }
                std::cout << std::endl;
            }
        }
        std::cout << "size: " << _n << std::endl;
    }

private:
    size_t _n;
    EllpackMatrix<T> _ellpack_data;
};
