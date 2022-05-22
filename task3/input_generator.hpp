#pragma once
#include "matrix.hpp"
#include "vector.hpp"
#include "formula_generators.hpp"
#include <cstdio>
#include <memory>

template <typename matrix_T, typename vector_T>
class InputGenerator {
public:
    InputGenerator(int n_x, int n_y, int n_z, int test_id) : _n_x(n_x), _n_y(n_y), _n_z(n_z), _test_id(test_id) {
        if (_test_id < 1 || _test_id > (int) cube_matrix_function_generators.size()) {
            throw "test id is invalid";
        }
    }
    
    SparseMatrix<matrix_T> get_A() {
        return SparseMatrix<matrix_T>(cube_matrix_function_generators[_test_id - 1].generate(_n_x, _n_y, _n_z));
    }

    DenseVector<vector_T> get_b() {
        return DenseVector<vector_T>(vector_function_generators[_test_id - 1].generate(_n_x * _n_y * _n_z));
    }

private:
    int _n_x;
    int _n_y;
    int _n_z;
    int _test_id;
};
