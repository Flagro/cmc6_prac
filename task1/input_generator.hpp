#pragma once
#include "matrix.hpp"
#include "vector.hpp"
#include "formula_generators.hpp"
#include <cstdio>
#include <memory>

template <typename matrix_T, typename vector_T>
class InputGenerator {
public:
    InputGenerator(int n, int test_id) : _n(n), _test_id(test_id) {
        if (_test_id < 1 || _test_id > (int) matrix_function_generators.size()) {
            throw "test id is invalid";
        }
    }
    std::unique_ptr<Matrix<matrix_T> > get_A() {
        std::vector<std::vector<matrix_T> > result_matrix(_n, std::vector<matrix_T>(_n, 0));
        for (int i = 0; i < _n; ++i) {
            for (int j = 0 ; j < _n; ++j) {
                result_matrix[i][j] = matrix_function_generators[_test_id - 1].get(i, j, _n);
            }
        }
        std::unique_ptr<Matrix<matrix_T> > result(new ArrayMatrix<matrix_T>(_n, result_matrix));
        return result;
    }

    std::unique_ptr<Vector<vector_T> > get_b() {
        std::vector<vector_T> result_vector(_n, 0);
        for (int i = 0; i < _n; ++i) {
            result_vector[i] = vector_function_generators[_test_id - 1].get(i, _n);
        }
        std::unique_ptr<Vector<vector_T> > result(new ArrayVector<vector_T>(_n, result_vector));
        return result;
    }

private:
    int _n;
    int _test_id;
};
