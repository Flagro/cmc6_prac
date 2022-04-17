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
        std::unique_ptr<Matrix<matrix_T> > result(new FormulaMatrix<matrix_T>(_n, matrix_function_generators[_test_id - 1]));
        return result;
    }

    std::unique_ptr<Vector<vector_T> > get_b() {
        std::unique_ptr<Vector<vector_T> > result(new FormulaVector<vector_T>(_n, vector_function_generators[_test_id - 1]));
        return result;
    }

private:
    int _n;
    int _test_id;
};
