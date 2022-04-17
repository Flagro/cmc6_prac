#pragma once
#include "matrix.hpp"
#include "vector.hpp"
#include "formula_generators.hpp"
#include <cstdio>
#include <memory>

template <typename matrix_T, typename vector_T>
class InputGenerator {
public:
    InputGenerator(int n, int test_id) : _n(n), _test_id(test_id) {}
    Matrix<matrix_T> get_A() {
        return Matrix<matrix_T>(_formula_generators[_test_id]->generate_A());
    }

    Vector<vector_T> get_b() {
        return Vector<vector_T>(_formula_generators[_test_id]->generate_b(), VectorType::Column);
    }

private:
    int _n;
    int _test_id;
};
