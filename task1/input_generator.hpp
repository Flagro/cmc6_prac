#pragma once
#include "matrix.hpp"
#include "vector.hpp"
#include "formula_generators.hpp"
#include <cstdio>
#include <memory>

template <typename matrix_T, typename vector_T>
class InputGenerator {
public:
    InputGenerator(int argc, char* argv[]) {
        if (argc < 4) {
            throw "not enough pragram parameters";
        }
        _n = strtol(argv[1], nullptr, 10);
        _threads_num = strtol(argv[2], nullptr, 10);
        _test_id = strtol(argv[3], nullptr, 10) - 1;
        _formula_generators.push_back(std::unique_ptr<BaseFormulaGenerator<matrix_T, vector_T> >
                                      (new FormulaGenerator1<matrix_T, vector_T>(_n)));
    }
    Matrix<matrix_T> get_A() {
        return Matrix<matrix_T>(_formula_generators[_test_id]->generate_A());
    }

    Vector<vector_T> get_b() {
        return Vector<vector_T>(_formula_generators[_test_id]->generate_b(), VectorType::Column);
    }

private:
    int _n;
    int _threads_num;
    int _test_id;
    std::vector<std::unique_ptr<BaseFormulaGenerator<matrix_T, vector_T> > > _formula_generators;
};
