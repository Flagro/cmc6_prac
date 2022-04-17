#pragma once
#include <vector>
#include <variant>
#include <memory>
#include "formula_generators.hpp"

template <typename T>
class Matrix {
public:
    Matrix(const BaseMatrixGenerator* matrix_generator_ptr) {
        _contain_type = MatrixContainType::FormulaGenerator;
        _n = matrix_generator_ptr->n;
        if (!_n) {
            throw "passed empty matrix";
        }
        _values = values;
    }

    size_t size() const {
        return _n;
    }

    void print(size_t print_count) const {
        if (2 * print_count >= _n) {
            for (size_t i = 0; i < _n; ++i) {
                for (size_t j = 0; j < _n; ++j) {
                    std::cout << _values[i][j];
                    if (j + 1 != _n) {
                        std::cout << ", ";
                    }
                }
                std::cout << std::endl;
            }
        } else {
            for (size_t i = 0; i < print_count; ++i) {
                for (size_t j = 0; j < print_count; ++j) {
                    std::cout << _values[i][j] << ", ";
                }
                std::cout << ". . ., ";
                for (size_t j = _n - print_count - 1; j < _n; ++j) {
                    std::cout << _values[i][j];
                    if (j + 1 != _n) {
                        std::cout << ", ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << ". . ., ";
            for (size_t i = _n - print_count - 1; i < _n; ++i) {
                for (size_t j = 0; j < print_count; ++j) {
                    std::cout << _values[i][j] << ", ";
                }
                std::cout << ". . ., ";
                for (size_t j = _n - print_count - 1; j < _n; ++j) {
                    std::cout << _values[i][j];
                    if (j + 1 != _n) {
                        std::cout << ", ";
                    }
                }
                std::cout << std::endl;
            }
        }
        std::cout << "size: " << _n << std::endl;
    }

private:
    std::vector<std::vector<T> > _values;
    size_t _n;
    MatrixContainType _contain_type;
};

template <typename T>
class FormulaMatrix {
public:
    Matrix(const BaseMatrixGenerator* matrix_generator_ptr) {
        _contain_type = MatrixContainType::FormulaGenerator;
        _n = matrix_generator_ptr->n;
        if (!_n) {
            throw "passed empty matrix";
        }
        _values = values;
    }

private:
        
};

template <typename T>
class ArrayMatrix {
public:
    Matrix(std::vector<std::vector<T> > values) {
        _contain_type = MatrixContainType::VectorContainer;
        if (!values.size()) {
            throw "passed empty matrix";
        }
        if (values.size() != values[0].size()) {
            throw "passed non square matrix";
        }
        _values = values;
        _n = values.size();
    }

private:
    std::vector<std::vector<T> > _values;
};
