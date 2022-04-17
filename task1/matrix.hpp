#pragma once
#include <vector>
#include <variant>
#include <memory>
#include "formula_generators.hpp"

template <typename T>
class Matrix {
public:
    Matrix(size_t n) : _n(n) {}

    virtual T get(size_t i, size_t j) const = 0;

    size_t size() const {
        return _n;
    }

    void print(size_t print_count) const {
        if (2 * print_count >= _n) {
            for (size_t i = 0; i < _n; ++i) {
                for (size_t j = 0; j < _n; ++j) {
                    std::cout << get(i, j);
                    if (j + 1 != _n) {
                        std::cout << ", ";
                    }
                }
                std::cout << std::endl;
            }
        } else {
            for (size_t i = 0; i < print_count; ++i) {
                for (size_t j = 0; j < print_count; ++j) {
                    std::cout << get(i, j) << ", ";
                }
                std::cout << ". . ., ";
                for (size_t j = _n - print_count - 1; j < _n; ++j) {
                    std::cout << get(i, j);
                    if (j + 1 != _n) {
                        std::cout << ", ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << ". . ., ";
            for (size_t i = _n - print_count - 1; i < _n; ++i) {
                for (size_t j = 0; j < print_count; ++j) {
                    std::cout << get(i, j) << ", ";
                }
                std::cout << ". . ., ";
                for (size_t j = _n - print_count - 1; j < _n; ++j) {
                    std::cout << get(i, j);
                    if (j + 1 != _n) {
                        std::cout << ", ";
                    }
                }
                std::cout << std::endl;
            }
        }
        std::cout << "size: " << _n << std::endl;
    }

protected:
    size_t _n;
};

template <typename T>
class FormulaMatrix : private Matrix<T> {
public:
    FormulaMatrix(size_t n,  const MatrixGenerator& matrix_generator) : Matrix(n) {
        if (!_n) {
            throw "passed empty matrix";
        }
        _matrix_generator = matrix_generator;
    }

    T get(size_t i, size_t j) {
        return _matrix_generator(i, j, _n);
    }

private:
    const MatrixGenerator& _matrix_generator;
};

template <typename T>
class ArrayMatrix : private Matrix<T> {
public:
    ArrayMatrix(size_t n, const std::vector<std::vector<T> >& values) {
        if (!values.size()) {
            throw "passed empty matrix";
        }
        if (values.size() != values[0].size()) {
            throw "passed non square matrix";
        }
        _values = values;
        _n = values.size();
    }

    T get(size_t i, size_t j) {
        return _values[i][j];
    }

    T & get(size_t i, size_t j) {
        return _values[i][j];
    }

private:
    std::vector<std::vector<T> > _values;
};
