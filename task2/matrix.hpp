#pragma once
#include <vector>
#include <memory>
#include "formula_generators.hpp"

template <typename T>
class Matrix {
public:
    Matrix(size_t n) : _n(n) {}
    
    virtual ~Matrix() {}

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

protected:
    size_t _n;
};

template <typename T>
class FormulaMatrix : public Matrix<T> {
public:
    FormulaMatrix(size_t n,  const MatrixGenerator<T>& matrix_generator) : Matrix<T>(n), 
            _matrix_generator(matrix_generator) {
        if (!n) {
            throw "passed empty matrix";
        }
    }

    T get(size_t i, size_t j) const {
        return _matrix_generator.get(i, j, this->_n);
    }

private:
    const MatrixGenerator<T>& _matrix_generator;
};

template <typename T>
class ArrayMatrix : public Matrix<T> {
public:
    ArrayMatrix(size_t n, const std::vector<std::vector<T> >& values) : Matrix<T>(n) {
        if (!values.size()) {
            throw "passed empty matrix";
        }
        if (values.size() != values[0].size()) {
            throw "passed non square matrix";
        }
        _values = values;
    }

    T get(size_t i, size_t j) const {
        return _values[i][j];
    }

    void set (size_t i, size_t j, T new_value) {
        _values[i][j] = new_value;
    }

private:
    std::vector<std::vector<T> > _values;
};
