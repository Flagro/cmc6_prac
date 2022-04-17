#pragma once
#include <vector>
#include <iostream>

enum class VectorType { Column, Row };

template <typename T>
class Vector {
public:
    Vector(size_t n) : _n(n) {}

    virtual T get(size_t i) const = 0;

    size_t size() const {
        return _n;
    }

    void print(size_t print_count) const {
        std::cout << "[";
        if (2 * print_count >= _n) {
            for (size_t i = 0; i < _n; ++i) {
                std::cout << get(i);
                if (i + 1 != _n) {
                    std::cout << ", ";
                }
            }
        } else {
            for (size_t i = 0; i < print_count; ++i) {
                std::cout << get(i) << ", ";
            }
            std::cout << ". . ., ";
            for (size_t i = _n - print_count - 1; i < _n; ++i) {
                std::cout << get(i);
                if (i + 1 != _n) {
                    std::cout << ", ";
                }
            }
        }
        std::cout << "], type: " << (_type == VectorType::Column ? "Column" : "Row") 
                  << ", size: " << _n << std::endl;
    }

protected:
    size_t _n;
};

template <typename T>
class FormulaVector : private Vector<T> {
public:
    FormulaVector(size_t n,  const VectorGenerator& vector_generator) : Vector(n) {
        if (!_n) {
            throw "passed empty matrix";
        }
        _vector_generator = vector_generator;
    }

    T get(size_t i) {
        return _vector_generator(i, _n);
    }

private:
    const MatrixGenerator& _vector_generator;
};

template <typename T>
class ArrayVector : private Vector<T> {
public:
    ArrayVector(size_t n, const std::vector<T>& values) : Vector(n)  {
        if (!values.size()) {
            throw "passed empty vector";
        }
        _values = values;
        _n = n;
    }

    T get(size_t i) {
        return _values[i];
    }

    T & get(size_t i) {
        return _values[i];
    }

private:
    std::vector<T> _values;
};
