#pragma once
#include <vector>
#include <iostream>

enum class VectorType { Column, Row };

template <typename T>
class Vector {
public:
    Vector(size_t n) : _n(n) {}

    virtual ~Vector() {}

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
        std::cout << "], size: " << _n << std::endl;
    }

protected:
    size_t _n;
};

template <typename T>
class FormulaVector : public Vector<T> {
public:
    FormulaVector(size_t n,  const VectorGenerator<T>& vector_generator) : Vector<T>(n), 
                _vector_generator(vector_generator) {
        if (!n) {
            throw "passed empty matrix";
        }
    }

    T get(size_t i) const {
        return _vector_generator.get(i, this->_n);
    }

private:
    const VectorGenerator<T>& _vector_generator;
};

template <typename T>
class ArrayVector : public Vector<T> {
public:
    ArrayVector(size_t n, const std::vector<T>& values) : Vector<T>(n)  {
        if (!values.size()) {
            throw "passed empty vector";
        }
        _values = values;
    }

    T get(size_t i) const {
        return _values[i];
    }

    void set(size_t i, T new_value) {
        _values[i] = new_value;
    }

private:
    std::vector<T> _values;
};
