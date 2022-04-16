#pragma once
#include <vector>
#include <iostream>

enum class VectorType { Column, Row };

template <typename T>
class Vector {
public:
    Vector(const std::vector<T>& values, const VectorType type) {
        _values = values;
        _type = type;
        _n = values.size();
    }

    size_t size() const {
        return _n;
    }

    void print(size_t print_count) const {
        std::cout << "[";
        if (2 * print_count >= _n) {
            for (size_t i = 0; i < _n; ++i) {
                std::cout << _values[i];
                if (i + 1 != _n) {
                    std::cout << ", ";
                }
            }
        } else {
            for (size_t i = 0; i < print_count; ++i) {
                std::cout << _values[i] << ", ";
            }
            std::cout << ". . ., ";
            for (size_t i = _n - print_count - 1; i < _n; ++i) {
                std::cout << _values[i];
                if (i + 1 != _n) {
                    std::cout << ", ";
                }
            }
        }
        std::cout << "], type: " << (_type == VectorType::Column ? "Column" : "Row") 
                  << ", size: " << _n << std::endl;
    }

private:
    std::vector<T> _values;
    VectorType _type;
    size_t _n;
};
