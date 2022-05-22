#pragma once
#include <vector>
#include <iostream>

template <typename T>
class DenseVector {
public:
    DenseVector(const std::vector<T>& x) : _n(x.size()), _vector_data(x) {}

    T get(size_t i) const {
        return _vector_data[i];
    }

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

private:
    size_t _n;
    std::vector<T> _vector_data;
};
