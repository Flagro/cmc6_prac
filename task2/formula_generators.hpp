#pragma once
#include <vector>
#include <cstdlib>
#include <functional>

template <typename matrix_T>
struct MatrixGenerator {
public:
    MatrixGenerator(const std::function<matrix_T(size_t, size_t, size_t)>& get_function, 
                    const std::string& representation_string): get(get_function), 
                                                               function_representation(representation_string) {}
    const std::function<matrix_T(size_t, size_t, size_t)> get;
    const std::string function_representation;
};

template <typename vector_T>
struct VectorGenerator {
public:
    VectorGenerator(const std::function<vector_T(size_t, size_t)>& get_function, 
                    const std::string& representation_string): get(get_function), 
                                                               function_representation(representation_string) {}
    const std::function<vector_T(size_t, size_t)> get;
    const std::string function_representation;
};

std::vector<MatrixGenerator<double> > matrix_function_generators = {
    {[](size_t i, size_t j, size_t n) { return (srand(i * 137 + j), rand() % 10027); }, "A[i, j] = (srand(i * 137 + j), rand() % 10027)"},
    {[](size_t i, size_t j, size_t n) { return (i * 37) % 5 + j * 13; }, "A[i, j] = (i * 37 + j * 13)"},
    {[](size_t i, size_t j, size_t n) { return ((i * 137) % 47 + j * 13 + 176584) % 97; }, "A[i, j] = ((i * 137) % 47 + j * 13 + 176584) % 97"}
};

std::vector<VectorGenerator<double> > vector_function_generators = {
    {[](size_t i, size_t n) { return (srand(i), rand() % 10027); }, "b[i] = (srand(i), rand() % 10027)"},
    {[](size_t i, size_t n) { return 11 * i + 17; }, "b[i] = (11 * i) + 17"},
    {[](size_t i, size_t n) { return ((159 * i) % 13 + 17) % 57; }, "b[i] = ((159 * i) % 13 + 17) % 57"}
};
