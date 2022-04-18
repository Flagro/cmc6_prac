#pragma once
#include <vector>
#include <functional>

template <typename matrix_T>
struct MatrixGenerator {
public:
    MatrixGenerator(const std::function<matrix_T(size_t, size_t, [[maybe_unused]] size_t)>& get_function, 
                    const std::string& representation_string): get(get_function), 
                                                               function_representation(representation_string) {}
    const std::function<matrix_T(size_t, size_t, size_t)> get;
    const std::string function_representation;
};

template <typename vector_T>
struct VectorGenerator {
public:
    VectorGenerator(const std::function<vector_T(size_t, [[maybe_unused]] size_t)>& get_function, 
                    const std::string& representation_string): get(get_function), 
                                                               function_representation(representation_string) {}
    const std::function<vector_T(size_t, size_t)> get;
    const std::string function_representation;
};

std::vector<MatrixGenerator<double> > matrix_function_generators = {
    {[](size_t i, size_t j, [[maybe_unused]] size_t n) { return rand() % 12381245823; }, "A[i, j] = TEST"},
    {[](size_t i, size_t j, [[maybe_unused]] size_t n) { return (i * 37) % 5 + j * 13; }, "A[i, j] = (i * 37 + j * 13)"},
    {[](size_t i, size_t j, [[maybe_unused]] size_t n) { return ((i * 137) % 47 + j * 13 + 176584) % 97; }, "A[i, j] = TEST"}
};

std::vector<VectorGenerator<double> > vector_function_generators = {
    {[](size_t i, [[maybe_unused]] size_t n) { return rand() % 12381245823; }, "b[i] = TEST"},
    {[](size_t i, [[maybe_unused]] size_t n) { return 11 * i + 17; }, "b[i] = (11 * i) + 17"},
    {[](size_t i, [[maybe_unused]] size_t n) { return ((159 * i) % 13 + 17) % 57; }, "b[i] = TEST"}
};
