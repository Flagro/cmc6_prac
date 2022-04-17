#include <vector>
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
class VectorGenerator {
public:
    VectorGenerator(const std::function<vector_T(size_t, size_t)>& get_function, 
                    const std::string& representation_string): get(get_function), 
                                                               function_representation(representation_string) {}
    const std::function<vector_T(size_t, size_t)> get;
    const std::string function_representation;
};

std::vector<MatrixGenerator<double> > matrix_function_generators = {
    {[](size_t i, size_t j, size_t n) { return (i * 1230) % j + 12; }, "A[i, j] = (i * 1230) % j + 12"}
};

std::vector<VectorGenerator<double> > vector_function_generators = {
    {[](size_t i, size_t n) { return (23123 * i) % 120; }, "b[i] = (23123 * i) % 120"}
};
