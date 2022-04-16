#include <vector>

template <typename matrix_T>
class BaseMatrixGenerator {
public:
    BaseMatrixGenerator(): n(0) {}
    BaseMatrixGenerator(std::size_t output_size): n(output_size) {}
    virtual ~BaseMatrixGenerator() = default;
    const std::size_t n;

    virtual get(size i, size_t j) = 0;
    virtual std::string get_representation() = 0;
};

template <typename vector_T>
class BaseVectorGenerator {
public:
    BaseVectorGenerator(): n(0) {}
    BaseVectorGenerator(std::size_t output_size): n(output_size) {}
    virtual ~BaseVectorGenerator() = default;
    const std::size_t n;

    virtual get(size i) = 0;
    virtual std::string get_representation() = 0;
};

template <typename matrix_T>
class MatrixGenerator1 {
public:
    MatrixGenerator1(): n(0) {}
    MatrixGenerator1(std::size_t output_size): n(output_size) {}
    const std::size_t n;

    get(size i, size_t j) {
        return (i * 1230) % j + 12;
    }

    std::string get_representation() {
        return "A[i, j] = (i * 1230) % j + 12";
    }
};

template <typename vector_T>
class VectorGenerator1 {
public:
    VectorGenerator1(): n(0) {}
    VectorGenerator1(std::size_t output_size): n(output_size) {}
    const std::size_t n;

    get(size i) {
        return (23123 * i) % 120;
    }
    std::string get_representation() {
        return "b[i] = (23123 * i) % 120";
    }
};
