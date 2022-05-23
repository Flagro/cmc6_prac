#include "operations.hpp"
#include "formula_generators.hpp"

int run_tests() {
    // dot product op testing
    static const double equal_answers_eps = 1e-7;
    static const size_t dot_product_test_n = 1000;
    static const size_t x_generator_id = 1;
    static const size_t y_generator_id = 2;
    DenseVector<double> X(vector_function_generators[x_generator_id].generate(dot_product_test_n));
    DenseVector<double> Y(vector_function_generators[y_generator_id].generate(dot_product_test_n));
    double elapsed_time;
    double result_parallel = dot_product<double>(X, Y, &elapsed_time);
    double result_sequential = dot_product_sequential<double>(X, Y, &elapsed_time);
    if (fabs(result_parallel - result_sequential) < equal_answers_eps) {
        throw "dot product testing failed (parallel version gives different result from sequential version)";
    }
    // SpMV op testing
    
    // axpby op testing

    return 0;
}
