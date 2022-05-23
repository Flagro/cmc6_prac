#include "operations.hpp"
#include "operations_safe.hpp"
#include "formula_generators.hpp"

void run_tests() {
    static const double equal_answers_eps = 1e-5;
    double elapsed_time;
    // dot product op testing
    static const size_t dot_product_test_n = 1000;
    static const size_t x_generator_id = 1;
    static const size_t y_generator_id = 2;
    DenseVector<double> X(vector_function_generators[x_generator_id].generate(dot_product_test_n));
    DenseVector<double> Y(vector_function_generators[y_generator_id].generate(dot_product_test_n));
    double result_parallel = dot_product<double>(X, Y, &elapsed_time);
    double result_sequential = dot_product_sequential<double>(X, Y, &elapsed_time);
    if (fabs(result_parallel - result_sequential) >= equal_answers_eps) {
        std::cout << result_parallel << " " << result_sequential << std::endl;
        throw "dot product testing failed (parallel version gives different result from sequential version)";
    }

    // SpMV op testing
    static const size_t spmv_test_n_x = 30, spmv_test_n_y = 40, spmv_test_n_z = 50;
    static const size_t a_generator_id = 0;
    static const size_t b_generator_id = 0;
    SparseMatrix<double> spmv_a(cube_matrix_function_generators[a_generator_id].generate(spmv_test_n_x, spmv_test_n_y, spmv_test_n_z));
    DenseVector<double> spmv_b(vector_function_generators[b_generator_id].generate(spmv_test_n_x * spmv_test_n_y * spmv_test_n_z));
    DenseVector<double> vector_result_parallel = SpMV<double>(spmv_a, spmv_b, &elapsed_time);
    DenseVector<double> vector_result_sequential = SpMV_sequential<double>(spmv_a, spmv_b, &elapsed_time);
    if (fabs(vector_result_parallel.norm_l2() - vector_result_sequential.norm_l2()) >= equal_answers_eps) {
        throw "SpMV testing failed (parallel version gives different result from sequential version)";
    }

    // axpby op testing
    static const double axbpy_a = 532.21;
    static const double axbpy_b = 0.7;
    DenseVector<double> axpby_result_parallel = axpby<double>(axbpy_a, X, axbpy_b, Y, &elapsed_time);
    DenseVector<double> axpby_result_sequential = axpby_sequential<double>(axbpy_a, X, axbpy_b, Y, &elapsed_time);
    if (fabs(axpby_result_parallel.norm_l2() - axpby_result_sequential.norm_l2()) >= equal_answers_eps) {
        throw "axpby testing failed (parallel version gives different result from sequential version)";
    }
}
