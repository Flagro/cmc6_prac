#include <iostream>
#include <iomanip>
#include "test.hpp"
#include "arguments_parser.hpp"
#include "input_generator.hpp"
#include "linear_system.hpp"
#include "result_printer.hpp"

using MatrixElementType = double;
using VectorElementType = double;
using ResultElementType = double;

enum { MATRIX_PRINT_COUNT = 5, VECTOR_PRINT_COUNT = 5, DOUBLE_PRINT_PRECISION = 7 };
enum { NO_CALC_RESUDUAL = -1, NO_CALC_ERROR = -1 };

int main(int argc, char *argv[]) {
    try {
        Parser parser(argc, argv);
        if (parser.run_tests) {
            try {
                run_tests();
                std::cout << "Tests successfully passed" << std::endl;
            } catch (char const* s) {
                std::cout << "OP Tests failed, reason: " << s << std::endl;
            }
        }

        InputGenerator<MatrixElementType, VectorElementType> input_generator(parser.n_x, parser.n_y, parser.n_z, parser.test_id);

        std::cout << "Generating matrix A..." << std::endl;
        auto A = input_generator.get_A();
        std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Generated matrix A:" << std::endl;
        A.print(MATRIX_PRINT_COUNT);
        std::cout << std::endl;

        std::cout << "Generating vector b..." << std::endl;
        auto b = input_generator.get_b();
        std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Generated vector b:" << std::endl;
        b.print(VECTOR_PRINT_COUNT);
        std::cout << std::endl;
        
        LinearSystem<MatrixElementType, VectorElementType, ResultElementType> linear_system(A, b);
        const auto cg_result = linear_system.solve_cg_method(parser.threads_num, parser.polus_used,
                                                             parser.calculate_bandwidth, 
                                                             parser.cg_max_iterations, parser.cg_epsilon,
                                                             parser.max_retries);

        cg_result.print(DOUBLE_PRINT_PRECISION);

        double residual = NO_CALC_RESUDUAL, error = NO_CALC_ERROR;
        if (parser.calc_residual) {
            residual = linear_system.calculate_residual(cg_result.x);
            std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Residual: " << residual << std::endl;
        }
        if (parser.calc_error) {
            error = linear_system.calculate_error(cg_result.x);
            std::cout << std::setprecision(DOUBLE_PRINT_PRECISION)<< "Error: " << error << std::endl;
        }

        print_results(parser, cg_result, residual, error);
    } catch (char const* s) {
        std::cout << "Error occured: " << s << std::endl;
    }
    return 0;
}
