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

enum { DEFAULT_TEST_ID = 1 };
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

        InputGenerator<MatrixElementType, VectorElementType> input_generator(parser.n_x, parser.n_y, parser.n_z, DEFAULT_TEST_ID);
        auto A = input_generator.get_A();
        auto b = input_generator.get_b();

        std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Generated matrix A:" << std::endl;
        A.print(MATRIX_PRINT_COUNT);
        std::cout << std::endl;

        std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Generated vector b:" << std::endl;
        b.print(VECTOR_PRINT_COUNT);
        std::cout << std::endl;

        LinearSystem<MatrixElementType, VectorElementType, ResultElementType> linear_system(A, b);
        double first_stage_elapsed_time, second_stage_elapsed_time;
        const auto x = linear_system.solve_cg_method(parser.threads_num, parser.polus_used,
                                                     &first_stage_elapsed_time,
                                                     &second_stage_elapsed_time,
                                                     parser.calculate_bandwidth, 100, 0.0001);

        std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Calculated vector x:" << std::endl;
        x.print(VECTOR_PRINT_COUNT);
        std::cout << std::endl;

        std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) 
                  << "Upper triangulation elapsed time: " << first_stage_elapsed_time << std::endl
                  << "Backward Gauss move elapsed time: " << second_stage_elapsed_time << std::endl;

        double residual = NO_CALC_RESUDUAL, error = NO_CALC_ERROR;
        if (parser.calc_residual) {
            residual = linear_system.calculate_residual(x);
            std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Residual: " << residual << std::endl;
        }
        if (parser.calc_error) {
            error = linear_system.calculate_error(x);
            std::cout << std::setprecision(DOUBLE_PRINT_PRECISION)<< "Error: " << error << std::endl;
        }

        //print_results(parser.n_x, parser.n_y, parser.n_z, DEFAULT_TEST_ID, 
        //              first_stage_elapsed_time, second_stage_elapsed_time, 
        //              residual, error, parser.threads_mode, parser.threads_num);
    } catch (char const* s) {
        std::cout << "Error occured: " << s << std::endl;
    }
    return 0;
}
