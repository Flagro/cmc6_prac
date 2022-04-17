#include <iostream>
#include <iomanip>
#include "input_generator.hpp"
#include "linear_system.hpp"
#include "result_printer.hpp"
#include "arguments_parser.hpp"

using MatrixElementType = double;
using VectorElementType = double;
using ResultElementType = double;

enum { MATRIX_PRINT_COUNT = 5, VECTOR_PRINT_COUNT = 5 };

int main(int argc, char *argv[]) {
    try {
        Parser parser(argc, argv);
        InputGenerator<MatrixElementType, VectorElementType> input_generator(parser.n, parser.test_id);
        const auto A = input_generator.get_A();
        const auto b = input_generator.get_b();

        std::cout << std::setprecision(7) << "Generated matrix A:" << std::endl;
        A.print(MATRIX_PRINT_COUNT);
        std::cout << std::setprecision(7) << "Generated vector b:" << std::endl;
        b.print(VECTOR_PRINT_COUNT);

        LinearSystem<MatrixElementType, VectorElementType, ResultElementType> linear_system(A, b);
        double first_stage_elapsed_time, second_stage_elapsed_time;
        const auto x = linear_system.solve_reflection_method(&first_stage_elapsed_time, 
                                                             &second_stage_elapsed_time);

        std::cout << std::setprecision(7) << "Calculated vector x:" << std::endl;
        x.print(VECTOR_PRINT_COUNT);

        std::cout << std::setprecision(7) << "Upper triangulation elapsed time: " << first_stage_elapsed_time << std::endl
                                          << "Backward Gauss move elapsed time: " << second_stage_elapsed_time << std::endl;

        double residual = linear_system.calculate_residual(x);
        double error = linear_system.calculate_error(x);

        std::cout << std::setprecision(7) << "Residual: " << residual << std::endl
                                          << "Error: " << error << std::endl;

        print_json_results();
    } catch(char const* s) {
        std::cout << "Error occured: " << s << std::endl;
    }
    return 0;
}
