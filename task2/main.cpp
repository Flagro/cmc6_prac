#include <iostream>
#include <iomanip>
#include <mpi.h>
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
        MPI_Init(NULL, NULL);
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        Parser parser(argc, argv);
        InputGenerator<MatrixElementType, VectorElementType> input_generator(parser.n, parser.test_id);
        const auto A_uptr = input_generator.get_A();
        const auto b_uptr = input_generator.get_b();

        if (!world_rank) {
            std::cout << std::setprecision(7) << "Generated matrix A:" << std::endl;
            A_uptr->print(MATRIX_PRINT_COUNT);
            std::cout << std::endl;
            std::cout << std::setprecision(7) << "Generated vector b:" << std::endl;
            b_uptr->print(VECTOR_PRINT_COUNT);
            std::cout << std::endl;
        }

        LinearSystem<MatrixElementType, VectorElementType, ResultElementType> linear_system(*A_uptr, *b_uptr);
        double first_stage_elapsed_time, second_stage_elapsed_time;
        const auto x_uptr = linear_system.solve_reflection_method(parser.threads_num, parser.polus_used, 
                                                                &first_stage_elapsed_time, 
                                                                &second_stage_elapsed_time);

        if (!world_rank) {
            std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) << "Calculated vector x:" << std::endl;
            x_uptr->print(VECTOR_PRINT_COUNT);
            std::cout << std::endl;
            std::cout << std::setprecision(DOUBLE_PRINT_PRECISION) 
                    << "Upper triangulation elapsed time: " << first_stage_elapsed_time << std::endl
                    << "Backward Gauss move elapsed time: " << second_stage_elapsed_time << std::endl;
            double residual = NO_CALC_RESUDUAL, error = NO_CALC_ERROR;
            if (parser.calc_residual) {
                residual = linear_system.calculate_residual(*x_uptr);
                std::cout << std::setprecision(7) << "Residual: " << residual << std::endl;
            }
            if (parser.calc_error) {
                error = linear_system.calculate_error(*x_uptr);
                std::cout << std::setprecision(7)<< "Error: " << error << std::endl;
            }

            print_results(parser.n, parser.test_id, first_stage_elapsed_time, second_stage_elapsed_time, 
                        residual, error, parser.threads_mode, parser.threads_num);
        }
    } catch(char const* s) {
        std::cout << "Error occured: " << s << std::endl;
    }
    MPI_Finalize();
    return 0;
}
