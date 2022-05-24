#pragma once
#include "arguments_parser.hpp"
#include "linear_system.hpp"
#include "omp.h"
#include <cstdio>

template <typename T>
void print_results(const Parser &parser, const CGSolverResult<T> &cg_result, double residual, double error) {
    int actual_threads_num = 0;
    #pragma omp parallel
    {
        #pragma omp single
        actual_threads_num = omp_get_num_threads();
    }

    int mode = parser.threads_mode;

    const char mode1_str[] = "Without_Threads_Anchorage";
    const char mode2_str[] = "1_Thread_per_Unit";
    const char mode3_str[] = "2_Threads_per_Unit";
    const char *mode_string;
    switch (mode) {
        case 1:
            mode_string = mode1_str;
            break;
        case 2:
            mode_string = mode2_str;
            break;
        case 3:
            mode_string = mode3_str;
            break;
        default:
            break;
    }

    FILE *my_f;
    my_f = fopen("./output.txt", "a");
    fprintf(my_f, "%d %d %d %d %d %d %d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %d %.10lf %.10lf %.10lf %s\n", 
            parser.n_x, parser.n_y, parser.n_z, 
            parser.test_id, actual_threads_num, parser.threads_num, 
            parser.cg_max_iterations, parser.cg_epsilon, 
            cg_result.total_cg_time, 
            cg_result.dot_product_elapsed_time, cg_result.spmv_elapsed_time, cg_result.axpby_elapsed_time,
            cg_result.dot_product_bandwidth, cg_result.spmv_bandwidth, cg_result.axpby_bandwidth,
            cg_result.iterations_taken, cg_result.cg_residual,
            residual, error, mode_string);
    fclose(my_f);
}
