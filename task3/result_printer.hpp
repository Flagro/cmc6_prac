#pragma once
#include "omp.h"
#include <cstdio>

void print_results(int n_x, int n_y, int n_z, int test_id, 
        double first_time, double second_time, double residual, double error, int mode, int passed_threads_num) {
    int threads_num = 0;
    #pragma omp parallel
    {
        #pragma omp single
        threads_num = omp_get_num_threads();
    }

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
    fprintf(my_f, "%d %d %d %.10lf %.10lf %.10lf %.10lf %.10lf %s %d\n", 
            n_x, n_y, n_z, test_id, threads_num, first_time, second_time, first_time + second_time, 
            residual, error, mode_string, passed_threads_num);
    fclose(my_f);
}
