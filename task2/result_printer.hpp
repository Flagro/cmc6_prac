#pragma once
#include <mpi.h>

void print_results(int n, int test_id, double first_time, double second_time, double residual, double error, int expected_processes_cnt) {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    FILE *my_f;
    my_f = fopen("./output.txt", "a");
    fprintf(my_f, "%d %d %d %.10lf %.10lf %.10lf %.10lf %.10lf %d\n",
            n, test_id, world_size, first_time, second_time, first_time + second_time, 
            residual, error, expected_processes_cnt);
    fclose(my_f);
}
