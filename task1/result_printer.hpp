#pragma once

void print_results(int n, int test_id, double first_time, double second_time, double residual, double error) {
    int threads_num = 0;
    #pragma omp parallel
    {
        #pragma omp single
        threads_num = omp_get_num_threads();
    }
    /*
    FILE *my_f;
    my_f = fopen("./output.txt", "a");
    fprintf(my_f, "%d %d %d %.10lf %.10lf %.10lf %.10lf %.10lf\n", n, test_id, threads_num, first_time, second_time, first_time + second_time, residual, error);
    fclose(my_f);
    */
}
