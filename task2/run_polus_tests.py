import os
import sys
import time


def get_mpi_bsub_input(matrix_size, mpi_processes, test_id):
    return (f"#BSUB -J \"OpenMPI_job_{matrix_size}_{mpi_processes}\"\n"
            f"#BSUB -n {mpi_processes}\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"mpiexec ./a.out {matrix_size} {mpi_processes} {test_id} --residual --polus\n")


def get_valid_bsub_inputs(matrix_size, threads_count, test_id):
    result = []
    if threads_count > 40:
        return result
    result.append(get_mpi_bsub_input(matrix_size, threads_count, test_id))
    return result


if len(sys.argv) > 1 and sys.argv[1] == "--polus":
    polus_used = True
else:
    polus_used = False

test_matrix_sizes = [100, 1000, 3000]
test_threads_counts = [1, 2, 4, 8, 10, 16, 32, 40]
valid_test_ids = [1, 2, 3]
test_test_ids = valid_test_ids[:1]
polus_retries_cnt = 2

for matrix_size in test_matrix_sizes:
    for threads_count in test_threads_counts:
        for test_id in test_test_ids:
            for cur_mode_input in get_valid_bsub_inputs(matrix_size, threads_count, test_id):
                if polus_used:
                    for i in range(polus_retries_cnt):
                        f = open("./OpenMP_job.lsf", "w")
                        f.write(cur_mode_input)
                        f.close()
                        os.system("bsub < OpenMP_job.lsf")
                        time.sleep(1)
                else:
                    print(cur_mode_input)
