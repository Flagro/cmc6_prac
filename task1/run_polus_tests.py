import os
import sys
import time


def get_bsub_input_simple(matrix_size, threads_count, test_id):
    core_units = threads_count // 8 + 1
    return (f"#BSUB -J \"OpenMP_job\"\n"
            f"#BSUB -n {core_units}\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"span[hosts=1]\"\n"
            f"OMP_NUM_THREADS={threads_count} ./a.out {matrix_size} {threads_count} {test_id} 1 --residual --polus\n")


def get_bsub_input_one_threads_per_unit(matrix_size, threads_count, test_id):
    return (f"#BSUB -J \"OpenMP_job\"\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n" 
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"affinity[core({threads_count})]\"\n"
            f"/polusfs/lsf/openmp/launchOpenMP.py ./a.out {matrix_size} {threads_count} {test_id} 2 --residual --polus\n")


def get_bsub_input_two_threads_per_unit(matrix_size, threads_count, test_id):
    core_units = (threads_count + 1) // 2
    return (f"#BSUB -J \"OpenMP_job\"\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"affinity[core({core_units})]\" OMP_NUM_THREADS={threads_count}\n"
            f"/polusfs/lsf/openmp/launchOpenMP.py ./a.out {matrix_size} {threads_count} {test_id} 3 --residual --polus\n")


def get_valid_bsub_inputs(matrix_size, threads_count, test_id):
    result = []
    if threads_count > 40:
        return result
    result.append(get_bsub_input_simple(matrix_size, threads_count, test_id))
    if threads_count <= 20:
        result.append(get_bsub_input_one_threads_per_unit(matrix_size, threads_count, test_id))
    if threads_count >= 2:
        result.append(get_bsub_input_two_threads_per_unit(matrix_size, threads_count, test_id))
    return result


if len(sys.argv) > 1 and sys.argv[1] == "--polus":
    polus_used = True
else:
    polus_used = False

test_matrix_sizes = [100, 1000, 3000, 5000]
test_threads_counts = [1, 2, 4, 8, 10, 16, 32, 40]
valid_test_ids = [1, 2, 3]
test_test_ids = valid_test_ids[:1]
polus_retries_cnt = 3

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
                        time.sleep(10)
                else:
                    print(cur_mode_input)
