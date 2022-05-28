import os
import sys
import time


def get_bsub_input_simple(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon):
    core_units = threads_count // 8 + 1
    epsilon_string = "{:.16f}".format(cg_epsilon)
    return (f"#BSUB -J \"OpenMP_job_{matrix_size[0]*matrix_size[1]*matrix_size[2]}_{threads_count}_{1}\"\n"
            f"#BSUB -n {core_units}\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"span[hosts=1]\"\n"
            f"OMP_NUM_THREADS={threads_count} ./a.out {matrix_size[0]} {matrix_size[1]} {matrix_size[2]} {threads_count} 1 {max_cg_iterations} {epsilon_string} --residual --polus --bandwidth\n")


def get_bsub_input_one_threads_per_unit(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon):
    epsilon_string = "{:.16f}".format(cg_epsilon)
    return (f"#BSUB -J \"OpenMP_job_{matrix_size[0]*matrix_size[1]*matrix_size[2]}_{threads_count}_{2}\"\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n" 
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"affinity[core({threads_count})]\"\n"
            f"/polusfs/lsf/openmp/launchOpenMP.py ./a.out {matrix_size[0]} {matrix_size[1]} {matrix_size[2]} {threads_count} 2 {max_cg_iterations} {epsilon_string} --residual --polus --bandwidth\n")


def get_bsub_input_two_threads_per_unit(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon):
    core_units = (threads_count + 1) // 2
    epsilon_string = "{:.16f}".format(cg_epsilon)
    return (f"#BSUB -J \"OpenMP_job_{matrix_size[0]*matrix_size[1]*matrix_size[2]}_{threads_count}_{3}\"\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"affinity[core({core_units})]\"\n"
            f"OMP_NUM_THREADS={threads_count}\n"
            f"/polusfs/lsf/openmp/launchOpenMP.py ./a.out {matrix_size[0]} {matrix_size[1]} {matrix_size[2]} {threads_count} 3 {max_cg_iterations} {epsilon_string} --residual --polus --bandwidth\n")


def get_valid_bsub_inputs(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon):
    result = []
    if threads_count > 40:
        return result
    result.append(get_bsub_input_simple(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon))
    if threads_count <= 20:
        result.append(get_bsub_input_one_threads_per_unit(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon))
    if threads_count >= 2:
        result.append(get_bsub_input_two_threads_per_unit(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon))
    return result


if len(sys.argv) > 1 and sys.argv[1] == "--polus":
    polus_used = True
else:
    polus_used = False

test_matrix_sizes = [tuple([10, 25, 40]), tuple([40, 50, 50]), tuple([100, 100, 100])]
test_threads_counts = [1, 2, 4, 8, 10, 16, 32, 40]
valid_test_ids = [1]
test_test_ids = valid_test_ids
max_cg_iterations = 200
cg_epsilon = 0.00000001
polus_retries_cnt = 5

for matrix_size in test_matrix_sizes:
    for threads_count in test_threads_counts:
        for test_id in test_test_ids:
            for cur_mode_input in get_valid_bsub_inputs(matrix_size, threads_count, test_id, max_cg_iterations, cg_epsilon):
                if polus_used:
                    for i in range(polus_retries_cnt):
                        f = open("./OpenMP_job.lsf", "w")
                        f.write(cur_mode_input)
                        f.close()
                        os.system("bsub < OpenMP_job.lsf")
                else:
                    print(cur_mode_input)
