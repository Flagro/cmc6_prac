def get_bsub_input(threads_count, parameters_string):
    if threads_count > 160:
        print("Warning: Polus can't run more than 160 omp threads, using 160 instead")
        threads_count = 160
    if threads_count > 40:
        print("Warning: Having more than 2 threads per unit is not optimal")
        core_units = (threads_count + 1) // 8
    else:
        core_units = (threads_count + 1) // 2
    if threads_count > 40:
        print("Polus doesn't have more than 20 units, running ")
    return (f"#BSUB -J \"OpenMP_job\"\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"affinity[core({core_units})]\" OMP_NUM_THREADS={threads_count}\n"
            f"/polusfs/lsf/openmp/launchOpenMP.py ./a.out {parameters_string}\n")


def get_bsub_input_slow(threads_count, parameters_string):
    if threads_count > 160:
        print("Warning: Polus can't run more than 160 omp threads, using 160 instead")
        threads_count = 160
    core_units = threads_count // 8 + 1
    return (f"#BSUB -n {core_units}\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"span[hosts=1]\"\n"
            f"OMP_NUM_THREADS={threads_count} ./a.out {parameters_string}\n")


def get_bsub_input_fast(threads_count, parameters_string):
    if threads_count > 20:
        threads_count = 20
    return (f"#BSUB -J \"OpenMP_job\"\n"
            f"#BSUB -W 0:15\n"
            f"#BSUB -o a.out.out\n"
            f"#BSUB -e a.out.err\n"
            f"#BSUB -R \"affinity[core({threads_count})]\"\n"
            f"/polusfs/lsf/openmp/launchOpenMP.py ./a.out\n {parameters_string}")


test_matrix_sizes = [2, 4, 10, 100, 1000, 2000, 3000]
test_threads_counts = [1, 2, 4, 8, 10, 16, 20, 32, 40]
valid_test_ids = [1, 2, 3]
test_test_ids = valid_test_ids[:1]

for matrix_size in test_matrix_sizes:
    for threads_count in test_threads_counts:
        for test_id in test_test_ids:
            f = open("./OpenMP_job.lsf", "w")
            f.write(get_bsub_input(threads_count, f"{matrix_size} {threads_count} {test_id} --residual --polus"))
            f.close()
            #os.system("bsub < OpenMP_job.lsf")
