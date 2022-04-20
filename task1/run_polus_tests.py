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
    return f"""#BSUB -J \"OpenMP_job\"
               #BSUB -W 0:15
               #BSUB -o a.out.out
               #BSUB -e a.out.err
               #BSUB -R “affinity[core({core_units})]” OMP_NUM_THREADS={threads_count}
               /polusfs/lsf/openmp/launchOpenMP.py ./a.out {parameters_string}"""



test_matrix_sizes = [2, 4, 10, 100, 1000, 2000, 3000]
test_threads_counts = [1, 2, 4, 8, 10, 16, 20, 32, 40]
valid_test_ids = [1, 2, 3]
test_test_ids = valid_test_ids[:1]

for matrix_size in test_matrix_sizes:
    for threads_count in test_threads_counts:
        for test_id in test_test_ids:
            print(get_bsub_input(threads_count, f"{matrix_size} {threads_count} {test_id} --residual --polus"))
