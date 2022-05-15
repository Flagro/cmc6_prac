#pragma once
#include <iostream>
#include <cmath>
#include <mpi.h>
#include "matrix.hpp"
#include "vector.hpp"
#include "cramer_rule.hpp"

template <typename matrix_T, typename vector_T, typename result_T>
class LinearSystem {
public:
    LinearSystem(const Matrix<matrix_T>& A, const Vector<vector_T>& b): _A(A), _b(b), _n(A.size()) {
        if (A.size() != b.size()) {
            throw "A and b does not have same dimensions";
        }
    }

    std::unique_ptr<Vector<result_T> > solve_reflection_method(double *first_stage_elapsed_time, 
                                                               double *second_stage_elapsed_time) {
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        size_t b_column_index = _n;
        size_t columns_per_process = _n / world_size + 1;

        std::vector<std::vector<result_T> > process_data(columns_per_process, std::vector<result_T>(_n));
        // q-th process should contain columns with ids: q, p + q, 2p + q, ...
        for (size_t i = 0; i < columns_per_process; ++i) {
            size_t cur_column_id = i * world_size + world_rank;
            if (cur_column_id > b_column_index) {
                break;
            } else if (cur_column_id == b_column_index) {
                for (size_t j = 0 ; j < _n; ++j) {
                    process_data[i][j] = _b.get(j);
                }
            } else {
                for (size_t j = 0 ; j < _n; ++j) {
                    process_data[i][j] = _A.get(j, cur_column_id);
                }
            }
        }

        std::vector<result_T> x(_n, 0);

        MPI_Barrier(MPI_COMM_WORLD);
        auto first_stage_elapsed_time_first = MPI_Wtime();

        for (size_t i = 0; i + 1 < _n; ++i) {
            // process that has i-th column:
            size_t cur_column_process_rank = i % world_size;
            if ((size_t) world_rank == cur_column_process_rank) {
                result_T s = 0;
                size_t ith_column_local_id = i / world_size;
                for (size_t j = i + 1; j < _n; ++j) {
                    s += process_data[ith_column_local_id][j] * process_data[ith_column_local_id][j];
                }
                result_T cur_a_norm = sqrt(process_data[ith_column_local_id][i] * process_data[ith_column_local_id][i] + s);
                result_T cur_x1 = process_data[ith_column_local_id][i] - cur_a_norm;
                result_T cur_x_norm = sqrt(cur_x1 * cur_x1 + s);
                x[i] = cur_x1 / cur_x_norm;

                for (size_t j = i + 1; j < _n; ++j) {
                    x[j] = process_data[ith_column_local_id][j] / cur_x_norm;
                }

                process_data[ith_column_local_id][i] = cur_a_norm;
                for (size_t j = i + 1; j < _n; ++j) {
                    process_data[ith_column_local_id][j] = 0;
                }
            }

            MPI_Bcast(&x[i], _n - i, MPI_DOUBLE, cur_column_process_rank, MPI_COMM_WORLD);

            // doing mxv (U(x)A_{*j}) = A_{*j} - 2x(A_{*j}, x) for each j-th column from A
            // starting from i + 1 that means that q-th process will only require ids from (i + 1 - q) / p
            size_t cur_data_start_id = (i + 1) / world_size;
            for (size_t j = cur_data_start_id; j < columns_per_process; ++j) {
                size_t cur_column_id = j * world_size + world_rank;
                if (cur_column_id < i + 1) {
                    continue;
                }
                if (cur_column_id > b_column_index) {
                    break;
                }
                result_T local_alpha = 0;
                for (size_t k = i; k < _n; ++k) {
                    local_alpha += process_data[j][k] * x[k];
                }
                for (size_t k = i; k < _n; ++k) {
                    process_data[j][k] -= 2 * x[k] * local_alpha;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        auto first_stage_elapsed_time_second = MPI_Wtime();
        if (!world_rank) {
            std::cout << "Finished triangulation" << std::endl;
        }

        std::vector<result_T> local_result_array(columns_per_process, 0);
        MPI_Barrier(MPI_COMM_WORLD);
        auto second_stage_elapsed_time_first = MPI_Wtime();

        for (int i = (int)_n - 1; i >= 0; --i) {
            result_T local_sum = 0;
            size_t cur_data_start_id = (i + 1 - world_rank) / world_size;
            for (size_t j = cur_data_start_id; j < columns_per_process; ++j) {
                size_t cur_column_id = j * world_size + world_rank;
                if (cur_column_id > b_column_index) {
                    break;
                } else if (cur_column_id == b_column_index) {
                    local_sum += process_data[j][i];
                } else {
                    local_sum -= process_data[j][i] * local_result_array[j];
                }
            }
            result_T global_sum = 0;
            size_t cur_column_process_rank = i % world_size;
            MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, cur_column_process_rank, MPI_COMM_WORLD);
            if ((size_t) world_rank == cur_column_process_rank) {
                size_t cur_ans_id = i / world_size;
                local_result_array[cur_ans_id] = global_sum / process_data[cur_ans_id][i];
            }
        }
        std::vector<result_T> result_array(_n, 0);
        for (size_t i = 0; i < _n; ++i) {
            size_t cur_column_process_rank = i % world_size;
            if ((size_t) world_rank == cur_column_process_rank) {
                result_array[i] = local_result_array[i / world_size];
            }
            MPI_Bcast(&result_array[i], 1, MPI_DOUBLE, cur_column_process_rank, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        auto second_stage_elapsed_time_second = MPI_Wtime();

        *first_stage_elapsed_time = first_stage_elapsed_time_second - first_stage_elapsed_time_first;
        *second_stage_elapsed_time = second_stage_elapsed_time_second - second_stage_elapsed_time_first;
        std::unique_ptr<Vector<result_T> > result_uptr(new ArrayVector<result_T>(_n, result_array));
        return result_uptr;
    }

    double calculate_residual(const Vector<result_T>& x) {
        std::vector<result_T> residual(_n, 0);
        for (size_t i = 0; i < _n; ++i) {
            for (size_t j = 0; j < _n; ++j) {
                residual[i] += _A.get(i, j) * x.get(j);
            }
            residual[i] -= _b.get(i);
        }
        result_T residual_norm = 0;
        for (size_t i = 0; i < _n; ++i) {
            residual_norm += residual[i] * residual[i];
        }
        return sqrt(residual_norm);
    }

    double calculate_error(const Vector<result_T>& x) {
        std::vector<std::vector<result_T> > result_matrix(_n, std::vector<result_T>(_n, 0));
        for (size_t i = 0; i < _n; ++i) {
            for (size_t j = 0 ; j < _n; ++j) {
                result_matrix[i][j] = _A.get(i, j);
            }
        }
        std::vector<result_T> result_vector(_n, 0);
        for (size_t i = 0; i < _n; ++i) {
            result_vector[i] = _b.get(i);
        }
        std::vector<std::pair<double, double> > almost_exact_result;
        try {
            almost_exact_result = solveCramer(result_matrix, result_vector);
        } catch(std::runtime_error &e) {
            std::cout << "Error occured, Not calculating the error" << std::endl;
            return 0;
        }
        std::cout << "Found an answer with cramer's rule:" << std::endl;
        for (size_t i = 0; i < _n; ++i) {
            std::cout << almost_exact_result[i].first / almost_exact_result[i].second << " ";
        }
        std::cout << std::endl;
        double error = 0;
        // comparing x_i and (a_1i / a_2i) means comparing (x_i * a_2i and a_1i)
        for (size_t i = 0; i < _n; ++i) {
            double cur_difference = (x.get(i) * almost_exact_result[i].second - almost_exact_result[i].first) / almost_exact_result[i].second;
            error += cur_difference * cur_difference;
        }
        return sqrt(error);
    }

private:
    const Matrix<matrix_T>& _A;
    const Vector<vector_T>& _b;
    const size_t _n;
};
