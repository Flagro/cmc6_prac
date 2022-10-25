#pragma once

enum { DEFAULT_TEST_ID = 1, DEFAULT_MAX_RETIRES = 100 };

class Parser {
public:
    Parser(int argc, char* argv[]) {
        if (argc < 8) {
            throw "not enough pragram parameters";
        }
        n_x = strtol(argv[1], nullptr, 10);
        n_y = strtol(argv[2], nullptr, 10);
        n_z = strtol(argv[3], nullptr, 10);
        threads_num = strtol(argv[4], nullptr, 10);
        threads_mode = strtol(argv[5], nullptr, 10);
        cg_max_iterations = strtol(argv[6], nullptr, 10);
        cg_epsilon = strtod(argv[7], nullptr);
        bool calc_residual_found = false;
        bool calc_error_found = false;
        bool polus_used_found = false;
        bool run_tests_found = false;
        bool calculate_bandwidth_found = false;
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--residual") {
                calc_residual_found = true;
            } else if (std::string(argv[i]) == "--error") {
                calc_error_found = true;
            } else if (std::string(argv[i]) == "--polus") {
                polus_used_found = true;
            } else if (std::string(argv[i]) == "--test") {
                run_tests_found = true;
            } else if (std::string(argv[i]) == "--bandwidth") {
                calculate_bandwidth_found = true;
            }
        }
        calc_residual = calc_residual_found;
        calc_error = calc_error_found;
        polus_used = polus_used_found;
        run_tests = run_tests_found;
        calculate_bandwidth = calculate_bandwidth_found;
        test_id = DEFAULT_TEST_ID;
        max_retries = DEFAULT_MAX_RETIRES;
    }

    int n_x, n_y, n_z;
    int threads_num;
    int test_id;
    int threads_mode;
    bool polus_used;
    bool calc_residual;
    bool calc_error;
    bool run_tests;
    bool calculate_bandwidth;
    int cg_max_iterations;
    double cg_epsilon;
    int max_retries;
};
