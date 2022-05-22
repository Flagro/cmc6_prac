#pragma once

class Parser {
public:
    Parser(int argc, char* argv[]) {
        if (argc < 6) {
            throw "not enough pragram parameters";
        }
        n_x = strtol(argv[1], nullptr, 10);
        n_y = strtol(argv[2], nullptr, 10);
        n_z = strtol(argv[3], nullptr, 10);
        threads_num = strtol(argv[4], nullptr, 10);
        threads_mode = strtol(argv[5], nullptr, 10);
        bool calc_residual_found = false;
        bool calc_error_found = false;
        bool polus_used_found = false;
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--residual") {
                calc_residual_found = true;
            } else if (std::string(argv[i]) == "--error") {
                calc_error_found = true;
            } else if (std::string(argv[i]) == "--polus") {
                polus_used_found = true;
            } 
        }
        calc_residual = calc_residual_found;
        calc_error = calc_error_found;
        polus_used = polus_used_found;
    }

    int n_x, n_y, n_z;
    int threads_num;
    int test_id;
    int threads_mode;
    bool polus_used;
    bool calc_residual;
    bool calc_error;
};
