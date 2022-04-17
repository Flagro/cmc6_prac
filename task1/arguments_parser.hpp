#pragma once

class Parser {
public:
    Parser(int argc, char* argv[]) {
        if (argc < 4) {
            throw "not enough pragram parameters";
        }
        n = strtol(argv[1], nullptr, 10);
        threads_num = strtol(argv[2], nullptr, 10);
        test_id = strtol(argv[3], nullptr, 10);
    }

    int n;
    int threads_num;
    int test_id;
};
