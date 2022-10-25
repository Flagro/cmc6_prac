#pragma once
#include <vector>

template <typename T>
struct EllpackMatrix {
    size_t m;
    size_t n;
    std::vector<std::vector<size_t> > ellpack_col;
    std::vector<std::vector<T> > ellpack_val;
};
