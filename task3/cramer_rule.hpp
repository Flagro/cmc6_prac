// This code was copied from Rosetta Code resource and was slightly modified

#include <algorithm>
#include <iostream>
#include <vector>
 
class SubMatrix {
    const std::vector<std::vector<double>> *source;
    std::vector<double> replaceColumn;
    const SubMatrix *prev;
    size_t sz;
    int colIndex = -1;
 
public:
    SubMatrix(const std::vector<std::vector<double>> &src, const std::vector<double> &rc) : source(&src), replaceColumn(rc), prev(nullptr), colIndex(-1) {
        sz = replaceColumn.size();
    }
 
    SubMatrix(const SubMatrix &p) : source(nullptr), prev(&p), colIndex(-1) {
        sz = p.size() - 1;
    }
 
    SubMatrix(const SubMatrix &p, int deletedColumnIndex) : source(nullptr), prev(&p), colIndex(deletedColumnIndex) {
        sz = p.size() - 1;
    }
 
    int columnIndex() const {
        return colIndex;
    }
    void columnIndex(int index) {
        colIndex = index;
    }
 
    size_t size() const {
        return sz;
    }
 
    double index(int row, int col) const {
        if (source != nullptr) {
            if (col == colIndex) {
                return replaceColumn[row];
            } else {
                return (*source)[row][col];
            }
        } else {
            if (col < colIndex) {
                return prev->index(row + 1, col);
            } else {
                return prev->index(row + 1, col + 1);
            }
        }
    }
 
    double det() const {
        if (sz == 1) {
            return index(0, 0);
        }
        if (sz == 2) {
            return index(0, 0) * index(1, 1) - index(0, 1) * index(1, 0);
        }
        SubMatrix m(*this);
        double det = 0.0;
        int sign = 1;
        for (size_t c = 0; c < sz; ++c) {
            m.columnIndex(c);
            double d = m.det();
            det += index(0, c) * d * sign;
            sign = -sign;
        }
        return det;
    }
};
 
std::vector<std::pair<double, double> > solve(SubMatrix &matrix) {
    double det = matrix.det();
    if (det == 0.0) {
        throw std::runtime_error("The determinant is zero.");
    }
 
    std::vector<std::pair<double, double> > answer(matrix.size());
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix.columnIndex(i);
        answer[i] = std::make_pair(matrix.det(), det);
    }
    return answer;
}
 
std::vector<std::pair<double, double> > solveCramer(const std::vector<std::vector<double>> &matrix, const std::vector<double> &column) {
    SubMatrix sm(matrix, column);
    return solve(sm);
}
