#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include "openmp_settings.h"

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows;
    size_t cols;

public:
    Matrix(
        size_t rows, 
        size_t cols, 
        const std::vector<std::vector<double>>& initialData = {}
    );
    
    int get_size() const;

    void print_matrix(const char* name = "") const;

    void set_element(size_t row, size_t col, double value);

    Matrix addiction(const Matrix& other) const;

    Matrix multiply(const Matrix& other) const;

    Matrix multiply(double scalar) const;

    Matrix square() const;

    Matrix cube() const;

    void fill_random(unsigned int seed, double minValue = 0, double maxValue = 1);

    static Matrix identity(size_t size);

    static double trace(const Matrix& matrix);
};

#endif // MATRIX_H