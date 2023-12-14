#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows;
    size_t cols;

public:
    Matrix(size_t rows, size_t cols, const std::vector<std::vector<double>>& initialData = {});

    void printMatrix(const char* name = "") const;

    void setElement(size_t row, size_t col, double value);

    Matrix addiction(const Matrix& other) const;

    Matrix multiply(const Matrix& other) const;

    Matrix multiply(double scalar) const;

    Matrix square() const;

    Matrix cube() const;

    static Matrix identity(size_t size);

    static double trace(const Matrix& matrix);
};

#endif // MATRIX_H