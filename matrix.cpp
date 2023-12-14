#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include "matrix.h"


Matrix::Matrix(size_t rows, size_t cols, const std::vector<std::vector<double>> &initialData) {
    if (!initialData.empty() && (initialData.size() != rows || initialData[0].size() != cols)) {
        std::cerr << "Некорректные размеры для инициализации матрицы" << std::endl;
        exit(1);
    }
    
    this->cols = cols;
    this->rows = rows;
    if (!initialData.empty()) {
        this->data = initialData;
    } else {
        this->data = std::vector<std::vector<double>> (rows, std::vector<double>(cols, 0.0));
    }
}

void Matrix::printMatrix(const char* name) const {
    printf("Matrix %s: \n", name);
    for (const auto &row : data) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n"; 
}

void Matrix::setElement(size_t row, size_t col, double value) {
    if (row < rows && col < cols) {
        data[row][col] = value;
    } else {
        std::cerr << "Неверные индексы элементов матрицы" << std::endl;
        exit(1);
    }
}

Matrix Matrix::addiction(const Matrix &other) const {

    if (rows != other.cols || cols != other.rows) {
        std::cerr << "Невозможно сложить матрицы: неправильные размеры" << std::endl;
        exit(1);
    }

    Matrix result(rows, cols);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }

    return result;
}

Matrix Matrix::multiply(const Matrix &other) const {

    if (rows != other.cols || cols != other.rows) {
        std::cerr << "Невозможно сложить матрицы: неправильные размеры" << std::endl;
        exit(1);
    }

    Matrix result(rows, cols);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < other.cols; ++j) {
            for (size_t k = 0; k < cols; ++k) {
                result.data[i][j] += data[i][k] * other.data[k][j]; 
            }
        }
    }

    return result;
}

Matrix Matrix::multiply(double scalar) const {
    Matrix result(rows, cols);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] * scalar;
        }
    }

    return result;
}  

Matrix Matrix::square() const {
    return multiply(*this);
}

Matrix Matrix::cube() const {
    return multiply(square());
}

Matrix Matrix::identity(size_t size) {
    Matrix identityMatrix(size, size);

    for (size_t i = 0; i < size; ++i) {
        identityMatrix.data[i][i] = 1.0;
    }

    return identityMatrix;
}

double Matrix::trace(const Matrix &matrix) {
    if (matrix.rows == matrix.cols) {
        double result = 0.0;
        for (size_t i = 0; i < matrix.rows; ++i) {
            result += matrix.data[i][i];
        }
        return result;
    } else {
        std::cerr << "Невозможно вычислить след матрицы: матрица не квадратная" << std::endl;
        exit(1);    
    }
}


// int main() {
//     Matrix matrix_B = Matrix(2,2, {{4, 2}, {4.0, 5.0}});
//     matrix_B.printMatrix("B");
//     Matrix matrix_B2 = matrix_B.square();
//     matrix_B2.printMatrix("B^2");
//     Matrix matrix_B3 = matrix_B.cube();
//     matrix_B3.printMatrix("B^3");
//     Matrix matrix_I = Matrix::identity(2);
//     matrix_I.printMatrix("I");
//     double trace_result = Matrix::trace(matrix_B.addiction(matrix_B2));
//     printf("След: Tr(B + B^2) = %f \n", trace_result);
//     Matrix matrix_A = matrix_B3.addiction(matrix_I.multiply(trace_result));
//     matrix_A.printMatrix("A");
//     return 0;
// };