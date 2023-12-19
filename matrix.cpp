#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include "matrix.h"
#include "openmp_settings.h"
#include "omp.h"
#include <random>

Matrix::Matrix(
    size_t rows, 
    size_t cols, 
    const std::vector<std::vector<double>> &initialData
    ) {
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

int Matrix::get_size() const {
    return rows;
}

void Matrix::print_matrix(const char* name) const {
    printf("Matrix %s: \n", name);
    for (const auto &row : data) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n"; 
}

void Matrix::set_element(size_t row, size_t col, double value) {
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
    size_t i,j;

    #pragma omp parallel for private(i,j) shared(result, other, data)
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
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

    double start_time = omp_get_wtime();
    size_t i,j,k;

    #pragma omp parallel for private(i,j,k) shared(result, other, data) 
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < other.cols; ++j) {
            for (k = 0; k < cols; ++k) {
                result.data[i][j] += data[i][k] * other.data[k][j]; 
            }
        }
    }
        
    return result;
}


Matrix Matrix::multiply(double scalar) const {
    Matrix result(rows, cols);
    size_t i,j,k;
    #pragma omp parallel for private(i,j,k) shared(result, data, scalar)
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] * scalar;
        }
    }

    return result;
}  

void Matrix::fill_random(unsigned int seed, double minValue, double maxValue) {
    // Установка seed
    srand(seed);
    #pragma omp parallel for collapse(2) shared(data)
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            double randomNumber = minValue + static_cast<double>(rand()) / RAND_MAX * (maxValue - minValue);
            #pragma omp critical
            {
                data[i][j] = randomNumber;
            }
        }
    }
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
        size_t i;
        #pragma omp parallel for private(i) reduction (+:result)
        for (i = 0; i < matrix.rows; ++i) {
            result += matrix.data[i][i];
        }
        return result;
    } else {
        throw std::invalid_argument( "Невозможно вычислить след матрицы: матрица не квадратная" );    
    }
}

