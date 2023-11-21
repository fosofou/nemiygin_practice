#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows;
    size_t cols;
public:
    // конструктор
    Matrix(size_t rows, size_t cols, const std::vector<std::vector<double>> &initialData = {}) {
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


    void printMatrix() const {
        for (const auto &row : data) {
            for (double elem : row) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "\n"; 
    }

    void setElement(size_t row, size_t col, double value) {
        if (row < rows && col < cols) {
            data[row][col] = value;
        } else {
            std::cerr << "Неверные индексы элементов матрицы" << std::endl;
            exit(1);
        }
    }


    Matrix addiction(const Matrix &other) const {

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

    Matrix multiply(const Matrix &other) const {

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

    Matrix multiply(double scalar) const {
        Matrix result(rows, cols);

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }

        return result;
    }  

    Matrix square() const {
        return multiply(*this);
    }

    Matrix cube() const {
        return multiply(square());
    }

    static Matrix identity(size_t size) {
        Matrix identityMatrix(size, size);

        for (size_t i = 0; i < size; ++i) {
            identityMatrix.data[i][i] = 1.0;
        }

        return identityMatrix;
    }

    static double trace(const Matrix &matrix) {
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
};

int main() {
    // std::vector<std::vector<double>> B = {{1.0, 2.0, },
    //                                       {4.0, 5.0, }};
    Matrix* matrix_B = new Matrix(2,2, {{4, 2}, {4.0, 5.0}});
    // matrix_B->setElement(0,0,1);
    printf("Матрица B \n");
    matrix_B->printMatrix();
    Matrix matrix_B2 = matrix_B->square();
    printf("Матрица B^2 \n");
    matrix_B2.printMatrix();
    Matrix matrix_B3 = matrix_B->cube();
    printf("Матрица B^3 \n");
    matrix_B3.printMatrix();
    Matrix matrix_I = Matrix::identity(2);
    Matrix multiplyB_B2 = matrix_B->addiction(matrix_B2);
    printf("сложение B + B^2 \n");
    multiplyB_B2.printMatrix();
    double trace_result = Matrix::trace(multiplyB_B2);
    printf("След %f \n", trace_result);
    Matrix multiplyTrace_I = matrix_I.multiply(trace_result);
    printf("Единичная матрица * след \n");
    multiplyTrace_I.printMatrix();
    Matrix matrix_A = matrix_B3.addiction(multiplyTrace_I);
    printf("Матрица А \n");
    matrix_A.printMatrix();
    return 0;
};