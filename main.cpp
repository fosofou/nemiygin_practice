#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "omp.h"
#include <random>
#include <ctime>

class OpenMPSettings {
private:
    int numThreads;

public:
    OpenMPSettings();

    int getNumThreads() const;

    void setNumThreads(int threads);
};

OpenMPSettings::OpenMPSettings() : numThreads(omp_get_max_threads()) {}

int OpenMPSettings::getNumThreads() const {
    return numThreads;
}

void OpenMPSettings::setNumThreads(int threads) {
    if (threads > 0) {
        numThreads = threads;
        // omp_set_num_threads(threads);
    }
}

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

Matrix calculate(Matrix B) {
    Matrix matrix_B2 = B.square();
    Matrix matrix_B3 = B.cube();
    Matrix matrix_I = Matrix::identity(B.get_size());
    double trace_result = Matrix::trace(B.addiction(matrix_B2));
    Matrix matrix_A = matrix_B3.addiction(matrix_I.multiply(trace_result));
    return matrix_A;
}

double calculate_with_timer(Matrix B) {
    double start = omp_get_wtime();
    Matrix result = calculate(B);
    double end = omp_get_wtime();
    return end - start;
}

double get_avg_time(int runs_count, Matrix B) {
    double sum_t = 0;
    for (size_t i = 0; i < runs_count; i++) {
        sum_t += calculate_with_timer(B);
    }

    return sum_t / runs_count;
}

int main(int argc, char* argv[]) {
        // current date/time based on current system
    time_t now = time(0);

    // convert now to string form
    char* dt = ctime(&now);
    std::cout << "The local date and time is: " << dt << std::endl;
    std::vector<int> matrix_size = { 100, 300, 500, 700 };

    for (size_t i = 0; i < matrix_size.size(); i++) {
        int threads = 1;
        double prev_time = 0;
        double cur_time = RAND_MAX;
        Matrix B(matrix_size[i], matrix_size[i]);
        B.fill_random(42, 0, 1000);
        std::cout << "Matrix size: " << matrix_size[i] << std::endl;
        do {
            omp_set_num_threads(threads);
            prev_time = cur_time;
            cur_time = get_avg_time(10, B);
            std::cout << "Threads: " << threads << ", Time: " << cur_time << std::endl;

            threads++;

        } while (cur_time < prev_time);
    }

    return 0;
}

