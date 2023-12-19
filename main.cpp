#include <iostream>
#include "matrix.h"
#include "openmp_settings.h"
#include "omp.h"
#include <vector>

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