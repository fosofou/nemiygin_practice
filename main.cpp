#include <iostream>
#include "matrix.h"

int main(int, char**){
    Matrix matrix_B = Matrix(2,2, {{4, 2}, {4.0, 5.0}});
    matrix_B.printMatrix("B");
    Matrix matrix_B2 = matrix_B.square();
    matrix_B2.printMatrix("B^2");
    Matrix matrix_B3 = matrix_B.cube();
    matrix_B3.printMatrix("B^3");
    Matrix matrix_I = Matrix::identity(2);
    matrix_I.printMatrix("I");
    double trace_result = Matrix::trace(matrix_B.addiction(matrix_B2));
    printf("След: Tr(B + B^2) = %f \n", trace_result);
    Matrix matrix_A = matrix_B3.addiction(matrix_I.multiply(trace_result));
    matrix_A.printMatrix("A");
    return 0;
}
