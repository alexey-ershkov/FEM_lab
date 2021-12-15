//
// Created by a.ershkov on 12.12.2021.
//

#ifndef MKE_MATRIX_H
#define MKE_MATRIX_H

#include <utility>
#include <vector>
#include <utility>
#include <exception>
#include <vector>

class Matrix {
public:
    Matrix(std::vector<std::vector<double>> init_matrix) : matrix(std::move(init_matrix)) {};

    Matrix operator*=(int value);

    Matrix operator*(int value);

    Matrix operator+(const Matrix &rhs_matrix);

    Matrix operator-(const Matrix &rhs_matrix);

    Matrix operator+=(const Matrix &rhs_matrix);

    Matrix operator-=(const Matrix &rhs_matrix);

    Matrix operator-();

    std::vector<std::vector<double>> getMatrix() {
        return matrix;
    }

    std::vector<double> solve_gauss(Matrix equal_to);

private:
    std::vector<std::vector<double>> matrix;
};


#endif //MKE_MATRIX_H
