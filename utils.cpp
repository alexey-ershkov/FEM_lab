//
// Created by a.ershkov on 12.12.2021.
//

#include "utils.h"
#include "constants.h"

Matrix assemble(Matrix to_assemble, int n) {
    auto matrix_vector = to_assemble.getMatrix();

    if (matrix_vector[0].size() == 1) {
        std::vector<std::vector<double>> res(n, {0.});
        for (int i = 0; i < n - 1; ++i) {
            res[i][0] += matrix_vector[0][0];
            res[i + 1][0] += matrix_vector[1][0];
        }

        res[0][0] = U_0;
        res[n - 1][0] += A * dU_dX_at_100;

        Matrix res_matrix(res);
        return res_matrix;
    }

    if (matrix_vector[0].size() == 2) {
        std::vector<std::vector<double>> res(n, std::vector<double>(n, 0));
        for (int i = 0; i < n - 1; ++i) {
            res[i][i] += matrix_vector[0][0];
            res[i][i + 1] += matrix_vector[0][1];
            res[i + 1][i] += matrix_vector[1][0];
            res[i + 1][i + 1] += matrix_vector[1][1];
        }

        res[0][0] = 1;
        res[1][0] = 0;
        res[0][1] = 0;

        Matrix res_matrix(res);
        return res_matrix;
    }

    throw std::range_error("Can't assemble matrix of size " + std::to_string(matrix_vector[0].size()));
}

std::vector<double> error_calc(const std::vector<double> &lhs, const std::vector<double> &rhs) {
    if (rhs.size() != lhs.size()) {
        throw std::range_error("Can't count error due to sizes");
    }


    std::vector<double> error;
    error.reserve(lhs.size());

    for (int i = 0; i < lhs.size(); ++i) {
        error.push_back(abs(lhs[i] - rhs[i]));
    }

    return error;
}

void print_vector(const std::vector<double>& vector) {
    for (auto elem : vector) {
        std::cout << elem << " ";
    }

    std::cout << std::endl;
}

std::vector<double> analytical_solve(int n, double L) {
    std::vector<double> res(n, 0);
    double x = START;
    for (auto &elem: res) {
        elem = 11 + 90 / (7 * exp(140)) - (90 * exp(7 * x / 5)) / (7 * exp(140)) + 20 * x / 7;
        x += L;
    }


    return res;
}

std::vector<double> linear_solve(int n, double L) {
    Matrix quad_matrix({{1. / L,  -1. / L},
                        {-1. / L, 1. / L}});
    Matrix elementary_matrix({{-1. / 2., 1. / 2.},
                              {-1. / 2., 1. / 2.}});
    Matrix equality_matrix({{L / 2.},
                            {L / 2.}});

    quad_matrix *= A;
    elementary_matrix *= B;
    equality_matrix *= D;

    quad_matrix += elementary_matrix;

    equality_matrix = assemble(equality_matrix, n);
    quad_matrix = assemble(quad_matrix, n);


    return quad_matrix.solve_gauss(equality_matrix);
}
