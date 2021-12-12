//
// Created by a.ershkov on 12.12.2021.
//

#include "utils.h"
#include "constants.h"
#include <fstream>

Matrix assemble(Matrix to_assemble, int n) {
    auto matrix_vector = to_assemble.getMatrix();

    if (matrix_vector[0].size() == 1 && matrix_vector.size() == 2) {
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

    if (matrix_vector[0].size() == 1 && matrix_vector.size() == 4) {
        std::vector<std::vector<double>> res(n * 3 + 1, {0.});
        for (int i = 0; i < n * 3; i += 3) {
            res[i][0] += matrix_vector[0][0];
            res[i + 1][0] += matrix_vector[1][0];
            res[i + 2][0] += matrix_vector[2][0];
            res[i + 3][0] += matrix_vector[3][0];
        }

        res[0][0] = U_0;
        res[n * 3][0] += A * dU_dX_at_100;

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
        res[0][1] = 0;

        Matrix res_matrix(res);
        return res_matrix;
    }

    if (matrix_vector[0].size() == 4) {
        std::vector<std::vector<double>> res(n * 3 + 1, std::vector<double>(n * 3 + 1, 0));
        for (int i = 0; i < n * 3; i += 3) {
            for (int j = 0; j < 4; ++j) {
                res[i + j][i] += matrix_vector[j][0];
                res[i + j][i + 1] += matrix_vector[j][1];
                res[i + j][i + 2] += matrix_vector[j][2];
                res[i + j][i + 3] += matrix_vector[j][3];
            }
        }

        res[0][0] = 1;
        res[0][1] = 0;
        res[0][2] = 0;
        res[0][3] = 0;

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

double rsme(const std::vector<double> &lhs, const std::vector<double> &rhs) {
    if (rhs.size() != lhs.size()) {
        throw std::range_error("Can't count error due to sizes");
    }

    double error = 0.0;
    for (int i = 0; i < lhs.size(); ++i) {
        error += (lhs[i] - rhs[i]) * (lhs[i] - rhs[i]);
    }

    return sqrt(error / (double) lhs.size());

}

void print_vector(const std::vector<double> &vector) {
    for (auto elem: vector) {
        std::cout << elem << " ";
    }

    std::cout << std::endl;
}

std::vector<double> analytical_solve(const std::vector<double> &xes) {
    std::vector<double> res;
    res.reserve(xes.size());
    for (auto x: xes) {
        res.push_back(11. + 450. / (49. * exp(140.)) - (450. * exp(7. * x / 5.)) / (49. * exp(140.)) + 20. * x / 7.);
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

std::vector<double> cubic_solve(int n, double L) {
    Matrix quad_matrix({{37. / (10. * L),   -189. / (40. * L), 27. / (20. * L),   -13 / (40. * L)},
                        {-189. / (40. * L), 54. / (5. * L),    -297. / (40. * L), 27. / (20. * L)},
                        {27. / (20. * L),   -297. / (40. * L), 54. / (5. * L),    -189. / (40. * L)},
                        {-13 / (40. * L),   27. / (20. * L),   -189. / (40. * L), 37. / (10. * L)}});

    Matrix elementary_matrix({{-1. / 2.,   57. / 80.,  -3. / 10.,  7. / 80.},
                              {-57. / 80., 0.,         81. / 80.,  -3. / 10,},
                              {3. / 10.,   -81. / 80., 0.,         57. / 80.},
                              {-7. / 80.,  3. / 10.,   -57. / 80., 1. / 2.}});

    Matrix equality_matrix({{L / 8.},
                            {3. * L / 8.},
                            {3. * L / 8.},
                            {L / 8.}});

    auto left = quad_matrix * A + elementary_matrix * B;
    equality_matrix *= D;

    left = assemble(left, n);
    equality_matrix = assemble(equality_matrix, n);

    auto res = left.solve_gauss(equality_matrix);
    std::vector<double> reduced_res(n);
    for (int i = 0, j = 0; i < n; ++i, j += 3) {
        reduced_res[i] = res[j];
    }

    return reduced_res;
};


void save_to_file(const std::string &filename, const std::vector<double> &xes, const std::vector<double> &res) {
    std::ofstream file;
    file.open("data/" + filename);
    for (int i = 0; i < xes.size(); ++i) {
        file << xes[i] << " " << res[i] << std::endl;
    }

    file.close();
}
