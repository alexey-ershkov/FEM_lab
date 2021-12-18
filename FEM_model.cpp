//
// Created by a.ershkov on 12.12.2021.
//

#include "FEM_model.h"

FEM_model::FEM_model(Type model_type, int fems) : left({}), right({}) {
    double L = (double) (END - START) / fems;
    type = model_type;

    if (type == Linear) {

        left = -Matrix({{1. / L,  -1. / L},
                        {-1. / L, 1. / L}}) * A + Matrix({{-1. / 2., 1. / 2.},
                                                          {-1. / 2., 1. / 2.}}) * B;
        right = -Matrix({{L / 2.},
                         {L / 2.}}) * D;
        nodes_num = fems + 1;

        left = assemble(left, nodes_num);
        right = assemble(right, nodes_num);

        return;
    }

    if (type == Cubic) {

        left = -Matrix({{37. / (10. * L),   -189. / (40. * L), 27. / (20. * L),   -13 / (40. * L)},
                        {-189. / (40. * L), 54. / (5. * L),    -297. / (40. * L), 27. / (20. * L)},
                        {27. / (20. * L),   -297. / (40. * L), 54. / (5. * L),    -189. / (40. * L)},
                        {-13 / (40. * L),   27. / (20. * L),   -189. / (40. * L), 37. / (10. * L)}}) * A +
               Matrix({{-1. / 2.,   57. / 80.,  -3. / 10.,  7. / 80.},
                       {-57. / 80., 0.,         81. / 80.,  -3. / 10,},
                       {3. / 10.,   -81. / 80., 0.,         57. / 80.},
                       {-7. / 80.,  3. / 10.,   -57. / 80., 1. / 2.}}) * B;

        right = -Matrix({{L / 8.},
                         {3. * L / 8.},
                         {3. * L / 8.},
                         {L / 8.}}) * D;
        nodes_num = 3 * fems + 1;

        left = assemble(left, fems);
        right = assemble(right, fems);

        return;
    }


    throw std::range_error("Can't init model");
}

std::vector<double> FEM_model::solve() {
    if (type == Linear) {
        return left.solve_gauss(right);
    }

    if (type == Cubic) {
        auto res = left.solve_gauss(right);
        std::vector<double> reduced_res;
        for (int i = 0; i < nodes_num; i += 3) {
            reduced_res.push_back(res[i]);
        }

        return reduced_res;
    }
}

Matrix FEM_model::assemble(Matrix to_assemble, int n) {
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

void FEM_model::print_matrix() {
    int i = 0;
    for (const auto &line: left.getMatrix()) {
        for (auto elem: line) {
            std::cout << elem << " ";
        }
        std::cout << right.getMatrix()[i][0] << std::endl;
        ++i;
    }
}
