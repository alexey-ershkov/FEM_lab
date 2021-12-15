//
// Created by a.ershkov on 12.12.2021.
//

#include "utils.h"
#include "FEM_model.h"
#include <fstream>

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


void save_to_file(const std::string &filename, const std::vector<double> &xes, const std::vector<double> &res) {
    std::ofstream file;
    file.open("data/" + filename);
    for (int i = 0; i < xes.size(); ++i) {
        file << xes[i] << " " << res[i] << std::endl;
    }

    file.close();
}

double find_linear_rsme_for_cubic_rsme(double rsme_wanted) {
    int fem_num = 1;
    double current_rsme = 100;

    while (rsme_wanted < current_rsme) {
        ++fem_num;
        std::vector<double> x_nodes(fem_num + 1);
        for (int i = 0; i <= fem_num; i++) {
            x_nodes[i] = START + (double) i * (END - START) / fem_num;
        }

        auto linear_model = FEM_model(Linear, fem_num);
        auto linear_res = linear_model.solve();

        auto analytical_res = analytical_solve(x_nodes);

        current_rsme = rsme(linear_res, analytical_res);

    }

    return fem_num;
}

double find_linear_error_for_cubic_error(double max_error) {
    int fem_num = 1;
    double current_error = 100;

    while (max_error < current_error) {
        ++fem_num;
        std::vector<double> x_nodes(fem_num + 1);
        for (int i = 0; i <= fem_num; i++) {
            x_nodes[i] = START + (double) i * (END - START) / fem_num;
        }

        auto linear_model = FEM_model(Linear, fem_num);
        auto linear_res = linear_model.solve();

        auto analytical_res = analytical_solve(x_nodes);

        auto linear_error = error_calc(linear_res, analytical_res);
        current_error = *std::max_element(linear_error.begin(), linear_error.end());

    }

    return fem_num;
}