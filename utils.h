//
// Created by a.ershkov on 12.12.2021.
//

#ifndef MKE_UTILS_H
#define MKE_UTILS_H

#include "Matrix.h"
#include <vector>
#include <exception>
#include <string>
#include <cmath>
#include <iostream>

std::vector<double> analytical_solve(const std::vector<double> &xes);

std::vector<double> error_calc(const std::vector<double> &lhs, const std::vector<double> &rhs);

double rsme(const std::vector<double> &lhs, const std::vector<double> &rhs);

void print_vector(const std::vector<double> &vector);

void save_to_file(const std::string &filename, const std::vector<double> &xes, const std::vector<double> &res);

double find_linear_rsme_for_cubic_rsme(double rsme_wanted);
double find_linear_error_for_cubic_error(double max_error);


#endif //MKE_UTILS_H
