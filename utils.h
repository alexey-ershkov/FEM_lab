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

Matrix assemble(Matrix to_assemble, int n);

std::vector<double> analytical_solve(const std::vector<double>& xes);

std::vector<double> linear_solve(int n, double L);

std::vector<double> cubic_solve(int n, double L);

std::vector<double> error_calc(const std::vector<double> &lhs, const std::vector<double> &rhs);

double rsme(const std::vector<double> &lhs, const std::vector<double> &rhs);

void print_vector(const std::vector<double> &vector);

void save_to_file(const std::string &filename, const std::vector<double> &xes, const std::vector<double> &res);

#endif //MKE_UTILS_H
