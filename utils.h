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

std::vector<double> analytical_solve(int n, double L);

std::vector<double> linear_solve(int n, double L);

std::vector<double> error_calc(const std::vector<double> &lhs, const std::vector<double> &rhs);

void print_vector(const std::vector<double> &vector);

#endif //MKE_UTILS_H
