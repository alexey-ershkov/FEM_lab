//
// Created by a.ershkov on 12.12.2021.
//

#ifndef MKE_FEM_MODEL_H
#define MKE_FEM_MODEL_H


#include "Matrix.h"
#include <string>
#include "constants.h"
#include <exception>
#include <iostream>

enum Type {
    Linear,
    Cubic
};

class FEM_model {
public:
    FEM_model(Type model_type, int fems);

    std::vector<double> solve();

    void print_matrix();

private:
    Matrix assemble(Matrix to_assemble, int n);

    int nodes_num;
    Type type;
    Matrix left;
    Matrix right;
};


#endif //MKE_FEM_MODEL_H
