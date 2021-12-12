#include "utils.h"
#include "constants.h"


int main() {
    int fems = 20;
    int nodes = fems + 1;
    double L = (END - START) / fems;

    std::vector<double> x_nodes(nodes);
    for (int i = 0; i < nodes; i++) {
        x_nodes[i] = START + i * (END - START) / fems;
    }


    system("rm -rf data/");
    system("mkdir data");

    auto analytical_res = analytical_solve(x_nodes);
    auto linear_res = linear_solve(fems+1, L);
    auto cubic_res = cubic_solve(fems+1, L);

    std::cout << rsme(linear_res, analytical_res) << std::endl;
    std::cout << rsme(cubic_res, analytical_res) << std::endl;

    save_to_file("linear", x_nodes, linear_res);
    save_to_file("cubic", x_nodes, cubic_res);

    return 0;
}

