#include "utils.h"
#include "constants.h"
#include "FEM_model.h"

int main() {
    system("rm -rf data/");
    system("mkdir data");
    system("rm -rf ../graphs/");
    system("mkdir ../graphs");

    auto fems = {5, 20, 40};
    double rsme_wanted = 0;
    double error_wanted = 0;
    for (auto fem_num: fems) {
        std::vector<double> x_nodes(fem_num + 1);
        for (int i = 0; i <= fem_num; i++) {
            x_nodes[i] = START + (double) i * (END - START) / fem_num;
        }

        auto linear_model = FEM_model(Linear, fem_num);
        auto linear_res = linear_model.solve();

        auto cubic_model = FEM_model(Cubic, fem_num);
        auto cubic_res = cubic_model.solve();

        auto analytical_res = analytical_solve(x_nodes);

        auto linear_error = error_calc(linear_res, analytical_res);
        auto cubic_error = error_calc(cubic_res, analytical_res);

        std::cout << "Linear RSME for " << fem_num << " fem_num is " << rsme(linear_res, analytical_res) << std::endl;
        std::cout << "Linear MAX ERROR for " << fem_num << " fem_num is "
                  << *std::max_element(linear_error.begin(), linear_error.end()) << std::endl;
        std::cout << "Cubic RSME for " << fem_num << " fem_num is " << rsme(cubic_res, analytical_res) << std::endl;
        std::cout << "Cubic MAX ERROR for " << fem_num << " fem_num is "
                  << *std::max_element(cubic_error.begin(), cubic_error.end()) << std::endl;

        if (fem_num == 20) {
            rsme_wanted = rsme(cubic_res, analytical_res);
            error_wanted = *std::max_element(cubic_error.begin(), cubic_error.end());
        }

        save_to_file("analytic" + std::to_string(fem_num), x_nodes, analytical_res);
        save_to_file("linear" + std::to_string(fem_num), x_nodes, linear_res);
        save_to_file("cubic" + std::to_string(fem_num), x_nodes, cubic_res);

        save_to_file("linear_error" + std::to_string(fem_num), x_nodes, linear_error);
        save_to_file("cubic_error" + std::to_string(fem_num), x_nodes, cubic_error);

        std::cout << std::endl;
    }

    system("gnuplot ../scripts/plot_solutions.p");
    system("gnuplot ../scripts/plot_errors.p");

    std::cout << "Wanted Linear for Cubic RSME " << rsme_wanted << " at fems 20 is "
              << find_linear_rsme_for_cubic_rsme(rsme_wanted) << std::endl;
    std::cout << "Wanted Linear for Cubic MAX ERROR " << error_wanted << " at fems 20 is "
              << find_linear_error_for_cubic_error(error_wanted) << std::endl;

    return 0;
}
