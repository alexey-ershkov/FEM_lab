#include "utils.h"
#include "constants.h"


int main() {
    int n = 20;
    double L = (END - START) / (n - 1.);

    auto analytical_res = analytical_solve(n, L);
    auto linear_res = linear_solve(n, L);
    auto linear_error = error_calc(analytical_res, linear_res);

    print_vector(linear_error);

    return 0;
}
