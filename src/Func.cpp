//
// Created by Yuanyuan Qi on 6/26/24.
//

#include "Func.h"

std::vector<std::pair<real,real> > helper;

real log_eps(real x, real eps, int s_n) {
    if (x < eps) {
        real return_val = log(eps);
        real iter_val = (eps - x) / eps;
        for (int i=1; i < s_n; i++) {
            return_val -= pow(iter_val , i) / i;
        }
        return return_val;
    }
    return log(x);
}

